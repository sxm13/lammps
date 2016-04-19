/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/*
  Generic pseudo random number generator (pRNG) wrapper for LAMMPS.
  Contributed by Axel Kohlmeyer <akohlmey@gmail.com>
  Park/Miller and Masaglia pRNG code by Steve Plimpton
  PCG pRNG code by Stefan Paquay and Maarten Sebregts
  Mersenne Twister pRNG code by Axel Kohlmeyer

  Basic design:
  - individual pRNGs are provided as static functions that store
    their state in an opaque void pointer.
  - the pRNGs has to provide a function to initialize its state from a seed
    and a function to return a uniformly distributed 0.0 < pRNG <= 1.0
  - the state is described in a struct and all structs have the same
    first section to access state flags common to all pRNGs
  - the wrapper class accesses the uniform pRNG through a function pointer,
    which is set during initialization
  - gaussian distributed random numbers are produces from either the polar
    method or the ziggurat method
  - default settings for current pRNG method and seed are kept in Update;
    the "random" command can be used to change those
  - generic serialization and deserialization

  How to add a new pRNG to this class:
  - add a new flag to the enumerator in the header (note they are 2^n)
  - add a typedef to a suitable struct for a rng_XXX_t type.
  - add static rng_XXX_init() and rng_XXX_uniform() functions
  - add the required else if hooks to call the pRNG specific functions
 */

#include <math.h>
#include <string.h>

#include "random.h"
#include "update.h"
#include "comm.h"
#include "error.h"

using namespace LAMMPS_NS;

/* common constants */

static const int zigg_layers = 256;
static const double ran_uint32_scale = 1.0 / (256.0*256.0*256.0*256.0);

// keyword names. must be in sync with enumerator in header.
static const char * const ran_names[] = {
    "(unknown)", "equal", "ziggurat", "polar",
    "park", "mars", "pcg", "mt", "" };

/* common elements of all pRNG data structs */

typedef struct {
  int style;
  int save;
  double second;
  double xlayers[zigg_layers];
  double ylayers[zigg_layers];
} ran_any_t;


/* ---------------------------------------------------------------------- */

/* marsaglia pRNG after RANMAR in F James, Comp Phys Comm, 60, 329 (1990) */

typedef struct {
  int style;
  int save;
  double second;
  double xlayers[zigg_layers];
  double ylayers[zigg_layers];
  int i97,j97;
  double c,cd,cm;
  double u[97+1];
} ran_mars_t;

static void ran_mars_init(int seed, void *ptr);
static double ran_mars_uniform(void *ptr);

/* ---------------------------------------------------------------------- */

/* Park/Miller pRNG. Do not use. Only provided for backward compatibility.  */

typedef struct {
  int style;
  int save;
  double second;
  double xlayers[zigg_layers];
  double ylayers[zigg_layers];
  int seed;
} ran_park_t;

static void ran_park_init(int seed, void *ptr);
static double ran_park_uniform(void *ptr);

/* ----------------------------------------------------------------------
   PCG random number generator

   Paper:
   - http://www.pcg-random.org/pdf/toms-oneill-pcg-family-v1.02.pdf
   - https://web.archive.org/web/20150817092619/http://www.pcg-random.org/pdf/toms-oneill-pcg-family-v1.02.pdf
   According to website submitted to ACM Transactions on Mathematical
   Software and as of August 2015 under review.

   Code contributed by Stefan Paquay and Maarten Sebregts
   (Eindhoven University of Technology).
------------------------------------------------------------------------ */

typedef struct {
  int style;
  int save;
  double second;
  double xlayers[zigg_layers];
  double ylayers[zigg_layers];
  uint64_t state;
  uint64_t inc;
} ran_pcg_t;

static void ran_pcg_init(int seed, int stream, void *ptr);
static double ran_pcg_uniform(void *ptr);

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Mersenne Twister (MT19937) pseudo random number generator:
   M. Matsumoto & T. Nishimura,
   ACM Transactions on Modeling and Computer Simulation,
   Vol. 8, No. 1, 1998, pp. 3-30.

   Uses the Marsaglia pRNG from above to generate the initial seeds
------------------------------------------------------------------------ */

static const int MT_N = 624;
typedef struct {
  int style;
  int save;
  double second;
  double xlayers[zigg_layers];
  double ylayers[zigg_layers];
  int idx;
  uint32_t m[MT_N];
} ran_mt_t;

static void ran_mt_init(int seed, void *ptr);
static double ran_mt_uniform(void *ptr);

/* ---------------------------------------------------------------------- */

Random::Random(LAMMPS *lmp, int seed) : Pointers(lmp)
{
  _uniform = 0;
  rng_style = update->rng_style;
  setup_rng(seed);
}

/* ---------------------------------------------------------------------- */

Random::Random(LAMMPS *lmp, int seed, int style) : Pointers(lmp)
{
  _uniform = 0;
  rng_style = style;
  setup_rng(seed);
}

/* ---------------------------------------------------------------------- */

Random::~Random()
{
  _uniform = 0;
  clear_rng();
  state = 0;
}

/* ----------------------------------------------------------------------
   perform common setup tasks for pRNGs
------------------------------------------------------------------------- */

void Random::setup_rng(int seed)
{
  int tmp_seed = seed;

  if (seed == 0) tmp_seed = update->rng_seed;
  if (seed < 0) error->one(FLERR,"Invalid seed for random number generator");

  if (rng_style & RNG_PARK) {
    _uniform = &ran_park_uniform;
    state = new ran_park_t;
    nbytes = sizeof(ran_park_t);
  } else if (rng_style & RNG_MARS) {
    _uniform = &ran_mars_uniform;
    state = new ran_mars_t;
    nbytes = sizeof(ran_mars_t);
  } else if (rng_style & RNG_PCG) {
    _uniform = &ran_pcg_uniform;
    state = new ran_pcg_t;
    nbytes = sizeof(ran_pcg_t);
  } else if (rng_style & RNG_MT) {
    _uniform = &ran_mt_uniform;
    state = new ran_mt_t;
    nbytes = sizeof(ran_mt_t);
  } else error->one(FLERR,"Unsupported pRNG style");

  // the first bytes of the state struct is the rng_style setting
  // copy it over as it is required when deserializing the pRNG.
  memcpy(state,&rng_style,sizeof(int));

  if (rng_style & RNG_EQUAL) init(tmp_seed);
  else init(tmp_seed + comm->me);

  // preparation for gaussian distributed random numbers
  ran_any_t &rng = *static_cast<ran_any_t *>(state);
  const double xn = 3.6541528853610088;
  const double A  = 0.00492867323399;

  rng.save = 0;
  rng.second = 0.0;

  rng.xlayers[0] = xn;
  rng.ylayers[0] = exp(-0.5*xn*xn);
  for (int i=1; i < zigg_layers; ++i) {
    rng.ylayers[i] = A/rng.xlayers[i-1] + rng.ylayers[i-1];
    rng.xlayers[i] = sqrt(-2.0*log(rng.ylayers[i]));
  }
  rng.xlayers[zigg_layers-1] = 0.0;
  rng.ylayers[zigg_layers-1] = 1.0;
}

/* ----------------------------------------------------------------------
   initialize the state of pRNG based on provided seed.
------------------------------------------------------------------------- */

void Random::init(int seed)
{
  if (rng_style & RNG_PARK) {
    ran_park_init(seed,state);
  } else if (rng_style & RNG_MARS) {
    ran_mars_init(seed,state);
  } else if (rng_style & RNG_PCG) {
    int stream;
    if (rng_style & RNG_EQUAL) stream = 0;
    else stream = comm->me;
    ran_pcg_init(seed,stream,state);
  } else if (rng_style & RNG_MT) {
    ran_mt_init(seed,state);
  } else error->one(FLERR,"Unsupported pRNG style");
}

/* ----------------------------------------------------------------------
   initialize the pRNG based on coordinate triple and provided seed.
   use hash function, treating user seed and coords as sequence of input bytes
   this is Jenkins One-at-a-time hash, see Wikipedia entry on hash tables
------------------------------------------------------------------------- */

void Random::init(int seed, double *coord)
{
  int i;

  char *str = (char *) &seed;
  int n = sizeof(int);

  unsigned int hash = 0;
  for (i = 0; i < n; i++) {
    hash += str[i];
    hash += (hash << 10);
    hash ^= (hash >> 6);
  }

  str = (char *) coord;
  n = 3 * sizeof(double);
  for (i = 0; i < n; i++) {
    hash += str[i];
    hash += (hash << 10);
    hash ^= (hash >> 6);
  }

  hash += (hash << 3);
  hash ^= (hash >> 11);
  hash += (hash << 15);

  // keep 31 bits of unsigned int as new seed
  // do not allow seed = 0, since will cause hang in gaussian()

  seed = hash & 0x7ffffff;
  if (!seed) seed = 1;
  init(seed);

  // warm up pRNG a little bit
  for (i = 0; i < 5; i++) uniform();
}

/* ----------------------------------------------------------------------
   perform common cleanup tasks for pRNGs
------------------------------------------------------------------------- */

void Random::clear_rng()
{
  if (rng_style & RNG_PARK) {
    ran_park_t *tmp = static_cast<ran_park_t *>(state);
    delete tmp;
  } else if (rng_style & RNG_MARS) {
    ran_mars_t *tmp = static_cast<ran_mars_t *>(state);
    delete tmp;
  } else if (rng_style & RNG_PCG) {
    ran_pcg_t *tmp = static_cast<ran_pcg_t *>(state);
    delete tmp;
  } else if (rng_style & RNG_MT) {
    ran_mt_t *tmp = static_cast<ran_mt_t *>(state);
    delete tmp;
  }
}

/* ---------------------------------------------------------------------- */

// convert rng flag to option name
const char *Random::rng2name(int flag) const
{
    for (int i = 0; 1<<i < RNG_LAST; ++i) {
        if (1<<i == flag) return ran_names[i];
    }
    return ran_names[0];
}

// convert option name to rng flag
int Random::name2rng(const char *option) const
{
    for (int i = 0; 1<<i < RNG_LAST; ++i) {
        if (strcmp(option,ran_names[i]) == 0) return 1<<i;
    }
    return 0;
}

// make flags and state of current pRNG available
int Random::get_state(char **ptr)
{
    if (ptr == 0) return 0;
    *ptr = (char *) state;
    return nbytes;
}

// set type flags and state of this pRNG to provided data
// free the current pRNG and its data, initialize a new one
// and then copy the entire state over.
int Random::set_state(char *ptr)
{
    if (ptr == 0) return 0;
    ran_any_t &rng = *(ran_any_t *)ptr;
    clear_rng();
    rng_style = rng.style;
    setup_rng(0);
    memcpy(state,ptr,nbytes);
    return 1;
}

/* -------------------------------------------------------------------------
   Generate gaussian distribution from uniform via Mersenne polar method
------------------------------------------------------------------------- */

double Random::gauss_polar()
{
  ran_any_t &rng = *static_cast<ran_any_t *>(state);
  double first,v1,v2,rsq,fac;

  if (!rng.save) {
    int again = 1;
    while (again) {
      v1 = 2.0*uniform()-1.0;
      v2 = 2.0*uniform()-1.0;
      rsq = v1*v1 + v2*v2;
      if (rsq < 1.0 && rsq != 0.0) again = 0;
    }
    fac = sqrt(-2.0*log(rsq)/rsq);
    rng.second = v1*fac;
    first = v2*fac;
    rng.save = 1;
  } else {
    first = rng.second;
    rng.save = 0;
  }
  return first;
}

/* -------------------------------------------------------------------------
   Generate gaussian distribution from uniform via ziggurat method
------------------------------------------------------------------------- */

double Random::gauss_zigg()
{
  ran_any_t &rng = *static_cast<ran_any_t *>(state);
  const double A  = 0.00492867323399;
  double x,y;

  int n = static_cast<int>(static_cast<double>(zigg_layers)*uniform());
  if (n == 0) {                 // special case: tail
    const double xbound = A / rng.ylayers[0];
    x = 2.0*xbound * uniform() - xbound;
    if (fabs(x) < rng.xlayers[0]) {
      return x;
    } else {
      double rv;
      do {
        const double x0 = rng.xlayers[0];
        x = -log(uniform()) / x0;
        y = -log(uniform());
        rv = (uniform() < 0.5) ? -x0-x : x0+x;
      } while (2*y < x*x);
      return rv;
    }
  } else if (n == zigg_layers-1) {  // special case: top
    const double xbound = rng.xlayers[n-1];
    x = 2.0*xbound * uniform() - xbound;
    const double delta = rng.ylayers[n]-rng.ylayers[n-1];
    y = delta*uniform() + rng.ylayers[n-1];

    if (y < exp(-0.5*x*x))
      return x;
    else
      return gauss_zigg();
  } else {                      // normal case
    const double xbound = rng.xlayers[n];
    x = 2.0*xbound * uniform() - xbound;
    if (fabs(x) < rng.xlayers[n+1])
      return x;

    const double delta = rng.ylayers[n+1]-rng.ylayers[n];
    y = delta*uniform() + rng.ylayers[n];
    if (y < exp(-0.5*x*x))
      return x;
    else
      return gauss_zigg();
  }
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

void ran_mars_init(int seed, void *ptr)
{
  int ij,kl,i,j,k,l,ii,jj,m;
  double s,t;

  ran_mars_t &state = *static_cast<ran_mars_t *>(ptr);

  // sanitize seed. assumes the calling code ensures it to be > 0
  seed %= 900000000;

  ij = (seed-1)/30082;
  kl = (seed-1) - 30082*ij;
  i = (ij/177) % 177 + 2;
  j = ij %177 + 2;
  k = (kl/169) % 178 + 1;
  l = kl % 169;
  for (ii = 1; ii <= 97; ii++) {
    s = 0.0;
    t = 0.5;
    for (jj = 1; jj <= 24; jj++) {
      m = ((i*j) % 179)*k % 179;
      i = j;
      j = k;
      k = m;
      l = (53*l+1) % 169;
      if ((l*m) % 64 >= 32) s = s + t;
      t = 0.5*t;
    }
    state.u[ii] = s;
  }
  state.c = 362436.0 / 16777216.0;
  state.cd = 7654321.0 / 16777216.0;
  state.cm = 16777213.0 / 16777216.0;
  state.i97 = 97;
  state.j97 = 33;
  ran_mars_uniform(ptr);
}

double ran_mars_uniform(void *ptr)
{
  ran_mars_t &state = *static_cast<ran_mars_t *>(ptr);
  double uni = state.u[state.i97] - state.u[state.j97];

  if (uni < 0.0) uni += 1.0;
  state.u[state.i97] = uni;
  state.i97--;
  if (state.i97 == 0) state.i97 = 97;
  state.j97--;
  if (state.j97 == 0) state.j97 = 97;
  state.c -= state.cd;
  if (state.c < 0.0) state.c += state.cm;
  uni -= state.c;
  if (uni < 0.0) uni += 1.0;
  return uni;
}

/* ---------------------------------------------------------------------- */

static const int IA = 16807;
static const int IM = 2147483647;
static const double AM = (1.0/IM);
static const int IQ = 127773;
static const int IR =  2836;

void ran_park_init(int seed, void *ptr)
{
  ran_park_t &state = *static_cast<ran_park_t *>(ptr);
  state.seed = seed;
}

double ran_park_uniform(void *ptr)
{
  ran_park_t &state = *static_cast<ran_park_t *>(ptr);
  const int k = state.seed/IQ;

  state.seed = IA*(state.seed-k*IQ) - IR*k;
  if (state.seed < 0) state.seed += IM;
  return AM*state.seed;
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   32-bits uniform RN
------------------------------------------------------------------------- */

static inline uint32_t pcg_rotate32(uint32_t v, uint32_t r)
{
  // rotate v by r bits
  // https://en.wikipedia.org/wiki/Bitwise_operation#Circular_rotates_in_C_and_C.2B.2B
  return (v >> r) | (v << ((-r) & 31));
}

static uint32_t ran_pcg_uniformi(void *ptr)
{
  ran_pcg_t &state = *static_cast<ran_pcg_t *>(ptr);
  uint64_t oldstate = state.state;

  // advance state, constant taken from implementation by Melissa O'Neill:
  // https://github.com/imneme/pcg-c/blob/master/include/pcg_variants.h#L253
  // (PCG_DEFAULT_MULTIPLIER_64)
  state.state = oldstate * 6364136223846793005ULL + state.inc;
  // calculate output, section 6.3, PCG-XSH-RR
  return pcg_rotate32((oldstate ^ (oldstate >> 18)) >> 27, oldstate >> 59);
}

void ran_pcg_init(int seed, int stream, void *ptr)
{
  ran_pcg_t &state = *static_cast<ran_pcg_t *>(ptr);
  state.state = seed;
  state.inc = 2*stream + 1;
  ran_pcg_uniformi(ptr);
  state.state += seed;
  ran_pcg_uniformi(ptr);
}

double ran_pcg_uniform(void *ptr)
{
  return ran_uint32_scale * static_cast<double>(ran_pcg_uniformi(ptr));
}

/* ---------------------------------------------------------------------- */

static const int MT_S=7;
static const int MT_U=11;
static const int MT_T=15;
static const int MT_L=18;
static const int MT_R=31;
static const int MT_M=397;
static const uint32_t MT_A = 0x9908B0DFU;
static const uint32_t MT_B = 0x9D2C5680U;
static const uint32_t MT_C = 0xEFC60000U;

static uint32_t ran_mt_randomize(void *ptr)
{
  ran_mt_t &state = *static_cast<ran_mt_t *>(ptr);
  uint32_t r;

  if (state.idx >= MT_N) {
    // fill the entire status array with new data in one sweep
    const uint32_t LMASK = (1LU << MT_R) - 1;  // Lower MT_R bits
    const uint32_t UMASK = 0xFFFFFFFF << MT_R; // Upper (32 - MT_R) bits
    static const uint32_t magic[2] = {0, MT_A};
    const int diff = MT_N-MT_M;
    int i;

    for (i=0; i < diff; ++i) {
      r = (state.m[i] & UMASK) | (state.m[i+1] & LMASK);
      state.m[i] = state.m[i+MT_M] ^ (r >> 1) ^ magic[r & 1];}

    for (i=diff; i < MT_N-1; ++i) {
      r = (state.m[i] & UMASK) | (state.m[i+1] & LMASK);
      state.m[i] = state.m[i-diff] ^ (r >> 1) ^ magic[r & 1];}

    r = (state.m[MT_N-1] & UMASK) | (state.m[0] & LMASK);
    state.m[MT_N-1] = state.m[MT_M-1] ^ (r >> 1) ^ magic[r & 1];
    state.idx = 0;
  }
  r = state.m[state.idx++];

  r ^=  r >> MT_U;
  r ^= (r << MT_S) & MT_B;
  r ^= (r << MT_T) & MT_C;
  r ^=  r >> MT_L;

  return r;
}

void ran_mt_init(int seed, void *ptr)
{
  ran_mt_t &state = *static_cast<ran_mt_t *>(ptr);
  const uint32_t f = 1812433253UL;
  int i;

  // start minimal initialization
  state.m[0] = seed;
  state.idx = MT_N-1;
  for (i=1; i < MT_N; ++i)
    state.m[i] = (f * (state.m[i-1] ^ (state.m[i-1] >> 30)) + i);

  // to seed the RNG some more we need a second RNG
  ran_mars_t *mars = new ran_mars_t;
  ran_mars_init(seed,mars);
  for (i=0; i < MT_N-1; ++i) {
    state.m[i+1] = (state.m[i+1]^((state.m[i]^(state.m[i]>>30))*1664525UL))
      + (uint32_t) (ran_mars_uniform(mars)* (1U<<31)) + i;
  }
  state.m[0] = state.m[MT_N-1];
  delete mars;

  for (i=0; i < MT_N-1; ++i) {
    state.m[i+1] = (state.m[i+1] ^ ((state.m[i] ^ (state.m[i] >> 30))
                                    * 1566083941UL))-i-1;
  }
  state.m[0] = 0x80000000UL;

  // randomize one more turn
  state.idx = 0;
  for (i=0; i < MT_N-1; ++i) ran_mt_randomize(&state);
}

double ran_mt_uniform(void *ptr)
{
  return ran_uint32_scale * static_cast<double>(ran_mt_randomize(ptr));
}

