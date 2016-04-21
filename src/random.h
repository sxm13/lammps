/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_RANDOM_H
#define LMP_RANDOM_H

#include "pointers.h"

namespace LAMMPS_NS {

class Random : protected Pointers {
 public:
  Random(class LAMMPS *, int);      // constructor with seed only
  Random(class LAMMPS *, int, int); // constructor with seed and style
  virtual ~Random();

  // return next uniformly distributed random number between 0.0 and 1.0
  double uniform() {
    if (_uniform) return (*_uniform)(state);
    else return 0.0;
  };

  // return next gaussian distributed random number
  double gaussian() {
    if (rng_style & RNG_ZIGG) return gauss_zigg();
    else return gauss_polar();
  };

  void init(int);
  void init(int, double *);

  int get_state(char **); // give access to pRNG state array and return length
  int set_state(char *);  // initialize pRNG from given state
  void read_state(const char *, int);  // read pRNG state from file
  void write_state(const char *, int); // write pRNG state to file

  enum {RNG_EQUAL=1<<0, // pRNG is initialized equally across all parallel tasks
        RNG_ZIGG =1<<1, // gaussian RNGs via ziggurat instead of marsaglia polar
        RNG_POLAR=1<<2, // gaussian RNGs via marsaglia polar instead of ziggurat
        RNG_GAUSS=(RNG_ZIGG|RNG_POLAR), // combined bits for gauss dist flags
        RNG_PARK=1<<3,  // Park/Miller pRNG (not recommanded)
        RNG_MARS=1<<4,  // RANMAR in F James, Comp Phys Comm, 60, 329 (1990)
        RNG_PCG=1<<5 ,  // PCG, a fast pRNG from http://www.pcg-random.org/
        RNG_MT=1<<6,    // portable Mersenne Twister MT19937 pRNG
        RNG_SARU=1<<7,  // SARU LCG RNG by S.P.Wolsky
        // combined bits for all pRNGs
        RNG_MASK=(RNG_PARK|RNG_MARS|RNG_PCG|RNG_MT|RNG_SARU),
        RNG_LAST=1<<7 };

  const char *rng2name(int) const;  // return name of pRNG from constant
  int name2rng(const char *) const; // find pRNG constant from name

 private:
  int rng_style;     // style settings for this pRNG instance
  int nbytes;        // length in bytes of the pRNG state array
  void *state;       // pRNG state data. cast to a struct for individual pRNGs.
  double (*_uniform)(void *);  // function pointer to lo-level pRNG call

 private:
  void setup_rng(int);  // do basic pRNG setup tasks
  void clear_rng();     // reset pRNG to pre-constructor state
  double gauss_polar();
  double gauss_zigg();
};

}

#endif

/* ERROR/WARNING messages:

E: Invalid seed for random # generator

The initial seed for this random number generator must be a positive
integer.

*/
