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

/* ----------------------------------------------------------------------
   Author: Amin Aramoon (Johns Hopkins U) aaramoo1@jhu.edu
------------------------------------------------------------------------- */

#ifndef LMP_TYPE_DETECTOR_H
#define LMP_TYPE_DETECTOR_H

#include <stdlib.h>
#include <vector>
#include <ctype.h>

namespace LAMMPS_NS {

class TypeDetector {

  struct Entity {
    int *types; // list of atom types. * -> 0
    int value;
    int level;

    Entity() {
      types = NULL;
      value = -1;
      level = 0;
    }

    virtual ~Entity() {
      delete[] types;
      types = NULL;
    }
  };

 public:
  enum Style {BOND,ANGLE,DIHEDRAL,IMPROPER};

  int type;
  int length;

 private:
  TypeDetector() {}

 public:
  TypeDetector(enum Style t) : type(t),length(0) {
    if (type == BOND)
      length = 2;
    else if (type == ANGLE)
      length = 3;
    else if (type == DIHEDRAL)
      length = 4;
    else if (type == IMPROPER)
      length = 4;
  }
  virtual ~TypeDetector() {
    for (std::vector<Entity *>::iterator it = values.begin();
         it != values.end(); ++it)
      delete (*it);
  }

  bool init(const char* input) {
    std::vector<char*> args;
    int len = strlen(input) + 1;
    char* s = new char[len];
    strcpy(s, input);
    int ntype = tokenize(s, args);
    if (ntype < 0) return false;

    int *t = new int[length];
    int val = 0;
    val = atoi(args[0]);
    for (int num = 0; num < ntype; ++num) {
      const int j = num*(length+1);
      if (!isdigit(args[j][0])) return false;
      val=atoi(args[j]);
      for (int i = 1; i <= length; i++) {
        if (strcmp(args[j+i], "*") == 0)  t[i-1] = 0;
        else if (isdigit(args[j+i][0])) t[i-1] = atoi(args[j+i]);
        else return false;   // syntax error in type pattern
      }
      this->set(t, val);
    }
    delete[] t;
    delete[] s;
    return true;
  }

  int get(const int *t1) {
    for (std::vector<Entity *>::iterator it = values.begin();
         it != values.end(); it++) {
      if (equals((*it)->types, t1))
        return (*it)->value;
    }
    return -1;
  }

  bool check_types(int nentitytypes, int natomtypes) {
    int entitymax = 0;
    int atommax = -1;
    for (std::vector<Entity *>::iterator it = values.begin();
         it != values.end(); it++) {
      if ((*it)->value > entitymax) entitymax = (*it)->value;
      int *types = (*it)->types;
      for (int i = 0; i < length; ++i)
        if (types[i] > atommax) atommax = types[i];
    }
    if (entitymax < 1 || entitymax > nentitytypes) return false;
    if (atommax < 0 || atommax > natomtypes) return false;
    return true;
  }

 private:

  void set(const int *t1, int value) {
    int i, p = 0;
    for (i = 0; i < length; ++i)
      if (t1[i] == 0) ++p;
    std::vector<Entity *>::iterator it = values.begin();
    bool equal = false;
    while (it != values.end()) {
      Entity *e = (*it);
      if (equals(e->types, t1) && e->level == p) {
        equal = true;
        break;
      } else if (e->level > p) {
        break;
      }
      it++;
    }
    if (equal)
      (*it)->value = value;
    else {
      Entity *ne = new Entity();
      ne->types = new int[length];
      for (int i = 0; i < length; ++i)
        ne->types[i] = t1[i];
      ne->value = value;
      ne->level = p;
      values.insert(it, ne);
    }
  }

  int equals(const int *types, const int *id) {
    for (int i = 0; i < length; ++i)
      if (types[i] == id[i] || types[i] == 0) return 1;

    return 0;
  }

  int tokenize(char* dup, std::vector<char*> & v) {
    char *next = dup;
    int num = 0;
    do {
      if (dup) ++num;
      next = strchr(next,';');
      if (next != NULL) {
        next[0] = '\0';
        ++next;
      }
      char *token = strtok(dup, " \t\r\n");
      int len = 0;
      while (token != NULL) {
        ++len;
        v.push_back(token);
        token = strtok(NULL, " \t\r\n");
      }
      if (len != length+1) return -1;
      dup = next;
    } while (next != NULL);
    if (num < 1) return -1;
    return num;
  }

  std::vector<Entity*> values;
};

}

#endif /* SRC_TYPE_DETECTOR_H */
