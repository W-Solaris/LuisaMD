#ifndef LJ_H_
#define LJ_H_

#include "types.h"

struct In {
  int nx, ny, nz;
  md_float t_request;
  md_float rho;
  int units;
  ForceStyle forcetype;
  md_float epsilon, sigma;
  char *datafile;
  int ntimes;
  md_float dt;
  int neigh_every;
  md_float force_cut;
  md_float neigh_cut;
  int thermo_nstat;
};

#endif