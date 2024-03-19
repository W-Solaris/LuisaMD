#ifndef LJ_H_
#define LJ_H_

#include "types.h"

struct In {
  int nx, ny, nz;
  float t_request;
  float rho;
  int units;
  ForceStyle forcetype;
  float epsilon, sigma;
  char *datafile;
  int ntimes;
  float dt;
  int neigh_every;
  float force_cut;
  float neigh_cut;
  int thermo_nstat;
};

#endif