
#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "atom.h"
#include "force.h"
#include "neighbor.h"
#include "thermo.h"

class Integrate {
public:
  struct TagInitialIntegrate {};
  struct TagFinalIntegrate {};

  float dt;
  float dtforce;
  int ntimes;
  int nlocal, nmax;
  Buffer<float3> x, v, f, xold;
  float mass;

  int sort_every;

  Integrate();
  ~Integrate();
  void setup();
  void initialIntegrate(Stream &stream, Device &device);
  void finalIntegrate(Stream &stream, Device &device);
  void run(Stream &stream, Device &device, Atom &, Force *, Neighbor &,
           Thermo &);
};
#endif