
#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "atom.h"
#include "comm.h"
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
  // Buffer<float3> x, v, f, xold;
  float mass;

  int sort_every;

  Shader<1> initial_integrate;
  Shader<1> final_integrate;

  Integrate();
  ~Integrate();
  void setup();
  void setup_shader(Device &device, Atom &atom, Force *force);
  void run(Stream &stream, Atom &, Force *, Comm &, Neighbor &, Thermo &);
};
#endif