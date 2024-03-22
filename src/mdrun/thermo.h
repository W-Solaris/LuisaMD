#ifndef THERMO_H
#define THERMO_H

enum units { LJ, METAL };
#include "atom.h"
#include "force.h"
#include "neighbor.h"
#include "types.h"

class Integrate;

class Thermo {
public:
  int nstat;
  int ntimes;
  int nlocal;
  Buffer<int> mstat;
  Buffer<int> steparr;
  Buffer<float> tmparr;
  Buffer<float> engarr;
  Buffer<float> prsarr;

  std::vector<int> h_steparr;
  std::vector<float> h_tmparr;
  std::vector<float> h_engarr;
  std::vector<float> h_prsarr;

  Thermo();
  ~Thermo();
  void setup(Stream &stream, Device &device, float, Integrate &integrate,
             Atom &atom, int);
  void temperature(Stream &stream, Device &device, Atom &);

  //   float energy(Atom &, Neighbor &, Force *);

  //   float pressure(float, Force *);
  void compute(Stream &stream, Device &device, int, Atom &, Neighbor &,
               Force *);

  Buffer<float3> v;
  float mass;

  Buffer<float> t_act; // length = 1
  float t_scale, e_scale, p_scale, dof_boltz, mvv2e;

private:
  float rho;
};

#endif