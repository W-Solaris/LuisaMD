#ifndef COMM_H
#define COMM_H

#include "atom.h"

class Comm {
public:
  Comm();
  ~Comm();
  int setup(Stream &, Device &, float, Atom &);
  void setup_shaders(Device &, Atom &);
  void communicate(Stream &);
  void borders(Stream &);

public:
  int me;    // my proc ID
  int nswap; // # of swaps to perform

  Buffer<float> pbc_flagx; // PBC correction in x for this swap
  Buffer<float> pbc_flagy; // same in y
  Buffer<float> pbc_flagz; // same in z

  int sendnum[6], recvnum[6]; // # of atoms to send/recv in each swap

  int firstrecv[6]; // where to put 1st recv atom in each swap
  int *maxsendlist;

  Buffer<int> ghost_num;
  Image<int> ghost_atoms;

  int need[3]; // how many procs away needed in each dim
  std::vector<float> h_slablo, h_slabhi;
  Buffer<float> slablo, slabhi;
  Shader<1> communicate_shader;
  Shader<1> borders_shader;
  Shader<1> nghost_shader;
};

#endif
