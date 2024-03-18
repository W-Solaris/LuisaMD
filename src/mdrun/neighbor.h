#ifndef NEIGHBOR_H
#define NEIGHBOR_H

#include "atom.h"
class Neighbor {
public:
  typedef int value_type;

  int every;                // re-neighbor every this often
  int nbinx, nbiny, nbinz;  // # of global bins
  md_float cutneigh;        // neighbor cutoff
  Buffer<float> cutneighsq; // neighbor cutoff squared
  md_float cutneighsq_stack[MAX_STACK_TYPES * MAX_STACK_TYPES];
  int ncalls;         // # of times build has been called
  int max_totalneigh; // largest # of neighbors ever stored

  Buffer<int> numneigh; // # of neighbors for each atom
  Image<int> neighbors; // array of neighbors of each atom
  int maxneighs;        // max number of neighbors per atom
  int halfneigh;
  int team_neigh_build;

  int ghost_newton;
  int count;
  Neighbor(Device device, int ntypes_);
  ~Neighbor();

  int setup(Stream &stream, Device device,
            Atom &); // setup bins based on box and cutoff
  void build(Stream &stream, Device device, Atom &); // create neighbor list

  //   Atom is going to call binatoms etc for sorting
  //   void binatoms(Atom &atom, int count = -1); // bin all atoms

  //   Buffer<int> bincount; // ptr to 1st atom in each bin
  //   Image<int> bins;      // ptr to next atom in each bin
  //   Buffer<int> bin_has_local;
  //   Buffer<int> bin_list;

  int mbins; // binning parameters
  int mbinx, mbiny, mbinz;
  int atoms_per_bin;
  int shared_mem_size;

private:
  md_float xprd, yprd, zprd; // box size

  int nmax;   // max size of atom arrays in neighbor
  int ntypes; // number of atom types

  int nstencil;        // # of bins in stencil
  Buffer<int> stencil; // stencil list of bin offsets

  int mbinxlo, mbinylo, mbinzlo;
  int nextx, nexty, nextz;
  md_float binsizex, binsizey, binsizez;
  md_float bininvx, bininvy, bininvz;

  int resize;

  //   md_float bindist(int, int, int); // distance between binx

  //   int coord2bin(md_float, md_float,
  //                 md_float) const; // mapping atom coord to a bin

  //   x_rnd_view_type x;
  //   int_1d_rnd_view_type type;
  int nlocal;
  Buffer<int> new_maxneighs;
  std::vector<int> h_new_maxneighs;
};

#endif