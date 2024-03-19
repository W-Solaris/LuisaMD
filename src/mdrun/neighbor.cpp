#include "neighbor.h"
#include "stdio.h"
#include "stdlib.h"

#define FACTOR 0.999
#define SMALL 1.0e-6

Neighbor::Neighbor(Device &device, int ntypes_) {
  ncalls = 0;
  ntypes = ntypes_;
  max_totalneigh = 0;
  maxneighs = 100;
  nmax = 0;
  atoms_per_bin = 8;
  halfneigh = 0;
  ghost_newton = 1;

  cutneighsq = device.create_buffer<float>(ntypes * ntypes);
  new_maxneighs = device.create_buffer<int>(1);
  team_neigh_build = 0;

  shared_mem_size = 0;
}

Neighbor::~Neighbor() {}

void Neighbor::build(Stream &stream, Device &device, Atom &atom) {
  ncalls++;
  nlocal = atom.nlocal;
  const int nall = atom.nlocal + atom.nghost;
  /* extend atom arrays if necessary */

  if (nall > nmax) {
    nmax = nall;

    numneigh = device.create_buffer<int>(nmax);
    maxneighs = nlocal;  // only used for no bins: all atoms as neighbors
    neighbors = device.create_image<int>(PixelStorage::BYTE4,
                                         make_uint2(nmax, maxneighs));
  }

  count = 0;
  ntypes = atom.ntypes;

  Kernel1D neighbor_init_kernel = [&]() noexcept {
    auto i = dispatch_x();
    Float xtmp = atom.x->read(i)[0];
    Float ytmp = atom.x->read(i)[1];
    Float ztmp = atom.x->read(i)[2];
    Int type_i = atom.type->read(i);

    UInt n = 0;

    $for(j, nlocal) {
      Float delx = xtmp - atom.x->read(j)[0];
      Float dely = xtmp - atom.x->read(j)[1];
      Float delz = xtmp - atom.x->read(j)[2];
      Int type_j = atom.type->read(j);
      Float dist = delx * delx + dely * dely + delz * delz;

      $if(dist <= cutneighsq->read(type_i * ntypes + type_j)) {
        UInt2 coord = $if(n < maxneighs) {
          neighbors->write(make_uint2(i, n), j);
        };
        n += 1;
      };
    };

    numneigh->write(i, n);
    $if(n >= maxneighs) {
      $if(n >= new_maxneighs->read(0)) { new_maxneighs->write(0, n); };
    };
  };

  auto neighbor_init = device.compile(neighbor_init_kernel);
  stream << neighbor_init().dispatch(nlocal) << synchronize();
}

int Neighbor::setup(Stream &stream, Device &device, Atom &atom) {
  int i, j, k, nmax;
  float coord;
  int mbinxhi, mbinyhi, mbinzhi;
  std::vector<float> h_cutneighsq(ntypes * ntypes);

  for (int i = 0; i < ntypes * ntypes; i++) {
    h_cutneighsq[i] = cutneigh * cutneigh;
    if (i < MAX_STACK_TYPES * MAX_STACK_TYPES)
      cutneighsq_stack[i] = cutneigh * cutneigh;
  }

  stream << cutneighsq.copy_from(h_cutneighsq.data());

  xprd = atom.box.xlen;
  yprd = atom.box.ylen;
  zprd = atom.box.zlen;

  // TODO: bins implementation
  // buffers initialized
  mbins = 1;
  atoms_per_bin = 1;
  bincount = device.create_buffer<int>(mbins);
  bin_has_local = device.create_buffer<int>(mbins);
  bin_list = device.create_buffer<int>(mbins);
  bins = device.create_image<int>(PixelStorage::BYTE4, mbins, atoms_per_bin);

  return 0;
}