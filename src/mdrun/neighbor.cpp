#include "luisa/dsl/stmt.h"
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

  cutneighsq = device.create_buffer<float>(ntypes_ * ntypes_);
  new_maxneighs = device.create_buffer<int>(1);
  team_neigh_build = 0;

  shared_mem_size = 0;
}

Neighbor::~Neighbor() {}

int Neighbor::setup(Stream &stream, Device &device, Atom &atom) {
  int i, j, k, nmax;
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

  nlocal = atom.nlocal;
  nmax = nlocal;
  /* extend atom arrays if necessary */

  numneigh = device.create_buffer<int>(nmax);
  maxneighs =
      nlocal - 1;  // only used for no bins: all other atoms as neighbors
  neighbors =
      device.create_image<int>(PixelStorage::INT1, make_uint2(nmax, maxneighs));

  count = 0;
  ntypes = atom.ntypes;
  // // TODO: bins implementation
  // // buffers initialized
  // mbins = 1;
  // atoms_per_bin = 1;
  // bincount = device.create_buffer<int>(mbins);
  // bin_has_local = device.create_buffer<int>(mbins);
  // bin_list = device.create_buffer<int>(mbins);
  // bins = device.create_image<int>(PixelStorage::INT1, mbins, atoms_per_bin);

  return 0;
}

void Neighbor::setup_shader(Device &device, Atom &atom) {
  Kernel1D neighbor_init_kernel = [&]() noexcept {
    auto i = dispatch_x();
    Int nall = atom.nlocal + atom.nghost->read(0);
    Float3 tmpi = atom.x->read(i);
    Float xtmp = tmpi[0];
    Float ytmp = tmpi[1];
    Float ztmp = tmpi[2];
    Int type_i = atom.type->read(i);

    UInt n = 0;

    $for(j, nall) {
      $if(j == i) { $continue; };

      Float3 tmpj = atom.x->read(j);
      Float delx = xtmp - tmpj[0];
      Float dely = ytmp - tmpj[1];
      Float delz = ztmp - tmpj[2];
      Int type_j = atom.type->read(j);
      Float dist = delx * delx + dely * dely + delz * delz;

      $if(dist <= cutneighsq->read(type_i * ntypes + type_j)) {
        $if(n < maxneighs) {
          neighbors->write(make_uint2(i, n), make_int4(j));
          n += 1;
        };
      };
    };
    numneigh->write(i, n.cast<int>());
    // $if(i == 2000) { device_log("num neighbors: {}", n); };

    // $if(n >= maxneighs) {
    //   $if(n >= new_maxneighs->read(0)) {
    //     new_maxneighs->write(0, n.cast<int>());
    //   };
    // };
  };

  neighbor_init = device.compile(neighbor_init_kernel);
}

void Neighbor::build(Stream &stream) {
  ncalls++;
  stream << neighbor_init().dispatch(nlocal) << synchronize();
}