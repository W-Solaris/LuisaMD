#include "comm.h"
#include "luisa/runtime/stream.h"
#include "stdio.h"
#include "stdlib.h"

#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 100
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

Comm::Comm() {}

Comm::~Comm() {}

/* setup spatial-decomposition communication patterns */

int Comm::setup(Stream &stream, Device &device, float cutneigh, Atom &atom) {
  int i;
  float lo, hi;
  int ineed, idim;

  need[0] = 1;
  need[1] = 1;
  need[2] = 1;

  ghost_num = device.create_buffer<int>(6);
  ghost_atoms =
      device.create_image<int>(PixelStorage::INT1, make_uint2(6, atom.nlocal));
  // to save 3 dims * 2 directions
  slablo = device.create_buffer<float>(6);
  slabhi = device.create_buffer<float>(6);
  pbc_flagx = device.create_buffer<float>(6);
  pbc_flagy = device.create_buffer<float>(6);
  pbc_flagz = device.create_buffer<float>(6);

  int maxswap = 2 * (need[0] + need[1] + need[2]);

  int iswap = 0;
  nswap = 0;

  std::vector<float> h_pbc_flagx = {1, -1, 0, 0, 0, 0};
  std::vector<float> h_pbc_flagy = {0, 0, 1, -1, 0, 0};
  std::vector<float> h_pbc_flagz = {0, 0, 0, 0, 1, -1};

  for (idim = 0; idim < 3; idim++) {
    for (ineed = 0; ineed < 2 * need[idim]; ineed++) {
      if (ineed % 2 == 0) {
        lo = 0;
        hi = cutneigh;
      } else {
        float total = 0;
        if (idim == 0)
          total = atom.box.xlen;
        if (idim == 1)
          total = atom.box.ylen;
        if (idim == 2)
          total = atom.box.zlen;
        lo = total - cutneigh;
        hi = total;
      }

      h_slablo.push_back(lo);
      h_slabhi.push_back(hi);
      nswap++;
    }
  }
  stream << slablo.copy_from(h_slablo.data())
         << slabhi.copy_from(h_slabhi.data());
  stream << pbc_flagx.copy_from(h_pbc_flagx.data())
         << pbc_flagy.copy_from(h_pbc_flagy.data())
         << pbc_flagz.copy_from(h_pbc_flagz.data()) << synchronize();

  return 0;
}

void Comm::setup_shaders(Device &device, Atom &atom) {
  Kernel1D communicate_kernel = [&]() noexcept {
    auto i = dispatch_x();
    Int offset = 0; // offset of atom.x to store ghost atoms
    $for(j, i) { offset += ghost_num->read(j); };

    $for(k, ghost_num->read(i)) {
      Int ghost_id = ghost_atoms->read(make_uint2(i, k.cast<uint>()))[0];
      Float3 ghostx = atom.x->read(ghost_id);
      ghostx.x += pbc_flagx->read(i) * atom.box.xlen;
      ghostx.y += pbc_flagy->read(i) * atom.box.ylen;
      ghostx.z += pbc_flagz->read(i) * atom.box.zlen;
      atom.x->write(atom.nlocal + offset + k, ghostx);
      atom.type->write(atom.nlocal + offset + k, atom.type->read(ghost_id));
      // $if(atom.nlocal + offset + k == 4690) {
      //   device_log("ghost_id: {}", ghost_id);
      // };
    };
  };

  Kernel1D borders_kernel = [&]() noexcept {
    auto i = dispatch_x();
    Float lo = slablo->read(i);
    Float hi = slabhi->read(i);
    Int nghost_i = 0;
    Int idim = (i - (i % 2)) / 2;
    $for(k, atom.nlocal) {
      // $if(i == 4 & k == 24) {
      //   device_log("low: {} ", lo);
      //   device_log("high: {} ", hi);
      //   device_log("idim: {} ", idim);
      //   device_log("z: {} ", atom.x->read(k)[2]);
      // };
      $if(lo <= atom.x->read(k)[idim] & atom.x->read(k)[idim] <= hi) {
        ghost_atoms->write(make_uint2(i, nghost_i.cast<uint>()), make_int4(k));
        nghost_i += 1;
        // $if(k == 24) { device_log("have!!!"); };
      };
    };
    ghost_num->write(i, nghost_i);
  };

  Kernel1D nghost_kernel = [&]() noexcept {
    Int num = 0;
    $for(i, 6) { num += ghost_num->read(i); };
    atom.nghost->write(0, num);
  };

  communicate_shader = device.compile(communicate_kernel);
  borders_shader = device.compile(borders_kernel);
  nghost_shader = device.compile(nghost_kernel);
}

void Comm::communicate(Stream &stream) {
  stream << communicate_shader().dispatch(6) << synchronize();
}

void Comm::borders(Stream &stream) {
  stream << borders_shader().dispatch(6) << synchronize();
  stream << nghost_shader().dispatch(1) << synchronize();
  stream << communicate_shader().dispatch(6) << synchronize();
}
