#include "thermo.h"
#include "force_lj.h"
#include "integrate.h"
#include "luisa/dsl/stmt.h"
#include "stdio.h"
#include "stdlib.h"

Thermo::Thermo() {}
Thermo::~Thermo() {}

void Thermo::setup(Stream &stream, Device &device, float rho_in,
                   Integrate &integrate, Atom &atom, int units) {
  rho = rho_in;
  ntimes = integrate.ntimes;
  nlocal = atom.nlocal;

  if (nstat == 0)
    maxstat = 2;
  else
    maxstat = ntimes / nstat + 2;

  steparr = device.create_buffer<int>(maxstat);
  tmparr = device.create_buffer<float>(maxstat);
  engarr = device.create_buffer<float>(maxstat);
  prsarr = device.create_buffer<float>(maxstat);

  mstat = device.create_buffer<int>(1);
  t_act = device.create_buffer<float>(1);

  if (units == LJ) {
    mvv2e = 1.0;
    dof_boltz = atom.natoms * 3 - 3;
    t_scale = mvv2e / dof_boltz;
    p_scale = 1.0 / 3.0 / atom.box.xlen / atom.box.ylen / atom.box.zlen;
    e_scale = 0.5;

  } else if (units == METAL) {
    mvv2e = 1.036427e-04;
    dof_boltz = (atom.natoms * 3 - 3) * 8.617343e-05;
    t_scale = mvv2e / dof_boltz;
    p_scale = 1.602176e+06 / 3 / atom.box.xlen / atom.box.ylen / atom.box.zlen;
    e_scale = 524287.985533; // 16.0;
    integrate.dtforce /= mvv2e;
  }
}

void Thermo::setup_shader(Device &device, Atom &atom, Neighbor &neighbor,
                          Force *force) {
  Kernel1D reset_kernel = [&]() noexcept { t_act->write(0, 0.f); };

  Kernel1D temperature_kernel = [&]() noexcept {
    Float t = 0.f;
    Int max = 0;
    Float max_t = 0;
    Int target = 2000;
    $for(i, atom.nlocal) {
      Float3 v_ = atom.v->read(i);
      $if(i == target) {
        device_log("velocity: ({}, {}, {})", v_[0], v_[1], v_[2]);
      };
      t += (v_[0] * v_[0] + v_[1] * v_[1] + v_[2] * v_[2]) * atom.mass;
      // Float dt = (v_[0] * v_[0] + v_[1] * v_[1] + v_[2] * v_[2]) * atom.mass;
      // $if(dt > max_t) {
      //   max_t = dt;
      //   max = i;
      // };
    };
    // device_log("max_i: {}", max);
    // device_log("max_t: {}", max_t);

    t_act->write(0, t);
    // device_log("t: {}", t);
  };

  Kernel1D e_p_kernel = [&](Int iflag) noexcept {
    Int natoms_ = atom.natoms;

    Int istep = iflag;
    Int mstat_;
    $if(iflag == 0) { mstat_ = 0; }
    $else { mstat_ = mstat->read(0); };

    $if(iflag == -1) { istep = ntimes; };

    Float t_ = t_act->read(0) * t_scale;
    Float e_ = force->eng_vdwl->read(0) * e_scale / natoms_.cast<float>();
    Float p_ = (t_ * dof_boltz + force->virial->read(0)) * p_scale;

    // device_log("natoms: {}", natoms_);
    // device_log("eng_vdwl: {}", force->eng_vdwl->read(0));
    tmparr->write(mstat_, t_);
    engarr->write(mstat_, e_);
    prsarr->write(mstat_, p_);
    steparr->write(mstat_, istep);

    mstat_ += 1;
    mstat->write(0, mstat_);
    // device_log("t: {}", t_);
  };

  auto o = luisa::compute::ShaderOption{.enable_fast_math = false};
  reset_shader = device.compile(reset_kernel, o);
  temperature_shader = device.compile(temperature_kernel, o);
  e_p_shader = device.compile(e_p_kernel, o);
}

void Thermo::compute(Stream &stream, int iflag) {
  float t, eng, p;

  if (iflag > 0 && iflag % nstat)
    return;

  if (iflag == -1 && nstat > 0 && ntimes % nstat == 0)
    return;

  temperature(stream);

  stream << e_p_shader(iflag).dispatch(1) << synchronize();

  // TODO: output
}

void Thermo::temperature(Stream &stream) {
  stream << reset_shader().dispatch(1) << synchronize();
  stream << temperature_shader().dispatch(1) << synchronize();
}

void Thermo::output(Stream &stream, Device &device) {
  // copy info from device
  std::vector<int> h_mstat(1);
  // return;
  stream << mstat.copy_to(h_mstat.data());
  std::vector<int> h_steparr(maxstat);
  stream << steparr.copy_to(h_steparr.data());
  std::vector<float> h_tmparr(maxstat);
  stream << tmparr.copy_to(h_tmparr.data());
  std::vector<float> h_engarr(maxstat);
  stream << engarr.copy_to(h_engarr.data());
  std::vector<float> h_prsarr(maxstat);
  stream << prsarr.copy_to(h_prsarr.data());
  // printf("mstat: %i\n", h_mstat[0]);
  for (int i = 0; i < h_mstat[0]; i++)
    fprintf(stdout, "%d %.12f %.12f %.12f\n", h_steparr[i], h_tmparr[i],
            h_engarr[i], h_prsarr[i]);
}