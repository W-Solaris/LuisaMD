#include <luisa/dsl/sugar.h>
#include <luisa/luisa-compute.h>

#include "lj.h"
#include "luisa/core/logging.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

using namespace luisa;
using namespace luisa::compute;

int input(In &, const char *);

int main(int argc, char **argv) {
  In in;
  in.datafile = NULL;
  int num_steps = -1;   // number of timesteps (if -1 use value from lj.in)
  int system_size = -1; // size of the system (if -1 use value from lj.in)
  int nx = -1;
  int ny = -1;
  int nz = -1;

  int screen_yaml = 0; // print yaml output to screen also
  int yaml_output = 0; // print yaml output
  int halfneigh = 1; // 1: use half neighborlist; 0: use full neighborlist; -1:
                     // use original miniMD version half neighborlist force
  int teams = 1;
  char *input_file = NULL;
  int ghost_newton = 1;
  int sort = -1;
  int ntypes = 8;

  for (int i = 0; i < argc; i++) {
    if ((strcmp(argv[i], "-i") == 0) ||
        (strcmp(argv[i], "--input_file") == 0)) {
      input_file = argv[++i];
      continue;
    }
  }

  int error = 0;

  if (input_file == NULL)
    error = input(in, "in.lj.miniMD");
  else
    error = input(in, input_file);

  if (error) {
    exit(0);
  }

  srand(5413);

  for (int i = 0; i < argc; i++) {
    if ((strcmp(argv[i], "-n") == 0) || (strcmp(argv[i], "--nsteps") == 0)) {
      num_steps = atoi(argv[++i]);
      continue;
    }

    if ((strcmp(argv[i], "-s") == 0) || (strcmp(argv[i], "--size") == 0)) {
      system_size = atoi(argv[++i]);
      continue;
    }

    if ((strcmp(argv[i], "-nx") == 0)) {
      nx = atoi(argv[++i]);
      continue;
    }

    if ((strcmp(argv[i], "-ny") == 0)) {
      ny = atoi(argv[++i]);
      continue;
    }

    if ((strcmp(argv[i], "-nz") == 0)) {
      nz = atoi(argv[++i]);
      continue;
    }

    if ((strcmp(argv[i], "--ntypes") == 0)) {
      ntypes = atoi(argv[++i]);
      continue;
    }

    if ((strcmp(argv[i], "--sort") == 0)) {
      sort = atoi(argv[++i]);
      continue;
    }

    if ((strcmp(argv[i], "-o") == 0) ||
        (strcmp(argv[i], "--yaml_output") == 0)) {
      yaml_output = atoi(argv[++i]);
      continue;
    }

    if ((strcmp(argv[i], "--yaml_screen") == 0)) {
      screen_yaml = 1;
      continue;
    }

    if ((strcmp(argv[i], "-f") == 0) || (strcmp(argv[i], "--data_file") == 0)) {
      if (in.datafile == NULL)
        in.datafile = new char[1000];

      strcpy(in.datafile, argv[++i]);
      continue;
    }

    if ((strcmp(argv[i], "-u") == 0) || (strcmp(argv[i], "--units") == 0)) {
      in.units = strcmp(argv[++i], "metal") == 0 ? 1 : 0;
      continue;
    }

    if ((strcmp(argv[i], "-p") == 0) || (strcmp(argv[i], "--force") == 0)) {
      in.forcetype = strcmp(argv[++i], "eam") == 0 ? FORCEEAM : FORCELJ;
      continue;
    }

    if ((strcmp(argv[i], "-gn") == 0) ||
        (strcmp(argv[i], "--ghost_newton") == 0)) {
      ghost_newton = atoi(argv[++i]);
      continue;
    }

    if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "--help") == 0)) {
      printf("Commandline Options:\n");
      printf("\n  Execution configuration:\n");
      printf(
          "\t-gn / --ghost_newton <int>:   set usage of newtons third law for "
          "ghost atoms\n"
          "\t                                (only applicable with half "
          "neighborlists)\n");
      printf("\n  Simulation setup:\n");
      printf(
          "\t-i / --input_file <string>:   set input file to be used (default: "
          "in.lj.miniMD)\n");
      printf("\t-n / --nsteps <int>:          set number of timesteps for "
             "simulation\n");
      printf("\t-s / --size <int>:            set linear dimension of "
             "systembox\n");
      printf(
          "\t-nx/-ny/-nz <int>:            set linear dimension of systembox "
          "in x/y/z direction\n");
      printf("\t-b / --neigh_bins <int>:      set linear dimension of neighbor "
             "bin grid\n");
      printf(
          "\t-u / --units <string>:        set units (lj or metal), see LAMMPS "
          "documentation\n");
      printf("\t-p / --force <string>:        set interaction model (lj or "
             "eam)\n");
      printf(
          "\t-f / --data_file <string>:    read configuration from LAMMPS data "
          "file\n");
      printf("\t--ntypes <int>:               set number of atom types for "
             "simulation (default 8)\n");

      printf("\t--sort <n>:                   resort atoms (simple bins) every "
             "<n> steps (default: use reneigh frequency; never=0)");
      printf(
          "\t-o / --yaml_output <int>:     level of yaml output (default 1)\n");
      printf(
          "\t--yaml_screen:                write yaml output also to screen\n");
      printf("\t-h / --help:                  display this help message\n\n");
      printf("---------------------------------------------------------\n\n");

      exit(0);
    }

    Atom atom(ntypes);
    Neighbor neighbor(ntypes);
    Integrate integrate;
    Thermo thermo;
    Timer timer;

    Force *force;

    if (in.forcetype == FORCELJ)
      force = (Force *)new ForceLJ(ntypes);

    if (in.forcetype == FORCELJ) {
      Buffer<float> d_epsilon = device.create_buffer<float>(ntypes * ntypes);
      std::vector<float> h_epsilon(ntypes * ntypes);
      force->epsilon = d_epsilon;
      force->epsilon_scalar = in.epsilon;

      Buffer<float> d_sigma6 = device.create_buffer<float>(ntypes * ntypes);
      std::vector<float> h_sigma6(ntypes * ntypes);
      force->sigma6 = d_sigma6;

      Buffer<float> d_sigma = device.create_buffer<float>(ntypes * ntypes);
      std::vector<float> h_sigma(ntypes * ntypes);
      force->sigma = d_sigma;
      force->sigma_scalar = in.sigma;

      for (int i = 0; i < ntypes * ntypes; i++) {
        h_epsilon[i] = in.epsilon;
        h_sigma[i] = in.sigma;
        h_sigma6[i] =
            in.sigma * in.sigma * in.sigma * in.sigma * in.sigma * in.sigma;
        if (i < MAX_STACK_TYPES * MAX_STACK_TYPES) {
          force->epsilon_s[i] = h_epsilon[i];
          force->sigma6_s[i] = h_sigma6[i];
        }
      }
      Kokkos::deep_copy(d_epsilon, h_epsilon);
      Kokkos::deep_copy(d_sigma6, h_sigma6);
      Kokkos::deep_copy(d_sigma, h_sigma);
    }
  }

  Context context{argv[0]};
  Device device = context.create_device("cuda");
  Stream stream = device.create_stream(StreamTag::COMPUTE);
}