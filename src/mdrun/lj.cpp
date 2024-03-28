#include <luisa/dsl/sugar.h>
#include <luisa/luisa-compute.h>

#include "Integrate.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "force_lj.h"
#include "lj.h"
#include "luisa/core/logging.h"
#include "neighbor.h"
#include "setup.cpp"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "thermo.h"
#include "variant.h"

using namespace luisa;
using namespace luisa::compute;

int input(In &, const char *);
void output(Stream &stream, Device &device, In &, Atom &, Force *, Neighbor &,
            Thermo &, Integrate &, int);
int main(int argc, char **argv) {
  luisa::log_level_verbose();
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

  int neighbor_size = -1;
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
    if ((strcmp(argv[i], "-b") == 0) ||
        (strcmp(argv[i], "--neigh_bins") == 0)) {
      neighbor_size = atoi(argv[++i]);
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
  }

  Context context{argv[0]};
  Device device = context.create_device("cuda");
  Stream stream = device.create_stream(StreamTag::COMPUTE);

  Atom atom(device, ntypes);
  Neighbor neighbor(device, ntypes);
  Integrate integrate;
  Thermo thermo;
  Comm comm;

  Force *force;

  if (in.forcetype == FORCELJ)
    force = (Force *)new ForceLJ(stream, device, in, ntypes);
  neighbor.ghost_newton = ghost_newton;
  if (num_steps > 0)
    in.ntimes = num_steps;
  if (system_size > 0) {
    in.nx = system_size;
    in.ny = system_size;
    in.nz = system_size;
  }
  if (nx > 0) {
    in.nx = nx;
    if (ny > 0)
      in.ny = ny;
    else if (system_size < 0)
      in.ny = nx;

    if (nz > 0)
      in.nz = nz;
    else if (system_size < 0)
      in.nz = nx;
  }
  if (neighbor_size > 0) {
    neighbor.nbinx = neighbor_size;
    neighbor.nbiny = neighbor_size;
    neighbor.nbinz = neighbor_size;
  }
  if (neighbor_size < 0 && in.datafile == NULL) {
    float neighscale = 0.6;
    neighbor.nbinx = neighscale * in.nx;
    neighbor.nbiny = neighscale * in.ny;
    neighbor.nbinz = neighscale * in.nz;
  }

  if (neighbor_size < 0 && in.datafile)
    neighbor.nbinx = -1;

  if (neighbor.nbinx == 0)
    neighbor.nbinx = 1;

  if (neighbor.nbiny == 0)
    neighbor.nbiny = 1;

  if (neighbor.nbinz == 0)
    neighbor.nbinz = 1;

  integrate.ntimes = in.ntimes;
  integrate.dt = in.dt;
  integrate.sort_every = sort > 0 ? sort : (sort < 0 ? in.neigh_every : 0);
  neighbor.every = in.neigh_every;
  neighbor.cutneigh = in.neigh_cut;
  force->cutforce = in.force_cut;
  thermo.nstat = in.thermo_nstat;
  printf("# Create System:\n");
  if (in.datafile) {
    read_lammps_data(stream, device, atom, neighbor, integrate, thermo,
                     in.datafile, in.units);
    float volume = atom.box.xlen * atom.box.ylen * atom.box.zlen;
    in.rho = 1.0 * atom.natoms / volume;
    force->setup(stream, device, atom);

    if (in.forcetype == FORCEEAM)
      atom.mass = force->mass;
  } else {
    create_box(atom, in.nx, in.ny, in.nz, in.rho);

    integrate.setup();

    if (in.forcetype == FORCEEAM)
      atom.mass = force->mass;

    create_atoms(atom, in.nx, in.ny, in.nz, in.rho);
    thermo.setup(stream, device, in.rho, integrate, atom, in.units);

    create_velocity(stream, device, in.t_request, atom, thermo, neighbor,
                    force);
    comm.setup(stream, device, neighbor.cutneigh, atom);
    comm.setup_shaders(device, atom);

    neighbor.setup(stream, device, atom);
    neighbor.setup_shader(device, atom);

    force->setup(stream, device, atom);
  }
  printf("# Done .... \n");
  fprintf(stdout, "# " VARIANT_STRING " output ...\n");
  fprintf(stdout, "# Run Settings: \n");
  fprintf(stdout, "\t# Inputfile: %s\n",
          input_file == 0 ? "in.lj.miniMD" : input_file);
  fprintf(stdout, "\t# Datafile: %s\n", in.datafile ? in.datafile : "None");
  fprintf(stdout, "# Physics Settings: \n");
  fprintf(stdout, "\t# ForceStyle: %s\n",
          in.forcetype == FORCELJ ? "LJ" : "EAM");
  fprintf(stdout, "\t# Force Parameters: %2.2lf %2.2lf\n", in.epsilon,
          in.sigma);
  fprintf(stdout, "\t# Units: %s\n", in.units == 0 ? "LJ" : "METAL");
  fprintf(stdout, "\t# Atoms: %i\n", atom.natoms);
  fprintf(stdout, "\t# Atom types: %i\n", atom.ntypes);
  fprintf(stdout,
          "\t# System size: %2.2lf %2.2lf %2.2lf (unit cells: %i %i %i)\n",
          atom.box.xlen, atom.box.ylen, atom.box.zlen, in.nx, in.ny, in.nz);
  fprintf(stdout, "\t# Density: %lf\n", in.rho);
  fprintf(stdout, "\t# Force cutoff: %lf\n", force->cutforce);
  fprintf(stdout, "\t# Timestep size: %lf\n", integrate.dt);
  fprintf(stdout, "# Technical Settings: \n");
  fprintf(stdout, "\t# Team neighborlist construction: %i\n",
          neighbor.team_neigh_build);
  // fprintf(stdout, "\t# Neighbor bins: %i %i %i\n", neighbor.nbinx,
  //         neighbor.nbiny, neighbor.nbinz);
  // fprintf(stdout, "\t# Neighbor frequency: %i\n", neighbor.every);
  // fprintf(stdout, "\t# Sorting frequency: %i\n", integrate.sort_every);
  fprintf(stdout, "\t# Thermo frequency: %i\n", thermo.nstat);
  // fprintf(stdout, "\t# Use intrinsics: %i\n", force->use_sse);
  fprintf(stdout, "\t# Size of float: %i\n\n", (int)sizeof(float));

  if (sort > 0)
    atom.sort(neighbor);

  comm.borders(stream);

  force->evflag = 1;
  neighbor.build(stream);

  force->setup_shader(device, atom, neighbor);
  atom.setup_shader(device);
  integrate.setup_shader(device, atom, force);

  force->compute(stream);
  // std::vector<float> tt(1);
  // stream << force->eng_vdwl.copy_to(tt.data()) << synchronize();
  // LUISA_INFO("energy: {}", tt[0]);

  // delete force;
  // return 0;

  printf("# Starting dynamics ...\n");

  printf("# Timestep T U P Time\n");

  { thermo.compute(stream, 0); }
  // thermo.output(stream, device);
  // delete force;
  // return 0;

  // run the main logic

  integrate.run(stream, atom, force, comm, neighbor, thermo);

  int natoms = atom.nlocal;

  force->evflag = 1;

  force->compute(stream);

  thermo.compute(stream, -1);

  // output thermo!
  thermo.output(stream, device);

  printf("\n\n");
  printf("# Performance Summary:\n");
  printf("# MPI_proc OMP_threads nsteps natoms t_total t_force t_neigh t_comm "
         "t_other performance perf/thread grep_string t_extra\n");
  // printf("%i %i %i %i %lf %lf %lf %lf %lf %lf %lf PERF_SUMMARY %lf\n\n\n",
  //        nprocs, num_threads, integrate.ntimes, natoms,
  //        timer.array[TIME_TOTAL], timer.array[TIME_FORCE],
  //        timer.array[TIME_NEIGH], timer.array[TIME_COMM], time_other, 1.0 *
  //        natoms * integrate.ntimes / timer.array[TIME_TOTAL], 1.0 * natoms *
  //        integrate.ntimes / timer.array[TIME_TOTAL] / nprocs /
  //            num_threads,
  //        timer.array[TIME_TEST]);

  // if (yaml_output)
  //   output(stream, device, in, atom, force, neighbor, thermo, integrate,
  //          screen_yaml);

  delete force;

  return 0;
}