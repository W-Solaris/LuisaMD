#include "stdio.h"
#include "stdlib.h"

#include "lj.h"

int main(int argc, char **argv) {
  In in;
  in.datafile = NULL;
  int num_steps = -1;   // number of timesteps (if -1 use value from lj.in)
  int system_size = -1; // size of the system (if -1 use value from lj.in)
  int nx = -1;
  int ny = -1;
  int nz = -1;
  int check_safeexchange = 0; // if 1 complain if atom moves further than 1
                              // subdomain length between exchanges
  int do_safeexchange = 0; // if 1 use safe exchange mode [allows exchange over
                           // multiple subdomains]
  int use_sse = 0;         // setting for SSE variant of miniMD only
  int screen_yaml = 0;     // print yaml output to screen also
  int yaml_output = 0;     // print yaml output
  int halfneigh = 1; // 1: use half neighborlist; 0: use full neighborlist; -1:
                     // use original miniMD version half neighborlist force
  int teams = 1;
  int device = 0;
  int ngpu = 1;
  int skip_gpu = 99999;
  int neighbor_size = -1;
  char *input_file = NULL;
  int ghost_newton = 1;
  int sort = -1;
  int ntypes = 8;
  int team_neigh = 0;
}