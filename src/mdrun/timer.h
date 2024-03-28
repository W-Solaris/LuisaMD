#ifndef TIMER_H
#define TIMER_H
#define PREC_TIMER
#define TIME_TOTAL 0
#define TIME_COMM 1
#define TIME_FORCE 2
#define TIME_NEIGH 3
#define TIME_TEST 4
#define TIME_N 5

#include <ctime>

class Timer {
public:
  Timer();
  ~Timer();
  void stamp();
  void stamp(int);
  void stamp_extra_start();
  void stamp_extra_stop(int);
  void barrier_start(int);
  void barrier_stop(int);
  double *array;

private:
#ifdef PREC_TIMER
  timespec previous_time, previous_time_extra;
#endif
  double previous_time_d, previous_time_extra_d;
};

#endif
