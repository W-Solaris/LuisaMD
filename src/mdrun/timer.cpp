#include "timer.h"
#include "stdlib.h"
#include "time.h"
#include "types.h"

Timer::Timer() {
  array = new double[TIME_N];

  for (int i = 0; i < TIME_N; i++)
    array[i] = 0.0;
}

Timer::~Timer() { delete[] array; }

void Timer::stamp() { clock_gettime(CLOCK_REALTIME, &previous_time); }

void Timer::stamp(int which) {
  timespec current_time;
  clock_gettime(CLOCK_REALTIME, &current_time);
  array[which] +=
      (current_time.tv_sec - previous_time.tv_sec +
       1.0 * (current_time.tv_nsec - previous_time.tv_nsec) / 1000000000);
  previous_time = current_time;
}

void Timer::stamp_extra_start() {
  clock_gettime(CLOCK_REALTIME, &previous_time_extra);
}

void Timer::stamp_extra_stop(int which) {
  timespec current_time;
  clock_gettime(CLOCK_REALTIME, &current_time);
  array[which] +=
      (current_time.tv_sec - previous_time_extra.tv_sec +
       1.0 * (current_time.tv_nsec - previous_time_extra.tv_nsec) / 1000000000);
  previous_time_extra = current_time;
}

void Timer::barrier_start(int which) {
  MPI_Barrier(MPI_COMM_WORLD);
  timespec current_time;
  clock_gettime(CLOCK_REALTIME, &current_time);
  array[which] = current_time.tv_sec + 1.0e-9 * current_time.tv_nsec;
}

void Timer::barrier_stop(int which) {
  MPI_Barrier(MPI_COMM_WORLD);
  timespec current_time;
  clock_gettime(CLOCK_REALTIME, &current_time);
  array[which] =
      current_time.tv_sec + 1.0e-9 * current_time.tv_nsec - array[which];
}
