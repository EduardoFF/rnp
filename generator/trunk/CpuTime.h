#ifndef CPU_TIME_H
#define CPU_TIME_H

#include <stdlib.h>
#include <stdio.h>
#include <sys/times.h>
#include <unistd.h>

// for gettimeofday()
#include <sys/time.h>


class CpuTime
{
  double 	  cpu_time;   
  double	  total_time;
  struct tms      cpu_time_start, cpu_time_end;
  struct timeval  total_time_start, total_time_end;


public:

  CpuTime() {}
  
  void start();

  void end();

  void end(FILE *log_fp); 

  void end(FILE *log_fp, char *msg);

  double cpu_time_elapsed();
  double cpu_time_elapsed(FILE *log_fp, char *msg);

  double total_time_elapsed();

};

#endif
