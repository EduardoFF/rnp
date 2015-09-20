#include "CpuTime.h"
  
void CpuTime::start() 
{ 
  gettimeofday(&total_time_start, NULL); 
  times(&cpu_time_start); 
}
  
void CpuTime::end() 
{
  times(&cpu_time_end);
  
  gettimeofday(&total_time_end, NULL);
  
  cpu_time = (cpu_time_end.tms_utime - cpu_time_start.tms_utime) +
    (cpu_time_end.tms_stime - cpu_time_start.tms_stime);
  cpu_time /= (double)sysconf(_SC_CLK_TCK); 
  
  total_time = (total_time_end.tv_sec - total_time_start.tv_sec) +
    ((double)(total_time_end.tv_usec - total_time_start.tv_usec) / 1000000.0);
  
  /*
  printf("Cpu time elapsed: %f - Total time elapsed %f\n", 
	 cpu_time, total_time);
  fflush(stdout);
  */
}

void CpuTime::end(FILE *log_fp) 
{
  times(&cpu_time_end);
  
  gettimeofday(&total_time_end, NULL);
  
  cpu_time = (cpu_time_end.tms_utime - cpu_time_start.tms_utime) +
    (cpu_time_end.tms_stime - cpu_time_start.tms_stime);
  cpu_time /= (double)sysconf(_SC_CLK_TCK); 
  
  total_time = (total_time_end.tv_sec - total_time_start.tv_sec) +
    ((double)(total_time_end.tv_usec - total_time_start.tv_usec) / 1000000.0);
  
  /*
  printf("Cpu time elapsed: %f - Total time elapsed %f\n", 
	 cpu_time, total_time);
  fflush(stdout);
  */
  fprintf(log_fp, "Cpu time elapsed: %f - Total time elapsed %f\n", 
	  cpu_time, total_time);
  fflush(log_fp);
}

void CpuTime::end(FILE *log_fp, char *msg) 
{
  times(&cpu_time_end);
  
  gettimeofday(&total_time_end, NULL);
  
  cpu_time = (cpu_time_end.tms_utime - cpu_time_start.tms_utime) +
    (cpu_time_end.tms_stime - cpu_time_start.tms_stime);
  cpu_time /= (double)sysconf(_SC_CLK_TCK); 
  
  total_time = (total_time_end.tv_sec - total_time_start.tv_sec) +
    ((double)(total_time_end.tv_usec - total_time_start.tv_usec) / 1000000.0);
  
  fprintf(log_fp, "%s: %g\n", 
	  msg, cpu_time);
}

double CpuTime::cpu_time_elapsed(FILE *log_fp, char *msg) 
{
  times(&cpu_time_end);
  
  cpu_time = (cpu_time_end.tms_utime - cpu_time_start.tms_utime) +
    (cpu_time_end.tms_stime - cpu_time_start.tms_stime);
  cpu_time /= (double)sysconf(_SC_CLK_TCK); 

  fprintf(log_fp, "%s: Cpu time: %g\n",
	  msg, cpu_time);
  
  return cpu_time;
}

double CpuTime::cpu_time_elapsed()
{
  times(&cpu_time_end);
  
  cpu_time = (cpu_time_end.tms_utime - cpu_time_start.tms_utime) +
    (cpu_time_end.tms_stime - cpu_time_start.tms_stime);
  cpu_time /= (double)sysconf(_SC_CLK_TCK); 
  
  return cpu_time;

}

double CpuTime::total_time_elapsed() 
{
  gettimeofday(&total_time_end, NULL);
  
  total_time = (total_time_end.tv_sec - total_time_start.tv_sec) +
    ((double)(total_time_end.tv_usec - total_time_start.tv_usec) / 1000000.0);
  
  return total_time;
}

