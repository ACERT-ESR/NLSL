/*
 * Returns the delta time since the last call to dtime.
 *
 * calling sequence:
 * 	real time(2)
 * 	call dtime(time)
 * where:
 * 	the 2 element array time will receive the user and system
 * 	elapsed time since the last call to dtime, or since the start
 * 	of execution.
 *
 * This routine can be called as function, and returns the sum of
 * user and system times. The time_array argument must always be given.
 *
 * Original BSD-licensed routine replaced by WH 20020212 to fix
 * bugs and non-portable code. This routine should conform to POSIX.
 *
 */

#include <time.h>
#include <sys/times.h>
#include "fortrancall.h"

float FORTRAN(dtime)(float *dt)
{
	static time_t usr=0, sys=0;
	struct tms theclock;

	times(&theclock);
	dt[0] = (theclock.tms_utime - usr) / (float) CLK_TCK;
	dt[1] = (theclock.tms_stime - sys) / (float) CLK_TCK;
	usr = theclock.tms_utime;
	sys = theclock.tms_stime;
	return dt[0]+dt[1];
}
