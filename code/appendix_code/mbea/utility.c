// Copyright 2006-2014 Yun Zhang 
// Department of Computer Science
// University of Tennessee
//
// This file is part of MBEA.                  
//
// MBEA is free software: you can redistribute it and/or modify      
// it under the terms of the GNU General Public License as published by   
// the Free Software Foundation, either version 3 of the License, or   
// (at your option) any later version.   
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of   
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   
// GNU General Public License for more details.   
//
// You should have received a copy of the GNU General Public License   
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

/* APIs of utility for undirected graphs */
/* Implement functions defined in utility.h */
/* Graph Theory Team, Computer Science Department */
/* University of Tennessee, Knoxville */
/* Yun Zhang, yzhang@cs.utk.edu, January 20, 2005 */


#include <sys/time.h>  /* gettimeofday */
#include "utility.h"
 

/* Returns current time in seconds as a double value
   and the precision can be as accurate as microseconds
   (10^-6 second) 
 */
double get_cur_time() {
  struct timeval   tv;
  struct timezone  tz;
  double cur_time;
  gettimeofday(&tv, &tz);
  cur_time = tv.tv_sec + tv.tv_usec / 1000000.0;
  return cur_time;
}



