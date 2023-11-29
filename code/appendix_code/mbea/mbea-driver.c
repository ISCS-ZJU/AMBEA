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

/* Driver file for MBEA 
 *
 * Author: Yun Zhang
 * Created: April 2006
 *
 * Copyright 2005-2014
 * Department of Computer Science, University of Tennessee, Knoxville
 * Computer Science and Mathematics Division, Oak Ridge National Laboratory
 *
 * Last modified: March 2014
 */

#include <string.h>
#include <stdio.h>
#include "bigraph.h"
#include "utility.h"
#include "biclique.h"


extern int LLEAST, RLEAST;
extern int CACHE_SIZE;
extern int VERSION;
extern int PRINT;
#ifdef PERFORMANCE
extern long long node_num;
extern double time_check, time_expand, time_out, time_sort;
double time_in, stime;
#endif
FILE *fp;
char infn[100];
char *outfn;

void print_options(void)
{
  fprintf(stderr, "\n Options: \n");
  fprintf(stderr, "  -v [1|2]       algorithm version <default = 2>\n");
  fprintf(stderr, "                   1 - original MBEA\n");
  fprintf(stderr, "                   2 - improved MBEA\n");
  fprintf(stderr, "  -p             print out bicliques <default = no print out>\n");
  fprintf(stderr, "  -o <filename>  output filename to store bicliques if choose to print\n");
  fprintf(stderr, "                 <default = STDOUT>\n");
  fprintf(stderr, "  -l <value>     least number of left vertices in biclique <default = 1>\n");
  fprintf(stderr, "  -r <value>     least number of right vertices in biclqiue <default = 1>\n");
  fprintf(stderr, "\n");
}

void argument_parse(int argc, char **argv)
{
  int i;
  
  if (argc < 2) {
    fprintf(stderr, "Usage: %s Graph <options>\n", argv[0]);
    print_options();
	exit(1);
  }
  
  LLEAST = 1;
  RLEAST = 1;
  VERSION = 2;
  PRINT = 0;
  outfn = NULL;
  
  for (i = 2; i < argc; i++) {
	if (!strcmp(argv[i], "-o")) {
	  outfn = strdup(argv[++i]);
	}
	if (!strcmp(argv[i], "-l")) {
	  LLEAST = atoi(argv[++i]);
	}
	if (!strcmp(argv[i], "-r")) {
	  RLEAST = atoi(argv[++i]);
	}
	if (!strcmp(argv[i], "-v")) {
	  VERSION = atoi(argv[++i]);
	}
	if (!strcmp(argv[i], "-p")) {
	  PRINT = 1;
	}
  }

  strcpy(infn, argv[1]);
  if ((fp = fopen(infn, "r")) == NULL) {
    fprintf(stderr, "Can't open file %s\n", infn);
    exit(-1);
  }

  return;
}


void maximal_biclique(char *fn, BiGraph *G)
{
  FILE *fp1=NULL, *fp2;
  char fname[100];
#ifdef PERFORMANCE
  double utime;
#endif
  int n2 = G->_num_v2;
  vid_t cand[n2];
  int i;

  if (PRINT) { 
	if (fn != NULL) { fp1 = fopen64(fn, "w"); }
    else { fp1 = stdout; }
  }
  if (fn != NULL) {
    sprintf(fname, "%s.profile", fn);
    fp2 = fopen(fname, "w");
  }
  else { fp2 = stdout; }

#ifdef PERFORMANCE
  node_num = 0;
  time_check = 0.0;
  time_expand = 0.0;
  time_out = 0.0;
  time_sort = 0.0;
#endif

  for (i = 0; i < n2; i++) cand[i] = i;
  biclique_enumerate(fp1, fp2, G, cand, n2);

#ifdef PERFORMANCE
  utime = get_cur_time() - stime;
  fprintf(fp2, "\n");
  fprintf(fp2, "%lld tree nodes\n", node_num);
  fprintf(fp2, "%f readin\n", time_in);
  fprintf(fp2, "%f expand\n", time_expand);
  fprintf(fp2, "%f check\n", time_check);
  fprintf(fp2, "%f output\n", time_out);
  fprintf(fp2, "%f sorting\n", time_sort);
  fprintf(fp2, "%f total time\n", utime);
#endif

#ifdef VERBOSE 
  printf("find biclique in %f seconds\n", utime);
#endif

  if (fp1 != NULL) fclose(fp1);
  fclose(fp2);

}


int main(int argc, char  **argv)
{
  BiGraph *G = NULL;

  argument_parse(argc, argv);
  
#ifdef PERFORMANCE
  stime = get_cur_time();
#endif

  G = bigraph_edgelist_in(fp);
  fclose(fp);

#ifdef PERFORMANCE
  time_in = get_cur_time() - stime;
#endif

  maximal_biclique(outfn, G);

  bigraph_free(G);

  exit(0);
}

