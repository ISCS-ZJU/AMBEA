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

/* Enumerate biclique in bipartite graphs
 * Author: Yun Zhang
 * Created: September 2006
 * Last Update: August 2007
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utility.h"
#include "bigraph.h"
#include "biclique.h"

#define LMIN 2

/* Global Variables */
int LLEAST, RLEAST;
int VERSION;
int PRINT;

typedef unsigned long num_t;

#ifdef PERFORMANCE
long long node_num;
double time_check, time_expand, time_out, time_sort;
#endif



/* ------------------------------------------------------------- *
 * Function: biclique_profile_out()                              *
 *   Print out the profile, no. of left/right vertices, no. of   *
 *   edges, no.of biclique, size of vertex/edge maximum biclique *
 * ------------------------------------------------------------- */
void biclique_profile_out(FILE *fp, BiGraph *G, num_t *nclique)
{
  unsigned int n1 = G->_num_v1;
  unsigned int n2 = G->_num_v2;
  num_t num, sum=0;
  int ei=0, ej=0;
  int vi=0, vj=0;
  int i, j;
 
  fprintf(fp, "|Right Vertex|\t|Left Vertex|\tNumber\n");
  for (i = 1; i <= n2; i++) {
	for (j = 1; j <= n1; j++) {
	  num = nclique[(i-1)*n1+(j-1)];
	  if (num > 0) {
		fprintf(fp, "%d\t%d\t%lu\n", i, j, num);
		sum += num;
		if (i*j > ei*ej) { ei = i; ej = j; }
		if (i+j > vi+vj) { vi = i; vj = j; }
	  }
	}
  }
  
  fprintf(fp, "\n");
  fprintf(fp, "Number of left vertices     : %d\n", n1);
  fprintf(fp, "Number of right vertices    : %d\n", n2);
  fprintf(fp, "Number of edges             : %d\n", G->_num_edges);
  fprintf(fp, "Number of bicliques         : %lu\n", sum);
  fprintf(fp, "Size of edge max biclique   : (%d, %d)\n", ei, ej);
  fprintf(fp, "Size of vertex max biclique : (%d, %d)\n", vi, vj);
  fprintf(fp, "\n");
}



/* ------------------------------------------------------------- *
 * Function: biclique_out()                                      *
 *   Print out a biclique: left and right vertices on two lines  *
 * ------------------------------------------------------------- */
void biclique_out(FILE *fp, BiGraph *G, vid_t *right, \
		int nr, vid_t *left, int nl)
{
  int i;
#ifdef PERFORMANCE
  double utime = get_cur_time();
#endif
  for (i = 0; i < nr-1; i++) {
	fprintf(fp, "%s\t", G->_label_v2[right[i]]);
  }
  fprintf(fp, "%s\n", G->_label_v2[right[i]]);
  for (i = 0; i < nl-1; i++) {
	fprintf(fp, "%s\t", G->_label_v1[left[i]]);
  }
  fprintf(fp, "%s\n", G->_label_v1[left[i]]);
  fprintf(fp, "\n");
#ifdef PERFORMANCE
  time_out += get_cur_time() - utime;
#endif
}



/* ------------------------------------------------------------- *
 * Function: searchtreenode_out()                                *
 *   Print out the elements in biclique, candidates, and former  *
 *   candidates, in both left and right.                         *
 * ------------------------------------------------------------- */
void searchtreenode_out(BiGraph *G, vid_t *clique, vid_t *right, \
		vid_t *left, int nc, int ne, int ce, int nl) 
{
  int i;
  for (i = 0; i < nc; i++) printf(" %s", G->_label_v2[clique[i]]);
  printf("\t|");
  for (i = 0; i < ne; i++) printf(" %s", G->_label_v2[right[i]]);
  printf("\t|");
  for (i = ne; i < ce; i++) printf(" %s", G->_label_v2[right[i]]);
  printf("\t|");
  for (i = 0; i < nl; i++) printf(" %s", G->_label_v1[left[i]]);
  printf("\n");
}

void searchtreenode_out2(BiGraph *G, vid_t *clique, vid_t *right, \
		vid_t *left, vid_t *old_right, vid_t w, int nc, int ne, \
		int ce, int nl, int old_ne) 
{
  int i;
  for (i = 0; i < nc; i++) printf(" %s", G->_label_v2[clique[i]]);
  printf("\t|");
  for (i = 0; i < ne; i++) printf(" %s", G->_label_v2[right[i]]);
  printf(" %s", G->_label_v2[w]); 
  printf("\t|");
  for (i = old_ne+1; i < ce; i++) printf(" %s", G->_label_v2[old_right[i]]);
  printf("\t|");
  for (i = 0; i < nl; i++) printf(" %s", G->_label_v1[left[i]]);
  printf("\n");
}



/* --------------------------------------------------- *
 *  Function: biclique_find_basic()
 *    The basic version of biclique_find().
 * --------------------------------------------------- */
void biclique_find_basic(FILE *fp, BiGraph *G, num_t *nclique, \
	vid_t *clique, int nc, vid_t *left, int nl, vid_t *right, int ne, int ce)
{
  unsigned int n1 = G->_num_v1;
  vid_t new_left[nl];
  vid_t new_right[ce];
  vid_t u, v, w, j, k;
  int new_nc, new_nl, new_ne, new_ce;
  int count, is_maximal=1;
#ifdef PERFORMANCE
  double utime;
#endif

  while (ne < ce) {

	/* Choose one vertex from candidate set */
	v = right[ne];
	
#ifdef DEBUG
  searchtreenode_out(G, clique, right, left, nc, ne, ce, nl);
#endif

#ifdef PERFORMANCE
  node_num++;
  utime = get_cur_time();
#endif

    /* Set right vertices in clique */
    new_nc = nc;
	clique[new_nc++] = v;
	
    /* Set neighbors on left */
	new_nl = 0;
	for (j = 0; j < nl; j++) {
	  u = left[j];
	  if (bigraph_edge_exists(G, u, v)) new_left[new_nl++] = u;
	}
	
#ifdef PERFORMANCE
  time_expand += get_cur_time() - utime;
  utime = get_cur_time();
#endif

    /* Set right vertices in not */
	new_ne = 0;
	is_maximal = 1;
	for (j = 0; j < ne; j++) {
	  w = right[j];
	  count = 0;
	  for (k = 0; k < new_nl; k++) {
		u = new_left[k];
		if (bigraph_edge_exists(G, u, w)) count++;
	  }
	  if (count == new_nl) { is_maximal = 0; break; }
	  else if (count > 0) new_right[new_ne++] = w;
	}

    /* Stop this branch if it is not maximal */
	if (!is_maximal) { 
#ifdef DEBUG
  searchtreenode_out2(G, clique, new_right, new_left, right, w, \
		  new_nc, new_ne, new_ce, new_nl, ne);
#endif
      ne++; continue;
	}
	
#ifdef PERFORMANCE
  time_check += get_cur_time() - utime;
  utime = get_cur_time();
#endif

    /* Set right vertices in cand */
    new_ce = new_ne;
	for (j = ne+1; j < ce; j++) {
	  w = right[j];
	  count = 0;
	  for (k = 0; k < new_nl; k++) {
		u = new_left[k];
		if (bigraph_edge_exists(G, u, w)) count++;
	  }
	  if (count == new_nl) { 
		clique[new_nc++] = w; 
	  }
	  else if (count > 0)  { 
		new_right[new_ce++] = w;
	  }
	}
	
#ifdef PERFORMANCE
  time_expand += get_cur_time() - utime;
#endif

    /* Print out the found maximal biclique */
    if (new_nc >= RLEAST && new_nl >= LLEAST) {
	  nclique[(new_nc-1)*n1+(new_nl-1)]++;
	  if (PRINT) biclique_out(fp, G, clique, new_nc, new_left, new_nl);
    }
	
    /* Recursively find bicliques */
	if ((new_ne < new_ce) && (new_nc+(new_ce-new_ne) >= RLEAST)) {
	  biclique_find_basic(fp, G, nclique, clique, new_nc, \
		new_left, new_nl, new_right, new_ne, new_ce);
	}
#ifdef DEBUG
	else {
  searchtreenode_out(G, clique, new_right, new_left, \
		  new_nc, new_ne, new_ce, new_nl);
   }
#endif

	/* Move v to former candidate set */
	ne++;

  }
  
  return;
}


/* --------------------------------------------------- *
 *  Function: biclique_find_improve()
 *    The Improved Version of biclique_find().
 * --------------------------------------------------- */
void biclique_find_improve(FILE *fp, BiGraph *G, num_t *nclique, \
	vid_t *clique, int nc, vid_t *left, int nl, vid_t *right, \
	int ne, int ce)
{
  unsigned int n1 = G->_num_v1;
  vid_t new_left[nl];
  vid_t new_right[ce];
  vid_t u, v, w, j, k;
  int new_nc, new_nl, new_ne, new_ce;
  int count, is_maximal=1;
  int x, noc[ce-ne];
#ifdef PERFORMANCE
  double utime;
#endif
  // Improvement II
  // divide new_left to two parts: in L' and not in L'
  // ----------------------------
  // |  L'          |   not L'  |
  // ----------------------------
  // ^              ^           ^
  // new_nl ->          <- not_nl  
  // L'=new_left[0..new_nl-1], ~L'=new_left[new_nl,nl-1]
  int not_nl;  // position of not_left
  int nn;      // number of vertices directly go to not
  // End
 
  /* Same operations as v2 on each candidate in order */
  while (ne < ce) {
  
	/* Choose the next candidate in P */
	v = right[ne];
	
#ifdef DEBUG
  searchtreenode_out(G, clique, right, left, nc, ne, ce, nl);
#endif

#ifdef PERFORMANCE
  node_num++;
  utime = get_cur_time();
#endif

	/* Choose one vertex from candidate set */
	new_nc = nc;
	clique[new_nc++] = v;
	
    /* Set right vertices in clique */
	new_nl = 0; not_nl = nl;
	for (j = 0; j < nl; j++) {
	  u = left[j];
	  if (bigraph_edge_exists(G, u, v))
	    new_left[new_nl++] = u;
	  else new_left[not_nl--] = u;
	}
	
#ifdef PERFORMANCE
  time_expand += get_cur_time() - utime;
  utime = get_cur_time();
#endif

    /* Set right vertices in not */
	new_ne = 0;
	is_maximal = 1;
	for (j = 0; j < ne; j++) {
	  w = right[j];
	  count = 0;
	  for (k = 0; k < new_nl; k++) {
		u = new_left[k];
		if (bigraph_edge_exists(G, u, w)) count++;
	  }
	  if (count == new_nl) { is_maximal = 0; break; }
	  else if (count > 0) new_right[new_ne++] = w;
	}

    /* Stop this branch if it is not maximal */
	if (!is_maximal) { 
#ifdef DEBUG
  searchtreenode_out2(G, clique, new_right, new_left, right, w, \
		  new_nc, new_ne, new_ce, new_nl, ne);
#endif
      ne++; continue;
	}
	
#ifdef PERFORMANCE
  time_check += get_cur_time() - utime;
  utime = get_cur_time();
#endif

    /* Set right vertices in cand */
	new_ce = new_ne;
	nn = 1; // number of vertice will be put in not when backtracking
	for (j = ne+1; j < ce; j++) {
	  w = right[j];
      /* count the connections to L */
	  count = 0;
	  for (k = 0; k < new_nl; k++) {
		u = new_left[k];
		if (bigraph_edge_exists(G, u, w)) count++;
	  }
	  if (count == new_nl) {
		clique[new_nc++] = w;
		// Improvement II
		for (k = nl; k > not_nl; k--) {
			u = new_left[k];
			if (bigraph_edge_exists(G,u,w)) count++;
		}
		// switch such vertex with the one next to 
		// the last picked vertex to biclique
		if (count == new_nl) {
			right[j] = right[ne+nn]; 
			right[ne+nn] = w;
			nn++;
		}
	  }
	  else if (count > 0)  { 
	    /* Improvement I: Sorting candidates in new_right */
		x = new_ce-1;
		while (x >= new_ne && noc[x-new_ne] > count) {
		  noc[x+1-new_ne] = noc[x-new_ne];
		  new_right[x+1] = new_right[x];
		  x--;
	    }
	    noc[x+1-new_ne] = count;
	    new_right[x+1] = w;
		new_ce++;
	  }
	}

#ifdef PERFORMANCE
  time_expand += get_cur_time() - utime;
#endif

    /* Print out the found maximal biclique */
    if (new_nc >= RLEAST && new_nl >= LLEAST) {
	  nclique[(new_nc-1)*n1+(new_nl-1)]++;
	  if (PRINT) biclique_out(fp, G, clique, new_nc, new_left, new_nl);
    }
	
    /* Recursively find bicliques */
	if ((new_ne < new_ce) && (new_nc+(new_ce-new_ne) >= RLEAST)) {
	  biclique_find_improve(fp, G, nclique, clique, new_nc, \
		new_left, new_nl, new_right, new_ne, new_ce);
	}
	else {
#ifdef DEBUG
  searchtreenode_out(G, clique, new_right, new_left, \
		  new_nc, new_ne, new_ce, new_nl);
#endif
    }
	
	/* Move v and other qualified vertics to former candidate set */
	ne += nn;

  }
  
  return;
}



/* --------------------------------------------------- *
 *  Biclique Enumerating Main Function 
 * --------------------------------------------------- */
void biclique_enumerate(FILE *fp1, FILE *fp2, BiGraph *G, \
		vid_t *cand, int lcand)
{
  unsigned int n1 = G->_num_v1;
  unsigned int n2 = G->_num_v2;
  num_t *nclique;
  vid_t left[n1], right[n2], clique[n2];
  vid_t u, v;
  int noc[n2], tmpnoc, x;
  
  /* Initialization */
  nclique = (num_t *) calloc(n1*n2, sizeof(num_t));
  if (!nclique) { perror("malloc nclique\n"); exit(-1); }
  memset(clique, -1, n2*sizeof(vid_t));  // initially Clique is empty
  for (u = 0; u < n1; u++) left[u] = u;  // every left vertex is candidate 
  for (v = 0; v < n2; v++) right[v] = v; // every right vertex is candidate
 
  /* Call the recursive function to find maximal bicliques */
  if (VERSION == 1) {
    biclique_find_basic(fp1, G, nclique, clique, 0, left, n1, right, 0, n2);
  }
  else if (VERSION == 2) {
	/* Sort the candidate right vertices */
    memset(noc, 0, n2*sizeof(int));
	for (v = 0; v < n2; v++) {
	  tmpnoc = 0;
	  for (u = 0; u < n1; u++)
		if (bigraph_edge_exists(G, u, v)) tmpnoc++;
	  /* Sorting candidates in new_right */
	  x = v - 1;
	  while (x >= 0 && noc[x] > tmpnoc) {
		  noc[x+1] = noc[x];
		  right[x+1] = right[x];
		  x--;
	  }
	  noc[x+1] = tmpnoc;
	  right[x+1] = v;
    }	    
    biclique_find_improve(fp1, G, nclique, clique, 0, left, n1, right, 0, n2);
  }
 
  /* Print out the profile of maximal bicliques */
  biclique_profile_out(fp2, G, nclique);

  /* Free memory */
  free(nclique);

  return;
}

