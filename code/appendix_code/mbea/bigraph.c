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

/* API for Bit-Based Adjacency Matrix for bipartite graphs */
/* Implement functions defined in bigraph.h */
/* Graph Theory Team, Computer Science Department */
/* University of Tennessee, Knoxville */
/* Yun Zhang, yzhang@eecs.utk.edu, December 12, 2004 */
/* Last Updated: 12/18/2014 */

#include <string.h>
#include <search.h>
#include "bigraph.h"
#include "utility.h"

#define LINE_LENGTH 1000

/* Malloc and initialize a graph, returns a pointer to it */
BiGraph *bigraph_make(unsigned int num_v1, unsigned int num_v2)
{
  BiGraph *G;
  int i;
  int num_ints_v1 = bit_num_ints(num_v2);
  int num_ints_v2 = bit_num_ints(num_v1);
  
  G = (BiGraph *) malloc(sizeof(BiGraph));
  if (G == NULL) { perror("malloc"); exit(-1); }
  G->_num_v1 = num_v1;
  G->_num_v2 = num_v2;
  G->_num_edges = 0;
  G->_num_bytes_v1 = num_ints_v1 * sizeof(int);
  G->_num_bytes_v2 = num_ints_v2 * sizeof(int);

  G->_label_v1 = (char **) malloc(num_v1 * sizeof(char *));
  G->_label_v2 = (char **) malloc(num_v2 * sizeof(char *));
  if (G->_label_v1 == NULL || G->_label_v2 == NULL) {
	perror("bigraph_make: malloc label");
	exit(-1);
  }

  G->_neighbor_v1 = (unsigned int **) malloc(num_v1 * sizeof(unsigned int *));
  if (G->_neighbor_v1 == NULL) { perror("malloc"); exit(-1); }
  G->_neighbor_v1[0] = (unsigned int *) malloc(G->_num_bytes_v1 * num_v1);
  if (G->_neighbor_v1[0] == NULL) { perror("malloc"); exit(-1); }
  for (i = 0; i < num_v1; i++)
    G->_neighbor_v1[i] = G->_neighbor_v1[0] + i * num_ints_v1;

  G->_neighbor_v2 = (unsigned int **) malloc(num_v2 * sizeof(unsigned int *));
  if (G->_neighbor_v2 == NULL) { perror("malloc"); exit(-1); }
  G->_neighbor_v2[0] = (unsigned int *) malloc(G->_num_bytes_v2 * num_v2);
  if (G->_neighbor_v2[0] == NULL) { perror("malloc"); exit(-1); }
  for (i = 0; i < num_v2; i++)
    G->_neighbor_v2[i] = G->_neighbor_v2[0] + i * num_ints_v2;
  
  G->_degree_v1 = (vid_t *) malloc(num_v1 * sizeof(vid_t));
  G->_degree_v2 = (vid_t *) malloc(num_v2 * sizeof(vid_t));
  if (!G->_degree_v1 | !G->_degree_v2) { perror("malloc degree"); exit(-1); }
  
  memset(G->_neighbor_v1[0], 0, G->_num_bytes_v1 * num_v1);
  memset(G->_neighbor_v2[0], 0, G->_num_bytes_v2 * num_v2);
  memset(G->_degree_v1, 0, num_v1 * sizeof(vid_t));
  memset(G->_degree_v2, 0, num_v2 * sizeof(vid_t));
 
  // added adjacency lists
  G->_adj_list_v1 = (vid_t **) malloc(num_v1 * sizeof(vid_t *));
  if (G->_adj_list_v1 == NULL) { perror("malloc"); exit(-1); }
  G->_adj_list_v2 = (vid_t **) malloc(num_v2 * sizeof(vid_t *));
  if (G->_adj_list_v2 == NULL) { perror("malloc"); exit(-1); }

//  fprintf(stderr,"bigraph_make() done\n");
  return G;
}


/* Create adjacency list after read in graph from a file */
void bigraph_adjlist_make(BiGraph *G)
{
//  fprintf(stderr,"bigraph_adjlist_make() start\n");
  unsigned int i, j, e;
  // malloc space
  G->_adj_list_v1[0] = (vid_t *) malloc(G->_num_edges * sizeof(vid_t));
  if (G->_adj_list_v1[0] == NULL) { perror("malloc"); exit(-1); }
  for (i = 1; i < G->_num_v1; i++) 
    G->_adj_list_v1[i] = G->_adj_list_v1[i-1] + G->_degree_v1[i-1];
  G->_adj_list_v2[0] = (vid_t *) malloc(G->_num_edges * sizeof(vid_t));
  if (G->_adj_list_v2[0] == NULL) { perror("malloc"); exit(-1); }
  for (i = 1; i < G->_num_v2; i++) 
    G->_adj_list_v2[i] = G->_adj_list_v2[i-1] + G->_degree_v2[i-1];
  // fill in neighbors' index
  for (i = 0; i < G->_num_v1; i++) {
    e = 0;
    for (j = 0; j < G->_num_v2; j++) {
      if (bigraph_edge_exists(G, i, j)) {
        G->_adj_list_v1[i][e++] = j;
      }
    }
  }
  for (j = 0; j < G->_num_v2; j++) {
    e = 0;
    for (i = 0; i < G->_num_v1; i++) {
      if (bigraph_edge_exists(G, i, j)) {
        G->_adj_list_v2[j][e++] = i;
      }
    }
  }
//  fprintf(stderr,"bigraph_adjlist_make() done\n");
}


/* Free the memory of a graph */
void bigraph_free(BiGraph *G)
{
  int i;
  if (G != NULL) {
    // adjacency matrix
    if (G->_neighbor_v1) {
      if (G->_neighbor_v1[0]) free(G->_neighbor_v1[0]);
      free(G->_neighbor_v1);
    }
    if (G->_neighbor_v2) {
      if (G->_neighbor_v2[0]) free(G->_neighbor_v2[0]);
      free(G->_neighbor_v2);
    }
	for (i = 0; i < G->_num_v1; i++) free(G->_label_v1[i]);
	for (i = 0; i < G->_num_v2; i++) free(G->_label_v2[i]);
    if (G->_degree_v1) free(G->_degree_v1);
    if (G->_degree_v2) free(G->_degree_v2);
    // added adjacency lists
    if (G->_adj_list_v1) {
      if (G->_adj_list_v1[0]) free(G->_adj_list_v1[0]);
      free(G->_adj_list_v1);
    }
    if (G->_adj_list_v2) {
      if (G->_adj_list_v2[0]) free(G->_adj_list_v2[0]);
      free(G->_adj_list_v2);
    }
	free(G);
  }
}

/** I/O functions for Graph **/

/* Read in an unweighted edge-list file, return a pointer to the graph */
BiGraph * bigraph_edgelist_in(FILE *fp)
{
  int u, v;
  unsigned int n1, n2, e;
  int k1=0, k2=0, edges=0, r;
  char word1[100], word2[100];
  BiGraph *G;
  int *id, *id1, *id2;
  ENTRY item;
  ENTRY *found_item;
  
  if (fscanf(fp, "%d %d %d", &n1, &n2, &e) != 3) {
    fprintf(stderr, "Bad file format: n1 n2 e incorrect\n");
	exit(1);
  }
  
  G = bigraph_make(n1, n2);
  
  /* create a hash table */
  (void) hcreate(n1+n2);
  id1 = (int *) malloc(n1 * sizeof(int));
  id2 = (int *) malloc(n2 * sizeof(int));
  
  while ((r = fscanf(fp, "%s\t%s", word1, word2)) != EOF) {
	if (r != 2) {
	  fprintf(stderr, "Bad file format: label1 label2 incorrect\n");
	  exit(1);
	}

    item.key = word1;
	if ((found_item = hsearch(item, FIND)) != NULL) {
		id = (int *) (found_item->data);
		u = *id;
	}
	else {
		u = k1; 
		G->_label_v1[k1++] = strdup(word1);
		item.key = G->_label_v1[u];
		id1[u] = u;
		item.data = (void *) (id1+u);
		(void) hsearch(item, ENTER);
	}

	item.key = word2;
	if ((found_item = hsearch(item, FIND)) != NULL) {
		id = (int *) (found_item->data);
		v = *id;
	}
	else {
		v = k2; 
		G->_label_v2[k2++] = strdup(word2);
		item.key = G->_label_v2[v];
		id2[v] = v;
		item.data = (void *) (id2+v);
		(void) hsearch(item, ENTER);
	}

	if (k1 > n1) {
	  fprintf(stderr, "Bad file format: too many left vertex labels\n");
	  exit(1);
	}
	if (k2 > n2) {
	  fprintf(stderr, "Bad file format: too many right vertex labels\n");
	  exit(1);
	}
	
    bigraph_add_edge(G, u, v);
	edges++;
  }
  
  if (edges != e) 
	fprintf(stderr, "edgelist_in: number of edges incorrect\n");
  if (k1 != n1) 
	fprintf(stderr, "edgelist_in: number of left vertices incorrect %d\n", k1);
  if (k2 != n2) 
	fprintf(stderr, "edgelist_in: number of right vertices incorrect %d\n", k2);

  (void) hdestroy();
  free(id1);
  free(id2);
 
  bigraph_adjlist_make(G);

  return G;
}


/* Read in a bigraph from a binary matrix format file */
BiGraph * bigraph_binarymatrix_in(FILE *fp)
{
  BiGraph *G;
  char line[LINE_LENGTH];
  char delims[] = " \t\n";
  char *a = NULL;
  int n1, n2, k1, k2, i, j;


  if (fgets(line, LINE_LENGTH, fp) == NULL) {
    fprintf(stderr,"ERROR: Input file is empty.\n");
	exit(1);
  }

  a = strtok(line, delims);
  n1 = atoi(a);
  a = strtok(NULL, delims);
  n2 = atoi(a);
  
  G = bigraph_make(n1, n2);
  
  if (fgets(line, LINE_LENGTH, fp) == NULL) {
    fprintf(stderr,"ERROR: Input file is invalid.\n");
    exit(1);
  }
  a = strtok(line, delims);
  k2 = 0;
  G->_label_v2[k2++] = strdup(a);
  while ((a = strtok(NULL, delims)) != NULL) {
    G->_label_v2[k2++] = strdup(a);
  }

  k1 = 0;
  while (fgets(line, LINE_LENGTH, fp) != NULL) {
    a = strtok(line, delims);
    G->_label_v1[k1] = strdup(a);
	i = 0;
    while ((a = strtok(NULL, delims)) != NULL) {
	  j = atoi(a);
      if (j==1) { bigraph_add_edge(G,k1,i); }
	  i++;
	}
	k1++;
  }

  if (k1 != n1) 
	fprintf(stderr, "binarymatrix_in: # left vertices incorret %d\n", k1);
  if (k2 != n2) 
	fprintf(stderr, "binarymatrix_in: # right vertices incorret %d\n", k2);

  return G;
}


/* Print out the edgelist of a bigraph to a file */
void bigraph_edgelist_out(FILE *fp, BiGraph *G)
{
  unsigned int n1 = bigraph_num_v1(G);
  unsigned int n2 = bigraph_num_v2(G);
  unsigned int e = bigraph_num_edges(G);
  unsigned int i, j;
  fprintf(fp, "EdgeList %d %d %d\n", n1, n2, e);
  for (i = 0; i < n1; i++) {
	for (j = 0; j < n2; j++) {
	  if (bigraph_edge_exists(G, i, j)) {
		fprintf(fp, "%s\t%s\n", G->_label_v1[i], G->_label_v2[j]);
	  }
	}
  }
  return;
}


/* Write out the degree of each vertex in a graph to a file pointer */
void bigraph_degreelist_out(FILE *fp, BiGraph *G)
{
  unsigned int i, n1, n2;
  n1 = bigraph_num_v1(G);
  n2 = bigraph_num_v2(G);
  fprintf(fp, "DegreeList %d %d\n", n1, n2);
  for (i = 0; i < n1; i++)
    fprintf(fp, "%s\t%d\n", G->_label_v1[i], (int)bigraph_degree_v1(G, i));
  for (i = 0; i < n2; i++)
    fprintf(fp, "%s\t%d\n", G->_label_v2[i], (int)bigraph_degree_v2(G, i));
  return;
}
 

/* Write out the adjacency list of each vertex to a file pointer */
void bigraph_adjlist_out(FILE *fp, BiGraph *G)
{
  unsigned int n1 = bigraph_num_v1(G);
  unsigned int n2 = bigraph_num_v2(G);
  unsigned int i, j, v;
  fprintf(fp, "AdjList %d %d\n", n1, n2);
  for (i = 0; i < n1; i++) {
    fprintf(fp, "%s :", G->_label_v1[i]);
    for (j = 0; j < G->_degree_v1[i]; j++) {
      v = G->_adj_list_v1[i][j];
      fprintf(fp, " %s", G->_label_v2[v]);
    }
    fprintf(fp, "\n");
  }
  return;
}

