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

/* API for Bit-Based Adjacency Matrix for undirected bipartite graphs */
/* Graph Theory Team, Computer Science Department */ 
/* University of Tennessee, Knoxville */
/* Yun Zhang, yzhang@cs.utk.edu, September 2006 */
/* Updated: 12/17/2013 */

#ifndef __BIGRAPH_H
#define __BIGRAPH_H

#include <stdio.h>
#include <stdlib.h>
#include "utility.h"

/** Data structure **/

/* data structure of a bipartite graph */
typedef struct bipartite_graph_t {
  unsigned int _num_v1;  /* number of vertices in V1 */
  unsigned int _num_v2;  /* number of vertices in V2 */
  unsigned int _num_edges;       /* number of active edges */
  unsigned int _num_bytes_v1;    /* number of bytes for V1 */
  unsigned int _num_bytes_v2;    /* number of bytes for V2 */
  char ** _label_v1;     /* labels of vertices in V1 */
  char ** _label_v2;     /* labels of vertices in V2 */
  unsigned int **_neighbor_v1;       /* bit-based adjacency matrix for V1 */
  unsigned int **_neighbor_v2;       /* bit-based adjacency matrix for V2 */
  vid_t *_degree_v1;   /* number of edges each vertex have (<65536) */
  vid_t *_degree_v2;   /* number of edges each vertex have (<65536) */
  // added adjacency list 
  vid_t **_adj_list_v1;     /* adjacency list for V1 */
  vid_t **_adj_list_v2;     /* adjacency list for V2 */
} BiGraph;



/* Macros for BiGraph */

/* returns the total number of vertices in V1 of graph */
#define bigraph_num_v1(g)  (g->_num_v1)

/* returns the number of active vertices in V2 of graph */
#define bigraph_num_v2(g)  (g->_num_v2)

/* returns the number of active edges in graph */
#define bigraph_num_edges(g)  (g->_num_edges)

/* check to see if an edge exists or not */
#define bigraph_edge_exists(g, u, v)  (IS_SET(g->_neighbor_v2[v], u))

/* return the degree of a vertex in V1 */
#define bigraph_degree_v1(g, u)  (g->_degree_v1[u])

/* return the degree of a vertex in V1 */
#define bigraph_degree_v2(g, u)  (g->_degree_v2[u])

/* return the pointer to the bit-based neighbors of vertex u */
#define bigraph_neighbors(g, u)  (g->_neighbor_v2[u])

/* return the pointer to the bit-based neighbors of vertex u */
#define bigraph_neighbor_v1(g, u)  (g->_neighbor_v1[u])

/* add an edge to graph */
#define bigraph_add_edge(g, u, v)  { \
  if (!IS_SET(g->_neighbor_v2[v], u)) {\
    SET_BIT(g->_neighbor_v2[v], u); SET_BIT(g->_neighbor_v1[u], v);\
    g->_num_edges++; \
    g->_degree_v1[u]++; g->_degree_v2[v]++; }\
}	

/* delete an edge from graph */
#define bigraph_delete_edge(g, u, v) { \
  if (IS_SET(g->_neighbor_v2[v], u)) {\
    DEL_BIT(g->_neighbor_v2[v], u); DEL_BIT(g->_neighbor_v1[u], v); \
    g->_degree_v1[u]--; g->_degree_v2[v]--; \
    g->_num_edges--; } \
}


/** Function prototypes **/

/* Malloc a graph, initialize it and returns a pointer to it */
BiGraph *bigraph_make(unsigned int num_v1, unsigned int num_v2);

/* Free the memory of a graph */
void bigraph_free(BiGraph *G);

/* Read in a graph from an unweighted edge-list format file */
BiGraph *bigraph_edgelist_in(FILE *fp);

/* Read in a graph from a binary matrix format file */
BiGraph *bigraph_binarymatrix_in(FILE *fp);

/* Write out the degree of each vertex in a graph to a file pointer */
void bigraph_degreelist_out(FILE *fp, BiGraph *G);

/* Create adjacency lists, this is after read in a graph from a file */
void bigraph_adjlist_make(BiGraph *G);

/* Write out the adjacency list to a file */
void bigraph_adjlist_out(FILE *fp, BiGraph *G);


#endif  /* __BIGRAPH_H */

