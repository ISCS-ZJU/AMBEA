#ifndef MBEA_BIGRAPH_H
#define MBEA_BIGRAPH_H

#include <cstdlib>
#include <iostream>
#include <vector>

#include "Utility.h"

typedef std::vector<int> VertexSet;
struct Node {
  int vc;
  VertexSet L;
  VertexSet R;
  Node(int vc_in) : vc(vc_in) {}
  Node(int vc_in, VertexSet L_in) : vc(vc_in), L(L_in) {}
};

enum OrderEnum { LInc, LDec, RInc, RDec, UC, LRInc, Rand };

typedef std::vector<Node *> NodeSet;

typedef unsigned int BITSET_T;
class BiGraph {
public:
  BiGraph() = delete;
  BiGraph(const char *filename);
  BiGraph(const BiGraph &graph);
  ~BiGraph();

  void Reorder(OrderEnum op = LInc);
  void Transpose();

  int GetLSize();
  int GetRSize();
  int GetEdges();

  std::vector<int> &NeighborsL(int v);
  std::vector<int> NeighborsL(std::vector<int> &v_set);
  std::vector<int> &NeighborsR(int v);
  // bool EdgeExist(int l, int r);

  void PrintProfile();
  void PrintTotalGraph();
  void ConvertIdVector(std::vector<int> &l_vec, std::vector<int> &r_vec);

  void Prune1H(int min_neighbors_of_l, int min_neighbors_of_r);
  void Prune2H(int min_neighbors_of_l, int min_neighbors_of_r);
  void Prune2HOpt(int min_neighbors_of_l, int min_neighbors_of_r);
  void PopOrigGraph();
  
protected:
  void CalGraphStatistics();

  // int l_size_;
  // int r_size_;
  // int l_valid_;
  // int r_valid_;
  int min_neighbors_of_l_;
  int min_neighbors_of_r_;
  int l_degree_;
  int r_degree_;

  unsigned long long edges_;
  // logic id -> real id
  std::vector<int> l_map_;
  std::vector<int> r_map_;
  // all the adjlists stores logic id
  std::vector<std::vector<int>> l_adj_lists_;
  std::vector<std::vector<int>> r_adj_lists_;

  // graph buffer
  int min_neighbors_of_l_buf_;
  int min_neighbors_of_r_buf_;
  int edges_buf_;

  std::vector<std::vector<int>> l_adj_lists_buf_;
  std::vector<std::vector<int>> r_adj_lists_buf_;
};

class Biclique {
public:
  Biclique();
  Biclique(const Biclique &biclique);
  void Print(FILE *fp = stdout);
  void CompareAndSet(std::vector<int> &left_in, std::vector<int> &right_in,
                     BiGraph *graph);
  void CompareAndSet(std::vector<std::pair<int, int>> &left_in,
                     std::vector<int> &right_in, BiGraph *graph);
  void CompareAndSet(int L_size, int R_size);
  void Reset();
  int GetEdges();

private:
  std::vector<int> left_nodes_, right_nodes_;
  int edges_;
};

#endif