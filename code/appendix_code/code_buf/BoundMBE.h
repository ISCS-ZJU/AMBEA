#ifndef MBEA_BOUNDMBE_H
#define MBEA_BOUNDMBE_H

#include "BicliqueFinder.h"
#include "Utility.h"

typedef std::vector<std::pair<int, int>> VertexExtSet;

class BoundMBEFinder : public BicliqueFinder {
public:
  BoundMBEFinder() = delete;
  BoundMBEFinder(BiGraph *graph_in, const char *name = "BoundMBEFinder");
  void Execute(int min_l_size = 1, int min_r_size = 1);

protected:
  void biclique_find(VertexExtSet &L, VertexSet &R, VertexExtSet &C);
  void trim_node(VertexExtSet &L, VertexSet &R, VertexExtSet &C);
  void trim_node_2(VertexExtSet &L, VertexSet &R, VertexExtSet &C);

  void trim_node_max_edge(VertexExtSet &L, VertexSet &R, VertexExtSet &C);









  bool maximality_check(VertexExtSet &L, VertexSet &R);
  void print_node(VertexExtSet &L, VertexSet &R, VertexExtSet &C);
  double node_trim_time_;
  double intersect_time_;
};

#endif
