#ifndef MBEA_QUASI_H
#define MBEA_QUASI_H

#include "BicliqueFinder.h"
#include "Utility.h"
#include <bitset>
#include <unordered_map>

enum StateEnum {P, Q, D};
struct VertexExtNode {
  int v;
  int LN;
  StateEnum state;
};

typedef std::vector<VertexExtNode> VertexSetExt;

class QuasiMBEFinder : public BicliqueFinder {
public:
  QuasiMBEFinder() = delete;
  QuasiMBEFinder(BiGraph *graph_in, const char *name = "QuasiMBEFinder");
  void Execute(int min_l_size = 1, int min_r_size = 1);
  void setMiu(double miu_in);

 private:
  void biclique_find(VertexSet &L, VertexSet &R, VertexSetExt &C);
  void updateMB(int l_mb_in, int r_mb_in);
  double miu_;
  int l_mb_, r_mb_;
};

#endif