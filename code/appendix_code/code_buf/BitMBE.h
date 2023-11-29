#ifndef MBEA_BITMBE_H
#define MBEA_BITMBE_H

#include "BicliqueFinder.h"
#include "Utility.h"

class BitMBEFinder : public BicliqueFinder {
public:
  BitMBEFinder() = delete;
  BitMBEFinder(BiGraph *graph_in, const char *name = "BitMBEFinder");
  void Execute(int min_l_size = 1, int min_r_size = 1);

protected:
  void biclique_find(const VertexSet &L, const VertexSet &R,
                     std::vector<std::pair<int, int>> &P_Q, // with LN size
                     int p_start);
  void biclique_find(BITSET_T L_bs, std::vector<BITSET_T> P_Q_bs, int p_start);
  void biclique_find_2(BITSET_T L_bs, int L_start, int max_L,
                       std::vector<BITSET_T> P_Q_bs, int p_start);
  double bitset_time_;
};

#endif