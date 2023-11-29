#ifndef MBEA_COMMBE_H
#define MBEA_COMMBE_H

#include "BicliqueFinder.h"
#include "Utility.h"
#include <bitset>
#include <unordered_map>

#define N 64
typedef std::bitset<N> bitset_t;

class ComMBEFinder : public BicliqueFinder {
public:
  ComMBEFinder() = delete;
  ComMBEFinder(BiGraph *graph_in, const char *name = "ComMBEFinder");
  void Execute(int min_l_size = 1, int min_r_size = 1);

protected:
  void biclique_find(VertexSet &L, VertexSet &R, std::vector<Node *> &C);

  double r_clock_, iter_clock_, gen_cand_clock_, r_clock_mb_;
  int cur_v_;
private:
  void setup_bs(std::vector<int> &R);
  void biclique_find_bs(bitset_t R_bs, int R_start);
  void biclique_find_bs(std::vector<int> &L, bitset_t R_bs,
                        std::vector<int> &C_rev);

  std::vector<bitset_t> L_bitsets_;
  std::vector<std::vector<int>> R_neighbors_bs_;
  int R_size_bs_;

  std::vector<bool> l_valid;
  std::vector<bool> r_valid;
};

class ComMBEFinder2: public BicliqueFinder{
public:
  ComMBEFinder2() = delete;
  ComMBEFinder2(BiGraph *graph_in, const char *name = "ComMBEFinder2");
  void Execute(int min_l_size = 1, int min_r_size = 1);
protected:
 void biclique_find(VertexSet &L, VertexSet &R, std::vector<Node *> &C);
 double r_clock_;
 std::vector<VertexSet> global_2d_buf_;
};

class ComMBEFinder3 : public BicliqueFinder {
public:
  ComMBEFinder3() = delete;
  ComMBEFinder3(BiGraph *graph_in, const char *name = "ComMBEFinder3_bs_64");
  void Execute(int min_l_size = 1, int min_r_size = 1);

protected:
  void biclique_find(VertexSet &L, VertexSet &R,
                     std::vector<std::pair<int, VertexSet>> &C);
  void biclique_find(bitset_t &L, VertexSet &R,
                     std::vector<std::pair<int, bitset_t>> &Q_C, int c_start);

  std::vector<VertexSet> global_2d_buf_;
  std::vector<int> idx_buf_;
};

#endif