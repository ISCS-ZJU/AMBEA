#ifndef AGGRESSIVE_FINDER_H
#define AGGRESSIVE_FINDER_H
#include "BicliqueFinder.h"

/**
 * @brief 
 * final version: AmbeaFinder
 * parallel version: ParAmbeaFinder
 * disable merging prune: AmbeaFinderIntra
 * disable AE-tree / pure merging prune: IntraFinder (active)
 * naive version: AmbeaFinderNaive
 * pure passive merging prune: IntraFinder (passive)
 */




class AmbeaFinderNaive : public BicliqueFinder {
public:
  AmbeaFinderNaive() = delete;
  AmbeaFinderNaive(BiGraph *graph, const char *name = "AmbeaFinderNaive");
  void Execute(int min_l_size = 1, int min_r_size = 1);

protected:
  void biclique_find(const VertexSet &L, const VertexSet &R, const VertexSet &C,
                     int v);
};

class AmbeaFinderInter : public BicliqueFinder {
public:
  AmbeaFinderInter() = delete;
  AmbeaFinderInter(BiGraph *graph, const char *name = "AmbeaFinderInter");
  void Execute(int min_l_size = 1, int min_r_size = 1);

protected:
  void biclique_find(const VertexSet &L, const VertexSet &R,
                     const std::vector<std::pair<int, int>> &C, int v);
};

class IntraFinder : public BicliqueFinder {
public:
  IntraFinder() = delete;
  IntraFinder(BiGraph *graph, bool mode = true, const char *name = "IntraFinder");
  void Execute(int min_l_size = 1, int min_r_size = 1);

 protected:
  VertexSet Partition(const VertexSet &C, const VertexSet &L_prime, int v = -1);
  void biclique_find(const VertexSet &L, const VertexSet &R,
                     const VertexSet &C);
  void biclique_find_passive(const VertexSet &L, const VertexSet &R,
                             VertexSet &C);
  bool active_node_;
};

class AmbeaFinder : public BicliqueFinder {
public:
  AmbeaFinder() = delete;
  AmbeaFinder(BiGraph *graph_in, const char *name = "AmbeaFinder");
  void Execute(int min_l_size = 1, int min_r_size = 1);
  void SetOrder(OrderEnum order);

 protected:
  void biclique_find(const std::vector<int> &L, const std::vector<int> &R,
                     std::vector<std::pair<int, int>> &C, int vp, int vc_id);
  inline void Partition(const std::vector<int> &L,
                        const std::vector<int> &L_current,
                        std::vector<int> &L_prime,
                        std::vector<std::pair<int, int>> &C,
                        std::vector<std::pair<int, int>> &C_prime,
                        std::vector<int> &cand_ids, int vp, int vc);
  inline void PartitionHash(const std::vector<int> &L,
                            const std::vector<int> &L_current,
                            std::vector<int> &L_prime,
                            std::vector<std::pair<int, int>> &C,
                            std::vector<int> &R_add,
                            std::vector<std::pair<int, int>> &C_prime,
                            std::vector<int> &cand_ids, int vp, int vc);
  inline void PartitionBitset(const std::vector<int> &L,
                              const std::vector<int> &L_current,
                              std::vector<int> &L_prime,
                              std::vector<std::pair<int, int>> &C,
                              std::vector<int> &R_add,
                              std::vector<std::pair<int, int>> &C_prime,
                              std::vector<int> &cand_ids, int vp, int vc);
  std::vector<int> cand_rev_vec_;
  std::vector<bool> quick_check_bitset_;
  double my_clock_;
  OrderEnum order_;
};
class ParAmbeaFinder : public AmbeaFinder {
public:
  ParAmbeaFinder() = delete;
  ParAmbeaFinder(BiGraph *graph_in, int thread_num = 64,
                     const char *name = "ParAmbeaFinder");
  void Execute(int min_l_size = 1, int min_r_size = 1);

private:
  void biclique_find(const std::vector<int> &L, const std::vector<int> &R,
                     std::vector<std::pair<int, int>> C, int vp, int vc_id);
  void parallel_biclique_find(const std::vector<int> &L,
                              const std::vector<int> &R,
                              std::vector<std::pair<int, int>> C, int vp,
                              int vc_id);
  int thread_num;
};

#endif