#ifndef MBEA_BICLIQUE_FINDER_FAST_H
#define MBEA_BICLIQUE_FINDER_FAST_H

#include "BicliqueFinder.h"
#include "Utility.h"
#include <atomic>
#include <map>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/parallel_for_each.h>
#include <stack>
#include <tbb/global_control.h>
#include <unordered_map>

class MmbeaIntraFinderFast : public BicliqueFinder {
public:
  MmbeaIntraFinderFast() = delete;
  MmbeaIntraFinderFast(BiGraph *graph_in,
                       const char *name = "MmbeaIntraFinderFast");
  void Execute(int min_l_size = 1, int min_r_size = 1);

private:
  void biclique_find(const std::vector<int> &L, const std::vector<int> &R,
                     const std::vector<int> &C, int vc_id);
};

class MmbeaFinderFast : public BicliqueFinder {
public:
  MmbeaFinderFast() = delete;
  MmbeaFinderFast(BiGraph *graph_in, const char *name = "MmbeaFinderFast");
  void Execute(int min_l_size = 1, int min_r_size = 1);

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
};
class ParMmbeaFinderFast : public MmbeaFinderFast {
public:
  ParMmbeaFinderFast() = delete;
  ParMmbeaFinderFast(BiGraph *graph_in, int thread_num = 64,
                     const char *name = "ParMmbeaFinderFastBitset");
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

typedef unsigned int BITSET_T;
// suppose enumerate all MBE
class MbeaFinderBitset : public BicliqueFinder {
public:
  MbeaFinderBitset() = delete;
  MbeaFinderBitset(BiGraph *graph_in, const char *name = "MbeaFinderBitset");
  void Execute(int min_l_size = 1, int min_r_size = 1);

private:
  void biclique_find(const std::vector<int> &L, const std::vector<int> &R,
                     const std::vector<int> &C, int vc);
  void biclique_find(BITSET_T parent_L, BITSET_T current_L, int vvc,
                     std::vector<BITSET_T> &c_bs);
  double my_clock_;
};

class MbeaFinderAuto : public BicliqueFinder {
public:
  MbeaFinderAuto() = delete;
  MbeaFinderAuto(BiGraph *graph_in, const char *name = "MbeaFinderAuto");
  void Execute(int min_l_size = 1, int min_r_size = 1);

protected:
  void biclique_find(const std::vector<int> &L, int vc, int vc_id);
  void biclique_find(BITSET_T parent_L, BITSET_T current_L, int vvc,
                     std::vector<BITSET_T> &c_bs);
  int max_c_size_;
  double my_clock_;
};
class ParMbeaFinderAuto : public MbeaFinderAuto {
public:
  ParMbeaFinderAuto() = delete;
  ParMbeaFinderAuto(BiGraph *graph_in, int thread_num = 64,
                    const char *name = "ParMbeaFinderAuto");
  void Execute(int min_l_size = 1, int min_r_size = 1);

private:
  void parallel_biclique_find(const std::vector<int> &L, int vc, int vc_id);
  void parallel_biclique_find(BITSET_T parent_L, BITSET_T current_L, int vvc,
                              std::vector<BITSET_T> &c_bs);
  int thread_num;
};

#endif