#ifndef MBEA_BICLIQUE_FINDER_H
#define MBEA_BICLIQUE_FINDER_H

#include "BiGraph.h"
#include <atomic>
#include <cstring>
#include <map>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/parallel_for_each.h>
#include <queue>
#include <set>
#include <stack>
#include <tbb/global_control.h>
#include <unordered_map>

#define COMPUTE_LEVEL

struct CExtNode {
  int r_id;
  int nc;
  std::vector<int> r_cands;
  CExtNode(int r_id, int nc) : r_id(r_id), nc(nc) {
    r_cands.emplace_back(r_id);
  }
  CExtNode(int r_id, int nc, std::vector<int> r_cands)
      : r_id(r_id), nc(nc), r_cands(r_cands) {}
  CExtNode() : r_id(-1), nc(0), r_cands(0) {}
};

struct PNode {
  int start;
  int size;
  int nc;
  PNode(int start, int size, int nc) : start(start), size(size), nc(nc) {}
};

class BicliqueFinder {
public:
  BicliqueFinder() = delete;
  BicliqueFinder(BiGraph *graph_in, const char *name);
  virtual void Execute(int min_l_size = 1, int min_r_size = 1) = 0;
  void PrintResult(char *fn = nullptr);
  void SetMemory(double memory_usage);
  void SetThreads(int threads);
  friend class MEBFinder;

 protected:
  BiGraph *graph_;
  Biclique maximum_biclique_;
  char *finder_name_;
  std::atomic<long long int> processing_nodes_, maximal_nodes_;
  int min_l_size_, min_r_size_;
  double exe_time_, start_time_;
  bool is_transposed_;
  double memory_usage_;
  int threads_;
#ifdef COMPUTE_LEVEL
  int max_level_;
  long long int level_accumulation_;
  int cur_level_;
#endif

  inline void Partition(const std::vector<int> &L_prime,
                        const std::vector<CExtNode> &C,
                        std::vector<int> &reordered_map, std::vector<PNode> &P);
  inline void Partition(const std::vector<int> &L_prime,
                        const std::vector<std::pair<int, int>> &C,
                        std::vector<int> &reordered_map, std::vector<PNode> &P);
  inline void Partition(const std::vector<int> &L_prime,
                        const std::vector<int> &C,
                        std::vector<int> &reordered_map,
                        std::vector<std::pair<int, int>> &P);
  void setup(int min_l_size, int min_r_size);
  void finish();

  /**
   * @brief
   *
   * @param X stands for R
   * @param GamaX stands for L
   * @param tailX stands for C
   */
  void MineLMBC(std::vector<int> X, std::vector<int> GamaX,
                std::vector<int> tailX);
};

#endif
