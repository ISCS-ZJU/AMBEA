#ifndef BASELINE_FINDER_H
#define BASELINE_FINDER_H
#include "BicliqueFinder.h"

/**
 * @brief 
 * include optimized code for
 * MineLMBC (Dawak 06)
 * MBEA, iMBEA (BMC 14)
 * FMBE (HiPC 19)
 * PMBE (IJCAI 20)
 */

class MbeaFinder : public BicliqueFinder {
public:
  MbeaFinder() = delete;
  MbeaFinder(BiGraph *graph_in, const char *name = "MbeaFinder");
  void Execute(int min_l_size = 1, int min_r_size = 1);

private:
  void biclique_find(std::vector<int> L, std::vector<int> R, std::vector<int> P,
                     std::vector<int> Q);
};

class ImbeaFinder : public BicliqueFinder {
public:
  ImbeaFinder() = delete;
  ImbeaFinder(BiGraph *graph_in, const char *name = "ImbeaFinder");
  void Execute(int min_l_size = 1, int min_r_size = 1);

private:
  void biclique_find(std::vector<int> L, std::vector<int> R, std::vector<int> P,
                     std::vector<int> Q);
};

class MineLMBCFinder : public BicliqueFinder {
public:
  MineLMBCFinder() = delete;
  MineLMBCFinder(BiGraph *graph_in, const char *name = "MineLMBCFinder");
  void Execute(int min_l_size = 1, int min_r_size = 1);
};

class FmbeFinder : public BicliqueFinder {
public:
  FmbeFinder() = delete;
  FmbeFinder(BiGraph *graph_in, const char *name = "FmbeFinder");
  void Execute(int min_l_size = 1, int min_r_size = 1);
};

class PmbeBiGraph : public BiGraph {
public:
  PmbeBiGraph() = delete;
  PmbeBiGraph(const BiGraph &graph);
  std::pair<int, int> GetRangeIndex(int v);

private:
  std::vector<std::pair<int, int>> range_index_;
};

class PmbeFinder : public BicliqueFinder { // disable vertex Q technique
public:
  PmbeFinder() = delete;
  PmbeFinder(BiGraph *graph_in, const char *name = "PmbeFinder");
  void Execute(int min_l_size = 1, int min_r_size = 1);

private:
  void biclique_find(std::vector<int> L, std::vector<int> R,
                     std::vector<int> C);
  PmbeBiGraph *pgraph_;
};

#endif