#pragma once

#include "BicliqueFinder.h"
#include <queue>

class PmbeBiGraph : public BiGraph {
public:
  PmbeBiGraph() = delete;
  PmbeBiGraph(const BiGraph &graph);
  std::pair<int, int> GetRangeIndex(int v);

private:
  std::vector<std::pair<int, int>> range_index_;
};

class PmbeFinder : public BicliqueFinder {
public:
  PmbeFinder() = delete;
  PmbeFinder(BiGraph *graph_in, const char *name = "PmbeFinder");
  void Execute(int min_l_size = 1, int min_r_size = 1);

private:
  void biclique_find(std::vector<int> L, std::vector<int> R, std::vector<int> P,
                     std::vector<int> Q);
  PmbeBiGraph *pgraph_;
};

class PmbeFinderV2 : public BicliqueFinder { // disable vertex Q technique
public:
  PmbeFinderV2() = delete;
  PmbeFinderV2(BiGraph *graph_in, const char *name = "PmbeFinderV2");
  void Execute(int min_l_size = 1, int min_r_size = 1);

private:
  void biclique_find(std::vector<int> L, std::vector<int> R,
                     std::vector<int> C);
  PmbeBiGraph *pgraph_;
};

class PmbeAdvFinder : public BicliqueFinder { // disable vertex Q technique
public:
  PmbeAdvFinder() = delete;
  PmbeAdvFinder(BiGraph *graph_in, const char *name = "PmbeAdvFinder");
  void Execute(int min_l_size = 1, int min_r_size = 1);

private:
  void biclique_find(std::vector<int> L, std::vector<int> R,
                     std::vector<int> C);
  PmbeBiGraph *pgraph_;
};
