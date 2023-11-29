#ifndef MBEA_BICLIQUE_FINDER_H
#define MBEA_BICLIQUE_FINDER_H

#include "BiGraph.h"
#include <atomic>
#include <cstring>
#include <map>
#include <set>
#include <stack>

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
  friend class MEBFinder;

 protected:
  BiGraph *graph_;
  Biclique maximum_biclique_;
  char *finder_name_;
  std::atomic<long long int> processing_nodes_, maximal_nodes_;
  int min_l_size_, min_r_size_;
  double exe_time_, start_time_;
  bool is_transposed_;
  int max_level_;

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

class MbeaFinder : public BicliqueFinder {
public:
  MbeaFinder() = delete;
  MbeaFinder(BiGraph *graph_in, const char *name = "MbeaFinder");
  void Execute(int min_l_size = 1, int min_r_size = 1);

private:
  void biclique_find(std::vector<int> L, std::vector<int> R, std::vector<int> P,
                     std::vector<int> Q);
};

class MbeaAdvFinder : public BicliqueFinder {
public:
  MbeaAdvFinder() = delete;
  MbeaAdvFinder(BiGraph *graph_in, const char *name = "MbeaAdvFinder");
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

class ImbeaAdvFinder : public BicliqueFinder {
public:
  ImbeaAdvFinder() = delete;
  ImbeaAdvFinder(BiGraph *graph_in, const char *name = "ImbeaAdvFinder");
  void Execute(int min_l_size = 1, int min_r_size = 1);

private:
  void biclique_find(std::vector<int> L, std::vector<int> R, std::vector<int> P,
                     std::vector<int> Q, int level = 0);
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

class MmbeaIntraFinder : public BicliqueFinder {
public:
  MmbeaIntraFinder() = delete;
  MmbeaIntraFinder(BiGraph *graph_in, const char *name = "MmbeaIntraFinder");
  void Execute(int min_l_size = 1, int min_r_size = 1);

private:
  void biclique_find(const std::vector<int> &L, const std::vector<int> &R,
                     std::vector<CExtNode> &C, int vc);
};

class MmbeaIntraFinderV2
    : public BicliqueFinder { // only keep leading cand vertex
public:
  MmbeaIntraFinderV2() = delete;
  MmbeaIntraFinderV2(BiGraph *graph_in,
                     const char *name = "MmbeaIntraFinderV2");
  void Execute(int min_l_size = 1, int min_r_size = 1);

private:
  void biclique_find(const std::vector<int> &L, const std::vector<int> &R,
                     std::vector<int> &C, int vc);
};

class MmbeaInterFinder : public BicliqueFinder {
public:
  MmbeaInterFinder() = delete;
  MmbeaInterFinder(BiGraph *graph_in, const char *name = "MmbeaInterFinder");
  void Execute(int min_l_size = 1, int min_r_size = 1);

private:
  void biclique_find(const std::vector<int> &L, const std::vector<int> &R,
                     const std::vector<std::pair<int, int>> &C, int vp, int vc);
};

class MmbeaInterAdvFinder : public BicliqueFinder {
public:
  MmbeaInterAdvFinder() = delete;
  MmbeaInterAdvFinder(BiGraph *graph_in,
                      const char *name = "MmbeaInterAdvFinder");
  void Execute(int min_l_size = 1, int min_r_size = 1);

private:
  void biclique_find(const std::vector<int> &L, const std::vector<int> &R,
                     const std::vector<std::pair<int, int>> &C, int vp, int vc);
};

// partition and merge the vertex set to one vertex

class MmbeaFinder : public BicliqueFinder {
public:
  MmbeaFinder() = delete;
  MmbeaFinder(BiGraph *graph_in, const char *name = "MmbeaFinder");
  void Execute(int min_l_size = 1, int min_r_size = 1);

private:
  void biclique_find(const std::vector<int> &L, const std::vector<int> &R,
                     std::vector<std::pair<int, int>> &C, int vp, int vc);
};

struct CNode {
  int r_id, nc, size;
  CNode(int r_id, int nc, int size) : r_id(r_id), nc(nc), size(size) {}
};

class MmbeaFinderV2 : public BicliqueFinder {
public:
  MmbeaFinderV2() = delete;
  MmbeaFinderV2(BiGraph *graph_in, const char *name = "MmbeaFinderV2");
  void Execute(int min_l_size = 1, int min_r_size = 1);

private:
  void biclique_find(const std::vector<int> &L, const std::vector<int> &R,
                     std::vector<CNode> &C, int vp, int vc);
  void Partition(const std::vector<int> &L_prime, const std::vector<CNode> &C,
                 std::vector<int> &reordered_map, std::vector<PNode> &P);
};

class MEBFinder {
public:
  MEBFinder() = delete;
  MEBFinder(BiGraph *graph);
  void Execute(int ctrl = 0);
  void PrintResult(char *fn = nullptr);

private:
  Biclique maximum_biclique_;
  BiGraph *orig_graph_;
  double exe_time_;
};

class LevelFinder : public BicliqueFinder {
public:
  LevelFinder() = delete;
  LevelFinder(BiGraph *graph, const char *name = "LevelFinder");
  void Execute(int min_l_size = 1, int min_r_size = 1);

protected:
  double recover_time_;
  std::vector<std::pair<int, int>> L_; // u_id, cur_level
  std::vector<int> R_;
  std::vector<std::pair<std::pair<int, int>, int>> C_; // v_id, nc, lock_level
  std::vector<int> cid_stack_;

  int FindNext(int last_cid);
  bool Push(int c_id);
  void Pop();
};

#define CACHE_SIZE 0x10000

class LevelFinderCache : public LevelFinder {
public:
  LevelFinderCache() = delete;
  LevelFinderCache(BiGraph *graph, const char *name = "LevelFinderCache");
  void Execute(int min_l_size = 1, int min_r_size = 1);

private:
  std::vector<int> cache_;
  bool CalUsingCache(int cur_cid);
};

class AMBEAFinderNaive : public BicliqueFinder {
public:
  AMBEAFinderNaive() = delete;
  AMBEAFinderNaive(BiGraph *graph, const char *name = "AMBEAFinderNaive");
  void Execute(int min_l_size = 1, int min_r_size = 1);

protected:
  void biclique_find(const VertexSet &L, const VertexSet &R, const VertexSet &C,
                     int v);
};

class AMBEAFinderInter : public BicliqueFinder {
public:
  AMBEAFinderInter() = delete;
  AMBEAFinderInter(BiGraph *graph, const char *name = "AMBEAFinderInter");
  void Execute(int min_l_size = 1, int min_r_size = 1);

protected:
  void biclique_find(const VertexSet &L, const VertexSet &R,
                     const std::vector<std::pair<int, int>> &C, int v);
};

class AMBEAFinderIntra : public BicliqueFinder {
public:
  AMBEAFinderIntra() = delete;
  AMBEAFinderIntra(BiGraph *graph, const char *name = "AMBEAFinderIntra");
  void Execute(int min_l_size = 1, int min_r_size = 1);

protected:
  VertexSet Partition(const VertexSet &C, const VertexSet &L_prime, int v = -1);

private:
  void biclique_find(const VertexSet &L, const VertexSet &R, const VertexSet &C,
                     int v);
};

class IntraFinder : public AMBEAFinderIntra {
public:
  IntraFinder() = delete;
  IntraFinder(BiGraph *graph, const char *name = "IntraFinder");
  void Execute(int min_l_size = 1, int min_r_size = 1);

protected:
  void biclique_find(const VertexSet &L, const VertexSet &R,
                     const VertexSet &C);
  void biclique_find_passive(const VertexSet &L, const VertexSet &R,
                             VertexSet &C);
};

class AMBEAFinder : public AMBEAFinderIntra {
public:
  AMBEAFinder() = delete;
  AMBEAFinder(BiGraph *graph, const char *name = "AMBEAFinder");
  void Execute(int min_l_size = 1, int min_r_size = 1);

protected:
  void biclique_find(const VertexSet &L, const VertexSet &R,
                     const VertexSet &C, int v);
  double my_clock_;
};

#endif
