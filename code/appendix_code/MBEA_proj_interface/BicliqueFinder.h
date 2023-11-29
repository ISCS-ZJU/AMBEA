#ifndef MBEA_BICLIQUE_FINDER_H
#define MBEA_BICLIQUE_FINDER_H

#include <assert.h>
#include "BiGraph.h"

class BicliqueFinder {
public:
  BicliqueFinder() = delete;
  BicliqueFinder(BiGraph* graph_in);
  virtual void Execute(int ctrl = 0) = 0;
  void PrintResult(char* fn, double time_find);
  void SetBound(int l_bound, int r_bound);
  const Biclique& GetMaximumBiclique();
  friend class MaximumEdgeBicliqueFinder;

protected:
  BiGraph* graph_;
  Biclique maximum_biclique_;
  int processing_nodes_, maximal_nodes_;
  int l_bound_, r_bound_;
};

class BasicMbeaFinder : public BicliqueFinder {
public:
  BasicMbeaFinder() = delete;
  BasicMbeaFinder(BiGraph* graph_in);
  void Execute(int ctrl = 0);

private:
  void biclique_find(std::vector<int> L, std::vector<int> R, std::vector<int> P,
    std::vector<int> Q);
};

class ImproveMbeaFinder : public BicliqueFinder {
public:
  ImproveMbeaFinder() = delete;
  ImproveMbeaFinder(BiGraph* graph_in);
  void Execute(int ctrl = 0);

private:
  void biclique_find(std::vector<int> L, std::vector<int> R, std::vector<int> P,
    std::vector<int> Q);
};

#define MMBEA_FULL_MODE
#ifdef MMBEA_FULL_MODE
typedef std::vector<int> RType;
#else
typedef int RType;
#endif  // MMBEA_FULL_MODE

class MMbeaFinder : public BicliqueFinder {
public:
  MMbeaFinder() = delete;
  MMbeaFinder(BiGraph* graph_in);
  void Execute(int ctrl = 0);

private:
  struct CExtNode {
    int r_id;
    RType R;
    int l_cnt;
    CExtNode() {}
    CExtNode(int r_id, RType R, int l_cnt) : r_id(r_id), R(R), l_cnt(l_cnt) {}
  };

  void biclique_find(const std::vector<int>& L, const RType& R,
    const std::vector<CExtNode>& C, int parent_r_id, int r_id);
  void biclique_find_v2(const std::vector<int>& L, const RType& R,
    const std::vector<CExtNode>& C, int parent_r_id,
    int r_id);  // only with inter-pruning technique
  void finder_init(int ctrl = 0);   // ctrl = 0 inter+intra; ctrl = 1 inter
};

class MMbeaIntraFinder : public BicliqueFinder {
public:
  MMbeaIntraFinder() = delete;
  MMbeaIntraFinder(BiGraph* graph_in);
  void Execute(int ctrl = 0);

private:
  struct CExtNode {
    bool is_visited;
    std::vector<int> r_cands;
    CExtNode() :is_visited(false) {
      r_cands.clear();
    }
  };
  void biclique_find(const std::vector<int>& L, const std::vector<int>& R,
    std::vector<CExtNode>& C);
};

class PMbeaFinder : public BicliqueFinder {
public:
  PMbeaFinder() = delete;
  PMbeaFinder(PMbeaBiGraph* graph_in);
  void Execute(int ctrl = 0);

private:
  void biclique_find(std::vector<int> L, std::vector<int> R, std::vector<int> P,
    std::vector<int> Q);
};

class MaximumEdgeBicliqueFinder {
public:
  MaximumEdgeBicliqueFinder() = delete;
  MaximumEdgeBicliqueFinder(BiGraph* graph);
  void Execute(int ctrl = 0);
  void Print();

private:
  int max_edges;
  std::vector<int> max_l, max_r;
  BiGraph* orig_graph_;
};



#endif  // !MBEA_BICLIQUE_FINDER_H
