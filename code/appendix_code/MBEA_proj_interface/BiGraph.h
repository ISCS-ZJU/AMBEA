#ifndef MBEA_BIGRAPH_H
#define MBEA_BIGRAPH_H

#include <map>
#include <bitset>

#include "utility.h"
#include <queue>
const int kMaxSize = 0x40000;

class BiGraph {
 public:
  BiGraph();
  BiGraph(const BiGraph& graph);
  virtual ~BiGraph();
  virtual void LoadGraphFromAdjListFile(const char* filename,
                                        int omit_nums = 0, int l_bound = 0);
  virtual void LoadGraphFromBiGraph(BiGraph* graph, int l_bound = 0,
                                    int r_bound = 0);

  virtual void Reorder();
  void Prune(int l_bound, int r_bound);

  std::vector<int>& GetNeighborsOfLVertex(int v);
  std::vector<int>& GetNeighborsOfRVertex(int v);
  int GetLSize();
  int GetRSize();
  int GetLBound();
  int GetRBound();

  virtual bool EdgeExist(int l, int r);
  inline virtual void EdgeErase(int l, int r);

  void PrintProfile();
  void PrintTotalGraph();
  void ConvertIdVector(std::vector<int>& l_vec, std::vector<int>& r_vec);

  friend class MaximumEdgeBicliqueFinder;
 protected:
  virtual void Transpose();
  inline void Clear();
  void PruneCore(BiGraph* graph, int l_bound, int r_bound);

  bool transpose_flag_ = false;
  int l_bound_ = 0;
  int r_bound_ = 0;
  int edges_ = 0;

  // logic id -> real id
  std::vector<int> l_map_;  
  std::vector<int> r_map_;
  // all the adj list stores logic ids
  std::vector<std::vector<int>> l_adj_list_;  
  std::vector<std::vector<int>> r_adj_list_;

};



class BitBiGraph :virtual public BiGraph{
public:
  BitBiGraph();
  virtual ~BitBiGraph();
  void LoadGraphFromAdjListFile(const char* filename, int omit_nums = 0,
                                        int l_bound = 0);
  void LoadGraphFromBiGraph(BiGraph* graph, int l_bound = 0,
                                    int r_bound = 0);
  void Reorder();
  inline bool EdgeExist(int l, int r);
  inline void EdgeErase(int l, int r);

protected:
  void Transpose();
  //std::vector<std::bitset<kMaxSize>> bitset_;
  std::vector<MyBitset> bitset_;
  inline void BuildBitset();
};

class PMbeaBiGraph : virtual public BiGraph {
public:
  virtual ~PMbeaBiGraph();
  virtual void LoadGraphFromAdjListFile(const char* filename, int omit_nums = 0,
                                int l_bound = 0);
  virtual void LoadGraphFromBiGraph(BiGraph* graph, int l_bound = 0, int r_bound = 0);
  std::vector<std::pair<int, int>> range_index_;
private:
  void RevTopReorder();
};

class PMbeaBitBiGraph : public PMbeaBiGraph, BitBiGraph {
public:
  virtual~ PMbeaBitBiGraph();
  void LoadGraphFromAdjListFile(const char* filename, int omit_nums = 0,
                                int l_bound = 0);
  void LoadGraphFromBiGraph(BiGraph* graph, int l_bound = 0, int r_bound = 0);
  
  void Reorder();
  bool EdgeExist(int l, int r);
  void EdgeErase(int l, int r);
protected:
  void Transpose();
};


struct Biclique {
  std::vector<int> left_nodes, right_nodes;
  int edges;
  BiGraph* graph;
  void Print(FILE* fp = NULL);
};

#endif  // !MBEA_BIGRAPH_H
