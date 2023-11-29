#include "BiGraph.h"
#include <assert.h>

namespace {
std::vector<int> getLine(FILE* fp, int omit_num) {
  std::vector<int> list;
  char c;

  do {
    int item = 0, pos = 0;
    c = getc(fp);
    while ((c >= '0') && (c <= '9')) {
      item *= 10;
      item += int(c) - int('0');
      c = getc(fp);
      pos++;
    }
    if (pos) {
      if (omit_num > 0)
        --omit_num;
      else
        list.push_back(item);
    }
  } while (c != '\n' && !feof(fp));
  return list;
}
}  // namespace

BiGraph::BiGraph() {}

BiGraph::BiGraph(const BiGraph& graph) {
  l_bound_ = graph.l_bound_;
  r_bound_ = graph.r_bound_;
  edges_ = graph.edges_;
  transpose_flag_ = graph.transpose_flag_;
  l_map_ = graph.l_map_;
  r_map_ = graph.r_map_;

  for (auto& list : graph.l_adj_list_) {
    l_adj_list_.emplace_back(list);
  }
  for (auto& list : graph.r_adj_list_) {
    r_adj_list_.emplace_back(list);
  }
}

BiGraph::~BiGraph() { Clear(); }

void BiGraph::LoadGraphFromAdjListFile(const char* filename, int omit_nums,
                                       int l_bound) {
  printf("Building bigraph with file \"%s\" ... \n", filename);
  FILE* fp = fopen(filename, "rt");

  if (fp == nullptr) {
    printf("Could not open file %s. Bigraph building failed!\n", filename);
    exit(1);
  }

  Clear();
  l_bound_ = l_bound;
  r_bound_ = 0;
  edges_ = 0;
  transpose_flag_ = false;

  std::vector<int> neighbor_list;
  std::map<int, std::vector<int>> r_map_tmp;
  int left_id = 0;

  while (!feof(fp)) {
    neighbor_list = move(getLine(fp, omit_nums));
    if (neighbor_list.size() <= l_bound_) {
      if (neighbor_list.size() != 0) left_id++;
      continue;
    }
    edges_ += neighbor_list.size();
    int v_l_id = l_map_.size();
    l_map_.emplace_back(left_id++);
    for (int r_id : neighbor_list) r_map_tmp[r_id].emplace_back(v_l_id);
  }
  int r_size = r_map_tmp.size();
  r_map_.resize(r_size);
  r_adj_list_.resize(r_size);
  l_adj_list_.resize(l_map_.size());
  int v_r_id = 0;
  for (auto& p : r_map_tmp) {
    r_map_[v_r_id] = p.first;
    r_adj_list_[v_r_id] = std::move(p.second);
    for (auto v_l_id : r_adj_list_[v_r_id])
      l_adj_list_[v_l_id].emplace_back(v_r_id);
    v_r_id++;
  }
  r_map_tmp.clear();
  if (l_map_.size() < r_map_.size()) Transpose();
  printf("Bigraph building finished!\n\n");
}

void BiGraph::LoadGraphFromBiGraph(BiGraph* graph, int l_bound, int r_bound) {
  // printf("Building bigraph with graph ... \n");
  if (graph->transpose_flag_) std::swap(graph->l_bound_, graph->r_bound_);
  l_bound = std::max(l_bound, graph->l_bound_);
  r_bound = std::max(r_bound, graph->r_bound_);
  if (graph->transpose_flag_) std::swap(graph->l_bound_, graph->r_bound_);
  transpose_flag_ = graph->transpose_flag_;
  PruneCore(graph, l_bound, r_bound);
  if (l_map_.size() < r_map_.size()) Transpose();
  // printf("Bigraph building finished!\n\n");
}

void BiGraph::Transpose() {
  transpose_flag_ = ~transpose_flag_;
  std::swap(l_bound_, r_bound_);
  std::swap(l_map_, r_map_);
  std::swap(l_adj_list_, r_adj_list_);
}

void BiGraph::Reorder() {
  std::vector<int> reorder_map;
  reorder_map.resize(r_map_.size());
  for (int i = 0; i < r_map_.size(); i++) reorder_map[i] = i;
  std::sort(reorder_map.begin(), reorder_map.end(),
            [=](int id0, int id1) -> bool {
              return r_adj_list_[id0].size() < r_adj_list_[id1].size();
            });
  for (int i = 0; i < l_adj_list_.size(); i++) l_adj_list_[i].clear();
  std::vector<int> n_r_map;
  std::vector<std::vector<int>> n_r_adj_list;
  n_r_map.resize(r_map_.size());
  n_r_adj_list.resize(r_map_.size());
  for (int i = 0; i < r_map_.size(); i++) {
    n_r_map[i] = r_map_[reorder_map[i]];
    n_r_adj_list[i] = std::move(r_adj_list_[reorder_map[i]]);
    for (int l_id : n_r_adj_list[i]) {
      l_adj_list_[l_id].emplace_back(i);
    }
  }
  for (int i = 0; i < r_map_.size(); i++) r_adj_list_[i].clear();
  r_adj_list_.clear();
  r_map_.clear();
  std::swap(n_r_adj_list, r_adj_list_);
  std::swap(n_r_map, r_map_);

  // std::vector<int> reorder_map;
  // reorder_map.resize(l_map_.size());
  // for (int i = 0; i < reorder_map.size(); i++) {
  //  reorder_map[i] = i;
  //}

  // std::sort(reorder_map.begin(), reorder_map.end(),
  //          [=](int id0, int id1) -> bool {
  //            return l_adj_list_[id0].size() < l_adj_list_[id1].size();
  //          });
  // for (int i = 0; i < r_map_.size(); i++) r_adj_list_[i].clear();

  // std::vector<int> n_l_map;
  // std::vector<std::vector<int>> n_l_adj_list;
  // n_l_map.resize(l_map_.size());
  // n_l_adj_list.resize(l_map_.size());

  // for (int i = 0; i < l_map_.size(); i++) {
  //  n_l_map[i] = l_map_[reorder_map[i]];
  //  n_l_adj_list[i] = std::move(l_adj_list_[reorder_map[i]]);
  //  for (int r_id : n_l_adj_list[i]) {
  //    r_adj_list_[r_id].emplace_back(i);
  //  }
  //}
  // for (int i = 0; i < l_map_.size(); i++) l_adj_list_.clear();
  // l_adj_list_.clear();
  // l_map_.clear();

  // std::swap(l_map_, n_l_map);
  // std::swap(l_adj_list_, n_l_adj_list);
}

void BiGraph::Prune(int l_bound, int r_bound) {
  PruneCore(this, l_bound, r_bound);
}

std::vector<int>& BiGraph::GetNeighborsOfLVertex(int v) {
  return l_adj_list_[v];
}

std::vector<int>& BiGraph::GetNeighborsOfRVertex(int v) {
  return r_adj_list_[v];
}

int BiGraph::GetLSize() { return l_map_.size(); }

int BiGraph::GetRSize() { return r_map_.size(); }

int BiGraph::GetLBound() { return l_bound_; }

int BiGraph::GetRBound() { return r_bound_; }

bool BiGraph::EdgeExist(int l, int r) {
  return  // l >= 0 && l < GetLSize() && r >= 0 && r < GetRSize() &&
      std::binary_search(l_adj_list_[l].begin(), l_adj_list_[l].end(), r);
}

void BiGraph::EdgeErase(int l, int r) {
  if (l >= 0 && l < l_adj_list_.size() && r >= 0 && r < r_adj_list_.size()) {
    auto& l_list = l_adj_list_[l];
    auto& r_list = r_adj_list_[r];
    l_list.erase(find(l_list.begin(), l_list.end(), r));
    r_list.erase(find(r_list.begin(), r_list.end(), l));
  }
}

void BiGraph::PrintProfile() {
  printf("Bigraph Profile:\n");
  printf("Number of left vertices      : %d\n", (int)l_map_.size());
  printf("Number of right vertices     : %d\n", (int)r_map_.size());
  printf("Number of edges              : %d\n", edges_);
  printf("Min left degree bound        : %d\n", l_bound_ + 1);
  printf("Min right degree bound       : %d\n\n", r_bound_ + 1);
}

void BiGraph::PrintTotalGraph() {
  for (int i = 0; i < l_map_.size(); i++) {
    printf("[%d] ", l_map_[i]);
    for (int j = 0; j < l_adj_list_[i].size(); j++) {
      printf(" %d", r_map_[l_adj_list_[i][j]]);
    }
    printf("\n");
  }
  printf("\n");
}

void BiGraph::ConvertIdVector(std::vector<int>& l_vec,
                              std::vector<int>& r_vec) {
  for (int i = 0; i < l_vec.size(); i++) {
    if (l_vec[i] >= 0 && l_vec[i] < l_map_.size())
      l_vec[i] = l_map_[l_vec[i]];
    else
      l_vec[i] = -1;
  }
  std::sort(l_vec.begin(), l_vec.end());
  for (int i = 0; i < r_vec.size(); i++) {
    if (r_vec[i] >= 0 && r_vec[i] < r_map_.size())
      r_vec[i] = r_map_[r_vec[i]];
    else
      r_vec[i] = -1;
  }
  std::sort(r_vec.begin(), r_vec.end());
  if (transpose_flag_) {
    std::swap(l_vec, r_vec);
  }
}

inline void BiGraph::Clear() {
  l_map_.clear();
  r_map_.clear();
  for (auto vec : l_adj_list_) vec.clear();
  l_adj_list_.clear();
  for (auto vec : r_adj_list_) vec.clear();
  r_adj_list_.clear();
}

void BiGraph::PruneCore(BiGraph* graph, int l_bound, int r_bound) {
  if (transpose_flag_) {
    std::swap(l_bound, r_bound);
  }
  if (graph == this && l_bound <= l_bound_ && r_bound <= r_bound_) return;
  if (graph != this && l_bound <= graph->l_bound_ &&
      r_bound <= graph->r_bound_) {
    *this = *graph;
    return;
  }
  std::vector<int> l_degree(graph->l_map_.size());
  std::vector<int> r_degree(graph->r_map_.size());

  std::vector<int> l_pruned_seed, r_pruned_seed;

  for (int i = 0; i < graph->l_map_.size(); i++) {
    l_degree[i] = graph->l_adj_list_[i].size();
    if (l_degree[i] <= r_bound) l_pruned_seed.emplace_back(i);
  }

  for (int i = 0; i < graph->r_map_.size(); i++) {
    r_degree[i] = graph->r_adj_list_[i].size();
    if (r_degree[i] <= l_bound) r_pruned_seed.emplace_back(i);
  }

  while (!l_pruned_seed.empty() || !r_pruned_seed.empty()) {
    for (int l_id : l_pruned_seed) {
      for (int r_id : graph->l_adj_list_[l_id]) {
        r_degree[r_id]--;
        if (r_degree[r_id] == l_bound) r_pruned_seed.emplace_back(r_id);
      }
    }
    l_pruned_seed.clear();
    for (int r_id : r_pruned_seed) {
      for (int l_id : graph->r_adj_list_[r_id]) {
        l_degree[l_id]--;
        if (l_degree[l_id] == r_bound) l_pruned_seed.emplace_back(l_id);
      }
    }
    r_pruned_seed.clear();
  }

  std::vector<int> n_l_map;
  std::vector<int> n_r_map;
  std::vector<std::vector<int>> n_l_adj_list;
  std::vector<std::vector<int>> n_r_adj_list;

  int left_id = 0, right_id = 0;
  int edges = 0;

  for (int i = 0; i < l_degree.size(); i++) {
    if (l_degree[i] > r_bound) {
      edges += l_degree[i];
      l_degree[i] = left_id++;
      n_l_map.emplace_back(graph->l_map_[i]);
    } else
      l_degree[i] = -1;
  }
  for (int i = 0; i < r_degree.size(); i++) {
    if (r_degree[i] > l_bound) {
      r_degree[i] = right_id++;
      n_r_map.emplace_back(graph->r_map_[i]);
    } else
      r_degree[i] = -1;
  }
  n_l_adj_list.resize(n_l_map.size());
  n_r_adj_list.resize(n_r_map.size());
  for (int i = 0; i < l_degree.size(); i++) {
    if (l_degree[i] >= 0) {
      for (int r_id : graph->l_adj_list_[i]) {
        if (r_degree[r_id] >= 0) {
          n_l_adj_list[l_degree[i]].emplace_back(r_degree[r_id]);
          n_r_adj_list[r_degree[r_id]].emplace_back(l_degree[i]);
        }
      }
      std::sort(n_l_adj_list[l_degree[i]].begin(),
                n_l_adj_list[l_degree[i]].end());
    }
  }
  printf("[Pruning Core] Orignal edges:%d\t; New edges:%d (lb:%d, rb:%d)\n",
         graph->edges_, edges, l_bound + 1, r_bound + 1);
  Clear();
  edges_ = edges;
  l_bound_ = l_bound;
  r_bound_ = r_bound;
  std::swap(l_map_, n_l_map);
  std::swap(r_map_, n_r_map);
  std::swap(l_adj_list_, n_l_adj_list);
  std::swap(r_adj_list_, n_r_adj_list);
  if (l_map_.size() < r_map_.size()) Transpose();
}

void Biclique::Print(FILE* fp) {
  if (graph != NULL) {
    graph->ConvertIdVector(left_nodes, right_nodes);
    graph = NULL;
  }
  fprintf(fp, "L nodes[%d]:", (int)left_nodes.size());
  for (int node : left_nodes) fprintf(fp, "%d ", node);
  fprintf(fp, "\n");
  int r_size = left_nodes.size() == 0 ? 0 : edges / left_nodes.size();

  fprintf(fp, "R nodes:[%d]:", r_size);
  for (int node : right_nodes) fprintf(fp, "%d ", node);
  fprintf(fp, "\n\n");
}

BitBiGraph::BitBiGraph() : BiGraph() {}

BitBiGraph::~BitBiGraph() {
  Clear();
  bitset_.clear();
}

void BitBiGraph::LoadGraphFromAdjListFile(const char* filename, int omit_nums,
                                          int l_bound) {
  BiGraph::LoadGraphFromAdjListFile(filename, omit_nums, l_bound);
  BuildBitset();
}

void BitBiGraph::LoadGraphFromBiGraph(BiGraph* graph, int l_bound,
                                      int r_bound) {
  BiGraph::LoadGraphFromBiGraph(graph, l_bound, r_bound);
  BuildBitset();
}

void BitBiGraph::Transpose() {
  BiGraph::Transpose();
  BuildBitset();
}

void BitBiGraph::Reorder() {
  BiGraph::Reorder();
  BuildBitset();
}

inline bool BitBiGraph::EdgeExist(int l, int r) { return bitset_[l].test(r); }

inline void BitBiGraph::EdgeErase(int l, int r) {
  BiGraph::EdgeErase(l, r);
  bitset_[l].reset(r);
}

inline void BitBiGraph::BuildBitset() {
  bitset_.clear();
  const int left_size = l_map_.size();
  const int right_size = r_map_.size();
  bitset_.resize(left_size);
  for (int l = 0; l < left_size; l++) {
    for (int r : l_adj_list_[l]) {
      bitset_[l].resize(right_size);
      bitset_[l].set(r);
    }
  }
}

PMbeaBiGraph::~PMbeaBiGraph() { 
  Clear();
 }

void PMbeaBiGraph::LoadGraphFromAdjListFile(const char* filename, int omit_nums,
                                            int l_bound) {
  BiGraph::LoadGraphFromAdjListFile(filename, omit_nums, l_bound);
  RevTopReorder();
}

void PMbeaBiGraph::LoadGraphFromBiGraph(BiGraph* graph, int l_bound,
                                        int r_bound) {
  BiGraph::LoadGraphFromBiGraph(graph, l_bound, r_bound);
  RevTopReorder();
}

void PMbeaBiGraph::RevTopReorder() {
  const int r_size = r_map_.size();
  std::vector<int> in_degree(r_size, 0);
  std::vector<std::vector<int>> child_r_nodes(r_size);
  std::vector<int> rename_map(r_size);
  typedef std::pair<int, int> IDPAIR;
  std::vector<int> ready_buf;
  auto rid_cmp = [=](int id0, int id1) -> bool {
    return r_adj_list_[id0].size() > r_adj_list_[id1].size();
  };
  std::queue<int> ready_q;
  range_index_.clear();
  range_index_.resize(r_size);
  for (int i = 0; i < r_size; i++) {
    for (int j = i + 1; j < r_size; j++) {
      int contain_res = seq_contain(r_adj_list_[i], r_adj_list_[j]);
      if (contain_res < 0) {  //  i contains j
        in_degree[j]++;
        child_r_nodes[i].emplace_back(j);
      } else if (contain_res > 0) {
        in_degree[i]++;
        child_r_nodes[j].emplace_back(i);
      }
    }
  }
  int scan_id = r_size;
  for (int i = 0; i < r_size; i++) {
    if (in_degree[i] == 0)  // ready_q.push(i);
      ready_buf.emplace_back(i);
  }
  std::sort(ready_buf.begin(), ready_buf.end(), rid_cmp);
  while (!ready_buf.empty()) {
    ready_q.push(ready_buf.back());
    ready_buf.pop_back();
  }
  while (!ready_q.empty()) {
    int r_id = ready_q.front();
    ready_q.pop();
    rename_map[--scan_id] = r_id;
    int start_point, end_point;
    start_point = scan_id - ready_q.size();
    end_point = start_point - 1;
    for (int neighbor : child_r_nodes[r_id]) {
      if (--in_degree[neighbor] == 0) {
        start_point--;
        ready_buf.emplace_back(neighbor);
        // ready_q.push(neighbor);
      }
    }
    std::sort(ready_buf.begin(), ready_buf.end(), rid_cmp);
    while (!ready_buf.empty()) {
      ready_q.push(ready_buf.back());
      ready_buf.pop_back();
    }
    range_index_[scan_id].first = start_point;
    range_index_[scan_id].second = end_point;
  }

  std::vector<int> n_r_map(r_size);
  std::vector<std::vector<int>> n_r_adj_list(r_size);
  std::vector<std::vector<int>> n_l_adj_list(l_map_.size());
  for (int i = 0; i < r_size; i++) {
    int r_id = rename_map[i];
    for (int l : r_adj_list_[r_id]) n_l_adj_list[l].emplace_back(i);
    n_r_map[i] = r_map_[r_id];
    n_r_adj_list[i] = std::move(r_adj_list_[r_id]);
  }
  std::swap(n_r_map, r_map_);
  std::swap(n_r_adj_list, r_adj_list_);
  std::swap(n_l_adj_list, l_adj_list_);
  n_l_adj_list.clear();
  n_r_adj_list.clear();
}

PMbeaBitBiGraph::~PMbeaBitBiGraph() { Clear(); bitset_.clear(); }

void PMbeaBitBiGraph::LoadGraphFromAdjListFile(const char* filename,
                                               int omit_nums, int l_bound) {
  PMbeaBiGraph::LoadGraphFromAdjListFile(filename, omit_nums, l_bound);
  BitBiGraph::BuildBitset();
}

void PMbeaBitBiGraph::LoadGraphFromBiGraph(BiGraph* graph, int l_bound,
                                           int r_bound) {
  PMbeaBiGraph::LoadGraphFromBiGraph(graph, l_bound, r_bound);
  BitBiGraph::BuildBitset();
}

void PMbeaBitBiGraph::Transpose() { BitBiGraph::Transpose(); }

void PMbeaBitBiGraph::Reorder() { BitBiGraph::Reorder(); }

bool PMbeaBitBiGraph::EdgeExist(int l, int r) {
  return BitBiGraph::EdgeExist(l, r);
}

void PMbeaBitBiGraph::EdgeErase(int l, int r) { BitBiGraph::EdgeErase(l, r); }
