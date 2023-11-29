#include "BiGraph.h"

#include <algorithm>
#include <cstdlib>
#include <queue>
#include <unordered_map>

namespace {
/**
 * @brief read a line from file to a vector and update the maximum item id
 * @param fp
 * @param max_item
 * @return
 */
std::vector<int> getLine(FILE *fp, int &max_item) {
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
      list.push_back(item);
      max_item = std::max(max_item, item);
    }
  } while (c != '\n' && !feof(fp));
  // std::sort(list.begin(), list.end());
  return list;
}
} // namespace

BiGraph::BiGraph(const char *filename)
    : min_neighbors_of_l_(1), min_neighbors_of_r_(1), l_map_(0), r_map_(0) {
  edges_ = 0;
  // printf("graph name:%s\n", filename);
  // printf("Building bigraph with file \"%s\" ... \n", filename);
  FILE *fp = fopen(filename, "rt");
  if (fp == nullptr) {
    printf("Could not open file %s. Bigraph building failed!\n", filename);
    exit(1);
  }
  std::vector<int> neighbor_list;
  int r_size = 0;
  while (!feof(fp)) {
    neighbor_list = std::move(getLine(fp, r_size));
    edges_ += neighbor_list.size();
    if (!neighbor_list.empty())
      l_adj_lists_.emplace_back(std::move(neighbor_list));
  }
  r_size++;
  r_adj_lists_.resize(r_size);
  for (int i = 0; i < l_adj_lists_.size(); i++) {
    for (int id : l_adj_lists_[i]) {
      r_adj_lists_[id].emplace_back(i);
    }
  }
  l_map_.resize(GetLSize());
  for (int i = 0; i < GetLSize(); i++)
    l_map_[i] = i;
  int r_valid = 0;
  r_map_.resize(r_size);
  for (int i = 0; i < r_size; i++) {
    if (!r_adj_lists_[i].empty()) {
      if (r_valid != i)
        r_adj_lists_[r_valid] = std::move(r_adj_lists_[i]);
      r_map_[r_valid++] = i;
    }
  }
  r_map_.resize(r_valid);
  r_adj_lists_.resize(r_valid);
  if (r_valid != GetRSize()) {
    for (int i = 0; i < GetLSize(); i++)
      l_adj_lists_[i].clear();
    for (int i = 0; i < GetRSize(); i++) {
      for (int id : r_adj_lists_[i])
        l_adj_lists_[id].emplace_back(i);
    }
  }
  // r_size = r_valid;
  l_degree_ = 0;
  for (int i = 0; i < GetLSize(); i++)
    l_degree_ = std::max(l_degree_, (int)l_adj_lists_[i].size());
  r_degree_ = 0;
  for (int i = 0; i < GetRSize(); i++)
    r_degree_ = std::max(r_degree_, (int)r_adj_lists_[i].size());
}

BiGraph::BiGraph(const BiGraph &graph) {
  min_neighbors_of_l_ = graph.min_neighbors_of_l_;
  min_neighbors_of_r_ = graph.min_neighbors_of_r_;
  l_degree_ = graph.l_degree_;
  r_degree_ = graph.r_degree_;
  l_map_ = graph.l_map_;
  r_map_ = graph.r_map_;
  edges_ = graph.edges_;
  for (auto vec : graph.l_adj_lists_) {
    l_adj_lists_.emplace_back(vec);
  }
  for (auto vec : graph.r_adj_lists_) {
    r_adj_lists_.emplace_back(vec);
  }
}

BiGraph::~BiGraph() {
  for (auto &vec : l_adj_lists_)
    vec.clear();
  l_adj_lists_.clear();
  for (auto &vec : r_adj_lists_)
    vec.clear();
  r_adj_lists_.clear();
  l_map_.clear();
  r_map_.clear();
}

void BiGraph::Reorder(OrderEnum op) {
  if (op == Rand) return;
  // std::vector<int> L_rename(GetLSize());
  // std::vector<int> L_rename_rev(GetLSize());
  // std::vector<int> R_rename(GetRSize());
  // for (int i = 0; i < L_rename.size(); i++)
  //   L_rename[i] = i;
  // for (int i = 0; i < R_rename.size(); i++)
  //   R_rename[i] = i;
  // switch (op)
  // {
  // case LInc:
  //     std::sort(
  //         L_rename.begin(), L_rename.end(), [=](int id0, int id1) -> bool {
  //           return l_adj_lists_[id0].size() < l_adj_lists_[id1].size() ||
  //                  (l_adj_lists_[id0].size() == l_adj_lists_[id1].size() &&
  //                   l_adj_lists_[id0] < l_adj_lists_[id1]);
  //         });
  //     break;
  // case LDec:
  //   std::sort(L_rename.begin(), L_rename.end(), [=](int id0, int id1) -> bool
  //   {
  //     return l_adj_lists_[id0].size() > l_adj_lists_[id1].size() ||
  //            (l_adj_lists_[id0].size() == l_adj_lists_[id1].size() &&
  //             l_adj_lists_[id0] < l_adj_lists_[id1]);
  //   });
  //   break;
  // case RInc:
  //   std::sort(R_rename.begin(), R_rename.end(), [=](int id0, int id1) -> bool
  //   {
  //     return r_adj_lists_[id0].size() < r_adj_lists_[id1].size() ||
  //            (r_adj_lists_[id0].size() == r_adj_lists_[id1].size() &&
  //             r_adj_lists_[id0] < r_adj_lists_[id1]);
  //   });
  //   break;
  // case RDec:
  //   std::sort(R_rename.begin(), R_rename.end(), [=](int id0, int id1) -> bool
  //   {
  //     return r_adj_lists_[id0].size() > r_adj_lists_[id1].size() ||
  //            (r_adj_lists_[id0].size() == r_adj_lists_[id1].size() &&
  //             r_adj_lists_[id0] < r_adj_lists_[id1]);
  //   });
  //   break;
  // case UC:
  //   break;

  // case LRInc:
  //   std::sort(L_rename.begin(), L_rename.end(), [=](int id0, int id1) -> bool {
  //     return l_adj_lists_[id0].size() < l_adj_lists_[id1].size() ||
  //            (l_adj_lists_[id0].size() == l_adj_lists_[id1].size() &&
  //             l_adj_lists_[id0] < l_adj_lists_[id1]);
  //   });
  //   std::sort(R_rename.begin(), R_rename.end(), [=](int id0, int id1) -> bool {
  //     return r_adj_lists_[id0].size() < r_adj_lists_[id1].size() ||
  //            (r_adj_lists_[id0].size() == r_adj_lists_[id1].size() &&
  //             r_adj_lists_[id0] < r_adj_lists_[id1]);
  //   });
  //   break;

  // default:
  //   break;
  // }

  // std::vector<std::vector<int>> l_adj_lists_new(GetLSize());
  // std::vector<std::vector<int>> r_adj_lists_new(GetRSize());

  // for (int i = 0; i < L_rename.size(); i++)
  //   L_rename_rev[L_rename[i]] = i;
  // for (int r_id = 0; r_id < GetRSize(); r_id++){
  //   int orig_r_id = R_rename[r_id];
  //   for(int orig_l_id: r_adj_lists_[orig_r_id]){
  //     int l_id = L_rename_rev[orig_l_id];
  //     l_adj_lists_new[l_id].emplace_back(r_id);
  //   }
  // }
  // for (int l_id = 0; l_id < GetLSize(); l_id++){
  //   for(int r_id:l_adj_lists_new[l_id])
  //     r_adj_lists_new[r_id].emplace_back(l_id);
  // }

  // std::swap(l_adj_lists_, l_adj_lists_new);
  // std::swap(r_adj_lists_, r_adj_lists_new);

  // for (auto &list : l_adj_lists_new)
  //   list.clear();
  // l_adj_lists_new.clear();
  // for (auto &list : r_adj_lists_new)
  //   list.clear();
  // r_adj_lists_new.clear();

  int size = GetLSize();
  switch (op) {
  case 2:
  case 3:
  case 4:
    size = GetRSize();
    break;
  default:
    break;
  }

  std::vector<int> reorder_map(size);
  for (int i = 0; i < size; i++)
    reorder_map[i] = i;

  switch (op) {
  case 0:
    std::sort(reorder_map.begin(), reorder_map.end(),
              [=](int id0, int id1) -> bool {
                return l_adj_lists_[id0].size() < l_adj_lists_[id1].size();
              });
    break;
  case 1:
    std::sort(reorder_map.begin(), reorder_map.end(),
              [=](int id0, int id1) -> bool {
                return l_adj_lists_[id0].size() > l_adj_lists_[id1].size();
              });
    break;
  case 2:
    Transpose();
    std::sort(reorder_map.begin(), reorder_map.end(),
              [=](int id0, int id1) -> bool {
                return l_adj_lists_[id0].size() < l_adj_lists_[id1].size() ||
                       (l_adj_lists_[id0].size() == l_adj_lists_[id1].size() &&
                        l_adj_lists_[id0] < l_adj_lists_[id1]);
              });
    break;
  case 3:
    Transpose();
    std::sort(reorder_map.begin(), reorder_map.end(),
              [=](int id0, int id1) -> bool {
                return l_adj_lists_[id0].size() > l_adj_lists_[id1].size() ||
                       (l_adj_lists_[id0].size() == l_adj_lists_[id1].size() &&
                        l_adj_lists_[id0] < l_adj_lists_[id1]);
              });
    break;
  case 4: {

    Transpose();
    std::vector<std::vector<int>> N2(size);
    std::vector<int> NON2(size);
    auto start = get_cur_time();
    /*for(int i = 0; i < size; i++)
    {
      std::vector<int> two_hop_neighbors;
      for(auto neighbor: l_adj_lists_[i])
      {
        two_hop_neighbors = seq_union(two_hop_neighbors,
    r_adj_lists_[neighbor]);
      }
      N2[i] = two_hop_neighbors;
      NON2[i] = two_hop_neighbors.size();
    }*/
    auto computeExa2HopN = [&]() {
      std::vector<int> bfs(size, -1);
      for (int i = 0; i < size; i++) {
        int n2 = 0;
        for (auto &v : l_adj_lists_[i]) { // v is in R
          for (auto &u : r_adj_lists_[v]) {
            if (u != i) {
              if (bfs[u] != i) {
                bfs[u] = i;
                n2++;
              }
            }
          }
        }
        NON2[i] = n2;
      }
    };
    computeExa2HopN();

    int maxDeg2 = 0;
    for (auto &deg2 : NON2) {
      if (deg2 > maxDeg2) {
        maxDeg2 = deg2;
      }
    }

    int *buck = new int[maxDeg2 + 1]();
    for (auto &deg2 : NON2) {
      buck[deg2]++;
    }

    int pos = 0;
    for (int i = 0; i < maxDeg2 + 1; i++) {
      int size = buck[i];
      buck[i] = pos;
      pos = size + pos;
    }
    // buck store the starting position

    // first sort by 2-hop neighbours
    int *sortedV = new int[size];
    int *vLoc = new int[size]; // location of a vertex in a the sorted array.

    for (int i = 0; i < size; i++) {
      sortedV[buck[NON2[i]]] = i;
      vLoc[i] = buck[NON2[i]];
      buck[NON2[i]]++;
    }
    // after the above loop, given a degree buck[degree] is the starting
    // position for vertices with degree+1
    for (int i = maxDeg2; i > 0; i--) {
      buck[i] = buck[i - 1];
    }
    buck[0] = 0;

    std::vector<int> bfs(size, -1);
    for (int i = 0; i < size; i++) {
      int u = sortedV[i];
      reorder_map[i] = u;
      for (auto &v : l_adj_lists_[u]) { // v is in R
        for (auto &u_prime : r_adj_lists_[v]) {
          if (bfs[u_prime] != u) {
            bfs[u_prime] = u;
            if (NON2[u_prime] > NON2[u]) {
              int deg_u_prime = NON2[u_prime];
              int p_u_prime = vLoc[u_prime];
              int positionw = buck[deg_u_prime];
              int w = sortedV[positionw];
              if (u_prime != w) {
                vLoc[u_prime] = positionw;
                sortedV[p_u_prime] = w;
                vLoc[w] = p_u_prime; // the previous u here almost kill me.
                sortedV[positionw] = u_prime;
              }
              buck[deg_u_prime]++;
              NON2[u_prime]--;
            }
          }
        }
      }
    }

    std::cout << "ucorder time: " << get_cur_time() - start << std::endl;
    delete[] vLoc;
    delete[] sortedV;
    delete[] buck;
  } break;
  default:
    break;
  }

  std::vector<int> n_l_map(GetLSize());
  std::vector<std::vector<int>> n_l_adj_list(GetLSize());

  for (int i = 0; i < GetLSize(); i++)
    n_l_map[i] = l_map_[reorder_map[i]];
  std::swap(l_map_, n_l_map);
  n_l_map.clear();
  for (int i = 0; i < GetLSize(); i++)
    n_l_adj_list[i] = std::move(l_adj_lists_[reorder_map[i]]);
  std::swap(l_adj_lists_, n_l_adj_list);
  n_l_adj_list.clear();

  for (int i = 0; i < GetRSize(); i++)
    r_adj_lists_[i].clear();
  for (int i = 0; i < GetLSize(); i++) {
    for (int r_id : l_adj_lists_[i]) {
      r_adj_lists_[r_id].emplace_back(i);
    }
  }
  switch (op) {
  case 2:
  case 3:
  case 4:
    Transpose();
    break;
  default:
    break;
  }
}

void BiGraph::Transpose() {
  std::swap(l_degree_, r_degree_);
  std::swap(min_neighbors_of_l_, min_neighbors_of_r_);
  std::swap(l_adj_lists_, r_adj_lists_);
  std::swap(l_map_, r_map_);
}

int BiGraph::GetLSize() { return l_adj_lists_.size(); }

int BiGraph::GetRSize() { return r_adj_lists_.size(); }

int BiGraph::GetEdges() { return edges_; }

std::vector<int> &BiGraph::NeighborsL(int v) { return l_adj_lists_[v]; }

std::vector<int> BiGraph::NeighborsL(std::vector<int> &v_set) {
  auto res = l_adj_lists_[v_set[0]];
  for (int i = 1; i < v_set.size(); i++)
    seq_intersect_local(res, l_adj_lists_[v_set[i]]);
  // for (int i = 1; i < v_set.size(); i++)
  //   res = seq_intersect(res, l_adj_lists_[v_set[i]]);
  return res;
}

std::vector<int> &BiGraph::NeighborsR(int v) { return r_adj_lists_[v]; }

// bool BiGraph::EdgeExist(int l, int r) {
//   return std::binary_search(l_adj_lists_[l].begin(), l_adj_lists_[l].end(),
//   r);
// }

void BiGraph::CalGraphStatistics() {
  // l_valid_ = r_valid_ = l_degree_ = r_degree_ = edges_ = 0;
  l_degree_ = r_degree_ = edges_ = 0;
  for (int i = 0; i < GetLSize(); i++) {
    if (!l_adj_lists_[i].empty()) {
      // l_valid_++;
      l_degree_ = std::max(l_degree_, (int)l_adj_lists_[i].size());
      edges_ += l_adj_lists_[i].size();
    }
  }
  for (int i = 0; i < GetRSize(); i++) {
    if (!r_adj_lists_[i].empty()) {
      // r_valid_++;
      r_degree_ = std::max(r_degree_, (int)r_adj_lists_[i].size());
    }
  }
}

void BiGraph::PrintProfile() {
  CalGraphStatistics();
  printf("Bigraph Profile:\n");
  printf("Number of left vertices      : %d\n", GetLSize());
  printf("Number of right vertices     : %d\n", GetRSize());
  printf("Maximum left degree          : %d\n", l_degree_);
  printf("Maximum right degree         : %d\n", r_degree_);
  printf("Number of edges              : %llu\n", edges_);
  printf("Min left degree bound        : %d\n", min_neighbors_of_l_);
  printf("Min right degree bound       : %d\n\n", min_neighbors_of_r_);
}

void BiGraph::PrintTotalGraph() {
  for (int i = 0; i < GetLSize(); i++) {
    printf("[%d] ", l_map_[i]);
    for (int id : l_adj_lists_[i])
      printf(" %d", r_map_[id]);
    printf("\n");
  }
  printf("\n");
}

void BiGraph::ConvertIdVector(std::vector<int> &l_vec,
                              std::vector<int> &r_vec) {
  for (int i = 0; i < l_vec.size(); i++)
    l_vec[i] = l_map_[l_vec[i]];
  for (int i = 0; i < r_vec.size(); i++)
    r_vec[i] = r_map_[r_vec[i]];
}

void BiGraph::Prune1H(int min_neighbors_of_l, int min_neighbors_of_r) {
  double start = get_cur_time();
  min_neighbors_of_l_ = std::max(min_neighbors_of_l_, min_neighbors_of_l);
  min_neighbors_of_r_ = std::max(min_neighbors_of_r_, min_neighbors_of_r);
  int orig_edges = edges_;
  edges_ = 0;

  std::vector<int> l_degree(GetLSize());
  std::vector<int> r_degree(GetRSize());
  std::vector<int> l_pruned_seed, r_pruned_seed;
  for (int i = 0; i < GetLSize(); i++) {
    l_degree[i] = l_adj_lists_[i].size();
    if (l_degree[i] < min_neighbors_of_r_)
      l_pruned_seed.emplace_back(i);
  }
  for (int i = 0; i < GetRSize(); i++) {
    r_degree[i] = r_adj_lists_[i].size();
    if (r_degree[i] < min_neighbors_of_l_)
      r_pruned_seed.emplace_back(i);
  }

  while (!l_pruned_seed.empty() || !r_pruned_seed.empty()) {
    for (int l_id : l_pruned_seed) {
      for (int r_id : l_adj_lists_[l_id]) {
        r_degree[r_id]--;
        if (r_degree[r_id] == min_neighbors_of_l - 1)
          r_pruned_seed.emplace_back(r_id);
      }
    }
    l_pruned_seed.clear();
    for (int r_id : r_pruned_seed) {
      for (int l_id : r_adj_lists_[r_id]) {
        l_degree[l_id]--;
        if (l_degree[l_id] == min_neighbors_of_r - 1)
          l_pruned_seed.emplace_back(l_id);
      }
    }
    r_pruned_seed.clear();
  }
  std::vector<int> L_rename(l_adj_lists_.size(), -1);
  std::vector<int> R_rename(r_adj_lists_.size(), -1);
  int l_valid = 0, r_valid = 0;
  for (int i = 0; i < l_adj_lists_.size(); i++)
    if (l_degree[i] >= min_neighbors_of_r_) {
      l_map_[l_valid] = l_map_[i];
      L_rename[i] = l_valid++;
    }
  for (int i = 0; i < r_adj_lists_.size(); i++)
    if (r_degree[i] >= min_neighbors_of_l_) {
      r_map_[r_valid] = r_map_[i];
      R_rename[i] = r_valid++;
    }
  l_map_.resize(l_valid);
  r_map_.resize(r_valid);
  auto l_adj_lists_buf_ = std::move(l_adj_lists_);
  l_adj_lists_.clear();
  l_adj_lists_.resize(l_valid);
  for (auto &list : r_adj_lists_)
    list.clear();
  r_adj_lists_.resize(r_valid);

  for (int i = 0; i < l_adj_lists_buf_.size(); i++) {
    if (l_degree[i] >= min_neighbors_of_r_) {
      int l_id = L_rename[i];
      for (int orig_rid : l_adj_lists_buf_[i]) {
        int r_id = R_rename[orig_rid];
        if (r_id >= 0) {
          l_adj_lists_[l_id].emplace_back(r_id);
          r_adj_lists_[r_id].emplace_back(l_id);
          edges_++;
        }
      }
    }
  }
  printf("[1H]lb %d, rb %d, edges %llu Time:%lfs\n", min_neighbors_of_l_,
         min_neighbors_of_r_, edges_, get_cur_time()-start);
}

void BiGraph::Prune2H(int min_neighbors_of_l, int min_neighbors_of_r) {
  double start = get_cur_time();
  min_neighbors_of_l_ = std::max(min_neighbors_of_l_, min_neighbors_of_l);
  min_neighbors_of_r_ = std::max(min_neighbors_of_r_, min_neighbors_of_r);
  int orig_edges = edges_;
  edges_ = 0;

  std::vector<int> Index(std::max(GetLSize(), GetRSize()), -1);
  std::vector<int> L_rename(l_adj_lists_.size(), -1);
  std::vector<int> R_rename(r_adj_lists_.size(), -1);
  int l_valid = 0;
  int r_valid = 0;

  for (int l = 0; l < GetLSize(); l++) {
    std::vector<std::pair<int, int>> id_cnt_vec;
    for (int r : l_adj_lists_[l]) {
      for (int ll : r_adj_lists_[r]) {
        if (Index[ll] < 0) {
          Index[ll] = id_cnt_vec.size();
          id_cnt_vec.emplace_back(std::make_pair(ll, 1));
        } else {
          id_cnt_vec[Index[ll]].second++;
        }
      }
    }
    int valid_neighbors_cnt = 0;
    for (auto &p : id_cnt_vec) {
      if (p.second >= min_neighbors_of_r)
        valid_neighbors_cnt++;
      Index[p.first] = -1;
    }
    if (valid_neighbors_cnt >= min_neighbors_of_l_){
      l_map_[l_valid] = l_map_[l];
      L_rename[l] = l_valid++;
    } else {
      l_adj_lists_[l].clear();
    }
  }

  for (int r = 0; r < GetRSize(); r++) {
    std::vector<std::pair<int, int>> id_cnt_vec;
    for (int l : r_adj_lists_[r]) {
      for (int rr : l_adj_lists_[l]) {
        if (Index[rr] < 0) {
          Index[rr] = id_cnt_vec.size();
          id_cnt_vec.emplace_back(std::make_pair(rr, 1));
        } else {
          id_cnt_vec[Index[rr]].second++;
        }
      }
    }
    int valid_neighbors_cnt = 0;
    for (auto &p : id_cnt_vec) {
      if (p.second >= min_neighbors_of_l)
        valid_neighbors_cnt++;
      Index[p.first] = -1;
    }
    if (valid_neighbors_cnt >= min_neighbors_of_r_) {
      r_map_[r_valid] = r_map_[r];
      R_rename[r] = r_valid++;
    } else {
      r_adj_lists_[r].clear();
    }
  }

  for (auto &list : r_adj_lists_) list.clear();
  r_adj_lists_.resize(r_valid);
  for (int i = 0; i < l_adj_lists_.size(); i++){
    int l_id = L_rename[i];
    if (l_id >= 0) {
      auto neighbors = std::move(l_adj_lists_[i]);
      l_adj_lists_[i].clear();
      for(int r_orig:neighbors){
        int r_id = R_rename[r_orig];
        if(r_id >= 0){
          l_adj_lists_[l_id].emplace_back(r_id);
          r_adj_lists_[r_id].emplace_back(l_id);
          edges_++;
        }
      }
    }
  }
  l_adj_lists_.resize(l_valid);
  printf("[2H]lb %d, rb %d, edges %llu Time:%lfs\n", min_neighbors_of_l_,
         min_neighbors_of_r_, edges_, get_cur_time() - start);
}

// void BiGraph::Prune2H(int min_neighbors_of_l, int min_neighbors_of_r) {
//   min_neighbors_of_l_buf_ = min_neighbors_of_l_;
//   min_neighbors_of_r_buf_ = min_neighbors_of_r_;
//   edges_buf_ = edges_;
//   min_neighbors_of_l_ = std::max(min_neighbors_of_l,
//   min_neighbors_of_l_buf_); min_neighbors_of_r_ =
//   std::max(min_neighbors_of_r, min_neighbors_of_r_buf_); edges_ = 0;
//   std::vector<int> L_prime, R_prime;
//   for (int i = 0; i < GetLSize(); i++) {
//     if (l_adj_lists_[i].size() < min_neighbors_of_r_)
//       continue;
//     std::unordered_map<int, int> nc_count_map;
//     int _2H_neighbors = 0;
//     for (int r : l_adj_lists_[i]) {
//       for (int new_l : r_adj_lists_[r]) {
//         if (++nc_count_map[new_l] == min_neighbors_of_r_)
//           if (++_2H_neighbors == min_neighbors_of_l_)
//             break;
//         if (_2H_neighbors == min_neighbors_of_l_)
//           break;
//       }
//       if (_2H_neighbors == min_neighbors_of_l_)
//         break;
//     }
//     if (_2H_neighbors == min_neighbors_of_l_)
//       L_prime.emplace_back(i);
//   }

//   for (int i = 0; i < GetRSize(); i++) {
//     if (r_adj_lists_[i].size() < min_neighbors_of_l_)
//       continue;
//     std::unordered_map<int, int> nc_count_map;
//     int _2H_neighbors = 0;
//     auto neighbors = seq_intersect(L_prime, r_adj_lists_[i]);
//     for (int l : neighbors) {
//       for (int new_r : l_adj_lists_[l]) {
//         if (++nc_count_map[new_r] == min_neighbors_of_l_)
//           if (++_2H_neighbors == min_neighbors_of_r_)
//             break;
//         if (_2H_neighbors == min_neighbors_of_r_)
//           break;
//       }
//       if (_2H_neighbors == min_neighbors_of_r_)
//         break;
//     }
//     if (_2H_neighbors == min_neighbors_of_r_)
//       R_prime.emplace_back(i);
//   }
//   l_adj_lists_buf_ = std::move(l_adj_lists_);
//   r_adj_lists_buf_ = std::move(r_adj_lists_);
//   l_adj_lists_ = std::vector<std::vector<int>>(GetLSize(),
//   std::vector<int>(0)); r_adj_lists_ =
//   std::vector<std::vector<int>>(GetRSize(), std::vector<int>(0)); for
//   (int l : L_prime) {
//     l_adj_lists_[l] = seq_intersect(R_prime, l_adj_lists_buf_[l]);
//     edges_ += l_adj_lists_[l].size();
//   }
//   for (int r : R_prime) {
//     r_adj_lists_[r] = seq_intersect(L_prime, r_adj_lists_buf_[r]);
//   }
//   printf("[2H]lb %d, rb %d, edges %llu\n", min_neighbors_of_l_,
//          min_neighbors_of_r_, edges_);
// }

void BiGraph::Prune2HOpt(int min_neighbors_of_l, int min_neighbors_of_r) {
  min_neighbors_of_l_buf_ = min_neighbors_of_l_;
  min_neighbors_of_r_buf_ = min_neighbors_of_r_;
  edges_buf_ = edges_;
  min_neighbors_of_l_ = std::max(min_neighbors_of_l, min_neighbors_of_l_buf_);
  min_neighbors_of_r_ = std::max(min_neighbors_of_r, min_neighbors_of_r_buf_);
  edges_ = 0;

  // pre-processing, remove some low-degree vertices
  std::vector<bool> L_valid(l_adj_lists_.size(), true);
  std::vector<bool> R_valid(r_adj_lists_.size(), true);
  for (int i = 0; i < l_adj_lists_.size(); i++)
    if (l_adj_lists_.size() < min_neighbors_of_r)
      L_valid[i] = false;
  for (int i = 0; i < r_adj_lists_.size(); i++)
    if (r_adj_lists_.size() < min_neighbors_of_l)
      R_valid[i] = false;

  std::vector<int> L_degree(l_adj_lists_.size(), 0);
  std::vector<int> R_degree(r_adj_lists_.size(), 0);
  for (int i = 0; i < l_adj_lists_.size(); i++) {
    int degree = 0;
    if (L_valid[i]) {
      for (int r : l_adj_lists_[i])
        if (R_valid[r])
          degree++;
    }
    if (degree >= min_neighbors_of_r)
      L_degree[i] = degree;
    else
      L_valid[i] = false;
  }
  for (int i = 0; i < r_adj_lists_.size(); i++) {
    int degree = 0;
    if (R_valid[i]) {
      for (int l : r_adj_lists_[i])
        if (L_valid[l])
          degree++;
    }
    if (degree >= min_neighbors_of_l)
      R_degree[i] = degree;
    else
      R_valid[i] = false;
  }

  // processing L stream, using score
  std::vector<int> L_score(l_adj_lists_.size(), 0);
  std::vector<int> L_process_stream;
  for (int i = 0; i < l_adj_lists_.size(); i++) {
    if (L_valid[i]) {
      int score = 0;
      for (int r : l_adj_lists_[i])
        if (R_valid[r])
          score += R_degree[r];
      L_score[i] = score;
      L_process_stream.emplace_back(i);
    }
  }

  std::sort(L_process_stream.begin(), L_process_stream.end(),
            [=](int v0, int v1) -> bool { return L_score[v0] < L_score[v1]; });
  for (int l_processing : L_process_stream) {
    std::unordered_map<int, int> nc_count_map;
    int _2H_neighbors = 0;
    for (int r : l_adj_lists_[l_processing]) {
      if (R_valid[r]) {
        for (int new_l : r_adj_lists_[r]) {
          if (L_valid[new_l]) {
            if (++nc_count_map[new_l] == min_neighbors_of_r_)
              if (++_2H_neighbors == min_neighbors_of_l_)
                break;
            if (_2H_neighbors == min_neighbors_of_l_)
              break;
          }
        }
        if (_2H_neighbors == min_neighbors_of_l_)
          break;
      }
    }
    if (_2H_neighbors < min_neighbors_of_l_)
      L_valid[l_processing] = false;
  }

  // processing R stream, using score
  std::vector<int> R_score(r_adj_lists_.size(), 0);
  std::vector<int> R_process_stream;
  for (int i = 0; i < r_adj_lists_.size(); i++) {
    if (R_valid[i]) {
      int score = 0;
      for (int l : r_adj_lists_[i])
        if (L_valid[l])
          score += L_degree[l];
      R_score[i] = score;
      R_process_stream.emplace_back(i);
    }
  }

  std::sort(R_process_stream.begin(), R_process_stream.end(),
            [=](int v0, int v1) -> bool { return R_score[v0] < R_score[v1]; });
  for (int r_processing : R_process_stream) {
    std::unordered_map<int, int> nc_count_map;
    int _2H_neighbors = 0;
    for (int l : r_adj_lists_[r_processing]) {
      if (L_valid[l]) {
        for (int new_r : l_adj_lists_[l]) {
          if (R_valid[new_r]) {
            if (++nc_count_map[new_r] == min_neighbors_of_l_)
              if (++_2H_neighbors == min_neighbors_of_r_)
                break;
            if (_2H_neighbors == min_neighbors_of_r_)
              break;
          }
        }
        if (_2H_neighbors == min_neighbors_of_r_)
          break;
      }
    }
    if (_2H_neighbors < min_neighbors_of_r_)
      R_valid[r_processing] = false;
  }
  l_adj_lists_buf_ = std::move(l_adj_lists_);
  r_adj_lists_buf_ = std::move(r_adj_lists_);
  l_adj_lists_ = std::vector<std::vector<int>>(GetLSize(), std::vector<int>(0));
  r_adj_lists_ = std::vector<std::vector<int>>(GetRSize(), std::vector<int>(0));

  for (int l = 0; l < l_adj_lists_.size(); l++) {
    if (L_valid[l]) {
      for (int r : l_adj_lists_buf_[l]) {
        if (R_valid[r]) {
          l_adj_lists_[l].emplace_back(r);
          r_adj_lists_[r].emplace_back(l);
          edges_++;
        }
      }
    }
  }

  printf("[2H Opt]lb %d, rb %d, edges %llu\n", min_neighbors_of_l_,
         min_neighbors_of_r_, edges_);
}

void BiGraph::PopOrigGraph() {
  min_neighbors_of_l_ = min_neighbors_of_l_buf_;
  min_neighbors_of_r_ = min_neighbors_of_r_buf_;
  l_adj_lists_ = std::move(l_adj_lists_buf_);
  r_adj_lists_ = std::move(r_adj_lists_buf_);
  edges_ = edges_buf_;
}

void Biclique::Print(FILE *fp) {
  fprintf(fp, "Max edges: %d (%d %d)\n", GetEdges(), (int)left_nodes_.size(),
          (int)right_nodes_.size());
  // std::sort(left_nodes_.begin(), left_nodes_.end());
  // std::sort(right_nodes_.begin(), right_nodes_.end());

  // fprintf(fp, "L nodes[%d]: ", (int)left_nodes_.size());
  // for (int node : left_nodes_) fprintf(fp, "%d ", node);
  // fprintf(fp, "\n");
  // fprintf(fp, "R nodes:[%d]: ", (int)right_nodes_.size());
  // for (int node : right_nodes_) fprintf(fp, "%d ", node);
  // fprintf(fp, "\n\n");
}

Biclique::Biclique(const Biclique &biclique) {
  left_nodes_ = biclique.left_nodes_;
  right_nodes_ = biclique.right_nodes_;
  edges_ = biclique.edges_;
}

Biclique::Biclique() {
  left_nodes_.clear();
  right_nodes_.clear();
  edges_ = 0;
}

void Biclique::CompareAndSet(std::vector<int> &left_in,
                             std::vector<int> &right_in, BiGraph *graph) {
  // Biclique b;
  // b.left_nodes = left_in;
  // b.right_nodes = right_in;
  // graph->ConvertIdVector(b.left_nodes, b.right_nodes);
  // b.Print();

  if (left_in.size() * right_in.size() >
      left_nodes_.size() * right_nodes_.size()) {
    left_nodes_ = left_in;
    right_nodes_ = right_in;
    graph->ConvertIdVector(left_nodes_, right_nodes_);
    edges_ = left_nodes_.size() * right_nodes_.size();
  }
}

void Biclique::CompareAndSet(std::vector<std::pair<int, int>> &left_in,
                             std::vector<int> &right_in, BiGraph *graph) {
  if (left_in.size() * right_in.size() >
      left_nodes_.size() * right_nodes_.size()) {
    left_nodes_.clear();
    for (auto &lp : left_in) left_nodes_.emplace_back(lp.first);
    right_nodes_ = right_in;
    graph->ConvertIdVector(left_nodes_, right_nodes_);
    edges_ = left_nodes_.size() * right_nodes_.size();
  }
}

void Biclique::CompareAndSet(int L_size, int R_size) {
  left_nodes_.clear();
  right_nodes_.clear();
  edges_ = L_size * R_size;
}

void Biclique::Reset() {
  edges_ = 0;
  left_nodes_.clear();
  right_nodes_.clear();
}

int Biclique::GetEdges() { return edges_; }
