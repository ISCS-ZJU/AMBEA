#include "BicliqueFinderFast.h"


MmbeaIntraFinderFast::MmbeaIntraFinderFast(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {}

void MmbeaIntraFinderFast::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);
  std::vector<int> L, R;
  std::vector<int> C;
  for (int i = 0; i < graph_->GetLSize(); i++)
    if (graph_->NeighborsL(i).size() >= min_r_size)
      L.emplace_back(i);
  for (int i = 0; i < graph_->GetRSize(); i++)
    if (graph_->NeighborsR(i).size() >= min_l_size) {
      if (i == 0 || graph_->NeighborsR(i) != graph_->NeighborsR(i - 1))
        C.emplace_back(i);
    }
  for (int i = 0; i < C.size(); i++)
    biclique_find(L, R, C, i);

  finish();
}

void MmbeaIntraFinderFast::biclique_find(const std::vector<int> &L,
                                         const std::vector<int> &R,
                                         const std::vector<int> &C, int vc_id) {
  int vc = C[vc_id];
  std::vector<int> L_prime = seq_intersect(L, graph_->NeighborsR(vc));
  if (L_prime.size() < min_l_size_)
    return;
  std::vector<int> R_prime = graph_->NeighborsL(L_prime);
  processing_nodes_++;

  // maximality check
  auto R_add = seq_except(R_prime, R);
  if (R_add[0] < vc)
    return;
  if (R_prime.size() >= min_r_size_) {
    maximal_nodes_++;
    maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
  }

  // generate C_prime
  std::vector<int> group_array(C.size() - vc_id - 1, 0);
  std::unordered_map<int, int> cand_rev_index;
  auto R_add_iter = R_add.begin();
  for (int i = vc_id + 1; i < C.size(); i++) {
    while (R_add_iter != R_add.end() && *R_add_iter < C[i])
      R_add_iter++;
    if (R_add_iter == R_add.end() || *R_add_iter > C[i])
      cand_rev_index[C[i]] = i - vc_id - 1;
  }
  GroupHelper g_helper;
  SetOpCounterAdd(L_prime.size());
  for (int l : L_prime) {
    for (int v : graph_->NeighborsL(l)) {
      if (cand_rev_index.find(v) != cand_rev_index.end()) {
        int c_id = cand_rev_index[v];
        group_array[c_id] = g_helper.GetCurGroupId(group_array[c_id], l);
      }
    }
  }
  std::vector<int> C_prime;
  for (int i = vc_id + 1; i < C.size(); i++) {
    int cur_id = i - vc_id - 1;
    if (group_array[cur_id] > 0) {
      int nc = g_helper.GetNc(group_array[cur_id]);
      if (nc >= min_l_size_ && g_helper.IsFirst(group_array[cur_id]))
        C_prime.emplace_back(C[i]);
    }
  }

  for (int i = 0; i < C_prime.size(); i++) {
    biclique_find(L_prime, R_prime, C_prime, i);
  }
}

MmbeaFinderFast::MmbeaFinderFast(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {
  cand_rev_vec_.resize(graph_->GetRSize());
  quick_check_bitset_ = std::vector<bool>(graph_->GetRSize(), false);
}

void MmbeaFinderFast::Execute(int min_l_size, int min_r_size) {
  is_transposed_ = min_l_size < min_r_size;
  if (is_transposed_) {
    graph_->Transpose();
    std::swap(min_l_size, min_r_size);
  }
  setup(min_l_size, min_r_size);
  // graph_->Reorder(4);
  graph_->Reorder(RInc);
  std::vector<int> L_all;

  for (int v = 0; v < graph_->GetRSize(); v++) {
    int vp = -1;
    if (graph_->NeighborsR(v).size() < min_l_size_)
      continue;
    if (v > 0 && graph_->NeighborsR(v - 1) == graph_->NeighborsR(v))
      continue;
    processing_nodes_++;

    std::vector<int> L, R;
    std::vector<std::pair<int, int>> C; // c_id, |Nc|
    L = graph_->NeighborsR(v);
    R = graph_->NeighborsL(L);
    std::vector<int> L_current;
    if (R[0] != v) {
      vp = R[0];
      L_current = graph_->NeighborsR(R[0]);
      for (int i = 1; R[i] < v && L_current.size() != L.size(); i++) {
        vp = R[i];
        L_current = seq_intersect(L_current, graph_->NeighborsR(R[i]));
      }
      if (L_current.size() == L.size()) {
        continue;
      }
    }

    std::map<int, int> c_map; // r_id, nc
    for (int w : L) {
      for (int y : graph_->NeighborsL(w)) {
        if (y >= v)
          c_map[y]++;
      }
    }

    for (auto c_node : c_map) {
      if (c_node.second == L.size()) {
        // R.emplace_back(c_node.first);
      } else if (c_node.second >= min_l_size_) {
        C.emplace_back(std::make_pair(c_node.first,
                                      graph_->NeighborsR(c_node.first).size()));
      }
    }

    if (R.size() >= min_r_size_) {
      maximal_nodes_++;
      maximum_biclique_.CompareAndSet(L, R, graph_);
    }

    if (C.size() == 0)
      continue;

    std::vector<std::pair<int, int>> C_prime;
    std::vector<int> cand_ids;

    Partition(L_all, L_current, L, C, C_prime, cand_ids, vp, v);
    // PartitionHash(L, L, L, C, R, C_prime, cand_ids, vp, v);
    // PartitionBitset(L, L, L, C, R, C_prime, cand_ids, vp, v);

    for (int cand_id : cand_ids) {
      biclique_find(L, R, C_prime, v, cand_id);
    }
  }
  if (is_transposed_)
    graph_->Transpose();
  finish();
}

void MmbeaFinderFast::biclique_find(const std::vector<int> &L,
                                    const std::vector<int> &R,
                                    std::vector<std::pair<int, int>> &C, int vp,
                                    int vc_id) {
  int vc = C[vc_id].first;
  std::vector<int> L_prime, R_prime;
  std::vector<std::pair<int, int>> C_prime;

  L_prime = seq_intersect(L, graph_->NeighborsR(vc));

  if (L_prime.size() < min_l_size_)
    return;
  processing_nodes_++;
  R_prime = graph_->NeighborsL(L_prime);

  // maximality check
  std::vector<int> R_add = seq_except(R_prime, R);
  if (R_add[0] < vp)
    return;
  std::vector<int> L_current = L;
  // for (int i = 0, j = 0; R_add[i] < vc;) {
  //   if (R_add[i] == C[j].first) {
  //     L_current = seq_intersect(L_current, graph_->NeighborsR(R_add[i]));
  //     if (L_current.size() == L_prime.size())
  //       return;
  //     i++;
  //     j++;
  //   } else if (R_add[i] < C[j].first)
  //     i++;
  //   else
  //     j++;
  // }

  for (int i = 0; R_add[i] < vc; i++) {
    L_current = seq_intersect(L_current, graph_->NeighborsR(R_add[i]));
    if (L_current.size() == L_prime.size())
      return;
  }

  if (R_prime.size() >= min_r_size_) {
    if (maximal_nodes_ % 1000000 == 0) {
      printf("7  %d %lfs\n", (int)maximal_nodes_, get_cur_time() - start_time_);
    }
    maximal_nodes_++;
    maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
  }

  std::vector<int> cand_ids;
  Partition(L, L_current, L_prime, C, C_prime, cand_ids, vp, vc);
  // PartitionHash(L, L_current, L_prime, C, R_add, C_prime, cand_ids, vp, vc);
  // PartitionBitset(L, L_current, L_prime, C, R_add, C_prime, cand_ids, vp,
  // vc);

  for (int cand_id : cand_ids) {
    biclique_find(L_prime, R_prime, C_prime, vc, cand_id);
  }
}

inline void MmbeaFinderFast::Partition(
    const std::vector<int> &L, const std::vector<int> &L_current,
    std::vector<int> &L_prime, std::vector<std::pair<int, int>> &C,
    std::vector<std::pair<int, int>> &C_prime, std::vector<int> &cand_ids,
    int vp, int vc) {
  SetOpCounterAdd(L_prime.size());
  GroupHelper g_helper;
  std::vector<int> group_array(C.size(), 0);

  for (int l : L_prime) {
    int c_id = 0;
    for (int v : graph_->NeighborsL(l)) {
      while (c_id < C.size() && C[c_id].first < v)
        c_id++;
      if (c_id == C.size())
        break;
      if (C[c_id].first == v) {
        group_array[c_id] = g_helper.GetCurGroupId(group_array[c_id], l);
      }
    }
  }

  for (int i = 0; i < C.size(); i++) {
    int nc = g_helper.GetNc(group_array[i]);
    if (nc >= min_l_size_ && nc < L_prime.size() &&
        g_helper.IsFirst(group_array[i])) {
      if (C[i].first > vc && C[i].second != nc) {
        if (L.size() == L_current.size() ||
            nc != seq_intersect_cnt(L_current, graph_->NeighborsR(C[i].first)))
          cand_ids.emplace_back(C_prime.size());
      }
      C_prime.emplace_back(C[i].first, nc);
    }
    if (nc == C[i].second && C[i].first < vp)
      C[i].first = -1;
  }
}

inline void MmbeaFinderFast::PartitionBitset(
    const std::vector<int> &L, const std::vector<int> &L_current,
    std::vector<int> &L_prime, std::vector<std::pair<int, int>> &C,
    std::vector<int> &R_add, std::vector<std::pair<int, int>> &C_prime,
    std::vector<int> &cand_ids, int vp, int vc) {
  SetOpCounterAdd(L_prime.size());
  GroupHelper g_helper;
  std::vector<int> group_array(C.size(), 0);
  // std::unordered_map<int, int> cand_rev_index;

  auto R_add_iter = R_add.begin();

  for (int i = 0; i < C.size(); i++) {
    if (C[i].first < 0)
      continue;
    while (R_add_iter != R_add.end() && *R_add_iter < C[i].first)
      R_add_iter++;
    if ((R_add_iter == R_add.end() || *R_add_iter > C[i].first)) {
      cand_rev_vec_[C[i].first] = i;
      quick_check_bitset_[C[i].first] = true;
    }
  }

  for (int l : L_prime) {
    for (int v : graph_->NeighborsL(l)) {
      if (quick_check_bitset_[v]) {
        int c_id = cand_rev_vec_[v];
        group_array[c_id] = g_helper.GetCurGroupId(group_array[c_id], l);
      }
    }
  }

  std::vector<int> cand_buf_for_quick_check;

  for (int i = 0; i < C.size(); i++) {
    if (C[i].first < 0)
      continue;
    quick_check_bitset_[C[i].first] = false;
    int nc = g_helper.GetNc(group_array[i]);
    if (nc >= min_l_size_ && nc < L_prime.size() &&
        g_helper.IsFirst(group_array[i])) {
      if (C[i].first > vc && C[i].second != nc) {
        if (L.size() == L_current.size())
          cand_ids.emplace_back(C_prime.size());
        else {
          cand_rev_vec_[C[i].first] = C_prime.size();
          quick_check_bitset_[C[i].first] = true;
          cand_buf_for_quick_check.emplace_back(C[i].first);
        }
      }
      C_prime.emplace_back(C[i].first, nc);
    }
    if (nc == C[i].second && C[i].first < vp)
      C[i].first = -1;
  }

  if (L.size() != L_current.size()) {
    auto L_add = seq_except(L_current, L_prime);
    SetOpCounterAdd(L_add.size());
    for (int l : L_add) {
      for (int v : graph_->NeighborsL(l)) {
        if (quick_check_bitset_[v]) {
          cand_ids.emplace_back(cand_rev_vec_[v]);
          quick_check_bitset_[v] = false;
        }
      }
    }
  }
  for (int v : cand_buf_for_quick_check)
    quick_check_bitset_[v] = false;
}

inline void MmbeaFinderFast::PartitionHash(
    const std::vector<int> &L, const std::vector<int> &L_current,
    std::vector<int> &L_prime, std::vector<std::pair<int, int>> &C,
    std::vector<int> &R_add, std::vector<std::pair<int, int>> &C_prime,
    std::vector<int> &cand_ids, int vp, int vc) {
  SetOpCounterAdd(L_prime.size());
  GroupHelper g_helper;
  std::vector<int> group_array(C.size(), 0);
  std::unordered_map<int, int> cand_rev_index;

  auto R_add_iter = R_add.begin();

  for (int i = 0; i < C.size(); i++) {
    if (C[i].first < 0)
      continue;
    while (R_add_iter != R_add.end() && *R_add_iter < C[i].first)
      R_add_iter++;
    if ((R_add_iter == R_add.end() || *R_add_iter > C[i].first))
      cand_rev_index[C[i].first] = i;
  }

  for (int l : L_prime) {
    for (int v : graph_->NeighborsL(l)) {
      if (cand_rev_index.find(v) != cand_rev_index.end()) {
        int c_id = cand_rev_index[v];
        group_array[c_id] = g_helper.GetCurGroupId(group_array[c_id], l);
      }
    }
  }

  for (int i = 0; i < C.size(); i++) {
    if (C[i].first < 0)
      continue;
    int nc = g_helper.GetNc(group_array[i]);
    if (nc >= min_l_size_ && nc < L_prime.size() &&
        g_helper.IsFirst(group_array[i])) {
      if (C[i].first > vc && C[i].second != nc) {
        if (L.size() == L_current.size())
          cand_ids.emplace_back(C_prime.size());
        else
          cand_rev_index[C[i].first] = -1 - C_prime.size();
      }
      C_prime.emplace_back(C[i].first, nc);
    }
    if (nc == C[i].second && C[i].first < vp)
      C[i].first = -1;
  }

  if (L.size() != L_current.size()) {
    auto L_add = seq_except(L_current, L_prime);
    SetOpCounterAdd(L_add.size());
    for (int l : L_add) {
      for (int v : graph_->NeighborsL(l)) {
        if (cand_rev_index.find(v) != cand_rev_index.end() &&
            cand_rev_index[v] < 0) {
          cand_ids.emplace_back(-cand_rev_index[v] - 1);
          cand_rev_index[v] = 0;
        }
      }
    }
  }
}

ParMmbeaFinderFast::ParMmbeaFinderFast(BiGraph *graph_in, int thread_num,
                                       const char *name)
    : MmbeaFinderFast(graph_in, name), thread_num(thread_num) {}

void ParMmbeaFinderFast::Execute(int min_l_size, int min_r_size) {
  auto mp = tbb::global_control::max_allowed_parallelism;
  tbb::global_control gc(mp, thread_num);
  is_transposed_ = min_l_size < min_r_size;
  if (is_transposed_) {
    graph_->Transpose();
    std::swap(min_l_size, min_r_size);
  }
  setup(min_l_size, min_r_size);
  graph_->Reorder(RInc);
  tbb::parallel_for(0, graph_->GetRSize(), [&](int v) {
    if (graph_->NeighborsR(v).size() < min_l_size_)
      return;
    if (v > 0 && graph_->NeighborsR(v - 1) == graph_->NeighborsR(v))
      return;
    processing_nodes_++;

    std::vector<int> L, R;
    std::vector<std::pair<int, int>> C; // c_id, |Nc|
    L = graph_->NeighborsR(v);

    std::map<int, int> c_map; // r_id, nc
    for (int w : L) {
      for (int y : graph_->NeighborsL(w)) {
        if (y >= v)
          c_map[y]++;
      }
    }
    // tbb::concurrent_hash_map<int, int> c_map;
    // tbb::parallel_for_each (L.begin(), L.end(), [&](int w) {
    //  for (int y : graph_->NeighborsL(w)) {
    //    if (y >= v)
    //    {
    //      tbb::concurrent_hash_map<int, int>::accessor a;
    //      c_map.insert(a, y);
    //      a->second ++;
    //      a.release();
    //    }
    //  }
    //});

    for (auto c_node : c_map) {
      if (c_node.second == L.size()) {
        R.emplace_back(c_node.first);
      } else if (c_node.second >= min_l_size_) {
        C.emplace_back(std::make_pair(c_node.first, c_node.second));
      }
    }

    if (R.size() >= min_r_size_) {
      maximal_nodes_++;
      // maximum_biclique_.CompareAndSet(L, R, graph_);
    }

    if (C.size() == 0)
      return;

    std::vector<int> group_array(C.size(), 0);
    std::vector<std::pair<int, int>> C_prime;
    std::vector<int> cand_ids;

    Partition(L, L, L, C, C_prime, cand_ids, -1, v);
    // PartitionHash(L, L, L, C, R, C_prime, cand_ids, -1, v);
    // PartitionBitset(L, L, L, C, R, C_prime, cand_ids, -1, v);

    tbb::parallel_for(0, (int)C_prime.size(), [&](int i) {
      parallel_biclique_find(L, R, C_prime, v, i);
    });
  });
  if (is_transposed_)
    graph_->Transpose();
  finish();
}
// void ParMmbeaFinderFast::biclique_find(const std::vector<int>& L,
//                                    const std::vector<int>& R,
//                                    std::vector<std::pair<int, int>>& C, int
//                                    vp, int vc_id) {
//  int vc = C[vc_id].first;
//  std::vector<int> L_prime, R_prime;
//  std::vector<std::pair<int, int>> C_prime;
//
//  L_prime = seq_intersect(L, graph_->NeighborsR(vc));
//
//  if (L_prime.size() < min_l_size_) return;
//  processing_nodes_++;
//  R_prime = graph_->NeighborsL(L_prime);
//
//  // maximality check
//  std::vector<int> R_add = seq_except(R_prime, R);
//  if (R_add[0] < vp) return;
//  std::vector<int> L_current = L;
//  for (int i = 0; R_add[i] < vc; i++) {
//    L_current = seq_intersect(L_current, graph_->NeighborsR(R_add[i]));
//    if (L_current.size() == L_prime.size()) return;
//  }
//  if (R_prime.size() >= min_r_size_) {
//    maximal_nodes_++;
//    maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
//  }
//
//  std::vector<int> cand_ids;
//  Partition(L, L_current, L_prime, C, C_prime, cand_ids, vp, vc);
//  // PartitionHash(L, L_current, L_prime, C, R_add, C_prime, cand_ids, vp,
//  vc);
//  // PartitionBitset(L, L_current, L_prime, C, R_add, C_prime, cand_ids, vp,
//  // vc);
//
//  for(int cand_id: cand_ids ){
//    biclique_find(L_prime, R_prime, C_prime, vc, cand_id);
//  };
//}
void ParMmbeaFinderFast::parallel_biclique_find(
    const std::vector<int> &L, const std::vector<int> &R,
    std::vector<std::pair<int, int>> C, int vp, int vc_id) {
  int vc = C[vc_id].first;
  std::vector<int> L_prime, R_prime;
  std::vector<std::pair<int, int>> C_prime;

  L_prime = seq_intersect(L, graph_->NeighborsR(vc));

  if (L_prime.size() < min_l_size_)
    return;
  processing_nodes_++;
  R_prime = graph_->NeighborsL(L_prime);

  // maximality check
  std::vector<int> R_add = seq_except(R_prime, R);
  if (R_add[0] < vp)
    return;
  std::vector<int> L_current = L;
  for (int i = 0; R_add[i] < vc; i++) {
    L_current = seq_intersect(L_current, graph_->NeighborsR(R_add[i]));
    if (L_current.size() == L_prime.size())
      return;
  }
  if (R_prime.size() >= min_r_size_) {
    maximal_nodes_++;
    // maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
  }

  std::vector<int> cand_ids;
  Partition(L, L_current, L_prime, C, C_prime, cand_ids, vp, vc);
  // PartitionHash(L, L_current, L_prime, C, R_add, C_prime, cand_ids, vp, vc);
  // PartitionBitset(L, L_current, L_prime, C, R_add, C_prime, cand_ids, vp,
  // vc);

  tbb::parallel_for_each(cand_ids.begin(), cand_ids.end(), [&](int cand_id) {
    if (L_prime.size() > 100)
      parallel_biclique_find(L_prime, R_prime, C_prime, vc, cand_id);
    else
      biclique_find(L_prime, R_prime, C_prime, vc, cand_id);
  });
}
void ParMmbeaFinderFast::biclique_find(const std::vector<int> &L,
                                       const std::vector<int> &R,
                                       std::vector<std::pair<int, int>> C,
                                       int vp, int vc_id) {
  int vc = C[vc_id].first;
  std::vector<int> L_prime, R_prime;
  std::vector<std::pair<int, int>> C_prime;

  L_prime = seq_intersect(L, graph_->NeighborsR(vc));

  if (L_prime.size() < min_l_size_)
    return;
  processing_nodes_++;
  R_prime = graph_->NeighborsL(L_prime);

  // maximality check
  std::vector<int> R_add = seq_except(R_prime, R);
  if (R_add[0] < vp)
    return;
  std::vector<int> L_current = L;
  for (int i = 0; R_add[i] < vc; i++) {
    L_current = seq_intersect(L_current, graph_->NeighborsR(R_add[i]));
    if (L_current.size() == L_prime.size())
      return;
  }
  if (R_prime.size() >= min_r_size_) {
    maximal_nodes_++;
    // maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
  }

  std::vector<int> cand_ids;
  Partition(L, L_current, L_prime, C, C_prime, cand_ids, vp, vc);
  // PartitionHash(L, L_current, L_prime, C, R_add, C_prime, cand_ids, vp, vc);
  // PartitionBitset(L, L_current, L_prime, C, R_add, C_prime, cand_ids, vp,
  // vc);

  for (int cand_id : cand_ids) {
    biclique_find(L_prime, R_prime, C_prime, vc, cand_id);
  }
}

MbeaFinderBitset::MbeaFinderBitset(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {}

void MbeaFinderBitset::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);
  my_clock_ = 0;
  graph_->Reorder(RInc);

  for (int v = 0; v < graph_->GetRSize(); v++) {
    if (v > 0 && graph_->NeighborsR(v - 1) == graph_->NeighborsR(v))
      continue;
    processing_nodes_++;
    maximal_nodes_++;
    std::vector<int> &L = graph_->NeighborsR(v);
    if (L.size() <= 32) {
      double start = get_cur_time();
      std::map<int, BITSET_T> c_bs_map;
      for (int i = 0; i < L.size(); i++) {
        int l = L[i];
        for (int w : graph_->NeighborsL(l)) {
          c_bs_map[w] |= (1 << i);
        }
      }
      std::vector<BITSET_T> c_bs_array;
      int vvc = -1;
      BITSET_T l_bs = 0xffffffff >> (32 - L.size());
      for (auto &p : c_bs_map) {
        if (p.second != l_bs) {
          if (vvc < 0 && p.first > v)
            vvc = c_bs_array.size();
          c_bs_array.emplace_back(p.second);
        }
      }
      for (int i = vvc; i < c_bs_array.size(); i++) {
        biclique_find(~0, c_bs_array[i], i, c_bs_array);
      }
      my_clock_ += get_cur_time() - start;
    } else {
      std::vector<int> R, C;
      std::map<int, int> c_nc_map;
      R.emplace_back(v);
      for (int l : L) {
        for (int r : graph_->NeighborsL(l)) {
          if (r > v)
            c_nc_map[r]++;
        }
      }
      for (auto &p : c_nc_map) {
        if (p.second == L.size())
          R.emplace_back(p.first);
        else
          C.emplace_back(p.first);
      }
      for (int i = 0; i < C.size(); i++) {
        biclique_find(L, R, C, i);
      }
    }
  }
  printf("bitset cal time:%lfs\n", my_clock_);

  finish();
}

void MbeaFinderBitset::biclique_find(const std::vector<int> &L,
                                     const std::vector<int> &R,
                                     const std::vector<int> &C, int vc) {
  std::vector<int> L_prime = seq_intersect(L, graph_->NeighborsR(C[vc]));
  processing_nodes_++;
  std::vector<int> R_prime = graph_->NeighborsL(L_prime);
  std::vector<int> R_add = seq_except(R_prime, R);
  if (!R_add.empty() && R_add[0] == C[vc]) {
    maximal_nodes_++;
    if (L_prime.size() <= 32) {
      double start = get_cur_time();
      std::map<int, BITSET_T> c_bs_map;
      for (int i = 0; i < L_prime.size(); i++) {
        int l = L_prime[i];
        for (int w : graph_->NeighborsL(l)) {
          c_bs_map[w] |= (1 << i);
        }
      }
      std::vector<BITSET_T> c_bs_array;
      int vvc = -1;
      BITSET_T l_bs = 0xffffffff >> (32 - L_prime.size());
      for (auto &p : c_bs_map) {
        if (p.second != l_bs) {
          if (vvc < 0 && p.first > C[vc])
            vvc = c_bs_array.size();
          c_bs_array.emplace_back(p.second);
        }
      }
      for (int i = vvc; i < c_bs_array.size(); i++) {
        biclique_find(~0, c_bs_array[i], i, c_bs_array);
      }
      my_clock_ += get_cur_time() - start;
    } else {
      std::vector<int> C_prime;
      for (int i = vc; i < C.size(); i++) {
        int nc = seq_intersect_cnt(L_prime, graph_->NeighborsR(C[i]));
        if (nc != 0 && nc != L_prime.size())
          C_prime.emplace_back(C[i]);
      }
      for (int i = 0; i < C_prime.size(); i++) {
        biclique_find(L_prime, R_prime, C_prime, i);
      }
    }
  }
}
void MbeaFinderBitset::biclique_find(BITSET_T parent_L, BITSET_T current_L,
                                     int vvc, std::vector<BITSET_T> &c_bs) {
  processing_nodes_++;
  // maximality check
  for (int i = 0; i < vvc; i++) {
    BITSET_T nc = current_L & c_bs[i];
    if (nc == current_L) {
      BITSET_T pnc = parent_L & c_bs[i];
      if (pnc != parent_L)
        return;
    }
  }

  maximal_nodes_++;
  for (int i = vvc; i < c_bs.size(); i++) {
    BITSET_T nc = current_L & c_bs[i];
    if (nc != current_L && nc != 0)
      biclique_find(current_L, nc, i, c_bs);
  }
}
MbeaFinderAuto::MbeaFinderAuto(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name), my_clock_(0) {}

void MbeaFinderAuto::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);
  graph_->Reorder(RInc);
  max_c_size_ = 0;

  for (int v = 0; v < graph_->GetRSize(); v++) {
    if (graph_->NeighborsR(v).size() < min_l_size_)
      continue;
    if (v > 0 && graph_->NeighborsR(v - 1) == graph_->NeighborsR(v))
      continue;
    processing_nodes_++;
    biclique_find(graph_->NeighborsR(v), v, 0);
  }
  finish();
  printf("max c size:%d \n", max_c_size_);
  printf("time of bitset:%lfs \n", my_clock_);
}

void MbeaFinderAuto::biclique_find(const std::vector<int> &L, int vc,
                                   int vc_id) {
  std::vector<int> L_prime = seq_intersect(L, graph_->NeighborsR(vc));
  if (L_prime.size() < min_l_size_)
    return;
  processing_nodes_++;

  if (L_prime.size() > 32) {
    std::map<int, int> c_nc_size_map;
    for (int l : L_prime) {
      for (int r : graph_->NeighborsL(l)) {
        c_nc_size_map[r]++;
      }
    }
    max_c_size_ = std::max(max_c_size_, (int)c_nc_size_map.size());
    int R_prime_size = 0;
    for (auto &p : c_nc_size_map) {
      if (p.second == L_prime.size()) {
        if (p.first == vc && R_prime_size != vc_id)
          return; // non-maximal
        R_prime_size++;
      } else if (p.first > vc) {
        biclique_find(L_prime, p.first, R_prime_size);
      }
    }
    if (R_prime_size >= min_r_size_)
      maximal_nodes_++;
  } else {
    double start = get_cur_time();
    std::map<int, BITSET_T> c_bs_map;
    for (int i = 0; i < L_prime.size(); i++) {
      int l = L_prime[i];
      for (int w : graph_->NeighborsL(l)) {
        c_bs_map[w] |= (1 << i);
      }
    }
    std::vector<BITSET_T> c_bs;
    int vvc = -1;
    BITSET_T l_bs = 0xffffffff >> (32 - L_prime.size());
    int R_prime_size = 0;
    for (auto &p : c_bs_map) {
      if (p.second == l_bs) {
        if (p.first == vc && R_prime_size != vc_id) {
          my_clock_ += get_cur_time() - start;
          return; // non-maximal
        }
        R_prime_size++;
      } else {
        if (vvc < 0 && p.first > vc)
          vvc = c_bs.size();
        c_bs.emplace_back(p.second);
      }
    }
    if (vvc >= 0)
      for (int i = vvc; i < c_bs.size(); i++) {
        biclique_find(~0, c_bs[i], i, c_bs);
      }
    my_clock_ += get_cur_time() - start;
    maximal_nodes_++;
  }
}

void MbeaFinderAuto::biclique_find(BITSET_T parent_L, BITSET_T current_L,
                                   int vvc, std::vector<BITSET_T> &c_bs) {
  processing_nodes_++;
  // maximality check
  for (int i = 0; i < vvc; i++) {
    BITSET_T nc = current_L & c_bs[i];
    if (nc == current_L) {
      BITSET_T pnc = parent_L & c_bs[i];
      if (pnc != parent_L)
        return;
    }
  }

  maximal_nodes_++;
  for (int i = vvc; i < c_bs.size(); i++) {
    BITSET_T nc = current_L & c_bs[i];
    if (nc != current_L && nc != 0)
      biclique_find(current_L, nc, i, c_bs);
  }
}
ParMbeaFinderAuto::ParMbeaFinderAuto(BiGraph *graph_in, int thread_num,
                                     const char *name)
    : MbeaFinderAuto(graph_in, name), thread_num(thread_num) {}
void ParMbeaFinderAuto::Execute(int min_l_size, int min_r_size) {
  auto mp = tbb::global_control::max_allowed_parallelism;
  tbb::global_control gc(mp, thread_num);
  setup(min_l_size, min_r_size);
  graph_->Reorder(RInc);
  max_c_size_ = 0;

  tbb::parallel_for(0, graph_->GetRSize(), [&](int v) {
    if (graph_->NeighborsR(v).size() < min_l_size_)
      return;
    if (v > 0 && graph_->NeighborsR(v - 1) == graph_->NeighborsR(v))
      return;
    processing_nodes_++;
    parallel_biclique_find(graph_->NeighborsR(v), v, 0);
  });
  finish();
  printf("max c size:%d \n", max_c_size_);
  printf("time of bitset:%lfs \n", my_clock_);
}
void ParMbeaFinderAuto::parallel_biclique_find(BITSET_T parent_L,
                                               BITSET_T current_L, int vvc,
                                               std::vector<BITSET_T> &c_bs) {
  processing_nodes_++;
  // maximality check
  for (int i = 0; i < vvc; i++) {
    BITSET_T nc = current_L & c_bs[i];
    if (nc == current_L) {
      BITSET_T pnc = parent_L & c_bs[i];
      if (pnc != parent_L)
        return;
    }
  }

  maximal_nodes_++;
  tbb::parallel_for(vvc, (int)c_bs.size(), [&](int i) {
    BITSET_T nc = current_L & c_bs[i];
    if (nc != current_L && nc != 0) {
      if (c_bs.size() < 1000)
        biclique_find(current_L, nc, i, c_bs);
      else
        parallel_biclique_find(current_L, nc, i, c_bs);
    }
  });
}

void ParMbeaFinderAuto::parallel_biclique_find(const std::vector<int> &L,
                                               int vc, int vc_id) {
  std::vector<int> L_prime = seq_intersect(L, graph_->NeighborsR(vc));
  if (L_prime.size() < min_l_size_)
    return;
  processing_nodes_++;

  if (L_prime.size() > 32) {
    std::map<int, int> c_nc_size_map;
    std::map<int, int> new_vc_id_map;
    for (int l : L_prime) {
      for (int r : graph_->NeighborsL(l)) {
        c_nc_size_map[r]++;
      }
    }
    max_c_size_ = std::max(max_c_size_, (int)c_nc_size_map.size());
    int R_prime_size = 0;
    for (auto &p : c_nc_size_map) {
      if (p.second == L_prime.size()) {
        if (p.first == vc && R_prime_size != vc_id)
          return; // non-maximal
        R_prime_size++;
      }
      new_vc_id_map[p.first] = R_prime_size;
    }
    tbb::parallel_for_each(
        c_nc_size_map.begin(), c_nc_size_map.end(),
        [&](std::pair<int, int> tp) {
          if (tp.second == L_prime.size() || tp.first <= vc)
            return;
          if (L_prime.size() > 1000) {
            parallel_biclique_find(L_prime, tp.first, new_vc_id_map[tp.first]);
          } else {
            biclique_find(L_prime, tp.first, new_vc_id_map[tp.first]);
          }
        });
    if (R_prime_size >= min_r_size_)
      maximal_nodes_++;
  } else {
    double start = get_cur_time();
    std::map<int, BITSET_T> c_bs_map;
    for (int i = 0; i < L_prime.size(); i++) {
      int l = L_prime[i];
      for (int w : graph_->NeighborsL(l)) {
        c_bs_map[w] |= (1 << i);
      }
    }
    std::vector<BITSET_T> c_bs;
    int vvc = -1;
    BITSET_T l_bs = 0xffffffff >> (32 - L_prime.size());
    int R_prime_size = 0;
    for (auto &p : c_bs_map) {
      if (p.second == l_bs) {
        if (p.first == vc && R_prime_size != vc_id) {
          my_clock_ += get_cur_time() - start;
          return; // non-maximal
        }
        R_prime_size++;
      } else {
        if (vvc < 0 && p.first > vc)
          vvc = c_bs.size();
        c_bs.emplace_back(p.second);
      }
    }
    if (vvc >= 0)
      for (int i = vvc; i < c_bs.size(); i++) {
        biclique_find(~0, c_bs[i], i, c_bs);
      }
    my_clock_ += get_cur_time() - start;
    maximal_nodes_++;
  }
}
