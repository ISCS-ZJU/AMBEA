#include "AggressiveFinder.h"

AmbeaFinderNaive::AmbeaFinderNaive(BiGraph *graph, const char *name)
    : BicliqueFinder(graph, name) {}

void AmbeaFinderNaive::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);
  graph_->Reorder(RInc);

  for (int v = 0; v < graph_->GetRSize(); v++) {
    if (graph_->NeighborsR(v).size() < min_l_size_)
      continue;
    if (v > 0 && graph_->NeighborsR(v - 1) == graph_->NeighborsR(v))
      continue;
    processing_nodes_++;
    VertexSet &L = graph_->NeighborsR(v);
    VertexSet R, C;
    std::map<int, int> c_map;

    for (int l : L) {
      for (int r : graph_->NeighborsL(l)) {
        if (r > v)
          c_map[r]++;
      }
    }
    R.emplace_back(v);
    for (auto c_node : c_map) {
      if (c_node.second == L.size())
        R.emplace_back(c_node.first);
      else if (c_node.second >= min_l_size_)
        C.emplace_back(c_node.first);
    }
    if (R.size() >= min_r_size_) {
      maximal_nodes_++;
      maximum_biclique_.CompareAndSet(L, R, graph_);
#ifdef COMPUTE_LEVEL
      level_accumulation_ += cur_level_;
      max_level_ = std::max(cur_level_, max_level_);
#endif    
    }
    for (int vc : C) {
#ifdef COMPUTE_LEVEL
      cur_level_++;
#endif
      biclique_find(L, R, C, vc);
#ifdef COMPUTE_LEVEL
      cur_level_--;
#endif      
    }
  }
  finish();
}
void AmbeaFinderNaive::biclique_find(const VertexSet &L, const VertexSet &R,
                                     const VertexSet &C, int v) {
  VertexSet L_prime, R_prime, C_prime, L_remove;
  processing_nodes_++;
  seq_intersect_diff(L, graph_->NeighborsR(v), L_prime, L_remove);
  R_prime = graph_->NeighborsL(L_prime);
  int idx = std::lower_bound(R.begin(), R.end(), C[0]) - R.begin();
  if (R_prime[idx] < C[0])
    return;

  for (int vc : C) {
    int ln = seq_intersect_cnt(L_prime, graph_->NeighborsR(vc));
    if (ln == L_prime.size()) {
      if (vc < v) {
        seq_intersect_local(L_remove, graph_->NeighborsR(vc));
        if (L_remove.empty())
          return;
      }
    } else if (ln >= min_l_size_ && vc > v)
      C_prime.emplace_back(vc);
  }

  if (R_prime.size() >= min_r_size_) {
    maximal_nodes_++;
    maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
#ifdef COMPUTE_LEVEL
    level_accumulation_ += cur_level_;
    max_level_ = std::max(cur_level_, max_level_);
#endif  
  }
  for (int vc : C_prime) {
    if (seq_intersect_cnt(graph_->NeighborsR(vc), L_remove) > 0){
#ifdef COMPUTE_LEVEL
      cur_level_++;
#endif
      biclique_find(L_prime, R_prime, C_prime, vc);
#ifdef COMPUTE_LEVEL
      cur_level_--;
#endif
    } else
      processing_nodes_++;
  }
}

AmbeaFinderInter::AmbeaFinderInter(BiGraph *graph, const char *name)
    : BicliqueFinder(graph, name) {}

void AmbeaFinderInter::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);
  graph_->Reorder(RInc);

  for (int v = 0; v < graph_->GetRSize(); v++) {
    if (graph_->NeighborsR(v).size() < min_l_size_)
      continue;
    if (v > 0 && graph_->NeighborsR(v - 1) == graph_->NeighborsR(v))
      continue;
    processing_nodes_++;
    VertexSet &L = graph_->NeighborsR(v);
    VertexSet R;
    std::vector<std::pair<int, int>> C;
    std::map<int, int> c_map;

    for (int l : L) {
      for (int r : graph_->NeighborsL(l)) {
        if (r > v)
          c_map[r]++;
      }
    }
    R.emplace_back(v);
    for (auto c_node : c_map) {
      if (c_node.second == L.size())
        R.emplace_back(c_node.first);
      else if (c_node.second >= min_l_size_)
        C.emplace_back(c_node);
    }
    if (R.size() >= min_r_size_) {
      maximal_nodes_++;
      maximum_biclique_.CompareAndSet(L, R, graph_);
#ifdef COMPUTE_LEVEL
      level_accumulation_ += cur_level_;
      max_level_ = std::max(cur_level_, max_level_);
#endif
    }
    for (auto c : C) {
#ifdef COMPUTE_LEVEL
      cur_level_++;
#endif
      biclique_find(L, R, C, c.first);
#ifdef COMPUTE_LEVEL
      cur_level_--;
#endif
    }
  }
  finish();
}

void AmbeaFinderInter::biclique_find(const VertexSet &L, const VertexSet &R,
                                     const std::vector<std::pair<int, int>> &C,
                                     int v) {
  VertexSet L_prime, R_prime, L_remove, C_cand;
  std::vector<std::pair<int, int>> C_prime;
  processing_nodes_++;
  seq_intersect_diff(L, graph_->NeighborsR(v), L_prime, L_remove);
  R_prime = graph_->NeighborsL(L_prime);

  int idx = std::lower_bound(R.begin(), R.end(), C[0].first) - R.begin();
  if (R_prime[idx] < C[0].first)
    return;

  for (auto c : C) {
    int ln = seq_intersect_cnt(L_prime, graph_->NeighborsR(c.first));
    if (ln == L_prime.size()) {
      if (c.first < v) {
        seq_intersect_local(L_remove, graph_->NeighborsR(c.first));
        if (L_remove.empty())
          return;
      }
    } else if (ln >= min_l_size_ && c.first > v) {
      if (ln != c.second)
        C_cand.emplace_back(c.first);
      C_prime.emplace_back(std::make_pair(c.first, ln));
    }
  }
  if (R_prime.size() >= min_r_size_) {
    maximal_nodes_++;
    maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
#ifdef COMPUTE_LEVEL
    level_accumulation_ += cur_level_;
    max_level_ = std::max(cur_level_, max_level_);
#endif
  }
  for (int vc : C_cand) {
    if (L_prime.size() + L_remove.size() == L.size() ||
        seq_intersect_cnt(graph_->NeighborsR(vc), L_remove) > 0){
#ifdef COMPUTE_LEVEL
      cur_level_++;
#endif      
      biclique_find(L_prime, R_prime, C_prime, vc);
#ifdef COMPUTE_LEVEL
      cur_level_--;
#endif    
    } else
      processing_nodes_++;
  }
}

VertexSet IntraFinder::Partition(const VertexSet &C, const VertexSet &L_prime,
                                 int v) {
  GroupHelper g_helper;
  std::vector<int> group_array(C.size(), 0);
  VertexSet C_prime;

  for (int l : L_prime) {
    int c_id = 0;
    for (int v : graph_->NeighborsL(l)) {
      while (c_id < C.size() && C[c_id] < v)
        c_id++;
      if (c_id == C.size())
        break;
      if (C[c_id] == v) {
        group_array[c_id] = g_helper.GetCurGroupId(group_array[c_id], l);
      }
    }
  }

  for (int i = 0; i < C.size(); i++) {
    int nc = g_helper.GetNc(group_array[i]);
    if (nc >= min_l_size_ && nc < L_prime.size() &&
        g_helper.IsFirst(group_array[i]) && C[i] > v) {
      C_prime.emplace_back(C[i]);
    }
  }

  return C_prime;
}

IntraFinder::IntraFinder(BiGraph *graph, bool mode, const char *name)
    : BicliqueFinder(graph, name), active_node_(mode) {}

void IntraFinder::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);
  graph_->Reorder(RInc);
  for (int v = 0; v < graph_->GetRSize(); v++) {
    if (graph_->NeighborsR(v).size() < min_l_size_)
      continue;
    if (v > 0 && graph_->NeighborsR(v - 1) == graph_->NeighborsR(v))
      continue;
    processing_nodes_++;
    VertexSet &L = graph_->NeighborsR(v);
    VertexSet R, C;
    std::map<int, int> c_map;

    for (int l : L) {
      for (int r : graph_->NeighborsL(l)) {
        if (r > v)
          c_map[r]++;
      }
    }
    R.emplace_back(v);
    for (auto c_node : c_map) {
      if (c_node.second == L.size())
        R.emplace_back(c_node.first);
      else if (c_node.second >= min_l_size_)
        C.emplace_back(c_node.first);
    }
    if (R.size() >= min_r_size_) {
      maximal_nodes_++;
      maximum_biclique_.CompareAndSet(L, R, graph_);
#ifdef COMPUTE_LEVEL
      level_accumulation_ += cur_level_;
      max_level_ = std::max(cur_level_, max_level_);
#endif    
    }
#ifdef COMPUTE_LEVEL
    cur_level_++;
#endif
    if (active_node_) {
      C = Partition(C, L);
      biclique_find(L, R, C);
    } else
      biclique_find_passive(L, R, C);
#ifdef COMPUTE_LEVEL
    cur_level_--;
#endif
  }
  finish();
}

void IntraFinder::biclique_find(const VertexSet &L, const VertexSet &R,
                                const VertexSet &C) {
  for (int i = 0; i < C.size(); i++) {
    processing_nodes_++;
    VertexSet L_prime, R_prime, C_prime;
    L_prime = seq_intersect(L, graph_->NeighborsR(C[i]));
    R_prime = graph_->NeighborsL(L_prime);
    VertexSet R_add = seq_except(R_prime, R);
    if (R_add[0] == C[i]) {
      if (R_prime.size() >= min_r_size_) {
        maximal_nodes_++;
        maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
#ifdef COMPUTE_LEVEL
        level_accumulation_ += cur_level_;
        max_level_ = std::max(cur_level_, max_level_);
#endif
      }
      C_prime = Partition(C, L_prime, C[i]);
#ifdef COMPUTE_LEVEL
      cur_level_++;
#endif
      biclique_find(L_prime, R_prime, C_prime);
#ifdef COMPUTE_LEVEL
      cur_level_--;
#endif
    }
  }
}

void IntraFinder::biclique_find_passive(const VertexSet &L, const VertexSet &R,
                                        VertexSet &C) {
  for (int i = 0; i < C.size(); i++) {
    if (C[i] < 0)
      continue;
    processing_nodes_++;
    VertexSet L_prime, R_prime, C_prime, L_remove;
    seq_intersect_diff(L, graph_->NeighborsR(C[i]), L_prime, L_remove);
    R_prime = graph_->NeighborsL(L_prime);
    VertexSet R_add = seq_except(R_prime, R);
    if (R_add[0] == C[i]) {
      if (R_prime.size() >= min_r_size_) {
        maximal_nodes_++;
        maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
#ifdef COMPUTE_LEVEL
        level_accumulation_ += cur_level_;
        max_level_ = std::max(cur_level_, max_level_);
#endif
      }
      for (int j = i + 1; j < C.size(); j++) {
        if (C[j] < 0)
          continue;
        int ln = seq_intersect_cnt(L_prime, graph_->NeighborsR(C[j]));
        if (ln == L_prime.size() &&
            seq_intersect_cnt(L_remove, graph_->NeighborsR(C[j])) == 0)
          C[j] = -1;
        else if (ln >= min_l_size_ && ln < L_prime.size())
          C_prime.emplace_back(C[j]);
      }
#ifdef COMPUTE_LEVEL
      cur_level_++;
#endif
      biclique_find(L_prime, R_prime, C_prime);
#ifdef COMPUTE_LEVEL
      cur_level_--;
#endif
    }
  }
}

AmbeaFinder::AmbeaFinder(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {
  order_ = RInc;
  cand_rev_vec_.resize(graph_->GetRSize());
  quick_check_bitset_ = std::vector<bool>(graph_->GetRSize(), false);
}

void AmbeaFinder::Execute(int min_l_size, int min_r_size) {
  // is_transposed_ = min_l_size < min_r_size;
  // if (is_transposed_) {
  //   graph_->Transpose();
  //   std::swap(min_l_size, min_r_size);
  // }
  setup(min_l_size, min_r_size);
  graph_->Reorder(order_);
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
#ifdef COMPUTE_LEVEL
      level_accumulation_ += cur_level_;
      max_level_ = std::max(cur_level_, max_level_);
#endif
    }

    if (C.size() == 0)
      continue;

    std::vector<std::pair<int, int>> C_prime;
    std::vector<int> cand_ids;

    Partition(L_all, L_current, L, C, C_prime, cand_ids, vp, v);

    for (int cand_id : cand_ids) {
#ifdef COMPUTE_LEVEL
      cur_level_++;
#endif
      biclique_find(L, R, C_prime, v, cand_id);
#ifdef COMPUTE_LEVEL
      cur_level_--;
#endif
    }
  }
  // if (is_transposed_)
  //   graph_->Transpose();
  finish();
}

void AmbeaFinder::SetOrder(OrderEnum order) { 
  order_ = order; 
  switch (order)
  {
  case RInc:
    break;
  case UC:
    strcpy(finder_name_, "AmbeaUC");
    break;
  case Rand:
    strcpy(finder_name_, "AmbeaRand");
    break;
  default:
    strcpy(finder_name_, "AmbeaError");
    break;
  }  
}

void AmbeaFinder::biclique_find(const std::vector<int> &L,
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

  for (int i = 0; R_add[i] < vc; i++) {
    L_current = seq_intersect(L_current, graph_->NeighborsR(R_add[i]));
    if (L_current.size() == L_prime.size())
      return;
  }

  if (R_prime.size() >= min_r_size_) {
    maximal_nodes_++;
    maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
#ifdef COMPUTE_LEVEL
    level_accumulation_ += cur_level_;
    max_level_ = std::max(cur_level_, max_level_);
#endif
  }

  std::vector<int> cand_ids;
  Partition(L, L_current, L_prime, C, C_prime, cand_ids, vp, vc);

  for (int cand_id : cand_ids) {
#ifdef COMPUTE_LEVEL
    cur_level_++;
#endif
    biclique_find(L_prime, R_prime, C_prime, vc, cand_id);
#ifdef COMPUTE_LEVEL
    cur_level_--;
#endif
  }
}

inline void AmbeaFinder::Partition(
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
      // if(C[i].first > vc)
      C_prime.emplace_back(C[i].first, nc);
    }
    if (nc == C[i].second && C[i].first < vp)
      C[i].first = -1;
  }
}

inline void AmbeaFinder::PartitionBitset(
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

inline void AmbeaFinder::PartitionHash(
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

ParAmbeaFinder::ParAmbeaFinder(BiGraph *graph_in, int thread_num,
                                       const char *name)
    : AmbeaFinder(graph_in, name), thread_num(thread_num) {}

void ParAmbeaFinder::Execute(int min_l_size, int min_r_size) {
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

void ParAmbeaFinder::parallel_biclique_find(
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

  tbb::parallel_for_each(cand_ids.begin(), cand_ids.end(), [&](int cand_id) {
    if (L_prime.size() > 100)
      parallel_biclique_find(L_prime, R_prime, C_prime, vc, cand_id);
    else
      biclique_find(L_prime, R_prime, C_prime, vc, cand_id);
  });
}
void ParAmbeaFinder::biclique_find(const std::vector<int> &L,
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

  for (int cand_id : cand_ids) {
    biclique_find(L_prime, R_prime, C_prime, vc, cand_id);
  }
}