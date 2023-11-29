#include "BaselineFinder.h"

MbeaFinder::MbeaFinder(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {}

void MbeaFinder::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);

  for (int v = 0; v < graph_->GetRSize(); v++) {
    if (graph_->NeighborsR(v).size() >= min_l_size_) {
      processing_nodes_++;
      std::vector<int> L, R, P, Q;
      L = graph_->NeighborsR(v);

      bool is_maximal = true;

      std::map<int, int> c_map;

      for (int w : L) {
        for (int y : graph_->NeighborsL(w)) {
          if (++c_map[y] == L.size() && y < v) {
            is_maximal = false;
            break;
          }
        }
        if (!is_maximal)
          break;
      }
      if (!is_maximal)
        continue;
      for (auto c_node : c_map) {
        if (c_node.second == L.size())
          R.emplace_back(c_node.first);
        else if (c_node.second >= min_l_size_) {
          if (c_node.first < v)
            Q.emplace_back(c_node.first);
          else
            P.emplace_back(c_node.first);
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
      for (int i = 0; i < P.size() / 2; i++) {
        std::swap(P[i], P[P.size() - 1 - i]);
      }
#ifdef COMPUTE_LEVEL
      cur_level_++;
#endif
      biclique_find(L, R, P, Q);
#ifdef COMPUTE_LEVEL
      cur_level_--;
#endif
    }
  }
  finish();
}

void MbeaFinder::biclique_find(std::vector<int> L, std::vector<int> R,
                                  std::vector<int> P, std::vector<int> Q) {
  while (!P.empty()) {
    std::vector<int> L_prime, R_prime, P_prime, Q_prime;
    bool is_maximal = true;

    int x = P.back();
    P.pop_back();
    L_prime = std::move(seq_intersect(L, graph_->NeighborsR(x)));
    R_prime = R;
    R_prime.emplace_back(x);

    if (L_prime.size() < min_l_size_) {
      Q.emplace_back(x);
      continue;
    }
    processing_nodes_++;

    for (int q : Q) {
      int Nc = seq_intersect_cnt(L_prime, graph_->NeighborsR(q));
      if (Nc == L_prime.size()) {
        is_maximal = false;
        break;
      } else if (Nc >= min_l_size_)
        Q_prime.emplace_back(q);
    }
    if (is_maximal) {
      for (int p : P) {
        int Nc = seq_intersect_cnt(L_prime, graph_->NeighborsR(p));
        if (Nc == L_prime.size())
          R_prime.emplace_back(p);
        else if (Nc >= min_l_size_)
          P_prime.emplace_back(p);
      }

      if (!P_prime.empty() && (R_prime.size() + P_prime.size() >= min_r_size_)){
#ifdef COMPUTE_LEVEL
        cur_level_++;
#endif
        biclique_find(L_prime, R_prime, P_prime, Q_prime);
#ifdef COMPUTE_LEVEL
        cur_level_--;
#endif
      }
      if (R_prime.size() >= min_r_size_) {
        maximal_nodes_++;
        maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
#ifdef COMPUTE_LEVEL
        level_accumulation_ += cur_level_;
        max_level_ = std::max(cur_level_, max_level_);
#endif      
      }
    }
    Q.emplace_back(x);
  }
}

ImbeaFinder::ImbeaFinder(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {}

void ImbeaFinder::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);
  graph_->Reorder(RInc);
  for (int v = 0; v < graph_->GetRSize(); v++) {
    if (graph_->NeighborsR(v).size() >= min_l_size_) {
      processing_nodes_++;
      std::vector<int> L, R, P, Q;
      L = graph_->NeighborsR(v);

      bool is_maximal = true;

      std::map<int, int> c_map;

      for (int w : L) {
        for (int y : graph_->NeighborsL(w)) {
          if (++c_map[y] == L.size() && y < v) {
            is_maximal = false;
            break;
          }
        }
        if (!is_maximal)
          break;
      }
      if (!is_maximal)
        continue;
      std::vector<std::pair<int, int>> P_pairs;

      for (auto c_node : c_map) {
        if (c_node.second == L.size())
          R.emplace_back(c_node.first);
        else if (c_node.second >= min_l_size_) {
          if (c_node.first < v)
            Q.emplace_back(c_node.first);
          else
            P_pairs.emplace_back(c_node);
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
      std::sort(P_pairs.begin(), P_pairs.end(),
                [&](std::pair<int, int> x0, std::pair<int, int> x1) -> bool {
                  return x0.second > x1.second ||
                         (x0.second == x1.second && x0.first > x1.first);
                });
      for (auto p : P_pairs)
        P.emplace_back(p.first);
#ifdef COMPUTE_LEVEL
      cur_level_++;
#endif      
      biclique_find(L, R, P, Q);
#ifdef COMPUTE_LEVEL
      cur_level_--;
#endif      
    }
  }

  finish();
}

void ImbeaFinder::biclique_find(std::vector<int> L, std::vector<int> R,
                                   std::vector<int> P, std::vector<int> Q) {
  while (!P.empty()) {
    std::vector<int> L_prime, R_prime, P_prime, Q_prime, L_comple, C;
    bool is_maximal = true;

    int x = P.back();
    P.pop_back();
    L_prime = std::move(seq_intersect(L, graph_->NeighborsR(x)));
    L_comple = std::move(seq_except(L, L_prime));
    R_prime = R;
    R_prime.emplace_back(x);
    C.emplace_back(x);

    if (L_prime.size() < min_l_size_) {
      Q.emplace_back(x);
      continue;
    }

    processing_nodes_++;

    for (int q : Q) {
      int Nc = seq_intersect_cnt(L_prime, graph_->NeighborsR(q));
      if (Nc == L_prime.size()) {
        is_maximal = false;
        break;
      } else if (Nc >= min_l_size_)
        Q_prime.emplace_back(q);
    }

    if (is_maximal) {
      std::vector<std::pair<int, int>> p_prime_with_nc;
      for (int p : P) {
        int Nc = seq_intersect_cnt(L_prime, graph_->NeighborsR(p));
        if (Nc == L_prime.size()) {
          R_prime.emplace_back(p);
          int Nc_comple = seq_intersect_cnt(L_comple, graph_->NeighborsR(p));
          if (Nc_comple == 0)
            C.emplace_back(p);
        } else if (Nc >= min_l_size_) {
          p_prime_with_nc.emplace_back(std::make_pair(p, Nc));
        }
      }
      if (!p_prime_with_nc.empty() &&
          (R_prime.size() + p_prime_with_nc.size() >= min_r_size_)) {
        std::sort(p_prime_with_nc.begin(), p_prime_with_nc.end(),
                  [&](std::pair<int, int> x0, std::pair<int, int> x1) -> bool {
                    return x0.second > x1.second;
                  });
        P_prime.resize(p_prime_with_nc.size());
        for (int i = 0; i < p_prime_with_nc.size(); i++) {
          P_prime[i] = p_prime_with_nc[i].first;
        }
#ifdef COMPUTE_LEVEL
        cur_level_++;
#endif        
        biclique_find(L_prime, R_prime, P_prime, Q_prime);
#ifdef COMPUTE_LEVEL
        cur_level_--;
#endif        
      }

      if (R_prime.size() >= min_r_size_) {
        maximal_nodes_++;
        maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
#ifdef COMPUTE_LEVEL
        level_accumulation_ += cur_level_;
        max_level_ = std::max(cur_level_, max_level_);
#endif
      }
      P = std::move(seq_except(P, C));
      Q.insert(Q.end(), C.begin(), C.end());
    }
  }
}

MineLMBCFinder::MineLMBCFinder(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {}

void MineLMBCFinder::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);

  std::vector<int> X, GamaX, tailX;
  for (int i = 0; i < graph_->GetLSize(); i++)
    if (graph_->NeighborsL(i).size() >= min_r_size)
      GamaX.emplace_back(i);
  for (int i = 0; i < graph_->GetRSize(); i++)
    if (graph_->NeighborsR(i).size() >= min_l_size)
      tailX.emplace_back(i);
#ifdef COMPUTE_LEVEL
  cur_level_++;
#endif
  MineLMBC(X, GamaX, tailX);
#ifdef COMPUTE_LEVEL
  cur_level_--;
#endif
  finish();
}


FmbeFinder::FmbeFinder(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {}

void FmbeFinder::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);
  for (int i = 0; i < graph_->GetRSize(); i++) {
    if (graph_->NeighborsR(i).size() >= min_l_size_) {
      processing_nodes_++;
      std::vector<int> X, GamaX, tailX; // R, L, P
      std::set<int> tailX_set;
      X.emplace_back(i);
      GamaX = graph_->NeighborsR(i);
      for (int w : GamaX) {
        for (int y : graph_->NeighborsL(w)) {
          if (graph_->NeighborsR(y).size() > graph_->NeighborsR(i).size() ||
              (graph_->NeighborsR(y).size() == graph_->NeighborsR(i).size() &&
               y > i))
            tailX_set.insert(y);
        }
      }
      for (int element : tailX_set)
        tailX.emplace_back(element);
#ifdef COMPUTE_LEVEL
      cur_level_++;
#endif
      MineLMBC(X, GamaX, tailX);
#ifdef COMPUTE_LEVEL
      cur_level_--;
#endif

      if (min_l_size_ == 1 && graph_->NeighborsL(GamaX).size() == 1) {
        maximal_nodes_++;
        maximum_biclique_.CompareAndSet(GamaX, X, graph_);
#ifdef COMPUTE_LEVEL
        level_accumulation_ += cur_level_;
        max_level_ = std::max(cur_level_, max_level_);
#endif
      }
    }
  }
  finish();
}

PmbeBiGraph::PmbeBiGraph(const BiGraph &graph) : BiGraph(graph) {
  const int r_size = r_map_.size();
  std::vector<int> in_degree(r_size, 0);
  std::vector<std::vector<int>> child_r_nodes(r_size);
  std::vector<int> rename_map(r_size);
  typedef std::pair<int, int> IDPAIR;
  std::vector<int> ready_buf;
  auto rid_cmp = [=](int id0, int id1) -> bool {
    return r_adj_lists_[id0].size() > r_adj_lists_[id1].size();
  };
  std::queue<int> ready_q;
  range_index_.clear();
  range_index_.resize(r_size);
  for (int i = 0; i < r_size; i++) {
    for (int j = i + 1; j < r_size; j++) {
      int nc = seq_intersect_cnt(r_adj_lists_[i], r_adj_lists_[j]);
      if (nc == r_adj_lists_[j].size()) { //  i contains j
        in_degree[j]++;
        child_r_nodes[i].emplace_back(j);
      } else if (nc == r_adj_lists_[i].size()) {
        in_degree[i]++;
        child_r_nodes[j].emplace_back(i);
      }
    }
  }

  int scan_id = r_size;
  for (int i = 0; i < r_size; i++) {
    if (in_degree[i] == 0) // ready_q.push(i);
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
    for (int l : r_adj_lists_[r_id])
      n_l_adj_list[l].emplace_back(i);
    n_r_map[i] = r_map_[r_id];
    n_r_adj_list[i] = std::move(r_adj_lists_[r_id]);
  }
  std::swap(n_r_map, r_map_);
  std::swap(n_r_adj_list, r_adj_lists_);
  std::swap(n_l_adj_list, l_adj_lists_);
  n_l_adj_list.clear();
  n_r_adj_list.clear();
}

std::pair<int, int> PmbeBiGraph::GetRangeIndex(int v) {
  return range_index_[v];
}


PmbeFinder::PmbeFinder(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {
  double start = get_cur_time();
  pgraph_ = new PmbeBiGraph(*graph_in);
  //printf("RevOrder and CDAG construction time: %lfs\n", get_cur_time() - start);
}

void PmbeFinder::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);

  for (int v = 0; v < pgraph_->GetRSize(); v++) {
    if (pgraph_->NeighborsR(v).size() >= min_l_size_) {
      processing_nodes_++;
      std::vector<int> L, R, C;
      L = pgraph_->NeighborsR(v);
      R = pgraph_->NeighborsL(L);
      if (R[0] < v)
        continue;
      if (R.size() >= min_r_size_) {
        maximal_nodes_++;
        maximum_biclique_.CompareAndSet(L, R, pgraph_);
#ifdef COMPUTE_LEVEL
        level_accumulation_ += cur_level_;
        max_level_ = std::max(cur_level_, max_level_);
#endif
      }
      std::map<int, int> c_map;
      for (int w : L) {
        for (int y : pgraph_->NeighborsL(w)) {
          if (y > v)
            c_map[y]++;
        }
      }
      for (auto c_node : c_map)
        if (c_node.second != L.size())
          C.emplace_back(c_node.first);
#ifdef COMPUTE_LEVEL
      cur_level_++;
#endif
      biclique_find(L, R, C);
#ifdef COMPUTE_LEVEL
      cur_level_--;
#endif
    }
  }

  // std::vector<int> L, R, C;
  // for (int i = 0; i < pgraph_->GetLSize(); i++)
  //   if (pgraph_->NeighborsL(i).size() >= min_r_size_)
  //     L.emplace_back(i);
  // for (int i = 0; i < pgraph_->GetRSize(); i++)
  //   if (pgraph_->NeighborsR(i).size() >= min_l_size_)
  //     C.emplace_back(i);
  // biclique_find(L, R, C);
  finish();
}

void PmbeFinder::biclique_find(std::vector<int> L, std::vector<int> R,
                                  std::vector<int> C) {
  std::queue<std::pair<int, int>> ex_q;
  while (!C.empty()) {
    std::vector<int> L_prime, R_prime, C_prime;
    bool is_maximal = true;
    int x = C.back();
    auto p = pgraph_->GetRangeIndex(x);
    if (p.first <= p.second)
      ex_q.push(p);
    while (!ex_q.empty() && ex_q.front().first > x)
      ex_q.pop();
    if (!ex_q.empty() && ex_q.front().second >= x) {
      C.pop_back();
      continue;
    }

    L_prime = seq_intersect(L, pgraph_->NeighborsR(x));

    if (L_prime.size() < min_l_size_) {
      C.pop_back();
      continue;
    }
    processing_nodes_++;

    // maximality check
    R_prime = pgraph_->NeighborsL(L_prime);
    auto R_add = seq_except(R_prime, R);
    if (seq_intersect_cnt(R_add, C) != R_add.size())
      is_maximal = false;

    // generate new node
    if (is_maximal) {
      for (int i = 0; i < C.size(); i++) {
        int neighbor_cnt =
                seq_intersect_cnt(L_prime, pgraph_->NeighborsR(C[i])),
            v = C[i];
        if (v == x)
          continue;
        if (neighbor_cnt > 0 && neighbor_cnt < L_prime.size())
          C_prime.push_back(v);
      }

      if ((!C_prime.empty()) &&
          (R_prime.size() + C_prime.size() >= min_r_size_)) {
#ifdef COMPUTE_LEVEL
        cur_level_++;
#endif
        biclique_find(L_prime, R_prime, C_prime);
#ifdef COMPUTE_LEVEL
        cur_level_--;
#endif
      }

      if (R_prime.size() >= min_r_size_) {
        maximal_nodes_++;
        maximum_biclique_.CompareAndSet(L_prime, R_prime, pgraph_);
#ifdef COMPUTE_LEVEL
        level_accumulation_ += cur_level_;
        max_level_ = std::max(cur_level_, max_level_);
#endif
      }
    }
    C.pop_back();
  }
}
