#include "BoundMBE.h"

BoundMBEFinder::BoundMBEFinder(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {}

void BoundMBEFinder::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);
  node_trim_time_ = 0;
  intersect_time_ = 0;
  graph_->Reorder(RInc);

  // VertexSet l_degrees;
  // for (int u = 0; u < graph_->GetLSize(); u++) {
  //   l_degrees.emplace_back(graph_->NeighborsL(u).size());
  // }
  for (int v = 0; v < graph_->GetRSize(); v++) {
    // if (v % 1000 == 0) printf("%d %lf\n", v, get_cur_time() - start_time_);

    VertexExtSet L, C;
    VertexSet R;
    // for (int u : graph_->NeighborsR(v)) {
    //   if (l_degrees[u] >= min_r_size_)
    //     L.emplace_back(u, l_degrees[u]);
    //   l_degrees[u]--;
    // }
    if (graph_->NeighborsR(v).size() < min_l_size_)
      continue;
    if (v > 0 && graph_->NeighborsR(v - 1) == graph_->NeighborsR(v))
      continue;
    processing_nodes_++;

    for (int u : graph_->NeighborsR(v)) {
      auto &l_neighbors = graph_->NeighborsL(u);
      if (l_neighbors.size() >= min_r_size_ &&
          l_neighbors[l_neighbors.size() - min_r_size_] >= v)
        L.emplace_back(std::make_pair(u, 0));
    }
    R.emplace_back(v);

    std::map<int, int> c_map;
    for (auto &lp : L) {
      auto &l_neighbors = graph_->NeighborsL(lp.first);
      for (int i = l_neighbors.size() - 1; l_neighbors[i] > v; i--) {
        c_map[l_neighbors[i]]++;
      }
    }

    for (auto &cp : c_map) {
      if (cp.second == L.size())
        R.emplace_back(cp.first);
      else if (cp.second >= min_l_size_)
        C.emplace_back(std::move(cp));
    }
    trim_node_2(L, R, C);
    if (L.size() < min_l_size_)
      continue;

    if (L.size() == graph_->NeighborsR(v).size() || maximality_check(L, R)) {
      if (R.size() >= min_r_size_) {
        maximal_nodes_++;
        maximum_biclique_.CompareAndSet(L, R, graph_);
      }
      biclique_find(L, R, C);
    }
  }

  finish();
  printf("node_trim_time: %lfs\n", node_trim_time_);
  printf("intersect_time: %lfs\n", intersect_time_);
}

void BoundMBEFinder::biclique_find(VertexExtSet &L, VertexSet &R,
                                   VertexExtSet &C) {

  VertexSet LN(C.size(), 0);
  while (L.size() >= min_l_size_ && R.size() + C.size() >= min_r_size_ &&
         !C.empty()) {
    processing_nodes_++;

    VertexExtSet L_prime, C_prime;
    VertexSet R_prime = R;

    int index = 0;
    for (int i = 1; i < C.size(); i++) {
      if (C[index].second > C[i].second)
        index = i;
    }

    auto &c_neighbors = graph_->NeighborsR(C[index].first);
    

    VertexSet LL;
    for (auto p : L) LL.emplace_back(p.first);
    double start = get_cur_time();
    auto x = seq_intersect(LL, graph_->NeighborsR(C[index].first));

    intersect_time_ += get_cur_time() - start;
    for (int li = 0, ci = 0; li < L.size() && ci < c_neighbors.size();) {
      if (L[li].first < c_neighbors[ci])
        li++;
      else if (L[li].first > c_neighbors[ci])
        ci++;
      else {
        L_prime.emplace_back(L[li].first, R_prime.size());
        auto &l_neighbors = graph_->NeighborsL(L[li].first);
        for (int lj = l_neighbors.size() - 1, cj = C.size() - 1;
             lj >= 0 && cj >= 0;) {
          if (l_neighbors[lj] > C[cj].first)
            lj--;
          else if (l_neighbors[lj] < C[cj].first)
            cj--;
          else {
            LN[cj]++;
            // L_prime.back().second++;
            lj--;
            cj--;
          }
        }
        li++;
        ci++;
      }
    }

    for (int ci = 0; ci < C.size(); ci++) {
      if (LN[ci] == L_prime.size())
        R_prime.emplace_back(C[ci].first);
      else if (LN[ci] >= min_l_size_)
        C_prime.emplace_back(std::make_pair(C[ci].first, LN[ci]));
    }

    trim_node(L_prime, R_prime, C_prime);

    if (L_prime.size() >= min_l_size_ && maximality_check(L_prime, R_prime)) {
      if (R_prime.size() >= min_r_size_) {
        maximal_nodes_++;
        
        maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
      }
      biclique_find(L_prime, R_prime, C_prime);
    }

    int r_size = R.size();

    for (int ci = 0; ci < C.size(); ci++) {
      if (LN[ci] == C[ci].second) {
        C[ci].second = 0;
      }
      LN[ci] = 0;
    }

    trim_node(L, R, C);

    if (R.size() != r_size && L.size() >= min_l_size_) {
      if (maximality_check(L, R)) {
        if (R.size() >= min_r_size_) {
          maximal_nodes_++;
          maximum_biclique_.CompareAndSet(L, R, graph_);
        }
      } else {
        break;
      }
    }
  }
}

void BoundMBEFinder::trim_node(VertexExtSet &L, VertexSet &R, VertexExtSet &C) {

  if (min_l_size_ == 1 && min_r_size_ == 1) {
    int c_valid = 0;
    for (int i = 0; i < C.size(); i++) {
      if (C[i].second == L.size())
        R.emplace_back(C[i].first);
      else if (C[i].second >= min_l_size_) {
        if (c_valid != i)
          C[c_valid] = C[i];
        c_valid++;
      }
    }
    C.resize(c_valid);
    return;
  }
  bool trim_flag = false;
  double start = get_cur_time();
  // initialization
  int l_valid, c_valid = 0;
  for (int i = 0; i < C.size(); i++) {
    if (C[i].second == L.size())
      R.emplace_back(C[i].first);
    else if (C[i].second >= min_l_size_) {
      if (c_valid != i)
        C[c_valid] = C[i];
      c_valid++;
    }
  }
  C.resize(c_valid);
  for (int i = 0; i < L.size(); i++) {
    auto &l_neighbors = graph_->NeighborsL(L[i].first);
    int LN = R.size();
    for (int li = 0, ci = 0; li < l_neighbors.size() && ci < C.size();) {
      if (l_neighbors[li] < C[ci].first)
        li++;
      else if (l_neighbors[li] > C[ci].first)
        ci++;
      else {
        LN++;
        li++;
        ci++;
      }
    }
    if (LN < min_r_size_)
      trim_flag = true;
    L[i].second = LN;
  }

  while (trim_flag) {
    trim_flag = false;
    l_valid = 0;
    for (int i = 0; i < L.size(); i++) {
      if (L[i].second >= min_r_size_) {
        if (l_valid != i)
          L[l_valid] = L[i];
        l_valid++;
      } else {
        auto &l_neighbors = graph_->NeighborsL(L[i].first);
        for (int li = 0, ci = 0; li < l_neighbors.size() && ci < C.size();) {
          if (l_neighbors[li] > C[ci].first)
            ci++;
          else if (l_neighbors[li] < C[ci].first)
            li++;
          else {
            C[ci].second--;
            li++;
            ci++;
          }
        }
      }
    }
    L.resize(l_valid);
    if (l_valid < min_l_size_)
      return;

    c_valid = 0;
    for (int i = 0; i < C.size(); i++) {
      if (C[i].second == L.size())
        R.emplace_back(C[i].first);
      else if (C[i].second >= min_l_size_) {
        if (c_valid != i)
          C[c_valid] = C[i];
        c_valid++;
      } else {
        auto &c_neighbors = graph_->NeighborsR(C[i].first);
        for (int li = 0, ci = 0; li < L.size() && ci < c_neighbors.size();) {
          if (L[li].first < c_neighbors[ci])
            li++;
          else if (L[li].first > c_neighbors[ci])
            ci++;
          else {
            if (--L[li].second < min_r_size_)
              trim_flag = true;
            li++;
            ci++;
          }
        }
      }
    }
    C.resize(c_valid);

    if (c_valid + R.size() < min_r_size_)
      return;
  }
  node_trim_time_ += get_cur_time() - start;
}

void BoundMBEFinder::trim_node_2(VertexExtSet &L, VertexSet &R,
                                 VertexExtSet &C) {
  if (min_l_size_ == 1 && min_r_size_ == 1)
    return;
  double start = get_cur_time();
  
    // initialization
  int c_valid = 0;
  for (int i = 0; i < C.size(); i++) {
    if (C[i].second == L.size())
      R.emplace_back(C[i].first);
    else if (C[i].second >= min_l_size_) {
      if (c_valid != i)
        C[c_valid] = C[i];
      c_valid++;
    }
  }
  C.resize(c_valid);

  std::vector<VertexSet> localNeighborsL(L.size());
  for (int i = 0; i < L.size(); i++) {
    auto &l_neighbors = graph_->NeighborsL(L[i].first);
    for (int li = 0, ci = 0; li < l_neighbors.size() && ci < C.size();) {
      if (l_neighbors[li] > C[ci].first)
        ci++;
      else if (l_neighbors[li] < C[ci].first)
        li++;
      else {
        localNeighborsL[i].emplace_back(ci);
        li++;
        ci++;
      }
    }
  }
  VertexSet LN(C.size(), 0);
  for (int i = 0; i < C.size(); i++) {
    auto &c_neighbors = graph_->NeighborsR(C[i].first);
    for (int li = 0, ci = 0; li < L.size() && ci < c_neighbors.size();) {
      if (L[li].first > c_neighbors[ci])
        ci++;
      else if (L[li].first < c_neighbors[ci])
        li++;
      else {
        for (int id : localNeighborsL[li])
          LN[id]++;
        li++;
        ci++;
      }
    }
    int cnt = 0;
    for (int j = 0; j < C.size(); j++) {
      if (LN[j] + R.size() >= min_r_size_)
        cnt++;
      LN[j] = 0;
    }
    if (R.size() + cnt < min_r_size_)
      C[i].second = 0;
  }

  node_trim_time_ += get_cur_time() - start;
  trim_node(L, R, C);
}

void BoundMBEFinder::trim_node_max_edge(VertexExtSet &L, VertexSet &R,
                                        VertexExtSet &C) {
  if (L.size() == 0) return;
  bool trim_flag = false;
  double start = get_cur_time();
  int bound = (maximum_biclique_.GetEdges() - 1) / (R.size() + C.size()) + 1;
  int cur_min_l = std::max(min_l_size_, bound);
  int l_valid, c_valid = 0;
  for (int i = 0; i < C.size(); i++) {
    if (C[i].second == L.size())
      R.emplace_back(C[i].first);
    else if (C[i].second >= cur_min_l) {
      if (c_valid != i)
        C[c_valid] = C[i];
      c_valid++;
    }
  }
  C.resize(c_valid);
  bound = (maximum_biclique_.GetEdges() - 1) / (R.size() + C.size()) + 1;
  cur_min_l = std::max(min_l_size_, bound);
  bound = (maximum_biclique_.GetEdges() - 1) / L.size() + 1;
  int cur_min_r = std::max(min_r_size_, bound);
  if (cur_min_l == 1 && cur_min_r == 1) return;

  for (int i = 0; i < L.size(); i++) {
    auto &l_neighbors = graph_->NeighborsL(L[i].first);
    int LN = R.size();
    for (int li = 0, ci = 0; li < l_neighbors.size() && ci < C.size();) {
      if (l_neighbors[li] < C[ci].first)
        li++;
      else if (l_neighbors[li] > C[ci].first)
        ci++;
      else {
        LN++;
        li++;
        ci++;
      }
    }
    if (LN < cur_min_r)
      trim_flag = true;
    L[i].second = LN;
  }

  while (trim_flag) {
    trim_flag = false;
    l_valid = 0;
    for (int i = 0; i < L.size(); i++) {
      if (L[i].second >= cur_min_r) {
        if (l_valid != i)
          L[l_valid] = L[i];
        l_valid++;
      } else {
        auto &l_neighbors = graph_->NeighborsL(L[i].first);
        for (int li = 0, ci = 0; li < l_neighbors.size() && ci < C.size();) {
          if (l_neighbors[li] > C[ci].first)
            ci++;
          else if (l_neighbors[li] < C[ci].first)
            li++;
          else {
            C[ci].second--;
            li++;
            ci++;
          }
        }
      }
    }
    if(L.size() != l_valid && l_valid > 0){
      bound = (maximum_biclique_.GetEdges() - 1) / l_valid + 1;
      cur_min_r = std::max(cur_min_r, bound);
    } else
      L.resize(l_valid);
    if (l_valid < cur_min_l)
      return;

    c_valid = 0;
    for (int i = 0; i < C.size(); i++) {
      if (C[i].second == L.size())
        R.emplace_back(C[i].first);
      else if (C[i].second >= cur_min_l) {
        if (c_valid != i)
          C[c_valid] = C[i];
        c_valid++;
      } else {
        auto &c_neighbors = graph_->NeighborsR(C[i].first);
        for (int li = 0, ci = 0; li < L.size() && ci < c_neighbors.size();) {
          if (L[li].first < c_neighbors[ci])
            li++;
          else if (L[li].first > c_neighbors[ci])
            ci++;
          else {
            if (--L[li].second < cur_min_r)
              trim_flag = true;
            li++;
            ci++;
          }
        }
      }
    }
    C.resize(c_valid);

    if (c_valid + R.size() < cur_min_r)
      return;
  }
  node_trim_time_ += get_cur_time() - start;
}

bool BoundMBEFinder::maximality_check(VertexExtSet &L, VertexSet &R) {
  VertexSet R_tmp = graph_->NeighborsL(L[0].first);
  if (L.size() < min_l_size_)
    return false;

  for (int i = 1; i < L.size(); i++) {
    R_tmp = seq_intersect(R_tmp, graph_->NeighborsL(L[i].first));
  }
  if (R_tmp.size() > R.size()) {
    return false;
  } else {
    R = R_tmp;
    return true;
  }
}

void BoundMBEFinder::print_node(VertexExtSet &L, VertexSet &R,
                                VertexExtSet &C) {
  printf("L(%ld):", L.size());
  for (auto l : L) {
    printf("[%d %d]", l.first, l.second);
  }
  printf("\nR(%ld):", R.size());

  for (auto r : R)
    printf("[%d]", r);
  printf("\nC(%ld):", C.size());

  for (auto c : C) {
    printf("[%d %d]", c.first, c.second);
  }
  printf("\n\n");
}
