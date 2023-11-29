#include "QuasiMBE.h"

QuasiMBEFinder::QuasiMBEFinder(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {
  miu_ = 1;
}

void QuasiMBEFinder::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);
  graph_->Reorder(RDec);

  std::vector<StateEnum> v_states(graph_->GetRSize(), P);
  for (int v = 0; v < graph_->GetRSize(); v++) {
    if (v_states[v] != P || graph_->NeighborsR(v).size() < min_l_size_)
      continue;
    processing_nodes_++;

    VertexSet &L = graph_->NeighborsR(v);
    VertexSet R;
    VertexSetExt C;
    int R_ext_cnt = 0;

    std::map<int, int> c_map;
    for (int l : L) {
      for (int r : graph_->NeighborsL(l)) {
        c_map[r]++;
      }
    }
    for (auto &c_node : c_map) {
      if (c_node.second == L.size()) {
        R.emplace_back(c_node.first);
        v_states[c_node.first] = (c_node.first == v) ? Q : D;
      } else if (c_node.second >= min_l_size_ && v_states[c_node.first] != D) {
        VertexExtNode node;
        node.v = c_node.first;
        node.LN = c_node.second;
        node.state = v_states[c_node.first];
        if (node.LN >= graph_->NeighborsR(v).size() * miu_) R_ext_cnt++;

        C.emplace_back(node);
        if (graph_->NeighborsR(c_node.first).size() == c_node.second)
          v_states[c_node.first] = D;
      }
    }
    if (R.size() >= min_r_size_) {
      maximal_nodes_++;
      updateMB(L.size(), R.size() + R_ext_cnt);
    }

    biclique_find(L, R, C);
  }
  // VertexSet L, R;
  // VertexSetExt C;

  // for (int i = 0; i < graph_->GetLSize(); i++)
  //   L.emplace_back(i);
  // for (int i = 0; i < graph_->GetRSize(); i++) {
  //   VertexExtNode node;
  //   node.v = i;
  //   node.LN = graph_->NeighborsR(i).size();
  //   node.state = P;
  //   C.emplace_back(std::move(node));
  // }
  // biclique_find(L, R, C);

  // graph_->Reorder(2);

  // for (int v = 0; v < graph_->GetRSize(); v++) {
  //   if (graph_->NeighborsR(v).size() < min_l_size_)
  //     continue;
  //   if (v > 0 && graph_->NeighborsR(v - 1) == graph_->NeighborsR(v))
  //     continue;
  //   processing_nodes_++;

  //   VertexSet &L = graph_->NeighborsR(v);
  //   VertexSet R, R_ext, P, Q;

  //   std::map<int, int> c_map;
  //   for (int w : L) {
  //     for (int y : graph_->NeighborsL(w)) {
  //       c_map[y]++;
  //     }
  //   }

  //   std::vector<std::pair<int, int>> p_with_nc;
  //   for (auto &c_node : c_map) {
  //     if (c_node.second == L.size()) {
  //       R.emplace_back(c_node.first);
  //     } else if (c_node.second >= min_l_size_) {
  //       if (c_node.first < v) {
  //         Q.emplace_back(c_node.first);
  //       } else if (c_node.second >= L.size() * miu_) {
  //         R_ext.emplace_back(c_node.first);
  //       } else {
  //         p_with_nc.emplace_back(std::move(c_node));
  //       }
  //     }
  //   }
  //   if(R.size() + R_ext.size() >= min_r_size_){
  //     maximal_nodes_++;
  //     updateMB(L.size(), R.size() + R_ext.size());
  //   }
  //   std::sort(p_with_nc.begin(), p_with_nc.end(),
  //             [&](std::pair<int, int> x0, std::pair<int, int> x1) -> bool {
  //               return x0.second < x1.second ||
  //                      (x0.second == x1.second && x0.first > x1.first);
  //             });
  //   for(auto& p:p_with_nc){
  //     P.emplace_back(p.first);
  //   }
  //   biclique_find(L, R, R_ext, P, Q);
  // }
  finish();
  printf("%d %d\n", l_mb_, r_mb_);
}

void QuasiMBEFinder::setMiu(double miu_in) { miu_ = miu_in; }

void QuasiMBEFinder::biclique_find(VertexSet &L, VertexSet &R,
                                   VertexSetExt &C) {
  auto update_LN = [=](int l, std::vector<int> &LN) {
    auto &l_neighbors = graph_->NeighborsL(l);
    for (int i = 0, j = 0; i < l_neighbors.size() && j < C.size();) {
      if (l_neighbors[i] > C[j].v)
        j++;
      else if (l_neighbors[i] < C[j].v)
        i++;
      else {
        LN[j++]++;
        i++;
      }
    }
  };

  for (int i = 0; i < C.size(); i++) {
    if (C[i].state == P && C[i].LN < L.size() * miu_) {
      processing_nodes_++;
      std::vector<int> LN(C.size(), 0);
      VertexSet L_prime = seq_intersect(L, graph_->NeighborsR(C[i].v));
      for (int l : L_prime)
        update_LN(l, LN);

      // pivot check
      int pivot_id = i;
      for (int j = i + 1; j < C.size(); j++) {
        if (LN[i] == LN[j] && C[j].state == P && C[j].LN > C[pivot_id].LN &&
            C[j].LN < L.size() * miu_)
          pivot_id = j;
      }
      if (pivot_id != i) {
        VertexSet L_prime_orig =
            seq_intersect(L, graph_->NeighborsR(C[pivot_id].v));
        std::swap(L_prime, L_prime_orig);
        VertexSet L_add = seq_except(L_prime, L_prime_orig);
        for (int l : L_add)
          update_LN(l, LN);
      }

      // processing
      bool is_maximal = true;
      for (int j = 0; j < C.size(); j++) {
        if (LN[j] == LN[pivot_id] && C[j].state != P) {
          is_maximal = false;
          break;
        }
      }
      if (is_maximal) {
        VertexSet R_prime = R;
        VertexSetExt C_prime;
        int R_ext_cnt = 0;
        for (int j = 0; j < C.size(); j++) {
          if (LN[j] == LN[pivot_id])
            R_prime.emplace_back(C[j].v);
          else if (LN[j] >= min_l_size_ && C[j].state != D) {
            VertexExtNode node = C[j];
            node.LN = LN[j];
            if (node.state == P && node.LN >= L_prime.size()*miu_)
              R_ext_cnt++;
            C_prime.emplace_back(std::move(node));
          }
        }
        if (R_prime.size() >= min_r_size_) {
          if (maximal_nodes_ % 100000000 == 0) {
            printf("%llu %lfs\n", (long long)maximal_nodes_,
                   get_cur_time() - start_time_);
          }
          maximal_nodes_++;
          updateMB(L_prime.size(), R_prime.size() + R_ext_cnt);
        }
        biclique_find(L_prime, R_prime, C_prime);
      }
      // after processing
      for (int j = 0; j < C.size(); j++)
        if (LN[j] == C[j].LN)
          C[j].state = D;
      C[pivot_id].state = Q;
    }
  }
}

void QuasiMBEFinder::updateMB(int l_mb_in, int r_mb_in) {
  if (l_mb_in * r_mb_in > l_mb_ * r_mb_) {
    l_mb_ = l_mb_in;
    r_mb_ = r_mb_in;
  }
}
