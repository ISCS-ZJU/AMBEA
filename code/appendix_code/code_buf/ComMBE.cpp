#include "ComMBE.h"
#include <map>
#include <unordered_set>

ComMBEFinder::ComMBEFinder(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {}

void ComMBEFinder::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);
  r_clock_ = 0;
  iter_clock_ = 0;
  gen_cand_clock_ = 0;
  r_clock_mb_ = 0;
  graph_ = new BiGraph(*graph_);

  graph_->Reorder(RInc);
  l_valid = std::vector<bool>(graph_->GetLSize(), false);
  r_valid = std::vector<bool>(graph_->GetRSize(), false);

  // std::vector<int> R;

  // for (int v = std::max(0, graph_->GetRSize() - N); v < graph_->GetRSize();
  // v++)
  //   R.emplace_back(v);
  // setup_bs(R);
  // bitset_t R_bs;
  // printf("start iteration %lfs %d %d\n", get_cur_time() - start_time_,
  //        R_size_bs_, (int)L_bitsets_.size());
  // std::vector<int> L, C_rev;
  // for (int i = 0; i < L_bitsets_.size(); i++)
  //   L.emplace_back(i);
  // for (int i = R_size_bs_ - 1; i >= 0; i--)
  //   C_rev.emplace_back(i);
  // biclique_find_bs(L, R_bs, C_rev);

  for (int v = 0; v < graph_->GetRSize(); v++) {
    cur_v_ = v;
    if (graph_->NeighborsR(v).size() < min_l_size_)
      continue;
    if (v > 0 && graph_->NeighborsR(v - 1) == graph_->NeighborsR(v))
      continue;
    processing_nodes_++;

    VertexSet &L = graph_->NeighborsR(v);
    VertexSet R;
    std::vector<Node *> C;
    std::map<int, std::vector<int>> c_map;
    for (int l : L) {
      auto &l_neighbors = graph_->NeighborsL(l);
      for (int i = l_neighbors.size() - 1; l_neighbors[i] > v; i--) {
        c_map[l_neighbors[i]].emplace_back(l);
      }
    }
    R.emplace_back(v);
    for (auto &c_node : c_map) {

      if (c_node.second.size() == L.size())
        R.emplace_back(c_node.first);
      else if (c_node.second.size() >= min_l_size) {
        C.emplace_back(new Node(c_node.first, std::move(c_node.second)));
      }
    }

    if (R.size() >= min_r_size_) {
      if (maximal_nodes_ % 10000000 == 0) {
        printf("maximal biclique:%lld time:%lfs   %d/%d\n",
               (long long)maximal_nodes_, get_cur_time() - start_time_, cur_v_,
               graph_->GetRSize());
      }
      for (int l : L)
        l_valid[l] = true;
      for (int r : R)
        r_valid[r] = true;
      maximal_nodes_++;
      maximum_biclique_.CompareAndSet(L, R, graph_);
    }
    double start = get_cur_time();
    biclique_find(L, R, C);
    iter_clock_ += get_cur_time() - start;
    for (auto ptr : C)
      delete ptr;
  }

  finish();
  int l_valid_cnt = 0, r_valid_cnt = 0;
  for (auto flag : l_valid)
    if (flag)
      l_valid_cnt++;
  for (auto flag : r_valid)
    if (flag)
      r_valid_cnt++;
  printf("%d/%zu %d/%zu \n", l_valid_cnt, l_valid.size(), r_valid_cnt,
         r_valid.size());
  printf("Time for computing R: %lf s\n", r_clock_);
  printf("Time for computing R (non maximal): %lf s\n", r_clock_mb_);
  printf("Time for iteration: %lf s\n", iter_clock_);
  printf("Time for generate candidates: %lf s\n", gen_cand_clock_);
}

void ComMBEFinder::biclique_find(VertexSet &L, VertexSet &R,
                                 std::vector<Node *> &C) {
  std::vector<int> LN(C.size(), 0);
  std::vector<Node *> ref_node(C.size(), nullptr);
  std::vector<Node *> create_node_ptrs;

  // suppose that C[i].L is OK.
  for (int i = 0; i < C.size(); i++) {
    if (!C[i]->R.empty())
      continue;
    processing_nodes_++;

    double R_time = get_cur_time();
    if (ref_node[i] == nullptr) {
      C[i]->R = graph_->NeighborsL(C[i]->L);
    } else {
      auto L_add = seq_except(C[i]->L, ref_node[i]->L);
      C[i]->R = seq_intersect(ref_node[i]->R, graph_->NeighborsL(L_add));
    }

    R_time = get_cur_time() - R_time;
    r_clock_ += R_time;

    int diff_element = first_diff_element(C[i]->R, R);
    if (diff_element == C[i]->vc) { // is_maximal
      if (C[i]->R.size() >= min_r_size_) {
        if (maximal_nodes_ % 10000000 == 0) {
          printf("maximal biclique:%lld time:%lfs   %d/%d\n",
                 (long long)maximal_nodes_, get_cur_time() - start_time_,
                 cur_v_, graph_->GetRSize());
        }
        for (int l : C[i]->L)
          l_valid[l] = true;
        for (int r : C[i]->R)
          r_valid[r] = true;
        maximal_nodes_++;
        maximum_biclique_.CompareAndSet(C[i]->L, C[i]->R, graph_);
      }
      std::vector<Node *> C_prime;
      double start = get_cur_time();
      if (C[i]->L.size() * 2 < (C.size() - i)) {
        for (int l : C[i]->L) {
          auto &l_neighbors = graph_->NeighborsL(l);
          for (int j = l_neighbors.size() - 1, k = C.size() - 1;
               j >= 0 && k >= 0 && l_neighbors[j] > C[i]->vc;) {
            if (l_neighbors[j] > C[k]->vc) {
              j--;
            } else if (l_neighbors[j] < C[k]->vc) {
              k--;
            } else {
              LN[k--]++;
              j--;
            }
          }
        }
      } else {
        for (int j = i + 1; j < C.size(); j++) {
          LN[j] = seq_intersect_cnt(C[i]->L, C[j]->L);
        }
      }

      for (int j = i + 1; j < C.size(); j++) {
        int ln = LN[j];
        LN[j] = 0;
        if (!C[j]->R.empty())
          continue;
        if (ln >= min_l_size_ && ln != C[i]->L.size()) {
          if (ln == C[j]->L.size()) {
            C_prime.emplace_back(C[j]);
          } else {
            Node *node = new Node(C[j]->vc);
            node->L = seq_intersect(C[i]->L, C[j]->L);
            C_prime.emplace_back(node);
            create_node_ptrs.emplace_back(node);
          }
        }
        if (ln == C[i]->L.size() && ln == C[j]->L.size() && C[j]->R.empty()) {
          C[j]->R.emplace_back(-1); //= C[i]->R;
        } else if (ln == C[i]->L.size() && ln > 2 &&
                   (ref_node[j] == nullptr || ref_node[j]->L.size() < ln)) {
          ref_node[j] = C[i];
        }
      }
      gen_cand_clock_ += get_cur_time() - start;
      if (C[i]->R.size() + C_prime.size() >= min_r_size_)
        biclique_find(C[i]->L, C[i]->R, C_prime);
    } else {
      double start = get_cur_time();
      auto superL = seq_intersect(L, graph_->NeighborsR(diff_element));
      if (superL.size() * 2 < C.size() - i) {
        for (int l : superL) { // C[i]->L) {
          auto &l_neighbors = graph_->NeighborsL(l);
          for (int j = l_neighbors.size() - 1, k = C.size() - 1;
               j >= 0 && k >= 0 && l_neighbors[j] > C[i]->vc;) {
            if (l_neighbors[j] > C[k]->vc) {
              j--;
            } else if (l_neighbors[j] < C[k]->vc) {
              k--;
            } else {
              LN[k--]++;
              j--;
            }
          }
        }
      } else {
        for (int j = i + 1; j < C.size(); j++) {
          LN[j] = seq_intersect_cnt(superL, C[j]->L);
        }
      }

      for (int j = i + 1; j < C.size(); j++) {
        int ln = LN[j];
        LN[j] = 0;
        if (!C[j]->R.empty())
          continue;
        if (ln == C[j]->L.size())
          C[j]->R.emplace_back(-1);
      }

      gen_cand_clock_ += get_cur_time() - start;
      r_clock_mb_ += R_time;
    }
  }
  for (auto ptr : create_node_ptrs)
    delete ptr;
}

void ComMBEFinder::setup_bs(std::vector<int> &R) {
  if (R.size() > N)
    printf("Error\n");
  R_size_bs_ = 0;
  std::map<int, bitset_t> mm;
  for (int i = 0; i < R.size(); i++) {
    if (i != 0 && graph_->NeighborsR(R[i]) == graph_->NeighborsR(R[i - 1]))
      continue;
    R_size_bs_++;
    for (int l : graph_->NeighborsR(R[i])) {
      mm[l].set(i);
    }
  }
  R_neighbors_bs_.clear();
  R_neighbors_bs_.resize(R_size_bs_);
  std::vector<bitset_t> L_bitsets_tmp;
  for (auto &p : mm)
    L_bitsets_tmp.emplace_back(std::move(p.second));
  printf("Ordering start\n");
  std::sort(L_bitsets_tmp.begin(), L_bitsets_tmp.end(),
            [&](bitset_t bs0, bitset_t bs1) -> bool {
              auto v0 = bs0.to_ullong();
              auto v1 = bs1.to_ullong();
              return v0 < v1;
              // bitset_t diff = bs0 ^ bs1;
              // if (diff.count() == 0)
              //   return false;
              // for (int i = 0; i < N; i = i * 2)
              //   diff |= i;
              // diff ^= diff >> 1;
              // return (bs0 & diff).count() != 0;
            });
  printf("Ordering end\n");
  for (int i = 0; i < L_bitsets_tmp.size(); i++) {
    if (i != 0 && (L_bitsets_tmp[i - 1] == L_bitsets_tmp[i]))
      continue;
    for (int j = 0; j < R_size_bs_; j++) {
      if (L_bitsets_tmp[i][j])
        R_neighbors_bs_[j].emplace_back(L_bitsets_.size());
    }
    L_bitsets_.emplace_back(std::move(L_bitsets_tmp[i]));
  }

  // for(auto &p: mm){
  //   for (int i = 0; i < R_size_bs_; i++) {
  //     if (p.second[i]) R_neighbors_bs_[i].emplace_back(L_bitsets_.size());
  //   }
  //   L_bitsets_.emplace_back(std::move(p.second));
  // }
}

void ComMBEFinder::biclique_find_bs(bitset_t R_bs, int R_start) {
  for (int i = R_start; i < R_size_bs_; i++) {
    if (R_bs[i] == 1)
      continue;
    bitset_t R_bs_current = R_bs, R_bs_next;
    int r_cnt = 0;
    for (int l : R_neighbors_bs_[i]) {
      if ((L_bitsets_[l] & R_bs_current) == R_bs_current) {
        if (r_cnt == 0)
          R_bs_next = L_bitsets_[l];
        else
          R_bs_next &= L_bitsets_[l];
        r_cnt++;
      }
    }
    processing_nodes_++;
    if (r_cnt != 0 && ((R_bs_next & (~R_bs_current)) << (N - i)).count() == 0) {
      if (maximal_nodes_ % 100000000 == 0) {
        printf("%llu %lfs\n", (long long)maximal_nodes_,
               get_cur_time() - start_time_);
      }
      maximal_nodes_++;
      biclique_find_bs(R_bs_next, i + 1);
    }
  }
}

void ComMBEFinder::biclique_find_bs(std::vector<int> &L, bitset_t R_bs,
                                    std::vector<int> &C_rev) {
  std::vector<int> C_prime_rev;
  for (int i = 0; i < C_rev.size(); i++) {
    if (R_bs.count() == 0)
      cur_v_ = i;
    if (R_bs[C_rev[i]] == 1)
      continue;

    std::vector<int> L_prime = seq_intersect(L, R_neighbors_bs_[C_rev[i]]);
    if (L_prime.empty())
      continue;
    processing_nodes_++;
    bitset_t R_bs_next = L_bitsets_[L_prime[0]];

    for (int j = 1; j < L_prime.size(); j++) {
      R_bs_next &= L_bitsets_[L_prime[j]];
    }

    if (((R_bs_next & (~R_bs)) << (N - C_rev.back())).count() != 0)
      continue;
    if (((R_bs_next & (~R_bs)) << (N - C_rev[i])).count() == 0) {
      if (maximal_nodes_ % 100000000 == 0) {
        printf("%llu %lfs (%d/%d)\n", (long long)maximal_nodes_,
               get_cur_time() - start_time_, cur_v_, R_size_bs_);
      }
      maximal_nodes_++;
      biclique_find_bs(L_prime, R_bs_next, C_prime_rev);
    }
    C_prime_rev.emplace_back(C_rev[i]);
  }
}

ComMBEFinder2::ComMBEFinder2(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {}

void ComMBEFinder2::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);
  r_clock_ = 0;

  graph_ = new BiGraph(*graph_);
  graph_->Reorder(RInc);

  global_2d_buf_ = std::vector<VertexSet>(graph_->GetRSize(), VertexSet());

  for (int v = 0; v < graph_->GetRSize(); v++) {
    if (graph_->NeighborsR(v).size() < min_l_size_)
      continue;
    if (v > 0 && graph_->NeighborsR(v - 1) == graph_->NeighborsR(v))
      continue;
    processing_nodes_++;

    VertexSet &L = graph_->NeighborsR(v);
    VertexSet R;
    std::vector<Node *> C;
    std::map<int, std::vector<int>> c_map;

    for (int l : L) {
      auto &l_neighbors = graph_->NeighborsL(l);
      for (int i = l_neighbors.size() - 1; l_neighbors[i] > v; i--) {
        c_map[l_neighbors[i]].emplace_back(l);
      }
    }
    R.emplace_back(v);
    for (auto &c_node : c_map) {

      if (c_node.second.size() == L.size())
        R.emplace_back(c_node.first);
      else if (c_node.second.size() >= min_l_size) {
        C.emplace_back(new Node(c_node.first, std::move(c_node.second)));
      }
    }

    if (R.size() >= min_r_size_) {
      maximal_nodes_++;
      maximum_biclique_.CompareAndSet(L, R, graph_);
    }
    double start = get_cur_time();
    biclique_find(L, R, C);
    r_clock_ += get_cur_time() - start;

    for (auto ptr : C)
      delete ptr;
  }

  printf("%lfs\n", r_clock_);
  finish();
}

void ComMBEFinder2::biclique_find(VertexSet &L, VertexSet &R,
                                  std::vector<Node *> &C) {
  // std::vector<VertexSet> global_2d_buf_(C.size(), std::vector<int>());
  std::vector<Node *> create_node_ptrs;
  VertexSet Q;

  // suppose that C[i].L is OK.
  for (int i = 0; i < C.size(); i++) {

    if (!C[i]->R.empty() || C[i]->L.size() < min_l_size_)
      continue;
    processing_nodes_++;

    int last_l = -1;
    for (int j = 0; j < C[i]->L.size(); j++) { // compute
      int l = C[i]->L[j];
      auto &l_neighbors = graph_->NeighborsL(l);
      int li = 0;

      if (j == 0) { // generate vertex set Q
        int ri = 0;
        while (l_neighbors[li] < C[i]->vc) {
          if (ri >= R.size() || l_neighbors[li] < R[ri]) {
            Q.emplace_back(l_neighbors[li]);
            li++;
          } else {
            li++;
            ri++;
          }
        }
      } else if (!Q.empty()) { // maximality check
        int q_valid = 0;
        for (int qi = 0; li < l_neighbors.size() && qi < Q.size();) {
          if (l_neighbors[li] > Q[qi])
            qi++;
          else if (l_neighbors[li] < Q[qi])
            li++;
          else {
            Q[q_valid++] = Q[qi];
            qi++;
            li++;
          }
        }
        if (q_valid == 0)
          last_l = l;
        Q.resize(q_valid);
      }

      for (int ci = i + 1; li < l_neighbors.size() && ci < C.size();) {
        if (!C[ci]->R.empty() || l_neighbors[li] > C[ci]->vc) {
          ci++;
        } else if (l_neighbors[li] < C[ci]->vc) {
          li++;
        } else {
          global_2d_buf_[ci].emplace_back(l);
          li++;
          ci++;
        }
      }
    }

    if (Q.empty()) {
      std::vector<Node *> C_prime;
      int ri = 0;
      while (ri < R.size() && R[ri] < C[i]->vc) {
        C[i]->R.emplace_back(R[ri++]);
      }
      C[i]->R.emplace_back(C[i]->vc);

      for (int ci = i + 1; ci < C.size(); ci++) {
        if (!C[ci]->R.empty())
          continue;
        while (ri < R.size() && R[ri] < C[ci]->vc) {
          C[i]->R.emplace_back(R[ri++]);
        }

        if (global_2d_buf_[ci].size() == C[i]->L.size()) {
          C[i]->R.emplace_back(C[ci]->vc);
          if (global_2d_buf_[ci].size() == C[ci]->L.size())
            C[ci]->R.emplace_back(-1);
        } else if (global_2d_buf_[ci].size() >= min_l_size_ &&
                   global_2d_buf_[ci].back() >= last_l) {
          if (global_2d_buf_[ci].size() == C[ci]->L.size())
            C_prime.emplace_back(C[ci]);
          else {
            Node *node = new Node(C[ci]->vc, global_2d_buf_[ci]);
            create_node_ptrs.emplace_back(node);
            C_prime.emplace_back(node);
          }
        }
        global_2d_buf_[ci].clear();
      }
      while (ri < R.size()) {
        C[i]->R.emplace_back(R[ri]);
        ri++;
      }

      if (C[i]->R.size() >= min_r_size_) {
        maximal_nodes_++;
        maximum_biclique_.CompareAndSet(C[i]->L, C[i]->R, graph_);
      }
      biclique_find(C[i]->L, C[i]->R, C_prime);
    } else {
      Q.clear();
      C[i]->R.emplace_back(-1);
      for (int ci = i + 1; ci < C.size(); ci++) {
        if (global_2d_buf_[ci].size() == C[ci]->L.size())
          C[ci]->R.emplace_back(-1);
        global_2d_buf_[ci].clear();
      }
    }
  }

  for (auto ptr : create_node_ptrs)
    delete ptr;
}

ComMBEFinder3::ComMBEFinder3(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {}

void ComMBEFinder3::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);
  graph_ = new BiGraph(*graph_);
  graph_->Reorder(RInc);
  global_2d_buf_.resize(graph_->GetRSize());
  idx_buf_ = std::move(std::vector<int>(graph_->GetRSize(), -1));

  for (int v = 0; v < graph_->GetRSize(); v++) {
    if (graph_->NeighborsR(v).size() < min_l_size_)
      continue;
    if (v > 0 && graph_->NeighborsR(v - 1) == graph_->NeighborsR(v))
      continue;
    processing_nodes_++;

    if (graph_->NeighborsR(v).size() > N) {

      VertexSet &L = graph_->NeighborsR(v);
      VertexSet R(1, v);
      std::vector<std::pair<int, VertexSet>> C;

      for (int l : L) {
        auto &l_neighbors = graph_->NeighborsL(l);
        for (int i = l_neighbors.size() - 1; l_neighbors[i] > v; i--) {
          int r = l_neighbors[i];
          if (idx_buf_[r] < 0) {
            idx_buf_[r] = C.size();
            C.emplace_back(
                std::make_pair(r, std::move(std::vector<int>(1, l))));
          } else {
            C[idx_buf_[r]].second.emplace_back(l);
          }
        }
      }
      std::sort(C.begin(), C.end(),
                [=](std::pair<int, VertexSet> &c0,
                    std::pair<int, VertexSet> &c1) -> bool {
                  return c0.first < c1.first;
                });

      int c_valid = 0;
      for (int i = 0; i < C.size(); i++) {
        idx_buf_[C[i].first] = -1;
        if (C[i].second.size() == L.size()) {
          R.emplace_back(C[i].first);
        } else if (C[i].second.size() >= min_l_size_) {
          if (c_valid != i)
            C[c_valid] = std::move(C[i]);
          c_valid++;
        }
      }
      C.resize(c_valid);

      if (R.size() >= min_r_size_) {
        maximal_nodes_++;
        maximum_biclique_.CompareAndSet(L, R, graph_);
      }

      biclique_find(L, R, C);
    } else {
      bitset_t L((~0LLU) >> (N - graph_->NeighborsR(v).size()));
      VertexSet R;
      std::vector<std::pair<int, bitset_t>> Q_C;
      int c_start;

      for (int li = 0; li < graph_->NeighborsR(v).size(); li++) {
        auto &l_neighbors = graph_->NeighborsL(graph_->NeighborsR(v)[li]);
        for (int r : l_neighbors) {
          if (idx_buf_[r] < 0) {
            idx_buf_[r] = Q_C.size();
            Q_C.emplace_back(std::make_pair(r, bitset_t(1LLU << li)));
          } else {
            Q_C[idx_buf_[r]].second.set(li);
          }
        }
      }
      std::sort(Q_C.begin(), Q_C.end(),
                [=](std::pair<int, bitset_t> &c0, std::pair<int, bitset_t> &c1)
                    -> bool { return c0.first < c1.first; });
      int Q_C_valid = 0;
      for (int i = 0; i < Q_C.size(); i++) {
        idx_buf_[Q_C[i].first] = -1;

        if (Q_C[i].second == L) {
          if (Q_C[i].first == v)
            c_start = Q_C_valid;
          R.emplace_back(Q_C[i].first);
        } else if (Q_C[i].second.count() >= min_l_size_) {
          if (Q_C_valid != i)
            Q_C[Q_C_valid] = std::move(Q_C[i]);
          Q_C_valid++;
        }
      }
      Q_C.resize(Q_C_valid);

      if (R.size() >= min_r_size_) {
        maximal_nodes_++;
        maximum_biclique_.CompareAndSet(L.count(), R.size());
      }
      biclique_find(L, R, Q_C, c_start);
    }
  }

  finish();
}

void ComMBEFinder3::biclique_find(VertexSet &L, VertexSet &R,
                                  std::vector<std::pair<int, VertexSet>> &C) {
  VertexSet Q;
  std::vector<int> id_array(C.size(), 0);
  for (int i = 0; i < C.size(); i++) {
    if (C[i].first < 0)
      continue;
    processing_nodes_++;

    if (C[i].second.size() > N) {

      VertexSet L_prime = std::move(C[i].second);
      VertexSet R_prime;
      std::vector<std::pair<int, VertexSet>> C_prime;
      int last_l = -1;

      for (int l : L_prime) {
        auto &l_neighbors = graph_->NeighborsL(l);
        int li = 0;

        if (L_prime[0] == l) { // initialize vertex set Q
          int ri = 0;
          while (l_neighbors[li] < C[i].first) {
            if (ri >= R.size() || l_neighbors[li] < R[ri]) {
              Q.emplace_back(l_neighbors[li]);
              li++;
            } else {
              li++;
              ri++;
            }
          }
        } else if (!Q.empty()) { // update vertex set Q
          int q_valid = 0;
          for (int qi = 0; li < l_neighbors.size() && qi < Q.size();) {
            if (l_neighbors[li] > Q[qi])
              qi++;
            else if (l_neighbors[li] < Q[qi])
              li++;
            else {
              Q[q_valid++] = Q[qi];
              qi++;
              li++;
            }
          }
          if (q_valid == 0)
            last_l = l;
          Q.resize(q_valid);
        }

        for (int ci = i + 1;
             li < l_neighbors.size() && ci < C.size();) { // update candidates
          if (l_neighbors[li] > C[ci].first) {
            ci++;
          } else if (l_neighbors[li] < C[ci].first) {
            li++;
          } else {
            // global_2d_buf_[ci].emplace_back(l);
            // id_array[ci] = -1;

            if (id_array[ci] < 0)
              global_2d_buf_[ci].emplace_back(l);
            else if (C[ci].second[id_array[ci]] == l)
              id_array[ci]++;
            else {
              for (int ii = 0; ii < id_array[ci]; ii++)
                global_2d_buf_[ci].emplace_back(C[ci].second[ii]);
              global_2d_buf_[ci].emplace_back(l);
              id_array[ci] = -1;
            }

            li++;
            ci++;
          }
        }
      }

      if (Q.empty()) {
        int ri = 0;
        while (ri < R.size() && R[ri] < C[i].first)
          R_prime.emplace_back(R[ri++]);
        R_prime.emplace_back(C[i].first);

        for (int ci = i + 1; ci < C.size(); ci++) {
          if (C[ci].first < 0)
            continue;
          while (ri < R.size() && R[ri] < C[ci].first)
            R_prime.emplace_back(R[ri++]);

          if (id_array[ci] == L_prime.size() ||
              global_2d_buf_[ci].size() == L_prime.size()) {
            R_prime.emplace_back(C[ci].first);
            if (C[ci].second.size() == L_prime.size()) {
              C[ci].first = -1;
              C[ci].second.clear();
            }
          } else if (C[ci].second.size() == id_array[ci]) {
            if (C[ci].second.back() >= last_l)
              C_prime.emplace_back(std::move(C[ci]));
            else
              C[ci].second.clear();
            C[ci].first = -1;
          } else if (id_array[ci] >= min_l_size_ &&
                     C[ci].second[id_array[ci] - 1] >= last_l) {
            VertexSet array(C[ci].second.begin(),
                            C[ci].second.begin() + id_array[ci]);
            C_prime.emplace_back(std::make_pair(C[ci].first, std::move(array)));
          } else if (global_2d_buf_[ci].size() >= min_l_size_ &&
                     global_2d_buf_[ci].back() >= last_l) {
            C_prime.emplace_back(C[ci].first, std::move(global_2d_buf_[ci]));
          }
          id_array[ci] = 0;
          global_2d_buf_[ci].clear();
        }

        while (ri < R.size()) {
          R_prime.emplace_back(R[ri]);
          ri++;
        }

        if (R_prime.size() >= min_r_size_) {
          maximal_nodes_++;
          maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
        }
        biclique_find(L_prime, R_prime, C_prime);
      } else {
        Q.clear();
        C[i].first = -1;
        C[i].second.clear();
        for (int ci = i + 1; ci < C.size(); ci++) {
          if (id_array[ci] == C[ci].second.size()) {
            C[ci].first = -1;
            C[ci].second.clear();
          }
          id_array[ci] = 0;
          global_2d_buf_[ci].clear();
        }
      }
    }
    else{
      bitset_t L_prime((~0LLU) >> (N - C[i].second.size()));
      VertexSet R_prime;
      std::vector<std::pair<int, bitset_t>> Q_C_prime;
      std::vector<bitset_t> C_bs(C.size(), 0);

      int c_start;

      for (int li = 0; li < C[i].second.size(); li++) { // compute local neighbors
        int l = C[i].second[li];
        auto &l_neighbors = graph_->NeighborsL(l);
        int ci = i + 1;

        for(int r: l_neighbors){
          if (r < C[i].first) {
            if (idx_buf_[r] < 0) {
              idx_buf_[r] = Q_C_prime.size();
              Q_C_prime.emplace_back(std::make_pair(r, bitset_t(1LLU << li)));
            } else {
              Q_C_prime[idx_buf_[r]].second.set(li);
            }
          } else {
            while (ci < C.size() && r > C[ci].first) ci++;
            if(ci >= C.size()) break;
            if (r == C[ci].first) 
              C_bs[ci].set(li);
          }
        }
      }

      int q_c_valid = 0;
      int r_size_before_vc = 0;

      for (int ci = 0; ci < Q_C_prime.size(); ci++) {
        idx_buf_[Q_C_prime[ci].first] = -1;
        if (Q_C_prime[ci].second == L_prime)
          r_size_before_vc++;
        else if (Q_C_prime[ci].second.count() >= min_l_size_) {
          if (q_c_valid != ci)
            Q_C_prime[q_c_valid] = std::move(Q_C_prime[ci]);
          q_c_valid++;
        }
      }

      if (R.size() >= r_size_before_vc &&
          R[r_size_before_vc -1] < C[i].first) { // maximal
        Q_C_prime.resize(q_c_valid);
        R_prime = R;
        R_prime.emplace_back(C[i].first);
        c_start = Q_C_prime.size();
        for (int ci = i + 1; ci < C.size(); ci++) {
          if (C_bs[ci] == L_prime)
            R_prime.emplace_back(C[ci].first);
          else if (C_bs[ci].count() >= min_l_size_)
            Q_C_prime.emplace_back(std::make_pair(C[ci].first, C_bs[ci]));
        }
        if (R_prime.size() >= min_r_size_) {
          maximal_nodes_++;
          maximum_biclique_.CompareAndSet(L_prime.count(), R_prime.size());
        }
        biclique_find(L_prime, R_prime, Q_C_prime, c_start);
      }
    }
  }
}

void ComMBEFinder3::biclique_find(bitset_t &L, VertexSet &R,
                                  std::vector<std::pair<int, bitset_t>> &Q_C,
                                  int c_start) {
  while (c_start < Q_C.size()) {

    processing_nodes_++;
    bitset_t L_prime = Q_C[c_start].second;
    VertexSet R_prime;
    std::vector<std::pair<int, bitset_t>> Q_C_prime;
    int c_start_prime;

    bool is_maximal = true;

    for (int i = 0; i < c_start; i++) {
      bitset_t bs = L_prime & Q_C[i].second;
      if (bs == L_prime) {
        if (L_prime != Q_C[i].second) {
          L_prime = Q_C[i].second;
          Q_C[c_start].second = L_prime;
        }
        is_maximal = false;
      } else if (bs.count() >= min_l_size_) {
        if (is_maximal)
          Q_C_prime.emplace_back(std::make_pair(Q_C[i].first, bs));
      }
      if (Q_C[i].second == bs)
        Q_C[i].second = 0LLU;
    }

    if (is_maximal) {
      R_prime = R;
      R_prime.emplace_back(Q_C[c_start].first);
      c_start_prime = Q_C_prime.size();
    }

    for (int i = c_start + 1; i < Q_C.size(); i++) {
      bitset_t bs = L_prime & Q_C[i].second;
      if (is_maximal && bs == L_prime)
        R_prime.emplace_back(Q_C[i].first);
      else if (is_maximal && bs.count() >= min_l_size_)
        Q_C_prime.emplace_back(std::make_pair(Q_C[i].first, bs));
      if (Q_C[i].second == bs)
        Q_C[i].second = 0LLU;
    }

    if (is_maximal) {
      if (R_prime.size() >= min_r_size_) {
        maximal_nodes_++;
        maximum_biclique_.CompareAndSet(L_prime.count(), R_prime.size());
      }
      biclique_find(L_prime, R_prime, Q_C_prime, c_start_prime);
    }
    int q_c_valid = 0;

    int c_start_next;
    for (int i = 0; i < Q_C.size(); i++) {
      if (Q_C[i].second != 0) {
        if (q_c_valid != i)
          Q_C[q_c_valid] = std::move(Q_C[i]);
        q_c_valid++;
      }
      if (i == c_start)
        c_start_next = q_c_valid;
    }
    Q_C.resize(q_c_valid);
    c_start = c_start_next;
  }
}
