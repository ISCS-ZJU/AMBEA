#include "BitMBE.h"
#define BIT_BOUND 32

BitMBEFinder::BitMBEFinder(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {}

void BitMBEFinder::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);
  bitset_time_ = 0;
  graph_->Reorder(RInc);
  for (int v = 0; v < graph_->GetRSize(); v++) {
    if (graph_->NeighborsR(v).size() < min_l_size_)
      continue;
    if (v > 0 && graph_->NeighborsR(v - 1) == graph_->NeighborsR(v))
      continue;
    processing_nodes_++;
    VertexSet L, R;
    std::vector<std::pair<int, int>> P_Q;
    int p_start = graph_->GetRSize();

    L = graph_->NeighborsR(v);

    if (L.size() > BIT_BOUND) {
      std::map<int, int> c_map;
      for (int w : L) {
        for (int y : graph_->NeighborsL(w)) {
          c_map[y]++;
        }
      }
      for (auto c_node : c_map) {
        if (c_node.second == L.size()) {
          R.emplace_back(c_node.first);
        } else if (c_node.second >= min_l_size_) {
          if (c_node.first > v && p_start > P_Q.size())
            p_start = P_Q.size();
          P_Q.emplace_back(c_node);
        }
      }

      if (R.size() >= min_r_size_) {
        maximal_nodes_++;
        maximum_biclique_.CompareAndSet(L, R, graph_);
      }

      biclique_find(L, R, P_Q, p_start);
    } else {
      double start = get_cur_time();
      BITSET_T L_bs = 0xffffffff >> (32 - L.size());
      std::map<int, int> c_bs_map;
      for (int i = 0; i < L.size(); i++) {
        int l = L[i];
        for (int w : graph_->NeighborsL(l)) {
          c_bs_map[w] |= (1 << i);
        }
      }

      std::vector<BITSET_T> P_Q_bs;
      int p_start;
      for (auto &p : c_bs_map) {
        if (p.first == v)
          p_start = P_Q_bs.size();
        else if (p.second != L_bs)
          P_Q_bs.emplace_back(p.second);
      }
      maximal_nodes_++;
      biclique_find(L_bs, P_Q_bs, p_start);
      bitset_time_ += get_cur_time() - start;
      // biclique_find_2(0, 0, L.size(), P_Q_bs, p_start);
    }
  }
  finish();
  printf("bitset time:%lfs\n", bitset_time_);
}

void BitMBEFinder::biclique_find(const VertexSet &L, const VertexSet &R,
                                 std::vector<std::pair<int, int>> &P_Q,
                                 int p_start) {

  std::vector<int> LN(P_Q.size(), 0);

  for (int pid = p_start; pid < P_Q.size(); pid++) {
    int v = P_Q[pid].first;
    if (v < 0)
      continue;
    processing_nodes_++;
    VertexSet L_prime = seq_intersect(L, graph_->NeighborsR(v));
    VertexSet R_prime;
    int p_start_prime = P_Q.size();

    if (L_prime.size() > BIT_BOUND) {
      std::vector<std::pair<int, int>> P_Q_prime;

      for (int l : L_prime) { // compute local neighbors
        auto &l_neighbors = graph_->NeighborsL(l);
        for (int i = 0, j = 0; i < P_Q.size() && j < l_neighbors.size();) {
          if (P_Q[i].first < l_neighbors[j])
            i++;
          else if (P_Q[i].first > l_neighbors[j])
            j++;
          else {
            LN[i++]++;
            j++;
          }
        }
      }

      bool is_maximal = true;
      for (int i = 0; i < pid && is_maximal; i++) { // maximality check
        if (LN[i] == L_prime.size())
          is_maximal = false;
      }

      if (is_maximal)
        R_prime = R;
      for (int i = 0; i < P_Q.size(); i++) {
        if (is_maximal) {
          if (LN[i] == L_prime.size())
            R_prime.emplace_back(P_Q[i].first);
          else if (LN[i] >= min_l_size_ && P_Q[i].first >= 0) {
            if (i > pid && p_start_prime > P_Q_prime.size())
              p_start_prime = P_Q_prime.size();
            P_Q_prime.emplace_back(std::make_pair(P_Q[i].first, LN[i]));
          }
        }
        if (i != pid && LN[i] == P_Q[i].second)
          P_Q[i].first = -1;
        LN[i] = 0;
      }

      if (is_maximal) {
        if (R_prime.size() >= min_r_size_) {
          maximal_nodes_++;
          maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
        }
        biclique_find(L_prime, R_prime, P_Q_prime, p_start_prime);
      }
    } else if(L_prime.size() >= min_l_size_){
      double start = get_cur_time();
      BITSET_T L_bs = 0xffffffff >> (32 - L_prime.size());
      for (int id = 0; id < L_prime.size(); id++) { // compute local neighbors
        int l = L_prime[id];
        auto &l_neighbors = graph_->NeighborsL(l);
        for (int i = 0, j = 0; i < P_Q.size() && j < l_neighbors.size();) {
          if (P_Q[i].first < l_neighbors[j])
            i++;
          else if (P_Q[i].first > l_neighbors[j])
            j++;
          else {
            LN[i++] |= (1 << id);
            j++;
          }
        }
      }
      bool is_maximal = true;
      for (int i = 0; i < pid && is_maximal; i++) { // maximality check
        if (LN[i] == L_bs)
          is_maximal = false;
      }

      if (is_maximal) {
        R_prime = R;
        std::vector<BITSET_T> P_Q_bs;
        for (int i = 0; i < P_Q.size(); i++) {
          if (LN[i] == L_bs)
            R_prime.emplace_back(P_Q[i].first);
          else if (bit_count(LN[i]) >= min_l_size_) {
            if (i > pid && p_start_prime > P_Q_bs.size())
              p_start_prime = P_Q_bs.size();
            P_Q_bs.emplace_back(LN[i]);
          }
        }
        if (R_prime.size() >= min_r_size_) {
          maximal_nodes_++;
          maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
        }
        biclique_find(L_bs, P_Q_bs, p_start_prime);
      }

      for (int i = 0; i < P_Q.size(); i++) {
        if (i != pid && bit_count(LN[i]) == P_Q[i].second) // prune technique
          P_Q[i].first = -1;
        LN[i] = 0;
      }
      bitset_time_ += get_cur_time() - start;
    }
  }
}

void BitMBEFinder::biclique_find(BITSET_T L_bs, std::vector<BITSET_T> P_Q_bs,
                                 int p_start) {
  for (int pid = p_start; pid < P_Q_bs.size(); pid++) {
    BITSET_T cur_L_bs = L_bs & P_Q_bs[pid];
    if (cur_L_bs == 0)
      continue;

    bool is_maximal = true;
    int p_start_prime = P_Q_bs.size();
    std::vector<BITSET_T> P_Q_bs_prime;
    processing_nodes_++;

    for (int i = 0; i < P_Q_bs.size() && is_maximal; i++) {
      if ((P_Q_bs[i] & cur_L_bs) == cur_L_bs) {
        if (i < pid)
          is_maximal = false;
      } else if ((P_Q_bs[i] & cur_L_bs) != 0) {
        if (i > pid && p_start_prime > P_Q_bs_prime.size())
          p_start_prime = P_Q_bs_prime.size();
        P_Q_bs_prime.emplace_back(P_Q_bs[i] & cur_L_bs);
      }
    }
    if (is_maximal) {
      biclique_find(cur_L_bs, P_Q_bs_prime, p_start_prime);
      maximal_nodes_++;
    }
    for (int i = 0; i < P_Q_bs.size(); i++) { // prune technique
      if ((P_Q_bs[i] & cur_L_bs) == P_Q_bs[i] && i != pid)
        P_Q_bs[i] = 0;
    }
  }
}

void BitMBEFinder::biclique_find_2(BITSET_T L_bs, int L_start, int max_L,
                                   std::vector<BITSET_T> P_Q_bs, int p_start) {
  for (int i = L_start; i < max_L; i++){
    if(L_bs & (1 << i)) continue;
    BITSET_T L_bs_current = L_bs | (1 << i);
    BITSET_T L_bs_next = 0xffffffff;
    bool is_maximal = true;
    int cnt = 0;
    for (int j = 0; j < P_Q_bs.size(); j++) {
      if ((P_Q_bs[j] & L_bs_current) == L_bs_current){
        if (j < p_start) 
          is_maximal = false;
        else
          cnt++;
        L_bs_next &= P_Q_bs[j];
      }
    }
    is_maximal = is_maximal && cnt > 0;
    if (((L_bs_next - L_bs_current) & ((1 << i) - 1)) == 0) {
      if (is_maximal) maximal_nodes_++;
      biclique_find_2(L_bs_next, i + 1, max_L, P_Q_bs, p_start);
    }
  }
}
