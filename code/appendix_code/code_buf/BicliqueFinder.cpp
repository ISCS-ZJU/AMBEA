#include "BicliqueFinder.h"

#include "BicliqueFinderFast.h"
#include <stack>
#include <unordered_map>


BicliqueFinder::BicliqueFinder(BiGraph *graph_in, const char *name) {
  finder_name_ = new char[strlen(name) + 1];
  start_time_ = 0;
  strcpy(finder_name_, name);
  graph_ = graph_in;
  processing_nodes_ = 0;
  maximal_nodes_ = 0;
  min_l_size_ = 1;
  min_r_size_ = 1;
  exe_time_ = 0;
  SetOpCounterReset();
}

void BicliqueFinder::PrintResult(char *fn) {
  FILE *fp = (fn == nullptr || strlen(fn) == 0) ? stdout : fopen(fn, "a+");
  if (fn != nullptr)
    fseek(fp, 0, SEEK_END);
  fprintf(fp, "Finder name: %s\n", finder_name_);
  fprintf(fp, "Total processing time: %lf seconds\n", exe_time_);
  fprintf(fp, "maximal nodes/processing nodes : %lld/%lld\n",
          maximal_nodes_.load(), processing_nodes_.load());
  fprintf(fp, "Max level: %d\n", max_level_);
  fprintf(fp, "min_l: %d\tmin_r: %d\n", min_l_size_, min_r_size_);
  maximum_biclique_.Print(fp);
  // long long int set_ops = GetSetOpCounter();
  // fprintf(fp, "Total set operations : %lld\n", set_ops);
  fprintf(fp, "\n");
  if (fn != NULL)
    fclose(fp);
}

MbeaFinder::MbeaFinder(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {}

void MbeaFinder::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);

  std::vector<int> L, R, P, Q;
  for (int i = 0; i < graph_->GetLSize(); i++)
    if (graph_->NeighborsL(i).size() >= min_r_size)
      L.emplace_back(i);
  for (int i = graph_->GetRSize() - 1; i >= 0; i--)
    if (graph_->NeighborsR(i).size() >= min_l_size)
      P.emplace_back(i);
  biclique_find(L, R, P, Q);
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

      if (!P_prime.empty() && (R_prime.size() + P_prime.size() >= min_r_size_))
        biclique_find(L_prime, R_prime, P_prime, Q_prime);
      if (R_prime.size() >= min_r_size_) {
        maximal_nodes_++;
        maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
      }
    }
    Q.emplace_back(x);
  }
}

MbeaAdvFinder::MbeaAdvFinder(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {}

void MbeaAdvFinder::Execute(int min_l_size, int min_r_size) {
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
      }
      for (int i = 0; i < P.size() / 2; i++) {
        std::swap(P[i], P[P.size() - 1 - i]);
      }
      biclique_find(L, R, P, Q);
    }
  }
  finish();
}

void MbeaAdvFinder::biclique_find(std::vector<int> L, std::vector<int> R,
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

      if (!P_prime.empty() && (R_prime.size() + P_prime.size() >= min_r_size_))
        biclique_find(L_prime, R_prime, P_prime, Q_prime);
      if (R_prime.size() >= min_r_size_) {
        maximal_nodes_++;
        maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
      }
    }
    Q.emplace_back(x);
  }
}

ImbeaFinder::ImbeaFinder(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {}

void ImbeaFinder::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);

  std::vector<int> L, R, P, Q;
  for (int i = 0; i < graph_->GetLSize(); i++)
    if (graph_->NeighborsL(i).size() >= min_r_size)
      L.emplace_back(i);
  for (int i = graph_->GetRSize() - 1; i >= 0; i--)
    if (graph_->NeighborsR(i).size() >= min_l_size)
      P.emplace_back(i);
  std::sort(P.begin(), P.end(), [=](int p0, int p1) -> bool {
    return graph_->NeighborsR(p0).size() > graph_->NeighborsR(p1).size();
  });

  biclique_find(L, R, P, Q);
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

        biclique_find(L_prime, R_prime, P_prime, Q_prime);
      }

      if (R_prime.size() >= min_r_size_) {
        maximal_nodes_++;
        maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
      }
      P = std::move(seq_except(P, C));
      Q.insert(Q.end(), C.begin(), C.end());
    }
  }
}

ImbeaAdvFinder::ImbeaAdvFinder(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {}

void ImbeaAdvFinder::Execute(int min_l_size, int min_r_size) {
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
      }
      std::sort(P_pairs.begin(), P_pairs.end(),
                [&](std::pair<int, int> x0, std::pair<int, int> x1) -> bool {
                  return x0.second > x1.second ||
                         (x0.second == x1.second && x0.first > x1.first);
                });
      for (auto p : P_pairs)
        P.emplace_back(p.first);
      biclique_find(L, R, P, Q);
    }
  }

  finish();
}

void ImbeaAdvFinder::biclique_find(std::vector<int> L, std::vector<int> R,
                                   std::vector<int> P, std::vector<int> Q,
                                   int level) {
  max_level_ = std::max(level, max_level_);
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

        biclique_find(L_prime, R_prime, P_prime, Q_prime, level + 1);
      }

      if (R_prime.size() >= min_r_size_) {
        maximal_nodes_++;
        maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
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

  MineLMBC(X, GamaX, tailX);
  finish();
}

void BicliqueFinder::MineLMBC(std::vector<int> X, std::vector<int> GamaX,
                              std::vector<int> tailX) {
  std::vector<std::pair<int, int>> tailx_with_nc;
  for (int v : tailX) {
    int Nc = seq_intersect_cnt(GamaX, graph_->NeighborsR(v));
    if (Nc >= min_l_size_)
      tailx_with_nc.emplace_back(std::make_pair(v, Nc));
  }
  if (X.size() + tailx_with_nc.size() < min_r_size_)
    return;
  std::sort(tailx_with_nc.begin(), tailx_with_nc.end(),
            [&](std::pair<int, int> x0, std::pair<int, int> x1) -> bool {
              return x0.second > x1.second ||
                     (x0.second == x1.second && x0.first > x1.first);
            });
  tailX.clear();
  tailX.resize(tailx_with_nc.size());
  for (int i = 0; i < tailx_with_nc.size(); i++)
    tailX[i] = tailx_with_nc[i].first;
  while (!tailX.empty()) {
    int v = tailX.back();
    std::vector<int> ordered_tailX = tailX;
    std::sort(ordered_tailX.begin(), ordered_tailX.end());
    if (X.size() + tailX.size() >= min_r_size_) {
      auto Nc = seq_intersect(GamaX, graph_->NeighborsR(v));
      std::vector<int> Y = graph_->NeighborsL(Nc);
      processing_nodes_++;
      auto Y_minus_X = seq_except(Y, X);
      if (seq_intersect_cnt(Y_minus_X, ordered_tailX) == Y_minus_X.size()) {
        if (Y.size() >= min_r_size_) {
          maximal_nodes_++;
          maximum_biclique_.CompareAndSet(Nc, Y, graph_);
        }
        MineLMBC(Y, Nc, seq_except(ordered_tailX, Y));
      }
    }
    tailX.pop_back();
  }
}

void BicliqueFinder::setup(int min_l_size, int min_r_size) {
  start_time_ = get_cur_time();
  processing_nodes_ = 0;
  maximal_nodes_ = 0;
  min_l_size_ = min_l_size;
  min_r_size_ = min_r_size;
  maximum_biclique_.Reset();
  max_level_ = 0;
}

void BicliqueFinder::finish() { exe_time_ = get_cur_time() - start_time_; }

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
      MineLMBC(X, GamaX, tailX);

      if (min_l_size_ == 1 && graph_->NeighborsL(GamaX).size() == 1) {
        maximal_nodes_++;
        maximum_biclique_.CompareAndSet(GamaX, X, graph_);
      }
    }
  }
  finish();
}

MmbeaIntraFinder::MmbeaIntraFinder(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {}

void MmbeaIntraFinder::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);
  std::vector<int> L, R;
  std::vector<CExtNode> C;
  for (int i = 0; i < graph_->GetLSize(); i++)
    if (graph_->NeighborsL(i).size() >= min_r_size)
      L.emplace_back(i);
  for (int i = 0; i < graph_->GetRSize(); i++)
    if (graph_->NeighborsR(i).size() >= min_l_size) {
      C.emplace_back(CExtNode(i, graph_->NeighborsR(i).size()));
    }
  for (auto &c : C)
    biclique_find(L, R, C, c.r_id);
  finish();
}

void MmbeaIntraFinder::biclique_find(const std::vector<int> &L,
                                     const std::vector<int> &R,
                                     std::vector<CExtNode> &C, int vc) {
  std::vector<int> L_prime = seq_intersect(L, graph_->NeighborsR(vc)),
                   R_prime = R;
  std::vector<CExtNode> C_prime;
  if (L_prime.size() < min_l_size_)
    return;
  processing_nodes_++;

  std::vector<int> reordered_map;
  std::vector<PNode> P;

  Partition(L_prime, C, reordered_map, P);

  // maximality check
  int min_r_id = reordered_map[0];
  for (int i = 1; i < P[0].size; i++)
    min_r_id = std::min(min_r_id, reordered_map[i]);
  if (C[min_r_id].r_id != vc)
    return; // non-maximal
  for (int i = 0; i < P[0].size; i++) {
    auto &R_tmp = C[reordered_map[i]].r_cands;
    R_prime.insert(R_prime.end(), R_tmp.begin(), R_tmp.end());
  }

  if (R_prime.size() >= min_r_size_) {
    maximal_nodes_++;
    maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
  }

  // generate C_prime
  for (int i = 1; i < P.size(); i++) {
    int start = P[i].start;
    int end = start + P[i].size;
    int nc = P[i].nc;

    if (nc < min_l_size_)
      continue;

    C_prime.emplace_back(CExtNode());
    auto &c = C_prime.back();
    c.r_id = C[reordered_map[start]].r_id;
    c.nc = nc;

    for (int j = start; j < end; j++) {
      auto &R_tmp = C[reordered_map[j]].r_cands;
      c.r_cands.insert(c.r_cands.end(), R_tmp.begin(), R_tmp.end());
      c.r_id = std::min(c.r_id, C[reordered_map[j]].r_id);
    }
  }
  std::sort(
      C_prime.begin(), C_prime.end(),
      [&](CExtNode c0, CExtNode c1) -> bool { return c0.r_id < c1.r_id; });
  std::vector<int> c_cands;
  for (auto &c : C_prime)
    if (c.r_id > vc)
      c_cands.emplace_back(c.r_id);
  for (int new_vc : c_cands)
    biclique_find(L_prime, R_prime, C_prime, new_vc);
}

MmbeaIntraFinderV2::MmbeaIntraFinderV2(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {}

void MmbeaIntraFinderV2::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);
  std::vector<int> L, R;
  std::vector<int> C;
  for (int i = 0; i < graph_->GetLSize(); i++)
    if (graph_->NeighborsL(i).size() >= min_r_size)
      L.emplace_back(i);
  for (int i = 0; i < graph_->GetRSize(); i++)
    if (graph_->NeighborsR(i).size() >= min_l_size) {
      C.emplace_back(i);
    }
  for (int c : C) {
    // printf("%d/%d\n", maximal_nodes_, processing_nodes_);
    biclique_find(L, R, C, c);
  }
  finish();
}

void MmbeaIntraFinderV2::biclique_find(const std::vector<int> &L,
                                       const std::vector<int> &R,
                                       std::vector<int> &C, int vc) {
  std::vector<int> L_prime = seq_intersect(L, graph_->NeighborsR(vc));
  std::vector<int> R_prime;
  std::vector<int> C_prime;

  if (L_prime.size() < min_l_size_ || L_prime.size() == L.size())
    return;
  processing_nodes_++;

  R_prime = graph_->NeighborsL(L_prime);
  auto R_add = seq_except(R_prime, R);
  if (R_add[0] != vc)
    return; // non-maximal
  if (R_prime.size() >= min_r_size_) {
    maximal_nodes_++;
    maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
  }
  auto C_buffer = seq_except_upper(C, R_add, vc);
  std::vector<std::pair<int, int>> P;
  std::vector<int> reordered_map;

  if (C_buffer.size() == 0)
    return;
  Partition(L_prime, C_buffer, reordered_map, P);

  C_prime.reserve(P.size());
  for (int i = 0; i < P.size(); i++) {
    int min_id = reordered_map[P[i].first];
    for (int j = P[i].first + 1; j < P[i].first + P[i].second; j++)
      min_id = std::min(min_id, reordered_map[j]);
    C_prime.emplace_back(C_buffer[min_id]);
  }
  std::sort(C_prime.begin(), C_prime.end());
  for (int c : C_prime) {
    biclique_find(L_prime, R_prime, C_prime, c);
  }
}

inline void BicliqueFinder::Partition(const std::vector<int> &L_prime,
                                      const std::vector<CExtNode> &C,
                                      std::vector<int> &reordered_map,
                                      std::vector<PNode> &P) {
  SetOpCounterAdd(L_prime.size());
  if (reordered_map.empty()) {
    reordered_map.reserve(C.size());
    for (int i = 0; i < C.size(); i++)
      if (C[i].r_id >= 0)
        reordered_map.emplace_back(i);
  }

  P.emplace_back(PNode(0, reordered_map.size(), 0));
  for (int l : L_prime) {
    std::vector<bool> edge_exist_bitset(C.size(), false);
    auto &l_neighbors = graph_->NeighborsL(l);
    for (int i = 0, j = 0; i < l_neighbors.size() && j < C.size();) {
      if (l_neighbors[i] < C[j].r_id)
        i++;
      else if (l_neighbors[i] > C[j].r_id)
        j++;
      else {
        edge_exist_bitset[j] = true;
        i++;
        j++;
      }
    }

    int size = P.size();
    for (int i = 0; i < size; i++) {
      auto &p = P[i];
      if (p.size == 1) {
        if (edge_exist_bitset[reordered_map[p.start]])
          p.nc++;
        continue;
      }
      int l_index = p.start, r_index = p.start + p.size - 1;

      while (l_index < r_index) {
        while (l_index <= r_index && edge_exist_bitset[reordered_map[l_index]])
          l_index++;
        while (l_index <= r_index && !edge_exist_bitset[reordered_map[r_index]])
          r_index--;
        if (l_index < r_index) {
          std::swap(reordered_map[l_index], reordered_map[r_index]);
        }
      }

      if (l_index > p.start) {
        p.nc++;
        if (l_index < p.start + p.size) {
          int next_size = p.start + p.size - l_index;
          p.size = l_index - p.start;
          P.emplace_back(PNode(l_index, next_size, p.nc - 1));
        }
      }
    }
  }
}

inline void BicliqueFinder::Partition(const std::vector<int> &L_prime,
                                      const std::vector<std::pair<int, int>> &C,
                                      std::vector<int> &reordered_map,
                                      std::vector<PNode> &P) {
  SetOpCounterAdd(L_prime.size());
  if (reordered_map.empty()) {
    reordered_map.reserve(C.size());
    for (int i = 0; i < C.size(); i++)
      if (C[i].first >= 0)
        reordered_map.emplace_back(i);
  }
  P.emplace_back(PNode(0, reordered_map.size(), 0));
  for (int l : L_prime) {
    std::vector<bool> edge_exist_bitset(C.size(), false);
    auto &l_neighbors = graph_->NeighborsL(l);
    for (int i = 0, j = 0; i < l_neighbors.size() && j < C.size();) {
      if (l_neighbors[i] < C[j].first)
        i++;
      else if (l_neighbors[i] > C[j].first)
        j++;
      else {
        edge_exist_bitset[j] = true;
        i++;
        j++;
      }
    }

    int size = P.size();
    for (int i = 0; i < size; i++) {
      auto &p = P[i];
      if (p.size == 1) {
        if (edge_exist_bitset[reordered_map[p.start]])
          p.nc++;
        continue;
      }
      int l_index = p.start, r_index = p.start + p.size - 1;

      while (l_index < r_index) {
        while (l_index <= r_index && edge_exist_bitset[reordered_map[l_index]])
          l_index++;
        while (l_index <= r_index && !edge_exist_bitset[reordered_map[r_index]])
          r_index--;
        if (l_index < r_index) {
          std::swap(reordered_map[l_index], reordered_map[r_index]);
        }
      }

      if (l_index > p.start) {
        p.nc++;
        if (l_index < p.start + p.size) {
          int next_size = p.start + p.size - l_index;
          p.size = l_index - p.start;
          P.emplace_back(PNode(l_index, next_size, p.nc - 1));
        }
      }
    }
  }
}

inline void BicliqueFinder::Partition(const std::vector<int> &L_prime,
                                      const std::vector<int> &C,
                                      std::vector<int> &reordered_map,
                                      std::vector<std::pair<int, int>> &P) {
  SetOpCounterAdd(L_prime.size());
  if (reordered_map.empty()) {
    reordered_map.reserve(C.size());
    for (int i = 0; i < C.size(); i++)
      reordered_map.emplace_back(i);
  }
  P.emplace_back(std::make_pair(0, reordered_map.size()));
  for (int l : L_prime) {
    std::vector<bool> edge_exist_bitset(C.size());
    auto &l_neighbors = graph_->NeighborsL(l);
    for (int i = 0, j = 0; i < l_neighbors.size() && j < C.size();) {
      if (l_neighbors[i] < C[j])
        i++;
      else if (l_neighbors[i] > C[j])
        j++;
      else {
        edge_exist_bitset[j] = true;
        i++;
        j++;
      }
    }
    int size = P.size();
    for (int i = 0; i < size; i++) {
      if (P[i].second == 1)
        continue;
      int l_index = P[i].first, r_index = P[i].first + P[i].second - 1;
      while (l_index < r_index) {
        while (l_index <= r_index && edge_exist_bitset[reordered_map[l_index]])
          l_index++;
        while (l_index <= r_index && !edge_exist_bitset[reordered_map[r_index]])
          r_index--;
        if (l_index < r_index)
          std::swap(reordered_map[l_index], reordered_map[r_index]);
      }

      if (l_index > P[i].first && l_index < P[i].first + P[i].second) {
        int next_size = P[i].first + P[i].second - l_index;
        P[i].second = l_index - P[i].first;
        P.emplace_back(std::make_pair(l_index, next_size));
      }
    }
  }
}

MEBFinder::MEBFinder(BiGraph *graph) : orig_graph_(graph), exe_time_(0) {}

void MEBFinder::Execute(int ctrl) {
  double start = get_cur_time();
  double my_clock = 0;
  maximum_biclique_.Reset();

  BiGraph *graph = orig_graph_;
  graph->Prune2HOpt(3, 3);
  BicliqueFinder *finder;

  int lb = 256; // graph->GetRSize() / 4;
  int rb = 3;

  while (maximum_biclique_.GetEdges() == 0) {
    lb = std::max(3, lb / 2);
    graph->Prune2HOpt(lb, rb);
    if (graph->GetEdges() > 0) {
      double dd = get_cur_time();
      switch (ctrl) {
      case 0:
        finder = new ImbeaFinder(graph);
        break;
      case 1:
        finder = new MmbeaFinderFast(graph);
        break;
      case 2:
        finder = new MmbeaFinderV2(graph);
        break;
      default:
        finder = new MineLMBCFinder(graph);
        break;
      }
      finder->Execute(lb, rb);
      if (finder->maximum_biclique_.GetEdges() > maximum_biclique_.GetEdges()) {
        maximum_biclique_ = finder->maximum_biclique_;
      }
      // printf("%d %d\n", lb, rb);
      // maximum_biclique_.Print();
      // finder->PrintResult();

      delete finder;
      my_clock += get_cur_time() - dd;
    }
    graph->PopOrigGraph();
  }
  // maximum_biclique_.Print();
  // printf("init time:%lf\n", get_cur_time() - start);

  while (lb > 3) {
    int max_edges = maximum_biclique_.GetEdges();
    rb = (max_edges - 1) / (lb - 1) + 1;
    lb = std::max(lb / 2, 3);
    graph->Prune2HOpt(lb, rb);
    while (graph->GetLSize() == 0 && lb > 3) {
      graph->PopOrigGraph();
      rb = (max_edges - 1) / (lb - 1) + 1;
      lb = std::max(lb / 2, 3);
      graph->Prune2HOpt(lb, rb);
    }

    if (lb == 1)
      break;
    double dd = get_cur_time();
    switch (ctrl) {
    case 0:
      finder = new ImbeaFinder(graph);
      break;
    case 1:
      finder = new MmbeaFinderFast(graph);
      break;
    case 2:
      finder = new MmbeaFinderV2(graph);
      break;
    default:
      finder = new MineLMBCFinder(graph);
      break;
    }

    finder->Execute(lb, rb);
    if (finder->maximum_biclique_.GetEdges() > maximum_biclique_.GetEdges()) {
      maximum_biclique_ = finder->maximum_biclique_;
    }
    // printf("%d %d\n", lb, rb);
    // maximum_biclique_.Print();
    // finder->PrintResult();

    delete finder;
    my_clock += get_cur_time() - dd;
    graph->PopOrigGraph();
  }
  maximum_biclique_.Print();
  printf("my clock:%lf\n", my_clock);
  exe_time_ = get_cur_time() - start;
}

void MEBFinder::PrintResult(char *fn) {
  FILE *fp = (fn == nullptr) ? stdout : fopen(fn, "a+");
  if (fn != nullptr)
    fseek(fp, 0, SEEK_END);
  fprintf(fp, "Total processing time: %lf seconds\n", exe_time_);
  maximum_biclique_.Print(fp);
  fprintf(fp, "\n");
  if (fn != NULL)
    fclose(fp);
}

MmbeaInterFinder::MmbeaInterFinder(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {}

void MmbeaInterFinder::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);

  graph_ = new BiGraph(*graph_);
  graph_->Reorder(RInc);

  for (int v = 0; v < graph_->GetRSize(); v++) {
    processing_nodes_++;

    std::vector<int> L, R;
    std::vector<std::pair<int, int>> C;
    std::map<int, int> c_map; // r_id, nc
    L = graph_->NeighborsR(v);
    R.emplace_back(v);
    for (int w : L) {
      for (int y : graph_->NeighborsL(w)) {
        if (y != v)
          c_map[y]++;
      }
    }
    bool is_maximal = true;
    for (auto c_node : c_map) {
      if (c_node.second == L.size()) {
        if (c_node.first < v) {
          is_maximal = false;
          break;
        }
        R.emplace_back(c_node.first);
      } else if (c_node.second >= min_l_size_) {
        C.emplace_back(c_node);
      }
    }
    if (is_maximal) {
      if (R.size() >= min_r_size_) {
        maximal_nodes_++;
        maximum_biclique_.CompareAndSet(L, R, graph_);
      }

      for (auto c : C)
        if (c.first > v)
          biclique_find(L, R, C, v, c.first);
    }
  }

  finish();
}

void MmbeaInterFinder::biclique_find(const std::vector<int> &L,
                                     const std::vector<int> &R,
                                     const std::vector<std::pair<int, int>> &C,
                                     int vp, int vc) {
  std::vector<int> L_prime = seq_intersect(L, graph_->NeighborsR(vc));
  std::vector<int> R_prime = R;
  std::vector<std::pair<int, int>> C_prime;

  if (L_prime.size() < min_l_size_)
    return;
  processing_nodes_++;
  std::vector<int> L_current = L;
  std::vector<int> c_cands;

  for (auto &c : C) {
    int nc = seq_intersect_cnt(L_prime, graph_->NeighborsR(c.first));
    if (nc < min_l_size_)
      continue;
    if (nc == L_prime.size()) {
      if (c.first < vp)
        return;
      if (c.first < vc) {
        L_current = seq_intersect(L_current, graph_->NeighborsR(c.first));
        if (L_current.size() == L_prime.size())
          return;
      }
      R_prime.emplace_back(c.first);
    } else {
      C_prime.emplace_back(std::make_pair(c.first, nc));
      if (nc != c.second && c.first > vc) {
        int nc_last = c.second;
        if (L_current.size() != L.size())
          nc_last = seq_intersect_cnt(L_current, graph_->NeighborsR(c.first));
        if (nc_last != nc)
          c_cands.emplace_back(c.first);
      }
    }
  }
  if (R_prime.size() >= min_r_size_) {
    maximal_nodes_++;
    maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
  }

  for (int new_vc : c_cands)
    biclique_find(L_prime, R_prime, C_prime, vc, new_vc);
}

MmbeaInterAdvFinder::MmbeaInterAdvFinder(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {}

void MmbeaInterAdvFinder::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);

  graph_ = new BiGraph(*graph_);
  graph_->Reorder(RInc);

  for (int v = 0; v < graph_->GetRSize(); v++) {
    processing_nodes_++;

    std::vector<int> L, R;
    std::vector<std::pair<int, int>> C;
    std::map<int, int> c_map; // r_id, nc
    L = graph_->NeighborsR(v);
    R.emplace_back(v);
    for (int w : L) {
      for (int y : graph_->NeighborsL(w)) {
        if (y != v)
          c_map[y]++;
      }
    }
    bool is_maximal = true;
    for (auto c_node : c_map) {
      if (c_node.second == L.size()) {
        if (c_node.first < v) {
          is_maximal = false;
          break;
        }
        R.emplace_back(c_node.first);
      } else if (c_node.second >= min_l_size_) {
        C.emplace_back(c_node);
      }
    }
    if (is_maximal) {
      if (R.size() >= min_r_size_) {
        maximal_nodes_++;
        maximum_biclique_.CompareAndSet(L, R, graph_);
      }

      for (auto c : C)
        if (c.first > v)
          biclique_find(L, R, C, v, c.first);
    }
  }

  finish();
}

void MmbeaInterAdvFinder::biclique_find(
    const std::vector<int> &L, const std::vector<int> &R,
    const std::vector<std::pair<int, int>> &C, int vp, int vc) {
  std::vector<int> L_prime = seq_intersect(L, graph_->NeighborsR(vc));
  std::vector<int> R_prime;
  std::vector<std::pair<int, int>> C_prime;

  if (L_prime.size() < min_l_size_)
    return;
  processing_nodes_++;
  std::vector<int> L_current = L;
  std::vector<int> c_cands;

  // maximality check
  R_prime = graph_->NeighborsL(L_prime);
  auto R_add = seq_except(R_prime, R);
  auto R_add_iter = R_add.begin();
  if (*R_add_iter <= vp)
    return;
  while (*R_add_iter < vc) {
    L_current = seq_intersect(L_current, graph_->NeighborsR(*R_add_iter));
    if (L_current.size() == L_prime.size())
      return;
    R_add_iter++;
  }

  for (auto &c : C) {
    if (c.first < vc)
      continue;
    if (R_add_iter != R_add.end() && c.first == *R_add_iter) {
      R_add_iter++;
      continue;
    }

    int nc = seq_intersect_cnt(L_prime, graph_->NeighborsR(c.first));
    if (nc < min_l_size_)
      continue;

    C_prime.emplace_back(std::make_pair(c.first, nc));
    if (nc != c.second) {
      int nc_last = c.second;
      if (L_current.size() != L.size())
        nc_last = seq_intersect_cnt(L_current, graph_->NeighborsR(c.first));
      if (nc_last != nc)
        c_cands.emplace_back(c.first);
    }
  }

  if (R_prime.size() >= min_r_size_) {
    maximal_nodes_++;
    maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
  }

  for (int new_vc : c_cands)
    biclique_find(L_prime, R_prime, C_prime, vc, new_vc);
}

MmbeaFinder::MmbeaFinder(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {}

void MmbeaFinder::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);

  graph_->Reorder(RInc);

  for (int v = 0; v < graph_->GetRSize(); v++) {
    if (graph_->NeighborsR(v).size() < min_l_size_)
      continue;
    if (v > 0 && graph_->NeighborsR(v - 1) == graph_->NeighborsR(v)) {
      continue;
    }
    processing_nodes_++;

    std::vector<int> L, R;
    std::vector<std::pair<int, int>> C;
    L = graph_->NeighborsR(v);

    std::map<int, int> c_map; // r_id, nc
    for (int w : L) {
      for (int y : graph_->NeighborsL(w)) {
        if (y >= v)
          c_map[y]++;
      }
    }

    for (auto c_node : c_map) {
      if (c_node.second == L.size()) {
        R.emplace_back(c_node.first);
      } else if (c_node.second >= min_l_size_) {
        C.emplace_back(std::make_pair(c_node.first, c_node.second));
      }
    }

    if (R.size() >= min_r_size_) {
      maximal_nodes_++;
      maximum_biclique_.CompareAndSet(L, R, graph_);
    }

    if (C.size() == 0)
      continue;
    std::vector<int> reordered_map;
    std::vector<PNode> P;

    Partition(L, C, reordered_map, P);
    std::vector<std::pair<int, int>> C_prime;
    C_prime.reserve(P.size());
    for (int i = 0; i < P.size(); i++) {
      if (P[i].nc < min_l_size_)
        continue;
      int min_index = reordered_map[P[i].start];
      for (int j = P[i].start + 1; j < P[i].start + P[i].size; j++)
        min_index = std::min(min_index, reordered_map[j]);
      C_prime.emplace_back(std::move(C[min_index]));
    }
    std::sort(C_prime.begin(), C_prime.end(),
              [&](std::pair<int, int> p0, std::pair<int, int> p1) -> bool {
                return p0.first < p1.first;
              });
    C.clear();
    for (const auto &c : C_prime) {
      biclique_find(L, R, C_prime, v, c.first);
    }
  }
  finish();
}

void MmbeaFinder::biclique_find(const std::vector<int> &L,
                                const std::vector<int> &R,
                                std::vector<std::pair<int, int>> &C, int vp,
                                int vc) {
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
  struct CBufferNode {
    int group_id;
    std::pair<int, int> &c;
    CBufferNode(int group_id, std::pair<int, int> &c)
        : group_id(group_id), c(c) {}
  };

  std::vector<CBufferNode> C_buffer;

  auto r_add_iter = R_add.begin();
  for (int i = 0; i < C.size(); i++) {
    if (C[i].first < 0)
      continue;
    while (r_add_iter != R_add.end() && *r_add_iter < C[i].first)
      r_add_iter++;
    if (r_add_iter != R_add.end() && *r_add_iter == C[i].first) {
      if (*r_add_iter < vc) {
        L_current = seq_intersect(L_current, graph_->NeighborsR(*r_add_iter));
        if (L_current.size() == L_prime.size())
          return;
      }
    } else
      C_buffer.emplace_back(CBufferNode(0, C[i]));
  }

  if (R_prime.size() >= min_r_size_) {
    maximal_nodes_++;
    maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
  }

  if (C_buffer.size() == 0)
    return;
  struct Node {
    int r_id;
    int nc;
    bool is_connect;
    int group_id;
    Node(int r_id, int nc, bool is_connect, int group_id)
        : r_id(r_id), nc(nc), is_connect(is_connect), group_id(group_id) {}
  };
  std::vector<Node> groups;
  groups.emplace_back(Node(C_buffer.front().c.first, 0, false, 0));

  for (int l : L_prime) {
    auto neighbor_iter = graph_->NeighborsL(l).begin();
    for (int i = 0; i < C_buffer.size(); i++) {
      auto &group = groups[C_buffer[i].group_id];
      int current_r_id = C_buffer[i].c.first;
      int current_group_id = C_buffer[i].group_id;
      while (neighbor_iter != graph_->NeighborsL(l).end() &&
             *neighbor_iter < current_r_id)
        neighbor_iter++;
      bool is_connect = (neighbor_iter != graph_->NeighborsL(l).end()) &&
                        (*neighbor_iter == current_r_id);
      if (group.r_id == current_r_id) {
        group.group_id = current_group_id;
        group.is_connect = is_connect;
        if (is_connect)
          group.nc++;
      } else if (group.is_connect != is_connect) {
        if (group.group_id != current_group_id)
          C_buffer[i].group_id = group.group_id;
        else {
          int nc = is_connect ? group.nc + 1 : group.nc - 1;
          int group_id = groups.size();
          group.group_id = group_id;
          C_buffer[i].group_id = group_id;
          groups.emplace_back(Node(current_r_id, nc, is_connect, group_id));
        }
      }
    }
  }

  std::vector<int> cand_r_ids;

  C_prime.reserve(groups.size());
  cand_r_ids.reserve(groups.size());

  auto c_buf_iter = C_buffer.begin();
  while (c_buf_iter != C_buffer.end()) {
    auto &group = groups[c_buf_iter->group_id];
    if (c_buf_iter->c.first == group.r_id && group.nc >= min_l_size_) {
      int r_id = c_buf_iter->c.first;
      C_prime.emplace_back(std::make_pair(r_id, group.nc));
      if (r_id > vc && c_buf_iter->c.second != group.nc) {
        if (L_current.size() == L.size() ||
            seq_intersect_cnt(L_current, graph_->NeighborsR(r_id)) > group.nc)
          cand_r_ids.emplace_back(r_id);
      }
      if (r_id < vp && c_buf_iter->c.second == group.nc)
        c_buf_iter->c.first = -1;
    }
    c_buf_iter++;
  }

  for (int this_r_id : cand_r_ids) {
    biclique_find(L_prime, R_prime, C_prime, vc, this_r_id);
  }
}

MmbeaFinderV2::MmbeaFinderV2(BiGraph *graph_in, const char *name)
    : BicliqueFinder(graph_in, name) {}

void MmbeaFinderV2::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);

  graph_->Reorder(RInc);
  for (int v = 0; v < graph_->GetRSize(); v++) {
    if (graph_->NeighborsR(v).size() < min_l_size_)
      continue;
    if (v > 0 && graph_->NeighborsR(v - 1) == graph_->NeighborsR(v))
      continue;
    processing_nodes_++;

    std::vector<int> L, R;
    std::vector<CNode> C;
    L = graph_->NeighborsR(v);

    std::map<int, int> c_map; // r_id, nc
    for (int w : L) {
      for (int y : graph_->NeighborsL(w)) {
        if (y >= v)
          c_map[y]++;
      }
    }
    for (auto c_node : c_map) {
      if (c_node.second == L.size()) {
        R.emplace_back(c_node.first);
      } else if (c_node.second >= min_l_size_) {
        C.emplace_back(CNode(c_node.first, c_node.second, 1));
      }
    }
    if (R.size() >= min_r_size_) {
      maximal_nodes_++;
      maximum_biclique_.CompareAndSet(L, R, graph_);
    }
    if (C.size() == 0)
      continue;
    std::vector<int> reordered_map;
    std::vector<PNode> P;

    Partition(L, C, reordered_map, P);
    std::vector<CNode> C_prime;
    C_prime.reserve(P.size());

    for (int i = 0; i < P.size(); i++) {
      if (P[i].nc < min_l_size_)
        continue;
      int min_index = reordered_map[P[i].start];
      int size = C[min_index].size;
      for (int j = P[i].start + 1; j < P[i].start + P[i].size; j++) {
        min_index = std::min(min_index, reordered_map[j]);
        size += C[reordered_map[j]].size;
      }
      C[min_index].size = size;
      C_prime.emplace_back(std::move(C[min_index]));
    }
    std::sort(C_prime.begin(), C_prime.end(),
              [&](CNode c0, CNode c1) -> bool { return c0.r_id < c1.r_id; });
    C.clear();
    int tail_size = 0;
    for (auto &c : C_prime)
      tail_size += c.size;

    if (R.size() + tail_size >= min_r_size_)
      for (auto &c : C_prime)
        biclique_find(L, R, C_prime, v, c.r_id);
  }

  finish();
}

void MmbeaFinderV2::biclique_find(const std::vector<int> &L,
                                  const std::vector<int> &R,
                                  std::vector<CNode> &C, int vp, int vc) {
  std::vector<int> L_prime, R_prime;
  std::vector<CNode> C_prime;

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

  std::vector<int> reordered_map;
  std::vector<PNode> P;

  auto r_add_iter = R_add.begin();

  for (int i = 0; i < C.size(); i++) {
    if (C[i].r_id < 0)
      continue;
    while (r_add_iter != R_add.end() && *r_add_iter < C[i].r_id)
      r_add_iter++;
    if (r_add_iter != R_add.end() && *r_add_iter == C[i].r_id) {
      if (*r_add_iter < vc) {
        L_current = seq_intersect(L_current, graph_->NeighborsR(*r_add_iter));
        if (L_current.size() == L_prime.size())
          return;
      }
    } else if (C[i].r_id > vc)
      reordered_map.emplace_back(i);
  }

  if (R_prime.size() >= min_r_size_) {
    maximal_nodes_++;
    maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
  }

  if (reordered_map.size() == 0)
    return;
  Partition(L_prime, C, reordered_map, P);

  std::vector<int> cand_r_ids;
  C_prime.reserve(P.size());

  for (int i = 0; i < P.size(); i++) {
    if (P[i].nc < min_l_size_)
      continue;
    int min_index = reordered_map[P[i].start];
    int size = C[min_index].size;
    for (int j = P[i].start + 1; j < P[i].start + P[i].size; j++) {
      min_index = std::min(min_index, reordered_map[j]);
      size += C[reordered_map[j]].size;
    }
    int r_id = C[min_index].r_id;
    C_prime.emplace_back(CNode(r_id, P[i].nc, size));
    if (C[min_index].nc != P[i].nc) {
      if (L_current.size() == L.size() ||
          seq_intersect_cnt(L_current, graph_->NeighborsR(r_id)) > P[i].nc)
        cand_r_ids.emplace_back(r_id);
    }
  }
  int tail_size = 0;
  for (auto &c : C_prime)
    tail_size += c.size;
  if (R_prime.size() + tail_size >= min_r_size_) {
    std::sort(C_prime.begin(), C_prime.end(),
              [&](CNode c0, CNode c1) -> bool { return c0.r_id < c1.r_id; });
    for (int c : cand_r_ids)
      biclique_find(L_prime, R_prime, C_prime, vc, c);
  }
}

void MmbeaFinderV2::Partition(const std::vector<int> &L_prime,
                              const std::vector<CNode> &C,
                              std::vector<int> &reordered_map,
                              std::vector<PNode> &P) {
  if (reordered_map.empty()) {
    reordered_map.reserve(C.size());
    for (int i = 0; i < C.size(); i++)
      if (C[i].r_id >= 0)
        reordered_map.emplace_back(i);
  }
  P.emplace_back(PNode(0, reordered_map.size(), 0));
  for (int l : L_prime) {
    std::vector<bool> edge_exist_bitset(C.size(), false);
    auto &l_neighbors = graph_->NeighborsL(l);
    for (int i = 0, j = 0; i < l_neighbors.size() && j < C.size();) {
      if (l_neighbors[i] < C[j].r_id)
        i++;
      else if (l_neighbors[i] > C[j].r_id)
        j++;
      else {
        edge_exist_bitset[j] = true;
        i++;
        j++;
      }
    }

    int size = P.size();
    for (int i = 0; i < size; i++) {
      auto &p = P[i];
      if (p.size == 1) {
        if (edge_exist_bitset[reordered_map[p.start]])
          p.nc++;
        continue;
      }
      int l_index = p.start, r_index = p.start + p.size - 1;

      while (l_index < r_index) {
        while (l_index <= r_index && edge_exist_bitset[reordered_map[l_index]])
          l_index++;
        while (l_index <= r_index && !edge_exist_bitset[reordered_map[r_index]])
          r_index--;
        if (l_index < r_index) {
          std::swap(reordered_map[l_index], reordered_map[r_index]);
        }
      }

      if (l_index > p.start) {
        p.nc++;
        if (l_index < p.start + p.size) {
          int next_size = p.start + p.size - l_index;
          p.size = l_index - p.start;
          P.emplace_back(PNode(l_index, next_size, p.nc - 1));
        }
      }
    }
  }
}

LevelFinder::LevelFinder(BiGraph *graph, const char *name)
    : BicliqueFinder(graph, name) {}

void LevelFinder::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);
  recover_time_ = 0;
  graph_->Reorder(RInc);

  for (int v = 0; v < graph_->GetRSize(); v++) {
    if (graph_->NeighborsR(v).size() < min_l_size_)
      continue;
    if (v > 0 && graph_->NeighborsR(v - 1) == graph_->NeighborsR(v))
      continue;
    maximal_nodes_++;
    processing_nodes_++;
    std::map<int, int> c_map;

    L_.clear();
    R_.clear();
    C_.clear();
    for (int l : graph_->NeighborsR(v)) {
      L_.emplace_back(std::make_pair(l, 0));
      for (int r : graph_->NeighborsL(l)) {
        c_map[r]++;
      }
    }
    int cid_start = -1;

    for (auto c_node : c_map) {
      if (c_node.second == L_.size())
        R_.emplace_back(c_node.first);
      else if (c_node.second >= min_l_size) {
        if (c_node.first < v)
          cid_start = C_.size();
        C_.emplace_back(std::make_pair(c_node, 0x7fffffff));
      }
    }

    cid_stack_.clear();
    cid_stack_.emplace_back(cid_start);
    while (!cid_stack_.empty()) {
      int last_cid = cid_stack_.back();
      cid_stack_.pop_back();
      int cur_cid = FindNext(last_cid);
      if (cur_cid != -1) {
        if (Push(cur_cid)) {
          if (C_[cur_cid].first.second <= min_l_size_)
            cur_cid = C_.size();
          cid_stack_.emplace_back(cur_cid);
        } else {
          Pop();
          cid_stack_.back() = cur_cid;
        }

      } else {
        if (!cid_stack_.empty()) {
          Pop();
        }
      }
    }
  }
  printf("recover time: %lfs\n", recover_time_);
  finish();
}

int LevelFinder::FindNext(int last_cid) {
  int level = cid_stack_.size();
  for (int i = last_cid + 1; i < C_.size(); i++) {
    // no lock
    if (C_[i].second > level && C_[i].first.second >= min_l_size_) {
      return i;
    }
  }
  return -1;
}

bool LevelFinder::Push(int c_id) {
  processing_nodes_++;
  int level = cid_stack_.size() + 1;
  int cur_l_size = C_[c_id].first.second;

  // get new L
  auto &r_neighbors = graph_->NeighborsR(C_[c_id].first.first);
  for (int i = 0, j = 0; i < L_.size() && j < r_neighbors.size();) {
    if (L_[i].first < r_neighbors[j] || L_[i].second != level - 1)
      i++;
    else if (L_[i].first > r_neighbors[j])
      j++;
    else {
      L_[i].second = level;
      i++;
      j++;
    }
  }

  int last_l_size = 0;
  if (cid_stack_.empty())
    last_l_size = L_.size();
  else
    last_l_size = C_[cid_stack_.back()].first.second;

  // update new C
  if (last_l_size <= (cur_l_size << 1)) {
    // if (true) {
    for (auto &lp : L_) {
      if (lp.second == level - 1) {
        auto &l_neighbors = graph_->NeighborsL(lp.first);
        for (int i = 0, j = 0; i < C_.size() && j < l_neighbors.size();) {
          if (C_[i].first.first < l_neighbors[j])
            i++;
          else if (C_[i].first.first > l_neighbors[j])
            j++;
          else {
            C_[i].first.second--;
            i++;
            j++;
          }
        }
      }
    }
  } else {
    for (auto &cp : C_)
      cp.first.second = 0;
    for (auto &lp : L_) {
      if (lp.second == level) {
        auto &l_neighbors = graph_->NeighborsL(lp.first);
        for (int i = 0, j = 0; i < C_.size() && j < l_neighbors.size();) {
          if (C_[i].first.first < l_neighbors[j])
            i++;
          else if (C_[i].first.first > l_neighbors[j])
            j++;
          else {
            C_[i].first.second++;
            i++;
            j++;
          }
        }
      }
    }
  }

  int R_size = R_.size();

  cid_stack_.emplace_back(c_id);
  for (auto &cp : C_) {
    if (cp.first.second == cur_l_size) {
      R_size++;
      if (cp.first.first < C_[c_id].first.first && cp.second >= level) {
        return false;
      } else if (cp.second >= level) {
        cp.second = level; // lock R set
      }
    }
  }

  // printf("L: ");
  // for (auto& lp : L_)
  //  if (lp.second == level) printf("%d ",lp.first);
  // printf("\nR: ");
  // for (int r : R_) printf("%d ", r);
  // for (auto& cp : C_)
  //  if (cp.first.second == cur_l_size) printf("%d ",cp.first.first);
  // printf("\n\n");

  max_level_ = std::max(max_level_, level);
  if (R_size >= min_r_size_)
    maximal_nodes_++;
  return true;
}

void LevelFinder::Pop() {
  double start = get_cur_time();
  int level = cid_stack_.size();
  int last_cid = cid_stack_.back();
  for (auto &cp : C_) {
    if (cp.second >= level)
      cp.second = level;
  }

  // recover L_
  for (auto &lp : L_) {
    if (lp.second == level - 1) {
      auto &l_neighbors = graph_->NeighborsL(lp.first);
      for (int i = 0, j = 0; i < C_.size() && j < l_neighbors.size();) {
        if (C_[i].first.first < l_neighbors[j])
          i++;
        else if (C_[i].first.first > l_neighbors[j])
          j++;
        else {
          if (C_[i].second >= level)
            C_[i].second = 0x7fffffff; // is active
          C_[i].first.second++;
          i++;
          j++;
        }
      }
    } else if (lp.second == level)
      lp.second = level - 1;
  }

  for (auto &cp : C_) {
    if (cp.second == level)
      cp.second = level - 1;
  }
  C_[last_cid].second = 0x7fffffff;
  recover_time_ += get_cur_time() - start;
}

LevelFinderCache::LevelFinderCache(BiGraph *graph, const char *name)
    : LevelFinder(graph, name) {}

void LevelFinderCache::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);
  cache_.resize(CACHE_SIZE / sizeof(int));
  graph_->Reorder(RInc);
  for (int v = 0; v < graph_->GetRSize(); v++) {
    if (graph_->NeighborsR(v).size() < min_l_size_)
      continue;
    if (v > 0 && graph_->NeighborsR(v - 1) == graph_->NeighborsR(v))
      continue;
    processing_nodes_++;
    std::map<int, int> c_map;

    L_.clear();
    R_.clear();
    C_.clear();
    for (int l : graph_->NeighborsR(v)) {
      L_.emplace_back(std::make_pair(l, 0));
      for (int r : graph_->NeighborsL(l)) {
        c_map[r]++;
      }
    }
    int cid_start = -1;

    for (auto c_node : c_map) {
      if (c_node.second == L_.size())
        R_.emplace_back(c_node.first);
      else if (c_node.second >= min_l_size) {
        if (c_node.first < v)
          cid_start = C_.size();
        C_.emplace_back(std::make_pair(c_node, 0x7fffffff));
      }
    }
    if (R_[0] < v)
      continue;
    cid_stack_.clear();
    maximal_nodes_++;
    if (CalUsingCache(-cid_start - 2))
      continue;

    cid_stack_.emplace_back(cid_start);
    while (!cid_stack_.empty()) {
      int last_cid = cid_stack_.back();
      cid_stack_.pop_back();
      int cur_cid = FindNext(last_cid);
      if (cur_cid != -1) {
        if (Push(cur_cid)) {
          if (CalUsingCache(cur_cid)) {
            Pop();
          } else {
            cid_stack_.emplace_back(cur_cid);
          }
        } else {
          Pop();
          cid_stack_.back() = cur_cid;
        }

      } else {
        if (!cid_stack_.empty()) {
          Pop();
        }
      }
    }
  }
  finish();
}

bool LevelFinderCache::CalUsingCache(int cur_cid) {
  int L_size = (cur_cid < 0) ? L_.size() : C_[cur_cid].first.second;
  int word_per_item = ((L_size - 1) >> 5) + 1;
  int item_size = 0;
  int start_item = -1;
  int level = cid_stack_.size();

  // if (word_per_item > 8)
  //   return false;

  auto CacheGet = [=](int item, int word) -> int {
    return cache_[item * word_per_item + word];
  };
  if (cur_cid < 0) {
    item_size = C_.size();
    start_item = -cur_cid - 1;

  } else {
    for (int i = 0; i < C_.size(); i++) {
      if (C_[i].second > level && C_[i].first.second >= min_l_size_) {
        if (start_item == -1 && i > cur_cid)
          start_item = item_size;
        item_size++;
      }
    }
  }

  if (item_size * word_per_item > cache_.size())
    return false;
  if (start_item != -1) { // there is some valid candidate vertex
    std::vector<int> cur_L;
    for (auto &lp : L_)
      if (lp.second == level)
        cur_L.emplace_back(lp.first);
    item_size = 0;
    // update cache
    for (auto &cp : C_) {
      if (cp.second > level && cp.first.second >= min_l_size_) {
        auto &r_neighbors = graph_->NeighborsR(cp.first.first);
        std::vector<int> val_bs(word_per_item, 0);
        for (int i = 0, j = 0; i < cur_L.size() && j < r_neighbors.size();) {
          if (cur_L[i] < r_neighbors[j])
            i++;
          else if (cur_L[i] > r_neighbors[j])
            j++;
          else {
            val_bs[i >> 5] |= (1 << (i & 0x1f));
            i++;
            j++;
          }
        }
        for (int k = 0; k < word_per_item; k++)
          cache_[item_size * word_per_item + k] = val_bs[k];
        item_size++;
      }
    }
    std::vector<int> items_stack;
    std::vector<int> id_stack;

    id_stack.emplace_back(start_item - 1);
    while (!id_stack.empty()) {
      int last_id = id_stack.back();
      id_stack.pop_back();
      std::vector<int> parent_item(word_per_item);
      std::vector<int> cur_item(word_per_item);
      for (int k = 0; k < word_per_item; k++) {
        int offset = items_stack.size() - word_per_item;
        parent_item[k] =
            id_stack.empty() ? 0xffffffff : items_stack[offset + k];
      }
      while (++last_id < item_size) {
        bool is_equal = true;
        bool is_zero = true;
        for (int k = 0; k < word_per_item; k++) {
          cur_item[k] = parent_item[k] & CacheGet(last_id, k);
          is_equal = is_equal && cur_item[k] == parent_item[k];
          is_zero = is_zero && cur_item[k] == 0;
        }
        if (!is_equal && !is_zero)
          break;
      }
      if (last_id < item_size) {
        bool is_maximal = true;
        for (int i = 0; i < last_id; i++) {
          bool parent_contained = true;
          bool cur_contained = true;
          for (int k = 0; k < word_per_item; k++) {
            int tmp_p = parent_item[k] & CacheGet(i, k);
            int tmp_cur = cur_item[k] & CacheGet(i, k);
            parent_contained = parent_contained && (parent_item[k] == tmp_p);
            cur_contained = cur_contained && (cur_item[k] == tmp_cur);
          }
          if (!parent_contained && cur_contained) {
            is_maximal = false;
            break;
          }
        }
        if (is_maximal) {
          maximal_nodes_++;
          for (int k = 0; k < word_per_item; k++)
            items_stack.emplace_back(cur_item[k]);
          id_stack.emplace_back(last_id);
        }
        id_stack.emplace_back(last_id);
        // else {
        //  id_stack.back() = last_id;
        //}
      } else {
        if (!id_stack.empty())
          for (int k = 0; k < word_per_item; k++)
            items_stack.pop_back();
        // id_stack.pop_back();
      }
    }
  }

  return true;
}

AMBEAFinderNaive::AMBEAFinderNaive(BiGraph *graph, const char *name)
    : BicliqueFinder(graph, name) {}

void AMBEAFinderNaive::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);
  graph_->Reorder(RInc);

  for (int v = 0; v < graph_->GetRSize(); v++){
    if (graph_->NeighborsR(v).size() < min_l_size_)
      continue;
    if (v > 0 && graph_->NeighborsR(v - 1) == graph_->NeighborsR(v))
      continue;
    processing_nodes_++;
    VertexSet &L = graph_->NeighborsR(v);
    VertexSet R, C;
    std::map<int, int> c_map;

    for(int l:L){
      for(int r:graph_->NeighborsL(l)){
        if(r > v)
          c_map[r]++;
      }
    }
    R.emplace_back(v);
    for(auto c_node:c_map){
      if (c_node.second == L.size()) 
        R.emplace_back(c_node.first);
      else if(c_node.second >= min_l_size_)
        C.emplace_back(c_node.first);
    }
    if (R.size() >= min_r_size_) {
      maximal_nodes_++;
      maximum_biclique_.CompareAndSet(L, R, graph_);
    }
    for(int vc : C){
      //if(!seq_except(graph_->NeighborsR(vc), L).empty()){
        biclique_find(L, R, C, vc);
      //}
    }
    //biclique_find(L, R, C, v);
  }
  finish();
}
void AMBEAFinderNaive::biclique_find(const VertexSet &L, const VertexSet &R,
                                     const VertexSet &C, int v){
  VertexSet L_prime, R_prime, C_prime, L_remove;
  processing_nodes_++;
  seq_intersect_diff(L, graph_->NeighborsR(v), L_prime, L_remove);
  R_prime = graph_->NeighborsL(L_prime);
  int idx = std::lower_bound(R.begin(), R.end(), C[0]) - R.begin();
  if (R_prime[idx] < C[0]) return;

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
  }
  for (int vc : C_prime) {
    if (seq_intersect_cnt(graph_->NeighborsR(vc), L_remove) > 0)
      biclique_find(L_prime, R_prime, C_prime, vc);
    else
      processing_nodes_++;
  }
}

AMBEAFinderInter::AMBEAFinderInter(BiGraph *graph, const char *name)
    : BicliqueFinder(graph, name) {}

void AMBEAFinderInter::Execute(int min_l_size, int min_r_size) {
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
    }
    for (auto c : C) {
      biclique_find(L, R, C, c.first);
    }
  }
  finish();
}

void AMBEAFinderInter::biclique_find(const VertexSet &L, const VertexSet &R,
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
  }
  for (int vc : C_cand) {
    if (L_prime.size() + L_remove.size() == L.size() ||
        seq_intersect_cnt(graph_->NeighborsR(vc), L_remove) > 0)
      biclique_find(L_prime, R_prime, C_prime, vc);
    else
      processing_nodes_++;
  }
}

AMBEAFinderIntra::AMBEAFinderIntra(BiGraph *graph, const char *name)
    : BicliqueFinder(graph, name) {}

void AMBEAFinderIntra::Execute(int min_l_size, int min_r_size) {
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
    }
    C = Partition(C, L);
    for (int vc : C) {
      biclique_find(L, R, C, vc);
    }
  }
  finish();
}

void AMBEAFinderIntra::biclique_find(const VertexSet &L, const VertexSet &R,
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
  C_prime = Partition(C_prime, L_prime);

  if (R_prime.size() >= min_r_size_) {
    maximal_nodes_++;
    maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
  }
  for (int vc : C_prime) {
    if (seq_intersect_cnt(graph_->NeighborsR(vc), L_remove) > 0)
      biclique_find(L_prime, R_prime, C_prime, vc);
    else
      processing_nodes_++;
  }
}

VertexSet AMBEAFinderIntra::Partition(const VertexSet &C,
                                      const VertexSet &L_prime, 
                                      int v) {
  GroupHelper g_helper;
  std::vector<int> group_array(C.size(), 0);
  VertexSet C_prime;

  for(int l:L_prime){
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

  for (int i = 0; i < C.size(); i++){
    int nc = g_helper.GetNc(group_array[i]);
    if (nc >= min_l_size_ && nc < L_prime.size() &&
        g_helper.IsFirst(group_array[i]) && C[i] > v){
      C_prime.emplace_back(C[i]);
    }
  }

  return C_prime;
}

IntraFinder::IntraFinder(BiGraph *graph, const char *name)
    : AMBEAFinderIntra(graph, name) {}

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
    }
    C = Partition(C, L);
    biclique_find(L, R, C);
  }
  finish();
}

void IntraFinder::biclique_find(const VertexSet &L, const VertexSet &R,
                                const VertexSet &C) {
  for (int i = 0; i < C.size(); i++){
    processing_nodes_++;
    VertexSet L_prime, R_prime, C_prime;
    L_prime = seq_intersect(L, graph_->NeighborsR(C[i]));
    R_prime = graph_->NeighborsL(L_prime);
    VertexSet R_add = seq_except(R_prime, R);
    if (R_add[0] == C[i]) {
      if (R_prime.size() >= min_r_size_) {
        maximal_nodes_++;
        maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
      }
      C_prime = Partition(C, L_prime, C[i]);

      // for (int j = i + 1; j < C.size(); j++) {
      //   int ln = seq_intersect_cnt(L_prime, graph_->NeighborsR(C[j]));
      //   if (ln >= min_l_size_ && ln < L_prime.size())
      //     C_prime.emplace_back(C[j]);
      // }
      // C_prime = Partition(C_prime, L_prime);
      biclique_find(L_prime, R_prime, C_prime);
    }
  }
}

void IntraFinder::biclique_find_passive(const VertexSet &L, const VertexSet &R,
                                        VertexSet &C) {
  for (int i = 0; i < C.size(); i++){
    if(C[i] < 0) continue;
    processing_nodes_++;
    VertexSet L_prime, R_prime, C_prime, L_remove;
    seq_intersect_diff(L, graph_->NeighborsR(C[i]), L_prime, L_remove);
    R_prime = graph_->NeighborsL(L_prime);
    VertexSet R_add = seq_except(R_prime, R);
    if(R_add[0] == C[i]){
      if (R_prime.size() >= min_r_size_) {
        maximal_nodes_++;
        maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
      }
      for (int j = i + 1; j < C.size(); j++) {
        if(C[j] < 0) continue;
        int ln = seq_intersect_cnt(L_prime, graph_->NeighborsR(C[j]));
        if (ln == L_prime.size() &&
            seq_intersect_cnt(L_remove, graph_->NeighborsR(C[j])) == 0)
          C[j] = -1;
        else if (ln >= min_l_size_ && ln < L_prime.size())
          C_prime.emplace_back(C[j]);
      }
      C_prime = Partition(C_prime, L_prime);
      biclique_find(L_prime, R_prime, C_prime);
    }
  }
}

AMBEAFinder::AMBEAFinder(BiGraph *graph, const char *name)
    : AMBEAFinderIntra(graph, name) {}

void AMBEAFinder::Execute(int min_l_size, int min_r_size) {
  setup(min_l_size, min_r_size);
  my_clock_ = 0;
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
    }
    C = Partition(C, L);
    for (int vc : C) {
      biclique_find(L, R, C, vc);
    }
  }
  finish();
  printf("%lfs\n", my_clock_);
}

void AMBEAFinder::biclique_find(const VertexSet &L, const VertexSet &R,
                                const VertexSet &C,  int v) {
  VertexSet L_prime, R_prime, C_prime, L_remove;
  processing_nodes_++;
  seq_intersect_diff(L, graph_->NeighborsR(v), L_prime, L_remove);
  R_prime = graph_->NeighborsL(L_prime);
  int idx = std::lower_bound(R.begin(), R.end(), C[0]) - R.begin();
  if (R_prime[idx] < C[0])
    return;
  
  int l_remove_size = 0;
  for (int i = 0; i < L_remove.size(); i++){
    auto l_neighbors = graph_->NeighborsL(L_remove[i]);
    int li = 0;
    for (int ri = 0; R_prime[ri] < v; ri++){
      while (li < l_neighbors.size() && l_neighbors[li] < R_prime[ri]) li++;
      if(li == l_neighbors.size() || l_neighbors[li] != R_prime[ri]){
        li = -1;
        break;
      }
    }
    if(li >= 0){
      if (l_remove_size != i) L_remove[l_remove_size] = L_remove[i];
      l_remove_size++;
    }
  }
  
  if (l_remove_size == 0) return;
  L_remove.resize(l_remove_size);
  if (R_prime.size() >= min_r_size_) {
    maximal_nodes_++;
    maximum_biclique_.CompareAndSet(L_prime, R_prime, graph_);
  }
  
  C_prime = Partition(C, L_prime, v);
  

  for (int vc : C_prime) {
    double ss = get_cur_time();
    int ln = seq_intersect_cnt(graph_->NeighborsR(vc), L_remove);
    my_clock_ += get_cur_time() - ss;
    if (ln > 0) 
      biclique_find(L_prime, R_prime, C_prime, vc);
  }
}