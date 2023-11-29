#include "BicliqueFinder.h"

namespace {
inline void RMerge(RType& r0, const RType& r1) {
#ifdef MMBEA_FULL_MODE
  r0.insert(r0.end(), r1.begin(), r1.end());
#else
  r0 += r1;
#endif  // MMBEA_FULL_MODE
}
}  // namespace

BicliqueFinder::BicliqueFinder(BiGraph* graph_in)
    : graph_(graph_in), processing_nodes_(0), maximal_nodes_(0) {
  l_bound_ = graph_in->GetLBound();
  r_bound_ = graph_in->GetRBound();
}

void BicliqueFinder::PrintResult(char* fn, double time_find) {
  FILE* fp = (fn == NULL) ? stdout : fopen(fn, "w");
  fprintf(fp, "[%s]\n", typeid(*this).name());
  fprintf(fp, "Total processing time: %lf seconds\n", time_find);
  fprintf(fp, "maximal nodes/processing nodes : %d/%d\n", maximal_nodes_,
          processing_nodes_);
  maximum_biclique_.Print(fp);
  if (fn != NULL) fclose(fp);
}

void BicliqueFinder::SetBound(int l_bound, int r_bound) {
  l_bound_ = l_bound;
  r_bound_ = r_bound;
}

const Biclique& BicliqueFinder::GetMaximumBiclique() {
  return maximum_biclique_;
}

BasicMbeaFinder::BasicMbeaFinder(BiGraph* graph_in)
    : BicliqueFinder(graph_in) {}

void BasicMbeaFinder::Execute(int ctrl) {
  maximum_biclique_.edges = 0;
  maximum_biclique_.graph = graph_;
  maximal_nodes_ = 0;
  processing_nodes_ = 0;
  std::vector<int> L, R, P, Q;
  for (int i = 0; i < graph_->GetLSize(); i++)
    if (graph_->GetNeighborsOfLVertex(i).size() > r_bound_) L.emplace_back(i);
  for (int i = 0; i < graph_->GetRSize(); i++)
    if (graph_->GetNeighborsOfRVertex(i).size() > l_bound_) P.emplace_back(i);
  biclique_find(L, R, P, Q);
}

void BasicMbeaFinder::biclique_find(std::vector<int> L, std::vector<int> R,
                                    std::vector<int> P, std::vector<int> Q) {
  while (!P.empty()) {
    std::vector<int> L_prime, R_prime, P_prime, Q_prime;
    bool is_maximal = true;
    int x = P.back();
    R_prime = R;
    R_prime.push_back(x);

    // for (int i = 0; i < L.size(); i++)
    //  if (graph_->EdgeExist(L[i], x)) L_prime.push_back(L[i]);

    L_prime = seq_intersect(L, graph_->GetNeighborsOfRVertex(x));
    if (L_prime.size() <= l_bound_) {
      P.pop_back();
      Q.push_back(x);
      continue;
    }
    processing_nodes_++;

    for (int i = 0; i < Q.size(); i++) {
      int neighbor_cnt = 0, v = Q[i];
      for (int node : L_prime) {
        if (graph_->EdgeExist(node, v)) neighbor_cnt++;
      }
      if (neighbor_cnt == L_prime.size()) {
        is_maximal = false;
        break;
      } else if (neighbor_cnt > 0) {
        Q_prime.push_back(v);
      }
    }

    if (is_maximal) {
      for (int i = 0; i < P.size(); i++) {
        int neighbor_cnt = 0, v = P[i];
        if (v == x) continue;
        for (int node : L_prime) {
          if (graph_->EdgeExist(node, v)) neighbor_cnt++;
        }
        if (neighbor_cnt == L_prime.size())
          R_prime.push_back(v);
        else if (neighbor_cnt > 0)
          P_prime.push_back(v);
      }

      if ((!P_prime.empty()) && (R_prime.size() + P_prime.size() > r_bound_))
        biclique_find(L_prime, R_prime, P_prime, Q_prime);

      if (R_prime.size() > r_bound_) {
        maximal_nodes_++;
        int biclique_edges = L_prime.size() * R_prime.size();
        if (biclique_edges > maximum_biclique_.edges) {
          maximum_biclique_.edges = biclique_edges;
          maximum_biclique_.left_nodes = std::move(L_prime);
          maximum_biclique_.right_nodes = std::move(R_prime);
        }
      }
    }
    P.pop_back();
    Q.push_back(x);
  }
}

ImproveMbeaFinder::ImproveMbeaFinder(BiGraph* graph_in)
    : BicliqueFinder(graph_in) {}

void ImproveMbeaFinder::Execute(int ctrl) {
  maximum_biclique_.edges = 0;
  maximum_biclique_.graph = graph_;
  maximal_nodes_ = 0;
  processing_nodes_ = 0;
  std::vector<int> L, R, P, Q;
  for (int i = 0; i < graph_->GetLSize(); i++)
    if (graph_->GetNeighborsOfLVertex(i).size() > r_bound_) L.emplace_back(i);
  for (int i = 0; i < graph_->GetRSize(); i++)
    if (graph_->GetNeighborsOfRVertex(i).size() > l_bound_) P.emplace_back(i);

  biclique_find(L, R, P, Q);
}

void ImproveMbeaFinder::biclique_find(std::vector<int> L, std::vector<int> R,
                                      std::vector<int> P, std::vector<int> Q) {
  while (!P.empty()) {
    std::vector<int> L_prime, R_prime, P_prime, Q_prime, L_comple, C;
    bool is_maximal = true;

    // Choose one vertex from candidate set
    int x = P.back();
    R_prime = R;
    R_prime.push_back(x);
    C.push_back(x);

    for (int i = 0; i < L.size(); i++)
      if (graph_->EdgeExist(L[i], x))
        L_prime.push_back(L[i]);
      else
        L_comple.push_back(L[i]);

    if (L_prime.size() <= l_bound_) {
      P.pop_back();
      Q.push_back(x);
      continue;
    }
    processing_nodes_++;

    for (int i = 0; i < Q.size(); i++) {
      int neighbor_cnt = 0, v = Q[i];
      for (int node : L_prime) {
        if (graph_->EdgeExist(node, v)) neighbor_cnt++;
      }
      if (neighbor_cnt == L_prime.size()) {
        is_maximal = false;
        break;
      } else if (neighbor_cnt > 0) {
        Q_prime.push_back(v);
      }
    }

    std::vector<std::pair<int, int>> vertex_neighbor_vec;

    if (is_maximal) {
      for (int i = 0; i < P.size(); i++) {
        int neighbor_cnt = 0, L_comple_cnt = 0, v = P[i];
        if (v == x) continue;
        for (int node : L_prime)
          if (graph_->EdgeExist(node, v)) neighbor_cnt++;

        if (neighbor_cnt == L_prime.size()) {
          R_prime.push_back(v);
          for (int node : L_comple)
            if (graph_->EdgeExist(node, v)) L_comple_cnt++;
          if (L_comple_cnt == 0) C.push_back(v);
        } else if (neighbor_cnt > 0) {
          vertex_neighbor_vec.push_back(std::pair<int, int>(v, neighbor_cnt));
        }
        // P_prime.push_back(v);
      }

      // Sorting candidates in P_prime by neighbor_cnt
      auto sort_by_degrees = [](const std::pair<int, int>& _left,
                                const std::pair<int, int>& _right) {
        return _left.second > _right.second;
      };
      std::sort(vertex_neighbor_vec.begin(), vertex_neighbor_vec.end(),
                sort_by_degrees);
      for (int i = 0; i < vertex_neighbor_vec.size(); i++)
        P_prime.push_back(vertex_neighbor_vec[i].first);

      if ((!P_prime.empty()) && (R_prime.size() + P_prime.size() > r_bound_))
        biclique_find(L_prime, R_prime, P_prime, Q_prime);
      if (R_prime.size() > r_bound_) {
        maximal_nodes_++;
        int biclique_edges = L_prime.size() * R_prime.size();
        if (biclique_edges > maximum_biclique_.edges) {
          maximum_biclique_.edges = biclique_edges;
          maximum_biclique_.left_nodes = L_prime;
          maximum_biclique_.right_nodes = R_prime;
        }
      }
    }

    for (auto node : C) P.erase(std::find(P.begin(), P.end(), node));
    Q.insert(Q.end(), C.begin(), C.end());
  }
}

MMbeaFinder::MMbeaFinder(BiGraph* graph_in) : BicliqueFinder(graph_in) {}

void MMbeaFinder::Execute(int ctrl) {
  processing_nodes_ = 0;
  maximal_nodes_ = 0;
  maximum_biclique_.edges = 0;
  maximum_biclique_.graph = graph_;
  finder_init(ctrl);
}

void MMbeaFinder::biclique_find(const std::vector<int>& L, const RType& R,
                                const std::vector<CExtNode>& C, int parent_r_id,
                                int r_id) {
  std::vector<int> L_prime =
      seq_intersect(L, graph_->GetNeighborsOfRVertex(r_id));
  RType R_prime = R;
  std::vector<CExtNode> C_prime;
  if (L_prime.size() <= l_bound_) return;
  std::vector<int> current_L = L;
  processing_nodes_++;

  std::vector<std::pair<int, const CExtNode*>> rid_c_list(C.size());
  // init rid_c_list
  for (int i = 0; i < C.size(); i++) {
    auto& c_node = C[i];
    rid_c_list[i] = std::make_pair(c_node.r_id, &C[i]);
  }

  struct CSetNode {
    int index;
    int cnt;
    int l_cnt;
    CSetNode(int index, int cnt, int l_cnt)
        : index(index), cnt(cnt), l_cnt(l_cnt) {}
  };
  std::vector<CSetNode> c_set_list;
  c_set_list.emplace_back(CSetNode(0, C.size(), 0));
  for (int l : L_prime) {
    int size = c_set_list.size();
    for (int i = 0; i < size; i++) {
      auto& c_set = c_set_list[i];
      if (c_set.cnt == 1) {
        if (graph_->EdgeExist(l, rid_c_list[c_set.index].first)) c_set.l_cnt++;
        continue;
      }
      int l_index = c_set.index, r_index = c_set.index + c_set.cnt - 1;
      while (l_index < r_index) {
        while (l_index <= r_index &&
               graph_->EdgeExist(l, rid_c_list[l_index].first))
          l_index++;
        while (l_index <= r_index &&
               !graph_->EdgeExist(l, rid_c_list[r_index].first))
          r_index--;
        if (l_index < r_index)
          std::swap(rid_c_list[l_index], rid_c_list[r_index]);
      }
      if (l_index == c_set.index)
        continue;  // all elements don't connect with l
      c_set.l_cnt++;
      if (l_index == c_set.index + c_set.cnt)
        continue;  // all elements connect with l
      int next_size = c_set.index + c_set.cnt - l_index;
      c_set.cnt = l_index - c_set.index;
      c_set_list.emplace_back(CSetNode(l_index, next_size, c_set.l_cnt - 1));
    }
  }
  typedef std::pair<int, const CExtNode*> RID_CLIST_T;
  std::sort(rid_c_list.begin(), rid_c_list.begin() + c_set_list[0].cnt,
            [&](RID_CLIST_T n0, RID_CLIST_T n1) -> bool {
              return n0.first < n1.first;
            });
  if (rid_c_list[0].first < parent_r_id) return;
  assert(c_set_list[0].l_cnt == L_prime.size());
  for (int i = 0; i < c_set_list[0].cnt; i++) {
    RMerge(R_prime, rid_c_list[i].second->R);
    if (rid_c_list[i].first < r_id) {
      current_L = seq_intersect(
          current_L, graph_->GetNeighborsOfRVertex(rid_c_list[i].first));
      if (L_prime.size() == current_L.size()) return;
    }
  }

  std::vector<int> cand_r_ids;
  // generate C_prime
  for (int i = 1; i < c_set_list.size(); i++) {
    if (c_set_list[i].l_cnt <= l_bound_) continue;
    int index = c_set_list[i].index;
    int size = c_set_list[i].cnt;

    C_prime.emplace_back(CExtNode());
    auto& c_node = C_prime.back();
    c_node.l_cnt = rid_c_list[index].second->l_cnt;
    c_node.R = rid_c_list[index].second->R;
    c_node.r_id = rid_c_list[index].second->r_id;

    // c_node.l_cnt = c_set_list[i].l_cnt;
    for (int j = 1 + index; j < size + index; j++) {
      auto* orig_c_node = rid_c_list[j].second;
      RMerge(c_node.R, orig_c_node->R);
      if (orig_c_node->r_id < c_node.r_id) {
        c_node.r_id = orig_c_node->r_id;
        c_node.l_cnt = orig_c_node->l_cnt;
      }
    }
    if (c_node.r_id > r_id && c_node.l_cnt > c_set_list[i].l_cnt) {
      if (current_L.size() == L.size() ||
          seq_intersect(current_L, graph_->GetNeighborsOfRVertex(c_node.r_id))
                  .size() > c_set_list[i].l_cnt)
        cand_r_ids.emplace_back(c_node.r_id);
    }
    c_node.l_cnt = c_set_list[i].l_cnt;
  }

  for (int this_r_id : cand_r_ids) {
    biclique_find(L_prime, R_prime, C_prime, r_id, this_r_id);
  }
  int l_size = L_prime.size(), r_size;
#ifdef MMBEA_FULL_MODE
  r_size = R_prime.size();
#else
  r_size = R_prime;
#endif  // MMBEA_FULL_MODE
  if (r_size > r_bound_) {
    maximal_nodes_++;
    int biclique_edges = l_size * r_size;
    if (biclique_edges > maximum_biclique_.edges) {
      maximum_biclique_.edges = biclique_edges;
      maximum_biclique_.left_nodes = std::move(L_prime);
#ifdef MMBEA_FULL_MODE
      maximum_biclique_.right_nodes = std::move(R_prime);
#endif
    }
  }
}

void MMbeaFinder::biclique_find_v2(const std::vector<int>& L, const RType& R,
                                   const std::vector<CExtNode>& C,
                                   int parent_r_id, int r_id) {
  std::vector<int> L_prime =
      seq_intersect(L, graph_->GetNeighborsOfRVertex(r_id));
  RType R_prime = R;
  std::vector<CExtNode> C_prime;
  std::vector<int> cand_r_ids;
  if (L_prime.size() <= l_bound_) return;
  std::vector<int> current_L = L;
  processing_nodes_++;

  for (int i = 0; i < C.size(); i++) {
    int c_id = C[i].r_id;
    int new_l_cnt =
        seq_intersect(L_prime, graph_->GetNeighborsOfRVertex(c_id)).size();
    if (new_l_cnt == L_prime.size()) {
      if (c_id <= parent_r_id)
        return;
      else {
        if (c_id < r_id) {
          current_L =
              seq_intersect(current_L, graph_->GetNeighborsOfRVertex(c_id));
          if (current_L.size() == L_prime.size()) return;
        }
        RMerge(R_prime, C[i].R);
      }
    } else if (new_l_cnt <= l_bound_)
      continue;
    else {
      C_prime.emplace_back(CExtNode(c_id, C[i].R, new_l_cnt));
      if (new_l_cnt != C[i].l_cnt && c_id > r_id) {
        if (current_L.size() == L.size() ||
            seq_intersect(current_L, graph_->GetNeighborsOfRVertex(c_id))
                    .size() != new_l_cnt)
          cand_r_ids.emplace_back(c_id);
      }
    }
  }

  for (int this_r_id : cand_r_ids) {
    biclique_find_v2(L_prime, R_prime, C_prime, r_id, this_r_id);
  }
  int l_size = L_prime.size(), r_size;
#ifdef MMBEA_FULL_MODE
  r_size = R_prime.size();
#else
  r_size = R_prime;
#endif  // MMBEA_FULL_MODE
  if (r_size > r_bound_) {
    maximal_nodes_++;
    int biclique_edges = l_size * r_size;
    if (biclique_edges > maximum_biclique_.edges) {
      maximum_biclique_.edges = biclique_edges;
      maximum_biclique_.left_nodes = std::move(L_prime);
#ifdef MMBEA_FULL_MODE
      maximum_biclique_.right_nodes = std::move(R_prime);
#endif
    }
  }
}

void MMbeaFinder::finder_init(int ctrl) {
  // double start, end;
  // start = get_cur_time();
  std::vector<int> L, R_vec;
  for (int i = 0; i < graph_->GetLSize(); i++) {
    if (graph_->GetNeighborsOfLVertex(i).size() > r_bound_) L.emplace_back(i);
  }
  for (int i = 0; i < graph_->GetRSize(); i++) {
    if (graph_->GetNeighborsOfRVertex(i).size() > l_bound_)
      R_vec.emplace_back(i);
  }
  // std::sort(R_vec.begin(), R_vec.end());
  RType R;
  std::vector<CExtNode> C;
#ifdef MMBEA_FULL_MODE
  R.clear();
#else
  R = 0;
#endif
  std::map<std::vector<int>, int> l_vector_map;

  for (int r_id : R_vec) {
    auto& l_vector = graph_->GetNeighborsOfRVertex(r_id);
    // if (l_vector.size() <= l_bound_) continue;
    RType R0;
#ifdef MMBEA_FULL_MODE
    R0.push_back(r_id);
#else
    R0 = 1;
#endif
    if (l_vector_map.find(l_vector) == l_vector_map.end()) {
      l_vector_map[l_vector] = C.size();
      C.emplace_back(CExtNode(r_id, R0, l_vector.size()));
    } else {
      int c_id = l_vector_map[l_vector];
      RMerge(C[c_id].R, R0);
    }
  }
  // end = get_cur_time();
  // printf("Init time : %lf\n", end - start);

  for (auto c_node : C) {
    if (ctrl == 0)
      biclique_find(L, R, C, -1, c_node.r_id);
    else
      biclique_find_v2(L, R, C, -1, c_node.r_id);
  }
}

PMbeaFinder::PMbeaFinder(PMbeaBiGraph* graph_in) : BicliqueFinder(graph_in) {}

void PMbeaFinder::Execute(int ctrl) {
  maximum_biclique_.edges = 0;
  maximum_biclique_.graph = graph_;
  maximal_nodes_ = 0;
  processing_nodes_ = 0;
  std::vector<int> L, R, P, Q;
  for (int i = 0; i < graph_->GetLSize(); i++)
    if (graph_->GetNeighborsOfLVertex(i).size() > r_bound_) L.emplace_back(i);
  for (int i = 0; i < graph_->GetRSize(); i++)
    if (graph_->GetNeighborsOfRVertex(i).size() > l_bound_) P.emplace_back(i);
  biclique_find(L, R, P, Q);
}

void PMbeaFinder::biclique_find(std::vector<int> L, std::vector<int> R,
                                std::vector<int> P, std::vector<int> Q) {
  std::queue<std::pair<int, int>> ex_q;

  while (!P.empty()) {
    std::vector<int> L_prime, R_prime, P_prime, Q_prime;
    bool is_maximal = true;
    int x = P.back();
    PMbeaBiGraph* p_graph = dynamic_cast<PMbeaBiGraph*>(graph_);
    auto& p = p_graph->range_index_[x];
    if (p.first <= p.second) ex_q.push(p);
    while (!ex_q.empty() && ex_q.front().first > x) ex_q.pop();
    if (!ex_q.empty() && ex_q.front().second >= x) {
      P.pop_back();
      Q.push_back(x);
      continue;
    }
    R_prime = R;
    R_prime.push_back(x);

    // for (int i = 0; i < L.size(); i++)
    //  if (graph_->EdgeExist(L[i], x)) L_prime.push_back(L[i]);
    L_prime = seq_intersect(L, graph_->GetNeighborsOfRVertex(x));

    if (L_prime.size() <= l_bound_) {
      P.pop_back();
      Q.push_back(x);
      continue;
    }
    processing_nodes_++;

    for (int i = 0; i < Q.size(); i++) {
      int neighbor_cnt = 0, v = Q[i];
      for (int node : L_prime) {
        if (graph_->EdgeExist(node, v)) neighbor_cnt++;
      }
      if (neighbor_cnt == L_prime.size()) {
        is_maximal = false;
        break;
      } else if (neighbor_cnt > 0) {
        Q_prime.push_back(v);
      }
    }

    if (is_maximal) {
      for (int i = 0; i < P.size(); i++) {
        int neighbor_cnt = 0, v = P[i];
        if (v == x) continue;
        for (int node : L_prime) {
          if (graph_->EdgeExist(node, v)) neighbor_cnt++;
        }
        if (neighbor_cnt == L_prime.size())
          R_prime.push_back(v);
        else if (neighbor_cnt > 0)
          P_prime.push_back(v);
      }

      if ((!P_prime.empty()) && (R_prime.size() + P_prime.size() > r_bound_))
        biclique_find(L_prime, R_prime, P_prime, Q_prime);

      if (R_prime.size() > r_bound_) {
        maximal_nodes_++;
        int biclique_edges = L_prime.size() * R_prime.size();
        if (biclique_edges > maximum_biclique_.edges) {
          maximum_biclique_.edges = biclique_edges;
          maximum_biclique_.left_nodes = std::move(L_prime);
          maximum_biclique_.right_nodes = std::move(R_prime);
        }
      }
    }
    P.pop_back();
    Q.push_back(x);
  }
}

MaximumEdgeBicliqueFinder::MaximumEdgeBicliqueFinder(BiGraph* graph)
    : orig_graph_(graph) {
  max_edges = 0;
  max_l.clear();
  max_r.clear();
}

void MaximumEdgeBicliqueFinder::Execute(int ctrl) {
  int lb = 1;
  int rb = orig_graph_->GetRBound() + 1;
  int last_lb, last_rb;
  if (orig_graph_->transpose_flag_) {
    for (int i = 0; i < orig_graph_->GetLSize(); i++)
      lb = std::max(lb, (int)(orig_graph_->GetNeighborsOfLVertex(i).size()));
  
  } else {
    for (int i = 0; i < orig_graph_->GetRSize(); i++)
      lb = std::max(lb, (int)(orig_graph_->GetNeighborsOfRVertex(i).size()));
  }
  max_edges = 0;
  last_rb = rb;
  last_lb = lb;
  lb = std::max(lb / 2, orig_graph_->GetLBound() + 1);
  
  BiGraph* graph;
  BicliqueFinder* finder;
  while (true) {
    if ((ctrl & 1) == 0)
      graph = new BitBiGraph();
    else
      graph = new BiGraph();
    graph->LoadGraphFromBiGraph(orig_graph_, lb - 1, rb - 1);
    if ((ctrl & 2) == 0)
      finder = new ImproveMbeaFinder(graph);  // new BasicMbeaFinder(graph);
    else {
      finder = new MMbeaFinder(graph);
    }
    finder->Execute();
    if (finder->maximum_biclique_.edges > max_edges) {
      max_edges = finder->maximum_biclique_.edges;
      printf("max_edges:%d\n", max_edges);
      max_l = finder->maximum_biclique_.left_nodes;
      max_r = finder->maximum_biclique_.right_nodes;
      graph->ConvertIdVector(max_l, max_r);
    }
    if (lb <= orig_graph_->GetLBound() + 1) break;
    last_lb = lb;
    last_rb = rb;
    rb = (max_edges - 1) / (lb - 1) + 1;
    lb = std::max(lb / 2, orig_graph_->GetLBound() + 1);
    //lb = (max_edges - 1) / (rb - 1) + 1;
    //rb = std::max(rb / 2, orig_graph_->GetRBound() + 1);
    delete graph;
    delete finder;
  }

}

void MaximumEdgeBicliqueFinder::Print() {
  printf("Maximum biclique edges: %d\n", max_edges);
  printf("L nodes[%d]:", (int)max_l.size());
  for (int l_id : max_l) printf("%d ", l_id);
  printf("\n");
  printf("R nodes[%d]:", (int)max_r.size());
  for (int r_id : max_r) printf("%d ", r_id);
  printf("\n\n");
}

MMbeaIntraFinder::MMbeaIntraFinder(BiGraph* graph_in)
    : BicliqueFinder(graph_in) {}

void MMbeaIntraFinder::Execute(int ctrl) {
  maximum_biclique_.edges = 0;
  maximum_biclique_.graph = graph_;
  maximal_nodes_ = 0;
  processing_nodes_ = 0;
  std::vector<int> L, R;
  std::vector<CExtNode> C;
  for (int i = 0; i < graph_->GetLSize(); i++)
    if (graph_->GetNeighborsOfLVertex(i).size() > r_bound_) L.emplace_back(i);
  for (int i = 0; i < graph_->GetRSize(); i++)
    if (graph_->GetNeighborsOfRVertex(i).size() > l_bound_) {
      C.emplace_back(CExtNode());
      C.back().r_cands.emplace_back(i);
    }
  biclique_find(L, R, C);

}

void MMbeaIntraFinder::biclique_find(const std::vector<int>& L,
                                     const std::vector<int>& R,
                                     std::vector<CExtNode>& C) {
  std::vector<CExtNode*> c_list(C.size());
  for (int i = 0; i < C.size(); i++)
    c_list[i] = &C[i];

  for (int i = 0; i < C.size(); i++) {
    if (C[i].is_visited) continue;
    std::vector<int> L_prime, R_prime;
    std::vector<CExtNode> C_prime;
    L_prime = seq_intersect(L, graph_->GetNeighborsOfRVertex(C[i].r_cands[0]));
    if (L_prime.size() > l_bound_) {
      processing_nodes_++;
      std::vector<std::pair<int, int>> c_set_list;  //(index, cnt)
      c_set_list.emplace_back(std::make_pair(0, C.size()));
      for (int l : L_prime) {
        int size = c_set_list.size();
        for (int j = 0; j < size; j++) {
          auto& c_set = c_set_list[j];
          int l_index = c_set.first, r_index = l_index + c_set.second - 1;
          while (l_index < r_index) {
            while (l_index <= r_index &&
                   graph_->EdgeExist(l, c_list[l_index]->r_cands[0]))
              l_index++;
            while (l_index <= r_index &&
                   !graph_->EdgeExist(l, c_list[r_index]->r_cands[0]))
              r_index--;
            if (l_index < r_index) 
              std::swap(c_list[l_index], c_list[r_index]);
          }
          if (l_index > c_set.first && l_index < c_set.first + c_set.second) {
            int next_size = c_set.first+c_set.second-l_index;
            c_set.second = l_index - c_set.first;
            c_set_list.emplace_back(std::make_pair(l_index, next_size));
          }
        }
      }
      // check_maximal
      bool is_maximal = true;
      for (int j = 0; j < c_set_list[0].second; j++)
        is_maximal = is_maximal && !c_list[j]->is_visited;
      if (is_maximal) {
        R_prime = R;
        for (int j = 0; j < c_set_list[0].second; j++) {
          auto& R_tmp = c_list[j]->r_cands;
          R_prime.insert(R_prime.end(), R_tmp.begin(), R_tmp.end());
        }
        for (int j = 1; j < c_set_list.size(); j++) {
          int index = c_set_list[j].first;
          int bound = index + c_set_list[j].second;
          bool is_visit = false;
          for (int k = index; k < bound; k++) 
            is_visit = is_visit || c_list[k]->is_visited;
          C_prime.emplace_back(CExtNode());
          auto& c_node = C_prime.back();
          if (!is_visit) {
            for (int k = index; k < bound; k++) {
              c_node.r_cands.insert(c_node.r_cands.end(),
                                    c_list[k]->r_cands.begin(),
                                    c_list[k]->r_cands.end());              
            }
          } else {
            c_node.is_visited = true;
            c_node.r_cands.emplace_back(c_list[index]->r_cands[0]);
          }
        }
        int cand_r_cnt = 0;
        for (auto& c_node:C_prime)
          if (!c_node.is_visited) cand_r_cnt += c_node.r_cands.size();
        if (cand_r_cnt > 0 && R_prime.size() + cand_r_cnt > r_bound_)
          biclique_find(L_prime, R_prime, C_prime);
        if (R_prime.size() > r_bound_) {
          maximal_nodes_++;
          int biclique_edges = L_prime.size() * R_prime.size();
          if (biclique_edges > maximum_biclique_.edges) {
            maximum_biclique_.edges = biclique_edges;
            maximum_biclique_.left_nodes = std::move(L_prime);
            maximum_biclique_.right_nodes = std::move(R_prime);
          }
        }
      }
    }
    C[i].is_visited = true;
  }
}
