#include "Utility.h"
#include <string.h>
#include <unordered_set>
#ifdef _MSC_VER
#include <psapi.h>
#include <windows.h>
#define fopen64 fopen
double get_cur_time() {
  LARGE_INTEGER nFreq;
  LARGE_INTEGER nTime;
  QueryPerformanceFrequency(&nFreq);
  QueryPerformanceCounter(&nTime);
  double time = (double)nTime.QuadPart / (double)nFreq.QuadPart;
  return time;
}

#else
#include <stdio.h>
#include <sys/time.h> /* gettimeofday */
double get_cur_time() {
  struct timeval tv;
  struct timezone tz;
  double cur_time;
  gettimeofday(&tv, &tz);
  cur_time = tv.tv_sec + tv.tv_usec / 1000000.0;
  return cur_time;
}
#define VMHWM_LINE 21
unsigned int get_proc_vmhwm(unsigned int pid) {

  char file_name[64] = {0};
  FILE *fd;
  char line_buff[512] = {0};
  sprintf(file_name, "/proc/%d/status", pid);

  fd = fopen(file_name, "r");
  if (nullptr == fd) {
    return 0;
  }

  char name[64];
  int vmhwm;
  for (int i = 0; i < VMHWM_LINE - 1; i++) {
    fgets(line_buff, sizeof(line_buff), fd);
  }

  while (true) {
    fgets(line_buff, sizeof(line_buff), fd);
    sscanf(line_buff, "%s %d", name, &vmhwm);
    if (!strcmp(name, "VmHWM:")) {
      break;
    }
  }

  fclose(fd);

  return vmhwm;
}
#endif

std::vector<int> seq_intersect(const std::vector<int> &v0,
                               const std::vector<int> &v1) {
  g_set_op_counter++;
  std::vector<int> res;

  if (v0.size() > 8 * v1.size()) {
    auto iter_begin = v0.begin();
    auto iter_end = v0.end();

    for(int v: v1){
      auto iter_aim = std::lower_bound(iter_begin, iter_end, v);
      if (iter_aim == v0.end()) break;
      else if(*iter_aim == v)
        res.emplace_back(v);
      iter_begin = iter_aim;
    }
  } else if (v0.size() * 8 < v1.size()) {
    auto iter_begin = v1.begin();
    auto iter_end = v1.end();

    for (int v : v0) {
      auto iter_aim = std::lower_bound(iter_begin, iter_end, v);
      if (iter_aim == v1.end())
        break;
      else if(*iter_aim == v) res.emplace_back(v);
      iter_begin = iter_aim;
    }
  } 
  else {
    for (auto it0 = v0.begin(), it1 = v1.begin();
         it0 != v0.end() && it1 != v1.end();) {
      if (*it0 == *it1) {
        res.emplace_back(*it0);
        it0++;
        it1++;
      } else if (*it0 > *it1)
        it1++;
      else
        it0++;
    }
  }
  return res;
}

void seq_intersect_local(std::vector<int> &v0, const std::vector<int> &v1) {
  int v0_valid = 0;
  auto it = v1.begin();
  for (int i = 0; i < v0.size() && it != v1.end(); i++) {
    while (it != v1.end() && *it < v0[i]) it++;
    if(it != v1.end() && *it == v0[i]){
      if (i != v0_valid) 
        v0[v0_valid] = v0[i];
      v0_valid++;
    }
  }
  v0.resize(v0_valid);
}

int seq_intersect_cnt(const std::vector<int> &v0, const std::vector<int> &v1) {
  g_set_op_counter++;
  int cnt = 0;
  for (auto it0 = v0.begin(), it1 = v1.begin();
       it0 != v0.end() && it1 != v1.end();) {
    if (*it0 == *it1) {
      cnt++;
      it0++;
      it1++;
    } else if (*it0 > *it1)
      it1++;
    else
      it0++;
  }
  return cnt;
}

std::vector<int> seq_except(const std::vector<int> &v0,
                            const std::vector<int> &v1) {
  g_set_op_counter++;
  std::vector<int> res;
  for (auto it0 = v0.begin(), it1 = v1.begin(); it0 != v0.end();) {
    if (it1 == v1.end() || *it0 < *it1) {
      res.emplace_back(*it0);
      it0++;
    } else if (*it0 == *it1) {
      it0++;
      it1++;
    } else {
      it1++;
    }
  }
  return res;
}

int first_diff_element(const std::vector<int> &v0, const std::vector<int> &v1) {
  for (int i = 0; i < std::min(v0.size(), v1.size()); i++) {
    if (v0[i] != v1[i])
      return v0[i];
  }
  return (v0.size() > v1.size()) ? v0[v1.size()] : -1;
}

std::vector<int> seq_union(const std::vector<int> &v0,
                           const std::vector<int> &v1) {
  g_set_op_counter++;
  std::vector<int> res;
  for (auto it0 = v0.begin(), it1 = v1.begin();
       it0 != v0.end() || it1 != v1.end();) {
    if (it1 == v1.end())
      res.emplace_back(*it0++);
    else if (it0 == v0.end())
      res.emplace_back(*it1++);
    else if (*it0 == *it1) {
      res.emplace_back(*it0);
      it0++;
      it1++;
    } else if (*it0 > *it1) {
      res.emplace_back(*it1);
      it1++;
    } else {
      res.emplace_back(*it0);
      it0++;
    }
  }
  return res;
}

std::vector<int> seq_intersect_upper(const std::vector<int> &v0,
                                     const std::vector<int> &v1, int bound) {
  g_set_op_counter++;
  std::vector<int> res;
  for (auto it0 = std::upper_bound(v0.begin(), v0.end(), bound),
            it1 = std::upper_bound(v1.begin(), v1.end(), bound);
       it0 != v0.end() && it1 != v1.end();) {
    if (*it0 == *it1) {
      res.emplace_back(*it0);
      it0++;
      it1++;
    } else if (*it0 > *it1)
      it1++;
    else
      it0++;
  }
  return res;
}

int seq_intersect_cnt_upper(const std::vector<int> &v0,
                            const std::vector<int> &v1, int bound) {
  g_set_op_counter++;
  int cnt = 0;
  for (auto it0 = std::upper_bound(v0.begin(), v0.end(), bound),
            it1 = std::upper_bound(v1.begin(), v1.end(), bound);
       it0 != v0.end() && it1 != v1.end();) {
    if (*it0 == *it1) {
      cnt++;
      it0++;
      it1++;
    } else if (*it0 > *it1)
      it1++;
    else
      it0++;
  }
  return cnt;
}

int seq_intersect_cnt_lower(const std::vector<int> &v0,
                            const std::vector<int> &v1, int bound) {
  g_set_op_counter++;
  int cnt = 0;
  for (auto it0 = v0.begin(), it1 = v1.begin();
       it0 != v0.end() && it1 != v1.end() && *it0 < bound && *it1 < bound;) {
    if (*it0 == *it1) {
      cnt++;
      it0++;
      it1++;
    } else if (*it0 > *it1)
      it1++;
    else
      it0++;
  }
  return cnt;
}

std::vector<int> seq_except_upper(const std::vector<int> &v0,
                                  const std::vector<int> &v1, int bound) {
  g_set_op_counter++;
  std::vector<int> res;
  for (auto it0 = std::upper_bound(v0.begin(), v0.end(), bound),
            it1 = std::upper_bound(v1.begin(), v1.end(), bound);
       it0 != v0.end();) {
    if (it1 == v1.end() || *it0 < *it1) {
      res.emplace_back(*it0);
      it0++;
    } else if (*it0 == *it1) {
      it0++;
      it1++;
    } else {
      it1++;
    }
  }
  return res;
}

std::vector<int> seq_union_upper(const std::vector<int> &v0,
                                 const std::vector<int> &v1, int bound) {
  g_set_op_counter++;
  std::vector<int> res;
  for (auto it0 = std::upper_bound(v0.begin(), v0.end(), bound),
            it1 = std::upper_bound(v1.begin(), v1.end(), bound);
       it0 != v0.end() || it1 != v1.end();) {
    if (it1 == v1.end())
      res.emplace_back(*it0++);
    else if (it0 == v0.end())
      res.emplace_back(*it1++);
    else if (*it0 == *it1) {
      res.emplace_back(*it0);
      it0++;
      it1++;
    } else if (*it0 > *it1) {
      res.emplace_back(*it1);
      it1++;
    } else {
      res.emplace_back(*it0);
      it0++;
    }
  }
  return res;
}

void seq_intersect_diff(const std::vector<int> &v0, const std::vector<int> &v1,
                        std::vector<int> &res_intersect,
                        std::vector<int> &res_diff) {
  for (auto it0 = v0.begin(), it1 = v1.begin(); it0 != v0.end();){
    if(it1 == v1.end() || *it0 < *it1){
      res_diff.emplace_back(*it0);
      it0++;
    } else if(*it0 == *it1){
      res_intersect.emplace_back(*it0);
      it0++;
      it1++;
    } else{
      it1++;
    }
  } 
}

void SetOpCounterReset() { g_set_op_counter = 0; }

void SetOpCounterAdd(int num) { g_set_op_counter += num; }

long long GetSetOpCounter() { return g_set_op_counter; }

unsigned bit_count(unsigned i) {
  i = (i & 0x55555555) + ((i >> 1) & 0x55555555);
  i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
  i = (i & 0x0F0F0F0F) + ((i >> 4) & 0x0F0F0F0F);
  i = (i * (0x01010101) >> 24);
  return i;
}

GroupHelper::GroupHelper() {
  groups_.emplace_back(GroupInfoNode());
  groups_.back().common_neighbor_size = 0;
  groups_.back().cur_neighbor = -1;
  groups_.back().next_group_id = -1;
  groups_.back().size = 0x7fffffff;
}

int GroupHelper::GetCurGroupId(int last_group_id, int cur_neighbor) {
  int ret_group_id;
  if (groups_[last_group_id].cur_neighbor ==
      cur_neighbor) { // if next group exist
    ret_group_id = groups_[last_group_id].next_group_id;
    groups_[ret_group_id].size++;
  } else {            // if not exist
    if (s_.empty()) { // get from new element
      ret_group_id = groups_.size();
      groups_.emplace_back(GroupInfoNode());
    } else { // get from stack
      ret_group_id = s_.top();
      s_.pop();
    }
    // update group info
    groups_[last_group_id].cur_neighbor = cur_neighbor;
    groups_[last_group_id].next_group_id = ret_group_id;

    // initialize next group info
    groups_[ret_group_id].common_neighbor_size =
        groups_[last_group_id].common_neighbor_size + 1;
    groups_[ret_group_id].cur_neighbor = -1;
    groups_[ret_group_id].next_group_id = -1;
    groups_[ret_group_id].size = 1;
  }
  if (--groups_[last_group_id].size == 0) {
    s_.push(last_group_id);
  }
  return ret_group_id;
}

int GroupHelper::GetNc(int group_id) {
  return groups_[group_id].common_neighbor_size;
}

bool GroupHelper::IsFirst(int group_id) {
  bool ret = groups_[group_id].cur_neighbor != VISIT_FLAG;
  groups_[group_id].cur_neighbor = VISIT_FLAG;
  return ret;
}