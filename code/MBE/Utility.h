#ifndef __UTILITY_H
#define __UTILITY_H
#include <algorithm>
#include <vector>
#include <stack>

/** Timer function **/

/* Returns current time in seconds as a double value
   and the precision can be as accurate as microseconds
   (10^-6 second)
 */
double get_cur_time();

/** Set operation **/

std::vector<int> seq_intersect(const std::vector<int>& v0,
  const std::vector<int>& v1);
void seq_intersect_local(std::vector<int> &v0, const std::vector<int> &v1);

int seq_intersect_cnt(const std::vector<int>& v0, const std::vector<int>& v1);
std::vector<int> seq_except(const std::vector<int>& v0,
  const std::vector<int>& v1);
int first_diff_element(const std::vector<int> &v0, const std::vector<int> &v1);

std::vector<int> seq_union(const std::vector<int>& v0,
  const std::vector<int>& v1);
std::vector<int> seq_intersect_upper(const std::vector<int>& v0,
                               const std::vector<int>& v1, int bound);
int seq_intersect_cnt_upper(const std::vector<int>& v0,
                            const std::vector<int>& v1, int bound);
int seq_intersect_cnt_lower(const std::vector<int> &v0,
                            const std::vector<int> &v1, int bound);

std::vector<int> seq_except_upper(const std::vector<int>& v0,
                                  const std::vector<int>& v1, int bound);
std::vector<int> seq_union_upper(const std::vector<int>& v0,
                                 const std::vector<int>& v1, int bound);

void seq_intersect_diff(const std::vector<int> &v0, const std::vector<int> &v1,
                        std::vector<int> &res_intersect,
                        std::vector<int> &res_diff);

/** Set operation helper **/
static long long g_set_op_counter;

void SetOpCounterReset();

void SetOpCounterAdd(int num);

long long GetSetOpCounter();

/** Memory util **/
unsigned int get_proc_vmhwm(unsigned int pid);

/** Bit operation helper **/
unsigned bit_count(unsigned i);

struct GroupInfoNode {
  int common_neighbor_size;
  int cur_neighbor; // check whether current v is
  int next_group_id;
  int size;
};

class GroupHelper {
public:
  GroupHelper();
  int GetCurGroupId(int last_group_id, int cur_neighbor);
  int GetNc(int group_id);
  bool IsFirst(int group_id);

private:
  std::vector<GroupInfoNode> groups_;
  std::stack<int> s_;
};
#define VISIT_FLAG 0x7ffffff

#endif /* __UTILITY_H */
