#ifndef __UTILITY_H
#define __UTILITY_H
#include <vector>
#include <algorithm>

/** Timer function **/

/* Returns current time in seconds as a double value
   and the precision can be as accurate as microseconds
   (10^-6 second)
 */
double get_cur_time();
std::vector<int> seq_intersect(const std::vector<int>& v0,
                               const std::vector<int>& v1);
int seq_contain(const std::vector<int>& v0, const std::vector<int>& v1);

class MyBitset {
public:
  bool test(int id);
  void set(int id);
  void reset(int id);
  void resize(int size); 
private:
  std::vector<bool> values_;
};

#endif /* __UTILITY_H */
