#include "utility.h"

/* Returns current time in seconds as a double value
   and the precision can be as accurate as microseconds
   (10^-6 second)
 */

#ifdef _MSC_VER
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
#include <sys/time.h> /* gettimeofday */
double get_cur_time() {
  struct timeval tv;
  struct timezone tz;
  double cur_time;
  gettimeofday(&tv, &tz);
  cur_time = tv.tv_sec + tv.tv_usec / 1000000.0;
  return cur_time;
}
#endif

std::vector<int> seq_intersect(const std::vector<int>& v0,
                               const std::vector<int>& v1) {
  std::vector<int> res;
  for (auto it0 = v0.begin(), it1 = v1.begin();
       it0 != v0.end() && it1 != v1.end();) {
    if (*it0 == *it1) {
      res.push_back(*it0);
      it0++;
      it1++;
    } else if (*it0 > *it1)
      it1++;
    else
      it0++;
  }
  return res;
}

int seq_contain(const std::vector<int>& v0, const std::vector<int>& v1) { // -1 v0 contains v1
  int mask = v0.size() > v1.size() ? -1 : (v0.size() == v1.size()) ? 0 : 1;
  auto it0 = v0.begin(), it1 = v1.begin();
  while (it0 != v0.end() && it1 != v1.end()) {
    int delta = *it0 - *it1;
    if (delta == 0) {
      it0++;
      it1++;
    } else if ((mask == 0) || ((mask ^ delta) < 0))
      return 0;
    else if (delta < 0)
      it0++;
    else
      it1++;
  }
  if ((it0 != v0.end() && mask > 0) || (it1 != v1.end() && mask < 0))
    return 0;
  if (mask == 0) mask = 1;
  return mask;
 }

bool MyBitset::test(int id) { 
  return values_[id]; 
}

void MyBitset::set(int id) { 
  values_[id] = true; 
}

void MyBitset::reset(int id) { 
  values_[id] = false; 
}

void MyBitset::resize(int size) { 
  values_.resize(size); 
}
