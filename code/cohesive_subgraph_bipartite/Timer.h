//
// Created by Luis on 1/8/2021.
//

#ifndef MBBP_TIMER_H
#define MBBP_TIMER_H
#include <cstdlib>
#include <sys/time.h>

class Timer {
public:
    Timer() { m_start = timestamp(); }
    void restart() { m_start = timestamp(); }
    long long elapsed() { return timestamp() - m_start; }

private:
    long long m_start;

    // Returns a timestamp ('now') in microseconds
    long long timestamp() {
        struct timeval tp;
        gettimeofday(&tp, nullptr);
        return ((long long)(tp.tv_sec))*1000000 + tp.tv_usec;
    }
};

#endif //MBBP_TIMER_H
