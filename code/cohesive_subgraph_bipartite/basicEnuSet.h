
#ifndef MBBP_BASICENUSET_H
#define MBBP_BASICENUSET_H

#include <vector>
#include "bgraph.h"
#include <unordered_set>
#include <algorithm>
#include "LCUtility.h"
#include "UniCore.h"
//#include <boost/unordered_set.hpp>
class basicEnuSet {
public:
    bgraph bg;
    basicEnuSet(bgraph &g);
    int ub;
    int maxSize;
    void basicEnu(vector<int> &,vector<int> &,vector<int> &,vector<int> &, vector<vector<int>> &, vector<vector<int>> &,int&);

    bool isVinAdj(int,vector<int> &);
    virtual ~basicEnuSet();
    bool coreP;
    int min(int a, int b);
    int upper_bound(vector<int>&,vector<int>&,vector<int>&,vector<int>&,int&);
    int upper_bound(vector<int>&,vector<int>&,vector<int>&,vector<int>&,int,bgraph&);
    void betterEnu();
};


#endif //MBBP_BASICENUSET_H
