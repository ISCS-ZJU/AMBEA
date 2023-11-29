
#ifndef MBBP_DENSEMBB_H
#define MBBP_DENSEMBB_H
#include <vector>
#include "bgraph.h"
#include "LCUtility.h"
using namespace std;

struct sInf{
    int vistis = 0;
    int flag = -1;
};


class denseMBB{

private:
public:
    int *La;
    int *llookupA;

    int *Ra;
    int *rlookupA;
    int maxsize=0;
    int cou= 0;
    bgraph& bg;

    denseMBB(bgraph &g):bg(g){
        maxsize = 0;
    };

    //denseMBB(const bgraph &);

    void recEnuA(int *La, int *Ra, int *llookupA , int *rlookupA,int &a, int &ca, int &xl, int &b, int &cb, int& xr, vector<vector<int>> &adjL, vector<vector<int>> &adjR);
    void recEnu(unordered_set<int> &,unordered_set<int> &,unordered_set<int> &,unordered_set<int> &,int);

    void computeInG(unordered_map<int,unordered_set<int>> &, unordered_set<int>&, vector<unordered_set<int>> &,
                    unordered_map<int,unordered_set<int>> &, unordered_set<int>&, vector<unordered_set<int>> &);

    void computeRInG(unordered_map<int,unordered_set<int>> &, unordered_map<int,unordered_set<int>> &,unordered_set<int>&); // compute reverse induced subgraph

    void startEnu();
    void indstart();

    void printIG(unordered_map<int,unordered_set<int>> &,unordered_map<int,unordered_set<int>> &);

    int checkPCondition();
    void dynamicMBB(vector<pair<int,int>>&,unordered_set<int> &, unordered_map<int,unordered_set<int>> &,unordered_set<int> &, unordered_map<int,unordered_set<int>> &);
    void bfsLR(vector<pair<int,int>>&, unordered_map<int,unordered_set<int>> &, unordered_map<int,unordered_set<int>> &);

    //void generateSizes(pair<int,int> &, vector<int> &);

    /*
     * Have found a good way to compute induced subgraphs
     *
     *
     */
    void algDenseMBB();

    void calcuateSizes(pair<int,int> &, vector<pair<int,int>> & );

};



#endif //MBBP_DENSEMBB_H

