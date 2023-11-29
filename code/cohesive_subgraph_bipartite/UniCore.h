

#ifndef MBBP_UNICORE_H
#define MBBP_UNICORE_H
#include "bgraph.h"

class UniCore {
public:
    //const bigraph
    int *deg= nullptr;
    int noV;
    int offset; //to make the bgraph unigraph
    int maxCore;
    UniCore();
    UniCore(const vector<vector<int>> &adjL,const vector<vector<int>> &adjR);
    void uniCoreDecomposition(const vector<vector<int>> &,const vector<vector<int>>&);

    int maxDeg(const vector<vector<int>> &,const vector<vector<int>>&);
    int maxDegOneSize(const vector<vector<int>>&);

    void countVerSamDeg(int *, const vector<vector<int>>&,const vector<vector<int>>&);

    void getVertexSetCore(vector<int>&,vector<int>&,int);
    int getMaxCore();
    virtual ~UniCore();
};


#endif //MBBP_UNICORE_H
