
#ifndef MBBP_LCUTILITY_H
#define MBBP_LCUTILITY_H

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <unordered_set>
#include "bgraph.h"
#include <stack>
#include <utility>
#include <stdlib.h>
class LCUtility {

};
void adjMtoAdj();
void mbsToMBB();
void coreFilter();
void eLtoAdjL();
void eLtoAdjM();
void cleanUp();
void cleanUpLJ();
void newEdgelist();
void newEdgeList(string, string, int);
void edgeFilter(string, string,int);
void computeIG(vector<int> &,bgraph &og, bgraph &ng);
void computeIG(vector<int> &,vector<int> &,bgraph &og, bgraph &ng);
//void readEdgeList(string &file);
void computeIG(vector<int> &,vector<int> &,bgraph &og, bgraph &ng,int);
void compute2IG(int v,bgraph &og, bgraph &ng,int);

struct vertexDeg{
    int v =0;
    int deg = 0;
};

struct vertexadj{
    int v =0;
    vector<int> list;
};

#endif //MBBP_LCUTILITY_H
