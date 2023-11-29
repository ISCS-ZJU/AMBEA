
#ifndef MBBP_BASICENU_H
#define MBBP_BASICENU_H

#include <vector>
#include "bgraph.h"
#include <unordered_set>
#include "UniCore.h"
//#include <boost/unordered_set.hpp>
using namespace std;

class basicEnu{

private:
public:
    /*
 *online indices for mbb algorithm
 */
    //left vertex set
    //vector<int> L;    //it is partitioned into A, CA,X

    int *La;
    //right vertex set
    //vector<int> R;  //it is partitioned into B, CB,X
    int *Ra;
    int *core;
    //left lookup table
   //vector<int> llookup; //using vertex id to get its location in L
    int *lookupLa;
    //right lookup table
    //vector<int> rlookup;  //using vertex id to get its location in R
    int *lookupRa;
    int maxsize;

    //constructor:
    basicEnu(const vector<vector<int>> &adjL,const vector<vector<int>> &adjR, int []);


    virtual ~basicEnu();
    //void recEnu( vector<int> & L,vector<int> & R, vector<int> & llookup, vector<int> & rlookup,int a, int ca, int xl, int b, int cb, int xr,const bgraph bg);
    void recEnuA(int *La, int *Ra, int *llookupA , int *rlookupA,int a, int ca, int xl, int b, int cb, int xr,const vector<vector<int>> adjL, const vector<vector<int>> adjR,int flag);

    void recEnuI(vector<int> R, int *La, int *lookupLa, int ca,  int xl,
                 int *Ra, int *lookupRa, int cb, int xr,
                 vector<vector<int>> &adjL, vector<vector<int>> &adjR,int flag);


    /*
     * Let design a simplified version of recEnuA
     */
    void srecEnu(int *La, int *Ra, int *llookupA , int *rlookupA,int a, int ca, int b, int cb,const vector<vector<int>> adjL, const vector<vector<int>> adjR);

    void recEnuAA(int *La, int *llookupA, int a, int ca, int xl,
                  int *Ra, int *rlookupA, int b, int cb, int xr,
                  const vector<vector<int>> adjL, const vector<vector<int>> adjR);




    int min(int a, int b);


};
#endif //MBBP_BASICENU_H


