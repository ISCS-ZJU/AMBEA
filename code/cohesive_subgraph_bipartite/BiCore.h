
#ifndef MBBP_BICORE_H
#define MBBP_BICORE_H
#include "bgraph.h"
#include "LCUtility.h"
#include <algorithm>


using namespace std;
struct V2N{
    int v = 0;
    int n1 = 0;
    int n2 = 0;
    int n3 = 0;
    int heapid = -1;
};


struct bidegCMP{
    bool operator()(const V2N  *v1, const V2N *v2) const{

        return  min(v1->n1,v1->n2)>min(v2->n1,v2->n2);

    }
};

struct BFST{
    int v=-1;
    bool isV = false;
    int count=0;
    int loc_in_P=-1;
};
//ooMBE(bgraph &bg):g(bg){  nomb = 0;};
class BiCore {
public:
    int nov=0;
    int offset = 0;
    int *bcn= nullptr;
    int *uni_l_core = nullptr;
    int maxE2Core = 0;
    //bgraph &bg;
    bgraph &bg;
    BiCore(bgraph &g):bg(g){offset = bg.adjL.size();};
    //BiCore(bgraph g){bg=g; offset = bg.adjL.size();};
    //BiCore();
    vector<int> biGegOrder;
    vector<int> uniOrder;
    vector<int> index_for_uniOrder;
    //vector<int>* biGegOrder;
    void bicoreDec(); // 1-hop and 2-hop

    void IncresingDegreeOrder(vector<vector<int>> &);
    
    void NoOrder(vector<vector<int>> &);

    void compute2HopN(vector<vector<int>> &,vector<vector<int>> &, vector<V2N> &, int);

    void update2HopDegree(int, vector<vector<int>> &,vector<vector<int>> &, vector<vector<V2N*>> &, vector<V2N> &, int);

    void write2Hop();
    void showDegeneracyOrder();
    void e2hopCore(vector<vector<int>> &, vector<vector<int>>&); // 2-hop only

    void computeExa2HopN(vector<vector<int>> &,vector<vector<int>> &, vector<V2N> &);

    void computeOrderedSubgraph(vector<int>&,vector<vector<int>>&,vector<vector<int>>&,vector<V2N>&);

    virtual ~BiCore();

};




#endif //MBBP_BICORE_H
