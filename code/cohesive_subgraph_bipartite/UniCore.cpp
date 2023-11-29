
#include "UniCore.h"

UniCore::UniCore(const vector<vector<int>> &adjL,const vector<vector<int>> &adjR){
    //assume vertex id (both sides) from 0 to ...
    maxCore = -1;
    offset = adjL.size(); //in uni-graph new id i' = i + offset
    noV = adjL.size()+adjR.size();
    deg = new int[noV]; // use to store the core number of each vertex

    uniCoreDecomposition(adjL,adjR);



}
/**
 *
 int maxC = 0;
    cout<<"****"<<endl;
    int count=0;
    int count2=0;
    for(int i=0;i<noV;i++){
        if(deg[i]>maxC){
            //if(deg[i]<200){
                maxC = deg[i];
           // }

        }
        if(deg[i]>=108){
            //cout<<i<<"has"<<deg[i]<<endl;
            count++;
        }
        if(deg[i]>=80){
            count2++;
        }

    }
    cout<<"vertices in 108:"<<count<<endl;
    cout<<"vertices in 80 core:"<<count2<<endl;
    cout<<"max core"<<maxC<<endl;
 */


int UniCore::getMaxCore() {
    if(maxCore<=-1){
        for(int i=0;i<noV;i++){
            if(deg[i]>maxCore){
                //if(deg[i]<200){
                maxCore = deg[i];
                // }

            }
        }

    }
    return maxCore;
}
void UniCore::getVertexSetCore(vector<int> &L, vector<int> &R, int th){
    for(int i=0;i<noV;i++){
        if(deg[i]>=th){
            if(i<offset){
                L.emplace_back(i);
            }else{
                R.emplace_back(i-offset);
            }
        }
    }
}

UniCore::~UniCore(){
    delete [] deg;
}

void UniCore::uniCoreDecomposition(const vector<vector<int>> &adjL,const vector<vector<int>> &adjR) {
    //step: get max degree
    int maxd = maxDeg(adjL,adjR);
    cout<<"max degree is "<< maxd<<endl;
    //cout<<"max degree"<<maxd<<endl;
    //step: for each degree, compute the number of vertices having the degree
    int *vsDeg = new int[maxd+1];
    //countVerSamDeg(vsDeg,adjL,adjR);
    for(int i=0;i<maxd+1;i++){
        vsDeg[i]=0;
    }
    for(auto &list: adjL){
        vsDeg[list.size()]=vsDeg[list.size()]+1;
    }

    for(auto &list:adjR){
        vsDeg[list.size()]=vsDeg[list.size()]+1;
    }

    //start the bucket sort
    //1: compute starting the position of vertices the same degree
    int *vlist = new int[noV]; // it is for storing the sorted vertex
    int pos =0;
    for(int i=0;i<maxd+1;i++){
        int tmepLoc = vsDeg[i];
        vsDeg[i] = pos;
        pos = tmepLoc+pos;
    }
    //cout<<"1"<<endl;
    //2:sort vertices and build an index for recording the position of each vertex in the sorted array
    int *vLoc = new int[noV];

    //L
    for(int i=0;i<adjL.size();i++){
        vLoc[i] = vsDeg[adjL[i].size()]; //adjL[i].size(): degree of i
        vlist[vLoc[i]] = i;
        vsDeg[adjL[i].size()]= vsDeg[adjL[i].size()]+1;
        deg[i] = adjL[i].size();
    }
    //cout<<"1.5"<<endl;
    //R
    for(int i=0;i<adjR.size();i++){
        vLoc[i+offset] = vsDeg[adjR[i].size()];
        vlist[vLoc[i+offset]] = i+offset;
        vsDeg[adjR[i].size()]= vsDeg[adjR[i].size()]+1;
        deg[i+offset] = adjR[i].size();
    }
    //cout<<"2"<<endl;
    // after the above loop, given a degree vsDeg[degree] is the starting position for vertices with degree+1
    for (int i=maxd;i>0;i--){
        vsDeg[i]=vsDeg[i-1];
    }
    vsDeg[0] = 0;
    //cout<<"3"<<endl;
    // start the core decomposition
    for (int i=0;i<noV;i++){

        int v = vlist[i];
        if(deg[v]>maxCore){
            maxCore = deg[v];
        }
        if(v<offset){//vertices in L
            vector<int> list = adjL[v];
            for(auto & u:list){ //note that u are in R
                if(deg[u+offset]>deg[v]){
                    int degu=deg[u+offset];
                    int positionu=vLoc[u+offset];
                    //get the vertex at the first position of degu+1
                    int positionw = vsDeg[degu];

                    int w = vlist[positionw];
                    //need to move if w!=u
                    if((u+offset)!=w){
                        vLoc[u+offset] = positionw;
                        vlist[positionu] =w;
                        vLoc[w] = positionu;
                        vlist[positionw] =u+offset;
                    }
                    vsDeg[degu]=vsDeg[degu]+1;
                    deg[u+offset] = deg[u+offset]-1;
                }
            }
        }else{
            //vertices in R
            vector<int> list =adjR[v-offset];
            for(auto & u:list){ //note that u are in L
                if(deg[u]>deg[v]){
                    int degu=deg[u];
                    int positionu=vLoc[u];
                    int positionw = vsDeg[degu];
                    int w = vlist[positionw];
                    if(u!=w){
                        vLoc[u] = positionw;
                        vlist[positionu] =w;
                        vLoc[w] = positionu;
                        vlist[positionw] =u;
                    }
                    vsDeg[degu]=vsDeg[degu]+1;
                    deg[u] = deg[u]-1;

                }

            }

        }

    }

    //clean garbage

    delete [] vLoc;
    delete [] vlist;
    delete[] vsDeg;

}

void UniCore::countVerSamDeg(int *vsDeg, const vector<vector<int>> &adjL, const vector<vector<int>>&adjR) {


}

int UniCore::maxDeg(const vector<vector<int>> &adjL,const vector<vector<int>> &adjR) {
    int maxL = maxDegOneSize(adjL);
    int maxR = maxDegOneSize(adjR);
    return maxL>maxR?maxL:maxR;
}

int UniCore::maxDegOneSize(const vector<vector<int>> &ajdList) {
    int max = 0;
    for(auto &list: ajdList){
        if(list.size()>max){
            max = list.size();
        }
    }
    return max;
}