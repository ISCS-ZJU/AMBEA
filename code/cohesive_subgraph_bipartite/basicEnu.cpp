
#include "basicEnu.h"

basicEnu::basicEnu(const vector<vector<int>> &adjL, const vector<vector<int>> &adjR, int c[]) {
    maxsize = 0;
    //ini L
    core=c;
    La = new int[adjL.size()];
    lookupLa = new int[adjL.size()];
    for (int i = 0; i < adjL.size(); i++) {
        La[i] = i;
        lookupLa[i] = i;
    }

    //ini R
    Ra = new int[adjR.size()];
    lookupRa = new int[adjR.size()];
    for (int i = 0; i < adjR.size(); i++) {

        Ra[i] = i;
        lookupRa[i] = i;
    }


    cout << "finish ini!" << endl;

    //ini R P X position
    //left side
    //a=-1,ca=0,xl=L.size()

    //0: starting point of enumerated vertex
    //a-1: end of the enumerated vertex
    //a:starting point of the partial result
    //ca-1:end of the partial result
    //ca:starting point of the candidate
    //xl-1: end of the candidate
    //xl: starting point of the other vertices
    //adjL.size()-1 end point of the other vertices

    /*
     *
     */
    //int xa = 0;
    int ca = 0;
    int xl = adjL.size();


    //right side
    //b=-1,cb=0,xr=R.size()
    //int b = 0;
    int cb = 0;
    int xr = adjR.size();
    //ini the search let's start from add a v in L to ca
    int flag =0;// 0 processing La  1: processing R
    cout << "start!" <<endl;
    //recEnuA(La, Ra, lookupLa, lookupRa, a, ca, xl, b, cb, xr, adjL, adjR,flag);

    //recEnuI(La,lookupLa,ca,xl,Ra,lookupRa,cb,xr,adjL,adjR,0);
}





void basicEnu::recEnuI(vector<int> R, int *La, int *lookupLa, int ca,  int xl,
                       int *Ra, int *lookupRa , int cb, int xr,
                        vector<vector<int>> &adjL,  vector<vector<int>> &adjR,int flag){







}

void basicEnu::srecEnu(int *La, int *Ra, int *llookupA, int *rlookupA, int a, int ca, int b, int cb,
                       const vector<vector<int>> adjL, const vector<vector<int>> adjR) {

}




//0: starting point of enumerated vertex
//a-1: end of the enumerated vertex
//a:starting point of the partial result
//ca-1:end of the partial result
//ca:starting point of the candidate
//xl-1: end of the partial result
//xl: starting point of the other vertices
//adjL.size()-1 end point of the other vertices
void basicEnu::recEnuA(int *La, int *Ra, int *lookupLa, int *lookupRa, int a, int ca, int xl, int b, int cb, int xr,
                       const vector<vector<int>> adjL, const vector<vector<int>> adjR, int flag) {

    int tsize = min(xl-a,xr-b);
    if(tsize<maxsize){
        return;
    }

    if(ca==xl||cb==xr){
        int tsize = min(xl-a,xr-b);
        if(tsize>maxsize){
            for(int i=a;i<xl;i++){
                cout<<La[i]<<",";
            }
            cout<<endl;
            for(int i=b;i<xr;i++){
                cout<<Ra[i]<<",";
            }
            cout<<endl;
            maxsize = tsize;
            cout<<maxsize<<endl;
            cout<<flag<<endl;
            cout<<"*******"<<endl;
        }
        return;
    }

    int v = La[ca];
    ca=ca+1;
    vector<int> adjv = adjL[v];

    int nxr = 0;




        if (cb == b) {//this means  B = \emptyset move adjv from cb to adjv.size()
            nxr = cb;
            for (int i = 0; i < adjv.size(); i++) {

                //could be double checked
                int tempV = Ra[nxr];
                Ra[nxr] = adjv[i];
                Ra[lookupRa[adjv[i]]] = tempV;

                lookupRa[tempV] = lookupRa[adjv[i]];
                lookupRa[adjv[i]] = nxr;
                nxr++;
            }

        } else { //this means that B contains some partial results the element
            //the vertices between cb and xr-1 should and adjv.size()
            unordered_set<int> intersection;
            unordered_set<int> hash;
            //compute intersection :: can be improved by sort merge
            for (int i = 0; i < adjv.size(); i++) {
                hash.emplace(adjv[i]);
            }
            for (int i = cb; i <= xr - 1; i++) {
                if (hash.find(Ra[i]) != hash.end()) {
                    intersection.emplace(Ra[i]);
                }
            }

            nxr = cb; //organize the intersected vertices from cb to  intersection.size()

            for (auto &v:intersection) {
                int tempV = Ra[nxr];
                Ra[nxr] = v;
                Ra[lookupRa[v]] = tempV;

                lookupRa[tempV] = lookupRa[v];
                lookupRa[v] = nxr;
                nxr++;
            }


        }

    //TODO check if intersection contains the commoneighbours of [a] to [ca-1]


    //move the interset
    //cout<<"2"<<endl;
    //branch
    //recEnuA(La,Ra,llookupA,rlookupA,a,ca,xl,b,cb,xr,adjL,adjR);
    int newflag=0;
    if(flag==0){
        newflag=1;
    }
    if(flag==1){
        newflag=0;
    }

    recEnuA(Ra,La,lookupRa,lookupLa,b,cb,nxr,a,ca,xl,adjR,adjL,newflag);
    //not select v
    /*
     *
     * move v to the [a] if necessary
     *
     */
    int lofv = lookupLa[v];
    if(lofv!=a){
        int tv = La[a];
        La[a] = v;
        La[lofv]=tv;
        lookupLa[tv]=lofv;
        lookupLa[v]=a;
    }
    a=a+1;
    //cout<<"4"<<endl;

    //need to move back some vertices

    recEnuA(La,Ra,lookupLa,lookupRa,a,ca,xl,b,cb,xr,adjL,adjR,flag);

    //recEnuA(La,Ra,llookupA,rlookupA,a,ca,xl,b,cb,xr,adjL,adjR,0);
    //recEnuA(Ra,La,rlookupA,llookupA,b,cb,xr,a,ca,xl,adjR,adjL,1);
}






int basicEnu::min(int a, int b) {
    if (a <= b) {
        return a;
    } else {
        return b;
    }
}





//0: starting point of enumerated vertex
//a-1: end of the enumerated vertex
//a:starting point of the partial result
//ca-1:end of the partial result
//ca:starting point of the candidate
//xl-1: end of the partial result
//xl: starting point of the other vertices
//adjL.size()-1 end point of the other vertices
void basicEnu::recEnuAA(int *La, int *llookupA, int a, int ca, int xl,
                        int *Ra, int *rlookupA, int b, int cb, int xr,
                        const vector<vector<int>> adjL, const vector<vector<int>> adjR ) {
    //balanced clique generation
    // ca==xl || cb==xr
    //if(ca==adjL.size()||cb==adjR.size()){
    if(ca==xl||cb==xr){
        cout<<endl;
        return;
    }


    //select v from L to ca
    ca = ca+1;

    /*
     * operate R
     * need to operate on b, cb
     * move N(v) to the end of B, adjust rlookupA, and b
     */

    //get the vertex moved to A,
    int v = La[ca-1];

    //cout<<"1"<<endl;
    //get neighbours of v
    vector<int> adjv = adjL[v];
    //move adjL[v]\cap Ca     //R shall reduce to N(v)
    int nomv=0;
    for(int i=0;i<adjv.size();i++ ){
        //adjv[i] is the vertex of u \in N(v)  and N(v) are in Ra
        //rlookupA[adjv[i]] u's location in Ra
        if(rlookupA[adjv[i]]>=cb){
            /*
             * switch u and Ra[cb]
             */
            int tempV=Ra[cb];
            Ra[cb]=adjv[i];
            Ra[rlookupA[adjv[i]]]=tempV;
            /*
             * update rlookupA
             */
            rlookupA[tempV] = rlookupA[adjv[i]];
            rlookupA[adjv[i]] = cb;
            //cb++ move the starting location of Cb afterward by 1;
            cb=cb+1; //potential problem
            nomv++;
        }
    }
    //branch
    //recEnuA(La,Ra,llookupA,rlookupA,a,ca,xl,b,cb,xr,adjL,adjR);
    //recEnuAA(Ra,La,rlookupA,llookupA,b,cb,a,ca,adjR,adjL);

    //not select v
    /*
     * increase a by 1
     *
     */
     a=a+1;
     /*
      * Cb should
      */
     cb=cb-nomv;
    //recEnuAA(La,Ra,llookupA,rlookupA,a,ca,b,cb,adjL,adjR);
}




basicEnu::~basicEnu(){
    delete [] La;
    delete [] lookupLa;
    delete [] Ra;
    delete [] lookupRa;
}





