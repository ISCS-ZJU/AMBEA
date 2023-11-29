

#include "basicEnuSet.h"
basicEnuSet::basicEnuSet(bgraph & b){

    bg = b;  //infact copy

    vector<int> A;
    vector<int> Ca;

    vector<int> B;
    vector<int> Cb;

    //coreP = true;
    //betterEnu();


      cout<<"start the algorithm"<<endl;
    //for(int i=0;i<b.adjL.size();i++){

    //    Ca.emplace_back(i);
   // }

    for(int i=b.adjL.size()-1;i>=0;i--){

        Ca.emplace_back(i);
    }

    //for(int i=0;i<b.adjR.size();i++){
//
 //       Cb.emplace_back(i);
   // }

    for(int i=b.adjR.size()-1;i>=0;i--){

        Cb.emplace_back(i);
    }
    maxSize = 0;
    int flag = 0;

    cout<<"start the algorithm 2"<<endl;


    basicEnu(A,Ca,B,Cb,b.adjL,b.adjR,flag);

}


void basicEnuSet::betterEnu(){
    vector<int> A;
    vector<int> Ca;

    vector<int> B;
    vector<int> Cb;

    for(int i=0;i<bg.adjR.size();i++){
        int v = i;
        bgraph nb;
        compute2IG(v,bg,nb,1);
        nb.sortAdj();

        //cout<<"new L:"<<nb.adjR.size()<<endl;
        //cout<<"new R:"<<nb.adjL.size()<<endl;
        //cout<<"edges:"<<nb.edgeNumber()<<endl;
        //cout<<"******"<<endl;
        for(int j=nb.adjR.size()-1;j>=0;j--){
            if(j!=i){
                Cb.emplace_back(j);
            }

        }
        //cout<<"3.5"<<endl;
        for(int j=0;j<nb.adjL.size();j++){

            Ca.emplace_back(j);
        }

        maxSize = 0;
        int flag = 1;
        B.emplace_back(i);
        //cout<<"n" <<endl;
        basicEnu(B,Cb,A,Ca,nb.adjR,nb.adjL,flag);
        //basicEnu(A,Ca,B,Cb,nb.adjL,nb.adjR,0);
    }

    /**
     *  for(int i=0;i<bg.adjL.size();i++){
        int v = i;
        cout<<"1" <<endl;
        //compute v 2-hop induced subgraph
        bgraph nb;

        compute2IG(v,bg,nb,0);
        cout<<"2" <<endl;
        nb.sortAdj();
        cout<<"3" <<endl;

        cout<<"new L:"<<nb.adjL.size()<<endl;
        cout<<"new R:"<<nb.adjR.size()<<endl;
        for(int j=0;j<nb.adjL.size();j++){
            if(j!=i){
                Ca.emplace_back(j);
            }

        }
        cout<<"3.5"<<endl;
        for(int j=0;j<nb.adjR.size();j++){

            Cb.emplace_back(j);
        }

        maxSize = 0;
        int flag = 0;
        A.emplace_back(i);
        cout<<"n" <<endl;
        basicEnu(A,Ca,B,Cb,nb.adjL,nb.adjR,0);

    }
     */


}

//adjust adjacency list to determine for induced sbugraph computation


void basicEnuSet::basicEnu(vector<int> &A, vector<int> &Ca,
                           vector<int> &B, vector<int> &Cb,
                           vector<vector<int>> &adjL, vector<vector<int>> &adjR, int & flag) {
    if(Ca.size() == 0 || Cb.size() == 0){
        int currentsize = A.size()+Ca.size();
        if(currentsize>B.size()+Cb.size()) {
            currentsize = B.size() + Cb.size();
        }

        //cout<<end<<endl;
        if(currentsize>maxSize){
            maxSize = currentsize;
            cout<<"*****"<<endl;
            cout<<maxSize<<endl;
            cout<<"*****"<<endl;
        }
        return;
        if(false){//debug purpose

            int offset1=0;
            int offset2=0;

            if(flag==0){
                offset1 =0;
                offset2 = adjL.size();
            }else{

                offset1 = adjL.size();
                offset2 =0;

            }


            vector<int> RA;
            RA.reserve(A.size()+Ca.size());
            RA.insert(RA.end(),A.begin(),A.end());
            RA.insert(RA.end(),Ca.begin(),Ca.end());


            vector<int> RB;
            RB.reserve(B.size()+Cb.size());
            RB.insert(RB.end(),B.begin(),B.end());
            RB.insert(RB.end(),Cb.begin(),Cb.end());
            cout<<"****"<<endl;
            for(auto &u:RA){
                cout<<u+offset1<<",";

            }
            cout<<endl;
            for(auto &u:RB){
                cout<<u+offset2<<",";

            }
            cout<<endl;
            cout<<"****"<<endl;
            return;
        }
    }
    if(min((A.size()+Ca.size()),(B.size() + Cb.size()))<=maxSize){
        return;
    }



    if(coreP){
        int ub = upper_bound(A,Ca,B,Cb, flag);

        //cout<<"losse:"<<min((A.size()+Ca.size()),(B.size() + Cb.size()))<<";tight"<<ub<<endl;
        if(ub<=maxSize){
            //cout<<"pruned"<<endl;
            return;
        }
    }



    //select a vertex from Ca
    int v = Ca.back();
    Ca.pop_back();
    A.emplace_back(v);


    //update Cb
    vector<int> newCb; //it should be adj(v) \cap Cb
    unordered_set<int> *adjofv;

    if(flag==0){
        //processing left
       // cout<<"start the algorithm 2.1:"<<v<<endl;
        //cout<<"start the algorithm 2.1:"<<bg.adjHL.size()<<endl;
        adjofv = & (bg.adjHL[v]);  //copy  refined !!
    }else{
        //processing r
        adjofv = & (bg.adjHR[v]);  //copy  refined !!
    }

    //cout<<"start the algorithm 3"<<endl;
    if(A.size()==1){
        //vector<int> adj = adjL[v];  //copy
        for(auto &u:adjL[v]){
            newCb.emplace_back(u);
        }
    }else{
        for(auto &u:Cb){
            if(adjofv->find(u)!=adjofv->end()){
                newCb.emplace_back(u);
            }
            /*
            bool in = isVinAdj(u,adjL[v]);
            if(in == true){
                newCb.emplace_back(u);
            }
             */
        }

    }

    //delete adjofv;

    //cout<<"start the algorithm 4"<<endl;
    //cout<<"start the algorithm 3"<<endl;

    int newflag = 0;
    if(flag == 0){
        newflag =1;
    }else{
        newflag=0;
    }


    basicEnu(B,newCb,A,Ca,adjR,adjL,newflag);


    //do not consider ad v in A
    int remoa =A.back();
    A.pop_back();
    if(remoa!=v){
        cout<<"><"<<endl;
    }
    basicEnu(A,Ca,B,Cb,adjL,adjR,flag);
}


int basicEnuSet::upper_bound(vector<int> &A,vector<int> &Ca,vector<int> &B,vector<int> &Cb, int & flag){
    //flag ==0 A->AdjL
    //flag ==1 A->AdjR

    vector<int> gA;
    for(auto&v:A){
        gA.emplace_back(v);
    }
    for(auto&v:Ca){
        gA.emplace_back(v);
    }

    vector<int> gB;
    for(auto&v:B){
        gB.emplace_back(v);
    }
    for(auto&v:Cb){
        gB.emplace_back(v);
    }
    bgraph nG;
    computeIG(gA,gB,bg,nG,flag);

    UniCore coreness(nG.adjL,nG.adjR);
    return coreness.getMaxCore();
}

bool basicEnuSet::isVinAdj(int u, vector<int> &list) {
    return binary_search(list.begin(), list.end(),u);
}



int basicEnuSet::min(int a, int b) {
    if (a <= b) {
        return a;
    } else {
        return b;
    }
}

basicEnuSet::~basicEnuSet(){


}
