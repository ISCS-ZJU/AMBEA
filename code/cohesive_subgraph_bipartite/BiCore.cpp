
#include <sys/time.h>
#include "BiCore.h"



void BiCore::bicoreDec() {
    nov = bg.adjL.size()+bg.adjR.size();
    bcn = new int [nov];
    //calculate 2-hop neighbours


    //int n2 = 0;


    //vector<V2N> v2degreesetL(bg.adjL.size());

    //vector<V2N> v2degreesetR(bg.adjL.size());
    vector<V2N> allV (nov); //index is the id of the vertex
    //cout<<"1"<<endl;
    compute2HopN(bg.adjL,bg.adjR,allV,0);

    //cout<<"2"<<endl;
    compute2HopN(bg.adjR,bg.adjL,allV,bg.adjL.size());


    int maxDeg2 = 0;
    for(auto &v: allV){
        if((v.n1+v.n2)>maxDeg2){
            maxDeg2 = v.n1+v.n2;
        }
    }

    int *buck = new int[maxDeg2+1]();
    for(auto &v: allV){
        ++buck[v.n1+v.n2];
    }


    int pos =0;
    for(int i=0;i<maxDeg2+1;i++){
        int tmepLoc = buck[i];
        buck[i] = pos;
        pos = tmepLoc+pos;
    }
    //buck store the starting position
    cout<<"max 2 degree is "<<maxDeg2<<endl;

    //first sort by 2-hop neighbours


    int *sortedV = new int [nov];
    int *vLoc = new int[nov]; // location of a vertex in a the sorted array.

    for(auto &v: allV){
        sortedV[buck[v.n1+v.n2]] = v.v;
        buck[v.n1+v.n2]++;
    }



    // after the above loop, given a degree vsDeg[degree] is the starting position for vertices with degree+1
    for (int i=maxDeg2;i>0;i--){
        buck[i]=buck[i-1];
    }
    buck[0] = 0;

    //collect vertices with the same 2-hop degree  allV
    //for each possible 2-hop degree, make a min-heap
    //if vertex with such 2-hop degree, make an empty-head
    vector<vector<V2N *>> heaps(maxDeg2+1);  //this ensure the V2N in every heap refers to the same V2N stored in allV
    for (int i=0;i<maxDeg2;i++){
        vector<V2N *> vl;
        //cout<<"bidegree is:"<<i<<endl;
        for(int j=buck[i];j<buck[i+1];j++){ // vertices with degree i
            allV[sortedV[j]].heapid = i;
            vl.emplace_back(&(allV[sortedV[j]]));
            //cout<<"vertex is "<< sortedV[j] <<" degrees are" <<allV[sortedV[j]].n1<<", "<<allV[sortedV[j]].n2<<", "<<allV[sortedV[j]].n1+allV[sortedV[j]].n2<<", and head id "<<allV[sortedV[j]].heapid<<endl;
        }

        make_heap(vl.begin(),vl.end(),bidegCMP());
        heaps[i]=vl;

    }

    //deal with the vertices with the maximum 2-hop degree
    vector<V2N*> vmaxdeg;
    //cout<<"bidegree is:"<<maxDeg2<<endl;
    for(int j = buck[maxDeg2];j<nov;j++){
        allV[sortedV[j]].heapid= maxDeg2;
        vmaxdeg.emplace_back(&(allV[sortedV[j]]));
        //cout<<"vertex is "<< sortedV[j] <<" degrees are" <<allV[sortedV[j]].n1<<", "<<allV[sortedV[j]].n2<<", "<<allV[sortedV[j]].n1+allV[sortedV[j]].n2<<", and head id "<<allV[sortedV[j]].heapid<<endl;
    }
    make_heap(vmaxdeg.begin(),vmaxdeg.end(),bidegCMP());
    heaps[maxDeg2]=vmaxdeg;

    for (int i=0;i<maxDeg2+1;i++){
        //get the heap
        vector<V2N*> *heap = &(heaps[i]);
        while(heap->size()){
            pop_heap(heap->begin(),heap->end(),bidegCMP());
            V2N *v =  heap->back();
            //loop on the two-hop neighbours
            //check the id of

            heap->pop_back();
            //
            if(v->heapid==i){
                //cout<<"vertex id is "<<v->v<<" and its bicoreness is "<<(v->n1+v->n2)<<" and its upper bound is "<<min(v->n1,v->n2)<<", 1-hop is "<<v->n1<<", and 2-hop is "<<v->n2<<endl;
                if((v->v)<(bg.adjL.size())){
                    //v is in L
                    update2HopDegree(v->v,bg.adjL,bg.adjR, heaps, allV, 0);
                }else{
                    // v is in R
                    update2HopDegree(v->v,bg.adjR,bg.adjL, heaps, allV, bg.adjL.size());
                }


                biGegOrder.emplace_back(v->v);
            }//else{

                //cout<<v->heapid<<" vs real id "<<i<<endl;
            //}



            //int heapSize = heap->size();



        }

    }

    delete [] vLoc;
    delete [] sortedV;
    delete [] buck;
    //show the sorted v debug purpose
    //for(int i=0;i<nov;i++){
    //    cout<<"vertex "<<sortedV[i]<<" has the 2-hop degree of "<<allV[sortedV[i]].n1+allV[sortedV[i]].n2<<endl;
    //}


    //L

    //for(auto & c:v2degreesetL){
    //    cout<<c.v<<","<<c.n1<<","<<c.n2<<endl;
    //}

//    cout<<"3"<<endl;
//    for(auto & c:v2degreeset){
//        delete c;
//    }
//    for(auto & c:v2degreeset){
//        cout<<c->v<<","<<c->n1<<c->n2<<endl;
//    }
}

void BiCore::update2HopDegree(int v, vector<vector<int>> &adjL, vector<vector<int>> &adjR, vector<vector<V2N*>> &heaps, vector<V2N> &allV, int offset) {
    //use relative relationship to simplify the code
    //pay extra attentions to the id changes!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //last step of the bicore decomposition
    vector<BFST> bfs (adjL.size());

    //v is vertex id passing in
    //we shall work on v-offset
    for(auto vertxL:adjL[v-offset]){
        //those are 1-hop neighbours of v
        int realvertxL;
        if(offset==0){
            //vertexL is actually in R, its id in allV is vertxL is vertxL+offset
            realvertxL = vertxL+adjL.size();
        }else{
            //vertexL is actually in L its id in allV is vertxL
            realvertxL = vertxL;
        }
        //now works on realvertxL
        //check the current total 2-hop degree of realvertxL
        if(( allV[realvertxL].n1+allV[realvertxL].n2)>(allV[v].n1+allV[v].n2) ){
            // the 1-hop neighbour for realvertxL shall be reduced by 1 (allV[realvertxL].n1+allV[realvertxL].n2)
            // put them into corresponding heap and update their heap id
            allV[realvertxL].n1 = allV[realvertxL].n1-1;
            int newHeapid = allV[realvertxL].heapid-1;
            //cout<<"new heap id "<<newHeapid<<endl;
            allV[realvertxL].heapid = newHeapid;
            heaps[newHeapid].push_back(&(allV[realvertxL]));
            push_heap(heaps[newHeapid].begin(),heaps[newHeapid].end(),bidegCMP());
        }

            //next, loop on the 2-hop neighbours of v
            //we work on the original id of vertxL
            //cout<<"test 3 " <<endl;
            for(auto vnn:adjR[vertxL]){
                //vnn is at the same side as v
                //Be careful !!!
                int newvnn;
                if(offset==0){
                    //vnn is in L
                    newvnn = vnn;
                }else{
                    //vnn is actually in R its id in allV
                    newvnn = vnn+offset;
                    //newvnn = vnn+adjL.size();
                }


                if(( allV[newvnn].n1+allV[newvnn].n2)>(allV[v].n1+allV[v].n2)){
                    if(vnn!=v-offset){
                        //cout<<"test 2 "<< vnn <<endl;

                        if(bfs[vnn].v!=v-offset){

                            //it is not visited yet
                            bfs[vnn].v = v-offset;

                            //cout<< "vnn is "<<vnn<<" and new vnn id is "<<newvnn<<endl;
                            //cout<<"size of allV is " <<sizeof(allV)<<endl;
                            //cout<< newvnn <<endl;
                            //cout<<allV[newvnn].v<<", "<<allV[newvnn].n1<<", "<<allV[newvnn].n2<<endl;
                            allV[newvnn].n2 = allV[newvnn].n2 - 1;
                            //cout<<"test"<<endl;
                            //cout<<"test!"<<endl;
                            int newHeapid = allV[newvnn].heapid;
                            newHeapid = newHeapid - 1;
                            //cout<<allV[newvnn].n2<<","<<allV[newvnn].n1<<endl;
                            allV[newvnn].heapid = newHeapid;
                            //cout<<newHeapid<<endl;
                            heaps[newHeapid].push_back(&(allV[newvnn]));
                            push_heap(heaps[newHeapid].begin(), heaps[newHeapid].end(),bidegCMP());


                        }
                    }
                }

                //if((allV[newvnn].n1+allV[newvnn].n2)
                //check the current degree of newvnn
                //cout<<"test 3 " <<endl;


            }





    }




}

void BiCore::compute2HopN(vector<vector<int>> &adjL,vector<vector<int>> &adjR, vector<V2N> & store, int offset) {
    vector<BFST> bfs (adjL.size());
    //cout<<"3"<<endl;
    for(int i=0;i<adjL.size();i++){
        int n1 = 0;
        int n2 = 0;
        n1 = adjL[i].size();
        vector<int> n2s;
        // index is the vertex

        for(auto &v:adjL[i]){ //v is in R
            for(auto &u: adjR[v]){
                if(u!=i){
                    if(bfs[u].v!=i){
                        bfs[u].v = i;
                        n2++;
                    }
                }
            }
           // vector<int> c(n2s.size()+adjR[v].size());
            //std::vector<int>::iterator it=set_union(n2s.begin(),n2s.end(),adjR[v].begin(),adjR[v].end(),c.begin());
            //c.resize(it-c.begin());
            //n2s = c;
        }
        //n2 = n2s.size()-1; //exclude i itself
        V2N u{i+offset,n1,n2, n1+n2};
        store[i+offset]=u;

        //cout<< "1-hop neighbours "<<n1<<endl;
        //cout<< "2-hop neighbours "<<n2<<endl;
        //cout<< "totoal 2-hop neighbours "<<n1+n2<<endl;

    }
}

void BiCore::write2Hop() {
    vector<V2N> v2degreesetL(bg.adjL.size());

    vector<V2N> v2degreesetR(bg.adjR.size());

    cout<<"1"<<endl;
    compute2HopN(bg.adjL,bg.adjR,v2degreesetL,0);

    cout<<"2"<<endl;
    compute2HopN(bg.adjR,bg.adjL,v2degreesetR,0);

    ofstream fout;
    fout.open("../2hopljL.txt");
    for(auto & c:v2degreesetL){
        fout<<c.v<<","<<c.n1<<","<<c.n2<<endl;
    }

    ofstream fout2;
    fout2.open("../2hopljR.txt");
    for(auto & c:v2degreesetR){
        fout2<<c.v<<","<<c.n1<<","<<c.n2<<endl;
    }


}

BiCore::~BiCore() {

   if(bcn!= nullptr) delete [] bcn;
    bcn=nullptr;
    if(uni_l_core!= nullptr)  delete [] uni_l_core;
    uni_l_core = nullptr;
   //cout<<"BiCore"<<endl;
}

void BiCore::showDegeneracyOrder() {
    if(biGegOrder.size()!=0){

        for(auto &v:biGegOrder){
            cout<<v<<",";
        }
        cout<<endl;

    }
}

void BiCore::NoOrder(vector<vector<int>> &adjL)
{
    index_for_uniOrder.resize(adjL.size());
    uniOrder.resize(adjL.size());
    for(int i = 0; i < uniOrder.size(); i++)
    {
        index_for_uniOrder[i] = i;
        uniOrder[i] = i;
    }
}

void BiCore::IncresingDegreeOrder(vector<vector<int>> &adjL)
{
    index_for_uniOrder.resize(adjL.size());
    uniOrder.resize(adjL.size());
    for(int i = 0; i < uniOrder.size(); i++)uniOrder[i] = i;
    sort(uniOrder.begin(), uniOrder.end(), [&](int a, int b){
         return (adjL[a].size() < adjL[b].size() || (adjL[a].size() == adjL[b].size() && adjL[a] < adjL[b]));
         });
    for(int i = 0; i < uniOrder.size(); i++)
    {
        index_for_uniOrder[uniOrder[i]] = i;
    }
}

//always compute e2hop core numbers for the vertices in adjL
void BiCore::e2hopCore(vector<vector<int>> &adjL,vector<vector<int>> &adjR) {
    auto get_cur_time = []()
    {
      struct timeval tv;
      struct timezone tz;
      double cur_time;
      gettimeofday(&tv, &tz);
      cur_time = tv.tv_sec + tv.tv_usec / 1000000.0;
      return cur_time;
    };
    auto start =get_cur_time();
        //compute two hop neighbours
    vector<V2N> allV (adjL.size()); //index is the id of the vertex
    uni_l_core = new int [adjL.size()];
    //cout<<"1"<<endl;
    //cout<<"start 2-hop neighbour"<<endl;
    computeExa2HopN(adjL,adjR,allV); //pass
    //cout<<adjR.size()<<endl;
    //cout<<adjL.size()<<endl;


    int maxDeg2 = 0;
    for(auto &v: allV){
        if(v.n2>maxDeg2){
            maxDeg2 = v.n2;
        }
    }

    int *buck = new int[maxDeg2+1]();
    for(auto &v: allV){
        buck[v.n2]++;
    }


    int pos =0;
    for(int i=0;i<maxDeg2+1;i++){
        int size = buck[i];
        buck[i] = pos;
        pos = size+pos;
    }
    //buck store the starting position
    //cout<<"max 2 degree is "<<maxDeg2<<endl;

    //first sort by 2-hop neighbours


    int *sortedV = new int [adjL.size()];
    int *vLoc = new int[adjL.size()]; // location of a vertex in a the sorted array.

    for(auto &v: allV){
        sortedV[buck[v.n2]] = v.v;
        vLoc[v.v] = buck[v.n2];
        buck[v.n2]++;
    }



    // after the above loop, given a degree buck[degree] is the starting position for vertices with degree+1
    for (int i=maxDeg2;i>0;i--){
        buck[i]=buck[i-1];
    }
    buck[0] = 0;

    int maxmax = 0;

    //debug 1
//    for(int i=0;i<maxDeg2+1;i++){
//        cout<<"2-hop degree is:"<<i<<", starting position is"<<buck[i]<<endl;
 //   }


    //for(int i=0;i<adjL.size();i++){
    //    int u = sortedV[i];
        //if(u==117) {
    //        cout << "vertex id is:" << u <<endl;// ", and 2 hop neighbour is:" << allV[u].n2 << ", vertex position is:" << i
                 //<< ", position in vlook is" << vLoc[u] << endl;
       // }
    //}




    //vector<int> order;
    index_for_uniOrder.resize(adjL.size());
    vector<BFST> bfs (adjL.size());
    for(int i=0;i<adjL.size();i++){
        int u = sortedV[i];
        uniOrder.emplace_back(u);
        index_for_uniOrder[u]=i;
        uni_l_core[u] = allV[u].n2;
        if(allV[u].n2>maxmax){
            maxmax = allV[u].n2;
        }
        //cout<<u<<endl;
        //if(i!=vLoc[u]){
        //    cout<<"why!!!"<<endl;
        //    continue;
        //}
        //loop on exact 2-hop neighbours of v
        //cout<<"vertex id is:"<<u<<", and 2 hop neighbour is:"<<allV[u].n2<<", vertex position is:"<<i<<", position in vlook is "<<vLoc[u]<<endl;

        for(auto &v:adjL[u]){ // v is in R
            //vector<BFST> bfs (adjL.size());
            for(auto &u_prime:adjR[v]){

                    if(bfs[u_prime].v!=u){
                        bfs[u_prime].v=u;
                        //do core job

                        if(allV[u_prime].n2>allV[u].n2){
                            //cout<<"***"<<"u_prime is"<<u_prime<<", 2-hop neighbour is"<<allV[u_prime].n2<<endl;
                            int deg_u_prime=allV[u_prime].n2;
                            int p_u_prime=vLoc[u_prime];
                            //get the vertex at the first position of degu+1
                            int positionw = buck[deg_u_prime];

                            int w = sortedV[positionw];
                            //cout<<"++++++++"<<"w is"<<w<<", w position is "<<positionw<<endl;
                            //need to move if w!=u
                            if(u_prime!=w){
                                vLoc[u_prime] = positionw;
                                sortedV[p_u_prime] = w;
                                vLoc[w] = p_u_prime;  //the previous u here almost kill me.
                                sortedV[positionw] = u_prime;
                            }
                            buck[deg_u_prime]++;
                            allV[u_prime].n2 = allV[u_prime].n2-1;
                        }
                    }

            }
        }
    }


    //report exact 2 hop core

    //for(int i=0;i<allV.size();i++){
    //    if(allV[i].n2==maxmax){
    //        cout<<"vertex:"<<i<<", exact 2-hop core number is:"<<allV[i].n2<<endl;
    //    }

    //
    // }
    //cout<<"max:"<<maxmax<<endl;
    //uniOrder = order;


    maxE2Core = maxmax;

    //computeOrderedSubgraph(biGegOrder,adjL,adjR,allV);

    // cout<<endl;
    // cout<<"ucorder time: "<<get_cur_time() - start<<endl;
    delete [] vLoc;
    delete [] sortedV;
    delete [] buck;

    //show the sorted v debug purpose
    //for(int i=0;i<nov;i++){
    //    cout<<"vertex "<<sortedV[i]<<" has the 2-hop degree of "<<allV[sortedV[i]].n1+allV[sortedV[i]].n2<<endl;
    //}


    //L

    //for(auto & c:v2degreesetL){
    //    cout<<c.v<<","<<c.n1<<","<<c.n2<<endl;
    //}

//    cout<<"3"<<endl;
//    for(auto & c:v2degreeset){
//        delete c;
//    }
//    for(auto & c:v2degreeset){
//        cout<<c->v<<","<<c->n1<<c->n2<<endl;
//    }

}

void BiCore::computeExa2HopN(vector<vector<int>> &adjL, vector<vector<int>> &adjR, vector<V2N> &store) {
    vector<BFST> bfs (adjL.size());
    //cout<<"3"<<endl;
    for(int i=0;i<adjL.size();i++){
        //int n1 = 0;
        int n2 = 0;
        //n1 = adjL[i].size();
        //vector<int> n2s;
        // index is the vertex

        for(auto &v:adjL[i]){ //v is in R
            for(auto &u: adjR[v]){
                if(u!=i){
                    if(bfs[u].v!=i){
                        bfs[u].v = i;
                        n2++;
                    }
                }
            }

        }
        V2N u{i,0,n2, n2};
        store[i]=u;

        //cout<< "1-hop neighbours "<<n1<<endl;
        //cout<< "2-hop neighbours "<<n2<<endl;
        //cout<< "totoal 2-hop neighbours "<<n1+n2<<endl;

    }
}



void BiCore::computeOrderedSubgraph(vector<int> & order, vector<vector<int>> &adjL, vector<vector<int>> &adjR,vector<V2N> &allV) {
    vector<int> vertexLookup(order.size());
    for(int i=0;i<order.size();i++){
        int u =  order[i];
        vertexLookup[u] = i;
    }

   // vector<pair<vector<int>,vector<int>>> sets;
    string filename = "../../sizes/sizesLG.txt";
    ofstream fout;

    //cout<<filename<<endl;
    fout.open(filename);
    vector<BFST> bfs(adjL.size());
    int lMax =0;
    int rMax = 0;
    for(int i=0;i<order.size();i++){
        //vector<pair<int,int>> edgeList; //the first vertex in the edge is in L, the second vertex in the edge is in R
        vector<int> L,R;
        int u = order[i];
        L.emplace_back(u);
        for(auto &v:adjL[u]){
            //edgeList.emplace_back(make_pair(u,v));
            R.emplace_back(v);
            for(auto &u_prime:adjR[v]){
                if(vertexLookup[u_prime]>i){
                    if(bfs[u_prime].v!=u) {
                        bfs[u_prime].v = u;
                        //edgeList.emplace_back(make_pair(u_prime, v));
                        L.emplace_back(u_prime);
                    }
                }
            }
        }
        //sets.emplace_back(make_pair(L,R));
        if(lMax<L.size()){
            lMax = L.size();
        }

        if(rMax<R.size()){
            rMax = R.size();
        }
        fout<<L.size()<<","<<R.size()<<endl;





    }
    cout<<lMax<<","<<rMax<<endl;

    //write the sets to file
    //cout<<"for test"<<endl;
    fout.close();
/**
    for(int i=0;i<sets.size();i++){
        ofstream fout;
        vector<int> L = sets[i].first;
        vector<int> R = sets[i].second;
        string filename = "../SLG/"+to_string(L[0])+".txt";
        //cout<<filename<<endl;
        fout.open(filename);



        for(auto &u:L){
            fout<<u<<",";
        }
        fout<<endl;
        for(auto &v:R){
            fout<<v<<",";
        }
        fout.close();
    }
**/
    //fout<<vL.size()<<","<<vR.size()<<endl;
    //for(int i=1;i<vL.size();i++){
    //    for(int j=1;j<vR.size();j++){
    //        if(j==vR.size()-1){
    //            fout<<adjM[i][j]<<"\n";
     //       }else{
     //           fout<<adjM[i][j]<<" ";
     //       }
     //   }
    //}
    //fout.close();
}


