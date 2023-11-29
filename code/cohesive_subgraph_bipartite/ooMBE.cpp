
#include "ooMBE.h"
#include "Utility.h"





ooMBE::~ooMBE() {
    //cout<<"ooMBE"<<endl;
    //cout<<++des_times<<endl;
    delete [] ln_for_vex;
    delete []  other_for_vex;
    delete [] two_hop;
    //delete [] pivots_arry;
    if(R!= nullptr) {delete [] R;R=nullptr;}

    //cout<<"$$$$$$1"<<endl;
    if(indexR!= nullptr){delete [] indexR;indexR=nullptr;}
    //cout<<"$$$$$$2"<<endl;
    if(adjIndexR!= nullptr) {delete [] adjIndexR;adjIndexR=nullptr;}
    if(L!= nullptr) {delete [] L;L = nullptr;}
    if(indexL!= nullptr) {delete [] indexL; indexL = nullptr;}
    if(adjIndexL!= nullptr){delete [] adjIndexL; adjIndexL = nullptr;}
    //cout<<"$$$$$$2"<<endl;
    if(counts!= nullptr){delete [] counts;counts = nullptr;}
    //cout<<"$$$$$$3"<<endl;

}

void ooMBE::adv_mbeStart() {
    // use an equation to determin which side to enumerate
    cout<<"call adv_mbeStart()"<<endl;
    if (g.max_degree_L <= g.max_degree_R) {
        // do Right
        //    g.switchLR();
        //    cout<<"switched!"<<endl;
    } else {
        //do left, but treat left as right
        g.switchLR();
        cout << "switched!" << endl;
    }




    //printVector(bcore.uniOrder);


    R = new int[g.adjR.size()];
    size_of_R = g.adjR.size();
    indexR = new int[g.adjR.size()];
    adjIndexR = new int[g.adjR.size()];

    L = new int[g.adjL.size()];
    size_of_L = g.adjL.size();
    indexL = new int[g.adjL.size()];
    adjIndexL = new int[g.adjL.size()];

    g_bfs.resize(size_of_R);
    isPovit.resize(size_of_R);

    fill(isPovit.begin(),isPovit.end(),false);

    counts = new int[g.adjR.size()];


    mbeDSini(g.adjL, g.adjR);
    //generate the unilateral order
    //TODO:this could be a potential problem
    BiCore bcore(g);
    bcore.e2hopCore(g.adjR, g.adjL);
    //end of problem

    if (maxPruning) {
        //Future work, further refine the unilateral order to maximize pruning
        //This optimisation is for sequential algorithm,
        vector<pair<int, int>> auD;
        int initUniCore = bcore.uni_l_core[bcore.uniOrder[1]];
        int count = 0;
        for (auto &v: bcore.uniOrder) {
            if (bcore.uni_l_core[v] == initUniCore) {
                count++;
            } else {
                auD.emplace_back(make_pair(initUniCore, count));
                initUniCore = bcore.uni_l_core[v];
                count = 1;
            }
        }

        auD.emplace_back(make_pair(initUniCore, count));

        int ini = 0;
        for (auto &pair: auD) { // adjust pairs for initial position of very unilateral core.
            int tmp = pair.second;
            pair.second = ini;
            ini = ini + tmp;
        }

        vector<vector<int>> vertex_with_same_uni_cores;
        for (int i = 0; i < auD.size(); i++) {
            pair<int, int> cp = auD[i];
            vector<int> ertex_with_same_uni_core;
            if (i < auD.size() - 1) {
                pair<int, int> np = auD[i + 1];
                for (int j = cp.second; j < np.second; j++) {
                    ertex_with_same_uni_core.emplace_back(bcore.uniOrder[j]);
                }
                vertex_with_same_uni_cores.emplace_back(ertex_with_same_uni_core);
            }
            if (i == auD.size() - 1) {
                for (int j = cp.second; j < bcore.uniOrder.size(); j++) {
                    ertex_with_same_uni_core.emplace_back(bcore.uniOrder[j]);
                }
                vertex_with_same_uni_cores.emplace_back(ertex_with_same_uni_core);
            }
        }
        bcore.uniOrder.clear();
        for (auto &vec: vertex_with_same_uni_cores) {
            std::sort(vec.begin(), vec.end(),
                      bind(&ooMBE::orderSort, this, std::placeholders::_1, std::placeholders::_2));
            bcore.uniOrder.insert(bcore.uniOrder.end(), vec.begin(), vec.end());
        }


        //cout<<"test tag"<<endl;//the above code has been tested
    }

    //printVector(bcore.uniOrder);
    //test
    /**
    if (bcore.uniOrder.size() != size_of_R) {
        cout << "Wrong!" << endl;
    }
     **/
    //adjust R, IndexR, but adjIndexR keeps the same
    vector<int> orderIndex(bcore.uniOrder.size());

    for (int i = 0; i < bcore.uniOrder.size(); i++) {
        int v = bcore.uniOrder[i];
        R[i] = v;
        indexR[v] = i;
        orderIndex[v] = i;
    }


    // let refine it now

    //compute order preserved batch pivots
    //refresh that give vertex v in R, indexR[v] returns the location of v in R, where R is an array. The same for L
    //parameters adjL, adjR, the order, adjIndexL, adjIndexR
    //note that the order is for R
    vector<int> pivots;
    vector<vector<int>> pivots_plus_neighbours;
    orderBatchPivots(pivots);


    int iP = 0;//g.adjR.size();
    int iQ = iP;
    int iB = g.adjR.size();


    int iA = 0;

    vector<pair<int, int>> RP;
    for (int i = iP; i < g.adjR.size(); i++) {
        RP.emplace_back(make_pair(i, 0));
    }


    vector<pair<int, int>> RQ, RB;

    vector<pair<int, int>> LA;
    for (int i = 0; i < g.adjL.size(); i++) {
        LA.emplace_back(make_pair(i, 0));
    }


    for (int i = 0; i < size_of_L; i++) {
        adjIndexL[i] = 0;
    }

    for (int i = 0; i < size_of_R; i++) {
        adjIndexR[i] = 0;
    }


    //cout<<"still fine!"<<endl;
    //for(int i=iP;i<size_of_R;i++){// When backtracking, R should be the same as here
    // This can be achieved by using the stored RP
    //for(auto &pa: RP){
    //for(int i=0;i<RP.size();i++){
    cout<<"max 2-hop core"<<bcore.maxE2Core<<endl;
    cout<<"breadth:"<<pivots.size()<<endl;
    int count =0;
    vector<BFST> bfs(size_of_R);

    int size_of_P;
    for (auto &v: pivots) {
        //cout<<"Progress:"<<++count<<endl;
        //move v to iB

        if(count==10000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==20000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==30000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==40000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==50000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==60000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==70000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==80000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==90000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==100000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==200000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==300000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==400000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==500000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==600000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==700000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==800000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==900000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==1000000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==1500000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==2000000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==2500000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==3100000){
            cout<<count<<":"<<nomb<<endl;
        }

        count++;
        //cout<<"Processing "<<v<<endl;
        int iP_prime = 0, iA_prime = 0, iB_prime = 0, iQ_prime = 0;
        iB_prime = iB - 1;
        if (indexR[v] != iB - 1) {
            int tmpv = R[iB - 1];
            R[iB - 1] = v;
            R[indexR[v]] = tmpv;
            indexR[tmpv] = indexR[v];
            indexR[v] = iB - 1;
        }

        vector<pair<int, int>> A_prime;
        vector<pair<int, int>> Q_prime;
        vector<pair<int, int>> P_prime;

        //move neighbours of v in L
        //   compute starting position of the neighbour
        int nV = g.adjR[v].size();
        iA_prime = size_of_L - nV;
        int cPos = size_of_L - 1;
        for (auto &u: g.adjR[v]) {
            if (indexL[u] != cPos) {
                int tmpu = L[cPos];
                L[cPos] = u;
                L[indexL[u]] = tmpu;
                //int tmpPos = indexL[u];
                indexL[tmpu] = indexL[u];
                indexL[u] = cPos;
            }
            cPos--;
        }
        //iQ_prime = iQ;
        //Put vertex in Q into the correct position

        //shrank P
        //let's do it with O(E)
        //*****  ******
        //the loop below runs in O(E)
        //vector<pair<int,int>> counts;

        //TODO problem starts here! Fixed!


        vector<int> v_2_hop;
        //fill(counts, counts + size_of_R, 0);
        for (int i = iA_prime; i < size_of_L; i++) {
            int u = L[i];
            //vector<int> neig_of_u = g.adjL[u];
            for(int j=0;j<g.adjL[u].size();j++){
                int v_prime = g.adjL[u][j];
            //for (auto &v_prime: neig_of_u) {
                if (orderIndex[v_prime] > orderIndex[v]) {
                    if(bfs[v_prime].v!=v){
                        bfs[v_prime].v=v;
                        bfs[v_prime].count=0;
                        v_2_hop.emplace_back(v_prime);
                    }
                    if (bfs[v_prime].v==v) {
                        bfs[v_prime].count++;
                    }

                }
            }
        }

        vector<int> s_p;
        for (auto &v_prime: v_2_hop) {
            //if(counts[v_prime]==g.adjR[v_prime].size()){ //this degree should be local degree
            if ( bfs[v_prime].count == size_of_L - iA_prime) {
                //move v_prime to B
                iB_prime = iB_prime - 1;
                if (indexR[v_prime] != iB_prime) {
                    //switch vertices at locations of iB_prime and indexR[v_prime]
                    int tmpv = R[iB_prime];
                    R[iB_prime] = v_prime;
                    R[indexR[v_prime]] = tmpv;
                    indexR[tmpv] = indexR[v_prime];
                    indexR[v_prime] = iB_prime;
                }
            } else {
                s_p.emplace_back(v_prime);
            }
            //reset
            bfs[v_prime].v=-1;
            bfs[v_prime].count=0;
        }

        nomb++;
level_accumulation_ += cur_level_;
max_level_ = std::max(cur_level_, max_level_);

        int iniPos = iB_prime - 1;
        for (auto &v_prime: s_p) {//align them to iB_prime
            if (indexR[v_prime] != iniPos) {
                int tmpv = R[iniPos];
                R[iniPos] = v_prime;
                R[indexR[v_prime]] = tmpv;
                indexR[tmpv] = indexR[v_prime];
                indexR[v_prime] = iniPos;
            }

            vector<int> localNeighbour;
            vector<int> other;
            for (int j = 0; j < g.adjR[v_prime].size(); j++) {  // adjust the neighbours of v_prime
                int u = g.adjR[v_prime][j];
                if (indexL[u] >= iA_prime) {
                    localNeighbour.emplace_back(u);
                }else{
                    other.emplace_back(u);
                }
            }

            int adj_ini = g.adjR[v_prime].size() - 1;
            for (auto & u:localNeighbour){
                g.adjR[v_prime][adj_ini] = u;
                adj_ini--;
            }
            adjIndexR[v_prime] = adj_ini + 1;
            for (auto & u:other){
                g.adjR[v_prime][adj_ini] = u;
                adj_ini--;
            }
            P_prime.emplace_back(make_pair(v_prime, adjIndexR[v_prime]));
            iniPos--;
        }
        iP_prime = ++iniPos;

        //Align old Q into the correct position
        int iniPos_Q = iP_prime - 1;
        for (auto &vqp: RQ) {
            int vq = vqp.first;
            if (indexR[vq] != iniPos_Q) {
                int tmpv = R[iniPos_Q];
                R[iniPos_Q] = vq;
                R[indexR[vq]] = tmpv;
                indexR[tmpv] = indexR[vq];
                indexR[vq] = iniPos_Q;
            }

            //adjust adj list for vq;
            /**
            int adj_ini = g.adjR[vq].size()-1;
            for(int j=0;j<g.adjR[vq].size();j++){
                int u = g.adjR[vq][j];
                if(j<=adj_ini){
                    if(indexL[u]>iA_prime){
                        int tmpu=  g.adjR[vq][adj_ini];
                        g.adjR[vq][adj_ini] = u;
                        g.adjR[vq][j] = tmpu;
                        adj_ini--;
                    }
                }
            }
            adjIndexR[vq] = adj_ini+1;
            Q_prime.emplace_back(make_pair(vq,adjIndexR[vq]));
             **/
            iniPos_Q--;
        }

        iQ_prime = iniPos_Q + 1;


        //now shrank Q to Q_prime
        vector<int> notQ;
        //for (int i = iQ_prime; i < iB_prime; i++) {
        for (int i = iQ_prime; i < iP_prime; i++) {
                int vq = R[i];
                int counts = 0;

                vector<int> localNeighbours;
                vector<int> other;
                for (int j = 0; j < g.adjR[vq].size(); j++) {
                    int u = g.adjR[vq][j];
                    if (indexL[u] >= iA_prime) {
                        counts++;
                        localNeighbours.emplace_back(u);
                    }else{
                        other.emplace_back(u);
                    }
                }

                    //vq.adj \cap A \ne emptyset
                    int adj_ini = g.adjR[vq].size() - 1;
                    for(auto u: localNeighbours){
                        g.adjR[vq][adj_ini]=u;
                        adj_ini--;
                    }
                    adjIndexR[vq] = adj_ini+1;
                    for(auto u: other){
                        g.adjR[vq][adj_ini]=u;
                        adj_ini--;
                    }
                    if(counts>0){
                        Q_prime.emplace_back(make_pair(vq,adjIndexR[vq]));
                    }



                if (counts == 0) {
                    notQ.emplace_back(vq);
                }
        }


        for (auto &vq: notQ) {
            int tmpv = R[iQ_prime];
            R[iQ_prime] = vq;
            R[indexR[vq]] = tmpv;
            indexR[tmpv] = indexR[vq];
            indexR[vq] = iQ_prime;
            iQ_prime = iQ_prime + 1;
        }

        //update adj index based on the current iQ, iP and iB
        for (int i = iA_prime; i < size_of_L; i++) {
            int u = L[i];

            vector<int> local_neighbours;
            vector<int> other;
            for (int j = 0; j < g.adjL[u].size(); j++) {
                int v_prime = g.adjL[u][j];
                if (indexR[v_prime] >= iQ_prime && indexR[v_prime] < iB_prime) {
                    local_neighbours.emplace_back(v_prime);
                }else{
                    other.emplace_back(v_prime);
                }
            }
            int adj_ini = g.adjL[u].size() - 1;
            for(auto &v:local_neighbours){
                g.adjL[u][adj_ini]=v;
                adj_ini--;
            }
            adjIndexL[u] = adj_ini + 1;
            for(auto &v:other){
                g.adjL[u][adj_ini]=v;
                adj_ini--;
            }
            A_prime.emplace_back(u, adjIndexL[u]);
        }

        // next is for the recursive call
//dvIMBEA(int iA,
//                     int iQ,
//                     int iP,
//                     int iB,
//                     vector<pair<int,int>> &A,
//                     vector<pair<int,int>> &P,
//                     vector<pair<int,int>> &Q)
        //if(P_prime.size()>0){
        /***
        cout<<"~~~~~~~~~~~"<<endl;
        for(int i=iB_prime;i<size_of_R;i++){
            cout<<R[i]<<",";
        }
        cout<<endl;
        //print L
        for(int i=iA_prime;i<size_of_L;i++){
            cout<<L[i]<<",";
        }
        cout<<endl;
        cout<<"~~~~~~~~~~~"<<endl;
        ***/
       // if(P_prime.size()>2501){
       //     cout<<"order does not work! "<<P_prime.size()<<endl;
       // }
        //if(Q_prime.size()==2){
        //    if(Q_prime[0].first==Q_prime[1].first){
        //        cout<<"Duplicated"<<endl;
        //    }
        //}



        if(P_prime.size()>0){
            //advIMBEA(iA_prime, iQ_prime, iP_prime, iB_prime, A_prime, P_prime, Q_prime);
            cur_level_++;
            advPIMBEA_local(iA_prime, iQ_prime, iP_prime, iB_prime, A_prime,
                            P_prime, Q_prime);
            cur_level_--;
        }

        if(true){
            int iB = iB - 1;
            //size_of_P = size_of_P-1;
            //iP = iB-;
            int last_pos_of_Q=count-1;
            //move v to position count
            /**
            if(indexR[v]!=last_pos_of_Q){
                int tmpv = R[last_pos_of_Q];
                R[last_pos_of_Q] = v;
                R[indexR[v]]=tmpv;
                indexR[tmpv] = indexR[v];
                indexR[v] = last_pos_of_Q;
            }
             **/
            iQ=0;
            RQ.emplace_back(make_pair(v,0));
            iP=last_pos_of_Q+1;
           //cout<<iP<<endl;

        }else{
            int ini_pos = iB - 1;
            //for(auto &pair:RP){
            int loc_of_v_in_RP;
            for (int i = 0; i < RP.size(); i++) {
                pair<int, int> p = RP[i];
                int bv = p.first;
                if (bv == v) {
                    loc_of_v_in_RP = i;
                }
                adjIndexR[bv] = p.second;
                if (indexR[bv] != ini_pos) {
                    int tmpv = R[ini_pos];
                    R[ini_pos] = bv;
                    R[indexR[bv]] = tmpv;
                    indexR[tmpv] = indexR[bv];
                    indexR[bv] = ini_pos;
                }
                ini_pos--;
            }
            //remove v form RP
            iP = ini_pos + 1;
            if(RP.size()>0){//this is always true
                if(RP.size()>=2){

                    pair<int, int> tmpp = RP[RP.size() - 1];
                    pair<int, int> tmpvp = RP[loc_of_v_in_RP];
                    RP[loc_of_v_in_RP] = tmpp;
                    RP.pop_back();
                    RQ.emplace_back(tmpvp);
                }else if (RP.size()==1) {
                    RQ.emplace_back(RP[0]);
                    RP.pop_back();
                }

                //refine array and indices
                int b_loc_of_v = indexR[v];
                R[b_loc_of_v] = R[iP];
                indexR[R[iP]] = b_loc_of_v;
                R[iP] = v;
                indexR[v] = iP;
                iP = iP + 1;
                int ini_pos_Q = iP - 1;
                for (auto &p: RQ) {
                    if (indexR[p.first] != ini_pos_Q) {
                        int tmpv = R[ini_pos_Q];
                        R[ini_pos_Q] = p.first;
                        R[indexR[p.first]] = tmpv;
                        indexR[tmpv] = indexR[p.first];
                        indexR[p.first] = ini_pos_Q;
                    }
                    ini_pos_Q--;
                    adjIndexR[p.first] = 0;
                }
                iQ = ini_pos_Q + 1;
                //if(iQ!=0){
                //    cout<<"Hehehehe"<<endl;
               // }
                //cout<<iQ<<endl;
                for(int i=iA_prime;i<size_of_L;i++){
                    int u = L[i];
                    adjIndexL[u]=0;
                }
                //fill(adjIndexL+(iA_prime-1), adjIndexL + size_of_L, 0);
            }

            cout<<iP<<endl;
        }

        //cout << "testint" << endl;
        //back track
        //align P,Q,B back first




        //for()
        //adj align back

        /***
        for(int j=i+1;j<RP.size();j++){
            int v_prime = RP[j].first;
            //int number_of_intersection = 0;
            int adj_v_start = adjIndexR[v_prime];
            int iniPos = g.adjR[v_prime].size()-1;
            for(int z=adj_v_start;z<g.adjR[v_prime].size();z++){
                int u = g.adjR[v_prime][z];
                if(indexL[u]>=iA_prime){
                    //move it to the iniPos
                    if(z<iniPos){
                        int tmpu =  g.adjR[v_prime][iniPos];
                        g.adjR[v_prime][iniPos] = u;
                        g.adjR[v_prime][z] = tmpu;

                    }
                    iniPos--;
                }
            }
            //iniPos++ is the adjIndexR[v_prime]
            int adj_v_prime_update = iniPos+1;
            adjIndexR[v_prime] = adj_v_prime_update;
            int intersection = g.adjR[v_prime].size()-adj_v_prime_update;
            if(intersection==g.adjR[v_prime].size()){
                //aggressive expansion
            }
        }
        ***/
        //cout<<"test"<<endl;

    }

    /**
     * vector<int> bfs(size_of_L,-1);
    for (auto & u:bcore.uniOrder){//loop 1
        cout<<u<<",";
        vector<pair<int,int>> A_prime;
        vector<pair<int,int>> Q_prime;
        vector<pair<int,int>> P_prime;
        //compute induced subgraph
        int loc_of_u=indexL[u];
        vertex_swap(loc_of_u,size_of_L-1,L,indexL);
        int iB_prime=size_of_L;
        //int iP_prime=
        //algin neighbours of u to A
        int ini_loc = size_of_R-1;
        for(auto &v:g.adjL[u]){
            int loc_of_v = indexR[v];
            if(loc_of_v!=ini_loc){
                vertex_swap(loc_of_v,ini_loc,R,indexR);
            }
            ini_loc--;
        }
        int iA_prime=ini_loc+1; //starting position of iA_prime
        //prepare P, i.e., vertices
        ini_loc = iB_prime-1;

        for(int loc=iA_prime;loc<size_of_R;loc++){
            int v = R[loc];
            for(auto &u_prime: g.adjR[v]){
                if(bcore.index_for_uniOrder[u_prime]>bcore.index_for_uniOrder[u]){// appears after u in the unilateral order
                    if(bfs[u_prime]!=u){
                        //first time it visits
                        bfs[u_prime] = u;
                        int loc_u_prime = indexL[u_prime];
                        if(loc_u_prime!=ini_loc){
                            vertex_swap(loc_u_prime,ini_loc,L,indexL);
                        }
                        ini_loc--;
                    }
                }

            }
        }
        int iP_prime = ini_loc+1;

        //adjust adjIndexR and adjust adjIndexL
        for(int loc=iA_prime;loc<size_of_R;loc++){
            int v_prime = R[loc];
            for(auto &u_prime:g.adjR[v_prime]){
                int loc_u_prime = indexL[u_prime];
                if(loc_u_prime>=iP_prime&&loc_u_prime<iB_prime){

                }
            }

        }

        for(int loc=iP_prime;loc<iB_prime;loc++){

        }


        int iB=size_of_L;


        //need

    }
    **/


    //if(R!= 0) {delete [] R;R=0;}


}
void ooMBE::mbeStart() {



    R = new int[g.adjR.size()];
    size_of_R = g.adjR.size();
    //fill(R,R+size_of_R,0);
    indexR = new int[g.adjR.size()];
    //fill(indexR,indexR+size_of_R,0);
    adjIndexR = new int[g.adjR.size()];
    //fill(adjIndexR,adjIndexR+size_of_R,0);

    L = new int[g.adjL.size()];
    size_of_L = g.adjL.size();
    //fill(L,L+size_of_L,0);
    indexL = new int[g.adjL.size()];
    //fill(indexL,indexL+size_of_L,0);
    adjIndexL = new int[g.adjL.size()];
    //fill(adjIndexL,adjIndexL+size_of_L,0);
    mbeDSini(g.adjR,g.adjL);

    //counts = new int[g.adjR.size()];


    BiCore bcore(g);
    bcore.e2hopCore(g.adjR, g.adjL);


    //
    int iP=0;//g.adjR.size();
    int iQ=0;
    int iB=g.adjR.size();


    int iA = 0;

    vector<pair<int,int>> RP;
    for(int i=0;i<g.adjR.size();i++){
        RP.emplace_back(make_pair(i,0));

    }


    vector<pair<int,int>> RQ,RB;

    vector<pair<int,int>> LA;
    for(int i=0;i<g.adjL.size();i++){
        LA.emplace_back(make_pair(i,0));
    }


    for(int i=0;i<size_of_L;i++){
        adjIndexL[i]=0;
    }

    for(int i=0;i<size_of_R;i++){
        adjIndexR[i]=0;
    }
    //iniLogical(g.adjR,g.adjL);
    advIMBEA(iA,iQ,iP,iB,LA,RP,RQ);
    //cout<<"******"<<nomb<<"******"<<endl;

}

void ooMBE::mbeDSini(vector<vector<int>> &adjL, vector<vector<int>> &adjR) {


    iniArrayValues(adjL,L,indexL,adjIndexL);

    iniArrayValues(adjR,R,indexR,adjIndexR);



    //cout<<adjR.size()<<endl;
    //cout<<adjL.size()<<endl;
}

// pair.first: vertex v, pair.second: loc_{v} for adjList of v
//TODO: several bugs to be fixed
void ooMBE::advIMBEA(int iA,
                     int iQ,
                     int iP,
                     int iB,
                     vector<pair<int,int>> &A,
                     vector<pair<int,int>> &P,
                     vector<pair<int,int>> &Q) {
    //cout<<"call recursive subproblem"<<endl;
    //special case:

    //vector<pair<int,int>> B_prime; is not necessary
    //vector<int> P(R+iP,R+iB-1);
    //for(auto &v:P){  //P could become P' later
    //this for loop can be changed iP to iB-1
    //for(int x=iP;x<iB;x++){
    for(int z=0;z<P.size();z++){
        pair<int,int> v= P[z];
        //int global_v = R[iP];
        //move v from the current location to iB+1 if necessary
        vector<pair<int,int>> A_prime;
        vector<pair<int,int>> Q_prime;
        vector<pair<int,int>> P_prime;
        int iQ_prime = iQ;


        int locV=indexR[v.first];
        /**
         * This is essentially moving a vertex to a known location
         * I
         */
        if(locV != iB-1){
            int vprime = R[iB-1];

            R[locV] = vprime;
            R[iB-1] = v.first;


            indexR[v.first] = iB-1;
            indexR[vprime] = locV;
            //cout<<"test"<<endl;

        }

        int iBprime = iB-1;

        //******refine A: works on set L and generate A' <- A\cap adj[v] ******
        //Since H is maintained, in H adj[v] must be a subset of A.
        //As such, adj[v]  \cap A is equivalent to align adj[v] to last few of positions of L

        int localAdjv = adjIndexR[v.first];
        int startLoc = g.adjL.size()-1;
        for(int i=localAdjv;i<g.adjR[v.first].size();i++){
            int u = g.adjR[v.first][i];
            int locU = indexL[u];
            if(indexL[u]!=startLoc){
                //move vertex to location of startLoc
                if(startLoc>size_of_L){
                    cout<<"Problem!"<<endl;
                }
                int uprime = L[startLoc];

                L[locU] = uprime;
                L[startLoc] =u;

                indexL[u]= startLoc;
                if(uprime>size_of_L){
                    cout<<" uprime  Problem!"<<endl;
                }
                indexL[uprime] = locU;
            }
            startLoc=startLoc-1;
        }


        int iAprime = startLoc+1;   //cout<<"haha test"<<endl; it is okay

        //*** adjust Q
        //first adjust Q  check maximality
        bool isMaximal = true;
  
        cout<<"call advimbea\n"<<endl;
        for(int i=iQ_prime;i<iP;i++){//loop on v\in Q
            //cout<<"--@@@---"<<endl;
            if(i<0){
                cout<<"730 problem <0"<<endl;
            }
            if(i>=size_of_R){
                cout<<"730 problem >size"<<endl;
            }
            int v = R[i];
            int insctv = 0;
            //loop on its local neighbours
            //aligan its neighbours from .. to ..
            int startLoc = g.adjR[v].size()-1;
            for(int j=adjIndexR[v];j<g.adjR[v].size();j++){//loop on the adj[v]
                    int u = g.adjR[v][j];
                    if(indexL[u]>=iAprime){
                        insctv++;
                    }
            }

            //TODO: update location for adjacency list of Q_prime
            /*
            if(insctv!=0){
                //startloc+1 is adj[v]
                Q_prime.emplace_back(make_pair(v,startLoc+1));
                adjIndexR[v]=startLoc+1;
                cout<<"--@@@---"<<endl;
            }
            */

            if(insctv==(size_of_L-iAprime)){
                //non maximal
                isMaximal = false;
                //cout<<"!!!!"<<endl;
                break;
            }

            if(insctv==0){
                //cout<<"--@@@---"<<endl;
                //move it out of Q
                int locV = indexR[v];
                if(locV!=iQ_prime){
                    int vprime = R[iQ_prime];

                    R[locV] = vprime;
                    R[iQ_prime] = v;

                    indexR[vprime] = locV;
                    indexR[v] = iQ_prime;
                }
                iQ_prime=iQ_prime+1;
            }
            //insctv=0;
        }//end of processing Q  Q' has been generated
        if(isMaximal==false){
          //cout<<"none maximal!"<<endl;
        }
        if(isMaximal){
            //now start to refine P according to the updated v, i.e., derive P_prime
            // Now vertices in P are located from iP to iBprime-1
            int iBprime_prime = iBprime;
            int iP_prime = iP;
            for(int i=iP_prime;i<iBprime_prime;i++){
                //compute induced subgraph
                int v = R[i];
                int local_adj_v_location = adjIndexR[v];
                int startLoc = g.adjR[v].size()-1;
                //int insctv = 0; insctv is not necessary, using startLoc is sufficient

                for(int j=local_adj_v_location;j<g.adjR[v].size();j++){//test if adj[v]\cap A_prime ne emptyset
                    int u = g.adjR[v][j];
                    //TODO: this loop should be refined and simplified  Done! figure crossed
                    //update
                    if(j<=startLoc){
                        if(indexL[u]>=iAprime){
                            //insctv++;
                            if(j==startLoc){
                                startLoc--;
                            }
                            if(j<startLoc){
                                //insctv++;
                                int uprime = g.adjR[v][startLoc];
                                while(indexL[uprime]>=iAprime&&uprime!=u){ //refine uprime if necessary
                                    //insctv++;
                                    startLoc--;
                                    uprime=g.adjR[v][startLoc];
                                }
                                //now
                                //switch u and uprime
                                if(j==startLoc){
                                    //insctv++;
                                    startLoc--;
                                    //break;
                                }else if(j<startLoc){
                                    g.adjR[v][startLoc]=u;
                                    g.adjR[v][j] = uprime;
                                    startLoc--;
                                }//else{
                                //    cout<<"wrong!"<<endl;
                                //}
                            }
                        }
                    }

                }
                //int test = size_of_L;

                int insctv = g.adjR[v].size()-(startLoc+1);
                if((insctv>0)&&(insctv<(size_of_L-iAprime))){
                    //adjust H' for v in P' and H' is used for the next recursion
                    int newadjv= startLoc+1;
                    P_prime.emplace_back(make_pair(v,newadjv));
                    //TODO: can be double checked
                    adjIndexR[v]=newadjv;
                    //cout<<"test hehehe---+++"<<endl;
                }
                if(insctv==0){
                    int locV = indexR[v];
                    //move v from P to Q, say, append v to the end of Q first and then move it to the beginning of Q
                    if(locV!=iP_prime){  // if v is not at the location iP_prime
                        //move v to the location of iP_prime
                        int v_prime = R[iP_prime];
                        R[locV] = v_prime;
                        R[iP_prime] = v;
                        indexR[v] = iP_prime;
                        indexR[v_prime] = locV;

                    }
                    iP_prime++;///Ahhhhhhhhhh
                    int loc_of_v = indexR[v];
                    if(loc_of_v!=iQ_prime){
                        int v_prime = R[iQ_prime];
                        R[loc_of_v] = v_prime;
                        R[iQ_prime] = v;
                        indexR[v_prime]=loc_of_v;
                        indexR[v]=iQ_prime;
                    }
                    iQ_prime++;///Ahhhhhhhhhhh
                }
                //aggressive
                if(insctv==(size_of_L-iAprime)){
                    int locV = indexR[v];
                    if(locV!=iBprime_prime-1){//if v is not at the location iBprime_prime-1
                        //move v to location iBprime_prime-1
                        int v_prime = R[iBprime_prime-1];
                        R[locV] = v_prime;
                        R[iBprime_prime-1] = v;
                        indexR[v]= iBprime_prime-1;
                        indexR[v_prime] = locV;
                        i--;//this is correct

                    }
                    iBprime_prime--;

                }
                //reset intersect
                //insctv=0;

            }
           // cout<<"test hehehe---"<<endl;
            //cout<<P_prime.size()<<endl;
            //report maximal biclique
            nomb++;
level_accumulation_ += cur_level_;
max_level_ = std::max(cur_level_, max_level_);
            //print maximal bicliques
            //print R
            /**
            cout<<"~~~~~~~~~~~"<<endl;
            for(int i=iBprime_prime;i<size_of_R;i++){
                cout<<R[i]<<",";
            }
            cout<<endl;
            //print L
            for(int i=iAprime;i<size_of_L;i++){
                cout<<L[i]<<",";
            }
            cout<<endl;
            cout<<"~~~~~~~~~~~"<<endl;
            **/
            if(P_prime.size()>0){
                //adjust adjlist for u  \in A  \in L
                for(int i=iAprime;i<size_of_L;i++){
                    int u = L[i];
                    int startLoc = g.adjL[u].size()-1;
                    //int intersec=0;
                    //TODO: logic adjustment
                    for(int j=adjIndexL[u];j<g.adjL[u].size();j++){
                        if(j<=startLoc){
                            int v = g.adjL[u][j];
                            if(indexR[v]>=iQ_prime && indexR[v]<=iBprime_prime-1){
                                //intersec++;
                                if(j==startLoc){
                                    startLoc--;
                                }
                                if(j<startLoc){
                                    int v_prime = g.adjL[u][startLoc];
                                    //if(j<startLoc){
                                    while(indexR[v_prime]>=iQ_prime && indexR[v_prime]<=iBprime_prime-1 &&v_prime!=v){// check if v_prime is also in P\cup Q if yes, keep it here
                                        startLoc--;
                                        //intersec++;
                                        v_prime = g.adjL[u][startLoc];
                                    }
                                    //if(startLoc>j){//could be removed
                                    if(j==startLoc){
                                        startLoc--;
                                    }else if(j<startLoc){
                                        g.adjL[u][startLoc] = v;
                                        g.adjL[u][j] = v_prime;
                                        startLoc--;
                                    }else{
                                        cout<<"wrong!"<<endl;
                                    }
                                }

                            }
                        }

                    }
                    //startLoc+1 is adjL for u now
                    //if(intersec>0){
                    A_prime.emplace_back(make_pair(u,startLoc+1));
                    adjIndexL[u]=startLoc+1;
                    //cout<<"test hehehe+++"<<endl;
                    //}

                }

                //adjust adj_list for Q
                //TODO: for all the adjlist adjustment, double-check is needed figure crossed
                //TODO: intersec is not necessary and should be simplified
                for(int i=iQ_prime;i<iP_prime;i++){
                    int v = R[iQ_prime];
                    int startLoc = g.adjR[v].size()-1;
                    //int intersec = 0;
                    for(int j=adjIndexR[v];j<g.adjR[v].size();j++){
                        int u = g.adjR[v][j];
                        if(j<=startLoc){
                            if(indexL[u]>=iAprime){ //v or u??
                                //intersec++;
                                if(j==startLoc){//do not need to move
                                    startLoc--;
                                }
                                if(j<startLoc){
                                    //TODO:5_Feb
                                    int u_prime = g.adjR[v][startLoc]; //find the first u_prime that is not in Aprime
                                    while(indexL[u_prime]>=iAprime&&u_prime!=u){
                                        startLoc--;
                                        //intersec++;
                                        u_prime = g.adjR[v][startLoc];
                                    }
                                    //if(startLoc>j){//
                                    if(j==startLoc){
                                        startLoc--;
                                    }else if(j<startLoc){
                                        g.adjR[v][startLoc]=u;
                                        g.adjR[v][j] = u_prime;
                                        startLoc--;
                                    }else{
                                        cout<<"wrong!"<<endl;
                                    }
                                }
                            }
                        }
                    }
                    //if(intersec>0){
                        Q_prime.emplace_back(make_pair(v,startLoc+1));
                        adjIndexR[v]=startLoc+1;
                    //}
                }//end of adjusting adjlists for Q
                //void ooMBE::advIMBEA(int iA, int iQ, int iP, int iB, vector<pair<int,int>> &A,vector<pair<int,int>> &P,vector<pair<int,int>> &Q) {
                if(false){
                    advIMBEA(iAprime,iQ_prime,iP_prime,iBprime_prime,A_prime,P_prime,Q_prime);
                }


            }


        }// is maximal finish
        //back tracking adjustment
        //use the stored  iA, iQ, iP, iB, A,P,Q to restore the logical partitions
        //align P to iP to iB-1 and adjust the location index
        //remove v from P
        //Add v to Q
        int loc_pointer = iB-1;
        //for(auto &v_adjLoc:P){
        int v_before_rec = P[z].first;
        if(v_before_rec!=v.first){
            cout<<"error mismatch!"<<endl;
        }
        adjIndexR[v.first]=P[z].second;
        for (int i=z+1;i<P.size();i++){//equivalent to moving v out
            auto v_adjLoc = P[i];
            int v=v_adjLoc.first;
            int adj_index_value = v_adjLoc.second;
            adjIndexR[v]=adj_index_value; //above okay

            if(loc_pointer>=size_of_R){
                cout<<"error!+++++ 1018"<<endl;
            }
            if(loc_pointer<0){
                cout<<"error!----- 1018"<<endl;
            }
            int v_prime = R[loc_pointer]; //the vertex at the pointer
            int loc_of_v = indexR[v];
            //if(loc_of_v<iP){ // loc_of_v may at iB, this vertex will be dealt with later
            //TODO:6_Feb
            //switch v and v_prime
            R[loc_of_v] = v_prime;
            R[loc_pointer] = v;
            //adjust locations for v and v_prime
            indexR[v] = loc_pointer;
            indexR[v_prime] = loc_of_v;
            loc_pointer--;
            //}
        }
        //iP=iP+1;
        iP=loc_pointer+1;
        //TODO:double check here!
        Q.emplace_back(v);
        loc_pointer = iP-1;
        for(auto &v_adjLoc:Q){
            int v = v_adjLoc.first;
            int adj_index_value = v_adjLoc.second;
            adjIndexR[v] = adj_index_value;
            if(loc_pointer>=size_of_R){
                cout<<"error++++! 1046"<<endl;
            }
            if(loc_pointer<0){
                cout<<"error----! 1046"<<endl;
            }
            int v_prime = R[loc_pointer];
            int loc_of_v = indexR[v];
            R[loc_of_v] = v_prime;
            R[loc_pointer] = v;
            //adjust locations for v and v_prime
            indexR[v] = loc_pointer;
            indexR[v_prime] = loc_of_v;
            loc_pointer--;
        }
        iQ = loc_pointer+1;
        //align A back to iA to the end of L
        loc_pointer=size_of_L-1;
        for(auto &u_adjLoc:A){
            adjIndexL[u_adjLoc.first] = u_adjLoc.second;

            int u=u_adjLoc.first;
            int loc_of_u = indexL[u];
            int u_prime = L[loc_pointer];
            //switch u and u_prime and update location information
            L[loc_pointer] = u;
            L[loc_of_u]=u_prime;
            indexL[u] = loc_pointer;
            indexL[u_prime]=loc_of_u;
            loc_pointer--;
        }

        // adjust adjlist back

        //TODO, initially, iQ is not well managed, need to refine!!!!! done !
        //move v to Q,
        //int loc_of_v = indexR[v.first];
        //if(loc_of_v!=iP){ // in most cases, this condition cannot be satisfied
        //    int v_prime = R[iP];
        //    R[iP] = v.first;
        //    R[loc_of_v]= v_prime;
        //    indexR[v.first]=iP;
        //    indexR[v_prime]=loc_of_v;
        //}
        //iP=iP+1;


        //1:recover adjlist

        //2:recover R

        //3:recover L

        //normal biclique operation

        //remove v from P to Q

    }

    //vector<int> PAdjLoc(R+iP,R+iB-1);


    //for(int i=)

    //cout<<"test";
    //vector<int> P()


}



void ooMBE::iniVector(vector<vector<int>> &adjX, vector<int> &X, vector<int> &indexX,vector<int> &adjIndexX) {


    for(int i=0;i<adjX.size();i++){
        X.emplace_back(i);
        indexX.emplace_back(i);
    }

    for(int i=0;i<X.size();i++){
        indexX[X[i]]=i;
        adjIndexX[X[i]]=0;   //initially, every vertex in the adjlist is considered.
    }



}

void ooMBE::iniArrayValues(vector<vector<int>> &adjX, int *X, int *indexX, int *adjIndexX) {

    for(int i=0;i<adjX.size();i++){
        X[i]=i;
        indexX[X[i]]=i;
        adjIndexX[i] = 0;
    }

}

//asume work on the right side




void ooMBE::orderBatchPivots_neighbours(vector<PivotsNeighbour> & pivots) {
    //vector<BFST> counts (g.adjR.size());
    vector<int> isDominated (g.adjR.size());
    vector<BFST> bfs (g.adjR.size());
    //vector<int> pivots;
    vector<int> pruned;
    int n_of_pivots=0;
    int pid_in_pvoits=0;
    for(int i=0;i<size_of_R;i++){
        int cV = R[i]; // the reference vertex
        //auto rightB = g.CommonNeighborsL(g.adjR[cV]);
        if(isDominated[cV]!=1){
            processing_nodes++;
            PivotsNeighbour p;
            p.vertex=cV;
           // n_of_pivots++;
            int locV = indexR[cV];
            vector<int> visitedVs;
            for(auto &u:g.adjR[cV]){ // the neighbours of cV, they are in L
                for(auto &v_prime:g.adjL[u]){
                    if(v_prime!=cV&&indexR[v_prime]>locV){ //in fact the first condition is not necessary. It's for readability.
                        //visitedVs.emplace_back(v_prime);
                        if(bfs[v_prime].v!=cV){
                            bfs[v_prime].v = cV;
                            bfs[v_prime].count=1;
                            visitedVs.emplace_back(v_prime);
                            //p.adj.emplace_back(v_prime);
                        }else{
                            bfs[v_prime].count++;
                        }
                    }else if (v_prime!=cV&&indexR[v_prime]<locV){
                        if(bfs[v_prime].v!=cV){
                            bfs[v_prime].v = cV;
                            if(isPovit[v_prime]==true){
                                p.Q.emplace_back(v_prime);
                            }
                        }
                    }
                }
            }

            //if()



            for(auto &vv: visitedVs){//vv: a visited vertex
                if(bfs[vv].count==g.adjR[vv].size()){
                    //cout<<"pruned"<<endl;
                    isDominated[vv]=1;
                    pruned.emplace_back(vv);

                }
                if(bfs[vv].count==g.adjR[cV].size()){//lixu
                     p.agg.emplace_back(vv);
                }
                if(bfs[vv].count>0&&(bfs[vv].count<g.adjR[cV].size())){
                    p.adj.emplace_back(vv);
                }
            }
            /*if(p.agg.size() != rightB.size() - 1)
            {
              std::cout<<cV<<":"<<std::endl;
              std::cout<<p.agg.size()<<" "<<rightB.size()<<std::endl;
              for(int rb: rightB)
              {
                if(rb != cV)
                {
                  std::cout<<rb<<": "<<bfs[rb].count<<" "<<g.adjR[cV].size()<<std::endl;
                }
              }
            }*/
            isPovit[p.vertex]=true;
            pivots.emplace_back(p);
        }

    }
/***
    int ini=0;
    //need to update indexR as well
    for(auto &v:pruned){
        R[ini]=v;
        indexR[v] = ini;
        ini++;
    }
    int start_of_p = ini;
    for(auto &v: pivots){
        R[ini]=v;
        indexR[v]=ini;
        ini++;
    }
    //rearrange R: put vertices
    //int iniP = size_of_R-1;


    //return locations

    //cout<<"test"<<endl; //the test for the above code finishes
  ****/
    //return start_of_p;
}

//work on R,
void ooMBE::orderBatchPivots(vector<int> & pivots) {
    vector<BFST> bfs (g.adjR.size());
    //vector<BFST> counts (g.adjR.size());
    vector<int> isDominated (g.adjR.size());
    //vector<int> pivots;
    vector<int> pruned;
    int n_of_pivots=0;
    for(int i=0;i<size_of_R;i++){
        int cV = R[i]; // the reference vertex
        if(isDominated[cV]!=1){
            pivots.emplace_back(cV);
            n_of_pivots++;
            int locV = indexR[cV];
            vector<int> visitedVs;
            for(auto &u:g.adjR[cV]){ // the neighbours of cV, they are in L
                for(auto &v_prime:g.adjL[u]){
                    if(v_prime!=cV&&indexR[v_prime]>locV){ //in fact the first condition is not necessary. It's for readability.
                        //visitedVs.emplace_back(v_prime);
                        if(bfs[v_prime].v!=cV){
                            bfs[v_prime].v = cV;
                            bfs[v_prime].count=1;
                            visitedVs.emplace_back(v_prime);
                        }else{
                            bfs[v_prime].count++;
                        }
                    }
                }
            }
            for(auto &vv: visitedVs){//vv: a visited vertex
                if(bfs[vv].count==g.adjR[vv].size()){
                    //cout<<"pruned"<<endl;
                    isDominated[vv]=1;
                    pruned.emplace_back(vv);
                }
            }
        }

    }
/***
    int ini=0;
    //need to update indexR as well
    for(auto &v:pruned){
        R[ini]=v;
        indexR[v] = ini;
        ini++;
    }
    int start_of_p = ini;
    for(auto &v: pivots){
        R[ini]=v;
        indexR[v]=ini;
        ini++;
    }
    //rearrange R: put vertices
    //int iniP = size_of_R-1;


    //return locations

    //cout<<"test"<<endl; //the test for the above code finishes
  ****/
    //return start_of_p;
}


void
ooMBE::adjustment_after_moving_a_vertex(int u,
                                        int *L,
                                        int *indexL,
                                        int *adjIndexL,
                                        int size_of_L,
                                        int *R,
                                        int *indexR,
                                        int *adjIndexR,
                                        int size_of_R) {

}

void ooMBE::vertex_swap(int loc1, int loc2, int* array ,int* index) {
    //update values
    int tmp = array[loc1];
    array[loc1] = array[loc2];
    array[loc2] = tmp;

    //update indices
    index[array[loc1]] = loc1;
    index[array[loc2]] = loc2;

}

void ooMBE::printVector(vector<int> o) {

    for(auto & v:o){
        cout<<v<<",";
    }
    cout<<endl;
}

void ooMBE::advPIMBEA(int iA,
                 int iQ,
                 int iP,
                 int iB,
                 vector<pair<int,int>> &A,
                 vector<pair<int,int>> &P,
                 vector<pair<int,int>> &Q) {
        //cout<<"call recursive subproblem"<<endl;
        //special case:

        //vector<pair<int,int>> B_prime; is not necessary
        //vector<int> P(R+iP,R+iB-1);
        //for(auto &v:P){  //P could become P' later
        //this for loop can be changed iP to iB-1
        //for(int x=iP;x<iB;x++){

    //vector<int> pivots;
    //vector<BFST> bfs(size_of_R);//using unordered
    //reset counting data structure

    //cout<<recursive_tracker<<endl;

    //it should be done with P

    for(int i=iQ;i<iB;i++){
        int v=R[i];
        g_bfs[v].v=-1;
        g_bfs[v].count=0;
        g_bfs[v].isV=false;
        isPovit[v]=true; //assuming all are pivots
    }
    for(int i=iQ;i<iB;i++){
        int v = R[i];
        //fill(counts, counts + size_of_R, 0);//this is O(n) but practically faster
        if(isPovit[v]==true){
            int adj_start=adjIndexR[v];
            vector<int> visited_2_hop;

            for(int j=adj_start;j<g.adjR[v].size();j++){
                int u=g.adjR[v][j];
                int ajd_u_start = adjIndexL[u];
                for(int k=ajd_u_start;k<g.adjL[u].size();k++){
                    int v_prime = g.adjL[u][k];
                    if(v!=v_prime){

                        if(indexR[v_prime]<iQ||indexR[v_prime]>=iB){
                            cout<<"wrong!"<<endl;
                            cout<<"loc if v_prime: "<<indexR[v_prime]<<endl;
                            cout<<"iQ: "<<iQ<<endl;
                            cout<<"iB: "<<iB<<endl;
                        }

                        if(g_bfs[v_prime].v!=v){
                            g_bfs[v_prime].v =v;
                            visited_2_hop.emplace_back(v_prime);
                        }
                        if(g_bfs[v_prime].v==v){
                            g_bfs[v_prime].count=g_bfs[v_prime].count+1;
                        }
                    }
                }
            }
            for(auto &v_prime:visited_2_hop){
                int loc_deg_v_prime = g.adjR[v_prime].size()-adjIndexR[v_prime];
                if(g_bfs[v_prime].count>g.adjR[v_prime].size()){
                    cout<<"2-hop based pruning very wrong!"<<endl;
                }
                if(g_bfs[v_prime].count==loc_deg_v_prime){
                    //v_prime is not a pivot
                    isPovit[v_prime]=false;
                    //cout<<"pruned!"<<endl;
                }
                g_bfs[v_prime].count =0;
            }
        }

    }

    //recursive_tracker++;

    vector<pair<int,int>> pivots;
    for(int i=0;i<P.size();i++){

        auto v_pair = P[i];
        int v = v_pair.first;
        if(isPovit[v]==true&&indexR[v]>=iP){ //the second condition is very important
            //Pivot p;
            //p.pair=v_pair;
            //p.loc = i;
            pivots.emplace_back(P[i]);
        }
    }

    //int pivot_starting = iniPos+1;

    //now, pivots should be aligned to the last few positions of P

    //TODO: neighbourhood adjustment has serious problems in this loop.
        for(int z=0;z<pivots.size();z++){
            pair<int,int> v= pivots[z];
            //int global_v = R[iP];
            //move v from the current location to iB+1 if necessary
            vector<pair<int,int>> A_prime;
            vector<pair<int,int>> Q_prime;
            vector<pair<int,int>> P_prime;
            int iQ_prime = iQ;


            int locV=indexR[v.first];
            /**
             * This is essentially moving a vertex to a known location
             * I
             * TODO:up to here dec
             */
            if(locV != iB-1){
                int vprime = R[iB-1];

                R[locV] = vprime;
                R[iB-1] = v.first;


                indexR[v.first] = iB-1;
                indexR[vprime] = locV;
                //cout<<"test"<<endl;

            }

            int iBprime = iB-1;

            //******refine A: works on set L and generate A' <- A\cap adj[v] ******
            //Since H is maintained, in H adj[v] must be a subset of A.
            //As such, adj[v]  \cap A is equivalent to align adj[v] to last few of positions of L

            int localAdjv = adjIndexR[v.first];
            int startLoc = g.adjL.size()-1;
            for(int i=localAdjv;i<g.adjR[v.first].size();i++){
                int u = g.adjR[v.first][i];
                int locU = indexL[u];
                if(indexL[u]!=startLoc){
                    //move vertex to location of startLoc
                    if(startLoc>size_of_L){
                        cout<<"Problem!"<<endl;
                    }
                    int uprime = L[startLoc];

                    L[locU] = uprime;
                    L[startLoc] =u;

                    indexL[u]= startLoc;
                    if(uprime>size_of_L){
                        cout<<" uprime  Problem!"<<endl;
                    }
                    indexL[uprime] = locU;
                }
                startLoc=startLoc-1;
            }


            int iAprime = startLoc+1;   //cout<<"haha test"<<endl; it is okay

            //*** adjust Q
            //first adjust Q  check maximality
            bool isMaximal = true;

            for(int i=iQ_prime;i<iP;i++){//loop on v\in Q
                //cout<<"--@@@---"<<endl;
                //if(i<0){
                //    cout<<"730 problem <0"<<endl;
                //}
                //if(i>=size_of_R){
                //    cout<<"730 problem >size"<<endl;
                //}
                int v = R[i];
                int insctv = 0;
                //loop on its local neighbours
                //aligan its neighbours from .. to ..

                //vector<int> localNeighbours_of_v;
                vector<int> localNeighbours_of_v;
                vector<int> other;
                for(int j=adjIndexR[v];j<g.adjR[v].size();j++){//loop on the adj[v]
                    int u = g.adjR[v][j];
                    if(indexL[u]>=iAprime){
                        insctv++;
                        localNeighbours_of_v.emplace_back(u);
                    }else{
                        other.emplace_back(u);
                    }
                }

                //TODO: update location for adjacency list of Q_prime
                /*
                if(insctv!=0){
                    //startloc+1 is adj[v]
                    Q_prime.emplace_back(make_pair(v,startLoc+1));
                    adjIndexR[v]=startLoc+1;
                    //cout<<"--@@@---"<<endl;
                }
                */
                int startLoc = g.adjR[v].size()-1;
                for(auto& u:localNeighbours_of_v){
                    g.adjR[v][startLoc] = u;
                    startLoc--;
                }
                adjIndexR[v]= startLoc+1;
                for(auto& u:other){
                    g.adjR[v][startLoc] = u;
                    startLoc--;
                }
                /**
                if(insctv==(size_of_L-iAprime)){
                    //non maximal
                    isMaximal = false;
                    cout<<"After refinement, this should be impossible !!!!"<<endl;
                //    break;
               }
                **/

                if(insctv==0){
                   // cout<<"--@@@---"<<endl;
                    //move it out of Q
                    int locV = indexR[v];
                    if(locV!=iQ_prime){
                        int vprime = R[iQ_prime];

                        R[locV] = vprime;
                        R[iQ_prime] = v;

                        indexR[vprime] = locV;
                        indexR[v] = iQ_prime;
                    }
                    iQ_prime=iQ_prime+1;
                }
                //insctv=0;
            }//end of processing Q  Q' has been generated
            //if(isMaximal==false){
                //cout<<"none maximal!"<<endl;
            //    cout<<"Should never gonna be here!"<<endl;
          // }
            if(isMaximal){
                //now start to refine P according to the updated v, i.e., derive P_prime
                // Now vertices in P are located from iP to iBprime-1
                int iBprime_prime = iBprime;
                int iP_prime = iP;
                for(int i=iP_prime;i<iBprime_prime;i++){
                    //compute induced subgraph
                    int v = R[i];
                    int local_adj_v_location = adjIndexR[v];

                    //int insctv = 0; insctv is not necessary, using startLoc is sufficient
                    int insctv=0;
                    vector<int> local_neighbours;
                    vector<int> other;
                    for(int j=local_adj_v_location;j<g.adjR[v].size();j++){
                        int u = g.adjR[v][j];
                        if(indexL[u]>=iAprime){
                            insctv++;
                            local_neighbours.emplace_back(u);
                        }else{
                            other.emplace_back(u);
                        }
                    }
                    int startLoc = g.adjR[v].size()-1;
                    for(auto& u:local_neighbours){
                        g.adjR[v][startLoc]=u;
                        startLoc--;
                    }
                    adjIndexR[v]=startLoc+1;
                    for(auto& u:other){
                        g.adjR[v][startLoc]=u;
                        startLoc--;
                    }
                    /***
                    for(int j=local_adj_v_location;j<g.adjR[v].size();j++){//test if adj[v]\cap A_prime ne emptyset
                        int u = g.adjR[v][j];
                        //TODO: this loop should be refined and simplified  Done! figure crossed
                        //update
                        if(j<=startLoc){
                            if(indexL[u]>=iAprime){
                                //insctv++;
                                if(j==startLoc){
                                    startLoc--;
                                }
                                if(j<startLoc){
                                    //insctv++;
                                    int uprime = g.adjR[v][startLoc];
                                    while(indexL[uprime]>=iAprime&&uprime!=u){ //refine uprime if necessary
                                        //insctv++;
                                        startLoc--;
                                        uprime=g.adjR[v][startLoc];
                                    }
                                    //now
                                    //switch u and uprime
                                    if(j==startLoc){
                                        //insctv++;
                                        startLoc--;
                                        //break;
                                    }else if(j<startLoc){
                                        g.adjR[v][startLoc]=u;
                                        g.adjR[v][j] = uprime;
                                        startLoc--;
                                    }//else{
                                    //    cout<<"wrong!"<<endl;
                                    //}
                                }
                            }
                        }

                    }
                     ***/
                    //int test = size_of_L;

                    //int insctv = g.adjR[v].size()-(startLoc+1);
                    if((insctv>0)&&(insctv<(size_of_L-iAprime))){
                        //adjust H' for v in P' and H' is used for the next recursion
                        //int newadjv= startLoc+1;
                        P_prime.emplace_back(make_pair(v,adjIndexR[v]));
                        //TODO: can be double checked
                        //djIndexR[v]=newadjv;
                        //cout<<"test hehehe---+++"<<endl;
                    }
                    if(insctv==0){
                        int locV = indexR[v];
                        //move v from P to Q, say, append v to the end of Q first and then move it to the beginning of Q
                        if(locV!=iP_prime){  // if v is not at the location iP_prime
                            //move v to the location of iP_prime
                            int v_prime = R[iP_prime];
                            R[locV] = v_prime;
                            R[iP_prime] = v;
                            indexR[v] = iP_prime;
                            indexR[v_prime] = locV;

                        }
                        iP_prime++;///Ahhhhhhhhhh
                        int loc_of_v = indexR[v];
                        if(loc_of_v!=iQ_prime){
                            int v_prime = R[iQ_prime];
                            R[loc_of_v] = v_prime;
                            R[iQ_prime] = v;
                            indexR[v_prime]=loc_of_v;
                            indexR[v]=iQ_prime;
                        }
                        iQ_prime++;///Ahhhhhhhhhhh
                    }
                    //aggressive
                    if(insctv==(size_of_L-iAprime)){
                        int locV = indexR[v];
                        if(locV!=iBprime_prime-1){//if v is not at the location iBprime_prime-1
                            //move v to location iBprime_prime-1
                            int v_prime = R[iBprime_prime-1];
                            R[locV] = v_prime;
                            R[iBprime_prime-1] = v;
                            indexR[v]= iBprime_prime-1;
                            indexR[v_prime] = locV;
                            i--;//this is correct

                        }
                        iBprime_prime--;

                    }
                    //reset intersect
                    //insctv=0;

                }
                // cout<<"test hehehe---"<<endl;
                //cout<<P_prime.size()<<endl;
                //report maximal biclique
                nomb++;
                level_accumulation_ += cur_level_;
                max_level_ = std::max(cur_level_, max_level_);
                //print maximal bicliques
                //print R
                /***
                cout<<"~~~~~~~~~~~"<<endl;
                for(int i=iBprime_prime;i<size_of_R;i++){
                    cout<<R[i]<<",";
                }
                cout<<endl;
                //print L
                for(int i=iAprime;i<size_of_L;i++){
                    cout<<L[i]<<",";
                }
                cout<<endl;
                cout<<"~~~~~~~~~~~"<<endl;
                ***/

                    //adjust adjlist for u  \in A  \in L
                    for(int i=iAprime;i<size_of_L;i++){
                        int u = L[i];
                        vector<int> local_neighbours;
                        vector<int> other;
                        for(int j=adjIndexL[u];j<g.adjL[u].size();j++){
                            int v = g.adjL[u][j];
                            if(indexR[v]>=iQ_prime && indexR[v]<iBprime_prime){
                                local_neighbours.emplace_back(v);
                            }else{
                                other.emplace_back(v);
                            }
                        }
                        int startLoc = g.adjL[u].size()-1;
                        for(auto &v:local_neighbours){
                            g.adjL[u][startLoc] = v;
                            startLoc--;
                        }
                        adjIndexL[u] = startLoc+1;
                        //int intersec=0;

                        for(auto &v:other){
                            g.adjL[u][startLoc] = v;
                            startLoc--;
                        }
                        A_prime.emplace_back(make_pair(u,adjIndexL[u]));
                        //adjIndexL[u]=startLoc+1;
                        //cout<<"test hehehe+++"<<endl;
                        //}

                    }

                    //adjust adj_list for Q
                    //TODO: for all the adjlist adjustment, double-check is needed figure crossed
                    //TODO: intersec is not necessary and should be simplified
                    for(int i=iQ_prime;i<iP_prime;i++){
                        int v = R[i];

                        Q_prime.emplace_back(make_pair(v,adjIndexR[v]));
                        //adjIndexR[v]=startLoc+1;
                        //}
                    }//end of adjusting adjlists for Q

                    //void ooMBE::advIMBEA(int iA, int iQ, int iP, int iB, vector<pair<int,int>> &A,vector<pair<int,int>> &P,vector<pair<int,int>> &Q) {
                    //if(P_prime.size()>2501){
                    //cout<<"order does not work! "<<P_prime.size()<<endl;

                     //}
                    //if(Q_prime.size()==2){
                    //    if(Q_prime[0].first==Q_prime[1].first){
                    //        cout<<"Duplicated"<<endl;
                    //    }
                    //}
                    if(P_prime.size()>0){
                    advPIMBEA(iAprime,iQ_prime,iP_prime,iBprime_prime,A_prime,P_prime,Q_prime);

                    //advIMBEA(iAprime,iQ_prime,iP_prime,iBprime_prime,A_prime,P_prime,Q_prime);

                }


            }// is maximal finish
            //back tracking adjustment
            //use the stored  iA, iQ, iP, iB, A,P,Q to restore the logical partitions
            //align P to iP to iB-1 and adjust the location index
            //remove v from P
            //Add v to Q


            //for(auto &v_adjLoc:P){
            int v_before_rec = pivots[z].first;
            //if(v_before_rec!=v.first){
            //    cout<<"error mismatch!"<<endl;
            //}
            adjIndexR[v.first]=pivots[z].second;
            /*
             * TODO:int i=z+1;i<P.size();i++ was wrong
             * vertices from z+1 to P.size()-1 are pivots only, all vertices in P except for  v.first
             * should be aligned back.
             */
            int loc_pointer = iB-1;
            int loc_of_to_be_remove=-1;
            for (int i=0;i<P.size();i++){//equivalent to moving v out/
                auto v_adjLoc = P[i];
                int v=v_adjLoc.first;
                int adj_index_value = v_adjLoc.second;
                adjIndexR[v]=adj_index_value; //above okay
                if(v==v_before_rec){
                    loc_of_to_be_remove=i;
                }
                if(v!=v_before_rec){
                    //if(loc_pointer>=size_of_R){
                    //    cout<<"error!+++++ 1018"<<endl;
                    //}
                    //if(loc_pointer<0){
                    //    cout<<"error!----- 1018"<<endl;
                    //}
                    int v_prime = R[loc_pointer]; //the vertex at the pointer
                    int loc_of_v = indexR[v];
                    //if(loc_of_v<iP){ // loc_of_v may at iB, this vertex will be dealt with later
                    //TODO:6_Feb
                    //switch v and v_prime
                    R[loc_of_v] = v_prime;
                    R[loc_pointer] = v;
                    //adjust locations for v and v_prime
                    indexR[v] = loc_pointer;
                    indexR[v_prime] = loc_of_v;
                    loc_pointer--;
                }
            }

            if(loc_of_to_be_remove!=P.size()-1){
                auto tmp = P[P.size()-1];
                P[P.size()-1]=v;
                P[loc_of_to_be_remove] = tmp;
            }

            P.pop_back();
            //iP=iP+1;
            iP=iP+1;
            //TODO:double check here!

            Q.emplace_back(v);
            loc_pointer = iP-1;
            for(auto &v_adjLoc:Q){
                int v = v_adjLoc.first;
                int adj_index_value = v_adjLoc.second;
                adjIndexR[v] = adj_index_value;
                //if(loc_pointer>=size_of_R){
                //    cout<<"error++++! 1046"<<endl;
                //}
                //if(loc_pointer<0){
                //    cout<<"error----! 1046"<<endl;
                //}
                int v_prime = R[loc_pointer];
                int loc_of_v = indexR[v];
                R[loc_of_v] = v_prime;
                R[loc_pointer] = v;
                //adjust locations for v and v_prime
                indexR[v] = loc_pointer;
                indexR[v_prime] = loc_of_v;
                loc_pointer--;
            }
            iQ = loc_pointer+1;
            //align A back to iA to the end of L
            loc_pointer=size_of_L-1;
            for(auto &u_adjLoc:A){
                adjIndexL[u_adjLoc.first] = u_adjLoc.second;

                int u=u_adjLoc.first;
                int loc_of_u = indexL[u];
                int u_prime = L[loc_pointer];
                //switch u and u_prime and update location information
                L[loc_pointer] = u;
                L[loc_of_u]=u_prime;
                indexL[u] = loc_pointer;
                indexL[u_prime]=loc_of_u;
                loc_pointer--;
            }
            int tmpiA= loc_pointer+1;
            //if(tmpiA!=iA){
            //    cout<<"back track wrong!"<<endl;
            //}
            // adjust adjlist back

            //TODO, initially, iQ is not well managed, need to refine!!!!! done !
            //move v to Q,
            //int loc_of_v = indexR[v.first];
            //if(loc_of_v!=iP){ // in most cases, this condition cannot be satisfied
            //    int v_prime = R[iP];
            //    R[iP] = v.first;
            //    R[loc_of_v]= v_prime;
            //    indexR[v.first]=iP;
            //    indexR[v_prime]=loc_of_v;
            //}
            //iP=iP+1;


            //1:recover adjlist

            //2:recover R

            //3:recover L

            //normal biclique operation

            //remove v from P to Q

        }

        //vector<int> PAdjLoc(R+iP,R+iB-1);


        //for(int i=)

        //cout<<"test";
        //vector<int> P()

    //recursive_tracker--;
    //cout<<recursive_tracker<<endl;
}

void ooMBE::advPIMBEA_local(int iA,
                            int iQ,
                            int iP,
                            int iB,
                            vector<pair<int,int>> &A,
                            vector<pair<int,int>> &P,
                            vector<pair<int,int>> &Q) {


    for(int i=iQ;i<iB;i++){
        int v=R[i];
        g_bfs[v].v=-1;
        g_bfs[v].count=0;
        g_bfs[v].isV=false;
        isPovit[v]=true; //assuming all are pivots
    }

    for(int i=iQ;i<iB;i++){
        int v = R[i];
        //fill(counts, counts + size_of_R, 0);//this is O(n) but practically faster
        if(isPovit[v]==true){
            int adj_start=adjIndexR[v];
            //vector<int> visited_2_hop;
            //int two_hop_ini = size_of_R-1;

            for(int j=adj_start;j<g.adjR[v].size();j++){
                int u=g.adjR[v][j];
                int ajd_u_start = adjIndexL[u];
                for(int k=ajd_u_start;k<g.adjL[u].size();k++){
                    int v_prime = g.adjL[u][k];
                    if(v!=v_prime){
                        /**
                        if(indexR[v_prime]<iQ||indexR[v_prime]>=iB){
                            cout<<"wrong!"<<endl;
                            cout<<"loc if v_prime: "<<indexR[v_prime]<<endl;
                            cout<<"iQ: "<<iQ<<endl;
                            cout<<"iB: "<<iB<<endl;
                        }
                        **/
                        if(g_bfs[v_prime].v!=v){
                            g_bfs[v_prime].v =v;
                            g_bfs[v_prime].count=0;
                            //visited_2_hop.emplace_back(v_prime);
                            //two_hop[two_hop_ini] = v_prime;
                            //two_hop_ini--;

                        }
                        if(g_bfs[v_prime].v==v){
                            g_bfs[v_prime].count=g_bfs[v_prime].count+1;
                        }
                        int loc_deg_v_prime = g.adjR[v_prime].size()-adjIndexR[v_prime];
                        if(g_bfs[v_prime].count==loc_deg_v_prime){
                            //v_prime is not a pivot
                            isPovit[v_prime]=false;
                            //cout<<"pruned!"<<endl;
                        }
                    }
                }
            }
            /**
            two_hop_ini=two_hop_ini+1;
            //for(auto &v_prime:visited_2_hop){
            for(int j=two_hop_ini;j<=size_of_R;j++){
                int v_prime = two_hop[j];
                int loc_deg_v_prime = g.adjR[v_prime].size()-adjIndexR[v_prime];
                if(g_bfs[v_prime].count>g.adjR[v_prime].size()){
                    cout<<"2-hop based pruning very wrong!"<<endl;
                }
                if(g_bfs[v_prime].count==loc_deg_v_prime){
                    //v_prime is not a pivot
                    isPovit[v_prime]=false;
                    //cout<<"pruned!"<<endl;
                }
                g_bfs[v_prime].count =0;
                g_bfs[v_prime].v=-1;
            }
             **/
        }

    }

    //recursive_tracker++;
    //int pivots_ini = size_of_R-1;
    vector<pair<int,int>> pivots;
    for(int i=0;i<P.size();i++){

        auto v_pair = P[i];
        int v = v_pair.first;
        if(isPovit[v]==true&&indexR[v]>=iP){ //the second condition is very important
            //Pivot p;
            //p.pair=v_pair;
            //p.loc = i;
            pivots.emplace_back(P[i]);
            //pivots_arry[pivots_ini]=P[i];
            //pivots_ini--;
        }
    }
    //pivots_ini=pivots_ini+1;
    //int pivot_starting = iniPos+1;

    //now, pivots should be aligned to the last few positions of P

    //TODO: neighbourhood adjustment has serious problems in this loop.
    for(int z=0;z<pivots.size();z++){
        pair<int,int> v= pivots[z];
        //int global_v = R[iP];
        //move v from the current location to iB+1 if necessary
        vector<pair<int,int>> A_prime;
        vector<pair<int,int>> Q_prime;
        vector<pair<int,int>> P_prime;
        A_prime.reserve(A.size());
        P_prime.reserve(P.size());
        Q_prime.reserve(Q.size());
        int iQ_prime = iQ;


        int locV=indexR[v.first];
        /**
         * This is essentially moving a vertex to a known location
         * I
         * TODO:up to here dec
         */
        if(locV != iB-1){
            int vprime = R[iB-1];

            R[locV] = vprime;
            R[iB-1] = v.first;


            indexR[v.first] = iB-1;
            indexR[vprime] = locV;
            //cout<<"test"<<endl;

        }

        int iB_prime = iB-1;

        //******refine A: works on set L and generate A' <- A\cap adj[v] ******
        //Since H is maintained, in H adj[v] must be a subset of A.
        //As such, adj[v]  \cap A is equivalent to align adj[v] to last few of positions of L

        int localAdjv = adjIndexR[v.first];
        int startLoc = g.adjL.size()-1;
        int local_degree=g.adjR[v.first].size()-localAdjv;
        /**
        if(local_degree>(size_of_L-iA)){
            cout<<"local A wrong"<<endl;
        }
         **/
        for(int i=localAdjv;i<g.adjR[v.first].size();i++){
            int u = g.adjR[v.first][i];
            int locU = indexL[u];

            if(indexL[u]!=startLoc){
                //move vertex to location of startLoc
                /**
                if(startLoc>size_of_L){
                    cout<<"Problem!"<<endl;
                }
                 **/
                int uprime = L[startLoc];

                L[locU] = uprime;
                L[startLoc] =u;

                indexL[u]= startLoc;
                /**
                if(uprime>size_of_L){
                    cout<<" uprime  Problem!"<<endl;
                }
                 **/
                indexL[uprime] = locU;
            }
            startLoc=startLoc-1;
        }


        int iA_prime = startLoc+1;   //cout<<"haha test"<<endl; it is okay
        //reset counts
        for(int i=iQ_prime;i<iB_prime;i++){
            int v = R[i];
            counts[v] = 0;
        }
        //vector<int> v_2_hop;
        int two_hop_ini = size_of_R-1;
        for (int i = iA_prime; i < size_of_L; i++) {
            int u = L[i];
            int start_adj = adjIndexL[u];
            for(int j=start_adj;j<g.adjL[u].size();j++){
                int v_prime = g.adjL[u][j];
                if(indexR[v_prime]>=iQ_prime && indexR[v_prime]<iB_prime){
                    if (counts[v_prime] == 0) {
                        //v_2_hop.emplace_back(v_prime);
                        two_hop[two_hop_ini] = v_prime;
                        two_hop_ini--;
                    }
                    counts[v_prime]++;
                }
            }
        }
        two_hop_ini=two_hop_ini+1;
        //vector<int> not_in_P;
        vector<int> tmp_P;
        //vector<int> not_in_Q;
        vector<int> tmp_Q;
        //for (auto &v_prime: v_2_hop) {
        for(int j=two_hop_ini;j<size_of_R;j++){
            int v_prime = two_hop[j];
            //if(counts[v_prime]==g.adjR[v_prime].size()){ //this degree should be local degree
            if (counts[v_prime] == size_of_L - iA_prime) {//aggressive expansion
                /**
                if(indexR[v_prime]<iP){
                    cout<<"Wrong!! It is impossible and should be pruned by batch-pivots methods!"<<endl;
                }
                 **/
                //move v_prime to B
                iB_prime = iB_prime - 1;
                if (indexR[v_prime] != iB_prime) {
                    //switch vertices at locations of iB_prime and indexR[v_prime]
                    int tmpv = R[iB_prime];
                    R[iB_prime] = v_prime;
                    R[indexR[v_prime]] = tmpv;
                    indexR[tmpv] = indexR[v_prime];
                    indexR[v_prime] = iB_prime;
                }
            }/*else if(counts[v_prime]==0){
                //the second condition is for securing iP==iB_prime, say, P' is an empty set
                if(indexR[v_prime]>=iP&&indexR[v_prime]<iB_prime){
                    not_in_P.emplace_back(v_prime);
                    cout<<"Will this happen>"<<endl;
                }else if(indexR[v_prime]>=iQ_prime && indexR[v_prime]<iP){
                    not_in_Q.emplace_back(v_prime);
                    cout<<"Will this happen>"<<endl;
                }
            }*/else if(counts[v_prime]>0 && (counts[v_prime] != (size_of_L - iA_prime))){
                if(indexR[v_prime]>=iP&&indexR[v_prime]<iB_prime){
                    tmp_P.emplace_back(v_prime);
                }else if(indexR[v_prime]>=iQ_prime && indexR[v_prime]<iP){
                    tmp_Q.emplace_back(v_prime);
                }
            }
        }

        nomb++;
        level_accumulation_ += cur_level_;
        max_level_ = std::max(cur_level_, max_level_);

        int iP_prime;
        //now align P,Q,not_in_P,not_in_Q from iB_prime-1
        int iniPos = iB_prime - 1;
        for (auto &v_prime: tmp_P) {//align them to iB_prime
            if (indexR[v_prime] != iniPos) {
                int tmpv = R[iniPos];
                R[iniPos] = v_prime;
                R[indexR[v_prime]] = tmpv;
                indexR[tmpv] = indexR[v_prime];
                indexR[v_prime] = iniPos;
            }

            //vector<int> localNeighbour;
            //vector<int> other;
            int ini_local = max_degree-1;
            int ini_other = max_degree-1;

            int start_adj = adjIndexR[v_prime];
            for (int j = start_adj; j < g.adjR[v_prime].size(); j++) {  // adjust the neighbours of v_prime
                int u = g.adjR[v_prime][j];
                if (indexL[u] >= iA_prime) {
                    ln_for_vex[ini_local] =u;
                    ini_local--;
                    //localNeighbour.emplace_back(u);
                }else{
                    //other.emplace_back(u);
                    other_for_vex[ini_other]=u;
                    ini_other--;
                }
            }

            ini_local=ini_local+1;
            ini_other=ini_other+1;

            int adj_ini = g.adjR[v_prime].size() - 1;
            for (int j=ini_local;j<max_degree;j++){
                int u = ln_for_vex[j];
                g.adjR[v_prime][adj_ini] = u;
                adj_ini--;
            }
            adjIndexR[v_prime] = adj_ini + 1;
            for (int j=ini_other;j<max_degree;j++){
                int u = other_for_vex[j];
                g.adjR[v_prime][adj_ini] = u;
                adj_ini--;
            }
            P_prime.emplace_back(make_pair(v_prime, adjIndexR[v_prime]));
            iniPos--;
        }

        iP_prime = iniPos+1;
        for (auto &v_prime: tmp_Q) {//align them to iB_prime
            if (indexR[v_prime] != iniPos) {
                int tmpv = R[iniPos];
                R[iniPos] = v_prime;
                R[indexR[v_prime]] = tmpv;
                indexR[tmpv] = indexR[v_prime];
                indexR[v_prime] = iniPos;
            }

            int ini_local = max_degree-1;
            int ini_other = max_degree-1;
            int local_adj = adjIndexR[v_prime];
            for (int j = local_adj; j < g.adjR[v_prime].size(); j++) {  // adjust the neighbours of v_prime
                int u = g.adjR[v_prime][j];
                if (indexL[u] >= iA_prime) {

                    ln_for_vex[ini_local] =u;
                    ini_local--;
                }else{
                    other_for_vex[ini_other]=u;
                    ini_other--;
                }
            }

            ini_local=ini_local+1;
            ini_other=ini_other+1;

            int adj_ini = g.adjR[v_prime].size() - 1;
            for (int j=ini_local;j<max_degree;j++){
                int u = ln_for_vex[j];
                g.adjR[v_prime][adj_ini] = u;
                adj_ini--;
            }
            adjIndexR[v_prime] = adj_ini + 1;
            for (int j=ini_other;j<max_degree;j++){
                int u = other_for_vex[j];
                g.adjR[v_prime][adj_ini] = u;
                adj_ini--;
            }
            Q_prime.emplace_back(make_pair(v_prime, adjIndexR[v_prime]));
            iniPos--;
        }
        iQ_prime = iniPos+1;
        /**
        for (auto &v_prime: not_in_Q) {//align them to iB_prime
            if (indexR[v_prime] != iniPos) {
                int tmpv = R[iniPos];
                R[iniPos] = v_prime;
                R[indexR[v_prime]] = tmpv;
                indexR[tmpv] = indexR[v_prime];
                indexR[v_prime] = iniPos;
            }
            //adjIndexR[v_prime] = g.adjR[v_prime].size();
            iniPos--;
        }
        for (auto &v_prime: not_in_P) {//align them to iB_prime
            if (indexR[v_prime] != iniPos) {
                int tmpv = R[iniPos];
                R[iniPos] = v_prime;
                R[indexR[v_prime]] = tmpv;
                indexR[tmpv] = indexR[v_prime];
                indexR[v_prime] = iniPos;
            }
            //adjIndexR[v_prime] = g.adjR[v_prime].size();
            //P_prime.emplace_back(make_pair(v_prime, adjIndexR[v_prime]));
            iniPos--;
        }
         **/
        //Align old Q into the correct position

        for (int i = iA_prime; i < size_of_L; i++) {
            int u = L[i];


            int ini_local = max_degree-1;
            int ini_other = max_degree-1;

            int local_adj = adjIndexL[u];
            for (int j = local_adj; j < g.adjL[u].size(); j++) {
                int v_prime = g.adjL[u][j];
                if (indexR[v_prime] >= iQ_prime && indexR[v_prime] < iB_prime) {
                    //local_neighbours.emplace_back(v_prime);
                    ln_for_vex[ini_local] = v_prime;
                    ini_local--;
                }else{
                    //other.emplace_back(v_prime);
                    other_for_vex[ini_other] = v_prime;
                    ini_other--;
                }
            }

            ini_local=ini_local+1;
            ini_other=ini_other+1;

            int adj_ini = g.adjL[u].size() - 1;
            for(int j=ini_local;j<max_degree;j++){
                int v= ln_for_vex[j];
                g.adjL[u][adj_ini]=v;
                adj_ini--;
            }
            adjIndexL[u] = adj_ini + 1;
            for(int j=ini_other;j<max_degree;j++){
                int v =other_for_vex[j];
                g.adjL[u][adj_ini]=v;
                adj_ini--;
            }
            A_prime.emplace_back(u, adjIndexL[u]);
        }

        // next is for the recursive call
//dvIMBEA(int iA,
//                     int iQ,
//                     int iP,
//                     int iB,
//                     vector<pair<int,int>> &A,
//                     vector<pair<int,int>> &P,
//                     vector<pair<int,int>> &Q)
        //if(P_prime.size()>0){
        /***
        cout<<"~~~~~~~~~~~"<<endl;
        for(int i=iB_prime;i<size_of_R;i++){
            cout<<R[i]<<",";
        }
        cout<<endl;
        //print L
        for(int i=iA_prime;i<size_of_L;i++){
            cout<<L[i]<<",";
        }
        cout<<endl;
        cout<<"~~~~~~~~~~~"<<endl;
        ***/
        // if(P_prime.size()>2501){
        //     cout<<"order does not work! "<<P_prime.size()<<endl;
        // }
        //if(Q_prime.size()==2){
        //    if(Q_prime[0].first==Q_prime[1].first){
        //        cout<<"Duplicated"<<endl;
        //    }
        //}
        A_prime.resize(size_of_L-iA_prime);
        P_prime.resize(iB_prime-iP_prime);
        Q_prime.resize(iP_prime-iQ_prime);
        if(P_prime.size()>0){
            //advIMBEA(iA_prime, iQ_prime, iP_prime, iB_prime, A_prime, P_prime, Q_prime);
            cur_level_++;
            advPIMBEA_local(iA_prime, iQ_prime, iP_prime, iB_prime, A_prime,
                            P_prime, Q_prime);
            cur_level_--;
        }
        //A strong local O(E) implementation for
        //*** adjust Q
        //first adjust Q  check maximality



        //back tracking adjustment
        //use the stored  iA, iQ, iP, iB, A,P,Q to restore the logical partitions
        //align P to iP to iB-1 and adjust the location index
        //remove v from P
        //Add v to Q


        //for(auto &v_adjLoc:P){
        int v_before_rec = pivots[z].first;
        //if(v_before_rec!=v.first){
        //    cout<<"error mismatch!"<<endl;
        //}
        adjIndexR[v.first]=pivots[z].second;
        /*
         * TODO:int i=z+1;i<P.size();i++ was wrong
         * vertices from z+1 to P.size()-1 are pivots only, all vertices in P except for  v.first
         * should be aligned back.
         */
        int loc_pointer = iB-1;
        int loc_of_to_be_remove=-1;
        for (int i=0;i<P.size();i++){//equivalent to moving v out/
            auto v_adjLoc = P[i];
            int v=v_adjLoc.first;
            int adj_index_value = v_adjLoc.second;
            adjIndexR[v]=adj_index_value; //above okay
            if(v==v_before_rec){
                loc_of_to_be_remove=i;
            }
            if(v!=v_before_rec){
                //if(loc_pointer>=size_of_R){
                //    cout<<"error!+++++ 1018"<<endl;
                //}
                //if(loc_pointer<0){
                //    cout<<"error!----- 1018"<<endl;
                //}
                int v_prime = R[loc_pointer]; //the vertex at the pointer
                int loc_of_v = indexR[v];
                //if(loc_of_v<iP){ // loc_of_v may at iB, this vertex will be dealt with later
                //TODO:6_Feb
                //switch v and v_prime
                R[loc_of_v] = v_prime;
                R[loc_pointer] = v;
                //adjust locations for v and v_prime
                indexR[v] = loc_pointer;
                indexR[v_prime] = loc_of_v;
                loc_pointer--;
            }
        }

        if(loc_of_to_be_remove!=P.size()-1){
            auto tmp = P[P.size()-1];
            P[P.size()-1]=v;
            P[loc_of_to_be_remove] = tmp;
        }

        P.pop_back();
        //iP=iP+1;
        iP=iP+1;
        //TODO:double check here!

        Q.emplace_back(v);
        loc_pointer = iP-1;
        for(auto &v_adjLoc:Q){
            int v = v_adjLoc.first;
            int adj_index_value = v_adjLoc.second;
            adjIndexR[v] = adj_index_value;
            //if(loc_pointer>=size_of_R){
            //    cout<<"error++++! 1046"<<endl;
            //}
            //if(loc_pointer<0){
            //    cout<<"error----! 1046"<<endl;
            //}
            int v_prime = R[loc_pointer];
            int loc_of_v = indexR[v];
            R[loc_of_v] = v_prime;
            R[loc_pointer] = v;
            //adjust locations for v and v_prime
            indexR[v] = loc_pointer;
            indexR[v_prime] = loc_of_v;
            loc_pointer--;
        }
        iQ = loc_pointer+1;
        //align A back to iA to the end of L
        loc_pointer=size_of_L-1;
        for(auto &u_adjLoc:A){
            adjIndexL[u_adjLoc.first] = u_adjLoc.second;

            int u=u_adjLoc.first;
            int loc_of_u = indexL[u];
            int u_prime = L[loc_pointer];
            //switch u and u_prime and update location information
            L[loc_pointer] = u;
            L[loc_of_u]=u_prime;
            indexL[u] = loc_pointer;
            indexL[u_prime]=loc_of_u;
            loc_pointer--;
        }
        int tmpiA= loc_pointer+1;
        //if(tmpiA!=iA){
        //    cout<<"back track wrong!"<<endl;
        //}
        // adjust adjlist back

        //TODO, initially, iQ is not well managed, need to refine!!!!! done !
        //move v to Q,
        //int loc_of_v = indexR[v.first];
        //if(loc_of_v!=iP){ // in most cases, this condition cannot be satisfied
        //    int v_prime = R[iP];
        //    R[iP] = v.first;
        //    R[loc_of_v]= v_prime;
        //    indexR[v.first]=iP;
        //    indexR[v_prime]=loc_of_v;
        //}
        //iP=iP+1;


        //1:recover adjlist

        //2:recover R

        //3:recover L

        //normal biclique operation

        //remove v from P to Q

    }

    //vector<int> PAdjLoc(R+iP,R+iB-1);


    //for(int i=)

    //cout<<"test";
    //vector<int> P()

    //recursive_tracker--;
    //cout<<recursive_tracker<<endl;
}

void ooMBE::adv_mbeStart_reuse() {
    // use an equation to determin which side to enumerate

    cout<<"call adv_mbeStart()"<<endl;
    if (g.max_degree_L <= g.max_degree_R) {
        // do Right
        //    g.switchLR();
        //    cout<<"switched!"<<endl;
    } else {
        //do left, but treat left as right
        g.switchLR();
        cout << "switched!" << endl;
    }


    //g.switchLR();

    //printVector(bcore.uniOrder);


    R = new int[g.adjR.size()];
    size_of_R = g.adjR.size();
    indexR = new int[g.adjR.size()];
    adjIndexR = new int[g.adjR.size()];

    L = new int[g.adjL.size()];
    size_of_L = g.adjL.size();
    indexL = new int[g.adjL.size()];
    adjIndexL = new int[g.adjL.size()];


    two_hop = new int [size_of_R];

    max_degree =g.max_degree_R;
    if(max_degree<g.max_degree_L){
        max_degree=g.max_degree_L;
    }
    //int *ln_for_vex;
    //int *other_for_vex;
    ln_for_vex = new int[max_degree];
    other_for_vex = new int[max_degree];


    g_bfs.resize(size_of_R);
    isPovit.resize(size_of_R);

    fill(isPovit.begin(),isPovit.end(),false);

    counts = new int[g.adjR.size()];


    mbeDSini(g.adjL, g.adjR);
    vector<BFST> bfs(size_of_R);
    //generate the unilateral order
    //TODO:this could be a potential problem
    auto start = high_resolution_clock::now();
    BiCore bcore(g);
    bcore.e2hopCore(g.adjR, g.adjL);
    //end of problem

    if (maxPruning) {
        //Future work, further refine the unilateral order to maximize pruning
        //This optimisation is for sequential algorithm,
        vector<pair<int, int>> auD;
        int initUniCore = bcore.uni_l_core[bcore.uniOrder[1]];
        int count = 0;
        for (auto &v: bcore.uniOrder) {
            if (bcore.uni_l_core[v] == initUniCore) {
                count++;
            } else {
                auD.emplace_back(make_pair(initUniCore, count));
                initUniCore = bcore.uni_l_core[v];
                count = 1;
            }
        }

        auD.emplace_back(make_pair(initUniCore, count));

        int ini = 0;
        for (auto &pair: auD) { // adjust pairs for initial position of very unilateral core.
            int tmp = pair.second;
            pair.second = ini;
            ini = ini + tmp;
        }

        vector<vector<int>> vertex_with_same_uni_cores;
        for (int i = 0; i < auD.size(); i++) {
            pair<int, int> cp = auD[i];
            vector<int> ertex_with_same_uni_core;
            if (i < auD.size() - 1) {
                pair<int, int> np = auD[i + 1];
                for (int j = cp.second; j < np.second; j++) {
                    ertex_with_same_uni_core.emplace_back(bcore.uniOrder[j]);
                }
                vertex_with_same_uni_cores.emplace_back(ertex_with_same_uni_core);
            }
            if (i == auD.size() - 1) {
                for (int j = cp.second; j < bcore.uniOrder.size(); j++) {
                    ertex_with_same_uni_core.emplace_back(bcore.uniOrder[j]);
                }
                vertex_with_same_uni_cores.emplace_back(ertex_with_same_uni_core);
            }
        }
        bcore.uniOrder.clear();
        for (auto &vec: vertex_with_same_uni_cores) {
            std::sort(vec.begin(), vec.end(),
                      bind(&ooMBE::orderSort, this, std::placeholders::_1, std::placeholders::_2));
            bcore.uniOrder.insert(bcore.uniOrder.end(), vec.begin(), vec.end());
        }


        //cout<<"test tag"<<endl;//the above code has been tested
    }

    //printVector(bcore.uniOrder);
    //test
    //if (bcore.uniOrder.size() != size_of_R) {
    //    cout << "Wrong!" << endl;
    //}
    //adjust R, IndexR, but adjIndexR keeps the same
    vector<int> orderIndex(bcore.uniOrder.size());
    for (int i = 0; i < bcore.uniOrder.size(); i++) {
        int v = bcore.uniOrder[i];
        R[i] = v;
        indexR[v] = i;
        orderIndex[v] = i;
    }


    // let refine it now

    //compute order preserved batch pivots
    //refresh that give vertex v in R, indexR[v] returns the location of v in R, where R is an array. The same for L
    //parameters adjL, adjR, the order, adjIndexL, adjIndexR
    //note that the order is for R
    //vector<int> pivots;
    vector<PivotsNeighbour> pivots_plus_neighbours;
    //orderBatchPivots(pivots);
    orderBatchPivots_neighbours(pivots_plus_neighbours);


    int iP = 0;//g.adjR.size();
    int iQ = iP;
    int iB = g.adjR.size();


    int iA = 0;

    vector<pair<int, int>> RP;
    for (int i = iP; i < g.adjR.size(); i++) {
        RP.emplace_back(make_pair(i, 0));
    }


    vector<pair<int, int>> RQ, RB;
    RQ.reserve(pivots_plus_neighbours.size());
    vector<pair<int, int>> LA;
    for (int i = 0; i < g.adjL.size(); i++) {
        LA.emplace_back(make_pair(i, 0));
    }


    for (int i = 0; i < size_of_L; i++) {
        adjIndexL[i] = 0;
    }

    for (int i = 0; i < size_of_R; i++) {
        adjIndexR[i] = 0;
    }


    //cout<<"still fine!"<<endl;
    //for(int i=iP;i<size_of_R;i++){// When backtracking, R should be the same as here
    // This can be achieved by using the stored RP
    //for(auto &pa: RP){
    //for(int i=0;i<RP.size();i++){
    //cout<<"max 2-hop core"<<bcore.maxE2Core<<endl;
    //max2Core = bcore.maxE2Core;
    //pivots_arry = new pair<int,int> [size_of_R];
    cout<<"breadth:"<<pivots_plus_neighbours.size()<<endl;
    int count =0;


    int size_of_P;
    //for(int z=0;z<pivots_plus_neighbours.size();z++){

    for (auto &vn: pivots_plus_neighbours) {
        //cout<<"Progress:"<<++count<<endl;
        //move v to iB
        int v = vn.vertex;
        /**
        if(count==10000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==20000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==30000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==40000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==50000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==60000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==70000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==80000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==90000){
            cout<<count<<":"<<nomb<<endl;
        }else **/
        if(count==1000000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==1500000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==2000000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==2500000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==3100000){
            cout<<count<<":"<<nomb<<endl;
        }

        count++;
        //cout<<"Processing "<<v<<endl;
        int iP_prime = 0, iA_prime = 0, iB_prime = 0, iQ_prime = 0;
        iB_prime = iB - 1;
        if (indexR[v] != iB - 1) {
            int tmpv = R[iB - 1];
            R[iB - 1] = v;
            R[indexR[v]] = tmpv;
            indexR[tmpv] = indexR[v];
            indexR[v] = iB - 1;
        }

        vector<pair<int, int>> A_prime;
        vector<pair<int, int>> Q_prime;
        vector<pair<int, int>> P_prime;

        //move neighbours of v in L
        //   compute starting position of the neighbour
        int nV = g.adjR[v].size();
        iA_prime = size_of_L - nV;
        int cPos = size_of_L - 1;
        for (auto &u: g.adjR[v]) {
            if (indexL[u] != cPos) {
                int tmpu = L[cPos];
                L[cPos] = u;
                L[indexL[u]] = tmpu;
                //int tmpPos = indexL[u];
                indexL[tmpu] = indexL[u];
                indexL[u] = cPos;
            }
            cPos--;
        }
        //iQ_prime = iQ;
        //Put vertex in Q into the correct position

        //shrank P
        //let's do it with O(E)
        //*****  ******
        //the loop below runs in O(E)
        //vector<pair<int,int>> counts;

        //TODO problem starts here! Fixed!
        //align 2-hop neighbours starting from iB-1
        //int ini_pos = iB_prime-1;
        vector<int> agg;
        //vector<int> notP;
        for(auto &v_p:vn.adj){
            //vector<int> localNeighbour;
            //vector<int> other;
            int ini_local=max_degree-1;
            int ini_other=max_degree-1;
            for (int j = 0; j < g.adjR[v_p].size(); j++) {  // adjust the neighbours of v_prime
                int u = g.adjR[v_p][j];
                if (indexL[u] >= iA_prime) {
                    ln_for_vex[ini_local] = u;
                    ini_local--;
                    //localNeighbour.emplace_back(u);
                }else{
                    other_for_vex[ini_other] = u;
                    ini_other--;
                }
            }
            ini_local = ini_local+1;
            ini_other = ini_other+1;
            int adj_ini = g.adjR[v_p].size() - 1;
            //for (auto & u:localNeighbour){
            for(int j=ini_local;j<max_degree;j++){
                int u = ln_for_vex[j];
                g.adjR[v_p][adj_ini] = u;
                adj_ini--;
            }
            adjIndexR[v_p] = adj_ini + 1;
            for(int j=ini_other;j<max_degree;j++){
                int u = other_for_vex[j];
                g.adjR[v_p][adj_ini] = u;
                adj_ini--;
            }
            if((max_degree-ini_local)==size_of_L - iA_prime){
                agg.emplace_back(v_p);
            }else if((max_degree-ini_local)>0&&(max_degree-ini_local)!=(size_of_L - iA_prime)){
                P_prime.emplace_back(make_pair(v_p, adjIndexR[v_p]));
            }

        }
        //iP_prime = ini_pos+1;

        int ini_pos = iB_prime-1;
        for(auto &v_b:agg){
            if(indexR[v_b]!=ini_pos){
                int tmpv = R[ini_pos];
                R[ini_pos] = v_b;
                R[indexR[v_b]] = tmpv;
                indexR[tmpv] = indexR[v_b];
                indexR[v_b] = ini_pos;
            }
            ini_pos--;
        }
        //iP_prime = ini_pos+1;
        iB_prime = ini_pos+1;
        ini_pos = iB_prime-1;

        for(auto &v_pp:P_prime){
            int v_b = v_pp.first;
            if(indexR[v_b]!=ini_pos){
                int tmpv = R[ini_pos];
                R[ini_pos] = v_b;
                R[indexR[v_b]] = tmpv;
                indexR[tmpv] = indexR[v_b];
                indexR[v_b] = ini_pos;
            }
            ini_pos--;
        }
        iP_prime = ini_pos+1;
        /*cout<<"***********"<<endl;
        for(int i=iB_prime;i<size_of_R;i++){
            cout<<R[i]<<",";
        }
        cout<<endl;
        //print L
        for(int i=iA_prime;i<size_of_L;i++){
            cout<<L[i]<<",";
        }
        cout<<endl;
        cout<<"***********"<<endl;
        */
        nomb++;
        level_accumulation_ += cur_level_;
        max_level_ = std::max(cur_level_, max_level_);


        //Align old Q into the correct position
        int iniPos_Q = iP_prime - 1;
        for (auto &vq:vn.Q){
            //int vq = vqp.first;
            if (indexR[vq] != iniPos_Q) {
                int tmpv = R[iniPos_Q];
                R[iniPos_Q] = vq;
                R[indexR[vq]] = tmpv;
                indexR[tmpv] = indexR[vq];
                indexR[vq] = iniPos_Q;
            }

            //adjust adj list for vq;
            iniPos_Q--;
        }

        iQ_prime = iniPos_Q + 1;


        //now shrank Q to Q_prime
       // vector<int> notQ;
        //for (int i = iQ_prime; i < iB_prime; i++) {
        for (int i = iQ_prime; i < iP_prime; i++) {
            int vq = R[i];
            int counts = 0;
            int ini_local=max_degree-1;
            int ini_other=max_degree-1;

            //vector<int> localNeighbours;
            //vector<int> other;
            for (int j = 0; j < g.adjR[vq].size(); j++) {
                int u = g.adjR[vq][j];
                if (indexL[u] >= iA_prime) {
                    counts++;
                    //localNeighbours.emplace_back(u);
                    ln_for_vex[ini_local] = u;
                    ini_local--;
                }else{
                    //other.emplace_back(u);
                    other_for_vex[ini_other] = u;
                    ini_other--;
                }
            }

            ini_local=ini_local+1;
            ini_other=ini_other+1;
            //vq.adj \cap A \ne emptyset
            int adj_ini = g.adjR[vq].size() - 1;
            for(int j = ini_local;j<max_degree;j++){
                int u = ln_for_vex[j];
                g.adjR[vq][adj_ini]=u;
                adj_ini--;
            }
            adjIndexR[vq] = adj_ini+1;
            for(int j= ini_other;j<max_degree;j++){
                int u = other_for_vex[j];
                g.adjR[vq][adj_ini]=u;
                adj_ini--;
            }

            if(counts>0){
                Q_prime.emplace_back(make_pair(vq,adjIndexR[vq]));
            }
            /**
            if (counts == 0) {
                //notQ.emplace_back(vq);
                cout<<"should not happlen!"<<endl;

            }
            **/
        }

        /**
        for (auto &vq: notQ) {
            int tmpv = R[iQ_prime];
            R[iQ_prime] = vq;
            R[indexR[vq]] = tmpv;
            indexR[tmpv] = indexR[vq];
            indexR[vq] = iQ_prime;
            iQ_prime = iQ_prime + 1;
        }
        **/

        //update adj index based on the current iQ, iP and iB
        for (int i = iA_prime; i < size_of_L; i++) {
            int u = L[i];

            //vector<int> local_neighbours;
            //vector<int> other;
            int ini_local=max_degree-1;
            int ini_other=max_degree-1;

            for (int j = 0; j < g.adjL[u].size(); j++) {
                int v_prime = g.adjL[u][j];
                if (indexR[v_prime] >= iQ_prime && indexR[v_prime] < iB_prime) {
                    ln_for_vex[ini_local]=v_prime;
                    ini_local--;
                }else{
                    other_for_vex[ini_other]=v_prime;
                    ini_other--;
                }
            }
            ini_local=ini_local+1;
            ini_other=ini_other+1;
            int adj_ini = g.adjL[u].size() - 1;
            for(int j=ini_local;j<max_degree;j++){
                int v = ln_for_vex[j];
                g.adjL[u][adj_ini]=v;
                adj_ini--;
            }
            adjIndexL[u] = adj_ini + 1;
            for(int j=ini_other;j<max_degree;j++){
                int v = other_for_vex[j];
                g.adjL[u][adj_ini]=v;
                adj_ini--;
            }
            A_prime.emplace_back(u, adjIndexL[u]);
        }

        if(L[size_of_L-1]==9&&L[size_of_L-2]==6){
            cout<<endl;
        }

        if(P_prime.size()>0){
            //advIMBEA(iA_prime, iQ_prime, iP_prime, iB_prime, A_prime, P_prime, Q_prime);
            if(iB_prime>=size_of_R){
                cout<<"why!!!"<<endl;
            }
            cur_level_++;
            advPIMBEA_local(iA_prime, iQ_prime, iP_prime, iB_prime, A_prime,
                            P_prime, Q_prime);
            cur_level_--;
        }

            int iB = iB - 1;
            //size_of_P = size_of_P-1;
            //iP = iB-;
            int last_pos_of_Q=count-1;
            //move v to position count
            /**
            if(indexR[v]!=last_pos_of_Q){
                int tmpv = R[last_pos_of_Q];
                R[last_pos_of_Q] = v;
                R[indexR[v]]=tmpv;
                indexR[tmpv] = indexR[v];
                indexR[v] = last_pos_of_Q;
            }
             **/
            //iQ=0;
            //RQ.emplace_back(make_pair(v,0));
            iP=last_pos_of_Q+1;
            //cout<<iP<<endl;
        /**
        }else{
            int ini_pos = iB - 1;
            //for(auto &pair:RP){
            int loc_of_v_in_RP;
            for (int i = 0; i < RP.size(); i++) {
                pair<int, int> p = RP[i];
                int bv = p.first;
                if (bv == v) {
                    loc_of_v_in_RP = i;
                }
                adjIndexR[bv] = p.second;
                if (indexR[bv] != ini_pos) {
                    int tmpv = R[ini_pos];
                    R[ini_pos] = bv;
                    R[indexR[bv]] = tmpv;
                    indexR[tmpv] = indexR[bv];
                    indexR[bv] = ini_pos;
                }
                ini_pos--;
            }
            //remove v form RP
            iP = ini_pos + 1;
            if(RP.size()>0){//this is always true
                if(RP.size()>=2){

                    pair<int, int> tmpp = RP[RP.size() - 1];
                    pair<int, int> tmpvp = RP[loc_of_v_in_RP];
                    RP[loc_of_v_in_RP] = tmpp;
                    RP.pop_back();
                    RQ.emplace_back(tmpvp);
                }else if (RP.size()==1) {
                    RQ.emplace_back(RP[0]);
                    RP.pop_back();
                }

                //refine array and indices
                int b_loc_of_v = indexR[v];
                R[b_loc_of_v] = R[iP];
                indexR[R[iP]] = b_loc_of_v;
                R[iP] = v;
                indexR[v] = iP;
                iP = iP + 1;
                int ini_pos_Q = iP - 1;
                for (auto &p: RQ) {
                    if (indexR[p.first] != ini_pos_Q) {
                        int tmpv = R[ini_pos_Q];
                        R[ini_pos_Q] = p.first;
                        R[indexR[p.first]] = tmpv;
                        indexR[tmpv] = indexR[p.first];
                        indexR[p.first] = ini_pos_Q;
                    }
                    ini_pos_Q--;
                    adjIndexR[p.first] = 0;
                }
                iQ = ini_pos_Q + 1;
                //if(iQ!=0){
                //    cout<<"Hehehehe"<<endl;
                // }
                //cout<<iQ<<endl;
                for(int i=iA_prime;i<size_of_L;i++){
                    int u = L[i];
                    adjIndexL[u]=0;
                }
                //fill(adjIndexL+(iA_prime-1), adjIndexL + size_of_L, 0);
            }

            cout<<iP<<endl;
        }
        **/
        //cout << "testint" << endl;
        //back track
        //align P,Q,B back first




        //for()
        //adj align back

        /***
        for(int j=i+1;j<RP.size();j++){
            int v_prime = RP[j].first;
            //int number_of_intersection = 0;
            int adj_v_start = adjIndexR[v_prime];
            int iniPos = g.adjR[v_prime].size()-1;
            for(int z=adj_v_start;z<g.adjR[v_prime].size();z++){
                int u = g.adjR[v_prime][z];
                if(indexL[u]>=iA_prime){
                    //move it to the iniPos
                    if(z<iniPos){
                        int tmpu =  g.adjR[v_prime][iniPos];
                        g.adjR[v_prime][iniPos] = u;
                        g.adjR[v_prime][z] = tmpu;

                    }
                    iniPos--;
                }
            }
            //iniPos++ is the adjIndexR[v_prime]
            int adj_v_prime_update = iniPos+1;
            adjIndexR[v_prime] = adj_v_prime_update;
            int intersection = g.adjR[v_prime].size()-adj_v_prime_update;
            if(intersection==g.adjR[v_prime].size()){
                //aggressive expansion
            }
        }
        ***/
        //cout<<"test"<<endl;

    }

    /**
     * vector<int> bfs(size_of_L,-1);
    for (auto & u:bcore.uniOrder){//loop 1
        cout<<u<<",";
        vector<pair<int,int>> A_prime;
        vector<pair<int,int>> Q_prime;
        vector<pair<int,int>> P_prime;
        //compute induced subgraph
        int loc_of_u=indexL[u];
        vertex_swap(loc_of_u,size_of_L-1,L,indexL);
        int iB_prime=size_of_L;
        //int iP_prime=
        //algin neighbours of u to A
        int ini_loc = size_of_R-1;
        for(auto &v:g.adjL[u]){
            int loc_of_v = indexR[v];
            if(loc_of_v!=ini_loc){
                vertex_swap(loc_of_v,ini_loc,R,indexR);
            }
            ini_loc--;
        }
        int iA_prime=ini_loc+1; //starting position of iA_prime
        //prepare P, i.e., vertices
        ini_loc = iB_prime-1;

        for(int loc=iA_prime;loc<size_of_R;loc++){
            int v = R[loc];
            for(auto &u_prime: g.adjR[v]){
                if(bcore.index_for_uniOrder[u_prime]>bcore.index_for_uniOrder[u]){// appears after u in the unilateral order
                    if(bfs[u_prime]!=u){
                        //first time it visits
                        bfs[u_prime] = u;
                        int loc_u_prime = indexL[u_prime];
                        if(loc_u_prime!=ini_loc){
                            vertex_swap(loc_u_prime,ini_loc,L,indexL);
                        }
                        ini_loc--;
                    }
                }

            }
        }
        int iP_prime = ini_loc+1;

        //adjust adjIndexR and adjust adjIndexL
        for(int loc=iA_prime;loc<size_of_R;loc++){
            int v_prime = R[loc];
            for(auto &u_prime:g.adjR[v_prime]){
                int loc_u_prime = indexL[u_prime];
                if(loc_u_prime>=iP_prime&&loc_u_prime<iB_prime){

                }
            }

        }

        for(int loc=iP_prime;loc<iB_prime;loc++){

        }


        int iB=size_of_L;


        //need

    }
    **/


    //if(R!= 0) {delete [] R;R=0;}
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << duration.count() << endl;
}


void ooMBE::adv_order_PIMBEA_local(int iA,
                                   int iQ,
                                   int iP,
                                   int iB,
                                   vector<pair<int,int>> &A,
                                   vector<pair<int,int>> &P,
                                   vector<pair<int,int>> &Q) {
#ifdef EXIT_TIME
  if(getExcuteTime() >= EXIT_TIME)
  {
    printExcuteTime();
    PrintMemory();
    cout<<"processing_nodes: "<<processing_nodes<<endl;
    cout<<"maximal nodes: "<<nomb<<endl;
    exit(0);
  }
#endif // 
    for(int i=iQ;i<iB;i++){
        int v=R[i];
        g_bfs[v].v=-1;
        g_bfs[v].count=0;
        g_bfs[v].isV=false;
        isPovit[v]=true; //assuming all are pivots
    }

    vector<PivotsNeighbour> pivots;
    pivots.reserve(iB-iP);
    int nop=0;
    for(int i=iQ;i<iB;i++){
        int v = R[i];
        //fill(counts, counts + size_of_R, 0);//this is O(n) but practically faster
        if(isPovit[v]==true){
            processing_nodes++;
            int adj_start=adjIndexR[v];
            int loc_deg_v = g.adjR[v].size()-adj_start;
            PivotsNeighbour p;
            p.vertex=v;
            p.adj.reserve(iB-iP);
            p.Q.reserve(iP-iQ);
            p.agg.reserve(iB-iP);
            int noAdj=0;
            int noQ=0;
            int noAgg=0;
            //vector<int> visited_2_hop;
            int two_hop_ini = size_of_R-1;
            /*for(int j=0;j<adj_start;j++)
            {
              if(std::find(rightA.begin(), rightA.end(), g.adjR[v][j])!=rightA.end())
              {
                std::cout<<"error"<<std::endl;
                exit(0);
              }
            }*/
      auto start = high_resolution_clock::now();
            for(int j=adj_start;j<g.adjR[v].size();j++){
                int u=g.adjR[v][j];
                int ajd_u_start = adjIndexL[u];
                for(int k=ajd_u_start;k<g.adjL[u].size();k++){
                    int v_prime = g.adjL[u][k];

                    if(v!=v_prime&&indexR[v_prime]>indexR[v]){

                        if(g_bfs[v_prime].v!=v){
                            g_bfs[v_prime].v =v;
                            g_bfs[v_prime].count=0;
                            two_hop[two_hop_ini] = v_prime;
                            two_hop_ini--;
                            //p.adj.emplace_back(v_prime);
                            //noAdj++;
                        }
                        if(g_bfs[v_prime].v==v){
                            g_bfs[v_prime].count=g_bfs[v_prime].count+1;
                        }
                        int loc_deg_v_prime = g.adjR[v_prime].size()-adjIndexR[v_prime];
                        if(g_bfs[v_prime].count==loc_deg_v_prime){
                            //v_prime is not a pivot
                            isPovit[v_prime]=false;
                            //cout<<"pruned!"<<endl;
                        }
                        if(v!=v_prime&&g_bfs[v_prime].count==loc_deg_v){//lixu
                            p.agg.emplace_back(v_prime);
                            noAgg++;
                        }
                    }else if(v!=v_prime&&indexR[v_prime]<indexR[v]){
                        if(g_bfs[v_prime].v!=v){
                            g_bfs[v_prime].v =v;
                            if(isPovit[v_prime]==true){
                                p.Q.emplace_back(v_prime);
                                noQ++;
                            }
                        }
                        
                    }
                }
            }

            two_hop_ini = two_hop_ini+1;
            for(int j=two_hop_ini;j<size_of_R;j++){
                int v_2_hop = two_hop[j];
                //int loc_deg_v_prime = g.adjR[v_2_hop].size()-adjIndexR[v_2_hop];
                if(g_bfs[v_2_hop].count<loc_deg_v){
                    //if(g_bfs[v_2_hop].count==0){
                    //    cout<<"Wrong!"<<endl;
                    //}
                    p.adj.emplace_back(v_2_hop);
                    noAdj++;
                }
            }


            p.locAdj = adjIndexR[v];
            p.adj.resize(noAdj);
            p.Q.resize(noQ);
            p.agg.resize(noAgg);
            if(indexR[v]>=iP){

                if(p.adj.size()==0){
                    //cout<<"more pruning opp"<<endl;
#ifdef THRESHOLD
                if(g.adjR.size() - iB + noAgg + 1 < THRESHOLD) //lixu 
#endif
                    nomb++;
                    level_accumulation_ += cur_level_;
                    max_level_ = std::max(cur_level_, max_level_);
                }
                if(p.adj.size()>0){
                    pivots.emplace_back(p);
                    nop++;
                }
            }

          auto end = high_resolution_clock::now();
          collection_time += duration<double>(end - start);
        }

    }
    pivots.resize(nop);

        for(int z=0;z<pivots.size();z++){

            if(z>0){
                for (int i=0;i<P.size();i++){//equivalent to moving v out/
                    auto v_adjLoc = P[i];
                    int v=v_adjLoc.first;
                    int adj_index_value = v_adjLoc.second;
                    adjIndexR[v]=adj_index_value; //above okay
                }

                for(auto &v_adjLoc:Q){
                    int v = v_adjLoc.first;
                    int adj_index_value = v_adjLoc.second;
                    adjIndexR[v] = adj_index_value;
                }

                //align A back to iA to the end of L
                //int loc_pointer=size_of_L-1;
                for(auto &u_adjLoc:A){
                    adjIndexL[u_adjLoc.first] = u_adjLoc.second;
                }
            }

            PivotsNeighbour v = pivots[z];
            //int global_v = R[iP];
            //move v from the current location to iB+1 if necessary
            vector<pair<int,int>> A_prime;
            vector<pair<int,int>> Q_prime;
            vector<pair<int,int>> P_prime;
            A_prime.reserve(A.size());
            P_prime.reserve(P.size());
            Q_prime.reserve(Q.size());
            //int iQ_prime = iQ;


            int locV=indexR[v.vertex];
            /**
             * This is essentially moving a vertex to a known location
             * I
             * TODO:up to here dec
             */
            if(locV != iB-1){
                int vprime = R[iB-1];

                R[locV] = vprime;
                R[iB-1] = v.vertex;


                indexR[v.vertex] = iB-1;
                indexR[vprime] = locV;
                //cout<<"test"<<endl;

            }

            int iB_prime = iB-1;


            int localAdjv = adjIndexR[v.vertex];
            int startLoc = g.adjL.size()-1;
            int local_degree=g.adjR[v.vertex].size()-localAdjv;


            for(int i=localAdjv;i<g.adjR[v.vertex].size();i++){
                int u = g.adjR[v.vertex][i];
                int locU = indexL[u];

                if(indexL[u]!=startLoc){
                    //move vertex to location of startLoc
                    int uprime = L[startLoc];

                    L[locU] = uprime;
                    L[startLoc] =u;

                    indexL[u]= startLoc;

                    indexL[uprime] = locU;
                }
                startLoc=startLoc-1;
            }


            int iA_prime = startLoc+1;   //cout<<"haha test"<<endl; it is okay
            //reset counts

            //vector<int> agg;
            //agg.reserve(P.size());
            //int noagg=0;
            for(auto &v_n:v.adj){
                int l_degree_v_n = g.adjR[v_n].size()-adjIndexR[v_n];
                //if(l_degree_v_n==(size_of_L - iA_prime)){

                    //noagg++;
                //}else{
                    int ini_local=max_degree-1;
                    int ini_other=max_degree-1;

                    int adj_loc = adjIndexR[v_n];
                    for(int j=adj_loc;j<g.adjR[v_n].size();j++){
                        int u= g.adjR[v_n][j];
                        if (indexL[u] >= iA_prime) {
                            ln_for_vex[ini_local] = u;
                            ini_local--;
                            //localNeighbour.emplace_back(u);
                        }else{
                            other_for_vex[ini_other] = u;
                            ini_other--;
                        }
                    }
                    ini_local = ini_local+1;
                    ini_other = ini_other+1;


                    int adj_ini = g.adjR[v_n].size() - 1;
                    //for (auto & u:localNeighbour){
                    for(int j=ini_local;j<max_degree;j++){
                        int u = ln_for_vex[j];
                        g.adjR[v_n][adj_ini] = u;
                        adj_ini--;
                    }
                    adjIndexR[v_n] = adj_ini + 1;
                    for(int j=ini_other;j<max_degree;j++){
                        int u = other_for_vex[j];
                        g.adjR[v_n][adj_ini] = u;
                        adj_ini--;
                    }
                    //if((max_degree-ini_local)==size_of_L - iA_prime){
                    //    agg.emplace_back(v_n);
                    //    noagg++;
                        //cout<<"Must be something here!"<<endl;
                    //}else if((max_degree-ini_local)>0&&(max_degree-ini_local)<(size_of_L - iA_prime)){
                        P_prime.emplace_back(make_pair(v_n, adjIndexR[v_n]));
                    //}
               // }

            }

            int ini_pos = iB_prime-1;

            for(auto &v_b:v.agg){
                if(indexR[v_b]!=ini_pos){
                    int tmpv = R[ini_pos];
                    R[ini_pos] = v_b;
                    R[indexR[v_b]] = tmpv;
                    indexR[tmpv] = indexR[v_b];
                    indexR[v_b] = ini_pos;
                }
                ini_pos--;
            }
            //int iP_prime = ini_pos+1;
            iB_prime = ini_pos+1;
            ini_pos = iB_prime-1;



            for(auto &v_pp:P_prime){
                int v_b = v_pp.first;
                if(indexR[v_b]!=ini_pos){
                    int tmpv = R[ini_pos];
                    R[ini_pos] = v_b;
                    R[indexR[v_b]] = tmpv;
                    indexR[tmpv] = indexR[v_b];
                    indexR[v_b] = ini_pos;
                }
                ini_pos--;
            }
            int iP_prime = ini_pos+1;
#ifdef THRESHOLD
            if(g.adjR.size() - iB_prime >= THRESHOLD)continue;//lixu
#endif

            nomb++;
            level_accumulation_ += cur_level_;
            max_level_ = std::max(cur_level_, max_level_);

            if(P_prime.size()==0){
                continue;
            }




            int iniPos_Q = iP_prime - 1;
            for (auto &vq:v.Q){
                //int vq = vqp.first;
                if (indexR[vq] != iniPos_Q) {
                    int tmpv = R[iniPos_Q];
                    R[iniPos_Q] = vq;
                    R[indexR[vq]] = tmpv;
                    indexR[tmpv] = indexR[vq];
                    indexR[vq] = iniPos_Q;
                }

                //adjust adj list for vq;
                iniPos_Q--;
            }

            int iQ_prime = iniPos_Q + 1;

            for (int i = iQ_prime; i < iP_prime; i++) {
                int vq = R[i];
                int counts = 0;
                int ini_local=max_degree-1;
                int ini_other=max_degree-1;

                //vector<int> localNeighbours;
                //vector<int> other;
                int adj_loc = adjIndexR[vq];
                for (int j = adj_loc; j < g.adjR[vq].size(); j++) {
                    int u = g.adjR[vq][j];
                    if (indexL[u] >= iA_prime) {
                        counts++;
                        //localNeighbours.emplace_back(u);
                        ln_for_vex[ini_local] = u;
                        ini_local--;
                    }else{
                        //other.emplace_back(u);
                        other_for_vex[ini_other] = u;
                        ini_other--;
                    }
                }

                ini_local=ini_local+1;
                ini_other=ini_other+1;
                //vq.adj \cap A \ne emptyset
                int adj_ini = g.adjR[vq].size() - 1;
                for(int j = ini_local;j<max_degree;j++){
                    int u = ln_for_vex[j];
                    g.adjR[vq][adj_ini]=u;
                    adj_ini--;
                }
                adjIndexR[vq] = adj_ini+1;
                for(int j= ini_other;j<max_degree;j++){
                    int u = other_for_vex[j];
                    g.adjR[vq][adj_ini]=u;
                    adj_ini--;
                }

                if(counts>0){
                    Q_prime.emplace_back(make_pair(vq,adjIndexR[vq]));
                }

            }


            //update adj index based on the current iQ, iP and iB
            for (int i = iA_prime; i < size_of_L; i++) {
                int u = L[i];

                //vector<int> local_neighbours;
                //vector<int> other;
                int ini_local=max_degree-1;
                int ini_other=max_degree-1;

                int adj_loc = adjIndexL[u];
                for (int j = adj_loc; j < g.adjL[u].size(); j++) {
                    int v_prime = g.adjL[u][j];
                    if (indexR[v_prime] >= iQ_prime && indexR[v_prime] < iB_prime) {
                        ln_for_vex[ini_local]=v_prime;
                        ini_local--;
                    }else{
                        other_for_vex[ini_other]=v_prime;
                        ini_other--;
                    }
                }
                ini_local=ini_local+1;
                ini_other=ini_other+1;
                int adj_ini = g.adjL[u].size() - 1;
                for(int j=ini_local;j<max_degree;j++){
                    int v = ln_for_vex[j];
                    g.adjL[u][adj_ini]=v;
                    adj_ini--;
                }
                adjIndexL[u] = adj_ini + 1;
                for(int j=ini_other;j<max_degree;j++){
                    int v = other_for_vex[j];
                    g.adjL[u][adj_ini]=v;
                    adj_ini--;
                }
                A_prime.emplace_back(u, adjIndexL[u]);
            }

            A_prime.resize(size_of_L-iA_prime);
            P_prime.resize(iB_prime-iP_prime);
            Q_prime.resize(iP_prime-iQ_prime);

            if(P_prime.size()>0){
                cur_level_++;
                // advIMBEA(iA_prime, iQ_prime, iP_prime, iB_prime, A_prime,
                // P_prime, Q_prime);
                adv_order_PIMBEA_local(iA_prime, iQ_prime, iP_prime, iB_prime, A_prime, P_prime, Q_prime);
                cur_level_--;
            }

            // adjust adjlist back

            //TODO, initially, iQ is not well managed, need to refine!!!!! done !
            //move v to Q,
            //int loc_of_v = indexR[v.first];
            //if(loc_of_v!=iP){ // in most cases, this condition cannot be satisfied
            //    int v_prime = R[iP];
            //    R[iP] = v.first;
            //    R[loc_of_v]= v_prime;
            //    indexR[v.first]=iP;
            //    indexR[v_prime]=loc_of_v;
            //}
            //iP=iP+1;


            //1:recover adjlist

            //2:recover R

            //3:recover L

            //normal biclique operation

            //remove v from P to Q

        }



    //TODO: neighbourhood adjustment has serious problems in this loop.

    //vector<int> PAdjLoc(R+iP,R+iB-1);


    //for(int i=)

    //cout<<"test";
    //vector<int> P()

    //recursive_tracker--;
    //cout<<recursive_tracker<<endl;
}

void ooMBE::printExcuteTime()
{
  auto time_now = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(time_now - time_start);
  exe_time_ = duration.count() / 1000000.0;
//   cout << "Total processing time: " << duration.count() / 1000000.0 << "s"
//        << endl;
}

double ooMBE::getExcuteTime()
{
  auto time_now = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(time_now - time_start);
  return duration.count()/1000000.0;
}

void ooMBE::adv_mbeStart_reuse_full(int order) {
    cur_level_ = 1;
    max_level_ = 0;
    level_accumulation_ = 0;

    // use an equation to determin which side to enumerate
    // cout<<"call adv_mbeStart()"<<endl;
    if (g.max_degree_L <= g.max_degree_R) {//lixu
        // do Right
        //    g.switchLR();
        //    cout<<"switched!"<<endl;
    } else {
        //do left, but treat left as right
        g.switchLR();
        // cout << "switched!" << endl;
    }
    //g.switchLR();


    //printVector(bcore.uniOrder);


    R = new int[g.adjR.size()];
    size_of_R = g.adjR.size();
    indexR = new int[g.adjR.size()];
    adjIndexR = new int[g.adjR.size()];

    L = new int[g.adjL.size()];
    size_of_L = g.adjL.size();
    indexL = new int[g.adjL.size()];
    adjIndexL = new int[g.adjL.size()];

    // std::cout<<"enumerated vertices size: "<<size_of_R<<std::endl;

    two_hop = new int [size_of_R];

    max_degree =g.max_degree_R;
    if(max_degree<g.max_degree_L){
        max_degree=g.max_degree_L;
    }
    //int *ln_for_vex;
    //int *other_for_vex;
    ln_for_vex = new int[max_degree];
    other_for_vex = new int[max_degree];


    g_bfs.resize(size_of_R);
    isPovit.resize(size_of_R);

    fill(isPovit.begin(),isPovit.end(),false);

    counts = new int[g.adjR.size()];


    mbeDSini(g.adjL, g.adjR);
    vector<BFST> bfs(size_of_R);
    //generate the unilateral order
    //TODO:this could be a potential problem
    time_start = high_resolution_clock::now();
    BiCore bcore(g);

    if(order == 1){
        cout << "Inc";
        bcore.IncresingDegreeOrder(g.adjR);
    } else if(order == 2){
        cout << "Rand";
        bcore.NoOrder(g.adjR);
    } else{
        bcore.e2hopCore(g.adjR, g.adjL);
    }
    
    //bcore.IncresingDegreeOrder(g.adjR);
    //bcore.NoOrder(g.adjR);
    //end of problem
    // cout<<"core finish"<<endl;
    if (maxPruning) {
        //Future work, further refine the unilateral order to maximize pruning
        //This optimisation is for sequential algorithm,
        vector<pair<int, int>> auD;
        int initUniCore = bcore.uni_l_core[bcore.uniOrder[1]];
        int count = 0;
        for (auto &v: bcore.uniOrder) {
            if (bcore.uni_l_core[v] == initUniCore) {
                count++;
            } else {
                auD.emplace_back(make_pair(initUniCore, count));
                initUniCore = bcore.uni_l_core[v];
                count = 1;
            }
        }

        auD.emplace_back(make_pair(initUniCore, count));

        int ini = 0;
        for (auto &pair: auD) { // adjust pairs for initial position of very unilateral core.
            int tmp = pair.second;
            pair.second = ini;
            ini = ini + tmp;
        }

        vector<vector<int>> vertex_with_same_uni_cores;
        for (int i = 0; i < auD.size(); i++) {
            pair<int, int> cp = auD[i];
            vector<int> ertex_with_same_uni_core;
            if (i < auD.size() - 1) {
                pair<int, int> np = auD[i + 1];
                for (int j = cp.second; j < np.second; j++) {
                    ertex_with_same_uni_core.emplace_back(bcore.uniOrder[j]);
                }
                vertex_with_same_uni_cores.emplace_back(ertex_with_same_uni_core);
            }
            if (i == auD.size() - 1) {
                for (int j = cp.second; j < bcore.uniOrder.size(); j++) {
                    ertex_with_same_uni_core.emplace_back(bcore.uniOrder[j]);
                }
                vertex_with_same_uni_cores.emplace_back(ertex_with_same_uni_core);
            }
        }
        bcore.uniOrder.clear();
        for (auto &vec: vertex_with_same_uni_cores) {
            std::sort(vec.begin(), vec.end(),
                      bind(&ooMBE::orderSort, this, std::placeholders::_1, std::placeholders::_2));
            bcore.uniOrder.insert(bcore.uniOrder.end(), vec.begin(), vec.end());
        }


        //cout<<"test tag"<<endl;//the above code has been tested
    }

    //printVector(bcore.uniOrder);
    //test
    //if (bcore.uniOrder.size() != size_of_R) {
    //    cout << "Wrong!" << endl;
    //}
    //adjust R, IndexR, but adjIndexR keeps the same
    //vector<int> orderIndex(bcore.uniOrder.size());
    for (int i = 0; i < bcore.uniOrder.size(); i++) {
        int v = bcore.uniOrder[i];
        R[i] = v;
        indexR[v] = i;
        //orderIndex[v] = i;
    }


    // let refine it now

    //compute order preserved batch pivots
    //refresh that give vertex v in R, indexR[v] returns the location of v in R, where R is an array. The same for L
    //parameters adjL, adjR, the order, adjIndexL, adjIndexR
    //note that the order is for R
    //vector<int> pivots;
    vector<PivotsNeighbour> pivots_plus_neighbours;
    //orderBatchPivots(pivots);
    orderBatchPivots_neighbours(pivots_plus_neighbours);


    int iP = 0;//g.adjR.size();
    int iQ = iP;
    int iB = g.adjR.size();


    int iA = 0;

    vector<pair<int, int>> RP;
    for (int i = iP; i < g.adjR.size(); i++) {
        RP.emplace_back(make_pair(i, 0));
    }


    vector<pair<int, int>> RQ, RB;
    RQ.reserve(pivots_plus_neighbours.size());
    vector<pair<int, int>> LA;
    for (int i = 0; i < g.adjL.size(); i++) {
        LA.emplace_back(make_pair(i, 0));
    }


    for (int i = 0; i < size_of_L; i++) {
        adjIndexL[i] = 0;
    }

    for (int i = 0; i < size_of_R; i++) {
        adjIndexR[i] = 0;
    }


    //cout<<"still fine!"<<endl;
    //for(int i=iP;i<size_of_R;i++){// When backtracking, R should be the same as here
    // This can be achieved by using the stored RP
    //for(auto &pa: RP){
    //for(int i=0;i<RP.size();i++){
    //cout<<"max 2-hop core"<<bcore.maxE2Core<<endl;
    //max2Core = bcore.maxE2Core;
    //pivots_arry = new pair<int,int> [size_of_R];
    // cout<<"breadth: "<<pivots_plus_neighbours.size()<<endl;
    int count =0;


    int size_of_P;
    //for(int z=0;z<pivots_plus_neighbours.size();z++){
    for (int z=0;z<pivots_plus_neighbours.size();z++){
    //for (auto &vn: pivots_plus_neighbours) {
        //cout<<"Progress:"<<++count<<endl;
        //move v to iB
        PivotsNeighbour vn = pivots_plus_neighbours[z];
        int v = vn.vertex;
        //auto rightB = g.CommonNeighborsL(g.adjR[v]);
        /**
        if(count==10000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==20000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==30000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==40000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==50000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==60000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==70000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==80000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==90000){
            cout<<count<<":"<<nomb<<endl;
        }else
        if(count==3160000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==3150000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==3140000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==3130000){
            cout<<count<<":"<<nomb<<endl;
        }else if(count==3120000){
            cout<<count<<":"<<nomb<<endl;
        }else
        if(z==3244000) {
            cout << nomb << endl;
        }
        if(z==3245000){
            cout<<z<<":"<<nomb<<endl;
        }else if(z==3246000){
            cout<<z<<":"<<nomb<<endl;
        }else if(z==3247000){
            cout<<z<<":"<<nomb<<endl;
        }else if(z==3248000){
            cout<<z<<":"<<nomb<<endl;
        }else if(z==3249000){
            cout<<z<<":"<<nomb<<endl;
        }else if(z==3250000){
            cout<<z<<":"<<nomb<<endl;
        }
         **/
        /**
        if(count==3290001){
            cout<<"start"<<endl;
        }
        **/
        //if(z>3250000){
        //    cout<<vn.adj.size()<<", "<<z<<":"<<nomb<<endl;
        //}//else if(count==3200000){
        //    cout<<count<<":"<<nomb<<endl;
        //}
        //if()
        //count++;

        //cout<<"Processing "<<v<<endl;
        int iP_prime = 0, iA_prime = 0, iB_prime = 0, iQ_prime = 0;
        iB_prime = iB - 1;
        if (indexR[v] != iB - 1) {
            int tmpv = R[iB - 1];
            R[iB - 1] = v;
            R[indexR[v]] = tmpv;
            indexR[tmpv] = indexR[v];
            indexR[v] = iB - 1;
        }

        vector<pair<int, int>> A_prime;
        vector<pair<int, int>> Q_prime;
        vector<pair<int, int>> P_prime;

        //move neighbours of v in L
        //   compute starting position of the neighbour
        int nV = g.adjR[v].size();
        iA_prime = size_of_L - nV;
        int cPos = size_of_L - 1;
        for (auto &u: g.adjR[v]) {
            if (indexL[u] != cPos) {
                int tmpu = L[cPos];
                L[cPos] = u;
                L[indexL[u]] = tmpu;
                //int tmpPos = indexL[u];
                indexL[tmpu] = indexL[u];
                indexL[u] = cPos;
            }
            cPos--;
        }
        //iQ_prime = iQ;
        //Put vertex in Q into the correct position

        //shrank P
        //let's do it with O(E)
        //*****  ******
        //the loop below runs in O(E)
        //vector<pair<int,int>> counts;

        //TODO problem starts here! Fixed!
        //align 2-hop neighbours starting from iB-1
        //int ini_pos = iB_prime-1;
        ///vector<int> agg;
        //vector<int> notP;
        for(auto &v_p:vn.adj){
            //vector<int> localNeighbour;
            //vector<int> other;
        #ifdef EXIT_TIME
          if(getExcuteTime() >= EXIT_TIME)
          {
            printExcuteTime();
            PrintMemory();
            cout<<"processing_nodes: "<<processing_nodes<<endl;
            cout<<"maximal nodes: "<<nomb<<endl;
            exit(0);
          }
        #endif 
            int ini_local=max_degree-1;
            int ini_other=max_degree-1;
            for (int j = 0; j < g.adjR[v_p].size(); j++) {  // adjust the neighbours of v_prime
                int u = g.adjR[v_p][j];
                if (indexL[u] >= iA_prime) {
                    ln_for_vex[ini_local] = u;
                    ini_local--;
                    //localNeighbour.emplace_back(u);
                }else{
                    other_for_vex[ini_other] = u;
                    ini_other--;
                }
            }
            ini_local = ini_local+1;
            ini_other = ini_other+1;
            int adj_ini = g.adjR[v_p].size() - 1;
            //for (auto & u:localNeighbour){
            for(int j=ini_local;j<max_degree;j++){
                int u = ln_for_vex[j];
                g.adjR[v_p][adj_ini] = u;
                adj_ini--;
            }
            adjIndexR[v_p] = adj_ini + 1;
            for(int j=ini_other;j<max_degree;j++){
                int u = other_for_vex[j];
                g.adjR[v_p][adj_ini] = u;
                adj_ini--;
            }
            //if((max_degree-ini_local)==size_of_L - iA_prime){
                //agg.emplace_back(v_p);
                //cout<<"wrong"<<endl;
           // }else if((max_degree-ini_local)>0&&(max_degree-ini_local)!=(size_of_L - iA_prime)){
           P_prime.emplace_back(make_pair(v_p, adjIndexR[v_p]));
            //}

        }
        //iP_prime = ini_pos+1;

        int ini_pos = iB_prime-1;
        for(auto &v_b:vn.agg){//vgg: 新增的B
            if(indexR[v_b]!=ini_pos){
                int tmpv = R[ini_pos];
                R[ini_pos] = v_b;
                R[indexR[v_b]] = tmpv;
                indexR[tmpv] = indexR[v_b];
                indexR[v_b] = ini_pos;
            }
            ini_pos--;
        }
        //iP_prime = ini_pos+1;
        iB_prime = ini_pos+1;
        ini_pos = iB_prime-1;

        for(auto &v_pp:P_prime){
            int v_b = v_pp.first;
            if(indexR[v_b]!=ini_pos){
                int tmpv = R[ini_pos];
                R[ini_pos] = v_b;
                R[indexR[v_b]] = tmpv;
                indexR[tmpv] = indexR[v_b];
                indexR[v_b] = ini_pos;
            }
            ini_pos--;
        }
        iP_prime = ini_pos+1;
#ifdef THRESHOLD
        if(size_of_R - iB_prime >= THRESHOLD)continue;//lixu
#endif
        nomb++;
        level_accumulation_ += cur_level_;
        max_level_ = std::max(cur_level_, max_level_);
        //std::cout<<"----------"<<nomb<<"--------------"<<std::endl;
        //std::cout<<vn.agg.size()<<" "<<rightB.size()<<std::endl;

        //for(int i: rightB)std::cout<<i<<" ";
        //std::cout<<std::endl;
        //if(rightB.size()>=2)nomb--;

        if(P_prime.size()==0){
            continue;
        }

        //Align old Q into the correct position
        int iniPos_Q = iP_prime - 1;
        for (auto &vq:vn.Q){
            //int vq = vqp.first;
            if (indexR[vq] != iniPos_Q) {
                int tmpv = R[iniPos_Q];
                R[iniPos_Q] = vq;
                R[indexR[vq]] = tmpv;
                indexR[tmpv] = indexR[vq];
                indexR[vq] = iniPos_Q;
            }

            //adjust adj list for vq;
            iniPos_Q--;
        }

        iQ_prime = iniPos_Q + 1;


        //now shrank Q to Q_prime
        // vector<int> notQ;
        //for (int i = iQ_prime; i < iB_prime; i++) {
        for (int i = iQ_prime; i < iP_prime; i++) {
            int vq = R[i];
            int counts = 0;
            int ini_local=max_degree-1;
            int ini_other=max_degree-1;

            //vector<int> localNeighbours;
            //vector<int> other;
            for (int j = 0; j < g.adjR[vq].size(); j++) {
                int u = g.adjR[vq][j];
                if (indexL[u] >= iA_prime) {
                    counts++;
                    //localNeighbours.emplace_back(u);
                    ln_for_vex[ini_local] = u;
                    ini_local--;
                }else{
                    //other.emplace_back(u);
                    other_for_vex[ini_other] = u;
                    ini_other--;
                }
            }

            ini_local=ini_local+1;
            ini_other=ini_other+1;
            //vq.adj \cap A \ne emptyset
            int adj_ini = g.adjR[vq].size() - 1;
            for(int j = ini_local;j<max_degree;j++){
                int u = ln_for_vex[j];
                g.adjR[vq][adj_ini]=u;
                adj_ini--;
            }
            adjIndexR[vq] = adj_ini+1;
            for(int j= ini_other;j<max_degree;j++){
                int u = other_for_vex[j];
                g.adjR[vq][adj_ini]=u;
                adj_ini--;
            }

            if(counts>0){
                Q_prime.emplace_back(make_pair(vq,adjIndexR[vq]));
            }
            /**
            if (counts == 0) {
                //notQ.emplace_back(vq);
                cout<<"should not happlen!"<<endl;

            }
            **/
        }

        /**
        for (auto &vq: notQ) {
            int tmpv = R[iQ_prime];
            R[iQ_prime] = vq;
            R[indexR[vq]] = tmpv;
            indexR[tmpv] = indexR[vq];
            indexR[vq] = iQ_prime;
            iQ_prime = iQ_prime + 1;
        }
        **/

        //update adj index based on the current iQ, iP and iB
        for (int i = iA_prime; i < size_of_L; i++) {
            int u = L[i];

            //vector<int> local_neighbours;
            //vector<int> other;
            int ini_local=max_degree-1;
            int ini_other=max_degree-1;

            for (int j = 0; j < g.adjL[u].size(); j++) {
                int v_prime = g.adjL[u][j];
                if (indexR[v_prime] >= iQ_prime && indexR[v_prime] < iB_prime) {
                    ln_for_vex[ini_local]=v_prime;
                    ini_local--;
                }else{
                    other_for_vex[ini_other]=v_prime;
                    ini_other--;
                }
            }
            ini_local=ini_local+1;
            ini_other=ini_other+1;
            int adj_ini = g.adjL[u].size() - 1;
            for(int j=ini_local;j<max_degree;j++){
                int v = ln_for_vex[j];
                g.adjL[u][adj_ini]=v;
                adj_ini--;
            }
            adjIndexL[u] = adj_ini + 1;
            for(int j=ini_other;j<max_degree;j++){
                int v = other_for_vex[j];
                g.adjL[u][adj_ini]=v;
                adj_ini--;
            }
            A_prime.emplace_back(u, adjIndexL[u]);

            /*for(int j = 0; j < adjIndexL[u]; j++)
            {
              if(std::find(R+iQ_prime, R+iB_prime,g.adjL[u][j])!=R+iB_prime)
              {
                std::cout<<"error"<<std::endl;
              }
            }*/
        }



        if(P_prime.size()>0){
            cur_level_++;
            // advIMBEA(iA_prime, iQ_prime, iP_prime, iB_prime, A_prime,
            // P_prime, Q_prime);
            adv_order_PIMBEA_local(iA_prime, iQ_prime, iP_prime, iB_prime, A_prime, P_prime, Q_prime);
            cur_level_--;
        }

        //int iB = iB - 1;
        //size_of_P = size_of_P-1;
        //iP = iB-;
        //int last_pos_of_Q=count-1;
        //move v to position count
        /**
        if(indexR[v]!=last_pos_of_Q){
            int tmpv = R[last_pos_of_Q];
            R[last_pos_of_Q] = v;
            R[indexR[v]]=tmpv;
            indexR[tmpv] = indexR[v];
            indexR[v] = last_pos_of_Q;
        }
         **/
        //iQ=0;
        //RQ.emplace_back(make_pair(v,0));
        //iP=last_pos_of_Q+1;
        //cout<<iP<<endl;
        /**
        }else{
            int ini_pos = iB - 1;
            //for(auto &pair:RP){
            int loc_of_v_in_RP;
            for (int i = 0; i < RP.size(); i++) {
                pair<int, int> p = RP[i];
                int bv = p.first;
                if (bv == v) {
                    loc_of_v_in_RP = i;
                }
                adjIndexR[bv] = p.second;
                if (indexR[bv] != ini_pos) {
                    int tmpv = R[ini_pos];
                    R[ini_pos] = bv;
                    R[indexR[bv]] = tmpv;
                    indexR[tmpv] = indexR[bv];
                    indexR[bv] = ini_pos;
                }
                ini_pos--;
            }
            //remove v form RP
            iP = ini_pos + 1;
            if(RP.size()>0){//this is always true
                if(RP.size()>=2){

                    pair<int, int> tmpp = RP[RP.size() - 1];
                    pair<int, int> tmpvp = RP[loc_of_v_in_RP];
                    RP[loc_of_v_in_RP] = tmpp;
                    RP.pop_back();
                    RQ.emplace_back(tmpvp);
                }else if (RP.size()==1) {
                    RQ.emplace_back(RP[0]);
                    RP.pop_back();
                }

                //refine array and indices
                int b_loc_of_v = indexR[v];
                R[b_loc_of_v] = R[iP];
                indexR[R[iP]] = b_loc_of_v;
                R[iP] = v;
                indexR[v] = iP;
                iP = iP + 1;
                int ini_pos_Q = iP - 1;
                for (auto &p: RQ) {
                    if (indexR[p.first] != ini_pos_Q) {
                        int tmpv = R[ini_pos_Q];
                        R[ini_pos_Q] = p.first;
                        R[indexR[p.first]] = tmpv;
                        indexR[tmpv] = indexR[p.first];
                        indexR[p.first] = ini_pos_Q;
                    }
                    ini_pos_Q--;
                    adjIndexR[p.first] = 0;
                }
                iQ = ini_pos_Q + 1;
                //if(iQ!=0){
                //    cout<<"Hehehehe"<<endl;
                // }
                //cout<<iQ<<endl;
                for(int i=iA_prime;i<size_of_L;i++){
                    int u = L[i];
                    adjIndexL[u]=0;
                }
                //fill(adjIndexL+(iA_prime-1), adjIndexL + size_of_L, 0);
            }

            cout<<iP<<endl;
        }
        **/
        //cout << "testint" << endl;
        //back track
        //align P,Q,B back first




        //for()
        //adj align back

        /***
        for(int j=i+1;j<RP.size();j++){
            int v_prime = RP[j].first;
            //int number_of_intersection = 0;
            int adj_v_start = adjIndexR[v_prime];
            int iniPos = g.adjR[v_prime].size()-1;
            for(int z=adj_v_start;z<g.adjR[v_prime].size();z++){
                int u = g.adjR[v_prime][z];
                if(indexL[u]>=iA_prime){
                    //move it to the iniPos
                    if(z<iniPos){
                        int tmpu =  g.adjR[v_prime][iniPos];
                        g.adjR[v_prime][iniPos] = u;
                        g.adjR[v_prime][z] = tmpu;

                    }
                    iniPos--;
                }
            }
            //iniPos++ is the adjIndexR[v_prime]
            int adj_v_prime_update = iniPos+1;
            adjIndexR[v_prime] = adj_v_prime_update;
            int intersection = g.adjR[v_prime].size()-adj_v_prime_update;
            if(intersection==g.adjR[v_prime].size()){
                //aggressive expansion
            }
        }
        ***/
        //cout<<"test"<<endl;

    }

    /**
     * vector<int> bfs(size_of_L,-1);
    for (auto & u:bcore.uniOrder){//loop 1
        cout<<u<<",";
        vector<pair<int,int>> A_prime;
        vector<pair<int,int>> Q_prime;
        vector<pair<int,int>> P_prime;
        //compute induced subgraph
        int loc_of_u=indexL[u];
        vertex_swap(loc_of_u,size_of_L-1,L,indexL);
        int iB_prime=size_of_L;
        //int iP_prime=
        //algin neighbours of u to A
        int ini_loc = size_of_R-1;
        for(auto &v:g.adjL[u]){
            int loc_of_v = indexR[v];
            if(loc_of_v!=ini_loc){
                vertex_swap(loc_of_v,ini_loc,R,indexR);
            }
            ini_loc--;
        }
        int iA_prime=ini_loc+1; //starting position of iA_prime
        //prepare P, i.e., vertices
        ini_loc = iB_prime-1;

        for(int loc=iA_prime;loc<size_of_R;loc++){
            int v = R[loc];
            for(auto &u_prime: g.adjR[v]){
                if(bcore.index_for_uniOrder[u_prime]>bcore.index_for_uniOrder[u]){// appears after u in the unilateral order
                    if(bfs[u_prime]!=u){
                        //first time it visits
                        bfs[u_prime] = u;
                        int loc_u_prime = indexL[u_prime];
                        if(loc_u_prime!=ini_loc){
                            vertex_swap(loc_u_prime,ini_loc,L,indexL);
                        }
                        ini_loc--;
                    }
                }

            }
        }
        int iP_prime = ini_loc+1;

        //adjust adjIndexR and adjust adjIndexL
        for(int loc=iA_prime;loc<size_of_R;loc++){
            int v_prime = R[loc];
            for(auto &u_prime:g.adjR[v_prime]){
                int loc_u_prime = indexL[u_prime];
                if(loc_u_prime>=iP_prime&&loc_u_prime<iB_prime){

                }
            }

        }

        for(int loc=iP_prime;loc<iB_prime;loc++){

        }


        int iB=size_of_L;


        //need

    }
    **/

    //if(R!= 0) {delete [] R;R=0;}
    printExcuteTime();
    // std::cout<<"Collection time: "<<collection_time.count()<<endl;
}



