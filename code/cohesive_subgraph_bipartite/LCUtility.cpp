

#include "LCUtility.h"
#include <unordered_map>
using namespace std;

//livejournal offset 5204176
//l: 1 to 5204176 [5204177]
//r: 5204177 to 12693249 [7489074]

//all vertices are indexed from 1

//int adjM[1000000000];



//livejournal |L|3201203
//degree for L, index=vid, value=deg


void mbsToMBB(){

}

void computeIG(vector<int> & L,vector<int> & R, bgraph &og, bgraph &ng, int flag){
    if(flag==0){
        unordered_map<int,int> oldLtoNewL;
        vector<int> newLtoOldL(L.size());
        for(int i=0;i<L.size();i++){
            newLtoOldL[i] = L[i];
            oldLtoNewL.emplace(L[i],i);
        }
        //cout<<"hash L size"<<oldLtoNewL.size()<<endl;

        unordered_map<int,int> oldRtoNewR;
        vector<int> newRtoOldR(R.size());
        for(int i=0;i<R.size();i++){
            newRtoOldR[i] = R[i];
            oldRtoNewR.emplace(R[i],i);
        }

        vector<vector<int>> newAdjL(L.size());
        vector<vector<int>> newAdjR(R.size());

        for(auto& u:L){
            int n_id_u = oldLtoNewL[u];
            vector<int> *adjs_of_u = &(og.adjL[u]);//get u's adjlist from its original graph this is cpoy, need to avoid
            vector<int> new_adj_of_u;
            for(auto& v:*adjs_of_u){ // those  are R
                if(oldRtoNewR.find(v)!=oldRtoNewR.end()){ //
                    int new_id_v = oldRtoNewR[v];
                    new_adj_of_u.emplace_back(new_id_v);// This may not
                }
            }
            newAdjL[n_id_u] = new_adj_of_u;
        }


        for(auto&v:R){
            int new_id_v = oldRtoNewR[v];
            vector<int> *adjs_of_v = &(og.adjR[v]); //copy
            vector<int> new_adj_of_v;
            for(auto& u:*adjs_of_v){ //those vertices are in L
                if(oldLtoNewL.find(u)!=oldLtoNewL.end()){
                    int new_id_u = oldLtoNewL[u];
                    new_adj_of_v.emplace_back(new_id_u);
                }
            }
            newAdjR[new_id_v] = new_adj_of_v;
        }

        ng.adjL = newAdjL;
        ng.adjR = newAdjR;
        ng.novL = newAdjL.size();
        ng.novR = newAdjR.size();

    }else{
        //L is R
        unordered_map<int,int> oldLtoNewL;
        vector<int> newLtoOldL(L.size());
        for(int i=0;i<L.size();i++){
            newLtoOldL[i] = L[i];
            oldLtoNewL.emplace(L[i],i);
        }
        //cout<<"hash L size"<<oldLtoNewL.size()<<endl;

        unordered_map<int,int> oldRtoNewR;
        vector<int> newRtoOldR(R.size());
        for(int i=0;i<R.size();i++){
            newRtoOldR[i] = R[i];
            oldRtoNewR.emplace(R[i],i);
        }


        vector<vector<int>> newAdjL(L.size());
        vector<vector<int>> newAdjR(R.size());


        for(auto& u:L){
            int n_id_u = oldLtoNewL[u];
            vector<int> *adjs_of_u = &(og.adjR[u]);//get u's adjlist from its original graph
            vector<int> new_adj_of_u;
            for(auto& v:*adjs_of_u){ // those  are R
                if(oldRtoNewR.find(v)!=oldRtoNewR.end()){ //
                    int new_id_v = oldRtoNewR[v];
                    new_adj_of_u.emplace_back(new_id_v);// This may not
                }
            }
            newAdjL[n_id_u] = new_adj_of_u;
        }


        for(auto&v:R){
            int new_id_v = oldRtoNewR[v];
            vector<int> *adjs_of_v = &(og.adjL[v]);
            vector<int> new_adj_of_v;
            for(auto& u:*adjs_of_v){ //those vertices are in L
                if(oldLtoNewL.find(u)!=oldLtoNewL.end()){
                    int new_id_u = oldLtoNewL[u];
                    new_adj_of_v.emplace_back(new_id_u);
                }
            }
            newAdjR[new_id_v] = new_adj_of_v;
        }

        ng.adjL = newAdjR;
        ng.adjR = newAdjL;
        ng.novL = newAdjR.size();
        ng.novR = newAdjL.size();

    }

}


void compute2IG(int u,  bgraph &og, bgraph &ng, int flag){
    //compute 2-hop vertex set
    if(flag==0){
        //flag ==0 u is in L
        //1-hop vertex set
        vector<int> R = og.adjL[u];

        unordered_set<int> L;
        L.emplace(u);
        for(auto &v:R){
            vector<int> list = og.adjR[v];
            for(auto &ul:list){
                L.emplace(ul);
            }
        }
        vector<int> nL;
        nL.insert(nL.end(),L.begin(),L.end());
        computeIG(nL,R,og,ng,0);

    }else{
        // u is in R
        vector<int> L = og.adjR[u];
        unordered_set<int> R;
        R.emplace(u);
        for(auto &v:L){
            vector<int> list = og.adjL[v];
            for(auto &ul:list){
                R.emplace(ul);
            }
        }

        vector<int>nR;
        nR.insert(nR.end(),R.begin(),R.end());
        computeIG(L,nR,og,ng,0);
    }

    //flag ==1 u is in R
}

void computeIG(vector<int> & L,vector<int> & R, bgraph &og, bgraph &ng){
    //processing L
    /**
     * for(auto &v:L){
        cout<<v<<",";
    }
    cout<<endl;

    for(auto &v:R){
        cout<<v<<",";
    }
    cout<<endl;
     */


    unordered_map<int,int> oldLtoNewL;
    vector<int> newLtoOldL(L.size());
    for(int i=0;i<L.size();i++){
        newLtoOldL[i] = L[i];
        oldLtoNewL.emplace(L[i],i);
    }
    //cout<<"hash L size"<<oldLtoNewL.size()<<endl;

    unordered_map<int,int> oldRtoNewR;
    vector<int> newRtoOldR(R.size());
    for(int i=0;i<R.size();i++){
        newRtoOldR[i] = R[i];
        oldRtoNewR.emplace(R[i],i);
    }


    vector<vector<int>> newAdjL(L.size());
    vector<vector<int>> newAdjR(R.size());

    vector<unordered_set<int>> newAdjHL(L.size());
    vector<unordered_set<int>> newAdjHR(R.size());

    for(auto& u:L){
        int n_id_u = oldLtoNewL[u];
        vector<int> adjs_of_u = og.adjL[u];//get u's adjlist from its original graph
        vector<int> new_adj_of_u;
        unordered_set<int> new_adjH_of_u;

        for(auto& v:adjs_of_u){ // those  are R
            if(oldRtoNewR.find(v)!=oldRtoNewR.end()){ //
                int new_id_v = oldRtoNewR[v];
                new_adj_of_u.emplace_back(new_id_v);// This may not
                new_adjH_of_u.emplace(new_id_v);
            }
        }
        newAdjL[n_id_u] = new_adj_of_u;
        newAdjHL[n_id_u] = new_adjH_of_u;
    }


    for(auto&v:R){
        int new_id_v = oldRtoNewR[v];
        vector<int> adjs_of_v = og.adjR[v];
        vector<int> new_adj_of_v;
        unordered_set<int> new_adjH_of_u;
        for(auto& u:adjs_of_v){ //those vertices are in L
            if(oldLtoNewL.find(u)!=oldLtoNewL.end()){
                int new_id_u = oldLtoNewL[u];
                new_adj_of_v.emplace_back(new_id_u);
                new_adjH_of_u.emplace(new_id_u);
            }
        }
        newAdjR[new_id_v] = new_adj_of_v;
        newAdjHR[new_id_v] = new_adjH_of_u;
    }

    ng.adjL = newAdjL;
    ng.adjR = newAdjR;
    ng.adjHL = newAdjHL;
    ng.adjHR = newAdjHR;
    ng.novL = newAdjL.size();
    ng.novR = newAdjR.size();


}




void computeIG(vector<int> & V,bgraph &og, bgraph &ng){

    unordered_map<int,int> oldVtoNewV; //old id as index
    vector<int> newVtoOldV(V.size());
    int count=0;
    for(int i=0;i<V.size();i++){
        newVtoOldV[i] = V[i];
        oldVtoNewV.emplace(V[i],i);
        if(V[i]<og.offset){ //offset is number of vertices in |L| , id of L starts from 0
            count++;
        }
    }


    vector<vector<int>> newAdjL(count);
    vector<vector<int>> newAdjR(V.size()-count);
    for(auto& u:V){
        //access its adjList
        if(u<og.offset){ //vertices in L
            vector<int> list = og.adjL[u];
            vector<int> newList;
            for(auto &v:list){ // vertices in R

                if(oldVtoNewV.find(v+og.offset)!=oldVtoNewV.end()){
                    newList.emplace_back(oldVtoNewV[v]); //put its new id into new adj list
                }

            }
            newAdjL[oldVtoNewV[u]] = newList;
        }else{

            vector<int> list = og.adjR[u];



        }

    }


}

void coreFilter(){
    ifstream myfile ("../core_decomposition.txt");
    string line, v, cnum;
    int coreTh=170;
    int count=0;
    unordered_set<int> bigCore;
    if (myfile.is_open()){
        while ( getline(myfile, line) ) {
            stringstream sline(line);
            getline(sline, v, '\t');
            getline(sline, cnum, '\t');
            int icum = stoi(cnum);
            int iv = stoi(v);
            if(icum>=coreTh){
                count++;
                bigCore.emplace(iv);
            }
        }
    }
    myfile.close();
    myfile.open ("../eedges_livejournal_new_id.txt");
    string sL,sR;
    vector<Edge> edgeList;
    if (myfile.is_open()){
        while ( getline(myfile, line) ) {
            stringstream sline(line);
            getline(sline, sL, ',');
            getline(sline, sR, ',');
            int l = stoi(sL);
            int r = stoi(sR)+3201203;
            if(bigCore.find(l)!=bigCore.end() && bigCore.find(r)!=bigCore.end() ){

                struct Edge e;
                e.l = l;
                e.r = r-3201203;
                edgeList.push_back(e);

            }

            //degL[l]++;
            //degR[r+3201203]++;
        }
    }
    myfile.close();

    ofstream fout;
    fout.open("../livejournal_170_edges.txt");
    for(auto &e:edgeList){
        fout<<e.l<<" "<<e.r<<"\n";
    }
    fout.close();
    //next


}

void eLtoAdjL(){

    ifstream myfile ("../livejournal.edges");
    string line;
    string sL,sR;

    vector<vector<int>> adjList;
    for(int i=0;i<3201203+7489073;i++){
        vector<int> list;
        adjList.push_back(list);
    }
    cout<<"start to read"<<endl;
    if (myfile.is_open()){
        while ( getline(myfile, line) ) {
            stringstream sline(line);
            getline(sline, sL, ',');
            getline(sline, sR, ',');
            int l = stoi(sL);
            int r = stoi(sR);
            adjList[l].push_back(r);
            adjList[r+3201203].push_back(l);
            //degL[l]++;
            //degR[r+3201203]++;
        }
    }
    myfile.close();

    ofstream fout;
    fout.open("../livejournal_adjL_deg.txt");
    for(int i=0;i<3201203+7489073;i++){
        vector<int> list = adjList[i];
        fout<<i<<","<<list.size();
        for(int j=0;j<list.size();j++){
            fout<<","<<list[j];
        }
        fout<<"\n";
    }
    fout.close();
}



void edgeFilter(string r, string s,int offset){

    ifstream myfile(r);
    string line;
    string sni, soi;
    ofstream fout;
    fout.open(s);

    if (myfile.is_open()){
        while ( getline(myfile, line) ) {
            stringstream sline(line);
            getline(sline, sni, ',');
            getline(sline, soi, ',');
            int ni = stoi(sni);
            int oi = stoi(soi);
            if(ni<=offset&&oi<=offset){
                continue;
            }
            fout<<sni<<" "<<soi<<"\n";
        }

    }
    myfile.close();
    fout.close();
}

void newEdgeList(string r, string s, int offset){
    ifstream myfile(r);
    string line;
    string sni, soi;
    vector<Edge> edgeList;
    unordered_set<int> L, R;
    if (myfile.is_open()){
        //getline(myfile, line); //skip the first line
        while ( getline(myfile, line) ) {
            stringstream sline(line);
            getline(sline, sni, ' ');
            getline(sline, soi, ' ');
            int ni = stoi(sni);
            int oi = stoi(soi);
            Edge e;
            //if(ni>offset){
            //    cout<<"unexpected l!"<<endl;
            //    break;
            //}
            e.l = ni;
            L.emplace(e.l);

            //if(oi<=offset){
            //    cout<<"unexpected r!"<<endl;
            //    break;
           // }

            e.r = oi-offset;
            R.emplace(e.r);
            edgeList.emplace_back(e);
        }
    }
    myfile.close();

    //now vertices are from 1:n1 and 1:n2
    vector<vector<int>> adjL (L.size());
    vector<vector<int>> adjR (R.size());

    //could be dangerous
    int * revL = new int [L.size()+1]();
    int * revR = new int [R.size()+1]();


    int i=0;
    for(auto &u:L){
        revL[u]=i;
        i++;
    }

    int j=0;
    for(auto &v:R){
        revR[v] =j;
        j++;
    }


    for(auto& e:edgeList){
        int nl = revL[e.l];
        int nr = revR[e.r];
        adjL[nl].emplace_back(nr);
        adjR[nr].emplace_back(nl);
    }


    //shall check biclique
    //write to adjfile

    int newoffset = L.size();

    ofstream foutl,foutr,foutall;
    foutall.open("../klj.txt");

    foutl.open("../kljL.txt");
    for(int i=0;i<adjL.size();i++){
        foutl<<i;
        foutall<<i;
        for(auto & v: adjL[i]){
            foutl<<","<<v;
            foutall<<","<<v+newoffset;
        }
        foutl<<endl;
        foutall<<endl;
    }
    foutl.close();

    foutr.open("../kljR.txt");
    for(int j=0;j<adjR.size();j++){
        foutr<<j;
        foutall<<j+newoffset;
        for(auto & u: adjR[j]){
            foutr<<","<<u;
            foutall<<","<<u;
        }
        foutr<<endl;
        foutall<<endl;
    }

    foutr.close();


    foutall.close();
    delete []revL;
    delete []revR;

}


void newEdgelist(){
    //l:3201202
    int *lnto = new int[3201203];
    int *lotn = new int[3201203+1];
    //r:7489072
    int *rnto = new int[7489073];
    int *rotn = new int[7489073+1];

    //read left
    ifstream myfile ("../lvnto.txt");
    string line;
    string sni, soi;
    if (myfile.is_open()){
        while ( getline(myfile, line) ) {
            stringstream sline(line);
            getline(sline, sni, ',');
            getline(sline, soi, ',');
            int ni = stoi(sni);
            int oi = stoi(soi);
            lnto[ni]=oi;
        }

    }
    myfile.close();

    myfile.open("../lvotn.txt");
    if (myfile.is_open()){
        while ( getline(myfile, line) ) {
            stringstream sline(line);
            getline(sline, soi, ',');
            getline(sline, sni, ',');
            int oi = stoi(soi);
            int ni = stoi(sni);
            lotn[oi]=ni;
        }
    }
    myfile.close();

    myfile.open("../rvnto.txt");
    if (myfile.is_open()){
        while ( getline(myfile, line) ) {
            stringstream sline(line);
            getline(sline, sni, ',');
            getline(sline, soi, ',');
            int ni = stoi(sni);
            int oi = stoi(soi);
            rnto[ni]=oi;
        }

    }
    myfile.close();

    myfile.open("../rvotn.txt");
    if (myfile.is_open()){
        while ( getline(myfile, line) ) {
            stringstream sline(line);
            getline(sline, soi, ',');
            getline(sline, sni, ',');
            int oi = stoi(soi);
            int ni = stoi(sni);
            rotn[oi]=ni;
        }

    }
    myfile.close();

    myfile.open("../eedges_live_journal.txt");
    vector<Edge> edgeList;
    string ls, rs;
    if (myfile.is_open()){
        while ( getline(myfile, line) ) {
            stringstream sline(line);
            getline(sline, ls, ',');
            getline(sline, rs, ',');
            int l = stoi(ls);
            int r = stoi(rs)-5204176;
            struct Edge e;
            e.l = lotn[l];
            e.r = rotn[r];
            edgeList.push_back(e);
        }
    }
    myfile.close();

    //writ to a new file
    ofstream fout;
    fout.open("../livejournal_new_id_space.txt");
    for(auto& e: edgeList){
        fout<<e.l<<" "<<e.r<<"\n";
    }
    fout.close();
}


//live journal: group id: 5204177 to xx
//live journal: people id: 1 to 5204176

void cleanUpLJ(){
    string line;
    ifstream myfile ("../livejournal.edges");
    vector<Edge> edgeList;
    string ls, rs;
    unordered_set<int> L;
    unordered_set<int> R;
    if (myfile.is_open()){
        while ( getline(myfile, line) ) {
            stringstream sline(line);
            getline(sline, ls, ',');
            getline(sline, rs, ',');
            int l = stoi(ls);
            int r = stoi(rs);
            if(l>5204176){
                cout<<"wrong";
            }
            if(r>5204176){
                struct Edge e;
                e.l = l-1;
                e.r = r-5204177;
                L.emplace(e.l);
                R.emplace(e.r);
                edgeList.push_back(e);
            }

        }
    }
    myfile.close();
    cout<<"Number of vertices in L:"<<L.size()<<endl;
    cout<<"Number of vertices in R:"<<R.size()<<endl;
    vector<vector<int>> adj;
    for(int i=0;i<L.size()+R.size();i++){
        vector<int> list;
        adj.push_back(list);
    }
    for(auto&e:edgeList){
        int r = e.r+L.size();
        adj[e.l].push_back(r);
        adj[r].push_back(e.l);
    }


    ofstream fout;
    fout.open("../dataset_1_12/adjLJ.txt");
    fout<<L.size()<<","<<R.size()<<endl;
    for(int i=0;i<adj.size();i++){
        vector<int> list =adj[i];
        fout<<i;
        for(auto &v:list){
            fout<<","<<v;
        }
        fout<<endl;
    }
    fout.close();


}

void cleanUp(){
    string line;
    ifstream myfile ("../livejournal.edges");
    vector<Edge> edgeList;
    string ls, rs;
    cout<<"start to read"<<endl;
    int lsize = 5204177;
    int rsize= 7489074;
    //int *lidex = new int[lsize];
    //int *index = new int[rsize];

    //unordered_set <int> lvs;
    unordered_set <int> rvs;
    if (myfile.is_open()){
        while ( getline(myfile, line) ) {
            stringstream sline(line);
            getline(sline, ls, ',');
            getline(sline, rs, ',');
            int l = stoi(ls);
            int r = stoi(rs);
            struct Edge e;
            e.l = l;
            e.r = r;
            edgeList.push_back(e);
        }
    }
    myfile.close();
    cout<<"read finish"<<endl;
    for(auto& e: edgeList){
        rvs.emplace(e.r);
        //lvs.insert(e.l);
    }
    cout<<rvs.size();
    //int *lnto = new int[lvs.size()];  //new id as index
    int *rnto = new int[rvs.size()];  //new id as index
    //int *lotn = new int[lvs.size()+1];  //old id as index
    int *rotn = new int[rvs.size()+1];  //old id as index
    //int i=0;
    //for (const auto& elem: lvs) {
    //    lnto[i]=elem;
     //   lotn[elem] =i;
   //     i++;
   // }
    int i=0;
    for (const auto& elem: rvs) {
        rnto[i]=elem;
        rotn[elem] = i;
        i++;
    }

    ofstream fout;
    fout.open("../rvnto.txt");
    for(i=0;i<rvs.size();i++){
        fout<<i<<","<<rnto[i]<<"\n";
    }
    fout.close();
    ofstream fout2;
    fout2.open("../rvotn.txt");
    for(i=0;i<rvs.size()+1;i++){
        fout2<<i<<","<<rotn[i]<<"\n";
    }
    fout2.close();
}





void eLtoAdjM(){
    //auto adjM = new int [520417][748907];
    cout<<"finish alloc"<<endl;

    ifstream myfile ("../livejournal_170_edges.txt");
    string line;
    vector<Edge> edgeList;
    string ls, rs;
    cout<<"start to read"<<endl;
    //unordered_set<int> vL;
    unordered_map<int,int> vL;
    unordered_map<int,int> vR;
    int il=0, ir=0;
    if (myfile.is_open()){
        while ( getline(myfile, line) ) {
            stringstream sline(line);
            getline(sline, ls, ' ');
            getline(sline, rs, ' ');
            int l = stoi(ls);
            int r = stoi(rs);
            if(vL.find(l)==vL.end()){
                vL.emplace(l,il);
                il++;
            }
            if(vR.find(r)==vR.end()){
                vR.emplace(r,ir);
                ir++;
            }

            struct Edge e;
            e.l = vL[l];
            e.r = vR[r];
            edgeList.push_back(e);
        }
    }
    myfile.close();

    int **adjM = new int*[vL.size()];
    for(int i=0;i<vL.size();i++){
        adjM[i] = new int[vR.size()];

    }


    for(auto &e: edgeList){
        adjM[e.l][e.r]=1;
    }
    //auto adjM = new int[vL.size()][vR.size()];


    cout<<"start to write"<<endl;
    //create adj matrix
    ofstream fout;
    fout.open("../adjm_lj_170.txt");
    fout<<vL.size()<<","<<vR.size()<<endl;
    for(int i=1;i<vL.size();i++){
        for(int j=1;j<vR.size();j++){
            if(j==vR.size()-1){
                fout<<adjM[i][j]<<"\n";
            }else{
                fout<<adjM[i][j]<<" ";
            }
        }
    }
    fout.close();
}

void adjMtoAdj(){

    vector<vector<int>> adjL; //591
    vector<vector<int>> adjR; //22311
    for(int i=0;i<591;i++){
        vector<int>list;
        adjL.push_back(list);
    }
    for(int j=0;j<22311;j++){
        vector<int>list;
        adjR.push_back(list);
    }
    ifstream myfile ("../adjm_livejournal_170.txt");
    string line;
    if (myfile.is_open()){
        int l=0;
        while (getline(myfile, line)) {
            stringstream sline(line);
            string ele;
            int r=0;
            while(getline(sline, ele, ' ')){
                int e = stoi(ele);
                if(e==1){
                    adjL[l].push_back(r);
                    adjR[r].push_back(l);

                }
                r++;
            }
            l++;
        }
    }
    myfile.close();


    ofstream fout;
    fout.open("../ljadjL_170.txt");
    for(int i=0;i<591;i++){
        vector<int> list =adjL[i];
        fout<<i;
        for(auto& v:list){
            fout<<","<<v;
        }
        fout<<endl;
    }
    fout.close();

    fout.open("../ljadjR_170.txt");
    for(int i=0;i<22311;i++){
        vector<int> list =adjR[i];
        fout<<i;
        for(auto& v:list){
            fout<<","<<v;
        }
        fout<<endl;
    }


    fout.close();
}



void readEdgeList(string &file, bgraph &g) {
    ifstream myfile (file);
    string line,u,v;
    vector<pair<int,int>> edge;
    //unordered_set<int> vsL,vsR;
    int numL = 0;
    int numR = 0;
    if (myfile.is_open()){
        while ( getline(myfile, line) ) {
            stringstream sline(line);
            getline(sline, u, ' ');
            getline(sline, v, ' ');
            int ul = stoi(u)-1;
            int vr = stoi(v)-1;
            if(ul>numL){
                numL = ul;
            }
            if(vr>numR){
                numR=vr;
            }
            edge.emplace_back(make_pair(ul,vr));
        }
    }

    vector<vector<int>> adjL(numL+1);
    vector<vector<int>> adjR(numR+1);
    for(auto &e:edge){
        adjL[e.first].emplace_back(e.second);
        adjR[e.second].emplace_back(e.first);
    }
    g.adjL = adjL;
    g.adjR = adjR;

    myfile.close();
}
