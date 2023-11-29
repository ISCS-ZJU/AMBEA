//
// Created by Luis on 24/09/2020.
//

#include "bgraph.h"
#include <set>
#include "Utility.h"
using namespace std;

bgraph::bgraph() {
    novL = 0;
    novR = 0;
    noE = 0;
    offset = 0;
}

void bgraph::readAdj(string  path) {
    ifstream myfile(path);
    string line;
    int maxR = -1;
    if (myfile.is_open()) {
        int count = 0;
        while (getline(myfile, line)) {
            stringstream sline(line);
            string ele;
            vector<int> list;
            while(getline(sline, ele, ' ')) {
                int iele = stoi(ele);
                list.push_back(iele);
                maxR = maxR > iele ? maxR : iele;
                noE ++;
            }
            //adjL.push_back(list);
            adjL.emplace_back(list);
        }
    }
    adjR.resize(maxR + 1);
    for(int i = 0; i < adjL.size(); i++)
    {
      for(auto r: adjL[i])
      {
        adjR[r].push_back(i);
      }
    }
    //sortAdj();
    //checkAdj();
    novL = adjL.size();
    novR = adjR.size();
    max_degree_L = maxDegree(adjL);
    max_degree_R = maxDegree(adjR);
    myfile.close();
    // cout<<novL<<" "<<novR<<" "<<noE<<endl;
}


int bgraph::edgeNumber(){
    if(noE == 0){
        for(int i=0;i<adjL.size();i++){
            vector<int> list = adjL[i];
            for(int j=0; j<list.size();j++){
                noE++;
            }
        }
    }
    return noE;
}


void bgraph::checkAdj(){
    bool change = false;
    for(int i=0;i<adjL.size();i++){

        vector<int> list = adjL[i]; // nodes in R

        for(int j=0;j<list.size();j++){
            int u = list[j];
            if(!binary_search(adjR[u].begin(),adjR[u].end(),i)){
                adjR[u].emplace_back(i);
                sort(adjR[u].begin(),adjR[u].end());
                change = true;
            }
        }
    }
    cout<<change<<endl;
    for(int i=0;i<adjR.size();i++){
        vector<int> list = adjL[i]; // nodes in L
        for(int j=0;j<list.size();j++){
            int v = list[j];
            if(!binary_search(adjL[v].begin(),adjL[v].end(),i)){
                adjL[v].emplace_back(i);
                sort(adjL[v].begin(),adjL[v].end());
                change = true;
            }

        }


    }
    cout<<change<<endl;
}

std::vector<int> bgraph::CommonNeighborsL(const std::vector<int> &L)
{
  //if(L.size() == 0)std::cout<<"commonneighborsl error"<<std::endl;
  std::vector<int>result;
  for(int i = 0; i < L.size(); i++)
  {
    auto adj = adjL[L[i]];
    sort(adj.begin(), adj.end());
    if(i == 0)
    {
      result = adj;
      continue;
    }
    result = seq_intersect(result, adj);
  }
  return result;
}

std::vector<int> bgraph::CommonNeighborsR(const std::vector<int> &R)
{
  //if(R.size() == 0)std::cout<<"commonneighborsr error"<<std::endl;
  std::vector<int>result;
  for(int i = 0; i < R.size(); i++)
  {
    auto adj = adjR[R[i]];
    sort(adj.begin(), adj.end());
    if(i == 0)
    {
      result = adj;
      continue;
    }
    result = seq_intersect(result, adj);
  }
  return result;
}


void bgraph::readAdj(int lr,string s){
    //test = new int [1];
    ifstream myfile(s);
    string line;
    if (myfile.is_open()) {
        while (getline(myfile, line)) {
            stringstream sline(line);
            string ele;
            vector<int> list;
            //getline(sline, ele, ','); //skip the first vertex
            while(getline(sline, ele, ' ')) {
                list.push_back(stoi(ele));
            }
            if (lr == 0) {
                adjL.push_back(list);
                //novL=adjL.size();
            }
            if (lr == 1) {
                adjR.push_back(list);
                //novR=adjR.size();
            }

        }
    }

    //sortAdj();


    novL = adjL.size();
    novR = adjR.size();
    calculate_number_of_edges();
   max_degree_L = maxDegree(adjL);
    max_degree_R = maxDegree(adjR);


    myfile.close();
}






void bgraph::dispalyAdj() {
    if(adjL.size()!=0||adjR.size()!=0){
        cout<<"The adj list for L"<<endl;
        for(int i=0; i<adjL.size();i++){
            cout<<i<<":";
            for(int j=0; j<adjL[i].size(); j++){
                cout<<adjL[i][j]<<",";
            }
            cout<<endl;
        }
        cout<<"The adj list for R"<<endl;
        for(int i=0; i<adjR.size();i++){
            cout<<i<<":";
            for(int j=0; j<adjR[i].size(); j++){
                cout<<adjR[i][j]<<",";
            }
            cout<<endl;
        }
    }else{
        cout<<"This bipartite graph is an empty set!";
    }

}

void bgraph::sortAdj() {
    for(auto &list:adjL){
        sort(list.begin(),list.end());
    }
    for(auto &list:adjR){
        sort(list.begin(),list.end());
    }
}

void bgraph::writAdj(int flag, string r) {
    ofstream fout;
    fout.open(r);
    if(flag==0){
        for(int i=0;i<adjL.size();i++){
            fout<<i;
            for(auto &v:adjL[i]){
                fout<<","<<v;
            }
            fout<<endl;
        }
    }else{
        for(int j=0;j<adjR.size();j++){
            fout<<j;
            for(auto& u:adjR[j]){
                fout<<","<<u;
            }
            fout<<endl;
        }
    }
    fout.close();
}

void bgraph::iniHash() {
    //cout<<"!!!"<<endl;
    for(int i=0;i<adjL.size();i++){
        unordered_set<int> set;
        for(auto & v:adjL[i]){
            set.emplace(v);
        }
        adjHL.emplace_back(set);
    }
    //cout<<"!!!***"<<endl;

    for(int j=0;j<adjR.size();j++){
        unordered_set<int> set;
        for(auto &u: adjR[j]){
            set.emplace(u);
        }
        adjHR.emplace_back(set);
    }


}

void bgraph::checkBi() {
    bool f = true;
    unordered_set<int> R;
    for(auto & ns: adjL){
        for(auto &v: ns){

            R.emplace(v);
        }
    }

    if(R.size()!=adjR.size()){
        cout<<"C L ne R Wrong"<<endl;
        f=false;
    }

    unordered_set<int> L;
    for(auto & ns: adjR){
        for(auto &u: ns){

            L.emplace(u);
        }
    }

    if(L.size()!=adjL.size()){
        cout<<"C R:"<<L.size()<<endl;
        cout<<"L:"<<adjL.size()<<endl;
        cout<<"C R ne L Wrong"<<endl;
        f=false;

    }

    if(f){

        cout<<"is a bipartite graph !"<<endl;
    }

}
void bgraph::readEdgeList(string file) {
    ifstream myfile (file);
    string line,u,v;
    //unordered_set<int> vsL,vsR;
    int numL = 0;
    int numR = 0;
    if (myfile.is_open()){
        while ( getline(myfile, line) ) {
            stringstream sline(line);
            getline(sline, u, ' ');
            getline(sline, v, ' ');
            int ul = stoi(u);
            int vr = stoi(v);
            if(numL == 0)
            {
                numL = ul;
                numR = vr;
                adjL.resize(numL);
                adjR.resize(numR);
                continue;
            }
            noE ++;
            adjL[ul].emplace_back(vr);
            adjR[vr].emplace_back(ul);
        }
    }
    novL = adjL.size();
    novR = adjR.size();
    max_degree_L = maxDegree(adjL);
    max_degree_R = maxDegree(adjR);
    myfile.close();
    // cout<<adjL.size()<<" "<<adjR.size()<<" "<<noE<<endl; 
    // cout<<"read done"<<endl;
}

void bgraph::switchLR() {
    adjL.swap(adjR);
    int tmp = max_degree_L;
    max_degree_L = max_degree_R;
    max_degree_R = tmp;
    isSwitch=!isSwitch;
}

int bgraph::maxDegree(vector<vector<int>> &adj) {
    int max=0;
    for(auto &neighbour_list:adj){
        if(neighbour_list.size()>max){
            max = neighbour_list.size();
        }
    }
    return max;
}

void bgraph::calculate_number_of_edges() {
    if(noE==0){
        if(novL<novR){
            int m=0;
            for(auto &neighbours:adjL){
                m=m+neighbours.size();
            }

            noE=m;
        }else{
            int m=0;
            for(auto& neighbours:adjR){
                m+=neighbours.size();
            }
            noE = m;
        }

    }
}
