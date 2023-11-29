

#ifndef MBBP_BGRAPH_H
#define MBBP_BGRAPH_H
#include<vector>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
//#include <boost/unordered_set.hpp>
//#include <boost/unordered_map.hpp>
using namespace std;
class bgraph{
public:
    int des_count=0;
    int novL;
    int novR;
    int noE;
    int max_degree_L;
    int max_degree_R;
    //int *test;
    int offset;
    bool isSwitch=false;
    //adjList left;
    vector<vector<int>> adjL;
    vector<unordered_set<int>> adjHL;
    //adjList right;
    vector<vector<int>> adjR;
    vector<unordered_set<int>> adjHR;




    bgraph();
    void iniHash();
    void readAdj(int,string); //0 read ajdL, 1 read adjR
    void readAdj(string);
    void dispalyAdj();
    void sortAdj();
    void checkAdj();
    int edgeNumber();
    void writAdj(int, string);
    void checkBi();
    void readEdgeList(string);
    void switchLR();
    int maxDegree(vector<vector<int>>&);
    void calculate_number_of_edges();
    std::vector<int> CommonNeighborsL(const std::vector<int>&);
    std::vector<int> CommonNeighborsR(const std::vector<int>&);

};


struct Edge{
    int l;
    int r;

};
/*
 * open question: do we need compute A\cup B\cup CL\cup CR induced subgraph
 * If this induced subgraph is used just once then we do not
 * If this induced subgraph is used multiple times O(n)-times then we build it.
 * */


#endif //MBBP_BGRAPH_H

