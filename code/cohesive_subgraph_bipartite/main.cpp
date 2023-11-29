#include <iostream>
#include "bgraph.h"
#include "basicEnu.h"
#include "LCUtility.h"
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <unistd.h>
//#include "Graph.h"
#include <cstring>
#include "Utility.h"
#include "UniCore.h"
#include "basicEnuSet.h"
#include "denseMBB.h"
#include "BiCore.h"
#include "ooMBE.h"
using namespace std;

void testDataDense(bgraph &g);
void testDataSparse(bgraph &g);

double GetMemory(){
  int pid = getpid();
  unsigned int vmhwm = get_proc_vmhwm(pid);
  return vmhwm / 1000.0;
}



int main(int argc, char * argv[]) {


    //newEdgeList("../out.livejournal-groupmemberships","",0);
    //edgeFilter("../orkut.edges","../borkut.edge",3072441);


    //vector<int> test;
    //test.reserve(10);
    //test.emplace_back(1);
    //cleanUpLJ();
    //newEdgelist();
    //cleanUp();
    //eLtoAdjL();
   // Graph *graph = new Graph("..");
    //graph->core_decomposition(1, true);
    //mbsToMBB();

    //adjMtoAdj();

    //bg.readAdj(0,"../kljR.txt");
    //bg.readAdj(1,"../kljL.txt");
    //bg.readAdj(0,"../466lork.txt");
    //bg.readAdj(1,"../466rork.txt");

    bgraph bg;
    //bg.readAdj(1,"../rmarveladj.txt"); //no_reuse 206135
    //bg.readAdj(0,"../lmarveladj.txt");//reuse 206135
    /**
     * BX for bookcrossing
     * IM2 for IMDB large
     * LG for livejournal
     */
    string filename = argv[1];
    // cout<<filename<<endl;
    if(filename.find(".graph") != filename.npos)
    {
        bg.readEdgeList(filename);
    }
    else
    {
        bg.readAdj(filename);
    }
    bg.sortAdj();
    //bg.dispalyAdj();
    //bg.readAdj(1,"../testRadj.txt");

    //bg.switchLR();
    //ooMBE alg1(bg);
    //alg1.mbeStart();
    //cout<<"!!!!!!!! 1:"<<alg1.nomb<<endl;
    //BiCore b(bg);
    //b.e2hopCore(bg.adjR,bg.adjL);
    //cout<<b.maxE2Core<<endl;

    //BiCore b2(bg);
    //b2.e2hopCore(bg.adjL,bg.adjR);
    //cout<<b2.maxE2Core<<endl;
    ooMBE alg2(bg);
    //alg2.advIMBEA();
    int order = 0;
    if (argc > 2) order = atoi(argv[2]);
    alg2.adv_mbeStart_reuse_full(order);
    //alg2.mbeStart();
    //alg2.adv_mbeStart_reuse();
    //alg2.adv_mbeStart();

    double memory = GetMemory();
    printf("ooMBEAFinder, %lf, %lf, %d, %d\n", alg2.exe_time_, memory,
           alg2.nomb, alg2.processing_nodes
           );
    //printf("ooMBEAFinder, %lf, %lf, %d, %d, %u, %lf\n", alg2.exe_time_, memory,
    //       alg2.nomb, alg2.processing_nodes, alg2.max_level_,
    //       (double)alg2.level_accumulation_ / alg2.nomb);

    // PrintMemory();
    // cout<<"processing_nodes: "<<alg2.processing_nodes<<endl;
    // cout<<"Maximal nodes: "<<alg2.nomb<<endl;
    // cout<<"Maximum level: "<<alg2.max_level_<<endl;
    // cout<<"Average level: "<<(double)alg2.level_accumulation_ / alg2.nomb <<endl;

    //bg.readEdgeList("../LG");
    //bg.readEdgeList("../MC");
    //cout<<bg.adjL.size()<<endl;
    //cout<<bg.adjR.size()<<endl;
    //bg.readAdj(0,"../lorkutadj.txt");
    //bg.readAdj(1,"../rorkutadj.txt");
    //testDataSparse(bg);
    //testDataDense(bg);
    //bg.readAdj(1,"../coreLLj.txt");
    //bg.readAdj(0,"../coreRLj.txt");
    //bg.iniHash();
    //denseMBB alg(bg);
    //alg.startEnu();
    //BiCore bc(bg);
    //bc.e2hopCore(bg.adjL, bg.adjR);
    //bc.e2hopCore(bg.adjR, bg.adjL);

    //bc.e2hopCore(bg.adjR, bg.adjL);

    //bc.bicoreDec();

    //bc.showDegeneracyOrder();
    //bc.bicoreDec();
    //bc.write2Hop();
    //bg.offset = 3201203;
    //bg.readAdj("../");
    //bg.readAdj("../dataset_1_12/adjLJ.txt");
    //bg.readAdj(0,"../466lork.txt");
    //bg.readAdj(1,"../466rork.txt");
    //bg.readAdj(0,"../lorkutadj.txt");
    //bg.readAdj(1,"../rorkutadj.txt");
    //bg.checkBi();
    //bg.iniHash();
    //testDataSparse(bg);
    //bg.iniHash();
    //int noe = bg.edgeNumber();
    //cout<<"|L|"<<bg.adjL.size()<<endl;
    //cout<<"|R|"<<bg.adjR.size()<<endl;
    //cout<<"|W|"<<noe<<endl;

    //testDataDense(bg);
    //bg.checkAdj();
    //UniCore core(bg.adjL,bg.adjR);
    //cout<<"max core:"<<core.getMaxCore();
    //vector<int> L;
    //vector<int> R;
    //core.getVertexSetCore(L,R,core.getMaxCore());
    //cout<<L.size()<<endl;
    //cout<<R.size()<<endl;
    //bgraph ibg;
    //computeIG(L,R,bg,ibg);
    //ibg.writAdj(0,"../kcoreLLj.txt");
    //ibg.writAdj(1,"../kcoreRLj.txt");
    //basicEnuSet alg(bg);
    //denseMBB alg(bg);
    //alg.maxsize=500;
   // alg.indstart();
    //ibg.dispalyAdj();
    //cout<<bg.novL<<" "<<bg.novR<<endl;
    //bg.dispalyAdj();
    //bg.novL=3;
    //bg.novR=2;
    //
    //testDataSparse(bg);
    //bg.dispalyAdj();
    //bg.sortAdj();
    //bg.dispalyAdj();

    //basicEnuSet alg(bg);
    //alg.maxSize = 15;

    return 0;





    //bg.dispalyAdj();
}

void testDataDense(bgraph &bg){
    //ini leftadj
    bg.adjL.push_back({0,1,2,3});
    bg.adjL.push_back({0,1,4,2});
    bg.adjL.push_back({1,2,3,4});
    bg.adjL.push_back({0,2,3,4});
    bg.adjL.push_back({0,1,3,4});


    //bg.adjL.push_back({0});
    //bg.adjL.push_back({0,1});
    //bg.adjL.push_back({0});


    //ini rightadj
    bg.adjR.push_back({0,1,3,4});
    bg.adjR.push_back({0,1,2,4});
    bg.adjR.push_back({0,1,2,3});
    bg.adjR.push_back({0,2,3,4});
    bg.adjR.push_back({1,2,3,4});


   // bg.adjR.push_back({0,1,2});
    //bg.adjR.push_back({1});
}


void testDataSparse(bgraph &bg){
    //ini leftadj
    bg.adjL.push_back({0,6});
    bg.adjL.push_back({0,1});
    bg.adjL.push_back({1,2,3});
    bg.adjL.push_back({2,3,4});
    bg.adjL.push_back({2,3, 5});
    bg.adjL.push_back({1});


    //ini rightadj
    bg.adjR.push_back({0,1});
    bg.adjR.push_back({1,2,5});
    bg.adjR.push_back({2,3,4});
    bg.adjR.push_back({2,3,4});
    bg.adjR.push_back({3});
    bg.adjR.push_back({4});
    bg.adjR.push_back({0});
}

/**
 *
 *  string line;
    ifstream myfile ("../orkut.edges");
    vector<Edge> edgeList;
    string ls, rs;
    if (myfile.is_open())
    {
        while ( getline(myfile, line) ){

            istringstream sline( line );
            getline(sline,ls,',');
            getline(sline,rs,',');

            int l = stoi (ls);
            int r = stoi (rs);
            struct Edge e;
            e.l = l;
            e.r = r;
            if(e.r>=3072442){
                edgeList.push_back (e);
            }

        }
        cout<<edgeList.size() << endl;
        myfile.close();


    }
    else
    {
        cout << "Unable to open file";
    }


    ofstream fout;
    fout.open("../eedges_orkut.txt");
    if (fout.is_open()){
        for(auto& e: edgeList){
            fout<<e.l<<","<<e.r<<"\n";
        }
        fout.close();
    }

*/


/**
 *
 * 5
0,4,0,1,2,3
1,4,0,1,2,4
2,4,1,2,3,4
3,4,0,2,3,4
4,4,0,1,3,4
5
0,4,0,1,3,4
1,4,0,1,2,4
2,4,0,1,2,3
3,4,0,2,3,4
4,4,1,2,3,4
 *
 *
 * */
