#include "BicliqueFinder.h"

void MBC_star_test() { 
  double stime;
  BiGraph* G = new BiGraph();
  G->LoadGraphFromAdjListFile("mui_mmbea",3);
  //G->LoadGraphFromAdjListFile("lj_mmbea",3);
  //G->Prune(2,2);
  G->PrintProfile();
  MaximumEdgeBicliqueFinder finder(G);
  //stime = get_cur_time();
  //finder.Execute(1);
  //printf("processing time:%lf s\n", get_cur_time()-stime);
  //finder.Print();
  stime = get_cur_time();
  finder.Execute(1);
  printf("processing time:%lf s\n", get_cur_time() - stime);
  finder.Print();
}

void MBEA_test() {
  double stime, time_find;
  BiGraph* G = new BitBiGraph();
  BiGraph* mG = new BitBiGraph();
  PMbeaBiGraph* pG = new PMbeaBitBiGraph();
  G->LoadGraphFromAdjListFile("10_10_10.data", 3);
  //G->Prune(2,2);
  G->PrintProfile();
  mG->LoadGraphFromBiGraph(G);
  pG->LoadGraphFromBiGraph(G);
  mG->Reorder();
  
  BasicMbeaFinder finder(mG);
  stime = get_cur_time();
  finder.Execute();
  time_find = get_cur_time() - stime;
  finder.PrintResult(NULL, time_find);

  ImproveMbeaFinder finder2(mG);
  stime = get_cur_time();
  finder2.Execute();
  time_find = get_cur_time() - stime;
  finder2.PrintResult(NULL, time_find);

  MMbeaFinder finder3(mG);
  stime = get_cur_time();
  finder3.Execute();
  time_find = get_cur_time() - stime;
  finder3.PrintResult(NULL, time_find);

  PMbeaFinder finder4(pG);
  stime = get_cur_time();
  finder4.Execute();
  time_find = get_cur_time() - stime;
  finder4.PrintResult(NULL, time_find);

  MMbeaIntraFinder finder5(mG);
  stime = get_cur_time();
  finder5.Execute();
  time_find = get_cur_time() - stime;
  finder5.PrintResult(NULL, time_find);
}

int main() {
  //MBEA_test();
  //MBC_star_test();
  MBEA_test();
  return 0;
}
