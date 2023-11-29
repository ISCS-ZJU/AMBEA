#include <unistd.h>

#include <cstdlib>
#include <iostream>

#include "BaselineFinder.h"
#include "BicliqueFinder.h"
#include "AggressiveFinder.h"
// #include "BicliqueFinderFast.h"
// #include "BitMBE.h"
// #include "ComMBE.h"
// #include "BoundMBE.h"
// #include "PMBE.h"
// #include "QuasiMBE.h"

// void FinderTest(BicliqueFinder* finder = nullptr, char* fn = nullptr) {
//   if (finder == nullptr) return;
//   finder->Execute();
//   finder->PrintResult(fn);
//   delete finder;
// }

// void FinderTestBatch(BiGraph* G, char* fn = nullptr) {
//   FinderTest(new MbeaFinder(new BiGraph(*G)), fn);
//   FinderTest(new ImbeaFinder(new BiGraph(*G)), fn);
//   FinderTest(new MineLMBCFinder(new BiGraph(*G)), fn);
//   FinderTest(new FmbeFinder(new BiGraph(*G)), fn);
//   FinderTest(new MmbeaFinder(new BiGraph(*G)), fn);
//   FinderTest(new MmbeaInterFinder(new BiGraph(*G)), fn);
//   FinderTest(new MmbeaIntraFinder(new BiGraph(*G)), fn);
//   FinderTest(new PmbeFinder(new BiGraph(*G)), fn);
//   FinderTest(new MMbeaV2Finder(new BiGraph(*G)), fn);

// }

void PrintMemory(char *fn = nullptr) {
  FILE *fp = (fn == nullptr || strlen(fn) == 0) ? stdout : fopen(fn, "a+");
  int pid = getpid();
  unsigned int vmhwm = get_proc_vmhwm(pid);
  if (fn != nullptr)
    fseek(fp, 0, SEEK_END);
  fprintf(fp, "Memory Usage : %lfMB\n", vmhwm / 1000.0);
  if (fn != NULL)
    fclose(fp);
}

void ExpFinderTest(char *graph_name, int finder_sel, int thread_num, char *fn, int lb=1, int rb=1, double miu = 1, OrderEnum order = RInc) {
  BiGraph *G = new BiGraph(graph_name);
  if (lb != 1 || rb != 1) {
    G->Prune1H(lb, rb);
    //G->Prune2H(lb, rb);
  }

  if (lb < rb ||  G->GetLSize() < G->GetRSize()) {
     G->Transpose();
     std::swap(lb, rb);
  }
  AmbeaFinder *tmp_finder;
  G->PrintProfile();
  BicliqueFinder *finder;
  switch (finder_sel) {
    // baseline approach
    case 0:
      finder = new MineLMBCFinder(G);
      break;
    case 1:
      finder = new MbeaFinder(G);
      break;
    case 2:
      finder = new ImbeaFinder(G);
      break;
    case 3:
      finder = new FmbeFinder(G);
      break;
    case 4:
      finder = new PmbeFinder(G);
      break;

    // aggressive approach and variants
    case 5:
      tmp_finder = new AmbeaFinder(G);
      tmp_finder->SetOrder(order);
      finder = tmp_finder;
      break;
    case 6:
      finder = new ParAmbeaFinder(G, thread_num);
      break;
    case 7:
      finder = new AmbeaFinderNaive(G);
      break;
    case 8:
      finder = new AmbeaFinderInter(G);
      break;
    case 9:
      finder = new IntraFinder(G, true);
      break;
    case 10:
      finder = new IntraFinder(G, false, "IntraFinderPassive");
      break;
  }
  finder->Execute(lb, rb);
  if (fn != nullptr) {
    // FILE *fp = fopen(fn, "a+");
    // fprintf(fp, "Graph : %s\t Finder:%d\n", graph_name, finder_sel);
    // fclose(fp);
  }
  // PrintMemory(fn);
  finder->SetMemory(get_proc_vmhwm(getpid()) / 1000.0);
  finder->SetThreads(thread_num);
  finder->PrintResult(fn);
  delete finder;
  delete G;
}

int main(int argc, char *argv[]) {
  int opt;
  int thread_num = 1;
  double miu = 1;
  opterr = 0;
  char graph_name[80];
  char *output_fn = nullptr;
  int sel, lb=1, rb=1;
  OrderEnum order = RInc;
  while ((opt = getopt(argc, argv, "i:s:t:x:y:u:o:l")) != -1) {
    switch (opt) {
    case 'i':
      memcpy(graph_name, optarg, strlen(optarg) + 1);
      break;
    case 's':
      sel = atoi(optarg);
      break;
    case 'x':
      lb = atoi(optarg);
      break;
    case 'y':
      rb = atoi(optarg);
      break;
    case 'l':
      output_fn = new char[100];
      sprintf(output_fn, "%s-log.txt", graph_name);
      break;
    case 't':
      thread_num = atoi(optarg);
      break;
    case 'o':
      order = (OrderEnum) atoi(optarg);
      break;
    case 'u':
      miu = atof(optarg);
      break;
    }
  }
  ExpFinderTest(graph_name, sel, thread_num, output_fn, lb, rb, miu, order);
}
