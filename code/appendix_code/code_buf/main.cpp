#include <unistd.h>

#include <cstdlib>
#include <iostream>

#include "BicliqueFinder.h"
#include "BicliqueFinderFast.h"
#include "BitMBE.h"
#include "ComMBE.h"
#include "BoundMBE.h"
#include "PMBE.h"
#include "QuasiMBE.h"

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

void ExpFinderTest(char *graph_name, int finder_sel, int thread_num, char *fn, int lb=1, int rb=1, double miu = 1) {
  BiGraph *G = new BiGraph(graph_name);
  if (lb != 1 || rb != 1) {
    G->Prune1H(lb, rb);
    //G->Prune2H(lb, rb);
  }

  if (lb < rb ||  G->GetLSize() < G->GetRSize()) {
     G->Transpose();
     std::swap(lb, rb);
  }

  G->PrintProfile();
  BicliqueFinder *finder;
  switch (finder_sel) {
  case 0:
    finder = new MbeaFinder(G);
    break;
  case 1:
    finder = new ImbeaFinder(G);
    break;
  case 2:
    finder = new PmbeFinderV2(G);
    break;
  case 3:
    finder = new MineLMBCFinder(G);
    break;
  case 4:
    finder = new FmbeFinder(G);
    break;
  case 5:
    finder = new MmbeaInterAdvFinder(G);
    break;
  case 6:
    finder = new MmbeaIntraFinderFast(G);
    break;
  case 7:
    finder = new MmbeaFinderFast(G);
    break;
  case 8:
    finder = new ParMmbeaFinderFast(G, thread_num);
    break;
  case 9:
    finder = new MbeaFinderBitset(G);
    break;
  case 10:
    finder = new MbeaFinderAuto(G);
    break;
  case 11:
    finder = new ParMbeaFinderAuto(G, thread_num);
    break;
  case 12:
    finder = new MbeaAdvFinder(G);
    break;
  case 13:
    finder = new ImbeaAdvFinder(G);
    break;
  case 14:
    finder = new PmbeAdvFinder(G);
    break;
  case 15:
    finder = new LevelFinder(G);
    break;
  case 16:
    finder = new LevelFinderCache(G);
    break;
  case 17:
    finder = new BitMBEFinder(G);
    break;
  case 18:
    finder = new ComMBEFinder(G);
    break;
  case 19:
    finder = new BoundMBEFinder(G);
    break;
  case 20:
    finder = new ComMBEFinder2(G);
    break;
  case 21:
    finder = new ComMBEFinder3(G);
    break;
  case 22:
    finder = new AMBEAFinderNaive(G);
    break;
  case 23:
    finder = new AMBEAFinderInter(G);
    break;
  case 24:
    finder = new AMBEAFinderIntra(G);
    break;
  case 25:
    finder = new IntraFinder(G);
    break;
  case 26:
    finder = new AMBEAFinder(G);
    break;
    // case 19:
    //   QuasiMBEFinder* cur_finder = new QuasiMBEFinder(G);
    //   cur_finder->setMiu(miu);
    //   finder = cur_finder;
    //   break;
  }
  finder->Execute(lb, rb);
  if (fn != nullptr) {
    FILE *fp = fopen(fn, "a+");
    fprintf(fp, "Graph : %s\t Finder:%d\n", graph_name, finder_sel);
    fclose(fp);
  }
  PrintMemory(fn);
  finder->PrintResult(fn);
  delete finder;
  delete G;
}

int main(int argc, char *argv[]) {
  // for (int i = 1; i >= 0; i--) {
  //   printf("sel %d\n", i);
  //   BiGraph* G = new BiGraph("../datasets/db/Github.adj");
  //   MEBFinder* finder = new MEBFinder(G);
  //   finder->Execute(i);
  //   finder->PrintResult();
  //   delete finder;
  //   delete G;
  // }

  int opt;
  int thread_num = 1;
  double miu = 1;
  opterr = 0;
  char graph_name[80];
  char *output_fn = nullptr;
  int sel, lb=1, rb=1;
  while ((opt = getopt(argc, argv, "i:s:t:x:y:u:l")) != -1) {
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
    case 'u':
      miu = atof(optarg);
      break;
    }
  }
  ExpFinderTest(graph_name, sel, thread_num, output_fn, lb, rb, miu);
}
