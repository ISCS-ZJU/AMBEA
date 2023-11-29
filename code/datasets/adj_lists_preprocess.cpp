#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unistd.h>

using namespace std;

int main(int argc, char *argv[]) {
  int opt;
  string infile = "", outfile = "out.dat";
  int start_id = 1;

  while ((opt = getopt(argc, argv, "i:o:b:")) != -1) {
    switch (opt) {
    case 'i':
      infile = optarg;
      break;
    case 'o':
      outfile = optarg;
      break;
    case 'b':
      start_id = atoi(optarg);
      break;

    default:
      break;
    }
  }
  std::ifstream fin(infile);
  std::ofstream fout(outfile);
  if (!fin.is_open() || !fout.is_open()) {
    std::cout << "Open file error!" << std::endl;
    exit(1);
  }
  std::map<int, int> item_map;
  std::string line;
  while (std::getline(fin, line)) {
    std::stringstream ss(line);
    int item;
    while (ss >> item) {
      item_map[item] = 0;
    }
  }
  fin.close();

  for (auto &p : item_map)
    p.second = start_id++;
  fin.open(infile);
  while (std::getline(fin, line)) {
    std::stringstream ss(line);
    int item;
    ss >> item;
    fout << item_map[item];
    while (ss >> item)
      fout << " " << item_map[item];
    fout << std::endl;
  }

  fin.close();
  fout.close();
}