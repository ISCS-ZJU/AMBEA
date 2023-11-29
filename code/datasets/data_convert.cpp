#include <map>
#include <set>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <unistd.h>

void IBMGenFormatConvert(const std::string& infile, const std::string& outfile) {
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
    int row_id, row_id2, neighbors, item;
    ss >> row_id >> row_id2 >> neighbors;
    for (int i = 0; i < neighbors; i++) {
      ss >> item;
      item_map[item] = 0;
    }
  }
  fin.close();

  int new_id = 0;
  for (auto& p : item_map) p.second = new_id++;

  fin.open(infile);
  while (std::getline(fin, line)) {
    std::stringstream ss(line);
    int row_id, row_id2, neighbors, item;
    ss >> row_id >> row_id2 >> neighbors;
    if (neighbors == 0) continue;
    ss >> item;
    fout << item_map[item];
    for (int i = 1; i < neighbors; i++) {
      ss >> item;
      fout << " " << item_map[item];
    }
    fout << std::endl;
  }

  fin.close();
  fout.close();
}

void KonectFormatConvert(const std::string& infile,
                         const std::string& outfile) {
  std::ifstream fin(infile);
  std::ofstream fout(outfile);
  if (!fin.is_open() || !fout.is_open()) {
    std::cout << "Open file error!" << std::endl;
    exit(1);
  }

  std::map<int, std::set<int>> item_map;
  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty() || line[0] == '%') continue;
    std::stringstream ss(line);
    int src, dst;
    ss >> src >> dst;
    item_map[src].insert(dst-1); // adjust the start id from 1 to 0
  }
  fin.close();

  for (auto& p : item_map) {
    for (auto it = p.second.begin(); it != p.second.end(); it++) {
      if (it != p.second.begin()) fout << " ";
      fout << *it;
    }
    fout << std::endl;
  }
  fout.close();
}

int main(int argc, char *argv[]) {
  int opt;
  std::string infile, outfile;
  infile = "";
  const char *optstring = "i:";
  while ((opt = getopt(argc, argv, optstring)) != -1) {
    infile = optarg;
  }
  if(infile != ""){
    if(infile.find("out.") != -1){
      int pos = infile.find("out.");
      outfile = "db/" + infile.substr(pos + 4) + ".adj";
      std::cout << infile << " >> " << outfile << std::endl;
      KonectFormatConvert(infile, outfile);
    }
    else if(infile.find(".data") != -1){
      int pos = infile.find(".data");
      int start_pos = infile.rfind("/");
      outfile = "db/" + infile.substr(start_pos + 1, pos - start_pos - 1) + ".adj";
      std::cout << infile << " >> " << outfile << std::endl;
      IBMGenFormatConvert(infile, outfile);
    }
  }
  return 0;
}
