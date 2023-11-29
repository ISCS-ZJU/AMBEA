#pragma once
#include <unistd.h>
#include <cstring>
#include <stdio.h>
#include <vector>

#define VMHWM_LINE 21
unsigned int get_proc_vmhwm(unsigned int pid);
void PrintMemory(char* fn = nullptr);
std::vector<int> seq_intersect(const std::vector<int>& v0,
  const std::vector<int>& v1);
