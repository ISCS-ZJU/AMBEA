#include "Utility.h"
unsigned int get_proc_vmhwm(unsigned int pid){

	char file_name[64]={0};
	FILE *fd;
	char line_buff[512]={0};
	sprintf(file_name,"/proc/%d/status",pid);

	fd =fopen(file_name,"r");
	if(nullptr == fd){
		return 0;
	}

	char name[64];
	int vmhwm;
	for (int i=0; i<VMHWM_LINE-1;i++){
		fgets(line_buff,sizeof(line_buff),fd);
	}

    while(true)
    {
	    fgets(line_buff,sizeof(line_buff),fd);
	    sscanf(line_buff,"%s %d",name,&vmhwm);
        if(!strcmp(name, "VmHWM:"))
        {
            break;
        }
    }
    
    fclose(fd);

	return vmhwm;
}
void PrintMemory(char* fn) {
  FILE* fp = (fn == nullptr || strlen(fn) == 0) ? stdout : fopen(fn, "a+");
  int pid = getpid();
  unsigned int vmhwm = get_proc_vmhwm(pid);
  if (fn != nullptr) fseek(fp, 0, SEEK_END);
  fprintf(fp, "Memory Usage : %lfMB\n", vmhwm / 1000.0);
  if (fn != NULL) fclose(fp);
}
std::vector<int> seq_intersect(const std::vector<int>& v0,
                               const std::vector<int>& v1) {
  std::vector<int> res;
  for (auto it0 = v0.begin(), it1 = v1.begin();
       it0 != v0.end() && it1 != v1.end();) {
    if (*it0 == *it1) {
      res.emplace_back(*it0);
      it0++;
      it1++;
    } else if (*it0 > *it1)
      it1++;
    else
      it0++;
  }
  return res;
}


