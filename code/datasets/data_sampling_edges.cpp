#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <unistd.h>
#include <cstring>
#include <random>
std::vector<int> CreateRandomNums(int min,int max, int num, int seed)
{
	std::vector<int> res;
	res.clear();
	if (max - min + 1 < num)
	{
		return res;
	}
	srand(seed);
	for (auto i{0}; i < num; i++)
	{
		while (true)
		{
			auto temp{ rand() % (max + 1 - min) + min };
			auto iter{ find(res.begin(),res.end(),temp) };
			if (res.end() == iter)
			{
				res.push_back(temp);
				break;
			}
		}
	}
	return res;
}
double randf() {
     return (double)(rand()/(double)RAND_MAX);
}
std::vector<int>getLine(std::ifstream& fin, int &max_item)
{
    std::vector<int> list;
    char c;
    do 
    {
        int item = 0, pos = 0;
        c = fin.get();
        while((c>='0')&&(c<='9'))
        {
            item *= 10;
            item += int(c - '0');
            c = fin.get();
            pos++;
        }
        if(pos)
        {
            list.push_back(item);
            max_item = std::max(max_item, item);
        }
    }while(c != '\n' && !fin.eof());
    return list;
}
int main(int argc, char *argv[])
{
    int seed = 0;
    std::string infile = "", outfile = "";
    int opt;
    float sampling_rate = 0.2;
    while((opt = getopt(argc, argv, "i:o:s:r:"))!=-1)
    {
        switch(opt)
        {
        case 'i':
            infile = optarg;
            break;
        case 'o':
            outfile = optarg;
            break;
        case 's':
            seed = atoi(optarg);
            break;
        case 'r':
            sampling_rate = atof(optarg);
            break;
        }
    }
    srand(seed);
    if(infile == "")
    {
        std::cout << "Error: infile is empty" << std::endl;
        return 0;
    }
    if(outfile == "")
    {
        std::cout << "Error: outfile is empty" << std::endl;
        return 0;
    }
    std::ifstream fin(infile);
    std::ofstream fout(outfile);
    int max_l = -1, max_r = 0, edge_num = 0;
    std::vector<std::vector<int> > neighbor_list;
    while(!fin.eof())
    {
        neighbor_list.push_back(getLine(fin, max_r));
        max_l ++;
        edge_num += neighbor_list[max_l].size();
    }
    std::cout<<max_l<<" "<<max_r<<" "<<edge_num<<std::endl;
    std::vector<int> l_removed; 
    std::vector<int> counts_r(max_r+1, 0);
    std::vector<std::vector<int> >new_neighbor_list;
    int new_max_l = 0, new_max_r = 0, new_edge_num = 0;;
    for(int i=0; i<neighbor_list.size(); i++)
    {
        std::vector<int>temp;
        for(int j=0; j<neighbor_list[i].size(); j++)
        {
            if(randf()<sampling_rate)
            {
                temp.push_back(neighbor_list[i][j]);
            }
        }
        if(temp.size()>0)
        {
            new_edge_num += temp.size();
            new_neighbor_list.push_back(temp);
        }
    }
    for(int i=0; i<new_neighbor_list.size(); i++)
    {
        for(int pos = 0; pos<new_neighbor_list[i].size(); pos++)
        {
            counts_r[new_neighbor_list[i][pos]]++;
        }
    }
    std::vector<int> newPos_r(max_r + 1);
    for(int i=0, count=0; i<newPos_r.size(); i++)
    {
        newPos_r[i] = i - count;
        if(counts_r[i] == 0)
        {
            count++;
            newPos_r[i] = -1;
        }
        else
        {
            new_max_r = newPos_r[i];
        }
    }
    for(int i=0; i<new_neighbor_list.size(); i++)
    {
        new_max_l ++;
        bool flag = false;
        for(int pos = 0; pos<new_neighbor_list[i].size(); pos++)
        {
            int newPos = newPos_r[new_neighbor_list[i][pos]];
            if(newPos >= 0)
            {
                flag = true;
                fout<<newPos<<" ";
            }
        }
        if(flag)
        {
            fout<<"\n";
        }
    }
    std::cout<<new_max_l<<" "<<new_max_r<<" "<<new_edge_num<<std::endl;
    return 0;
}

