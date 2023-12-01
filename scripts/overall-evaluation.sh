#! /bin/bash
if [ ! -d "./bin" ]
then
  mkdir ./bin
fi

if [ ! -f "./bin/MBE_ALL" ]
then
  cd ./code || exit
  cd MBE || exit
  mkdir build
  cd build || exit
  cmake ..
  make 
  mv MBE_ALL ../../../bin/
  cd ../../../
fi

if [ ! -f "./bin/mbbp" ]
then
  cd ./code || exit
  cd cohesive_subgraph_bipartite || exit
  mkdir build
  cd build || exit
  cmake ..
  make 
  mv mbbp ../../../bin/
  cd ../../../
fi

dataset_names=(Unicode UCforum Writers MovieLens Teams ActorMovies Wikipedia)
dataset_num=${#dataset_names[@]}


# figure 6: running on a machine with 96-core CPUs and a GPU
result_file="./scripts/results.txt"
progress_file="./scripts/progress.txt"
cur_time=$(date "+%Y-%m-%d %H:%M:%S")
echo $cur_time "Start overall evaluation." | tee -a $progress_file
rm $result_file
for ((i=0;i<dataset_num;i++)) 
do
  dataset_name=${dataset_names[i]}
  printf "%s \n" "$dataset_name" >> $result_file
  cur_time=$(date "+%Y-%m-%d %H:%M:%S")
  printf "%s %s %s %s %s\n" "Finder name," "Time(s)," "Memory(MB)," "Maximal Biclique," "Check Nodes" >> $result_file 
  echo $cur_time "Running MBEA on dataset" ${dataset_name} | tee -a $progress_file
  ./bin/MBE_ALL -i ./datasets/${dataset_name}.adj -s 1 >> ${result_file}  
  cur_time=$(date "+%Y-%m-%d %H:%M:%S")
  echo $cur_time "Running iMBEA on dataset" ${dataset_name} | tee -a $progress_file
  ./bin/MBE_ALL -i ./datasets/${dataset_name}.adj -s 2 >> ${result_file}  
  cur_time=$(date "+%Y-%m-%d %H:%M:%S")
  echo $cur_time "Running FMBE on dataset" ${dataset_name} | tee -a $progress_file
  ./bin/MBE_ALL -i ./datasets/${dataset_name}.adj -s 3 >> ${result_file}  
  cur_time=$(date "+%Y-%m-%d %H:%M:%S")
  echo $cur_time "Running PMBE on dataset" ${dataset_name} | tee -a $progress_file
  ./bin/MBE_ALL -i ./datasets/${dataset_name}.adj -s 4 >> ${result_file}  
  cur_time=$(date "+%Y-%m-%d %H:%M:%S")
  echo $cur_time "Running ooMBEA on dataset" ${dataset_name} | tee -a $progress_file
  ./bin/mbbp ./datasets/${dataset_name}.adj >> ${result_file} 
  cur_time=$(date "+%Y-%m-%d %H:%M:%S")
  echo $cur_time "Running AMBEA on dataset" ${dataset_name} | tee -a $progress_file
  ./bin/MBE_ALL -i ./datasets/${dataset_name}.adj -s 5 >> ${result_file}
  echo >> $result_file
done
