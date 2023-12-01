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

