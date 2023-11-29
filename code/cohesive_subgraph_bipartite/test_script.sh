rm -rf build
mkdir build
cd build 
cmake ..
make 


./mbbp ../../datasets/db/DBLP.adj >> ../../datasets/db/DBLP.adj-log.txt &
./mbbp ../../datasets/db/ActorMovies.adj >> ../../datasets/db/ActorMovies.adj-log.txt &
./mbbp ../../datasets/db/Wikipedia.adj >> ../../datasets/db/Wikipedia.adj-log.txt &
./mbbp ../../datasets/db/Teams.adj >> ../../datasets/db/Teams.adj-log.txt &
./mbbp ../../datasets/db/YouTube.adj >> ../../datasets/db/YouTube.adj-log.txt &
./mbbp ../../datasets/db/DBLP.adj 1 >> ../../datasets/db/DBLP.adj-log.txt &
./mbbp ../../datasets/db/ActorMovies.adj 1 >> ../../datasets/db/ActorMovies.adj-log.txt &
./mbbp ../../datasets/db/Wikipedia.adj 1 >> ../../datasets/db/Wikipedia.adj-log.txt &
./mbbp ../../datasets/db/Teams.adj 1 >> ../../datasets/db/Teams.adj-log.txt &
./mbbp ../../datasets/db/YouTube.adj 1 >> ../../datasets/db/YouTube.adj-log.txt &
./mbbp ../../datasets/db/DBLP.adj 2 >> ../../datasets/db/DBLP.adj-log.txt &
./mbbp ../../datasets/db/ActorMovies.adj 2 >> ../../datasets/db/ActorMovies.adj-log.txt &
./mbbp ../../datasets/db/Wikipedia.adj 2 >> ../../datasets/db/Wikipedia.adj-log.txt &
./mbbp ../../datasets/db/Teams.adj 2 >> ../../datasets/db/Teams.adj-log.txt &
./mbbp ../../datasets/db/YouTube.adj 2 >> ../../datasets/db/YouTube.adj-log.txt &

# for graphfile in `ls ../../datasets/db/*.adj`
# do 
#   ./mbbp $graphfile >> "$graphfile-log.txt" &
#   sleep 1
# done