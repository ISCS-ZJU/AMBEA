#!/bin/bash

# "../datasets/db/ActorMovies.adj" 
# "../datasets/db/BookCrossing.adj" 
# "../datasets/db/DBLP.adj" 
# "../datasets/db/Github.adj" 
# "../datasets/db/IMDB.adj" 
# "../datasets/db/LiveJournal.adj" 
# "../datasets/db/StackOverflow.adj" 
# "../datasets/db/Teams.adj" 
# "../datasets/db/UCforum.adj" 
# "../datasets/db/Unicode.adj" 
# "../datasets/db/WebTrackers.adj" 
# "../datasets/db/Wikinews.adj" 
# "../datasets/db/Wikipedia.adj" 
# "../datasets/db/Writers.adj" 
# "../datasets/db/YouTube.adj" 
# ls ../datasets/db/movielens-10m_ti.adj

cmake .
make
# ./MBE_ALL -i ../datasets/db/Writers.adj -s 16

./MBE_ALL -i ../datasets/db/Writers.adj -s 5
# for graphfile in `ls ../datasets/db/Movi*.adj|grep -vE '/LiveJournal|/edit|/HHH|/T10'`
# do
# # graphfile="../datasets/db/YouTube.adj"
#   for i in 0 1 2 3 4 5 7 8 9 10
#   do
#      ./MBE_ALL -i $graphfile -s $i -l &
#      sleep 1
#   done
# done



# ./MBE_ALL -i ../datasets/db/Github.adj -s 19 -x 5 -y 5 -l&
# ./MBE_ALL -i ../datasets/db/Github.adj -s 19 -x 10 -y 5 -l&
# ./MBE_ALL -i ../datasets/db/Github.adj -s 19 -x 5 -y 10 -l&
# ./MBE_ALL -i ../datasets/db/Github.adj -s 19 -x 10 -y 10 -l&

# ./MBE_ALL -i ../datasets/db/Github.adj -s 18 -x 5 -y 5 -l&
# ./MBE_ALL -i ../datasets/db/Github.adj -s 18 -x 10 -y 5 -l&
# ./MBE_ALL -i ../datasets/db/Github.adj -s 18 -x 5 -y 10 -l&
# ./MBE_ALL -i ../datasets/db/Github.adj -s 18 -x 10 -y 10 -l&

# ./MBE_ALL -i ../datasets/db/Github.adj -s 13 -x 5 -y 5 -l&
# ./MBE_ALL -i ../datasets/db/Github.adj -s 13 -x 10 -y 5 -l&
# ./MBE_ALL -i ../datasets/db/Github.adj -s 13 -x 5 -y 10 -l&
# ./MBE_ALL -i ../datasets/db/Github.adj -s 13 -x 10 -y 10 -l&

# # ps aux|grep pz|grep MBE|awk '{print $2}'|xargs kill -9