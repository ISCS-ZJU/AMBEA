#!/bin/bash

g++ data_convert.cpp -O3 -o data_convert
# for file in `find . -name "out.*"`
# do
#   ./data_convert -i ${file}
# done

for file in `find . -name "*.data"`
do
  ./data_convert -i ${file}
done

rm data_convert
