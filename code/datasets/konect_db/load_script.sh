#!/bin/bash

for line in `cat url.txt`
do
	wget -c ${line}
done

for file in `ls *.bz2`
do
  tar -jxvf ${file}
done

for file in `find . -name "out.*"`
do
  mv ${file} ../raw_dataset/
done

# rm -rf !(url.txt|load_script.sh)
