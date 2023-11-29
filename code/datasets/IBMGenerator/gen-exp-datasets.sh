#!/bin/bash
make
# ./gen lit -nitems 10 -ntrans 10 -tlen 5 -ascii  -fname XD10I10T5
# ./gen lit -nitems 10 -ntrans 10 -tlen 10 -ascii  -fname XD10I10T10
# ./gen lit -nitems 10 -ntrans 10 -tlen 20 -ascii  -fname XD10I10T20
# ./gen lit -nitems 10 -ntrans 10 -tlen 40 -ascii  -fname XD10I10T40
# ./gen lit -nitems 10 -ntrans 10 -tlen 80 -ascii  -fname XD10I10T80

# ./gen lit -nitems 50 -ntrans 50 -tlen 5 -ascii  -fname XD50I50T5
# ./gen lit -nitems 100 -ntrans 100 -tlen 5 -ascii  -fname XD100I100T5
# ./gen lit -nitems 200 -ntrans 200 -tlen 5 -ascii  -fname XD200I200T5
# ./gen lit -nitems 500 -ntrans 500 -tlen 5 -ascii  -fname XD500I500T5
# ./gen lit -nitems 1000 -ntrans 1000 -tlen 5 -ascii  -fname XD1000I1000T5
# ./gen lit -nitems 2000 -ntrans 2000 -tlen 5 -ascii  -fname XD2000I2000T5
./gen lit -nitems 5 -ntrans 0.064 -tlen 512 -ascii  -fname HHH

mv *.data ../raw_dataset/
make clean
# 2W个item 1000W事务 20亿边 