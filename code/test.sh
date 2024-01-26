#!/bin/sh
# 

./ae5.py 4PruebaTroph8 --numGen=2000 -p \
 --verbose --NumberOfCells=6000 \
 --NumberOfRsrcsInEachCell=20 \
 --species="0;NumberOfItems=32000;Distribution=100r;1;NumberOfItems=32000;Distribution=100r;2;Distribution=100r;NumberOfItems=20;3;Distribution=100r"
