#! /bin/env bash

# if [ $# != 1 ]
# then
#     echo "Usage: ./grab_ar.sh <beam>"
# fi

fl='2014-05-14-17:05:42_0000032532221952.000000_01.ar'
psredit $fl >> 1.head
psredit  -c int:freq $fl >> freq.txt
pdv -t $fl | awk 'NR!=1{print $2,$3,$4}' > "1.raw.cube"

