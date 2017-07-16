#! /bin/env bash

if [ $# != 1 ]
then
    echo "Usage: ./grab_ar.sh <beam>"
fi

psredit -c int:freq $1".ar" >> freq.txt
pdv -t -j "B 500" $1".ar" | awk 'NR!=1{print $2,$3,$4}' > $1".4ms.cube"
