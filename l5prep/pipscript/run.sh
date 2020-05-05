#!/bin/sh

NB_PROC=$1

clear

TMP=$(expr $NB_PROC - 1)

make --silent clean
make --silent

mpirun -n $NB_PROC ./*.run

for i in `seq 0 $TMP`
do 
    FILE="f$i"
    echo "$FILE output"
    cat "./$FILE"
    echo
    echo
done

for i in `seq 0 $TMP`
do 
    FILE="f$i"
    rm "./$FILE" 
done

make --silent clean