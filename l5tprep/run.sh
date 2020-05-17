#!/bin/sh

NB_PROC=${1:-4}
TYPE=${2:-2}
DIM_A=${3:-8}
NB_COL_B=${4:-4}
clear

make --silent clean
make --silent DBG_ON_FD=1

OUTPUT=$(mpirun --oversubscribe -n $NB_PROC ./*.run $TYPE $DIM_A $NB_COL_B)

for i in *.log;
do 
    FILE="$i"
    echo "$FILE output"
    cat "./$FILE"
    echo
    echo
done

rm *.log

echo $OUTPUT

