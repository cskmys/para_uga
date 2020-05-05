#!/bin/sh

NB_PROC=${1:-8}
DIM_A=${2:-8}
NB_COL_B=${3:-4}
clear

make --silent clean
make --silent DBG_ON_FD=1

OUTPUT=$(mpirun --oversubscribe -n $NB_PROC ./*.run $DIM_A $NB_COL_B)

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

