#!/bin/bash

# Position arg 1 is txt file name in data folder
# Position arg 2 is number of alpha electrons
# Position arg 3 is number of beta electrons
echo "Position arg 1" $1
echo "Position arg 2" $2
echo "Position arg 3" $3

# Run the optimization  
./moleculeGeometryOptimizer data/$1 $2 $3 > data/$1.out

# Copy .xyz output to Jmol folder for animation
x=$1
y=${x%.txt}
z=${y##*/}.xyz
echo $z
cp ./data/$z /private/var/folders/4b/cqxrc9r92jb1jc2qq33tn9jr0000gn/T/hsperfdata_xperrylinn/$z
