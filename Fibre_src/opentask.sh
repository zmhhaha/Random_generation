#!/bin/sh
while ps -e | grep "$1" > /dev/null; do
  sleep 20
done
cd /home/node2/gitmaster/Random_generation/Fibre_src
cd $2
nohup ../bin/$1 1 0.6 100 100 100 0.2 19 > /dev/null 2>&1 &
let addnumber=$2+1
#echo $addnumber
cd /home/node2/gitmaster/Random_generation/Fibre_src
if ls | grep "$addnumber" > /dev/null; then
  nohup ./opentask.sh $1 $addnumber > /dev/null 2>&1 &
fi
