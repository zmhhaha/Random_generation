#!/bin/sh
while ps -e | grep "$1" > /dev/null; do
    sleep 20
done
cd ~/gitmaster/Random_generation/Fibre_src
cd $2
c1=`echo $2 |cut -c 1`;#echo -e "c1=$c1\t"
c2=`echo $2 |cut -c 2`;#echo -e "c2=$c2\t"
c3=`echo $2 |cut -c 3`;#echo -e "c3=$c3\t"
c4=`echo $2 |cut -c 4`;#echo -e "c4=$c4\t"
if [ $c3 -eq 0 ]; then
    c=$c4
else
    c=${c3}${c4}
fi
#echo -e "c=$c"
nohup ../bin/$1 1 $c1.$c2 100 100 100 0.2 $c &> $2.out &
let addnumber=$c+10
if [ $addnumber -lt 10 ]; then
    addnumber=${c1}${c2}0${addnumber}
else
    addnumber=${c1}${c2}${addnumber}
fi
#echo -e "dirname=$addnumber"
cd ~/gitmaster/Random_generation/Fibre_src
if ls | grep "$addnumber" > /dev/null; then
    nohup ./opentask2.sh $1 $addnumber > /dev/null 2>&1 &
fi
