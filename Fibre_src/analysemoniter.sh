#!/bin/sh
local=$PWD
#./analysemoniter.sh analyze 0601 5 50
while ps -e | grep "$1" > /dev/null; do
    sleep 20
done
cd /mnt/sdb1
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
#echo -e "nohup ../$1 FibreParameter.dat RawData.dat $c1.$c2 $c $3 $4 &> anarecord.log &"
nohup ../$1 FibreParameter.dat RawData.dat $c1.$c2 $c $3 $4 &> anarecord.log &
let addnumber=$c+10
if [ $addnumber -lt 10 ]; then
    addnumber=${c1}${c2}0${addnumber}
else
    addnumber=${c1}${c2}${addnumber}
fi
#echo -e "dirname=$addnumber"
cd /mnt/sdb1

if ls | grep "$addnumber" > /dev/null; then
    cd $local
    #echo -e "nohup ./analysemoniter.sh $1 $addnumber $3 $4 > /dev/null 2>&1 &"
    nohup ./analysemoniter.sh $1 $addnumber $3 $4 > /dev/null 2>&1 &
fi