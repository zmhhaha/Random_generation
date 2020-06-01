#!/bin/sh
#./test.sh analyze 0601 1 1 df
#./test.sh analyze 0601 6 pcdata.dat df
dir=/mnt/sdb1
direction=$5
cd $dir
cd opencv
if [ ! -d $direction ];then
    mkdir $direction
fi
cd $direction
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
nohup ../../$1 ../FibreParameter.dat ../Meananalyse.dat $c1.$c2 $c $3 $4 &> anarecord.log &
