#!/bin/sh
local=$PWD
aimdir=/mnt/sdb1
cd $aimdir
for i in {6,7,8,9};do
    #for j in $(seq -f '%02g' 1 1 20);do
    for j in {01,10,20};do
        
        c=0$i$j
        if ls | grep "$c" > /dev/null; then
            #echo $c
            #cd $c
            ##sleep 60
            ##while ps -e |grep analyze &>/dev/null; do
            ##    sleep 60
            ##done
            #cd dfv
            #k=pathall.dat
            #for k in $(ls |grep path);do
            #    echo $k
            #    sleep 60
            #    while ps -e |grep analyze &>/dev/null; do
            #        sleep 60
            #    done
            #    nohup $local/analyse.sh analyze $c 6 $k &> /dev/null &
            #done
            nohup $local/analyse.sh analyze $c 1 1 &> /dev/null &
        fi
        cd $aimdir
    done
    #done
done