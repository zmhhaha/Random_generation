#!/bin/sh
local=$PWD
cd /mnt/sdb1
for i in {9,8};do
    for j in $(seq -f '%02g' 1 1 20);do
        c=0$i$j
        if ls | grep "$c" > /dev/null; then
            nohup $local/analyse.sh analyze $c 2 &> /dev/null &
        fi
    done
done