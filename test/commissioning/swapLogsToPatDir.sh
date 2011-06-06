#!/bin/bash

for file in $(ls '/tmp/'`whoami`'/harvest_scripts/*.sh'); do
    echo $file
    sed 's/lxbatch_log/lxbatch_pat_log/' $file > /tmp/`whoami`/temp.sh
    mv /tmp/`whoami`/temp.sh $file
done