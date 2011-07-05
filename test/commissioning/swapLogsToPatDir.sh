#!/bin/bash

echo 'Changing log output directory to lxbatch_pat_log...'
for file in $(ls /tmp/`whoami`/harvest_scripts/analyze_PATTuple_*.sh); do
#    echo $file
    sed 's/lxbatch_log/lxbatch_pat_log/' $file > /tmp/`whoami`/temp.sh
    mv /tmp/`whoami`/temp.sh $file
done
echo "Done."