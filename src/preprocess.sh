#!/bin/bash

for i in {0..1023}
do
    j=`printf "%03x" $i`
    echo "sorting store/$j"
    sort store/$j -o store/$j
done
