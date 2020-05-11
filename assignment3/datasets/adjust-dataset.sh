#!/bin/tcsh
if (-f $4) then
  rm $4
endif

echo "Dataset $1 with new a_len = $2 and new b_len $3 stored at $4"
echo "$2 $3" > $4
awk 'NR==2' $1 | cut -c -$2 >> $4
awk 'NR==3' $1 | cut -c -$3 >> $4
