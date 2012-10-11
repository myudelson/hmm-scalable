#!/bin/bash

# arg1 source test file
# arg2 created train file
# arg3 target train file
# arg4 if 0 - ss, if 1 - kts
# arg5 if 0  - default, 1 - section
echo "Creating test file..."

tmp4=tmp4_`date +%H%M%S`.txt

# first perform same thing as with training
./create_train.sh $1 $tmp4 $4 $5 t

# then concatenate to train
cat $2 $tmp4 > $3

rm $tmp4

echo "Done"

