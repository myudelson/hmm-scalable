#!/bin/bash

# arg1 old submission file
# arg2 new prediction file (1st column is correct)

# # e.g.
# ./create_submit.sh algebra_2008_2009_submission.txt a89_kts_pred2.txt

echo "Creating submission file..."

tmp1=tmp1_`date +%H%M%S`.txt
tmp2=tmp2_`date +%H%M%S`.txt
tmp3=tmp3_`date +%H%M%S`.txt
tmp4=tmp4_`date +%H%M%S`.txt
tmp5=tmp5_`date +%H%M%S`.txt

mv $1 $tmp1 # save copy
sed 1d $tmp1 > $tmp2 # delete column names, now it's tmp2
rm $tmp1

awk -F'\t' 'BEGIN{OFS=""} {print $1}' $tmp2 > $tmp3 # just row numbers
rm $tmp2
awk -F'\t' 'BEGIN{OFS=""} {print $1}' $2 > $tmp4 # just corrects

# merge
paste $tmp3 $tmp4 > $tmp5 # under proper name
rm $tmp3 $tmp4

# add header
awk 'NR==1{print "Row	Correct First Attempt"}1' $tmp5 > $1
rm $tmp5

echo "Done"

# # After run this
# zip -b . kddcup bridge_to_algebra_2008_2009_submission.txt algebra_2008_2009_submission.txt
