#!/bin/bash

#
# Compute RMSE
#

# arg1 train file
# arg2 predict file
# arg3 column to check against
echo "Computing rmse..."


# - 
# - 
# - 


tmp_ans=tmp_ans_`date +%H%M%S`.txt # just the correct answer
tmp_pred=tmp_pred_`date +%H%M%S`.txt # just the necessary predict column
tmp_ans_pred=tmp_ans_pred_`date +%H%M%S`.txt # correct answer and prediction together


#
awk '{print ('"$3"'==$1)?"1":"0"}' $1 > $tmp_ans

#
awk '{print $'"$3"'}' $2 > $tmp_pred

#
paste $tmp_ans $tmp_pred > $tmp_ans_pred


rm $tmp_ans $tmp_pred

#
awk -F'\t' 'BEGIN {n=0;sum=0} {n++; sum+=($1-$2)^2} END{print sqrt(sum/n)}' $tmp_ans_pred

rm $tmp_ans_pred

echo "Done"

