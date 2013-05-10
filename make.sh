#!/bin/bash

#
#
#
#
#


#
# KDD Cup 2010
#

# # Create data files
#
# ./create_train.sh bridge_to_algebra_2008_2009_train.txt bta89_ss_train.txt ss 0
# ./create_train.sh algebra_2008_2009_train.txt a89_ss_train.txt ss 0
# ./create_train.sh algebra_2008_2009_train.txt a89_kts_train.txt kts 0
# ./create_train.sh bridge_to_algebra_2008_2009_train.txt bta89_kts_train.txt kts 0
#
# ./create_train.sh algebra_2008_2009_train.txt a89_Uss_train.txt ss u
# ./create_train.sh algebra_2008_2009_train.txt a89_Ukts_train.txt kts u
# ./create_train.sh bridge_to_algebra_2008_2009_train.txt bta89_Uss_train.txt ss u
# ./create_train.sh bridge_to_algebra_2008_2009_train.txt bta89_Ukts_train.txt kts u
#
# ./create_test.sh bridge_to_algebra_2008_2009_test.txt bta89_ss_train.txt bta89_ss_test.txt ss 0
# ./create_test.sh bridge_to_algebra_2008_2009_test.txt bta89_kts_train.txt bta89_kts_test.txt kts 0
# ./create_test.sh algebra_2008_2009_test.txt a89_ss_train.txt a89_ss_test.txt ss 0
# ./create_test.sh algebra_2008_2009_test.txt a89_kts_train.txt a89_kts_test.txt kts 0
#
# ./create_test.sh algebra_2008_2009_test.txt a89_Uss_train.txt a89_Uss_test.txt ss u
# ./create_test.sh algebra_2008_2009_test.txt a89_Ukts_train.txt a89_Ukts_test.txt kts u
# ./create_test.sh bridge_to_algebra_2008_2009_test.txt bta89_Uss_train.txt bta89_Uss_test.txt ss u
# ./create_test.sh bridge_to_algebra_2008_2009_test.txt bta89_Ukts_train.txt bta89_Ukts_test.txt kts u
#
# # create a pure self-test file (same as train, but no observations)
# awk -F'\t' 'BEGIN{OFS=""} {print ".\t",$2,"\t",$3,"\t",$4}' a89_kts_train.txt > a89_kts_test0.txt
# awk -F'\t' 'BEGIN{OFS=""} {print ".\t",$2,"\t",$3,"\t",$4}' a89_ss_train.txt > a89_ss_test0.txt
#
# awk -F'\t' 'BEGIN{OFS=""} {print ".\t",$2,"\t",$3,"\t",$4}' a89_Ukts_train.txt > a89_Ukts_test0.txt
# awk -F'\t' 'BEGIN{OFS=""} {print ".\t",$2,"\t",$3,"\t",$4}' a89_Uss_train.txt > a89_Uss_test0.txt
#
# awk -F'\t' 'BEGIN{OFS=""} {print ".\t",$2,"\t",$3,"\t",$4}' bta89_kts_train.txt > bta89_kts_test0.txt
# awk -F'\t' 'BEGIN{OFS=""} {print ".\t",$2,"\t",$3,"\t",$4}' bta89_ss_train.txt > bta89_ss_test0.txt
#
# awk -F'\t' 'BEGIN{OFS=""} {print ".\t",$2,"\t",$3,"\t",$4}' bta89_Ukts_train.txt > bta89_Ukts_test0.txt
# awk -F'\t' 'BEGIN{OFS=""} {print ".\t",$2,"\t",$3,"\t",$4}' bta89_Uss_train.txt > bta89_Uss_test0.txt

# #
# # explore time gaps between studens' skill attempts in order to break sequences further
# #

# # explore when time $7 is empty and when time $8 is empty n=0 both full, n1 $7 is empty, n2 $8 is emp, n3 $7,$8 emp
# awk -F'\t' 'BEGIN{OFS=""; n0=0; n1=0; n2=0; n3=0} {n0+=($7!="" && $8!=""); n1+=($7=="" && $8!=""); n2+=($7!="" && $8==""); n3+=($7=="" && $8=="")} END{print "ok=",n0,"; t1=",n1,"; t2=",n2,"; both=", n3}' algebra_2008_2009_train.txt
# # ok=8652539; t1=265516; t2=0; both=0 :: only $7 is empty
# awk -F'\t' 'BEGIN{OFS=""} {print $2,"\t",$18,"\t",($7!="")?$7:$8}' algebra_2008_2009_train.txt > algebra_2008_2009_train_tmp.txt
# purge
# sed 1d algebra_2008_2009_train_tmp.txt > algebra_2008_2009_train_tmp2.txt
# rm algebra_2008_2009_train_tmp.txt
# sort algebra_2008_2009_train_tmp2.txt > algebra_2008_2009_train_SRT.txt
# rm algebra_2008_2009_train_tmp2.txt
# head -n1874 algebra_2008_2009_train_SRT.txt | tail -n2 # line 1874 is the first with all 3 measurements
# sed '1,1873d' algebra_2008_2009_train_SRT.txt > algebra_2008_2009_train_SRT2.txt # now clip first 1873 lines off
# rm algebra_2008_2009_train_SRT.txt
# awk -F'\t' 'BEGIN{OFS=""; t_m1=""} {if(t_m1!=""){print $3-t_m1;} t_m1=$3; }' algebra_2008_2009_train_SRT2.txt

# # Running EM, GD, EM+GD, GD+EM on KDD A89 and comparing residuals
# ./trainhmm -s 2 -0 0.5,1.0,0.4,0.8,0.2 -l 0,0,1,0,0,0,0,0,0,0 -u 1,1,1,0,1,1,1,0.3,0.3,1 -t 0.01 -q 1 -p 1 a89_kts_train.txt a89_kts_GDmodel.txt a89_kts_GDpred0.txt
# ./trainhmm -s 3 -0 0.5,1.0,0.4,0.8,0.2 -l 0,0,1,0,0,0,0,0,0,0 -u 1,1,1,0,1,1,1,0.3,0.3,1 -t 0.01 -q 1 -p 1 a89_kts_train.txt a89_kts_EMmodel.txt a89_kts_EMpred0.txt
# ./trainhmm -s 4 -0 0.5,1.0,0.4,0.8,0.2 -l 0,0,1,0,0,0,0,0,0,0 -u 1,1,1,0,1,1,1,0.3,0.3,1 -t 0.01 -q 1 -p 1 a89_kts_train.txt a89_kts_GDEMmodel.txt a89_kts_GDEMpred0.txt
# ./trainhmm -s 5 -0 0.5,1.0,0.4,0.8,0.2 -l 0,0,1,0,0,0,0,0,0,0 -u 1,1,1,0,1,1,1,0.3,0.3,1 -t 0.01 -q 1 -p 1 a89_kts_train.txt a89_kts_EMGDmodel.txt a89_kts_EMGDpred0.txt
# merge 1st columns of *pred0.txt files into one
# awk '{print $1}' a89_kts_GDpred0.txt > a1.txt
# awk '{print $1}' a89_kts_EMpred0.txt > a2.txt
# awk '{print $1}' a89_kts_GDEMpred0.txt > a3.txt
# awk '{print $1}' a89_kts_EMGDpred0.txt > a4.txt
# paste a1.txt a2.txt a3.txt a4.txt > a89_kts_GD,EM,GDEM,EMGDmodel.txt
# rm a1.txt a2.txt a3.txt a4.txt

# # Repeated skill labels for one step
# # For example
# grep 'Enter original median-1~~Enter original median-1' algebra_2008_2009_train.txt >o.txt
# # Find cases when KC's are listed next to each other
# grep -P '[\t|~]([A-Za-z][A-Za-z0-1 ,\-]+)~~\1[\t|~]' algebra_2008_2009_train.txt >o.txt
# perl -i.bak -pe 's/([~|\t])([A-Za-z][A-Za-z0-1 ,\-]+)~~\2/\1\2s/g' algebra_2008_2009_train.txt
# mv algebra_2008_2009_train.txt algebra_2008_2009_trainRM.txt
# mv algebra_2008_2009_train.txt.bak algebra_2008_2009_train.txt
# grep -P '[\t|~]([A-Za-z][A-Za-z0-1 ,\-]+)~~\1[\t|~]' algebra_2008_2009_trainRM.txt >o.txt
#
# grep -P '[\t|~]([A-Za-z][A-Za-z0-1 ,\-]+)~~\1[\n|~]' a89_kts_train.txt >o.txt
# perl -i.bak -pe 's/([~|\t])([A-Za-z][A-Za-z0-1 ,\-]+)~~\2/\1\2/g' a89_kts_train.txt
# mv a89_kts_train.txt a89_kts_trainRM.txt
# mv a89_kts_train.txt.bak a89_kts_train.txt
# 921 -> 866 KCs

#
# DataShop
#

# create list of input files 'ls ~/Documents/hcii/DataShop/data/data4hmm/*.txt > in.txt'

# # User-stratified cross-validation RMSE gradient descent
# for f in `ls ~/Documents/hcii/DataShop/data/data4hmm/*`; do
# 	echo "Processing $f"
# 	./trainhmm -s 2 -0 0.5,1.0,0.4,0.8,0.2 -l 0,0,1,0,0,0,0,0,0,0 -u 1,1,1,0,1,1,1,0.3,0.3,1 -t 0.01 -q 1 -v 10 "$f" m.txt
# 	echo " "
# done

# # Fitting AIC,BIC,RMSE (gradient descent)
# for f in `ls ~/Documents/hcii/DataShop/data/data4hmm/*`; do
# 	echo "Processing $f"
# 	./trainhmm -s 6 -0 0.5,1.0,0.4,0.8,0.2 -l 0,0,1,0,0,0,0,0,0,0 -u 1,1,1,0,1,1,1,0.3,0.3,1 -t 0.01 -q 1 -m 1 "$f" m.txt
# 	echo " "
# done
# run as './make.sh > out.txt'

#
# compare to CLI
#

# ./create_train.sh algebra_2008_2009_train.txt a89_kts_train.txt kts 0
# ./create_train.sh algebra_2008_2009_train.txt a89_ukts_train.txt kts u
# ./create_train.sh algebra_2008_2009_train.txt a89_uskts_train.txt kts s
# ./create_train.sh bridge_to_algebra_2008_2009_train.txt b89_kts_train.txt kts 0
# ./create_train.sh bridge_to_algebra_2008_2009_train.txt b89_ukts_train.txt kts u
# ./create_train.sh bridge_to_algebra_2008_2009_train.txt b89_uskts_train.txt kts s

# ./trainhmm -s 1.2 -m 1 -d ~ a89_kts_train.txt a89_kts_s1.2delim_model.txt > a89_kts_s1.2delim_out.txt
# ./trainhmm -s 4.2 -m 1 -d ~ a89_kts_train.txt a89_kts_s4.2delim_model.txt > a89_kts_s4.2delim_out.txt
# ./trainhmm -s 5.2 -m 1 -d ~ a89_kts_train.txt a89_kts_s5.2delim_model.txt > a89_kts_s5.2delim_out.txt
# ./trainhmm -s 6.2 -m 1 -d ~ a89_kts_train.txt a89_kts_s6.2delim_model.txt > a89_kts_s6.2delim_out.txt
#
# ./trainhmm -s 1.2 -m 1 -d ~ a89_ukts_train.txt a89_ukts_s1.2delim_model.txt > a89_ukts_s1.2delim_out.txt
# ./trainhmm -s 4.2 -m 1 -d ~ a89_ukts_train.txt a89_ukts_s4.2delim_model.txt > a89_ukts_s4.2delim_out.txt
# ./trainhmm -s 5.2 -m 1 -d ~ a89_ukts_train.txt a89_ukts_s5.2delim_model.txt > a89_ukts_s5.2delim_out.txt
# ./trainhmm -s 6.2 -m 1 -d ~ a89_ukts_train.txt a89_ukts_s6.2delim_model.txt > a89_ukts_s6.2delim_out.txt
#
# ./trainhmm -s 1.2 -m 1 -d ~ a89_uskts_train.txt a89_uskts_s1.2delim_model.txt > a89_uskts_s1.2delim_out.txt
# ./trainhmm -s 4.2 -m 1 -d ~ a89_uskts_train.txt a89_uskts_s4.2delim_model.txt > a89_uskts_s4.2delim_out.txt
# ./trainhmm -s 5.2 -m 1 -d ~ a89_uskts_train.txt a89_uskts_s5.2delim_model.txt > a89_uskts_s5.2delim_out.txt
# ./trainhmm -s 6.2 -m 1 -d ~ a89_uskts_train.txt a89_uskts_s6.2delim_model.txt > a89_uskts_s6.2delim_out.txt
#
#
# ./trainhmm -s 1.2 -m 1 -d ~ b89_kts_train.txt b89_kts_s1.2delim_model.txt > b89_kts_s1.2delim_out.txt
# ./trainhmm -s 4.2 -m 1 -d ~ b89_kts_train.txt b89_kts_s4.2delim_model.txt > b89_kts_s4.2delim_out.txt
# ./trainhmm -s 5.2 -m 1 -d ~ b89_kts_train.txt b89_kts_s5.2delim_model.txt > b89_kts_s5.2delim_out.txt
# ./trainhmm -s 6.2 -m 1 -d ~ b89_kts_train.txt b89_kts_s6.2delim_model.txt > b89_kts_s6.2delim_out.txt
#
# ./trainhmm -s 1.2 -m 1 -d ~ b89_ukts_train.txt b89_ukts_s1.2delim_model.txt > b89_ukts_s1.2delim_out.txt
# ./trainhmm -s 4.2 -m 1 -d ~ b89_ukts_train.txt b89_ukts_s4.2delim_model.txt > b89_ukts_s4.2delim_out.txt
# ./trainhmm -s 5.2 -m 1 -d ~ b89_ukts_train.txt b89_ukts_s5.2delim_model.txt > b89_ukts_s5.2delim_out.txt
# ./trainhmm -s 6.2 -m 1 -d ~ b89_ukts_train.txt b89_ukts_s6.2delim_model.txt > b89_ukts_s6.2delim_out.txt
#
# ./trainhmm -s 1.2 -m 1 -d ~ b89_uskts_train.txt b89_uskts_s1.2delim_model.txt > b89_uskts_s1.2delim_out.txt
# ./trainhmm -s 4.2 -m 1 -d ~ b89_uskts_train.txt b89_uskts_s4.2delim_model.txt > b89_uskts_s4.2delim_out.txt
# ./trainhmm -s 5.2 -m 1 -d ~ b89_uskts_train.txt b89_uskts_s5.2delim_model.txt > b89_uskts_s5.2delim_out.txt
# ./trainhmm -s 6.2 -m 1 -d ~ b89_uskts_train.txt b89_uskts_s6.2delim_model.txt > b89_uskts_s6.2delim_out.txt

#
# unique Unit-Section-problem
#
awk -F'\t' 'BEGIN{OFS=""} {print $3}' a89_kts_train.txt > a89_kts_train__ITEM.txt # save just item
sort a89_kts_train__ITEM.txt > a89_kts_train__ITEM_srt.txt

