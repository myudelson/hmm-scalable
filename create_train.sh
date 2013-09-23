#!/bin/bash

# arg1 source file
# arg2 target file
# arg3 ss, kts
# arg4 if 0  - default, s - add unit and section to kc name '__' is separator lowercase, u - just unit
# arg5 0 default, anything else markind test file
# arg6 0 default, anything else add 5th column with POSIX time (seconds since 1970)
echo "Creating train file..."

tmp1=tmp1_`date +%H%M%S`.txt
tmp2=tmp2_`date +%H%M%S`.txt

# Columns in source file
#  2 - student id
#  3 - problem hierarchy - unit & section
#  4 - problem name
#  6 - step name
#  8 - first transaction time
# 14 - first time correct
# 20 - KTS skill
# 18 - SS skill

first_column='2-$14'
# if [ -n $5 ] || [ $5 != "0" ] # if not null and not 0, then it's test
if [ "$5" != "0" ]
then
	first_column='"."'
fi

time_awk1a=''
time_awk1b=''
time_awk2=''
#if [ -n $6 ] || [ $6!= "0" ] # if not null and not 0, then add time
if [ "$6" != "0" ]
then
	 time_awk1a='gsub(/:/," ",$8);gsub(/-/," ",$8);t1=mktime($8);'
	 time_awk1b=',"\t",t1'
	 time_awk2=',"\t",$5'
fi


first_awk=''
# construct first awk
if [ "$3" = "ss" ]
then #ss
	first_awk='$18'
elif [ "$3" = "kts" ]
then #kts
	first_awk='$20'
fi

middle_awk1=''
middle_awk2=''
second_awk=''
if [ "$4" = "0" ]
then #default
	second_awk='((lengt h($4)==0)?".":$4)' # if the KC field is empty use '.'
	# no unit
	middle_awk1=''
	middle_awk2='$3'
	first_awk='tolower('$first_awk')'
elif [ "$4" = "s" ]
then # Unit and Section
	second_awk='((substr($4,length($4),1)=="_")?".":$4)'
	# add prefix for column
# 	middle_awk1='split($3,a,", Section"); '
# 	middle_awk2='a[1]'
# 	first_awk='a[1],"__",'$first_awk
	middle_awk1='sub(/, Section /, "__", $3); sub(/Unit /, "", $3); gsub(/~~/,"~~"$3"__",'$first_awk'); '
	middle_awk2='$3' # sub(/, Section/, "__", $3)
	first_awk='tolower($3),"__",tolower('$first_awk')'
elif [ "$4" = "u" ]
then # Unit only
	second_awk='((substr($4,length($4),1)=="_")?".":$4)'
	# add prefix for column
# 	middle_awk1='split($3,a,", Section "); '
# 	middle_awk2='a[1]'
# 	first_awk='a[1],"__",'$first_awk
	middle_awk1='sub(/, Section /, "__", $3); sub(/Unit /, "", $3);  split($3,a,"__"); gsub(/~~/,"~~"a[1]"__",'$first_awk'); '
	middle_awk2='$3'
	first_awk='tolower(a[1]),"__",tolower('$first_awk')'
fi

# echo $first_column
# echo $first_awk
# echo $second_awk

# awk -F'\t' 'BEGIN{OFS=""} {'"$middle_awk1"'print '"$first_column"',"\t",$2,"\t",'"$middle_awk2"',"__",$4,"__",$6,"\t",'"$first_awk"'}' $1 > $tmp1
# removed step ($6) to reduce number of items

gawk -F'\t' 'BEGIN{OFS=""} {'"$middle_awk1"''"$time_awk1a"'print '"$first_column"',"\t",$2,"\t",'"$middle_awk2"',"__",$4,"\t",'"$first_awk"''"$time_awk1b"'}' $1 > $tmp1
sed 1d $tmp1 > $tmp2
rm $tmp1
awk -F'\t' 'BEGIN{OFS=""} {print $1,"\t",$2,"\t",$3,"\t",'"$second_awk"''"$time_awk2"'}' $tmp2 > $2

rm $tmp2

echo "Done"