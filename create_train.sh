#!/bin/bash

# arg1 source file
# arg2 target file
# arg3 ss, kts
# arg4 if 0  - default, u - unit
# arg5 anything, markind test file
echo "Creating train file..."

tmp1=tmp1_`date +%H%M%S`.txt
tmp2=tmp2_`date +%H%M%S`.txt

first_column='2-$14'
if [ -n $5 ] # if exists, then it's test
then
	first_column='"."'
fi


first_awk=''
# construct first awk
if [ "$3" = "ss" ]
then #ss
	first_awk='$18}'
elif [ "$3" = "kts" ]
then #kts
	first_awk='$20}'
fi

middle_awk1=''
middle_awk2=''
second_awk=''
if [ "$4" = "0" ]
then #default
	second_awk='((length($4)==0)?".":$4)}'
	# no unit
	middle_awk1=''
	middle_awk2='$3'
elif [ "$4" = "u" ]
then # Section
	second_awk='((substr($4,length($4),1)==";")?".":$4)}'
	# add prefix for column
	middle_awk1='split($3,a,", Section"); '
	middle_awk2='a[1]'
	first_awk='a[1],";",'$first_awk
fi

# echo $first_awk
# echo $second_awk

awk -F'\t' 'BEGIN{OFS=""} {'"$middle_awk1"'print '"$first_column"',"\t",$2,"\t",'"$middle_awk2"',";",$4,";",$6,"\t",'"$first_awk"'' $1 > $tmp1
sed 1d $tmp1 > $tmp2
rm $tmp1
awk -F'\t' 'BEGIN{OFS=""} {print $1,"\t",$2,"\t",$3,"\t",'"$second_awk"'' $tmp2 > $2

rm $tmp2

echo "Done"

