#!/bin/sh
 
#solvers="1.2 4.2 6.2 5.2 7.2"
solvers="1.2 4.2"
flag="0"
#cv="10,i 10,g 10,n"
cv="10,i 10,n"
dataset="a89"
#dataset="a89 b89"
#a89
kcmodel="kts" # ss
predict="1"
metrics="1"
#delimiters="# ~" # here # is a fake delimiter that is nowhere in our data # HARDCODED
regularization="0 4 2 1 0.5 0.25 0.125 0.0625 0.03125 0.015625"

#./trainhmm -s 1.2 -m 1 -p 1 -d ~ -v 10 data/a89_kts_train.txt > a89_kts_out_s1.2_v10.txt

for d in $dataset;
do
	for k in $kcmodel;
	do
		train_file=$d"_"$k"_train.txt"
		for s in $solvers; # solvers
		do
			for p in $predict; # predict
			do
                for m in $metrics; # metrics
                do
                    for f in $flag; # fit1skill
                    do
                        for c in $regularization; # regularization
                        do
#                            for d in $delimiters # delimiters
#                            do
                                for v in $cv;  # cros-validate
                                do
                                    model_file=$d"_"$k"_model_s"$s"_f"$f"_p"$p"_m"$m"_v"$v"_c"$c"_d~.txt"
                                    predict_file=$d"_"$k"_predict_s"$s"_f"$f"_p"$p"_m"$m"_v"$v"_c"$c"_d~.txt"
                                    out_file=$d"_"$k"_out_s"$s"_f"$f"_p"$p"_m"$m"_v"$v"_c"$c"_d~.txt"
                                    if [ "$v" = "0" ]
                                    then # no cross validation
                                        echo "./trainhmm -s "$s" -t 0.01 -q 0 -f "$f"         -p "$p" -m "$m" -c "$c" -d ~ "$train_file" "$model_file" "$predict_file" "$out_file
                                        ./trainhmm -s $s -t 0.01 -q 0 -f $f       -p $p -m $m -c $c -d ~ $train_file $model_file $predict_file &> $out_file 2>&1
                                    else # with cross validation
                                        echo "./trainhmm -s "$s" -t 0.01 -q 0 -f "$f" -v "$v" -p 1    -m "$m" -c "$c" -d ~ "$train_file" "$model_file" "$predict_file" "$out_file
                                        ./trainhmm -s $s -t 0.01 -q 0 -f $f -v $v -p 1  -m $m -c $c -d ~ $train_file $model_file $predict_file > $out_file 2>&1
                                    fi  
                                    gzip -9 $predict_file
                                    echo "..."
                                done # cros-validate
#                            done # delimiters
                        done # regularization
                    done # fit1skill
                done # metrics
			done # predict
		done # solvers
	done # kcmodel
done # dataset
