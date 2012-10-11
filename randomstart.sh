#!/bin/bash

echo "Starting..."

for i in {1..100}
do
	echo "Random run #$i"
	./trainhmm -t 0.01 -q 1 a89_ktt_train.txt a89_ktt_model$i.txt a89_ktt_predict$i.txt	
done

echo "Done"
