#!/bin/bash
for f in {0..9}
do
	echo " " >> summary.txt
	cat ex0$f.txt >> summary.txt
done
for f in {10..99}
do
	echo " " >> summary.txt
	cat ex$f.txt >> summary.txt
done
