#!/bin/bash
FILES=~/git_repo/BorderApolarity/RESULTS3_15/*
for f in $FILES; do
	line=$(head -n 1 $f)
	let "line += 1"
	num=$(< $f wc -l)
	if [[ $line -eq $num ]]; then
		:
	else
		echo $f
		echo $line
		echo $num
	fi
done
		

