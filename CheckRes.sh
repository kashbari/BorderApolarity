#!/bin/bash
FILES=/home/kashbari/git_repo/BorderApolarity/RESULTS3_15
for f in $FILES
do
	line=$(head -n 1 f)
	let "line += 1"
	num=$(wc -l f)
	if [line != num]
	then
		echo $f
	fi
done
		
