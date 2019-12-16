#!/bin/bash

for((n = 9; n <= 10; n++)); do
	echo "Running with n = $n"
	python3 sol2.py $n > "toplotzeta$n"
	echo "Points computed :D"
done
