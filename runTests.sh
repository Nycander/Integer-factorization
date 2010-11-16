#!/bin/sh

exec <tests.txt

while read line
do
	echo $line
	# TODO: Run program with input $line
	# TODO: Check prime number
	# TODO: Multiply output, check against input
done
