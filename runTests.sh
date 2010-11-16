#!/bin/sh

exec <tests.txt
fails=0
while read line
do
	# Run program with input $line
	output=$(echo "$line\\n^C" | ./factoring.exe)

	result=1
	resultstr="1"
	for i in $output
	do
		# TODO: Check prime number
		isPrime=$(echo $i | ./checkForPrime.bsh)
		if [ "$isPrime" = "$i is not prime" ]; then
			echo "\tFailed test '$line' because: $isPrime\n"
			break
		fi

		resultstr="$resultstr*$i"
		result=$(echo "$result*$i" | bc)
	done

	# TODO: Multiply output, check against input
	if [ $result != $line ]; then
		echo "\t Failed test '$line' because: ($resultstr=$result) != $line"
		fails=$(echo "$fails+1" | bc)
	fi
done

return $fails
