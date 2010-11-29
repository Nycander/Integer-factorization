#!/bin/sh

exec <tests.txt
fails=0
while read line
do
	# Run program with input $line
	ulimit -t 10

	output=`echo "$line\\n^C" | ./factoring.exe 1` 2> /dev/null

	if [ "$?" = "137" ]; then
		echo "\tFailed test '$line' because: TLE(10)"
		continue
	fi
	if [ "$output" = "fail" ]; then
		continue
	fi

	result=1
	resultstr="1"
	for i in $output; do
		# Check for prime number
		#isPrime=$(echo $i | ./checkForPrime.bsh)
		#if [ "$isPrime" = "$i is not prime" ]; then
		#	echo "\tFailed test '$line' because: $isPrime"
		#	break
		#fi

		resultstr="$resultstr*$i"
		result=$(echo "$result*$i" | bc)
	done

	# Check the multiplied output against the input value
	if [ $result != $line ]; then
		echo "\tFailed test '$line' because: ($resultstr=$result) != $line"
		fails=$(echo "$fails+1" | bc)
	fi
done

return $fails
