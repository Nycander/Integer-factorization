#!/bin/sh

echo "! Pulling latest revision from repository..."

git pull

echo "! Compiling... "

make -s

echo "! Running tests... "

# Run a bunch of tests
./runTests.sh

failed_tests=$?
number_of_tests=$(cat tests.txt | wc -l)

if [ $failed_tests = 0 ]; then
	echo "! Passed all $number_of_tests tests!"
	# Send to kattis
	yes | ./submit.py factoring.c > /dev/null
	echo "! Sent mail to kattis :D"
else
	echo "! Failed $failed_tests / $number_of_tests tests."
fi
