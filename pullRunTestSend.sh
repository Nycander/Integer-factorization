#!/bin/sh

echo "! Pulling latest revision from repository..."

git pull

echo "! Compiling... "

make -s

number_of_tests=0

echo "! Running tests... "

# TODO: Run a bunch of tests
./runTests.sh

failed_tests=$?

if [ $failed_tests = 0 ]; then
	echo "! Passed all tests! Sending to Kattis..."
	# TODO: Send to kattis
else
	echo "! Failed $failed_tests tests."
fi
