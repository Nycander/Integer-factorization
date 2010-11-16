#!/bin/sh

echo "Pulling latest revision from repository..."

# TODO: Pull from git

echo "Running make... "

# TODO: Run "make"

number_of_tests=0
failed_tests=1

echo "Running tests... "

# TODO: Run a bunch of tests

if [ $failed_tests = 0 ]; then
	echo "Passed all tests! Sending to Kattis..."
	# TODO: Send to kattis
else
	echo "Failed $failed_tests tests."
fi

echo "Exiting."
