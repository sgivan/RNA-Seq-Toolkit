#!/bin/bash
echo 'Running all tests.'

for dir in `ls -1d TestData_*`
do
    echo ""
    echo "$dir"
    cd $dir
    ./setup_and_test.sh
    ./test_result
    cd ..
done

