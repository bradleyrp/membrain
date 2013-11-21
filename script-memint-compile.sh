#!/bin/sh

#touch memint.cpp
python setup.py build
cp build/lib.linux-x86_64-2.7/memint.so ./
echo
echo "Beginning script.py ...";
echo
python script-memint-test.py;
echo
echo "Finished running script.py.";
echo
