#!/bin/bash

argstring=""
for var in "$@"
do
	argstring+="$var "
done
python -i membrainrunner.py $argstring
