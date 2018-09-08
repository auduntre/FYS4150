#!/bin/bash

echo "Compiling..."
make
echo "Done compiling"
echo "Running the calculations..."
./np
echo "Done"
echo "Making plots and tables..."
python3 uvplot.py
echo "Done"
echo "Printing table"
cat results/re_table.txt
echo "Done"
