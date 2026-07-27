#!/bin/bash
for i in $(seq 1 35)
do
   echo "Run number $i started"
   python run_30_12_25.py
   echo "Run $i completed"
done
