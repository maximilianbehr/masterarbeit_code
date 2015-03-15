#!/bin/sh


#echo "KARMAN REF 1 VERTICAL"
#python src/benchmarks/karman/run_vertical.py 1 karman_vertical_ref_1 > karman_vertical_ref_1.log &
#python src/benchmarks/karman/run_vertical.py 1 karman_vertical_ref_1


#echo "KARMAN REF 2 VERTICAL"
#python src/benchmarks/karman/run_vertical.py 2 karman_vertical_ref_2 > karman_vertical_ref_2.log &

echo "KARMAN REF 3 VERTICAL"
python -u src/benchmarks/karman/run_vertical.py 3 karman_vertical_ref_3 > karman_vertical_ref_3.log &

#echo "KARMAN REF 4 VERTICAL"
#python src/benchmarks/karman/run_vertical.py 4 karman_vertical_ref_4 > karman_vertical_ref_4.log &




