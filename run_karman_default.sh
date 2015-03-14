#!/bin/sh


echo "KARMAN REF 1 DEFAULT"
python src/benchmarks/karman/run_default.py 1 karman_ref_1 > karman_ref_1.log &

echo "KARMAN REF 2 DEFAULT"
python src/benchmarks/karman/run_default.py 2 karman_ref_2 > karman_ref_2.log &

echo "KARMAN REF 3 DEFAULT"
python src/benchmarks/karman/run_default.py 3 karman_ref_3 > karman_ref_3.log &

echo "KARMAN REF 4 DEFAULT"
python src/benchmarks/karman/run_default.py 4 karman_ref_4 > karman_ref_4.log &




