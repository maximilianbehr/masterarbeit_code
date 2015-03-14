#!/bin/sh


echo "KARMAN REF 1 BERNOULLI"
python src/benchmarks/karman/run_bernoulli_stable.py 1 karman_bern_ref_1 > karman_bern_ref_1.log &

echo "KARMAN REF 2 BERNOULLI"
python src/benchmarks/karman/run_bernoulli_stable.py 2 karman_bern_ref_2 > karman_bern_ref_2.log &

echo "KARMAN REF 3 BERNOULLI"
python src/benchmarks/karman/run_bernoulli_stable.py 3 karman_bern_ref_3 > karman_bern_ref_3.log &

echo "KARMAN REF 4 BERNOULLI"
python src/benchmarks/karman/run_bernoulli_stable.py 4 karman_bern_ref_4 > karman_bern_ref_4.log &




