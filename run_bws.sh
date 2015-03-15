#!/bin/sh

#echo "BWS REF 1 lam=-0.25 default"
#python -u src/benchmarks/bws/run_default.py 1 -0.25 bws_lam_m025_ref_1

#echo "BWS REF 1 lam=0.0 default"
#python src/benchmarks/bws/run_default.py 1 0.0 bws_lam_0_ref_1



echo "BWS REF 3 lam=-0.25 default"
python -u src/benchmarks/bws/run_default.py 3 -0.25 bws_lam_m025_ref_3 > bws_lam_m025_ref_3.log &

echo "BWS REF 3 lam=0.0 default"
python -u src/benchmarks/bws/run_default.py 3 0.0 bws_lam_0_ref_3 > bws_lam_0_ref_3.log &
