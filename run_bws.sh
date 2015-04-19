#!/bin/sh

#echo "BWS REF 1 lam=-0.25 default"
#python -u src/benchmarks/bws/run_default.py 1 -0.25 bws_lam_m025_ref_1

#echo "BWS REF 1 lam=0.0 default"
#python src/benchmarks/bws/run_default.py 1 0.0 bws_lam_0_ref_1

echo "BWS REF 2 lam=-0.25/0.0/0.25 default"
python -u src/benchmarks/bws/run_default.py 2 -0.25 bws_lam_m025_ref_2_RE1000 1000 > bws_lam_m025_ref_2_RE1000.log &
python -u src/benchmarks/bws/run_default.py 2 -0.25 bws_lam_m025_ref_2_RE1500 1500 > bws_lam_m025_ref_2_RE1500.log &
python -u src/benchmarks/bws/run_default.py 2 -0.25 bws_lam_m025_ref_2_RE2000 2000 > bws_lam_m025_ref_2_RE2000.log &
python -u src/benchmarks/bws/run_default.py 2 -0.25 bws_lam_m025_ref_2_RE2500 2300 > bws_lam_m025_ref_2_RE2300.log &
python -u src/benchmarks/bws/run_default.py 2 0.0 bws_lam_0_ref_2_RE1000 1000 > bws_lam_0_ref_2_RE1000.log &
python -u src/benchmarks/bws/run_default.py 2 0.0 bws_lam_0_ref_2_RE1500 1500 > bws_lam_0_ref_2_RE1500.log &
python -u src/benchmarks/bws/run_default.py 2 0.0 bws_lam_0_ref_2_RE2000 2000 > bws_lam_0_ref_2_RE2000.log &
python -u src/benchmarks/bws/run_default.py 2 0.0 bws_lam_0_ref_2_RE2500 2300 > bws_lam_0_ref_2_RE2300.log &
python -u src/benchmarks/bws/run_default.py 2 0.25 bws_lam_025_ref_2_RE1000 1000 > bws_lam_025_ref_2_RE1000.log &
python -u src/benchmarks/bws/run_default.py 2 0.25 bws_lam_025_ref_2_RE1500 1500 > bws_lam_025_ref_2_RE1500.log &
python -u src/benchmarks/bws/run_default.py 2 0.25 bws_lam_025_ref_2_RE2000 2000 > bws_lam_025_ref_2_RE2000.log &
python -u src/benchmarks/bws/run_default.py 2 0.25 bws_lam_025_ref_2_RE2500 2300 > bws_lam_025_ref_2_RE2300.log &


#echo "BWS REF 3 lam=-0.25 default"
#python -u src/benchmarks/bws/run_default.py 3 -0.25 bws_lam_m025_ref_3_RE1000 1000 > bws_lam_m025_ref_3_RE1000.log &
#python -u src/benchmarks/bws/run_default.py 3 -0.25 bws_lam_m025_ref_3_RE1500 1500 > bws_lam_m025_ref_3_RE1500.log &
#python -u src/benchmarks/bws/run_default.py 3 -0.25 bws_lam_m025_ref_3_RE2000 2000 > bws_lam_m025_ref_3_RE2000.log &
#python -u src/benchmarks/bws/run_default.py 3 -0.25 bws_lam_m025_ref_3_RE2500 2500 > bws_lam_m025_ref_3_RE2500.log &

#echo "BWS REF 3 lam=0.0 default"
#python -u src/benchmarks/bws/run_default.py 3 0.0 bws_lam_0_ref_3 > bws_lam_0_ref_3.log &
