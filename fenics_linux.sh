#!/bin/sh
source venv/bin/activate
export PYTHONPATH=`pwd`:${PYTHONPATH}
export LD_LIBRARY_PATH=/home/daniels/workspace/mess_debug/lib/:$LD_LIBRARY_PATH

