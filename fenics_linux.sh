#!/bin/sh
source venv/bin/activate
export PYTHONPATH=`pwd`:${PYTHONPATH}
export PYTHONPATH=/home/daniels/build/usr/local/lib/python2.7/site-packages:${PYTHONPATH}
export LD_LIBRARY_PATH=/home/daniels/build/mess_build_debug/lib:${LD_LIBRARY_PATH}

echo "Numpy in use:"
python -c "import numpy; print numpy.__file__; numpy.__config__.show()"

echo "Scipy in use:"
python -c "import scipy; print scipy.__file__; scipy.__config__.show()"

echo "Pycmess in use:"
python -c "import pycmess; print pycmess.__file__"
