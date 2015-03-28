#!/bin/bash



export PYTHONPATH=`pwd`:${PYTHONPATH}


echo "Numpy in use:"
python -c "import numpy; print numpy.__file__; numpy.__config__.show()"

echo "Scipy in use:"
python -c "import scipy; print scipy.__file__; scipy.__config__.show()"

echo "Pycmess in use:"
python -c "import pycmess; print pycmess.__file__"