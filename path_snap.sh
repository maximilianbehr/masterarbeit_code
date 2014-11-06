#!/bin/sh


#fenics conf file
FENICS_CONF="/Applications/FEniCS_snapshot.app/Contents/Resources/share/fenics/fenics.conf"

#source fenics
source $FENICS_CONF

#paraview
PVPYTHON="/Applications/paraview.app/Contents/bin"
PARAVIEW_PYTHONPATH="/Applications/paraview.app/Contents/Python/paraview"

export PATH=${PVPYTHON}:${PATH}
export PYTHONPATH=${PARAVIEWPYTHON}:${PYTHONPATH}



