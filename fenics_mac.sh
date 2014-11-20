#!/bin/sh



if [ $# -eq 0  -o $#  -gt 1 ]
  then
    echo "usage: . ./fenics_mac.sh 1.2.0|1.3.0|snapshot";
    return 1
fi


#fenics conf file
if [ $1 == "1.2.0" ]
    then
    FENICS_CONF="/Applications/FEniCS_1.2.0.app/Contents/Resources/share/fenics/fenics.conf"
elif [ $1 == "1.3.0" ]
    then
    FENICS_CONF="/Applications/FEniCS_1.3.0.app/Contents/Resources/share/fenics/fenics.conf"
elif [ $1 == "snapshot" ]
    then
    FENICS_CONF="/Applications/FEniCS_snapshot.app/Contents/Resources/share/fenics/fenics.conf"
else
    echo $1 unknown FEniCS version
    return 1
fi


#source fenics
source $FENICS_CONF

#paraview
PVPYTHON="/Applications/paraview.app/Contents/bin"
PARAVIEW_PYTHONPATH="/Applications/paraview.app/Contents/Python/paraview"

export PATH=${PVPYTHON}:${PATH}
export PYTHONPATH=${PARAVIEWPYTHON}:`pwd`:${PYTHONPATH}



