#!/bin/bash

#current dir from which script is executed and dir of script
CUR_DIR=`pwd`
SCRIPT_DIR=$(dirname $(readlink -f $0)) 
#PROC=$(nproc)
PROC=1

#set scratch dir and local usr dir for compiled libs
if [ `hostname` == "editha" ]
then
	SCRATCH_BEHR=/scratch/vol1/behr
elif [ `hostname` == "jack" ]
then
	SCRATCH_BEHR=/home/daniels/build
else
	SCRATCH_BEHR=/scratch/behr
fi
USR_LOCAL=${SCRATCH_BEHR}/usr/local
export PYTHONPATH=${USR_LOCAL}/lib/python2.7/site-packages:$PYTHONPATH
export LD_LIBRARY_PATH=${USR_LOCAL}/lib:${LD_LIBRARY_PATH}
export LIBRARY_PATH=${USR_LOCAL}/lib:${LIBRARY_PATH}

#print informations
echo "Current Dir 	= " ${CUR_DIR}
echo "Script Dir 	= " ${SCRIPT_DIR}
echo "Processes 	= " ${PROC}
echo "Scratch Dir 	= " ${SCRATCH_BEHR}
echo "usr/local 	= " ${USR_LOCAL}

mkdir -p ${USR_LOCAL} 
mkdir -p ${USR_LOCAL}/lib/python2.7/site-packages

#clear dir
function clear_dir()
{
	rm -rf ${SCRATCH_BEHR}
	mkdir -p ${USR_LOCAL} 
}

#checkout mess
function mess_co()
{   
    cd ${SCRATCH_BEHR};
    svn co http://svncsc.mpi-magdeburg.mpg.de/repos/mess/cmess/branches/hiwi_behr/;
}

#update cmess
function mess_up()
{   
    cd ${SCRATCH_BEHR};
	cd hiwi_behr;
	svn up;
}


#checkout mmess
function mmess_co()
{   
    cd ${SCRATCH_BEHR};
    svn co http://svncsc.mpi-magdeburg.mpg.de/repos/mess/mmess/trunk/ mmess;
}


#update mmess
function mmess_up()
{   
    cd ${SCRATCH_BEHR};
	cd mmess;
	svn up;
}


#build openblas
function build_openblas()
{
    cd ${SCRATCH_BEHR};
    wget http://github.com/xianyi/OpenBLAS/tarball/v0.2.14;
    tar xfv v0.2.14 ;
    rm v0.2.14;
    cd xianyi-OpenBLAS-2b0d8a8/; 
    make -j${PROC};
    make PREFIX=${USR_LOCAL} install;
    cd ${SCRATCH_BEHR}
}

#build SuiteSparse, .a and .so libs compiled with metis and openblas
function build_suitesparse()
{
    cd ${SCRATCH_BEHR};
    wget http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-4.4.4.tar.gz;
    tar xfv SuiteSparse-4.4.4.tar.gz;
    rm SuiteSparse-4.4.4.tar.gz;
    cd SuiteSparse/;
    
    #prepare metis and fix log2 naming bug
    wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.1.tar.gz;
    tar xfv metis-4.0.1.tar.gz;
    rm metis-4.0.1.tar.gz;
    sed -i "s/__log2/METIS__log2/g" metis-4.0/Lib/rename.h
    
    #prepare umfpack for so build
    cp ${SCRIPT_DIR}/GNUmakefile_umfpackLib 	UMFPACK/Lib/GNUmakefile
    cp ${SCRIPT_DIR}/GNUmakefile_amdLib 		AMD/Lib/GNUmakefile    
    cp ${SCRIPT_DIR}/Makefile_umf 				UMFPACK/Makefile
    cp ${SCRIPT_DIR}/Makefile_amd 				AMD/Makefile
    cp ${SCRIPT_DIR}/Makefile_metis 			metis-4.0/Makefile.in
	
	#set correct dirs to config file
	sed -i "s#\(INSTALL_LIB *= *\).*#\1${USR_LOCAL}/lib#g" 			${SCRIPT_DIR}/SuiteSparse_config.mk;  
	sed -i "s#\(INSTALL_INCLUDE *= *\).*#\1${USR_LOCAL}/include#g" 	${SCRIPT_DIR}/SuiteSparse_config.mk;
  	cp ${SCRIPT_DIR}/SuiteSparse_config.mk 							SuiteSparse_config/SuiteSparse_config.mk;
    
    #make and install
    make -j${PROC};
    make install;

    #install metis
    cp metis-4.0/libmetis.a ${USR_LOCAL}/lib
    cp metis-4.0/Lib/*.h ${USR_LOCAL}/include
    cd ${SCRATCH_BEHR}
}

#download and install nose
function build_nose()
{
    cd ${SCRATCH_BEHR}
    wget https://pypi.python.org/packages/source/n/nose/nose-1.3.4.tar.gz#md5=6ed7169887580ddc9a8e16048d38274d
    tar xfv nose-1.3.4.tar.gz
    rm nose-1.3.4.tar.gz
    cd nose-1.3.4/
    python setup.py build
    python setup.py install --prefix=${USR_LOCAL}
    cd ${SCRATCH_BEHR}		
}

#download and install build numpy
function build_numpy()
{
    cd ${SCRATCH_BEHR}
    wget http://sourceforge.net/projects/numpy/files/latest/download?source=files 
    tar xfv download\?source\=files 
    rm download\?source\=files 
    cd numpy-1.9.2/

	#set correct dirs to config file
	sed -i "s#\(library_dirs *= *\).*#\1${USR_LOCAL}/lib#g" 		${SCRIPT_DIR}/site.cfg;  
	sed -i "s#\(include_dirs *= *\).*#\1${USR_LOCAL}/include#g" 	${SCRIPT_DIR}/site.cfg;
  	cp ${SCRIPT_DIR}/site.cfg 										site.cfg;

    python setup.py build --fcompiler=gnu95
    python setup.py install --prefix=${USR_LOCAL}
    export PYTHONPATH=${USR_LOCAL}/lib/python2.7/site-packages:$PYTHONPATH
    cd ${SCRATCH_BEHR}
    python -c "import numpy;print numpy.__file__; numpy.test()"
}

#function build scipy
function build_scipy()
{
    cd ${SCRATCH_BEHR}
    wget http://sourceforge.net/projects/scipy/files/latest/download?source=files
    tar xfv download\?source\=files 
    rm download\?source\=files 
    cd scipy-0.15.1/
	
    export BLAS=${USR_LOCAL}/lib/libopenblas.a
    export LAPACK=${USR_LOCAL}/lib/libopenblas.a
    python setup.py build --fcompiler=gnu95
    python setup.py install --prefix=${USR_LOCAL}
	python -c "import scipy;print scipy.__file__; scipy.test()"
}

#function build scikit-umfpack
function build_scikitumfpack()
{
    cd ${SCRATCH_BEHR}
    wget https://pypi.python.org/packages/source/s/scikit-umfpack/scikit-umfpack-0.1.tar.gz#md5=22c9189115b52087097f03ff257453aa
    tar xfv scikit-umfpack-0.1.tar.gz 
    rm scikit-umfpack-0.1.tar.gz 
    cd scikit-umfpack-0.1/
	
    export BLAS=${USR_LOCAL}/lib/libopenblas.a
    export LAPACK=${USR_LOCAL}/lib/libopenblas.a
    python setup.py build 
    python setup.py install --prefix=${USR_LOCAL}
	
}

#function build mess
function build_mess()
{ 
	#debug build
    cd ${SCRATCH_BEHR}
	rm -rf mess_build_debug;
    mkdir mess_build_debug;
    cd mess_build_debug;
    cmake ../hiwi_behr -DDEBUG=ON -DPYTHON=ON -DBLA_VENDOR=OpenBLAS -DBLAS=${USR_LOCAL}/lib/libopenblas.so;
    make -j2; 
	cd python;
	python setup.2.7.py build;
	python setup.2.7.py install --prefix=${USR_LOCAL}
}

#clear_dir;
#mmess_co;
#mess_co;
#build_openblas;
#build_suitesparse;
#build_nose;
#build_numpy;
#build_scipy;
build_scikitumfpack;
#mess_up;
#mmess_up;
#build_mess;
