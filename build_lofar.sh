# General environment settings.
export J=64
export INSTALLDIR=$HOME/opt/lofar_200124

# CEP3
#export cmake=/home/dijkema/opt/cmake/bin/cmake
# Leiden
#module load gcc/8.1.0
#module load make/4.2
#module load cmake/3.9
#export PYTHON_VERSION=3.6
#export PYTHON_VERSION_NODOT=36
# Hamburg
export PYTHON_VERSION=3.5
export PYTHON_VERSION_NODOT=35

# General compile and build settings.
export make=`which make`
export cmake=`which cmake`
export CC=`which gcc`
export CXX=`which g++`
export CFLAGS="-D_GLIB_USE_CXX_ABI=1 -DBOOST_NO_CXX11_SCOPED_ENUMS"
export CXXFLAGS="-D_GLIB_USE_CXX_ABI=1 -DBOOST_NO_CXX11_SCOPED_ENUMS"
# --std=c++11 -W -Wall -Woverloaded-virtual -Wno-unknown-pragmas -D_GLIBCXX_USE_CXX11_ABI=${CXX_ABI} -O3 -march=haswell"

# Path to where the patch for python-casacore's setup is stored.
#export PYTHON_CASACORE_PATCH=$HOME/opt/src/patch_python-casacore.patch
#export PATCH_AOFLAGGER=/net/lofar1/data1/sweijen/software/LOFAR/aoflagger.patch

# Settings relevant to the installed software.
#export AOFLAGGER_VERSION=latest
export AOFLAGGER_VERSION=v2.14.0
export ARMADILLO_VERSION=8.600.0
export BOOST_DOT_VERSION=1.67.0
export BOOST_VERSION=1_67_0
#export CASACORE_VERSION=v2.4.1
export CASACORE_VERSION=latest
# Leave at latest, release versions crash for some reason.
export CFITSIO_VERSION=3.47
export DYSCO_VERSION=v1.2.0
export FFTW_VERSION=3.3.8
export HDF5_VERSION=1.8.21
export LOSOTO_VERSION=2.0
export OPENBLAS_VERSION=v0.3.7
export PYBDSF_VERSION=v1.8.12
export PYTHON_CASACORE_VERSION=v3.1.1
# Do not change, Armadillo wants this version of SuperLU.
export SUPERLU_VERSION=v5.2.1
export WSCLEAN_VERSION=latest
export WCSLIB_VERSION=latest

mkdir -p $INSTALLDIR

#######################################
# Build all external libraries first. #
#######################################

if [ ! -d $INSTALLDIR/fftw ]; then
    #
    # Install FFTW
    #
    echo Installing FFTW...
    mkdir -p ${INSTALLDIR}/fftw
    cd ${INSTALLDIR}/fftw && wget http://www.fftw.org/fftw-${FFTW_VERSION}.tar.gz
    cd ${INSTALLDIR}/fftw && tar xzf fftw*.tar.gz
    cd ${INSTALLDIR}/fftw/fftw*/ && ./configure --prefix=${INSTALLDIR}/fftw --enable-threads --enable-shared
    cd ${INSTALLDIR}/fftw/fftw*/ && $make -j $J
    cd ${INSTALLDIR}/fftw/fftw*/ && $make install
    cd ${INSTALLDIR}/fftw/fftw*/ && ./configure --prefix=${INSTALLDIR}/fftw --enable-threads --enable-shared --enable-float
    cd ${INSTALLDIR}/fftw/fftw*/ && $make -j $J
    cd ${INSTALLDIR}/fftw/fftw*/ && $make install
    echo Installed FFTW.
else
    echo FFTW already installed.
fi

if [ ! -d $INSTALLDIR/hdf5 ]; then
    #
    # Install HDF5
    #
    echo Installing HDF5...
    mkdir -p ${INSTALLDIR}/hdf5
    cd ${INSTALLDIR}/hdf5 && wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${HDF5_VERSION%.*}/hdf5-${HDF5_VERSION}/src/hdf5-${HDF5_VERSION}.tar.gz
    cd ${INSTALLDIR}/hdf5 && tar xzf hdf5*.tar.gz
    cd ${INSTALLDIR}/hdf5/hdf5*/ && ./configure --prefix=${INSTALLDIR}/hdf5 --enable-fortran --enable-cxx
    cd ${INSTALLDIR}/hdf5/hdf5*/ && $make -j $J
    cd ${INSTALLDIR}/hdf5/hdf5*/ && $make install
    echo Installed HDF5.
else
    echo HDF5 already installed.
fi

export CMAKE_PREFIX_PATH=${INSTALLDIR}/hdf5

if [ ! -d $INSTALLDIR/boost ]; then
    #
    # Install Boost.Python
    #
    echo Installing Boost.Python...
    echo "using python : ${PYTHON_VERSION} : /usr/bin/python3 : /usr/include/python${PYTHON_VERSION} : /usr/lib ;" > $HOME/user-config.jam
    mkdir -p $INSTALLDIR/boost/src
    cd $INSTALLDIR/boost/ && wget https://dl.bintray.com/boostorg/release/${BOOST_DOT_VERSION}/source/boost_${BOOST_VERSION}.tar.gz && tar xzf boost_${BOOST_VERSION}.tar.gz
    cd $INSTALLDIR/boost/boost_*/ && ./bootstrap.sh --prefix=$INSTALLDIR/boost --with-python-version=${PYTHON_VERSION} && ./b2 headers && ./b2 install toolset=gcc cxxflags=-std=c++11 --prefix=$INSTALLDIR/boost --with-atomic --with-chrono --with-date_time --with-filesystem --with-program_options --with-python --with-signals --with-test --with-thread -j $J
    # CHECK THE PYTHON VERSION HERE
    cd $INSTALLDIR/boost/lib/
    ln -s libboost_python${PYTHON_VERSION_NODOT}.a         libboost_python.a
    ln -s libboost_python${PYTHON_VERSION_NODOT}.so        libboost_python.so
    ln -s libboost_python${PYTHON_VERSION_NODOT}.so.1.67.0 libboost_python.so.1.67.0
    ln -s libboost_python${PYTHON_VERSION_NODOT}.a         libboost_python3.a
    ln -s libboost_python${PYTHON_VERSION_NODOT}.so        libboost_python3.so
    ln -s libboost_python${PYTHON_VERSION_NODOT}.so.1.67.0 libboost_python3.so.1.67.0
    ln -s libboost_numpy${PYTHON_VERSION_NODOT}.a          libboost_numpy.a
    ln -s libboost_numpy${PYTHON_VERSION_NODOT}.so         libboost_numpy.so
    ln -s libboost_numpy${PYTHON_VERSION_NODOT}.so.1.67.0  libboost_numpy.so.1.67.0
    ln -s libboost_numpy${PYTHON_VERSION_NODOT}.a          libboost_numpy3.a
    ln -s libboost_numpy${PYTHON_VERSION_NODOT}.so         libboost_numpy3.so
    ln -s libboost_numpy${PYTHON_VERSION_NODOT}.so.1.67.0  libboost_numpy3.so.1.67.0
    echo Installed Boost.Python.
else
    echo Boost.Python already installed.
fi

if [ ! -d $INSTALLDIR/openblas ]; then
    #
    # Install OpenBLAS
    #
    echo Installing OpenBLAS...
    mkdir -p $INSTALLDIR/openblas/
    cd $INSTALLDIR/openblas/ && git clone https://github.com/xianyi/OpenBLAS.git src && cd src && git checkout $OPENBLAS_VERSION
    cd $INSTALLDIR/openblas/src && $make -j $J 
    cd $INSTALLDIR/openblas/src && $make install PREFIX=$INSTALLDIR/openblas
    echo Installed OpenBLAS.
else
    echo OpenBlas already installed.
fi

export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:${INSTALLDIR}/openblas

if [ ! -d $INSTALLDIR/superlu ]; then
    #
    # Install SuperLU
    #
    echo Installing SuperLU...
    mkdir -p $INSTALLDIR/superlu/build
    cd $INSTALLDIR/superlu/ && git clone https://github.com/xiaoyeli/superlu.git src && cd src && git checkout $SUPERLU_VERSION
    cd $INSTALLDIR/superlu/build && $cmake -DCMAKE_CXX_FLAGS=-D_GLIBCXX_USE_CXX11_ABI=1 -DCMAKE_INSTALL_PREFIX=$INSTALLDIR/superlu -DUSE_XSDK_DEFAULTS=TRUE -Denable_blaslib=OFF -DBLAS_LIBRARY=$INSTALLDIR/openblas/lib/libopenblas.so ../src
    cd $INSTALLDIR/superlu/build && $make -j $J 
    cd $INSTALLDIR/superlu/build && $make install
    echo Installed SuperLU.
else
    echo SuperLu already installed.
fi

if [ ! -d $INSTALLDIR/armadillo ]; then
    #
    # Install Armadillo
    #
    echo Installing Armadillo...
    # check if superlu lib is .../lib or .../lib64
    mkdir -p $INSTALLDIR/armadillo/
    cd $INSTALLDIR/armadillo && wget http://sourceforge.net/projects/arma/files/armadillo-$ARMADILLO_VERSION.tar.xz 
    cd $INSTALLDIR/armadillo && tar xf armadillo-$ARMADILLO_VERSION.tar.xz
    cd $INSTALLDIR/armadillo/armadillo*/ && ./configure && $cmake -DCMAKE_CXX_FLAGS=-D_GLIBCXX_USE_CXX11_ABI=1 -DCMAKE_INSTALL_PREFIX:PATH=$INSTALLDIR/armadillo -Dopenblas_LIBRARY:FILEPATH=$INSTALLDIR/openblas/lib/libopenblas.so  -DSuperLU_INCLUDE_DIR:PATH=$INSTALLDIR/superlu/include -DSuperLU_LIBRARY:FILEPATH=$INSTALLDIR/superlu/lib/libsuperlu.so	
    cd $INSTALLDIR/armadillo/armadillo*/ && $make -j $J
    cd $INSTALLDIR/armadillo/armadillo*/ && $make install
    echo Installed Armadillo.
else
    echo Armadillo already installed.
fi

if [ ! -d $INSTALLDIR/cfitsio ]; then
    #
    # Install-cfitsio
    #
    echo Installing CFITSIO...
    mkdir -p ${INSTALLDIR}/cfitsio/build
    cd ${INSTALLDIR}/cfitsio && wget http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-${CFITSIO_VERSION}.tar.gz
    cd ${INSTALLDIR}/cfitsio && tar xf cfitsio-${CFITSIO_VERSION}.tar.gz
    cd ${INSTALLDIR}/cfitsio/build && $cmake -DCMAKE_CXX_FLAGS=-D_GLIBCXX_USE_CXX11_ABI=1 -DCMAKE_INSTALL_PREFIX=${INSTALLDIR}/cfitsio/ ../cfitsio-${CFITSIO_VERSION}
    cd ${INSTALLDIR}/cfitsio/build && $make -j $J
    cd ${INSTALLDIR}/cfitsio/build && $make install
    echo Installed CFITSIO.
else
    echo CFITSIO already installed.
fi

if [ ! -d $INSTALLDIR/wcslib ]; then
    #
    # Install-wcslib
    #
    echo Installing WCSLIB...
    mkdir ${INSTALLDIR}/wcslib
    if [ "${WCSLIB_VERSION}" = "latest" ]; then cd ${INSTALLDIR}/wcslib && wget --retry-connrefused ftp://anonymous@ftp.atnf.csiro.au/pub/software/wcslib/wcslib.tar.bz2 -O wcslib-latest.tar.bz2; fi
    if [ "${WCSLIB_VERSION}" != "latest" ]; then cd ${INSTALLDIR}/wcslib && wget --retry-connrefused ftp://anonymous@ftp.atnf.csiro.au/pub/software/wcslib/wcslib-${WCSLIB_VERSION}.tar.bz2; fi
    cd ${INSTALLDIR}/wcslib && tar xf wcslib-*.tar.bz2
    #cd ${INSTALLDIR} && mkdir wcslib && cd wcslib && svn checkout https://github.com/astropy/astropy/trunk/cextern/wcslib
    cd ${INSTALLDIR}/wcslib/wcslib* && ./configure --prefix=${INSTALLDIR}/wcslib --with-cfitsiolib=${INSTALLDIR}/cfitsio/lib/ --with-cfitsioinc=${INSTALLDIR}/cfitsio/include/ --without-pgplot
    cd ${INSTALLDIR}/wcslib/wcslib* && $make -j $J
    cd ${INSTALLDIR}/wcslib/wcslib* && $make install -j $J
    echo Installed WCSLIB.
else
    echo WCSLIB already installed.
fi


# Make library and header directories easily available for future CMakes.
export CMAKE_INCLUDE=$INSTALLDIR/fftw/include:$INSTALLDIR/armadillo/include:$INSTALLDIR/boost/include:$INSTALLDIR/cfitsio/include:$INSTALLDIR/openblas/include:$INSTALLDIR/superlu/include:$INSTALLDIR/wcslib/include:$INSTALLDIR/hdf5/include:$INSTALLDIR/aoflagger/include
export CMAKE_LIBRARY=$INSTALLDIR/fftw/lib:$INSTALLDIR/armadillo/lib:$INSTALLDIR/boost/lib:$INSTALLDIR/cfitsio/lib:$INSTALLDIR/openblas/lib:$INSTALLDIR/superlu/lib64:$INSTALLDIR/wcslib/lib:$INSTALLDIR/hdf5/lib
export CMAKE_PREFIX_PATH=$INSTALLDIR/armadillo:$INSTALLDIR/boost:$INSTALLDIR/casacore:$INSTALLDIR/cfitsio:$INSTALLDIR/dysco:$INSTALLDIR/idg:$INSTALLDIR/openblas:$INSTALLDIR/superlu:$INSTALLDIR/wcslib:$INSTALLDIR/hdf5:/usr/lib64
export CPATH=$INSTALLDIR/armadillo/lib:$INSTALLDIR/boost/lib:$INSTALLDIR/cfitsio/lib:$INSTALLDIR/openblas/lib:$INSTALLDIR/superlu/lib:$INSTALLDIR/wcslib/lib

#########################################
# Install main LOFAR software packages. #
#########################################

if [ ! -d $INSTALLDIR/casacore ]; then
    #
    # Install CASAcore
    #
    # TODO: why do i need to use the other cmake here? the new one doesn't recognize python numpy...
    echo Installing CASAcore...
    mkdir -p ${INSTALLDIR}/casacore/build
    mkdir -p ${INSTALLDIR}/casacore/data
    cd $INSTALLDIR/casacore && git clone https://github.com/casacore/casacore.git src
    if [ "${CASACORE_VERSION}" != "latest" ]; then cd ${INSTALLDIR}/casacore/src && git checkout tags/${CASACORE_VERSION}; fi
    cd ${INSTALLDIR}/casacore/data && wget --retry-connrefused ftp://anonymous@ftp.astron.nl/outgoing/Measures/WSRT_Measures.ztar
    cd ${INSTALLDIR}/casacore/data && tar xf WSRT_Measures.ztar
    cd ${INSTALLDIR}/casacore/build && $cmake -DCMAKE_CXX_FLAGS=${CXX_FLAGS} -DCMAKE_BUILD_TYPE=Release -DCMAKE_LIBRARY_PATH=$CMAKE_LIBRARY -DCMAKE_INCLUDE_PATH=$CMAKE_INCLUDE -DCMAKE_INSTALL_PREFIX=${INSTALLDIR}/casacore/ -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DDATA_DIR:PATH=${INSTALLDIR}/casacore/data -DBoost_DIR=$INSTALLDIR/boost -DBOOST_PYTHON3_LIBRARY_NAME:STRING=python -DFFTW3F_LIBRARY:FILEPATH=$INSTALLDIR/fftw/lib/libfftw3f.so -DFFTW3F_THREADS_LIBRARY:FILEPATH=$INSTALLDIR/fftw/lib/libfftw3f_threads.so -DFFTW3_INCLUDE_DIR:PATH=$INSTALLDIR/fftw/include -DFFTW3_LIBRARY:FILEPATH=$INSTALLDIR/fftw/lib/libfftw3.so -DFFTW3_THREADS_LIBRARY:FILEPATH=$INSTALLDIR/fftw/lib/libfftw3_threads.so -DHDF5_hdf5_LIBRARY:FILEPATH=$INSTALLDIR/hdf5/lib/libhdf5.so -DENABLE_TABLELOCKING=OFF -DUSE_OPENMP=True -DUSE_FFTW3=TRUE -DUSE_HDF5=True -DBUILD_PYTHON:BOOL=False -DBUILD_PYTHON3:BOOL=True ../src/
    cd ${INSTALLDIR}/casacore/build && $make -j ${J}
    cd ${INSTALLDIR}/casacore/build && $make install
    echo Installed CASAcore.
else
    echo CASAcore already installed.
fi

if [ ! -d $INSTALLDIR/python-casacore ]; then
    #
    # Install-python-casacore
    #
    # Finding libraries is broken, patch the setup to include the previously installed boost and casacore libraries.
    export CFLAGS="-std=c++11"
    mkdir ${INSTALLDIR}/python-casacore
    cd ${INSTALLDIR}/python-casacore && git clone https://github.com/casacore/python-casacore
    if [ "$PYTHON_CASACORE_VERSION" != "latest" ]; then cd ${INSTALLDIR}/python-casacore/python-casacore && git checkout tags/${PYTHON_CASACORE_VERSION}; fi
    # OLD: patch setup.py $PYTHON_CASACORE_PATCH
    mkdir -p ${INSTALLDIR}/python-casacore/lib64/python${PYTHON_VERSION}/site-packages/  ${INSTALLDIR}/python-casacore/lib/python${PYTHON_VERSION}/site-packages/
    export PYTHONPATH=$INSTALLDIR/python-casacore/lib/python${PYTHON_VERSION}/site-packages/:$INSTALLDIR/DP3/lib/python${PYTHON_VERSION}/site-packages/:$PYTHONPATH
    export PYTHONPATH=$INSTALLDIR/python-casacore/lib64/python${PYTHON_VERSION}/site-packages/:$INSTALLDIR/DP3/lib64/python${PYTHON_VERSION}/site-packages/:$PYTHONPATH
    cd ${INSTALLDIR}/python-casacore/python-casacore && ./setup.py build_ext --swig-cpp -I${INSTALLDIR}/wcslib/include:${INSTALLDIR}/casacore/include/:${INSTALLDIR}/cfitsio/include:${INSTALLDIR}/boost/include -L${INSTALLDIR}/wcslib/lib:${INSTALLDIR}/casacore/lib/:${INSTALLDIR}/cfitsio/lib/:${INSTALLDIR}/boost/lib:/usr/lib/python3.5 install --prefix ${INSTALLDIR}/python-casacore/

else
    echo Python-CASAcore already installed.
fi

if [ ! -d $INSTALLDIR/dysco ]; then
    #
    # Install Dysco
    #
    echo Installing Dysco...
    mkdir -p $INSTALLDIR/dysco/build
    cd $INSTALLDIR/dysco && git clone https://github.com/aroffringa/dysco.git src
    cd $INSTALLDIR/dysco/build && $cmake -DCMAKE_CXX_FLAGS=-D_GLIB_USE_CXX_ABI=1 -DCMAKE_INSTALL_PREFIX=$INSTALLDIR/dysco -DCASACORE_ROOT_DIR=$INSTALLDIR/casacore -DCMAKE_LIBRARY_PATH=$CMAKE_LIBRARY -DCMAKE_INCLUDE_PATH=$CMAKE_INCLUDE -DCMAKE_Boost_DIR:PATH=$INSTALLDIR/Boost/ ../src
    cd $INSTALLDIR/dysco/build && $make -j $J 
    cd $INSTALLDIR/dysco/build && $make install
    echo Installed Dysco.
else
    echo Dysco already installed.
fi

if [ ! -d $INSTALLDIR/aoflagger ]; then
    #
    # Install-aoflagger
    #
    # might have to comment these in CMakeLists.txt
    # > pkg_check_modules(GTKMM gtkmm-3.0>=3.0.0)
    # > pkg_check_modules(SIGCXX sigc++-2.0)
    mkdir -p ${INSTALLDIR}/aoflagger/build
    if [ "${AOFLAGGER_VERSION}" = "latest" ]; then cd ${INSTALLDIR}/aoflagger && git clone git://git.code.sf.net/p/aoflagger/code aoflagger && cd ${INSTALLDIR}/aoflagger/aoflagger; fi
    if [ "${AOFLAGGER_VERSION}" != "latest" ]; then cd ${INSTALLDIR}/aoflagger && git clone git://git.code.sf.net/p/aoflagger/code aoflagger && cd ${INSTALLDIR}/aoflagger/aoflagger && git checkout tags/${AOFLAGGER_VERSION}; fi
    cd ${INSTALLDIR}/aoflagger/build && $cmake -DFFTW3_INCLUDE_DIR:PATH=${INSTALLDIR}/fftw/include -DFFTW3_LIB:FILEPATH=${INSTALLDIR}/fftw/lib/libfftw3.so -DCASACORE_ROOT_DIR=$INSTALLDIR/casacore -DCMAKE_CXX_FLAGS=-D_GLIB_USE_CXX_ABI=1 -DCMAKE_INSTALL_PREFIX=${INSTALLDIR}/aoflagger/ -DBUILD_SHARED_LIBS=ON -DPORTABLE=True ../aoflagger
    cd ${INSTALLDIR}/aoflagger/build && make -j $J
    cd ${INSTALLDIR}/aoflagger/build && make install
else
    echo AOFlagger already installed.
fi

 
if [ ! -d $INSTALLDIR/LOFARBeam ]; then
    #
    # Install the standalone StationResponse libraries.
    #
    echo Installing LOFARBeam...
    mkdir -p $INSTALLDIR/LOFARBeam/build
    cd $INSTALLDIR/LOFARBeam && git clone https://github.com/lofar-astron/LOFARBeam.git src
    cd $INSTALLDIR/LOFARBeam/build && $cmake -DCMAKE_INSTALL_PREFIX=$INSTALLDIR/LOFARBeam ../src
    cd $INSTALLDIR/LOFARBeam/build && $make -j $J
    cd $INSTALLDIR/LOFARBeam/build && $make install
fi


if [ ! -d $INSTALLDIR/idg ]; then
    #
    # Install the Image Domain Gridder (IDG).
    # 
    echo Installing IDG.
    mkdir -p $INSTALLDIR/idg/build
    cd $INSTALLDIR/idg && git clone https://gitlab.com/astron-idg/idg.git src
    cd $INSTALLDIR/idg/build && $cmake -DCMAKE_INSTALL_PREFIX=$INSTALLDIR/idg -DBUILD_WITH_MKL=ON ../src
    cd $INSTALLDIR/idg/build && $make -j $J
    cd $INSTALLDIR/idg/build && $make install
else
    echo IDG already installed.
fi


if [ ! -d $INSTALLDIR/DP3 ]; then
    #
    # Install DP3.
    #
    echo Installing DP3...
    #export LD_LIBRARY_PATH=/net/lofar1/data1/sweijen/software/HDF5_1.8/lib:/net/lofar1/data1/sweijen/software/LOFAR/2018_11_05_DP3/superlu/lib64:$INSTALLDIR/LOFARBeam/lib:$LD_LIBRARY_PATH
    mkdir -p $INSTALLDIR/DP3/build
    git clone https://github.com/lofar-astron/DP3.git $INSTALLDIR/DP3/src
    cd $INSTALLDIR/DP3/build && $cmake -DCMAKE_CXX_FLAGS="-D_GLIB_USE_CXX_ABI=1 -DBOOST_NO_CXX11_SCOPED_ENUMS" -DCMAKE_INSTALL_PREFIX:PATH=$INSTALLDIR/DP3 -DHDF5_DIR:PATH=${INSTALLDIR}/hdf5/hdf5-${HDF5_VERSION} -DIDGAPI_INCLUDE_DIRS:PATH=$INSTALLDIR/idg/include -DIDGAPI_LIBRARIES:PATH=$INSTALLDIR/idg/lib/libidg-api.so -DLOFAR_STATION_RESPONSE_DIR:PATH=$INSTALLDIR/LOFARBeam/include -DLOFAR_STATION_RESPONSE_LIB:FILEPATH=$INSTALLDIR/LOFARBeam/lib/libstationresponse.so -DAOFLAGGER_INCLUDE_DIR:PATH=$INSTALLDIR/aoflagger/include -DAOFLAGGER_LIB:FILEPATH=$INSTALLDIR/aoflagger/lib/libaoflagger.so ../src
    cd $INSTALLDIR/DP3/build && $make -j $J
    cd $INSTALLDIR/DP3/build && $make install
    echo Installed DP3.
else
    echo DP3 already installed.
fi


if [ ! -d $INSTALLDIR/wsclean ]; then
    #
    # Install WSClean
    #
    # Sometimes there are FFT to update by hand in CMakeCache.txt
    echo Installing WSClean.
    mkdir -p $INSTALLDIR/wsclean/build
    if [ "$WSCLEAN_VERSION" != "latest" ]; then cd ${INSTALLDIR}/wsclean && wget http://downloads.sourceforge.net/project/wsclean/wsclean-${WSCLEAN_VERSION}/wsclean-${WSCLEAN_VERSION}.tar.bz2 && tar -xjf wsclean-${WSCLEAN_VERSION}.tar.bz2 && cd wsclean-${WSCLEAN_VERSION}; fi
    if [ "$WSCLEAN_VERSION" = "latest" ]; then cd ${INSTALLDIR}/wsclean && git clone git://git.code.sf.net/p/wsclean/code src && cd src/wsclean; fi
    cd $INSTALLDIR/wsclean/build && $cmake -DCMAKE_INSTALL_PREFIX=$INSTALLDIR/wsclean -DIDGAPI_INCLUDE_DIRS:PATH=$INSTALLDIR/idg/include -DIDGAPI_LIBRARIES:PATH=$INSTALLDIR/idg/lib/libidg-api.so -DLOFAR_STATION_RESPONSE_DIR:PATH=$INSTALLDIR/LOFARBeam/include -DLOFAR_STATION_RESPONSE_LIB:FILEPATH=$INSTALLDIR/LOFARBeam/lib/libstationresponse.so -DFFTW3_INCLUDE_DIR:PATH=${INSTALLDIR}/fftw/include -DFFTW3_LIB:FILEPATH=${INSTALLDIR}/fftw/lib/libfftw3.so -DFFTW3_THREADS_LIB:FILEPATH=${INSTALLDIR}/fftw/lib/libfftw3_threads.so -DFFTW3F_LIB:FILEPATH=${INSTALLDIR}/fftw/lib/libfftw3f.so -DFFTW3F_THREADS_LIB:FILEPATH=${INSTALLDIR}/fftw/lib/libfftw3f_threads.so ../src/wsclean
    cd $INSTALLDIR/wsclean/build && $make -j $J
    cd $INSTALLDIR/wsclean/build && $make install
else
    echo WSClean already installed.
fi

if [ ! -d $INSTALLDIR/pyBDSF ]; then
    echo Installing pyBDSF.
    #
    # Install pyBDSF.
    #
    mkdir -p $INSTALLDIR/pyBDSF/lib64/python$PYTHON_VERSION/site-packages
    export PYTHONPATH=$INSTALLDIR/pyBDSF/lib64/python$PYTHON_VERSION/site-packages:$PYTHONPATH
    cd $INSTALLDIR/pyBDSF && git clone https://github.com/lofar-astron/PyBDSF pyBDSF
    cd $INSTALLDIR/pyBDSF/pyBDSF && python setup.py install --prefix=$INSTALLDIR/pyBDSF
else
    echo pyBDSF already installed.
fi

if [ ! -d $INSTALLDIR/lsmtool ]; then
    echo Installing LSMTool.
    #
    # Install LSMTool.
    #
    mkdir -p $INSTALLDIR/lsmtool/lib/python$PYTHON_VERSION/site-packages
    export PYTHONPATH=$INSTALLDIR/lsmtool/lib/python$PYTHON_VERSION/site-packages:$PYTHONPATH
    cd $INSTALLDIR/lsmtool && git clone https://github.com/darafferty/LSMTool.git lsmtool
    cd $INSTALLDIR/lsmtool/lsmtool && python setup.py install --prefix=$INSTALLDIR/lsmtool
else
    echo LSMTool already installed.
fi

###############################
# Finish up the installation. #
###############################
echo "Installation directory contents:"
ls ${INSTALLDIR}
#
# init-lofar
#
echo export INSTALLDIR=$INSTALLDIR > $INSTALLDIR/init.sh

echo export PYTHONPATH=\$INSTALLDIR/python-casacore/lib/python${PYTHON_VERSION}/site-packages/:\$INSTALLDIR/DP3/lib/python${PYTHON_VERSION}/site-packages/:\$INSTALLDIR/pyBDSF/lib/python$PYTHON_VERSION/site-packages:\$PYTHONPATH  >> $INSTALLDIR/init.sh
echo export PYTHONPATH=\$INSTALLDIR/python-casacore/lib64/python${PYTHON_VERSION}/site-packages/:\$INSTALLDIR/DP3/lib64/python${PYTHON_VERSION}/site-packages/:\$INSTALLDIR/pyBDSF/lib64/python$PYTHON_VERSION/site-packages:\$PYTHONPATH  >> $INSTALLDIR/init.sh
echo export PATH=\$INSTALLDIR/aoflagger/bin:\$PATH  >> $INSTALLDIR/init.sh
echo export PATH=\$INSTALLDIR/casacore/bin:\$PATH  >> $INSTALLDIR/init.sh
echo export PATH=\$INSTALLDIR/DP3/bin:\$PATH  >> $INSTALLDIR/init.sh
echo export PATH=\$INSTALLDIR/dysco/bin:\$PATH  >> $INSTALLDIR/init.sh
echo export PATH=\$INSTALLDIR/wsclean/bin:\$PATH  >> $INSTALLDIR/init.sh
echo export PATH=\$INSTALLDIR/pyBDSF/bin:\$PATH  >> $INSTALLDIR/init.sh
echo export LD_LIBRARY_PATH=\$INSTALLDIR/aoflagger/lib:\$INSTALLDIR/armadillo/lib64:\$INSTALLDIR/armadillo/lib\$INSTALLDIR/boost/lib:\$INSTALLDIR/casacore/lib:\$INSTALLDIR/cfitsio/lib:\$INSTALLDIR/DP3/lib:\$INSTALLDIR/dysco/lib:\$INSTALLDIR/idg/lib:\$INSTALLDIR/LOFARBeam/lib:\$INSTALLDIR/superlu/lib64:\$INSTALLDIR/superlu/lib:\$INSTALLDIR/wcslib/:\$INSTALLDIR/hdf5/lib:\$LD_LIBRARY_PATH  >> $INSTALLDIR/init.sh

