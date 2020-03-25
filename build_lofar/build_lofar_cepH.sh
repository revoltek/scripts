#!/bin/bash

# General environment settings.
. /sw/profile/init.sh
module switch env env/2019Q4-cuda-gcc-openmpi
# ThOr: This provides python and python-config pointing
# to current python3 and python3-config.
module load pythonversion/3

echo "TODO: Ensure conistent use of one openblas variant (the one cblas and numpy use)"

# ThOr: Maybe not even necessary, at least with cmake builds.
#export LDFLAGS="$(rrz-build-linkpath -R) $LDFLAGS"

export J=64
# ThOr: Use command line argument for the prefix to ease testing.
# Not all variables shall be exported to children!
scriptdir=$(cd $(dirname $0) && pwd)
INSTALLDIR=${1:-$HOME/opt/lofar_200325}
# ThOr: I hate to hardcode things that can be fetched automatically.
# This takes whatever python is linked to python3.
PYTHON_VERSION=$(python3 --version | cut -f 2 -d ' ' | cut -f 1,2 -d .)
PYTHON_VERSION_NODOT=$(echo $PYTHON_VERSION|tr -d .)
export CC=`which gcc`
export CXX=`which g++`
export FC=`which gfortran`
# ThOr: removed C++-specific defines from CFLAGS
export CFLAGS="-march=native -O3 -fno-math-errno -ftree-vectorize"
export FCFLAGS="$CFLAGS"
export CXXFLAGS="-D_GLIB_USE_CXX_ABI=1 -DBOOST_NO_CXX11_SCOPED_ENUMS -march=native -O3 -fno-math-errno -ftree-vectorize"
# ThOr: Those should not be global, should they? Let's add -Wl,--as-needed,
# so that unnecessary libraries are dropped from the build.
# Trying to ensure openblas_openmp here, but that needs to be fixed
# in the packages, recognizing that a dependency is already using
# openblas and then not trying to use a different one:-/
export LDFLAGS="-lgsl -lcblas -lopenblas_openmp -Wl,--as-needed"

# Ensure that we know what include search path variable is really used.
# Only leave CPATH set.
unset INCLUDE
unset C_INCLUDE_PATH

# Settings relevant to the installed software.
# ThOr: I recommend settling on with or without v and add the v on git
# command instead.
export AOFLAGGER_VERSION=v2.14.0
export BOOST_DOT_VERSION=1.67.0
export BOOST_VERSION=1_67_0
export CASACORE_VERSION=latest
export CFITSIO_VERSION=3.47
export DYSCO_VERSION=v1.2.0
export FFTW_VERSION=3.3.8
export HDF5_VERSION=1.8.21
export PYBDSF_VERSION=v1.8.12
export PYTHON_CASACORE_VERSION=3.1.1
export WSCLEAN_VERSION=latest
export WCSLIB_VERSION=latest

# ThOr: These functions emulare what proper environment modules would do.

# Add all given paths to LD_RUN_PATH, staring from
# the vanilla environment setting.
# Also sync LIBRARY_PATH, as some seem to pick that up.
use_prefix()
{
  if test -z "$env_run_path"; then
    env_run_path=$LD_RUN_PATH
    env_path=$PATH
    env_cpath=$CPATH
    #env_ldflags=$LDFLAGS
  fi
  LD_RUN_PATH=$env_run_path
  PATH=$env_path
  CPATH=$env_cpath
  #LDFLAGS=$end_ldflags
  for d in "$@"
  do
    LD_RUN_PATH="$d/lib:$LD_RUN_PATH"
    PATH="$d/bin:$PATH"
    CPATH="$d/include:$CPATH"
    #LDFLAGS="-L$d/lib $LDFLAGS"
  done
}

# ThOr: Run cmake, with standard arguments.
run_cmake()
{
  cmake \
  -DCMAKE_C_COMPILER="$CC" \
  -DCMAKE_CXX_COMPILER="$CXX" \
  -DCMAKE_Fortran_COMPILER="$FC" \
  -DCMAKE_C_FLAGS="$CFLAGS" \
  -DCMAKE_CXX_FLAGS="$CXXFLAGS" \
  -DCMAKE_Fortran_FLAGS="$FCFLAGS" \
  -DCMAKE_SKIP_RPATH=ON \
  "$@"
}

# ThOr: Locate a libray with provided linker
# Arguments: name to grep for (libpng), then compiler/linker arguments.
# Example: find_libfile libpng $(pkg-config --libs libpng) to return
# /some/prefix/lib/libpng16.so.16
find_libfile()
{
  local libname=$1; shift # Grep ldd output for that.
  wrk=$(mktemp -d) &&
  (
    cd $wrk &&
    echo 'int main(){ return 0; }' > test.c &&
    g++ "$@" -o test.bin test.c &&
    ldd test.bin | grep "$libname" | cut -f 3 -d ' '
  )
  rm -rf "$wrk"
}

mkdir -p $INSTALLDIR

#########################################
# Install main LOFAR software packages. #
#########################################

if [ ! -e $INSTALLDIR/casacore/.done ]; then
    #
    # Install CASAcore
    #
    echo Installing CASAcore... &&
    mkdir -p $INSTALLDIR/casacore &&
    cd $INSTALLDIR/casacore &&
    rm -rf bin lib &&
    if test -e src && test -e data; then
        echo "Re-using existing sources."
    else
        rm -rf src build data &&
        mkdir -p build data &&
        git clone https://github.com/casacore/casacore.git src &&
        if [ "${CASACORE_VERSION}" != "latest" ]; then
            ( cd src &&
            git checkout tags/${CASACORE_VERSION} )
        fi
        ( cd data &&
          wget --retry-connrefused \
            ftp://anonymous@ftp.astron.nl/outgoing/Measures/WSRT_Measures.ztar &&
          tar xf WSRT_Measures.ztar )
    fi &&
    cd build &&
    use_prefix $INSTALLDIR/casacore
    run_cmake \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=${INSTALLDIR}/casacore \
      -DDATA_DIR:PATH=${INSTALLDIR}/casacore/data \
      -DENABLE_TABLELOCKING=OFF \
      -DUSE_OPENMP=ON \
      -DUSE_FFTW3=ON \
      -DUSE_HDF5=ON \
      -DBUILD_PYTHON=OFF \
      -DBUILD_PYTHON3=ON \
      ../src/ &&
    make -j $J &&
    make install &&
    # ThOr: Let's see if that is really needed ...
    #ln -s ${INSTALLDIR}/casacore/lib/libcasa_python3.so \
    #      ${INSTALLDIR}/casacore/lib/libcasa_python.so &&
    echo Installed CASAcore. &&
    touch $INSTALLDIR/casacore/.done
else
    echo CASAcore already installed.
fi || exit 1

prefix=$INSTALLDIR/python-casacore
pycasacorelib=$prefix/lib/python$PYTHON_VERSION/site-packages
# ThOr: Subsequent builds might need that, so export it always here.
export PYTHONPATH=$pycasacorelib:$PYTHONPATH

if [ ! -e $prefix/.done ]; then
    #
    # Install-python-casacore
    #
    # Finding libraries is broken, patch the setup to include the previously installed boost and casacore libraries.
    # note the the python executable needs to be the right one
    mkdir -p $prefix &&
    cd $prefix &&
    rm -rf lib bin src &&
    if test -e src.tar.gz; then
        # ThOr: Dunno how to clean, but we can extract archive to avoid re-dowload.
        echo "Using archived sources."
        tar -xf src.tar.gz
    else
        echo "Fresh download."
        git clone https://github.com/casacore/python-casacore src &&
        if [ "$PYTHON_CASACORE_VERSION" != "latest" ]; then
            ( cd src && git checkout tags/v$PYTHON_CASACORE_VERSION )
        fi &&
        tar -czf src.tar.gz src
    fi &&
    # ThOr: Note that ${VAR} serves no other purpose than to delimit
    # the name VAR itself. It does not imply any quoting. "$VAR" is careful,
    # ${VAR} is't. "${VAR}" would be. You should either do use the braces
    # where required only (${VAR}_suffix) or apply them alwas as a matter
    # of style. I myself don't like them where not needed.
    mkdir -p $pycasacorelib &&
    cd src &&
    # ThOr: with the pythonversion/3 module, the plain python does the correct
    # thing. But an alternative would be to call python3 setup.py instead of
    # ./setup.py .
    use_prefix $INSTALLDIR/casacore $prefix &&
    # This magically picks up some RPATH, not sure how, but does neither
    # include confgured casacore path or its own install prefix. But it seems
    # to honour LDFLAGS.
    LDFLAGS="$(rrz-build-linkpath -R) $LDFLAGS" \
        ./setup.py build_ext --swig-opts=-c++  \
        -I$INSTALLDIR/casacore/include -L$INSTALLDIR/casacore/lib \
        install --prefix $prefix &&
    echo Installed Python CASAcore. &&
    touch $prefix/.done
else
    echo Python-CASAcore already installed.
fi || exit 1

prefix=$INSTALLDIR/dysco
if [ ! -e $prefix/.done ]; then
    #
    # Install Dysco
    #
    echo Installing Dysco...
    mkdir -p $prefix &&
    cd       $prefix &&
    rm -rf lib bin build  &&
    if [ -e src ]; then
        echo "Using existing source."
    else
        git clone https://github.com/aroffringa/dysco.git src
    fi &&
    mkdir build &&
    cd    build &&
    use_prefix $INSTALLDIR/casacore $prefix &&
    run_cmake \
        -DCMAKE_INSTALL_PREFIX=$prefix \
        -DCASACORE_ROOT_DIR=$INSTALLDIR/casacore \
        ../src &&
    make -j $J &&
    make install &&
    echo Installed Dysco. &&
    touch $prefix/.done
else
    echo Dysco already installed.
fi || exit 1

prefix=$INSTALLDIR/aoflagger
if [ ! -e $prefix/.done ]; then
    echo "Installing aoflagger..."
    #
    # Install-aoflagger
    #
    mkdir -p $prefix &&
    cd       $prefix &&
    rm -rf build bin lib presrc&&
    if [ -e src ]; then
        echo "Using existing sources."
    else
        git clone git://git.code.sf.net/p/aoflagger/code presrc
        if [ "${AOFLAGGER_VERSION}" != "latest" ]; then
            ( cd presrc && git checkout tags/${AOFLAGGER_VERSION} )
        fi
        if ! [ grep -q Boost::python presrc/src/CMakeLists.txt ]; then
          echo "Patching aoflagger, mainly for Boost::python" 
          (cd presrc && patch -Np1 < "$scriptdir/aoflagger-buildfix.patch" )
        fi &&
        mv presrc src
    fi &&
    mkdir build &&
    cd    build &&
    # ThOr: aoflagger build finds system's libpng, but dependencies use
    # libpng16 from pkgsrc. So that semi-hardcoded fix stays here.
    # But the Boost+Python stuff got fixed via the patch above.
    use_prefix $INSTALLDIR/casacore $prefix &&
    run_cmake \
        -DCASACORE_ROOT_DIR=$INSTALLDIR/casacore \
        -DCMAKE_INSTALL_PREFIX=$prefix \
        -DPNG_LIBRARY_RELEASE:FILEPATH="$(find_libfile libpng $(pkg-config --libs libpng))" \
        -DBUILD_SHARED_LIBS=ON -DPORTABLE=ON \
        ../src &&
    make -j $J &&
    make install &&
    echo "Done with aoflagger." &&
    touch $prefix/.done
else
    echo AOFlagger already installed.
fi || exit 1

prefix=$INSTALLDIR/LOFARBeam
if [ ! -e $prefix/.done ]; then
    #
    # Install the standalone StationResponse libraries.
    #
    echo Installing LOFARBeam...
    mkdir -p $prefix &&
    cd       $prefix &&
    rm -rf build bin lib presrc &&
    if [ -e src ]; then
        echo "Using existing sources."
    else
        git clone https://github.com/lofar-astron/LOFARBeam.git presrc
        if ! [ grep -q Boost:: presrc/CMakeLists.txt ]; then
          echo "Patching LOFARBeam, mainly for Boost::python" 
          (cd presrc && patch -Np1 < "$scriptdir/LOFARBeam-buildfix.patch" )
        fi &&
        mv presrc src
    fi &&
    mkdir build &&
    cd    build &&
    use_prefix $INSTALLDIR/casacore $prefix &&
    run_cmake \
        -DCASACORE_ROOT_DIR=$INSTALLDIR/casacore \
        -DCMAKE_INSTALL_PREFIX=$prefix \
        ../src &&
    make -j $J &&
    make install &&
    echo Done with LOFARBeam. &&
    touch $prefix/.done
fi || exit 1

prefix=$INSTALLDIR/idg
if [ ! -e $prefix/.done ]; then
    #
    # Install the Image Domain Gridder (IDG).
    # 
    echo Installing IDG.
    mkdir -p $prefix &&
    cd       $prefix &&
    rm -rf bin lib presrc build &&    
    if test -d src; then
        echo "Using existing source."
    else
        git clone https://gitlab.com/astron-idg/idg.git presrc &&
        mv presrc src
    fi &&
    mkdir build &&
    cd    build &&
    use_prefix $prefix &&
    # ThOr: The FindMKL included (two copies!) tries to locate libmkl_rt and
    # sets paths from that one. Uses MKL_LIB as environment hint.
    module load mkl/2020.0.166 &&
    MKL_LIB=$LD_RUN_PATH run_cmake \
        -DCMAKE_INSTALL_PREFIX=$prefix \
        -DBUILD_WITH_MKL=ON \
        -DBUILD_LIB_CPU=ON \
        -DBUILD_LIB_CUDA=ON \
        -DCUDA_INCLUDE_DIR=$CUDA_HOME/include \
        -DCUDA_FFT_LIBRARY=$(find_libfile libcufft.so -lcufft) \
        -DCUDA_NVTX_LIBRARY=$(find_libfile libnvToolsExt.so -lnvToolsExt) \
        ../src &&
    make -j $J &&
    make install &&
    module unload mkl &&
    echo "Done with IDG." &&
    touch $prefix/.done
else
    echo IDG already installed.
fi || exit 1

prefix=$INSTALLDIR/DP3
if [ ! -e $prefix/.done ]; then
    echo Installing DP3...
    #
    # Install DP3.
    #
    mkdir -p $prefix
    cd       $prefix
    rm -rf build bin lib presrc &&
    if test -e src; then
        echo "Using existing sources."
    else
        git clone https://github.com/lofar-astron/DP3.git presrc &&
        if ! [ grep -q Boost::python presrc/src/CMakeLists.txt ]; then
            echo "Patching build, mainly for Boost::python" 
            (cd presrc && patch -Np1 < "$scriptdir/DP3-buildfix.patch" )
        fi &&
        mv presrc src
    fi &&
    mkdir build &&
    cd    build &&
    use_prefix \
      $INSTALLDIR/casacore \
      $INSTALLDIR/idg \
      $INSTALLDIR/LOFARBeam \
      $INSTALLDIR/aoflagger \
      $prefix &&
    # ThOr: It's braindead that I still have to provide things like
    # AOFLAGGER_INCLUDE_DIR and AOFLAGGER_LIB! You should consider
    # just locating the library via LD_RUN_PATH or LIBRARY_PATH, set
    # by use_prefix. Or install an aoflagger-config script, like NetCDF.
    # Or a pkg-config file. Anything standard!
    # Really … having to specify the full path to the .so file?! Why do
    # you love hurting yourself (or your fellow astrophysicists) that much?
    run_cmake \
        -DCMAKE_INSTALL_PREFIX:PATH=$prefix \
        -DCASACORE_ROOT_DIR=$INSTALLDIR/casacore \
        -DIDGAPI_INCLUDE_DIRS:PATH=$INSTALLDIR/idg/include \
        -DIDGAPI_LIBRARIES:PATH=$INSTALLDIR/idg/lib/libidg-api.so \
        -DLOFAR_STATION_RESPONSE_DIR:PATH=$INSTALLDIR/LOFARBeam/include \
        -DAOFLAGGER_INCLUDE_DIR=$INSTALLDIR/aoflagger/include \
        -DAOFLAGGER_LIB=$INSTALLDIR/aoflagger/lib/libaoflagger.so \
        ../src &&
    make -j $J &&
    make install &&
    echo Installed DP3. &&
    touch $prefix/.done
else
    echo DP3 already installed.
fi

prefix=$INSTALLDIR/wsclean
if [ ! -e $prefix/.done ]; then
    #
    # Install WSClean
    #
    echo Installing WSClean.
    mkdir -p $prefix &&
    cd       $prefix &&
    rm -rf build bin lib presrc &&
    if test -e src; then
        echo "Using existing sources."
    else
        if [ "$WSCLEAN_VERSION" != "latest" ]; then
            # ThOr: Untested ...
            wget http://downloads.sourceforge.net/project/wsclean/wsclean-${WSCLEAN_VERSION}/wsclean-${WSCLEAN_VERSION}.tar.bz2 &&
            tar -xf wsclean-${WSCLEAN_VERSION}.tar.bz2 &&
            mkdir presrc &&
            mv wsclean-${WSCLEAN_VERSION} presrc/wsclean
        else
            git clone git://git.code.sf.net/p/wsclean/code presrc &&
            mv presrc src
        fi
    fi &&
    mkdir build &&
    cd    build &&
    use_prefix $INSTALLDIR/casacore $INSTALLDIR/LOFARBeam $INSTALLDIR/idg &&
    LDFLAGS="$(rrz-build-linkpath -R) $LDFLAGS" cmake \
        -DCMAKE_INSTALL_PREFIX=$prefix \
        -DCASACORE_ROOT_DIR=$INSTALLDIR/casacore \
        -DIDGAPI_INCLUDE_DIRS=$INSTALLDIR/idg/include \
        -DIDGAPI_LIBRARIES=$INSTALLDIR/idg/lib/libidg-api.so \
        -DLOFAR_STATION_RESPONSE_INCLUDE_DIR=$INSTALLDIR/LOFARBeam/include \
        -DLOFAR_STATION_RESPONSE_LIB=$INSTALLDIR/LOFARBeam/lib/libstationresponse.so \
        ../src/wsclean &&
    make -j $J &&
    make install &&
    echo "Done with WSClean." &&
    touch $prefix/.done
else
    echo WSClean already installed.
fi

prefix=$INSTALLDIR/pyBDSF &&
pybdsflib=$prefix/lib/python$PYTHON_VERSION/site-packages
export PYTHONPATH=$pybdsflib:$PYTHONPATH &&
if [ ! -e $prefix/.done ]; then
    echo Installing pyBDSF.
    #
    # Install pyBDSF.
    #
    mkdir -p $prefix &&
    cd       $prefix &&
    rm -rf build bin lib presrc &&
    if test -e src; then
        echo "Using existing sources."
    else
        git clone https://github.com/lofar-astron/PyBDSF presrc &&
        mv presrc src
    fi &&
    mkdir -p $pybdsflib &&
    cd src &&
    LDFLAGS="$(rrz-build-linkpath -R) $LDFLAGS" \
        NPY_DISTUTILS_APPEND_FLAGS=1 \
        python setup.py install --prefix=$prefix &&
    echo "Done with pyBDSF." &&
    touch $prefix/.done
else
    echo pyBDSF already installed.
fi || exit 1

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
echo export PATH=\$INSTALLDIR/aoflagger/bin:\$PATH  >> $INSTALLDIR/init.sh
echo export PATH=\$INSTALLDIR/casacore/bin:\$PATH  >> $INSTALLDIR/init.sh
echo export PATH=\$INSTALLDIR/DP3/bin:\$PATH  >> $INSTALLDIR/init.sh
echo export PATH=\$INSTALLDIR/dysco/bin:\$PATH  >> $INSTALLDIR/init.sh
echo export PATH=\$INSTALLDIR/wsclean/bin:\$PATH  >> $INSTALLDIR/init.sh
echo export PATH=\$INSTALLDIR/pyBDSF/bin:\$PATH  >> $INSTALLDIR/init.sh
echo export LD_LIBRARY_PATH=\$INSTALLDIR/aoflagger/lib:\$INSTALLDIR/casacore/lib:\$INSTALLDIR/DP3/lib:\$INSTALLDIR/dysco/lib:\$INSTALLDIR/idg/lib:\$INSTALLDIR/LOFARBeam/lib:\$LD_LIBRARY_PATH  >> $INSTALLDIR/init.sh
