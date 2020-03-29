#!/bin/bash

# A heavily changed LOFAR build script based on the one Francesco gave me.
# I put my changes into the public domain, hoping that they are useful.
# -- Dr. Thomas Orgis <thomas.orgis@uni-hamburg.de>

# TODO:
# - some syntax cleanup, "$var" instead of $var, most of the times
# - maybe more structre to reduce the repetitive identical
#   fetch/patch/build code

mandatory_vars="\
lofar_prefix      Installation prefix to use. Packages are put into sub-directories of that.
casacore_version  Version of CASAcore.
pycasa_version    Version of CASAcore Python bindings.
dysco_version     Version of Dysco.
pybdsf_version    Version of PyBDSF.
wsclean_version   Version of WSClean.
CC                C compiler command to use.
CXX               C++ compiler command to use.
FC                Fortran compiler command to use.
CFLAGS            C compiler flags (maybe empty, but set).
CXXFLAGS          C++ compiler flags (maybe empty, but set).
FCFLAGS           Fortran compiler flags (maybe empty, but set).
BLAS_LIBS         Linker flags to get BLAS (e.g. -lopenblas).
CBLAS_LIBS        Linker flags to get CBLAS (e.g. -lcblas)
LAPACK_LIBS       Linker flags to get LAPACK (e.g. -lopenblas).
"

optional_vars="\
MAKE_JOBS    Number of parallel make jobs to run.
LIBRARY_PATH Path list for link-time library search.
LD_RUN_PATH  Path for run-time library search path fixed at build time.
CPATH        Path list for C/C++ include file search.
idg_mkl      Passed on to idg BUILD_WITH_MKL
"

config=$0.config

if ! test -e "$config"; then
  cat <<EOT

HALTIS - Hamburg LOFAR Toolbox Installation Script

Please create/link an appropriate config script named

  $config

. It should load the correct environment for the site, possibly
using environment modules or just setting variables for paths
and compilers directly. It also has to set up PATH so that the
desired python interpreter is directly available as the command
python (symlinks/wrappers to call python3 and python3-config should
be enough). Furthermore, GNU Make is expected as make command,
and the various dependencies are expected to be present, either
in the base OS or in some additional environment you prepared.

Uppercase variables may be exported, lowercase one should not.

Mandatory variables:

$mandatory_vars

Optional variables:

$optional_vars

Essential variables are exportet by the build script. You
should only export variables in addition to setting them if you
really want their value to affect the builds via the environment
(standard variables like PATH are expected to be exported already).

EOT
  exit 1
fi


echo "Loading $config ..."
. "$config" || exit 1

while read var help
do
  test -n "$var" || continue
  echo "Checking mandatory variable $var ..."
  is_set=$(eval printf %s "\${$var+set}")
  val=$(eval printf %s "\$$var")
  if test "$is_set" != set; then
    echo "Mandatory variable $var is not set." >&2
    exit 1
  fi
  if test -z "$val"; then
    case "$var" in
      *FLAGS) true ;;
      *)
        echo "Variable $var needs to be non-empty." >&2
        exit 1
      ;;
    esac
  fi
done <<< "$mandatory_vars"

case "$lofar_prefix" in
/*)
  echo "Good. Given prefix could be an absolute path."
;;
*)
  echo "Please specify absolute LOFAR_PATH!" >&2
  exit 1
;;
esac

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

# Write the command to be executed to a shell script for debugging
# and run it. Some minimal quoting is attempted.
store_and_run()
{
  local script=$1; shift
  export > "$script.env"
  cat <<EOT > "$script"
#!$SHELL
#
# This reconstructs the central configuraton/build command.
# To debug interactively, consider loading $script.env, too.
#
EOT
  for n in "$@"
  do
    case "$n" in
      *\'*)
        echo "Really a command/argument with verbatim ' in it? $n" >&2
        # This might be broken.
        printf '%s \\\n' "\"$n\""
      ;;
      *[^[:alnum:]_-/:+,.=]*|'')
        printf '%s \\\n' "'$n'"
      ;;
      *) # hopefully safe
        printf '%s \\\n' "$n"
      ;;
    esac
  done >> "$script"
  echo >> "$script"
  "$@"
}

# Run cmake, with standard arguments.
# One of those is PORTABLE=ON to avoid messing
# with our provided compiler flags.
run_cmake()
{
  store_and_run "$prefix/configure.sh" cmake \
    -DCMAKE_C_COMPILER="$CC" \
    -DCMAKE_CXX_COMPILER="$CXX" \
    -DCMAKE_Fortran_COMPILER="$FC" \
    -DCMAKE_C_FLAGS="$CFLAGS" \
    -DCMAKE_CXX_FLAGS="$CXXFLAGS" \
    -DCMAKE_Fortran_FLAGS="$FCFLAGS" \
    -DCMAKE_SKIP_RPATH=ON \
    -DPORTABLE=ON \
    "$@"
}

# Take full version, like 1.2.3 and return git tag
# name. Some have v1.2.3, some use v1.2.3 but v1.2
# for 1.2.0.
git_version_tag()
{
  local fullvers=$1                        # 1.2.3
  local mmvers=$(cut -f 1,2 -d . <<< "$1") # 1.2
  local mvers=$(cut -f 1 -d . <<< "$1")    # 1
  for v in "$fullvers" "$mmvers" "$mvers"
  do
    if git tag | grep -q "^v$v\$"; then
      printf "v%s\n" "$v"
      return
    fi
  done
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

begin_pkg()
{
  pkg=$1
  echo "$(date) $pkg"
  prefix="$lofar_prefix/$pkg"
  mkdir -p "$prefix"
}

mkdir -p "$lofar_prefix" &&
env_python_path=$PYTHON_PATH &&
scriptdir=$(cd $(dirname $0) && pwd) &&
config="$scriptdir/$(basename $config)" &&
pyver=$(python --version | cut -f 2 -d ' ' | cut -f 1,2 -d .) &&
test -n "$pyver" &&
case "$pyver" in
  2.*)
    build_python=ON
    build_python3=OFF
  ;;
  3.*)
    build_python=OFF
    build_python3=ON
  ;;
  *)
    echo "There's another Python? Please check if I am OK with that."
    exit 1
  ;;
esac
# Ensure that we know what include search path variable is really used.
# Only leave CPATH set.
unset INCLUDE &&
unset C_INCLUDE_PATH &&

#########################################
# Install main LOFAR software packages. #
#########################################

begin_pkg casacore
if [ ! -e $prefix/.done ]; then
    #
    # Install CASAcore
    #
    echo Installing CASAcore... &&
    cd $prefix &&
    rm -rf bin build include lib presrc &&
    if test -e src && test -e data; then
        echo "Re-using existing sources."
    else
        rm -rf src data &&
        git clone https://github.com/casacore/casacore.git presrc &&
        if [ "${casacore_version}" != "latest" ]; then
            ( cd presrc &&
            git checkout tags/v${casacore_version} )
        fi &&
        echo "Patching build (for BLAS, Boost, Python, usually) ..." && 
        (cd presrc && patch -Np1 < "$scriptdir/casacore-buildfix.patch" ) &&
        mkdir data &&
        ( cd data &&
          wget --retry-connrefused \
            ftp://anonymous@ftp.astron.nl/outgoing/Measures/WSRT_Measures.ztar &&
          tar xf WSRT_Measures.ztar ) &&
        mv presrc src
    fi &&
    mkdir build &&
    cd    build &&
    use_prefix $prefix &&
    run_cmake \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=$lofar_prefix/casacore \
      -DDATA_DIR:PATH=$lofar_prefix/casacore/data \
      -DENABLE_TABLELOCKING=OFF \
      -DUSE_OPENMP=ON \
      -DUSE_FFTW3=ON \
      -DUSE_HDF5=ON \
      -DBUILD_PYTHON=$build_python \
      -DBUILD_PYTHON3=$build_python3 \
      ../src/ &&
    make -j $J &&
    make install &&
    echo Installed CASAcore. &&
    touch $prefix/.done
else
    echo CASAcore already installed.
fi >> "$prefix/build.log" 2>&1 || exit 1

begin_pkg python-casacore
pycasacorelib="$prefix/lib/python$pyver/site-packages"
# ThOr: Subsequent builds might need that, so export it always here.
export PYTHONPATH=$pycasacorelib:$PYTHONPATH
if [ ! -e $prefix/.done ]; then
    #
    # Install-python-casacore
    #
    # Finding libraries is broken, patch the setup to include the previously installed boost and casacore libraries.
    # note the the python executable needs to be the right one
    cd $prefix &&
    rm -rf lib bin src pre &&
    if test -e src.tar.gz; then
        # ThOr: Dunno how to clean, but we can extract archive to avoid re-dowload.
        echo "Using archived sources." &&
        tar -xf src.tar.gz
    else
        echo "Fresh download."
        git clone https://github.com/casacore/python-casacore presrc &&
        if [ "$pycasa_version" != "latest" ]; then
            ( cd presrc && git checkout tags/v$pycasa_version )
        fi &&
        mv presrc src &&
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
    use_prefix $lofar_prefix/casacore $prefix &&
    # This magically picks up some RPATH, not sure how, but does neither
    # include confgured casacore path or its own install prefix. But it seems
    # to honour LDFLAGS.
    LDFLAGS="$(rrz-build-linkpath -R) $LDFLAGS" \
        store_and_run "$prefix/install.sh" \
        python setup.py build_ext --swig-opts=-c++  \
        -I$lofar_prefix/casacore/include -L$lofar_prefix/casacore/lib \
        install --prefix $prefix &&
    echo Installed Python CASAcore. &&
    touch $prefix/.done
else
    echo Python-CASAcore already installed.
fi >> "$prefix/build.log" 2>&1 || exit 1

begin_pkg dysco
if [ ! -e $prefix/.done ]; then
    #
    # Install Dysco
    #
    echo Installing Dysco...
    cd $prefix &&
    rm -rf lib bin build presrc &&
    if [ -e src ]; then
        echo "Using existing source."
    else
        git clone https://github.com/aroffringa/dysco.git presrc &&
        if [ "$dysco_version" != "latest" ]; then
            ( cd presrc &&
              git checkout tags/"$(git_version_tag "$dysco_version")" )
        fi &&
        echo "Patching build (for BLAS, Boost, Python, usually) ..." && 
        (cd presrc && patch -Np1 < "$scriptdir/dysco-buildfix.patch" ) &&
        mv presrc src
    fi &&
    mkdir build &&
    cd    build &&
    use_prefix $lofar_prefix/casacore $prefix &&
    run_cmake \
        -DCMAKE_INSTALL_PREFIX=$prefix \
        -DCASACORE_ROOT_DIR=$lofar_prefix/casacore \
        ../src &&
    make -j $J &&
    make install &&
    echo Installed Dysco. &&
    touch $prefix/.done
else
    echo Dysco already installed.
fi >> "$prefix/build.log" 2>&1 || exit 1

begin_pkg aoflagger
if [ ! -e $prefix/.done ]; then
    echo "Installing aoflagger..."
    #
    # Install-aoflagger
    #
    cd $prefix &&
    rm -rf build bin lib presrc&&
    if [ -e src ]; then
        echo "Using existing sources."
    else
        git clone git://git.code.sf.net/p/aoflagger/code presrc
        if [ "$aoflagger_version" != "latest" ]; then
            ( cd presrc && git checkout tags/v$aoflagger_version )
        fi &&
        echo "Patching build (for BLAS, Boost, Python, usually) ..." && 
        (cd presrc && patch -Np1 < "$scriptdir/aoflagger-buildfix.patch" ) &&
        mv presrc src
    fi &&
    mkdir build &&
    cd    build &&
    # ThOr: aoflagger build finds system's libpng, but dependencies use
    # libpng16 from pkgsrc. So that semi-hardcoded fix stays here.
    # But the Boost+Python stuff got fixed via the patch above.
    use_prefix $lofar_prefix/casacore $prefix &&
    run_cmake \
        -DCASACORE_ROOT_DIR=$lofar_prefix/casacore \
        -DCMAKE_INSTALL_PREFIX=$prefix \
        -DPNG_LIBRARY_RELEASE:FILEPATH="$(find_libfile libpng $(pkg-config --libs libpng))" \
        -DBUILD_SHARED_LIBS=ON \
        ../src &&
    make -j $J &&
    make install &&
    echo "Done with aoflagger." &&
    touch $prefix/.done
else
    echo AOFlagger already installed.
fi >> "$prefix/build.log" 2>&1 || exit 1

begin_pkg LOFARBeam
if [ ! -e $prefix/.done ]; then
    #
    # Install the standalone StationResponse libraries.
    #
    echo Installing LOFARBeam...
    cd $prefix &&
    rm -rf build bin lib presrc &&
    if [ -e src ]; then
        echo "Using existing sources."
    else
        git clone https://github.com/lofar-astron/LOFARBeam.git presrc
        echo "Patching build (for BLAS, Boost, Python, usually) ..." && 
        (cd presrc && patch -Np1 < "$scriptdir/LOFARBeam-buildfix.patch" ) &&
        mv presrc src
    fi &&
    mkdir build &&
    cd    build &&
    use_prefix $lofar_prefix/casacore $prefix &&
    run_cmake \
        -DCASACORE_ROOT_DIR=$lofar_prefix/casacore \
        -DCMAKE_INSTALL_PREFIX=$prefix \
        ../src &&
    make -j $J &&
    make install &&
    echo Done with LOFARBeam. &&
    touch $prefix/.done
fi >> "$prefix/build.log" 2>&1 || exit 1

begin_pkg idg
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
        echo "Patching build (for BLAS, Boost, Python, usually) ..." && 
        (cd presrc && patch -Np1 < "$scriptdir/idg-buildfix.patch" ) &&
        mv presrc src
    fi &&
    mkdir build &&
    cd    build &&
    use_prefix $prefix &&
    # ThOr: The FindMKL included (two copies!) tries to locate libmkl_rt and
    # sets paths from that one. Uses MKL_LIB as environment hint.
    MKL_LIB=$LD_RUN_PATH run_cmake \
        -DCMAKE_INSTALL_PREFIX=$prefix \
        -DBUILD_WITH_MKL=${idg_mkl:-OFF} \
        -DBUILD_LIB_CPU=ON \
        -DBUILD_LIB_CUDA=ON \
        -DCUDA_INCLUDE_DIR=$CUDA_HOME/include \
        -DCUDA_FFT_LIBRARY=$(find_libfile libcufft.so -lcufft) \
        -DCUDA_NVTX_LIBRARY=$(find_libfile libnvToolsExt.so -lnvToolsExt) \
        ../src &&
    make -j $J &&
    make install &&
    echo "Done with IDG." &&
    touch $prefix/.done
else
    echo IDG already installed.
fi >> "$prefix/build.log" 2>&1 || exit 1

begin_pkg DP3
export "PYTHONPATH=$prefix/lib/python$pyver:$PYTHONPATH" &&
if [ ! -e $prefix/.done ]; then
    echo Installing DP3...
    #
    # Install DP3.
    #
    cd $prefix &&
    rm -rf build bin lib presrc &&
    if test -e src; then
        echo "Using existing sources."
    else
        git clone https://github.com/lofar-astron/DP3.git presrc &&
        echo "Patching build (for BLAS, Boost, Python, usually) ..." && 
        (cd presrc && patch -Np1 < "$scriptdir/DP3-buildfix.patch" ) &&
        mv presrc src
    fi &&
    mkdir build &&
    cd    build &&
    use_prefix \
      $lofar_prefix/casacore \
      $lofar_prefix/idg \
      $lofar_prefix/LOFARBeam \
      $lofar_prefix/aoflagger \
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
        -DCASACORE_ROOT_DIR=$lofar_prefix/casacore \
        -DIDGAPI_INCLUDE_DIRS:PATH=$lofar_prefix/idg/include \
        -DIDGAPI_LIBRARIES:PATH=$lofar_prefix/idg/lib/libidg-api.so \
        -DLOFAR_STATION_RESPONSE_DIR:PATH=$lofar_prefix/LOFARBeam/include \
        -DLOFAR_STATION_RESPONSE_LIB=$lofar_prefix/LOFARBeam/lib/libstationresponse.so \
        -DAOFLAGGER_INCLUDE_DIR=$lofar_prefix/aoflagger/include \
        -DAOFLAGGER_LIB=$lofar_prefix/aoflagger/lib/libaoflagger.so \
        ../src &&
    make -j $J &&
    make install &&
    echo Installed DP3. &&
    touch $prefix/.done
else
    echo DP3 already installed.
fi >> "$prefix/build.log" 2>&1 || exit 1

begin_pkg wsclean
if [ ! -e $prefix/.done ]; then
    #
    # Install WSClean
    #
    # Sometimes there are FFT to update by hand in CMakeCache.txt
    # ThOr: Please no! ;-)
    echo Installing WSClean.
    cd $prefix &&
    rm -rf build bin lib presrc &&
    if test -e src; then
        echo "Using existing sources."
    else
        if [ "$wsclean_version" != "latest" ]; then
            # ThOr: Untested ...
            wget http://downloads.sourceforge.net/project/wsclean/wsclean-${wsclean_version}/wsclean-${wsclean_version}.tar.bz2 &&
            tar -xf wsclean-${wsclean_version}.tar.bz2 &&
            mkdir presrc &&
            mv wsclean-${wsclean_version} presrc/wsclean
        else
            git clone git://git.code.sf.net/p/wsclean/code presrc
        fi &&
        echo "Patching build (for BLAS, Boost, Python, usually) ..." &&
        (cd presrc && patch -Np1 < "$scriptdir/$pkg-buildfix.patch" ) &&
        mv presrc src
    fi &&
    mkdir build &&
    cd    build &&
    use_prefix $lofar_prefix/casacore $lofar_prefix/LOFARBeam $lofar_prefix/idg &&
    LDFLAGS="$(rrz-build-linkpath -R) $LDFLAGS" run_cmake \
        -DCMAKE_INSTALL_PREFIX=$prefix \
        -DCASACORE_ROOT_DIR=$lofar_prefix/casacore \
        -DIDGAPI_INCLUDE_DIRS=$lofar_prefix/idg/include \
        -DIDGAPI_LIBRARIES=$lofar_prefix/idg/lib/libidg-api.so \
        -DLOFAR_STATION_RESPONSE_INCLUDE_DIR=$lofar_prefix/LOFARBeam/include \
        -DLOFAR_STATION_RESPONSE_LIB=$lofar_prefix/LOFARBeam/lib/libstationresponse.so \
        ../src/wsclean &&
    make -j $J &&
    make install &&
    echo "Done with WSClean." &&
    touch $prefix/.done
else
    echo WSClean already installed.
fi >> "$prefix/build.log" 2>&1 || exit 1

begin_pkg pyBDSF
pybdsflib="$prefix/lib/python$pyver/site-packages"
export PYTHONPATH=$pybdsflib:$PYTHONPATH &&
if [ ! -e $prefix/.done ]; then
    echo Installing pyBDSF.
    #
    # Install pyBDSF.
    #
    cd $prefix &&
    rm -rf build bin lib presrc &&
    if test -e src; then
        echo "Using existing sources."
    else
        git clone https://github.com/lofar-astron/PyBDSF presrc &&
        if [ "$pybdsf_version" != "latest" ]; then
            ( cd presrc && git checkout "tags/v$pybdsf_version" )
        fi &&
        mv presrc src
    fi &&
    mkdir -p $pybdsflib &&
    cd src &&
    LDFLAGS="$(rrz-build-linkpath -R) $LDFLAGS" \
        NPY_DISTUTILS_APPEND_FLAGS=1 \
        store_and_run "$prefix/install.sh" \
        python setup.py install --prefix=$prefix &&
    echo "Done with pyBDSF." &&
    touch $prefix/.done
else
    echo pyBDSF already installed.
fi >> "$prefix/build.log" 2>&1 || exit 1

begin_pkg lsmtool
lsmlib="$prefix/lib/python$pyver/site-packages" &&
export PYTHONPATH=$lsmlib:$PYTHONPATH &&
if [ ! -e $prefix/.done ]; then
    echo Installing LSMTool.
    #
    # Install LSMTool.
    #
    cd $prefix &&
    rm -rf build bin lib presrc &&
    if test -e src; then
        echo "Using existing sources."
    else
        git clone https://github.com/darafferty/LSMTool.git presrc &&
        mv presrc src
    fi &&
    mkdir -p $lsmlib &&    
    cd src &&
    store_and_run "$prefix/install.sh" \
      python setup.py install --prefix=$prefix &&
    echo "Done with LSMTool." &&
    touch $prefix/.done
else
    echo LSMTool already installed.
fi >> "$prefix/build.log" 2>&1 || exit 1

echo "$(date) done with packages"

###############################
# Finish up the installation. #
###############################

echo
echo "Installation directory contents:"
ls "$lofar_prefix"

# Use the same environment used to building.

print_init()
{
  cat <<EOT
# LOFAR Tools built by $USER@$(hostname) on $(date)
#
# Sourcing this script shall re-create the build environment and
# add things from the LOFAR prefix to appropriate path variables.
# It does _not_ set LD_LIBRARY_PATH (unless the build environment
# does), as a proper build should use RPATH in the binaries.
# But there is LD_RUN_PATH / LIBRARY_PATH set for link-time usage.
# You can set LD_LIBRARY_PATH=\$LD_RUN_PATH as a hack if some binary
# fails to find its libs. But this would be a bug in the build
# script that should be fixed properly.
#
# Step 1: The build config with all variables used.
#
EOT
  cat $config
cat <<EOT
#
# Step 2: Paths to LOFAR tools in the prefix.
#
EOT
  for p in lsmtool pyBDSF wsclean dysco DP3\
         idg aoflagger python-casacore casacore
  do
    echo "# $p"
    pp="$lofar_prefix/$p"
    test -d "$pp/bin" &&
    echo "export PATH=\"\$lofar_prefix/$p/bin:\$PATH\""
    test -d "$pp/lib" &&
      cat <<EOT
export LD_RUN_PATH=\"\$lofar_prefix/$p/lib:\$LD_RUN_PATH\"
export LIBRARY_PATH=\"\$lofar_prefix/$p/lib:\$LIBRARY_PATH\"
EOT
    test -d "$pp/include" &&
      echo "export CPATH=\"\$lofar_prefix/$p/include:\$CPATH\""
    pylib="$p/lib/python$pyver/site-packages"
    test -d "$lofar_prefix/$pylib" &&
      	echo "export PYTHONPATH=\"\$lofar_prefix/$pylib:\$PYTHONPATH\""
  done
}

print_init > "$lofar_prefix/init.sh" &&
echo "Finished installation. You can source $lofar_prefix/init.sh and run."
