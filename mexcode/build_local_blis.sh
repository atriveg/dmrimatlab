#!/bin/bash

TOTAL_MEM_KB=$(cat /proc/meminfo | grep MemTotal | sed -E 's/MemTotal:[ ]+([0-9]+)[ ]+kB/\1/' )
MAX_THREADS_BY_RAM=$(expr ${TOTAL_MEM_KB} \/ 2000000)
NPROCS=$(nproc --all)
if [ "${MAX_THREADS_BY_RAM}" -gt "0" ]; then
    if [ "${MAX_THREADS_BY_RAM}" -gt "${NPROCS}" ]; then
        BTHREADS="${NPROCS}"
    else
        BTHREADS="${MAX_THREADS_BY_RAM}"
    fi
else
    BTHREADS=1
fi

echo "Using ${BTHREADS} threads to compile source code"

CURRDIR=$(pwd)
SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
INSTALLDIR="${SCRIPTDIR}/BLIS"
TMPDIR="/var/tmp/blis_build_$(date "+%y-%m-%d-%H-%M-%S")"

echo "Script dir: ${SCRIPTDIR}"
echo "Install dir: ${INSTALLDIR}"
echo "Temp dir: ${TMPDIR}"

mkdir -p "${INSTALLDIR}"

if [ "$?" -ne 0 ]
then
    echo "Could not create install dir. Abort..."
    exit 1
fi

mkdir -p "${TMPDIR}"

if [ "$?" -ne 0 ]
then
    echo "Could not create tmp dir. Abort..."
    exit 2
fi

cd "${TMPDIR}"

if [ "$?" -ne 0 ]
then
    echo "Could not cd to tmp dir. Abort..."
    exit 3
fi

git clone "https://github.com/flame/blis.git"

if [ "$?" -ne 0 ]
then
    echo "Could not clone FLAME/BLIS repo. Abort..."
    exit 4
fi

cd blis && git checkout 2.1

if [ "$?" -ne 0 ]
then
    echo "Could not checkout BLIS 2.1. Abort..."
    exit 5
fi

export CFLAGS="-fexceptions -fPIC -fno-omit-frame-pointer -pthread -fwrapv -O3 -DNDEBUG"

./configure --prefix="${INSTALLDIR}" \
    --disable-static \
    --enable-shared \
    --disable-lapack-compat \
    --enable-threading=single \
    --disable-cblas \
    --enable-amd-frame-tweaks \
    --blas-int-size=32 \
    auto

if [ "$?" -ne 0 ]
then
    echo "Could not configure FLAME/BLIS. Abort..."
    exit 6
fi

make -j "${BTHREADS}"

if [ "$?" -ne 0 ]
then
    echo "Could not build FLAME/BLIS. Abort..."
    exit 7
fi

make install

if [ "$?" -ne 0 ]
then
    echo "Could not install FLAME/BLIS. Abort..."
    exit 8
fi

cd "${TMPDIR}"

if [ "$?" -ne 0 ]
then
    echo "Could not cd to tmp dir. Abort..."
    exit 9
fi

git clone "https://github.com/Reference-LAPACK/lapack.git"

if [ "$?" -ne 0 ]
then
    echo "Could not clone LAPACK repo. Abort..."
    exit 10
fi

cd lapack

if [ "$?" -ne 0 ]
then
    echo "Could not cd to LAPACK folder. Abort..."
    exit 11
fi

mkdir -p build

if [ "$?" -ne 0 ]
then
    echo "Could not create build folder. Abort..."
    exit 12
fi

cd build

if [ "$?" -ne 0 ]
then
    echo "Could not cd to build folder. Abort..."
    exit 13
fi

cmake .. \
    -G "Unix Makefiles" \
    -DBUILD_COMPLEX=OFF \
    -DBUILD_SHARED_LIBS=ON \
    -DLAPACKE=OFF \
    -DCBLAS=OFF \
    -DBLAS_LIBRARIES="${INSTALLDIR}/lib/libblis.so" \
    -DCMAKE_INSTALL_PREFIX="${INSTALLDIR}" \
    -DCMAKE_BUILD_TYPE=Release

if [ "$?" -ne 0 ]
then
    echo "Could not configure LAPACK. Abort..."
    exit 14
fi

make -j "${BTHREADS}"

if [ "$?" -ne 0 ]
then
    echo "Could not build LAPACK. Abort..."
    exit 15
fi

make install

if [ "$?" -ne 0 ]
then
    echo "Could not installº LAPACK. Abort..."
    exit 16
fi

rm -fR "${TMPDIR}"

if [ "$?" -ne 0 ]
then
    echo "Could not clean up. Abort..."
    exit 17
fi

cd "${CURRDIR}"
exit 0
