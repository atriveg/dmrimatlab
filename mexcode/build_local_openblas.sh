#!/bin/bash

DIRSUFFIX="-$2"
LIBSUFFIX="_$2"

SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
INSTALLDIR="${SCRIPTDIR}/openblas$DIRSUFFIX"
TMPDIR="/var/tmp/openblas_build_$(date "+%y-%m-%d-%H-%M-%S")"

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

git clone "https://github.com/OpenMathLib/OpenBLAS.git"

if [ "$?" -ne 0 ]
then
    echo "Could not clone OpenBLAS repo. Abort..."
    exit 4
fi

cd OpenBLAS

if [ "$?" -ne 0 ]
then
    echo "Could not cd to OpenBLAS code folder. Abort..."
    exit 5
fi


export CXXFLAGS="-fexceptions -fPIC -fno-omit-frame-pointer -pthread -fwrapv -O3 -DNDEBUG"
export CFLAGS="-fexceptions -fPIC -fno-omit-frame-pointer -pthread -fwrapv -O3 -DNDEBUG"
export FFLAGS="-O3"
case "$1" in
  "single-thread")
    make -j DYNAMIC_ARCH=0 BINARY=64 INTERFACE=64 NO_AFFINITY=1 NO_WARMUP=1 USE_OPENMP=0 USE_THREAD=0 USE_LOCKING=1 LIBNAMESUFFIX="${LIBSUFFIX}"
    ;;
  *)
    make -j DYNAMIC_ARCH=0 BINARY=64 INTERFACE=64 NO_AFFINITY=1 NO_WARMUP=1 USE_OPENMP=0 USE_THREAD=0 USE_LOCKING=1 LIBNAMESUFFIX="${LIBSUFFIX}"
    ;;
esac


if [ "$?" -ne 0 ]
then
    echo "Could not build OpenBLAS. Abort..."
    exit 6
fi

make PREFIX="${INSTALLDIR}" LIBNAMESUFFIX="${LIBSUFFIX}" install

if [ "$?" -ne 0 ]
then
    echo "Could not locally install OpenBLAS. Abort..."
    exit 7
fi

cd "${INSTALLDIR}"

if [ "$?" -ne 0 ]
then
    echo "Could not locally cd to install dir. Abort..."
    exit 8
fi

rm -fR "${TMPDIR}" "${INSTALLDIR}/bin" "${INSTALLDIR}/lib/"*.a "${INSTALLDIR}/lib/"*.a.*

if [ "$?" -ne 0 ]
then
    echo "Could not clean up. Abort..."
    exit 9
fi

exit 0
