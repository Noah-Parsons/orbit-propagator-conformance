#!/bin/sh
# Compiles Basilisk's lambertSolver module source directly (no SWIG, no Conan) into two small
# command-line programs: lambert_orig from the unmodified source and lambert_fixed with
# gooding_second_starter.diff applied. Needs a Basilisk source tree (commit ea95ec4 was used),
# g++ and Eigen 3 headers. Usage: sh build.sh /path/to/basilisk
set -e
B="$1/src"
EIGEN="${EIGEN:-/usr/include/eigen3}"
CXX="g++ -std=c++17 -O2 -I$B -I$EIGEN"
$CXX -c driver.cpp "$B/architecture/_GeneralModuleFiles/sys_model.cpp" "$B/architecture/utilities/bskLogging.cpp" \
     "$B/architecture/utilities/avsEigenSupport.cpp" "$B/architecture/utilities/moduleIdGenerator/moduleIdGenerator.cpp"
gcc -O2 -I"$B" -c "$B/architecture/utilities/rigidBodyKinematics.c" "$B/architecture/utilities/linearAlgebra.c"
$CXX -c "$B/fswAlgorithms/orbitControl/lambertSolver/lambertSolver.cpp" -o lambertSolver_orig.o
cp "$B/fswAlgorithms/orbitControl/lambertSolver/lambertSolver.cpp" lambertSolver_fixed.cpp
patch lambertSolver_fixed.cpp < gooding_second_starter.diff
$CXX -c lambertSolver_fixed.cpp -o lambertSolver_fixed.o
COMMON="driver.o sys_model.o bskLogging.o avsEigenSupport.o moduleIdGenerator.o rigidBodyKinematics.o linearAlgebra.o"
g++ -o lambert_orig $COMMON lambertSolver_orig.o
g++ -o lambert_fixed $COMMON lambertSolver_fixed.o
echo "built lambert_orig and lambert_fixed"
