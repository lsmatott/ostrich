#!/bin/bash

mkdir -p ../bin

make GCC_DBG 2>&1 | tee my_make_GCC_DBG.log
if [ -f Ostrich ]; then
  mv Ostrich ../bin/OstrichDebug
fi

make GCC 2>&1 | tee my_make_GCC.log
if [ -f Ostrich ]; then
   mv Ostrich ../bin/Ostrich
fi

make GCC_MPI 2>&1 | tee my_make_GCC_MPI.log
if [ -f OstrichMpi ]; then
   mv OstrichMpi ../bin/OstrichMpi
fi
