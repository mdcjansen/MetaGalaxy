#!/bin/bash -l

. $(dirname $0)/env.sh

scons test BOOST_ROOT=$BOOST_ROOT SAMTOOLS_DIR=$SAMTOOLS_DIR

