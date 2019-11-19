#!/bin/bash -l

. $(dirname $0)/env.sh

scons install PREFIX=$PREFIX BOOST_ROOT=$BOOST_ROOT SAMTOOLS_DIR=$SAMTOOLS_DIR
cp .genepool_deploy/module_dependencies $PREFIX/.deps
