#!/bin/bash -l

set -ex
module purge
module load PrgEnv-gnu/4.9
module load binutils/2.28
module load scons
module load boost/1.59.0

