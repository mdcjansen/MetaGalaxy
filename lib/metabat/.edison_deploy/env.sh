#!/bin/bash -l

set -ex
module swap PrgEnv-intel PrgEnv-gnu/4.9
module load boost
module load binutils
module load scons 
module load samtools

