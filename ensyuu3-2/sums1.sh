#!/bin/bash
echo "# Parallel with 4 cores and static"
OMP_NUM_THREADS=4 OMP_SCHEDULE="static,1" ~/ensyuu/ensyuu3-2/sum-mp
