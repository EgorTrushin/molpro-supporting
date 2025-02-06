#!/bin/bash -l

rm *.out *.xml *.log

MOLPRO=/home/trushin/Molpro/molpro/bin/molpro.exe
NUM_OMP_THREADS=16
MEM=2000

echo Molpro Path: $MOLPRO
echo Number of OMP threads: $NUM_OMP_THREADS
echo Memory: $MEM m

export OMP_NUM_THREADS=$NUM_OMP_THREADS

mpiexec -n 1 $MOLPRO -m $MEM --no-xml-output < b_uksinv.test
mpiexec -n 1 $MOLPRO -m $MEM --no-xml-output < be_acfd_scrpa.test
mpiexec -n 1 $MOLPRO -m $MEM --no-xml-output < co_acfd_rirpa.test
mpiexec -n 1 $MOLPRO -m $MEM --no-xml-output < co_acfd_scexx.test
mpiexec -n 1 $MOLPRO -m $MEM --no-xml-output < co_ksinv.test
mpiexec -n 1 $MOLPRO -m $MEM --no-xml-output < li_acfd_uscrpa.test
mpiexec -n 1 $MOLPRO -m $MEM --no-xml-output < li_ksinv.test
mpiexec -n 1 $MOLPRO -m $MEM --no-xml-output < nh2_acfd_urirpa.test
mpiexec -n 1 $MOLPRO -m $MEM --no-xml-output < nh2_acfd_uscexx.test