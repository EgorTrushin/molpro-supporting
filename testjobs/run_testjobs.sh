#!/bin/bash -l

rm *.out *.xml *.log

MOLPRO=/home/trushin/Molpro/molpro-fork/bin/molpro
NUM_OMP_THREADS=16

echo Molpro Path: $MOLPRO
echo Number of OMP threads: $NUM_OMP_THREADS

export OMP_THREAD_LIMIT=$NUM_OMP_THREADS

$MOLPRO -t $NUM_OMP_THREADS --no-xml-output < co_acfd_rirpa.test
$MOLPRO -t $NUM_OMP_THREADS --no-xml-output < co_acfd_scexx.test
$MOLPRO -t $NUM_OMP_THREADS --no-xml-output < co_ksinv.test
$MOLPRO -t $NUM_OMP_THREADS --no-xml-output < nh2_acfd_urirpa.test
$MOLPRO -t $NUM_OMP_THREADS --no-xml-output < nh2_acfd_uscexx.test
