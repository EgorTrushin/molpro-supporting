#!/bin/bash -l

rm *.out *.xml *.log

MOLPRO=~/Molpro/molpro/bin/molpro
NUM_OMP_THREADS=16
MEM=2000

echo Molpro Path: $MOLPRO
echo Number of OMP threads: $NUM_OMP_THREADS
echo Memory: $MEM m

export OMP_NUM_THREADS=$NUM_OMP_THREADS

$MOLPRO -t $NUM_OMP_THREADS -m $MEM --no-xml-output < b_uksinv.test
$MOLPRO -t $NUM_OMP_THREADS -m $MEM --no-xml-output < be_acfd_scrpa.test
$MOLPRO -t $NUM_OMP_THREADS -m $MEM --no-xml-output < co_acfd_rirpa.test
$MOLPRO -t $NUM_OMP_THREADS -m $MEM --no-xml-output < co_acfd_scexx.test
$MOLPRO -t $NUM_OMP_THREADS -m $MEM --no-xml-output < co_ksinv.test
$MOLPRO -t $NUM_OMP_THREADS -m $MEM --no-xml-output < li_acfd_uscrpa.test
$MOLPRO -t $NUM_OMP_THREADS -m $MEM --no-xml-output < li_ksinv.test
$MOLPRO -t $NUM_OMP_THREADS -m $MEM --no-xml-output < nh2_acfd_urirpa.test
$MOLPRO -t $NUM_OMP_THREADS -m $MEM --no-xml-output < nh2_acfd_uscexx.test
