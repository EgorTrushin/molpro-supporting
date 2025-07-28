#!/bin/bash -l

rm *.out *.xml *.log

MOLPRO=$HOME/Molpro/molpro/bin/molpro
OMP_NUM_THREADS=16
MEM=2000

echo Molpro Path: $MOLPRO
echo Number of OMP threads: $OMP_NUM_THREADS
echo Memory: $MEM m

export OMP_NUM_THREADS=$OMP_NUM_THREADS

$MOLPRO -m $MEM --no-xml-output < b_uksinv.test
$MOLPRO -m $MEM --no-xml-output < be_acfd_scrpa.test
$MOLPRO -m $MEM --no-xml-output < co_acfd_dftoep.test
$MOLPRO -m $MEM --no-xml-output < co_acfd_rirpa.test
$MOLPRO -m $MEM --no-xml-output < co_acfd_scexx.test
$MOLPRO -m $MEM --no-xml-output < co_ksinv.test
$MOLPRO -m $MEM --no-xml-output < li_acfd_uscrpa.test
$MOLPRO -m $MEM --no-xml-output < li_ksinv.test
$MOLPRO -m $MEM --no-xml-output < n_acfd_udftoep.test
$MOLPRO -m $MEM --no-xml-output < nh2_acfd_urirpa.test
$MOLPRO -m $MEM --no-xml-output < nh2_acfd_uscexx.test
