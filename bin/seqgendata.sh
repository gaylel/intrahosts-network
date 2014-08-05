#!/bin/bash
tr=$1
outfile=$2
dparamfile=$3


mrstr=$(cat $dparamfile | grep mr)
arr=(${mrstr//=/ })
mr=${arr[1]%%,}

slstr=$(cat $dparamfile | grep seqlen)
arr=(${slstr//=/ })
sl=${arr[1]%%,}

seq-gen -mHKY -s$mr -l${sl} < $tr > $outfile