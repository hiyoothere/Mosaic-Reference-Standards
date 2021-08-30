#!/bin/bash


## Input Parameters
IN_BAM=$1
OUT_PUP=$2
REF=$3

## Pileup
smart_pileup() {
    local bam=$1
    local ref=$2
    local pup=$3
    if [[ -f $pup ]]; then
        return
    elif [[ $bam =~ "," ]]; then
        samtools mpileup -q "1" -Q "1" -R -f $ref -b $bam -o $pup
    else
        samtools mpileup -q "1" -Q "1" -R -f $ref -o $pup $bam
    fi
}

echo "** Input: ${IN_BAM}"
echo "** Output: ${OUT_PUP}"
smart_pileup $IN_BAM $REF $OUT_PUP
