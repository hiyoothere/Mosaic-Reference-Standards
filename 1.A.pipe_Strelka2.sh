#!/bin/bash

ID=$1
REF=$2
AlignedPath=$3
INTERVAL=$4

StrelkaPath=/opt/Yonsei/Strelka2/2.9.10/bin
AnalysisPath=/data/project/MRS/0.Genotype/4.analysis/Strelka

$StrelkaPath/configureStrelkaGermlineWorkflow.py \
--bam=$AlignedPath$ID'.RGadded.marked.realigned.fixed.recal.LA.bam' \
--referenceFasta=$REF \
--exome --disableSequenceErrorEstimation \
--callRegions=$INTERVAL'.gz' \
--runDir=$AnalysisPath/$ID

$AnalysisPath/$ID/runWorkflow.py -m local -j 20
