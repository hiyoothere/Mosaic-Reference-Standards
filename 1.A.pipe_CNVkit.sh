#!/bin/bash
#$ -S /bin/bash
#$ -cwd
DIR=$1
REF=$2
Inputpath=$3
Analysispath=$2
ID=$3


cnvkit.py coverage $Inputpath/$PRI \
$Analysispath/my_targets.bed \
-o $Analysispath/$ID.targetcoverage.cnn

cnvkit.py coverage $Inputpath/$PRI \
$Analysispath/antitargets.bed \
-o $Analysispath/$ID.antitargetcoverage.cnn


cnvkit.py reference \
$Analysispath/FASTP_RF-1_RF_il_WES.RGadded.marked.realigned.fixed.targetcoverage.cnn \
$Analysispath/FASTP_RF-2_RF_il_WES.RGadded.marked.realigned.fixed.targetcoverage.cnn \
$Analysispath/FASTP_RF-3_RF_il_WES.RGadded.marked.realigned.fixed.targetcoverage.cnn \
$Analysispath/FASTP_RF-4_RF_il_WES.RGadded.marked.realigned.fixed.targetcoverage.cnn \
$Analysispath/FASTP_RF-5_RF_il_WES.RGadded.marked.realigned.fixed.targetcoverage.cnn \
-f $REF \
-o $Analysispath/FOR_RF6.reference.cnn


cnvkit.py fix \
$Analysispath/FASTP_$ID'_RF_il_WES.RGadded.marked.realigned.fixed.targetcoverage.cnn' \
$Analysispath/FASTP_$ID'_RF_il_WES.antitargetcoverage.cnn' \
$Analysispath/FOR_$ID'.reference.cnn' \
-o $Analysispath/$ID'.cnr'

Analysispath=/data/project/MRS/0.Genotype/4.analysis/CNVkit 
cnvkit.py segment $Analysispath/$ID.cnr \
-o $Analysispath/$ID.cns

cnvkit.py scatter $Analysispath/$ID.cnr \
-s $Analysispath/$ID.cns \
-o $Analysispath/$ID.pdf \
--y-max 4 \
--y-min -4
