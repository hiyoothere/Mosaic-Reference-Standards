#!/bin/bash
#$ -cwd

ID=$1
REF=$2
INTERVAL=$3

BIN_VERSION="1.0.0"

sudo docker run \
-v "/home/hiyoothere/MRS/3.aligned":"/input" \
-v "/home/hiyoothere/MRS/DV":"/output" \
-v "/home/hiyoothere/reference":"/reference" \
google/deepvariant:1.0.0 \
/opt/deepvariant/bin/run_deepvariant \
--model_type=WES \
--ref=$REF \
--reads=/input/$ID'.RGadded.marked.realigned.fixed.recal.LA.bam' \
--output_vcf=/output/$ID'.Target.DV.vcf' \
--output_gvcf=/output/$ID'.Target.DV.gvcf' \
--regions=$INTERVAL \
--num_shards 4
