#!/bin/bash
#$ -cwd
#$ -S /bin/bash


# basic argument
DIR=$1
SAMPLE=$2
REF=$3
DataPath=$4
FASTQ1=$5
FASTQ2=$6
PICARD=$7
OPT=''


date

# 1. Actual Alignment. -I option to use illumina 1.3+ quailities. For the latest version, we don't need -I option. 

bwa mem -t 4 -M $REF $DataPath$FASTQ1 $DataPath$FASTQ2 | java  -jar $PICARD/picard.jar SortSam \
SORT_ORDER=coordinate \
INPUT=/dev/stdin \
OUTPUT=$AlignedPath$SAMPLE$OPT.bam \
VALIDATION_STRINGENCY=LENIENT \
CREATE_INDEX=true \
TMP_DIR=$DIR/out/tmp/

# 2.1. Add or replace read groups
java -jar $PICARD/picard.jar AddOrReplaceReadGroups \
INPUT=$AlignedPath$SAMPLE$OPT.bam \
OUTPUT=$AlignedPath$SAMPLE$OPT.RGadded.bam \
SORT_ORDER=coordinate \
RGLB='MRS' \
RGPL='il' \
RGPU='WES' \
RGSM=$SAMPLE \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=$DIR/out/tmp/

# 3. Marking PCR duplicates
java -jar $PICARD/picard.jar MarkDuplicates \
INPUT=$AlignedPath$SAMPLE$OPT.RGadded.bam \
OUTPUT=$AlignedPath$SAMPLE$OPT.RGadded.marked.bam \
METRICS_FILE=$AlignedPath$SAMPLE$OPT.metrics \
CREATE_INDEX=true \
REMOVE_DUPLICATES=true \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=$DIR/out/tmp/

# 3-1. LeftAlignIndel 
/data/project/MRS/Resource/gatk-4.1.5.0/gatk LeftAlignIndels \
-R $REF \
-I $AlignedPath$SAMPLE$OPT.RGadded.marked.bam \
-O $AlignedPath$SAMPLE$OPT.RGadded.marked.LA.bam 


# 4. Local realignment around indels
#step1. To create a table of possible indels
java  -jar /opt/Yonsei/GATK/3.8-1/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R $REF \
-I $AlignedPath$SAMPLE$OPT.RGadded.marked.LA.bam \
-o $AlignedPath$SAMPLE$OPT.RGadded.marked.LA.forIndelRealigner.list \


# 5. [step2] To realign reads around indels targets
java  -jar /opt/Yonsei/GATK/3.8-1/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R $REF \
-I $AlignedPath$SAMPLE$OPT.RGadded.marked.LA.bam \
-targetIntervals $AlignedPath$SAMPLE$OPT.RGadded.marked.LA.forIndelRealigner.list \
-o $AlignedPath$SAMPLE$OPT.RGadded.marked.LA.realigned.bam \


# 6. The mate information must be fixed#
java -jar $PICARD/picard.jar FixMateInformation \
INPUT=$AlignedPath$SAMPLE$OPT.RGadded.marked.LA.realigned.bam \
OUTPUT=$AlignedPath$SAMPLE$OPT.RGadded.marked.LA.realigned.fixed.bam \
SO=coordinate \
VALIDATION_STRINGENCY=LENIENT \
CREATE_INDEX=true \
TMP_DIR=$DIR/out/tmp/


#7. Base Quality score recalibration.
# 7.1) BaseRecalibrator_1st pass
gatk  BaseRecalibrator \
-R $REF \
-I $AlignedPath$SAMPLE$OPT.RGadded.marked.LA.realigned.fixed.bam \
--known-sites $dbSNP \
-O $AlignedPath$SAMPLE$OPT.recal_pass1.table \
--tmp-dir $DIR/out/tmp/

## 7.2) ApplyBQSR
gatk  ApplyBQSR \
-R $REF \
-I $AlignedPath$SAMPLE$OPT.RGadded.marked.LA.realigned.fixed.bam \
--bqsr-recal-file $AlignedPath$SAMPLE$OPT.recal_pass1.table \
-O $AlignedPath$SAMPLE$OPT.RGadded.marked.LA.realigned.fixed.pre_recal.bam \
--tmp-dir $DIR/out/tmp/

# 7.3) BaseRecalibrator_2st pass
gatk  BaseRecalibrator \
-R $REF \
-I $AlignedPath$SAMPLE$OPT.RGadded.marked.LA.realigned.fixed.pre_recal.bam \
--known-sites $dbSNP \
-O $AlignedPath$SAMPLE$OPT.recal_pass2.table \
--tmp-dir $DIR/out/tmp/

# 7.4) ApplyBQSR
gatk  ApplyBQSR \
-R $REF \
-I $AlignedPath$SAMPLE$OPT.RGadded.marked.LA.realigned.fixed.pre_recal.bam \
--bqsr-recal-file $AlignedPath$SAMPLE$OPT.recal_pass2.table \
-O $AlignedPath$SAMPLE$OPT.RGadded.marked.LA.realigned.fixed.recal.bam \
--tmp-dir $DIR/out/tmp/

# 8. Remove intermediate files
rm $AlignedPath$SAMPLE$OPT.bam
rm $AlignedPath$SAMPLE$OPT.RGadded.bam
rm $AlignedPath$SAMPLE$OPT.RGadded.marked.bam
rm $AlignedPath$SAMPLE$OPT.RGadded.marked.realigned.bam
rm $AlignedPath$SAMPLE$OPT.RGadded.marked.realigned.fixed.bam
rm $AlignedPath$SAMPLE$OPT.RGadded.marked.realigned.fixed.pre_recal.bam


date




