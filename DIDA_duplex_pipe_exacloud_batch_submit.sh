#!/bin/bash
# DS bash script
# Version 2.1
# 

export TERM=xterm

echo Duplex_pipe_executable: "$0"
echo Sample Name "$1"

export PATH=$PATH:/home/exacloud/lustre1/users/peto/mutation_pipeline/tools/bwa-0.7.12/
export PATH=$PATH:/home/exacloud/lustre1/users/peto/mutation_pipeline/tools/samtools-1.2/
# Step 1: Setup variables for run:
clear

# Set up error checking
# Stop on any error
set -e
# Stop on an error inside a pipeline
set -o pipefail
# Throw an error on calling an unassigned variable
set -u

#DEFAULTS
DSpath=''
alignRef=''
runIdentifier=''
read1in=seq1.fq
read2in=seq2.fq
iSize=-1
minMem=3
maxMem=1000
cutOff=0.9
nCutOff=.1
readLength=98
barcodeLength=22
spacerLength=6
filtersSet='n'
readTypes='dpm'
repFilt=9
readOut=1000000

#NONDEFAULTS

DSpath='/home/exacloud/lustre1/users/butlerti/DIDA_scripts/'
alignRef='/home/exacloud/lustre1/users/peto/mutation_pipeline/genomes/human/g1k/human_g1k_v37.fasta'
runIdentifier="$1"
read1in="$1"_R1.fastq
read2in="$1"_R2.fastq
readLength=123
mpileup_region=/home/exacloud/lustre1/users/butlerti/DIDA_scripts/mpileup_region_files/good_dida1+dida_2_mut_file.txt # file listing mutations of interest being assayed
trim_length=3 # number of potentially low quality bases to trim read the start and end of the reads
trim1_start=1
HS_metrics_region=/home/exacloud/lustre1/users/butlerti/DIDA_scripts/HS_metrics_region_files/DIDA_well_performing_panel1+panel2_bail_intervals.txt # region file listing regions captured

#FINAL_READ_LENGTH
readLength=$((readLength-barcodeLength-spacerLength))

#cycles to trim
trim2_start=$((readLength-trim_length))

#LOG_FILE_NAME
logFile=${runIdentifier}.partial.log.txt

#Export all variables
export DSpath
export alignRef
export runIdentifier
export read1in
export read2in
export iSize
export minMem
export maxMem
export cutOff
export nCutOff
export readLength
export barcodeLength
export spacerLength
export filtersSet
export readTypes
export repFilt
export readOut
export mpileup_region
export trim_length
export trim1_start
export trim2_start
export HS_metrics_region

# Print out options used to log file
touch $logFile
echo "Run identifier: " $runIdentifier | tee -a ${logFile}
echo "Program path: " $DSpath | tee -a ${logFile}
echo "Reference genome: " $alignRef | tee -a ${logFile}
echo "Barcode length: " $barcodeLength | tee -a ${logFile}
echo "Spacer length: " $spacerLength | tee -a ${logFile}
echo "Post-tag_to_header read length: " $readLength | tee -a ${logFile}
echo "Repetitive tag filter length: " $repFilt | tee -a ${logFile}
echo "Minimum family size: " $minMem | tee -a ${logFile}
echo "Maximum family size: " $maxMem | tee -a ${logFile}
echo "Consensus cutoff: " $cutOff | tee -a ${logFile}
echo "Consensus N cutoff: " $nCutOff | tee -a ${logFile}
echo "Read types: " $readTypes | tee -a ${logFile}
echo "Filters: " $filtersSet | tee -a ${logFile}
echo "mpileup region" $mpileup_region | tee -a ${logFile}
echo "trimming length" $trim_length | tee -a ${logFile}
echo "trimming cycles: "$trim1_start"-"$trim_length" , "$trim2_start"-"$readLength  | tee -a ${logFile}
echo "Region File for HS Metrics: " $HS_metrics_region | tee -a ${logFile}
echo "" | tee -a ${logFile}

#  Step 2: Run tag_to_header.py on imput files

echo "Starting Run" | tee -a ${logFile}
echo "tag_to_header starting"  | tee -a ${logFile}
date | tee -a ${logFile}
echo "" | tee -a ${logFile}

/opt/installed/python-2.7.8/bin/python2 ${DSpath}/tag_to_header.py --infile1 $read1in --infile2 $read2in --outfile1 ${runIdentifier}.seq1.fq.smi --outfile2 ${runIdentifier}.seq2.fq.smi --barcode_length $barcodeLength --spacer_length $spacerLength --read_out $readOut

# Step 3: Align sequences

echo "Aligning with BWA" | tee -a ${logFile}
date | tee -a ${logFile}

# bwa aln $alignRef ${runIdentifier}.seq1.fq.smi > ${runIdentifier}.seq1.aln
# bwa aln $alignRef ${runIdentifier}.seq2.fq.smi > ${runIdentifier}.seq2.aln
# bwa sampe -s $alignRef ${runIdentifier}.seq1.aln ${runIdentifier}.seq2.aln ${runIdentifier}.seq1.fq.smi ${runIdentifier}.seq2.fq.smi > ${runIdentifier}.pe.sam
/home/exacloud/lustre1/users/peto/mutation_pipeline/tools/bwa-0.7.12/bwa mem $alignRef -t 8 -Ma ${runIdentifier}.seq1.fq.smi ${runIdentifier}.seq2.fq.smi > ${runIdentifier}.pe.sam


# Step 4: Sort aligned sequences
echo "Sorting aligned sequences" | tee -a ${logFile}
date | tee -a ${logFile}

/home/exacloud/lustre1/users/peto/mutation_pipeline/tools/samtools-1.2/samtools view -Sbu ${runIdentifier}.pe.sam | samtools sort - ${runIdentifier}.pe.sort

# Step 5: Run Consensus Maker
echo "Starting Consensus Maker" | tee -a ${logFile}
date | tee -a ${logFile}

/opt/installed/python-2.7.8/bin/python2 ${DSpath}/ConsensusMaker_fq_output.py --infile ${runIdentifier}.pe.sort.bam --tagfile ${runIdentifier}.pe.tagcounts --outfile ${runIdentifier}.sscs.bam --minmem $minMem --maxmem $maxMem --readlength $readLength --cutoff $cutOff --Ncutoff $nCutOff --read_type $readTypes --filt $filtersSet --isize $iSize

#Stp 6: Create .fastq from sscs.bam + Align SSCS
echo "aligning sscs" | tee -a ${logFile}
date | tee -a ${logFile}

awk  'BEGIN { OFS = "\t" } ; { sub(/:1/, "", $1) }1' ${runIdentifier}.sscs.r1.fq | awk  'BEGIN { OFS = "\t" } ; { sub(/:2/, "", $1) }1' > ${runIdentifier}.sscs.r1.no_num.fq

awk  'BEGIN { OFS = "\t" } ; { sub(/:1/, "", $1) }1' ${runIdentifier}.sscs.r2.fq | awk  'BEGIN { OFS = "\t" } ; { sub(/:2/, "", $1) }1' > ${runIdentifier}.sscs.r2.no_num.fq

/home/exacloud/lustre1/users/peto/mutation_pipeline/tools/bwa-0.7.12/bwa mem $alignRef -t 8 -Ma ${runIdentifier}.sscs.r1.no_num.fq ${runIdentifier}.sscs.r2.no_num.fq > ${runIdentifier}.sscs.aln.sam

# Step 8: Sort SSCSs
echo "Sorting SSCSs" | tee -a ${logFile}
date | tee -a ${logFile}

/home/exacloud/lustre1/users/peto/mutation_pipeline/tools/samtools-1.2/samtools view -bu ${runIdentifier}.sscs.aln.sam | samtools sort - ${runIdentifier}.sscs.aln.sort

# Step 9: Add ReadGroup Info
echo "Adding Readgroups" | tee -a ${logFile}
date | tee -a ${logFile}

java -Xmx25g -jar /home/exacloud/lustre1/users/peto/mutation_pipeline/tools/picard-tools-1.119/AddOrReplaceReadGroups.jar INPUT=${runIdentifier}.sscs.aln.sort.bam OUTPUT=${runIdentifier}.sscs.aln.readgroups.bam RGLB=UW RGPL=Illumina RGPU=ATATAT RGSM=default

# Step 10: Index Sorted SSCS
echo "Indexing sorted SSCSs" | tee -a ${logFile}
date | tee -a ${logFile}

/home/exacloud/lustre1/users/peto/mutation_pipeline/tools/samtools-1.2/samtools index ${runIdentifier}.sscs.aln.readgroups.bam

# Step 11: Local Relignment of SSCS
echo "Local Realignment SSCSs" | tee -a ${logFile}
date | tee -a ${logFile}

java -Xmx25g -jar /home/exacloud/lustre1/users/peto/mutation_pipeline/tools/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $alignRef -I ${runIdentifier}.sscs.aln.readgroups.bam -o ${runIdentifier}.sscs.aln.sort.intervals

java -Xmx25g -jar /home/exacloud/lustre1/users/peto/mutation_pipeline/tools/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar -T IndelRealigner --filter_bases_not_stored -R $alignRef -I ${runIdentifier}.sscs.aln.readgroups.bam -targetIntervals ${runIdentifier}.sscs.aln.sort.intervals -o ${runIdentifier}.sscs.realign.bam

# Step 12: Trim Ends of Reads
echo "Trimming $trim_length from read ends" | tee -a ${logFile}
date | tee -a ${logFile}

java -Xmx25g -jar /home/exacloud/lustre1/users/peto/mutation_pipeline/tools/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar -T ClipReads -I ${runIdentifier}.sscs.realign.bam -o ${runIdentifier}.sscs.trimmed.bam -R $alignRef --cyclesToTrim "$trim1_start-$trim_length,$trim2_start-$readLength" --clipRepresentation WRITE_NS_Q0S

# Step 13: Clip Overlapping Reads
echo "Clipping Overlapping Reads" | tee -a ${logFile}
date | tee -a ${logFile}

/home/users/butlerti/bamUtil_1.0.13/bamUtil/bin/bam clipOverlap --in ${runIdentifier}.sscs.trimmed.bam --stats --out ${runIdentifier}.sscs.final.bam

# Step 13: Index Clipped Reads
echo "Indexing Overlapped Clipped Reads" | tee -a ${logFile}
date | tee -a ${logFile}

/home/exacloud/lustre1/users/peto/mutation_pipeline/tools/samtools-1.2/samtools index ${runIdentifier}.sscs.final.bam


# Step 13: Generate Mpileup For Muts of Interest
echo "generating mpileup" | tee -a ${logFile}
date | tee -a ${logFile}

perl /home/exacloud/lustre1/users/butlerti/DIDA_scripts/batch_mpileup.pl $mpileup_region ${runIdentifier}.sscs.clipped.bam | perl /home/exacloud/lustre1/users/butlerti/DIDA_scripts/mutation_calculator.pl - > ${runIdentifier}.mutation_stats

# Step 14: Calculate HS Metrics
echo "Calculating HS Metrics" | tee -a ${logFile}
date | tee -a ${logFile}

java -Xmx25g -jar /home/exacloud/lustre1/users/peto/mutation_pipeline/tools/picard-tools-1.119/CalculateHsMetrics.jar BI=$HS_metrics_region TI=$HS_metrics_region I=${runIdentifier}.pe.sort.bam O=${runIdentifier}.pe.HS_metrics

java -Xmx25g -jar /home/exacloud/lustre1/users/peto/mutation_pipeline/tools/picard-tools-1.119/CalculateHsMetrics.jar BI=$HS_metrics_region TI=$HS_metrics_region I=${runIdentifier}.sscs.final.bam O=${runIdentifier}.sscs.HS_metrics


# Step 15: Remove Intermediary Files
echo "Removing Intermediate Files" | tee -a ${logFile}
date | tee -a ${logFile}

rm ${runIdentifier}.seq*
rm ${runIdentifier}.sscs.trimmed*
rm ${runIdentifier}.sscs.r2*
rm ${runIdentifier}.sscs.r1*
rm ${runIdentifier}.sscs.fq
rm ${runIdentifier}.sscs.aln*
rm ${runIdentifier}.pe.sam
rm ${runIdentifier}.sscs.realign*



echo "Finised with run " $runIdentifier | tee -a ${logFile}
date | tee -a ${logFile}

#Generated on Fri Apr  3 14:04:34 2015 with the command: 
#/usr/bin/python/home/exacloud/lustre1/users/butlerti/DIDA/PE_BASH_MAKER.py --ref /home/exacloud/lustre1/users/peto/mutation_pipeline/genomes/human/g1k/human_g1k_v37.fasta --r1src 50ng_random_0.1_1.fastq --r2src 50ng_random_0.1_2.fastq --runIdentifier spellman_50ng_0.1 --rlength 81 

#Updated on Wed Oct 5 2016 to remove unused DCS steps in DIDA data and adding alignment steps of SSCS
