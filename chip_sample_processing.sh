#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j yes

SAMPLE_DIR=$1
i=$2
NUM_SAMPLES=$3
INS_DIR=$4
EXP=$5

echo ""
echo "============================"
echo "|   PROCESSING SAMPLE $i   |"
echo "============================"
echo ""

cd $SAMPLE_DIR

## Sample quality control and read mapping to reference genome
fastqc sample_$i.fastq.gz
bowtie2 -x ../../genome/index -U sample_$i.fastq.gz -S sample_$i.sam

## Generating sorted bam file
samtools sort -o sample_$i.bam sample_$i.sam
rm sample_$i.sam
rm *.fastq.gz
samtools index sample_$i.bam

## Communication with blackboard.
echo "Processing Sample $i done!" >> ../../results/blackboard.txt
NUM_PROC=$(wc -l ../../results/blackboard.txt | awk '{ print $1 }')

if [ $NUM_PROC -eq $NUM_SAMPLES ]
then
        qsub -o peaks -N peaks $INS_DIR/ChipSeqPipeline/peak_determination.sh $SAMPLE_DIR/../../results $NUM_SAMPLES $EXP
fi

echo ""
echo "   Sample $i processing DONE!!"
echo ""
