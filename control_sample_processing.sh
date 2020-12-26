#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j yes

SAMPLE_DIR=$1
i=$2
NUM_SAMPLES=$3
INS_DIR=$4
EXP=$5
BROAD=$6

echo ""
echo "============================="
echo "|   PROCESSING CONTROL $i   |"
echo "============================="
echo ""

cd $SAMPLE_DIR

## Sample quality control and read mapping to reference genome
fastqc control_$i.fastq.gz
bowtie2 -x ../../genome/index -U control_$i.fastq.gz -S control_$i.sam

## Generating sorted bam file
samtools sort -o control_$i.bam control_$i.sam
rm control_$i.sam
rm *.fastq.gz
samtools index control_$i.bam

echo ""
echo "   Control $i processing DONE!!"
echo ""

## Communication with blackboard.
echo "Processing Control $i done!" >> ../../results/blackboard.txt
NUM_PROC=$(wc -l ../../results/blackboard.txt | awk '{ print $1 }')
TOTAL_PROC=$((${NUM_SAMPLES}*2))

cd ../../results

if [ $NUM_PROC -eq $TOTAL_PROC ]
then
	echo ""
	echo "   Continuing with peak determination, you are almost done!"
	echo ""
	qsub -o peaks -N peaks $INS_DIR/ChipSeqPipeline/peak_determination.sh $SAMPLE_DIR/../../results $NUM_SAMPLES $EXP $BROAD $INS_DIR
fi
