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
PAIRED=$7
UPSTREAM=$8
DOWNSTREAM=$9
MOTIFLENGTH=${10}
MOTIFSIZE=${11}

echo ""
echo "============================="
echo "|   PROCESSING CONTROL $i   |"
echo "============================="
echo ""

cd $SAMPLE_DIR

## Sample quality control and read mapping to reference genome

if [ $PAIRED -eq 0 ]
then
	fastqc control_$i.fastq.gz
	bowtie2 -x ../../genome/index -U control_$i.fastq.gz -S control_$i.sam
else
	fastqc control_${1}_1.fastq.gz
	fastqc control_${1}_2.fastq.gz
	bowtie2 -x ../../genome/index -1 control_${i}_1.fastq.gz -2 control_${i}_2.fastq.gz -S control_$i.sam
fi

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
	qsub -o peaks -N peaks $INS_DIR/ChipSeqPipeline/peak_determination.sh $SAMPLE_DIR/../../results $NUM_SAMPLES $EXP $BROAD $INS_DIR $UPSTREAM $DOWNSTREAM $MOTIFLENGTH $MOTIFSIZE
fi
