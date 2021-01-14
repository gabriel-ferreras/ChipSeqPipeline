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
NUM_EXP=${12}
EXP_DESIGN=${13}

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
	fastqc control_${i}_1.fastq.gz
	fastqc control_${i}_2.fastq.gz
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

## Communication with blackboard_1.
echo "Processing Control $i done!" >> ../../results/blackboard_1.txt
NUM_PROC=$(wc -l ../../results/blackboard_1.txt | awk '{ print $1 }')
TOTAL_PROC=$((${NUM_SAMPLES}*2))

cd ../../results

if [ $NUM_PROC -eq $TOTAL_PROC ]
then
	j=1
	while [ $j -le $NUM_SAMPLES ]
	do
		echo ""
		echo "   Continuing with peak determination of sample $j, it is almost done!"
		echo ""
		qsub -o peaks_${j} -N peaks_${j} $INS_DIR/ChipSeqPipeline/peak_determination.sh $SAMPLE_DIR/../../results $NUM_SAMPLES $EXP $BROAD $INS_DIR $UPSTREAM $DOWNSTREAM $MOTIFLENGTH $MOTIFSIZE $NUM_EXP $EXP_DESIGN $j
		((j++))
fi
