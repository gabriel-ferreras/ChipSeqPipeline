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
echo "=========================="
echo "|   PROCESSING CHIP $i   |"
echo "=========================="
echo ""

cd $SAMPLE_DIR

## Sample quality control and read mapping to reference genome
if [ $PAIRED -eq 0 ]
then
        fastqc chip_$i.fastq.gz
        bowtie2 -x ../../genome/index -U chip_$i.fastq.gz -S chip_$i.sam
else
	fastqc chip_${i}_1.fastq.gz
        fastqc chip_${i}_2.fastq.gz
        bowtie2 -x ../../genome/index -1 chip_${i}_1.fastq.gz -2 chip_${i}_2.fastq.gz -S chip_$i.sam
fi

## Generating sorted bam file
samtools sort -o chip_$i.bam chip_$i.sam
rm chip_$i.sam
rm *.fastq.gz
samtools index chip_$i.bam

echo ""
echo "   Chip $i processing DONE!!"
echo ""

## Communication with blackboard.
echo "Processing Chip $i done!" >> ../../results/blackboard.txt
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
