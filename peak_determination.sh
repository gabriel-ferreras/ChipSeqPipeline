#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j yes

RES_DIR=$1
NUM_SAMPLES=$2
EXP=$3

## Accessing results folder
cd $RES_DIR

## Peak determination
echo ""
echo "=========================="
echo "|   PEAK DETERMINATION   |"
echo "=========================="
echo ""

NUM_REP=(($NUM_SAMPLES))
macs2 callpeak -t ../samples/sample_$i -c ../samples/sample_$i.bam -f BAM --outdir . -n $EXP

echo "Analysis complete :)"
