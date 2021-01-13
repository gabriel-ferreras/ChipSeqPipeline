#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j yes

RES_DIR=$1
NUM_SAMPLES=$2
EXP=$3
BROAD=$4
INS_DIR=$5
UPSTREAM=$6
DOWNSTREAM=$7
MOTIFLENGTH=$8
MOTIFSIZE=$9
NUM_EXP=${10}
EXP_DESIGN=${11}
j=${12}

## Accessing results folder
cd $RES_DIR

## Creating sample result folder
mkdir sample_${j}_result
cd sample_${j}_result

## Peak determination
echo ""
echo "=========================="
echo "|   PEAK DETERMINATION   |"
echo "=========================="
echo ""

if [ $BROAD -eq 0 ]
then
	echo ""
	echo " Determining NARROW peaks for sample $j"
	echo ""
        macs2 callpeak -t ../samples/chip_$j/chip_$j.bam -c ../samples/control_$j/control_$j.bam -f BAM --outdir . -n ${EXP}_sample_$j

else
	echo ""
	echo " Determining BROAD peaks for sample $j"
	echo ""
        macs2 callpeak -t ../samples/chip_$j/chip_$j.bam -c ../samples/control_$j/control_$j.bam -f BAM --outdir . -n ${EXP}_sample_$j --broad

fi

echo ""
echo "==============================================="
echo "|   INITIATING R SCRIPT FOR PEAK ANNOTATION   |"
echo "==============================================="
echo ""

if [ $BROAD -eq 0 ]
then
        echo ""
	echo " Annotating NARROW peaks for sample $j"
        echo ""
        Rscript ${INS_DIR}/ChipSeqPipeline/target_genes.R ${EXP}_sample_${j}_peaks.narrowPeak ${EXP}_sample_${j}_summits.bed $UPSTREAM $DOWNSTREAM ${EXP}_sample_${j}_peaks_targetgenes.txt ${EXP}_sample_${j}_summits_targetgenes.txt

else
	echo ""
        echo " Annotating BROAD peaks for sample $j"
	echo ""
        Rscript ${INS_DIR}/ChipSeqPipeline/target_genes.R ${EXP}_sample_${j}_peaks.broadPeak ${EXP}_sample_${j}_summits.bed $UPSTREAM $DOWNSTREAM ${EXP}_sample_${j}_peaks_targetgenes.txt ${EXP}_sample_${j}_summits_targetgenes.txt
fi

echo ""
echo "==========================="
echo "|   HOMER MOTIF FINDING   |"
echo "==========================="
echo ""

echo ""
echo " Finding motives for sample $j"
echo ""
mkdir motifs_sample_${j}
cd motifs_sample_${j}
findMotifsGenome.pl ../${EXP}_sample_${j}_summits.bed tair10 . -len $MOTIFLENGTH -size $MOTIFSIZE
cd ..
cd ..

echo ""
echo "================================"
echo "|  EXPERIMENT GLOBAL ANALYSIS  |"
echo "================================"
echo ""

##Communication with blackboard_2.
echo "Peak annotation and HOMER motif finding of sample $j done" >> blackboard_2.txt
NUM_PROC=$(wc -l ../../results/blackboard_2.txt | awk '{ print $1 }')

if [ NUM_PROC -eq NUM_SAMPLES ]
then
	
fi

echo ""
echo "Analysis complete!!"
echo ""
