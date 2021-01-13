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

## Creating sample folder
mkdir sample_$i

## Peak determination
echo ""
echo "=========================="
echo "|   PEAK DETERMINATION   |"
echo "=========================="
echo ""

if [ $BROAD -eq 0 ]
then
	i=1
	while [ $i -le $NUM_SAMPLES ]
	do
		echo ""
		echo " Determining NARROW peaks for sample $i"
		echo ""
        	macs2 callpeak -t ../samples/chip_$i/chip_$i.bam -c ../samples/control_$i/control_$i.bam -f BAM --outdir . -n ${EXP}_sample_$i
        	((i++))
	done

else
	i=1
	while [ $i -le $NUM_SAMPLES ]
	do
		echo ""
		echo " Determining BROAD peaks for sample $i"
		echo ""
        	macs2 callpeak -t ../samples/chip_$i/chip_$i.bam -c ../samples/control_$i/control_$i.bam -f BAM --outdir . -n ${EXP}_sample_$i --broad
        	((i++))
	done
fi

echo ""
echo "==============================================="
echo "|   INITIATING R SCRIPT FOR PEAK ANNOTATION   |"
echo "==============================================="
echo ""

if [ $BROAD -eq 0 ]
then
        i=1
        while [ $i -le $NUM_SAMPLES ]
        do
                echo ""
                echo " Annotating NARROW peaks for sample $i"
                echo ""
                Rscript ${INS_DIR}/ChipSeqPipeline/target_genes.R ${EXP}_sample_${i}_peaks.narrowPeak ${EXP}_sample_${i}_summits.bed $UPSTREAM $DOWNSTREAM ${EXP}_sample_${i}_peaks_targetgenes.txt ${EXP}_sample_${i}_summits_targetgenes.txt
                ((i++))
        done

else
        i=1
        while [ $i -le $NUM_SAMPLES ]
        do
                echo ""
                echo " Annotating BROAD peaks for sample $i"
                echo ""
                Rscript ${INS_DIR}/ChipSeqPipeline/target_genes.R ${EXP}_sample_${i}_peaks.broadPeak ${EXP}_sample_${i}_summits.bed $UPSTREAM $DOWNSTREAM ${EXP}_sample_${i}_peaks_targetgenes.txt ${EXP}_sample_${i}_summits_targetgenes.txt
                ((i++))
        done
fi

echo ""
echo "==========================="
echo "|   HOMER MOTIF FINDING   |"
echo "==========================="
echo ""

i=1
while [ $i -le $NUM_SAMPLES ]
do
	echo ""
	echo " Finding motives for sample $i"
	echo ""
	mkdir motifs_sample_${i}
	cd motifs_sample_${i}
	findMotifsGenome.pl ../${EXP}_sample_${i}_summits.bed tair10 . -len $MOTIFLENGTH -size $MOTIFSIZE
	((i++))
	cd ..
done

echo ""
echo "Analysis complete!!"
echo ""
