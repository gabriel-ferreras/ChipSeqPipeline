#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j yes

RES_DIR=$1
NUM_SAMPLES=$2
EXP=$3
BROAD=$4

## Accessing results folder
cd $RES_DIR

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
		echo " Determining NARROW peaks"
		echo ""
        	macs2 callpeak -t ../samples/chip_$i/chip_$i.bam -c ../samples/control_$i/control_$i.bam -f BAM --outdir . -n ${EXP}_sample_$i
        	((i++))
	done

else
	i=1
	while [ $i -le $NUM_SAMPLES ]
	do
		echo ""
		echo " Determining BROAD peaks"
		echo ""
        	macs2 callpeak -t ../samples/chip_$i/chip_$i.bam -c ../samples/control_$i/control_$i.bam -f BAM --outdir . -n ${EXP}_sample_$i --broad
        	((i++))
	done
fi

echo ""
echo "Analysis complete!!"
echo ""
