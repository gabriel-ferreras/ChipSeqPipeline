#! /bin/bash

#Print help message with usage if no parameters indicated.
if [ $# -ne 1 ]
then
    echo ""
    echo "    Usage: bash ChipSeqPipeline.sh <params_file>"
    echo ""
    echo "    params_file: input file with the parameters."
    echo "    Example of params_file in test_params."
    exit
fi

#Reading parameters file.
PARAMS=$1
echo ""
echo "======================"
echo "| LOADING PARAMETERS |"
echo "======================"
echo ""
INS_DIR=$(grep installation_directory: $PARAMS | awk '{ print $2 }')
echo "      Installation directory is "$INS_DIR

WORK_DIR=$(grep working_directory: $PARAMS | awk '{ print $2 }')
echo "      Working directory is "$WORK_DIR

ANALYSIS=$(grep analysis_name: $PARAMS | awk '{ print $2 }')
echo "      Analysis name = "$ANALYSIS

NUM_SAMPLES=$(grep number_samples: $PARAMS | awk '{ print $2 }')
echo "      Number samples = "$NUM_SAMPLES

BROAD=$(grep broad: $PARAMS | awk '{ print $2 }')
echo "      Broad peaks = "$BROAD

ANNOTATION=$(grep path_annotation: $PARAMS | awk '{ print $2 }')
echo "      Annotation in "$ANNOTATION

GENOME=$(grep path_genome: $PARAMS | awk '{ print $2 }')
echo "      Genome in "$GENOME

UPSTREAM=$(grep upstream_limit: $PARAMS | awk '{ print $2 }')
echo "      Upstream limit for promoter is "$UPSTREAM

DOWNSTREAM=$(grep downstream_limit: $PARAMS | awk '{ print $2 }')
echo "      Downstream limit for promoter is "$DOWNSTREAM

MOTIFLENGTH=$(grep motif_length: $PARAMS | awk '{ print $2 }')
echo "      Selected motif length is "$MOTIFLENGTH

MOTIFSIZE=$(grep motif_size: $PARAMS | awk '{ print $2 }')
echo "      Selected motif size is "$MOTIFSIZE

PAIRED=$(grep paired: $PARAMS | awk '{ print $2 }')
echo "      Paired-end lectures = "$PAIRED

NUM_EXP=$(grep number_experiments: $PARAMS | awk '{ print $2 }')
echo "      Number of experiments = "$NUM_EXP

EXP_DESIGN=$(grep experimental_design: $PARAMS | awk '{ print $2 }')
echo "      Experimental design =" $EXP_DESIGN

CHR=$(grep chromosome: $PARAMS | awk '{ print $2 }')
echo "      Chromosome selected for analysis is" $CHR

if [ $PAIRED -eq 0 ]
then
	CHIPS=()
	i=0
	while [ $i -lt $NUM_SAMPLES ]
	do
        		j=$(( $i + 1 ))
        		CHIPS[$i]=$(grep path_chip_$j: $PARAMS | awk '{ print $2 }')
        		echo "      Chip $j in ${CHIPS[$i]}"
        		((i++))
	done

	CONTROLS=()
	i=0
	while [ $i -lt $NUM_SAMPLES ]
		do
        		j=$(( $i + 1 ))
        		CONTROLS[$i]=$(grep path_control_$j: $PARAMS | awk '{ print $2 }')
        		echo "      Control $j in ${CONTROLS[$i]}"
        		((i++))
	done

else
	CHIPS=()
        i=0
        while [ $i -lt $NUM_SAMPLES ]
        do
                j=$(( $i + 1 ))
                CHIPS[$i]=$(grep path_chip_${j}_1: $PARAMS | awk '{ print $2 }')
                echo "      Chip $j paired 1 in ${CHIPS[$i]}"
		CHIPS[$j]=$(grep path_chip_${j}_2: $PARAMS | awk '{ print $2 }')
                echo "      Chip $j paired 2 in ${CHIPS[$j]}"
                i=$(( $i + 2 ))
        done

        CONTROLS=()
        i=0
        while [ $i -lt $NUM_SAMPLES ]
        do
                j=$(( $i + 1 ))
                CONTROLS[$i]=$(grep path_control_${j}_1: $PARAMS | awk '{ print $2 }')
                echo "      Control $j paired 1 in ${CONTROLS[$i]}"
		CONTROLS[$j]=$(grep path_control_${j}_2: $PARAMS | awk '{ print $2 }')
                echo "      Control $j paired 2 in ${CONTROLS[$j]}"
                i=$(( $i + 2 ))
        done
fi

#Preparing working workspace.
echo ""
echo "======================"
echo "| CREATING WORKSPACE |"
echo "======================"
echo ""
echo "      Making directories"

cd $WORK_DIR
mkdir $ANALYSIS
cd $ANALYSIS
mkdir results genome annotation samples
cd samples
i=1
while [ $i -le $NUM_SAMPLES ]
do
        mkdir chip_$i
	mkdir control_$i
        ((i++))
done
cd ..

#Copying the data.
echo ""
echo "      Copying input data from the indicated paths"
echo ""
cp $ANNOTATION $WORK_DIR/$ANALYSIS/annotation/annotation.gtf
cp $GENOME $WORK_DIR/$ANALYSIS/genome/genome.fa
if [ $PAIRED -eq 0 ]
then
	i=1
	while [ $i -le $NUM_SAMPLES ]
	do
        	j=$((i - 1))
        	cp ${CHIPS[j]} $WORK_DIR/$ANALYSIS/samples/chip_$i/chip_$i.fastq.gz
		cp ${CONTROLS[j]} $WORK_DIR/$ANALYSIS/samples/control_$i/control_$i.fastq.gz
        	((i++))
	done
else
	echo "UNPAIRED"
	i=1
        while [ $i -le $NUM_SAMPLES ]
        do
                j=$((i - 1))
                cp ${CHIPS[j]} $WORK_DIR/$ANALYSIS/samples/chip_$i/chip_${i}_1.fastq.gz
                cp ${CONTROLS[j]} $WORK_DIR/$ANALYSIS/samples/control_$i/control_${i}_1.fastq.gz
                cp ${CHIPS[i]} $WORK_DIR/$ANALYSIS/samples/chip_$i/chip_${i}_2.fastq.gz
                cp ${CONTROLS[i]} $WORK_DIR/$ANALYSIS/samples/control_$i/control_${i}_2.fastq.gz
		i=$(( $i + 2 ))
        done

fi

#Creating genome index.
echo ""
echo "========================="
echo "| CREATING GENOME INDEX |"
echo "========================="
echo ""
cd $WORK_DIR/$ANALYSIS/genome
bowtie2-build genome.fa index


#Processing samples.
echo ""
echo "======================="
echo "|  SAMPLE PROCESSING  |"
echo "======================="
echo ""

cd ../results

i=1
while [ $i -le $NUM_SAMPLES ]
do
        echo "Sent to processing chip $i"
	qsub -o chip_$i -N chip_$i $INS_DIR/ChipSeqPipeline/chip_sample_processing.sh $WORK_DIR/$ANALYSIS/samples/chip_$i $i $NUM_SAMPLES $INS_DIR $ANALYSIS $BROAD $PAIRED $UPSTREAM $DOWNSTREAM $MOTIFLENGTH $MOTIFSIZE $NUM_EXP $EXP_DESIGN $CHR
        echo "Sent to processing control $i"
	qsub -o control_$i -N control_$i $INS_DIR/ChipSeqPipeline/control_sample_processing.sh $WORK_DIR/$ANALYSIS/samples/control_$i $i $NUM_SAMPLES $INS_DIR $ANALYSIS $BROAD $PAIRED $UPSTREAM $DOWNSTREAM $MOTIFLENGTH $MOTIFSIZE $NUM_EXP $EXP_DESIGN $CHR
	((i++))
done
