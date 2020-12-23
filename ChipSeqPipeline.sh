#! /bin/bash

#Print help message with usage if no parameters indicated.
if [ $# -ne 1 ]
then
    echo ""
    echo "    Usage: bash ChipSeqPipeline.sh <params_file>"
    echo ""
    echo "    params_file: imput file with the parameters."
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

EXP=$(grep experiment_name: $PARAMS | awk '{ print $2 }')
echo "      Experiment name = "$EXP

NUM_SAMPLES=$(grep number_samples: $PARAMS | awk '{ print $2 }')
echo "      Number samples = "$NUM_SAMPLES

BROAD=$(grep broad: $PARAMS | awk '{ print $2 }')
echo "      Broad peaks = "$BROAD

ANNOTATION=$(grep path_annotation: $PARAMS | awk '{ print $2 }')
echo "      Annotation in "$ANNOTATION

GENOME=$(grep path_genome: $PARAMS | awk '{ print $2 }')
echo "      Genome in "$GENOME

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

#Preparing working workspace.
echo ""
echo "======================"
echo "| CREATING WORKSPACE |"
echo "======================"
echo ""
cd $WORK_DIR
mkdir $EXP
cd $EXP
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
cp $ANNOTATION $WORK_DIR/$EXP/annotation/annotation.gtf
cp $GENOME $WORK_DIR/$EXP/genome/genome.fa
i=1
while [ $i -le $NUM_SAMPLES ]
do
        j=$((i - 1))
        cp ${CHIPS[j]} $WORK_DIR/$EXP/samples/chip_$i/chip_$i.fastq.gz
	cp ${CONTROLS[j]} $WORK_DIR/$EXP/samples/control_$i/control_$i.fastq.gz
        ((i++))
done

#Creating genome index.
echo ""
echo "========================="
echo "| CREATING GENOME INDEX |"
echo "========================="
echo ""
cd $WORK_DIR/$EXP/genome
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
	qsub -o chip_$i -N chip_$i $INS_DIR/ChipSeqPipeline/chip_sample_processing.sh $WORK_DIR/$EXP/samples/chip_$i $i $NUM_SAMPLES $INS_DIR $EXP $BROAD
        echo "Sent to processing control $i"
	qsub -o control_$i -N control_$i $INS_DIR/ChipSeqPipeline/control_sample_processing.sh $WORK_DIR/$EXP/samples/control_$i $i $NUM_SAMPLES $INS_DIR $EXP $BROAD
	((i++))
done
