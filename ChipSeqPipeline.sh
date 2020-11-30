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

ANNOTATION=$(grep path_annotation: $PARAMS | awk '{ print $2 }')
echo "      Annotation in "$ANNOTATION

GENOME=$(grep path_genome: $PARAMS | awk '{ print $2 }')
echo "      Genome in "$GENOME

SAMPLES=()
i=0
while [ $i -lt $NUM_SAMPLES ]
do
        j=$(( $i + 1 ))
        SAMPLES[$i]=$(grep path_sample_$j: $PARAMS | awk '{ print $2 }')
        echo "      Sample $j in ${SAMPLES[$i]}"
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
        mkdir sample_$i
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
        cp ${SAMPLES[j]} $WORK_DIR/$EXP/samples/sample_$i/sample_$i.fastq.gz
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
        qsub -o sample_$i -N sample_$i $INS_DIR/ChipSeqPipeline/chip_sample_processing.sh $WORK_DIR/$EXP/samples/sample_$i $i $NUM_SAMPLES $INS_DIR $EXP
        ((i++))
done
