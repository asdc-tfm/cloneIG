#!/bin/bash

help="cloneIG.sh -i config_csv -o outdir -t total_threads -m mixcr_threads"
cores=$(($(nproc) - 1))
user_path=${PWD%/*} #parent dir

while [ -n "$1" ]; do # while there are arguments
    case "$1" in
    -i) input=$2
        shift;;
    -o) output=$2
        shift;;
    -t) threads=$2
        shift;;
    -m) mixcr_threads=$2
        shift;;
    -h) echo $help ;;
    --) shift
        break ;;
    *) echo "Option $1 not recognized" ;;
    esac
    shift
done

# exit if no config_csv
if [[ ! $input  ]]
then
    echo $help
    exit
fi

if [[ $threads -gt $cores ]]
then
    echo -e "$threads is higher than number of available threads: $cores."
    exit
fi

if [[ $threads == '' ]]
then
    threads=1
fi


if [[ $mixcr_threads -gt $threads ]]
then
    echo -e "$mixcr_threads is higher than number of selected threads: $threads."
    exit
fi

if [[ $mixcr_threads == '' ]]
then
    mixcr_threads=1
fi


currentDate=$(date)
echo -e "Starting analysis... $currentDate\n"


if [[ $output == '' ]]
then
    output=$PWD
fi

# create config.yaml
outfile=config.yaml

echo -e "# path to config csv" > $outfile
echo -e "samples: ${input}\n" >> $outfile

echo -e "# path to outdir" >> $outfile
echo -e "outdir: ${output}/results\n" >> $outfile

echo -e "# path to logdir" >> $outfile
echo -e "logdir: ${output}/logs\n" >> $outfile

echo -e "resources:" >> $outfile

echo -e "    default:" >> $outfile
echo -e "        threads: 1" >> $outfile
echo -e "        mem: 4096" >> $outfile
echo -e "        walltime: 10" >> $outfile

echo -e "    mixcr:" >> $outfile
echo -e "        threads: ${mixcr_threads}" >> $outfile
echo -e "        mem: 8192" >> $outfile
echo -e "        walltime: 1440" >> $outfile

echo -e "    vidjil:" >> $outfile
echo -e "        threads: 1" >> $outfile
echo -e "        mem: 8192" >> $outfile
echo -e "        walltime: 1440" >> $outfile 

# create config CSV file
echo "Patient;Condition;Sample;R1;R2" > $input
echo -e "MM129;BM;MM129_BM;${user_path}/data/MM129/MM129_R1_BM.fastq.gz;${user_path}/data/MM129/MM129_R2_BM.fastq.gz" >> $input
echo -e "MM129;PPC;MM129_PPC;${user_path}/data/MM129/MM129_R1_PPC.fastq.gz;${user_path}/data/MM129/MM129_R2_PPC.fastq.gz" >> $input
echo -e "MM129;cfDNA;MM129_cfDNA;${user_path}/data/MM129/MM129_R1_cfDNA.fastq.gz;${user_path}/data/MM129/MM129_R2_cfDNA.fastq.gz" >> $input

# snakemake execution
snakemake --use-conda -c $threads

currentDate=$(date)

echo -e "\nProcess finished...$currentDate"

exit
