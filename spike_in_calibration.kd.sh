#!/bin/bash
#####################################
#Rewritten in 4/2021 to account for library size as well as spike in normalization
#Rewritten on 2/11/2020 to work with the bash shell. Original scirpt can be foound
#https://raw.githubusercontent.com/Henikoff/Cut-and-Run/master/spike_in_calibration.csh.
#Boolean Control flow was modified.
#####################################

#4. Spike-in calibration
#   NOTE: spike-in calibration may not be appropriate for your sample. You DO need to do it
#   when you need to compare samples, such as comparing treatment to control. The primary
#   genome to spike-in genome ratio per cell is expected to be the same for all samples.
#   The per base pair calculation used here to make a track is:
#      scale * (primary_genome_mapped_count_at_bp)/(spike-in_genome_total_of_mapped_fragments)
#   You can use any scale multiplier and any program you like to do the calculation, the
#   example here uses bedtools genomecov.


if [ "$#" -ne 8 ]; then
    echo "USAGE spike_calibrate.csh genome.bed spike_genome.bed scale output(bg|bga|d) genome_chr_lens min_len max_len"
    echo "Spike-in calibration using bedtools genomecov"
    echo "scale is an arbitrary large number used as multiplier (e.g. 10000)"
    echo "If there is no spike-in, use spike_genome.bed = none"
    echo "min_len and max_len refer to the lengths of fragments in genome.bed to calibrate"
    echo "To calibrate all fragments use min_len = 1 and max_len = 1000"
    echo "Output will be placed in the current directory"
    exit 1
fi

genome_bed=${1}
spike_bed=${2}
scale=${3}
report=${4}
genome_len=${5}
min_len=${6}
max_len=${7}
output_name=${8}

echo "VARIABLES GIVEN"
echo ${genome_bed} ${spike_bed} ${scale} ${report} ${min_len} ${max_len} ${output_name}


echo "Checking if files are real and Non-empty"
if [ -s $genome_bed ] ; then
    echo "${genome_bed} found and not empty"
else
    echo "${genome_bed} was found to either not exist or be empty"
    exit 0
fi


if [ -s $genome_len ] ; then
    echo "$genome_len found and not empty"
else
    echo "${genome_len} was found to either not exist or be empty"
    exit 0
fi

echo "Generating Scale Factor"
echo "Output is in file $output in the current directory"
#MPR is 1000000/total danio reads so total danio reads is getting multiplied by spike in reads
if [ -s $genome_bed ] ; then
        read_count=`wc -l $genome_bed | awk '{print $1}'`
        echo "Danio read count is $read_count"
        MPR=`echo "1000000/$read_count" | bc -l`
        echo "MPR is $MPR"
else
              echo "No Genome File given"
              read_count=0
              MPR=$scale
fi

if [ -s $spike_bed ] ; then
        spike_count=`wc -l $spike_bed | awk '{print $1}'`
        echo "Spike in count is $spike_count"
        spike_norm=`echo "$scale / $spike_count" | bc -l`
        echo "Spike in normalization factor is $spike_norm"
        scale_factor=`echo "$MPR * $spike_norm" | bc -l`
        echo "The scale factor is ${scale_factor}"
else
        echo "No Spike File given"
        spike_count=0
        scale_factor=$scale
fi

#Select fragments within the length range, assumes fragment length is in column 6 of the bed file
echo "Generating Fragment Lengths"
cat $genome_bed | awk -v min=10 -v max=1000 '$7 >= min &&  $7 <= max {print $0}' > frag_length.temp.bed

#Use genomecov to compute spike-in calibrated genome coverage
#see http://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html
#The first position is start-1 and the last is end for bed files, -bg and -bga
#-bga prints zero intervals -bg doesn't, -d prints each bp starting from 1 (not 0)

#This is for IGV which doesn't recognize .bg or .bga files
echo "Generating BG or BGA files"
#if [$report == "bg" || $report == "bga"] ; then
#    output=${name}.${min_len}-${max_len}.bedgraph
#    #echo track type=bedGraph name=$name > $output
#fi

echo "Generating Genome Coverage"
bedtools genomecov -${report} -scale ${scale_factor} -i frag_length.temp.bed -g ${genome_len} > ${output_name}

echo "SCRIPT DONE"
#rm frag_length.temp.bed
