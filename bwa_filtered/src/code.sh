#!/bin/bash
#
sleep 2h
# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#
# Fetch and uncompress genome
#
mkdir genome
dx cat "$genomeindex_targz" | tar zxvf - -C genome  # => genome/<ref>, genome/<ref>.ann, genome/<ref>.bwt, etc.
genome_file=`ls genome/*.bwt`     # Locate a file called <ref>.bwt
genome_file="${genome_file%.bwt}" # Remove the bwt suffix to keep the <ref>

mkdir -p ~/out/reads_out_fastqgz/ ~/out/reads2_out_fastqgz/ ~/out/removal_list/

#
# Fetch reads
#
dx-download-all-inputs --except genomeindex_targz --parallel

#
# Set up options
#
opts=""
if [ "$all_alignments" == "true" ]; then
  opts="$opts -a"
fi
if [ "$mark_as_secondary" == "true" ]; then
  opts="$opts -M"
fi
if [ "$add_read_group" == "true" ]; then
  opts="$opts -R @RG\\tID:${read_group_id}\\tPL:${read_group_platform}\\tPU:${read_group_platform_unit}\\tLB:${read_group_library}\\tSM:${read_group_sample}"
fi
if [ "$advanced_options" != "" ]; then
  opts="$advanced_options"
fi

input="./in/reads_fastqgz/* ./in/reads2_fastqgz/*"

sudo chmod +x /usr/bin/msamtools
bwa mem -t `nproc` "$genome_file" $input $opts | msamtools -m filter -p "$percent_identity" -z 50 -l 50 -S -  | awk 'BEGIN {FS="\t"}; {print $1}' > pre_removallist.txt 

cat pre_removallist.txt | uniq > removallist.txt

name="$reads_fastqgz_prefix"
# Remove any _R1 / _1 suffixes
name="${name%_1}"
name="${name%_R1}"
#gunzip ~/in/reads_fastqgz/*
#gunzip ~/in/reads2_fastqgz/*
( zcat < ./in/reads_fastqgz/* | python remove_match_paired_read.py removallist.txt | gzip > out/reads_out_fastqgz/"${name}.1.fq.gz" ) &
forward_pid=$!
zcat < ./in/reads2_fastqgz/* | python remove_match_paired_read.py removallist.txt | gzip > out/reads2_out_fastqgz/"${name}.2.fq.gz"
wait "$forward_pid"
mv removallist.txt out/removal_list/
#cd ~/out/reads_out_fastqgz
#gzip "${name%}_1.fastq" &
#cd ~/out/reads2_out_fastqgz
#gzip "${name%}_2.fastq"
#wait 
dx-upload-all-outputs --parallel
