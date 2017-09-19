# BWA-MEM contamination removal

## What does this app do?

This app maps FASTQ reads (paired or unpaired) to a reference genome with the BWA-MEM algorithm, select reads with higher percent identity indicated by user, removed those reads from raw data, and return filtered FASTQ reads to user.

## What are typical use cases for this app?

This app can be used removed contaminants from certain species or reference genome. As part of the microbiome pipeline, this tool is used to removed contamination from human sample. It is suitable as the final filter after user remove low quality bases and adaptors from raw data.

The first step of this app is BWA (Burrows-Wheeler Alignment Tool) software package. BWA includes three algorithms (BWA-backtrack,
BWA-SW and BWA-MEM), but this app is specifically using the BWA-MEM algorithm. The BWA-MEM algorithm is suitable for read lengths ranging from 70bp to 1Mbp. Compared to BWA-backtrack, it has
better performance for 70-100bp Illumina reads; compared to BWA-SW, it is faster and more accurate if the reads
are of high-quality.

The second step is selecting reads with minimum percent identity indicated by user using msamtools. 


## What data are required for this app to run?

This app requires reads files in gzipped FASTQ format (`*.fastq.gz` or `*.fq.gz`), such as those typically produced by Illumina
instruments. 

The app also requires a reference genome sequence index. This must be a gzipped tar archive file (`*.bwa-index.tar.gz`) containing
all the sequence index files as previously output by the BWA indexer. (Indexing is an one-time operation that needs to be performed to a
reference genome sequence in order for it to be usable by BWA. If you have created your own BWA index outside of DNAnexus,
place all the index files in a gzipped tar archive and provide that as the input; if you have a reference genome sequence in FASTA
format, you can index it with the BWA FASTA Indexer app). Some pre-indexed genomes are also available as suggested inputs. See
also ['which human reference sequence should I use?'](https://answers.dnanexus.com/p/183/) on DNAnexus Answers.

If you will be using the mappings to perform variation calling, we encourage you to provide the correct read group information,
and in particular to enter a read group sample. (The default behavior is to add a read group named after the input file, assign
the value `ILLUMINA` to the read group platform, and the value `1` to the read group sample).

## What does this app output?

This app outputs filtered reads, as a paired-end FASTQ files (`*.fastq.gz` or `*.fq.gz`). 

## How does this app work?

This app performs the following steps:

- Mapping with `bwa mem`.
- Retrieve human contamination with  `msamtools -m filter`.
- Remove human read from raw data using custom script  

BWA manual at: http://bio-bwa.sourceforge.net/bwa.shtml
