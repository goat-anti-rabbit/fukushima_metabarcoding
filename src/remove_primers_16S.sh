#!/bin/bash

# Define paths
input_dir="data/raw/16S_reads"                # Input FASTQ files
output_dir="data/intermediate/16S_trimmed"    # Output directory for primer-removed files
primer_f="CCTACGGGNGGCWGCAG"                  # Replace with your forward primer sequence
primer_r="GACTACHVGGGTATCTAATCC"              # Replace with your reverse primer sequence

# Create output directory if it doesn't exist
mkdir -p $output_dir

# Loop through FASTQ pairs
for fq1 in ${input_dir}/*_1.fq.gz; do
  fq2=${fq1/_1.fq.gz/_2.fq.gz}  # Match the corresponding reverse read file
  sample_name=$(basename $fq1 _1.fq.gz)

  # Run Cutadapt
  cutadapt \
    -g ^$primer_f -G ^$primer_r \
    -o ${output_dir}/${sample_name}_1_trimmed.fq.gz \
    -p ${output_dir}/${sample_name}_2_trimmed.fq.gz \
    $fq1 $fq2
done
