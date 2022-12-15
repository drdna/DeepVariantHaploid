#!/bin/bash

#SBATCH --time 12:00:00
#SBATCH --job-name=bowtie2
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --partition=normal
#SBATCH --mem=100GB
#SBATCH --mail-type ALL
#SBATCH	-A coa_farman_uksr
#SBATCH --mail-type ALL
#SBATCH --mail-user farman@uky.edu

conda activate bioinfo

refpath=$1

readspath=$2

prefix=$3

outdir=$4

singularity -s exec --nv -B /usr/lib/locale/:/usr/lib/locale/ /share/singularity/images/ccs/deepvarient/ccs-deepvariant-1.2.0.sinf /opt/deepvariant/bin/run_deepvariant \
--model_type=WGS \
--ref=$refpath  \
--reads=$readspath \
--output_vcf=${prefix}.vcf.gz \
--intermediate_results_dir $outdir \
--num_shards=16

 
