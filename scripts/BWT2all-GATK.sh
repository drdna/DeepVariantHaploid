#!/bin/bash

#SBATCH --time 12:00:00
#SBATCH --job-name=bowtie2-GATK
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --partition=CAC48M192_L 
#SBATCH --mem=180GB
#SBATCH --mail-type ALL
#SBATCH	-A col_farman_uksr
#SBATCH --mail-type ALL
#SBATCH --mail-user farman@uky.edu

echo "SLURM_NODELIST: "$SLURM_NODELIST
echo "PWD :" $PWD
echo "index: $1"
echo "path to reads: $2/$3"


#### RUN BOWTIE2

# Usage: bowtie2 <index-dir> <reads-directory> <reads-prefix>

	source /project/farman_uksr/miniconda3/etc/profile.d/conda.sh

# run bowtie2-build

	conda activate py27

	module load ccs/conda/py2.7-samtools-1.10

# index outdir

	readsdir=$2
	index=$1
	indexdir=/project/farman_uksr/${index}_index
#reads 

	reads=$3

#bowtie2-build $ref $indexdir/${ref/.fasta/}

# create output directory

	outdir=${index}_${reads/*\//}_ALIGNall
	mkdir $outdir

# run bowtie2 alignment

if [ "$4" == "SE" ]
then
        bowtie2 -a --threads 16 --very-sensitive-local --phred33 --no-unal -x $indexdir/$index \
	 -U $readsdir/$reads_*\.f*q* 2>$outdir/alignment_summary.txt |  samtools view -bS - > $outdir/accepted_hits.bam
else
        echo "Running bowtie2"
        bowtie2 -a --threads 16 --very-sensitive-local --phred33 --no-unal -x $indexdir/$index \
 	-1 $readsdir/${reads}_*1.f*q* -2 $readsdir/${reads}_*2.f*q* 2>$outdir/alignment_summary.txt | samtools view -bS - > $outdir/accepted_hits.bam
fi


# deleted: --no-unal

#bowtie2 --threads 16 --very-sensitive-local --phred33 -x $indexdir/$index \
# -1 $readsdir/${reads}_*1.f*q* -2 $readsdir/${reads}_*2.f*q* | samtools view -bS - > $outdir/accepted_hits.bam

# this was causing premature failure:  sed 's/#0\/[1-2]//'

# sort bamfile and remove original

	samtools sort -@ 16 -m 8G $outdir/accepted_hits.bam -o $outdir/accepted_hits_sorted.bam
#	rm $outdir/accepted_hits.bam

#conda deactivate	

# add read group info to bam header and remove source file

	o=${reads/*\//}
	conda activate gatk
	java -jar /project/farman_uksr/miniconda3/share/picard-2.21.6-0/picard.jar AddOrReplaceReadGroups I=$outdir/accepted_hits_sorted.bam \
	O=$outdir/accepted_hits_sortedRG.bam RGID=$o RGSM=$o RGLB=$o RGPI=50 RGPL=illumina RGPU=unit1
	conda deactivate

	#rm $outdir/accepted_hits_sorted.bam

# index the bamfile

	conda activate py27
	samtools index $outdir/accepted_hits_sortedRG.bam
	conda deactivate


##### RUN GATK

	source /project/farman_uksr/miniconda3/etc/profile.d/conda.sh

# specify index information

	fasta=${indexdir}/$index.fasta


# give reads prefix (base ID)

	readsprefix=$3



# specify input/output directory

	indir=$outdir


# index reference fasta file

	conda activate py27

	samtools faidx $fasta

	conda deactivate


# activate gatk environment 

	conda activate gatk


# Create dict and fai for reference

# remove existing dictionary

	rm ${fasta/fasta/dict}

#create dictionary

	java -jar /project/farman_uksr/miniconda3/share/picard-2.21.6-0/picard.jar CreateSequenceDictionary R=$fasta O=${fasta/fasta/dict}

# Call haplotype

	/project/farman_uksr/gatk-4.1.4.1/gatk --java-options "-Xmx35g" HaplotypeCaller \
        	--native-pair-hmm-threads 16 \
		-R $fasta \
        	-ploidy 1 \
	        -I $indir/accepted_hits_sortedRG.bam \
        	--emit-ref-confidence GVCF \
	        -O $indir/${readsprefix}.vcf \

# Determine genotypes

	/project/farman_uksr/gatk-4.1.4.1/gatk --java-options "-Xmx35g" GenotypeGVCFs \
	 -R $fasta \
	 -V $indir/$readsprefix.vcf \
	 --output $indir/${readsprefix}_genotyped.vcf


# Filter SNPs only

	/project/farman_uksr/gatk-4.1.4.1/gatk SelectVariants \
	 -R $fasta \
	 -select-type SNP \
	 -V $indir/${readsprefix}_genotyped.vcf \
	 --output $indir/${readsprefix}_genotyped-snps.vcf


	vcftools --vcf $indir/${readsprefix}_genotyped.snp-only.filters.vcf \
	        --remove-filtered-all --recode \
	        --out $indir/${readsprefix}_genotyped.snp-only.filters.subset.vcf

	conda deactivate
