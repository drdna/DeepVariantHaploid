# DeepVariantHaploid - AI-based variant calling in haploid genomes
## Motivation
Current work in the lab necessitates the use of near-truth SNP datasets. Industry-standard SNP calling pipelines produce an unacceptable number of false calls in fungal genomes (5-30%) due to alignment problems caused by the mutagenic process known as Repeat-Induced Point Mutation (RIP). The goal of this study is to benchmark DeepVariant - an Artifical Intelligence-based variant caller - against [SmartSNPs](), a simple, genome-informed, human intelligence-based strategy that produces near-truth SNP calls.
## Methods
1. Generate BWT and fasta indexes for the reference genome using bowtie 2.2.5:
```bash
mkdir B71v5_index
bowtie2-build B71v5.fasta B71v5_index/B71v5; mv B71v5.fasta B71v5_index
samtools faidx B71v5_index/B71v5.fasta
```
2. Use the [BWT2all_GATK.sh](/scripts/BWT2all_GATK.sh) script (GATK 4.1.4.1) to align reads against the reference and report ALL alignments. Perform SNP calling using standard parameters used for variant calling in fungi:
```bash
sbatch $script/BWT2all_GATK.sh B71v5 /path/to/reads_directory ERR2188722
```
3. Run DeepVariant on the resulting .bam file:
```bash
sbatch $script/DeepVariant B71v5_index/B71v5.fasta B71v5_ERR2188722_ALIGNall DeepVariantCalls DeepVariantTempFiles
```
