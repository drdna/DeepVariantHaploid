# DeepVariantHaploid - AI-based variant calling in haploid genomes
## Motivation
Current work in the lab necessitates the use of near-truth SNP datasets. Industry-standard SNP calling pipelines produce an unacceptable number of false calls in fungal genomes (5-30%) due to alignment problems caused by the mutagenic process known as Repeat-Induced Point Mutation (RIP). The goal of this study is to benchmark DeepVariant - an Artifical Intelligence-based variant caller - against [SmartSNPs](https://github.com/drdna/SmartSNPs), a simple, genome-informed, human intelligence-based strategy that produces near-truth SNP datasets.
## Methods
1. Generate BWT and fasta indexes for the reference genome using bowtie 2.2.5:
```bash
mkdir B71v5_index
bowtie2-build B71v5.fasta B71v5_index/B71v5; mv B71v5.fasta B71v5_index
samtools faidx B71v5_index/B71v5.fasta
```
2. Use the [BWT2-GATK.sh](/scripts/BWT2-GATK.sh) script (GATK 4.1.4.1) to align reads against the reference and report BEST alignments. Perform SNP calling using standard parameters used for variant calling in fungi:
```bash
sbatch $script/BWT2-GATK.sh B71v5 /path/to/reads_directory ERR2188722
```
Resulting VCF file: [ERR2188722_best_genotyped-snps.vcf](/data/ERR2188722_best_genotyped-snps.vcf)

3. Use the [BWT2all-GATK.sh](/scripts/BWT2all-GATK.sh) script (GATK 4.1.4.1) to align reads against the reference and report ALL alignments. Perform SNP calling using standard parameters used for variant calling in fungi:
```bash
sbatch $script/BWT2all-GATK.sh B71v5 /path/to/reads_directory ERR2188722
```
4. Run [DeepVariant.sh](/scripts/DeepVariant.sh) on the resulting .bam files:
```bash
sbatch $script/DeepVariant B71v5_index/B71v5.fasta B71v5_ERR2188722_ALIGN/accepted_hits_sortedRG.bam DeepVariantBest DeepVariantBestTemp
sbatch $script/DeepVariant B71v5_index/B71v5.fasta B71v5_ERR2188722_ALIGNall/accepted_hits_sortedRG.bam DeepVariantAll DeepVariantAllTemp
```
5. Extract SNP calls:
```bash
zgrep PASS DeepVariantAllAlign.vcf.gz | awk '$1 ~ /^##/ || (length($4) == 1 && length($5) == 1' | gzip - > DeepVariantAll.vcf.gz
zgrep PASS DeepVariantBestAlign.vcf.gz | awk '$1 ~ /^##/ || (length($4) == 1 && length($5) == 1)' | gzip - > DeepVariantBest.vcf.gz
```
Resulting VCF files: [DeepVariantAll.vcf.gz](/data/DeepVariantAll.vcf.gz), [DeepVariantBest.vcf.gz](/data/DeepVariantBest.vcf.gz)

6. Visualize html reports:
### Variants found with "best" alignments
![DeepVariantBest_report.tiff](data/DeepVariantBest_report.tiff)
### Variants found with "all" alignments
![DeepVariantAll_report.tiff](data/DeepVariantAll_report.tiff)

## Generate "truth" dataset:
1. Use SmartSNPsV2.pl to filter "best" alignments VCF file:
```bash
gunzip DeepVariantAllAlign.vcf.gz
perl SmartSNPsV2.pl DeepVariantAllAlign.vcf B71v5_align/B71v5.B71v5_alignments 20 10
```
Runtime report is as follows: 
#NumRecords: 7980; Allowed: 406; Repeated: 4958; Non-repeat heterozygotes: 2209; Low coverage: 407.

2. Extract valid variants:
```bash
awk '$1 ~ /^##/ || (length($4) == 1 && length($5) == 1)' ERR2188722_genotyped-snps_SSfilter.vcf | grep -v FAIL > SmartSNPsTruthSet.vcf
```
### SmartSNPs "Truth" dataset:

[SmartSNPsTruthSet.vcf](/data/SmartSNPsTruthSet.vcf)

Manual filtering resulted in the removal of XX false variants identifiable because several occurred within very short distances of one another. These SNPs likely come from short blocks of non-orthologous sequence in the query genome.

### GATK reported 6,830 SNPs (best alignment) of which only 406 are valid (true variants), with 4,958 calls being rejected because they came from repeats in the reference genome, and 2,209 rejected as they were derived from repeats in the query genome.

## A total of 4,823 SNPs passed the DeepVariant filter for the best alignment dataset. Only 2,450 of these were shared with the GATK dataset.

## DeepVariant makes inexplicable errors in calling homozygotes versus heterozygotes:
### Hets called as homos (just a tiny set of examples):
chr1	2288874	.	C	T	27	PASS	.	GT:GQ:DP:AD:VAF:PL	1/1:25:106:16,89:0.839623:26,28,0

chr2	3203907	.	A	T	17.7	PASS	.	GT:GQ:DP:AD:VAF:PL	1/1:14:201:62,135:0.671642:17,16,0

chr3	1799382	.	A	T	7.4	PASS	.	GT:GQ:DP:AD:VAF:PL	1/1:4:206:40,66:0.320388:5,4,0

chr3	4495178	.	T	C	17.3	PASS	.	GT:GQ:DP:AD:VAF:PL	1/1:9:139:69,61:0.438849:16,8,0

chr3	7895961	.	A	G	9.2	PASS	.	GT:GQ:DP:AD:VAF:PL	1/1:9:146:57,87:0.59589:8,22,0

chr3	7895978	.	C	A	4.2	PASS	.	GT:GQ:DP:AD:VAF:PL	1/1:4:155:53,95:0.612903:1,12,0

### Homos called as hets (much less common):

chr2	912475	.	C	T	14.2	PASS	.	GT:GQ:DP:AD:VAF:PL	0/1:7:60:0,56:0.933333:13,0,7

chr2	4883063	.	G	T	17.7	PASS	.	GT:GQ:DP:AD:VAF:PL	0/1:10:131:0,131:1:17,0,10






