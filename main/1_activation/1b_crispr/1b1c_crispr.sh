#!/bin/bash
ref=/homes/liu3zhen/references/B73Ref4/genome/B73Ref4.fa
bowtieDB=~/references/B73Ref4/bowtie/B73Ref4
gtf=../1a_promoters/1a1o_B73Ref4.ensembl46.promoter.gtf
genelist=1b1i_regeneration.candidate.genes 

#from Changtian “No requirement of A or G for the first nucleotide of sgRNA because of ZmUbi is employing instead of U3 or U6 promoters”
# therefore, --pattern was change from the default (G[ATGC]{20}GG) to [ATGC]{21}GG
perl ../../../crisprDesign/crispr.scan.pl --fasta $ref \
	--bowtie_db $bowtieDB \
	--gtf $gtf \
	--glist $genelist \
	--feature promoter \
	--pattern [ATGC]{21}GG \
	--force
 
