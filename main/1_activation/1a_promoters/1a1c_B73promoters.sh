#!/bin/bash
gofbed=/homes/liu3zhen/references/B73Ref4/genemodel2/B73Ref4.ensembl46.gene.bed
clen=/homes/liu3zhen/references/B73Ref4/genome/B73Ref4.length
ref=/homes/liu3zhen/references/B73Ref4/genome/B73Ref4.fa
out=1a1o_B73Ref4.ensembl46
promoter_len=150
bedtools flank -i $gofbed -g $clen -l $promoter_len -r 0 -s > ${out}.promoter.bed
bedtools getfasta -fi $ref -bed ${out}.promoter.bed -fo ${out}.promoter.fasta -s -name+

# convert bed to gtf
#1	.	promoter	46229	46342	.	+	.	gene_id "Zm00001d027230";
awk '{ print $1"\t.\tpromoter\t"$2+1"\t"$3"\t.\t"$6"\t.\tgene_id \""$4"\";" }' ${out}.promoter.bed > ${out}.promoter.gtf

