#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=24G
#SBATCH --time=1-00:00:00

module load R
list=1i_gene20.1eb.select.txt
ext5=0
ext3=0
source_ref=/homes/liu3zhen/references/maizeCurGenomes/genomes/B73-5.0.fasta
source_gtf=/homes/liu3zhen/references/maizeCurGenomes/gtf/B73-5.0.gtf
te_gff=/homes/liu3zhen/references/NAM1.0/TEs/Zm-B73-REFERENCE-NAM-5.0.TE.gff3
cds_fas=/homes/liu3zhen/references/NAM1.0/verion.b.annotation/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.cds.fa
cdna_fas=/homes/liu3zhen/references/NAM1.0/verion.b.annotation/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.cdna.fa
protein_fas=/homes/liu3zhen/references/NAM1.0/verion.b.annotation/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.fa

# extraction of sequences:
for gene in `cat $list`; do
	echo $gene
	source_seq_dir=$gene
	perl /homes/liu3zhen/scripts2/homotools/geneseq \
		--fas $source_ref \
		--gtf $source_gtf \
		--othergff $te_gff \
		--gene $gene \
		--prefix ${source_seq_dir} \
		--cds $cds_fas \
		--cdna $cdna_fas \
		--prot $protein_fas \
		--ext5 $ext5 \
		--ext3 $ext3
done

