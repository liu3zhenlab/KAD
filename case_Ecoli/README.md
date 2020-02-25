### A Case study to illustrate KAD analysis
*- Evaluation of *E. coli* genomes with KAD analysis*  

In this case, two versions of *E. coli* genome aseemblies and Illumina read data were downloaded. Analyses for KAD profiling, KAD distribution, and KAD comparison were performed.

#### PART I: Data preparation for KAD analysis

Through this preparation, we will have two genome assemblies for *Escherichia coli* K-12 MG1655
1. U00096.1.fasta
2. U00096.3.fasta

and a set of clean paired-end reads
1. MG1655_1.fq.gz
2. MG1655_2.fq.gz

**data preparation**
1. genome assemblies
For this example, we used two genome assemblies of *Escherichia coli* K-12 MG1655, which were from Genbank accession U00096.1 and U00096.3.

2. Illumina read data
Illumina data for the same bacterial strain were downloaded from NCBI with [sra-tools](https://github.com/ncbi/sra-tools].
```
fasterq-dump --skip-technical --split-3 ERR022075
```

3. trimming of read data with [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
```
########### subject to change #############
trimmomaticDir=<trimmomatic directory>
trimmomatic=$trimmomaticDir/trimmomatic-0.38.jar
adaptor=$trimmomaticDir/adapters/TruSeq3-PE.fa
nthreads=4
read1=ERR022075_1.fastq
read2=ERR022075_2.fastq
prefix=MG1655
minlen=60
###########################################
java -jar $trimmomatic PE \
    -threads $nthreads \
    $read1 $read2 \
    -baseout $prefix \
    ILLUMINACLIP:$adaptor:3:20:10:1:true \
    LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:13 \
    MINLEN:$minlen
```

4. subsample reads with [seqtk](https://github.com/lh3/seqtk) because original read data have >700x sequencing depth. 
```
nreads=2000000
seqtk sample -s 100 MG1655_1P $nreads | gzip > MG1655_1.fq.gz
seqtk sample -s 100 MG1655_2P $nreads | gzip > MG1655_2.fq.gz
```

#### PART II: KAD profiling
KAD profiling is a core step for KAD analysis. It produces k-mer profiles for assembled sequences and reads, combines them into a merged table, determines KAD values, and groups each k-mer based on its KAD value. It generates KAD profiling plots for each assembly.

**Run script**
```
perl <path-to-kad>/KADprofile.pl \
     --read MG1655_1P.gz --read MG1655_2P.gz \
     --asm U00096.1.fasta --aid U00096.1 \
     --asm U00096.3.fasta --aid U00096.3  \
     --prefix MG1655
```
**Outputs in the output directory (MG1655)**  
├── figures  
│   ├── reads.kmer.abundance.profile.pdf  
│   ├── U00096.1.kad.profile.pdf  
│   └── U00096.3.kad.profile.pdf  
├── MG1655_0_reads.kmer.table.txt  
├── MG1655_1_U00096.1.kmer.table.txt  
├── MG1655_1_U00096.3.kmer.table.txt  
├── MG1655_2_merge.kmer.table.txt  
├── MG1655_3_readkmer.cmode  
├── MG1655_4_kad.txt  
├── MG1655_5_readkmer.counts.txt  
├── MG1655_6_bincount.txt  
├── MG1655_7_kad.stat.txt  
├── MG1655_X_Rmd.render.R  
└── report  
    └── MG1655_KADprofile.report.html

The major outputs include
1. KAD result in MG1655_4_kad.txt and here are columns:  
Kmer: k-mer sequence  
reads: abundance of the k-mer in reads  
U00096.1: copy number of the k-mer in the assembly of U00096.1  
U00096.3: copy number of the k-mer in the assembly of U00096.3  
U00096.1.KAD: the KAD value of the k-mer of U00096.1  
U00096.3.KAD: the KAD value of the k-mer of U00096.3  

2. KAD summary in MG1655_7_kad.stat.txt.  
3. A report in HTML format "MG1655_KADprofile.report.html" in the directory of "report"  
4. Editable PDF figures in the directory of "figures"

Here are the KAD profiling plot for the assembly U00096.1:  

<img src="https://github.com/liu3zhenlab/KAD/blob/master/case_Ecoli/figures/U00096.1.kad.profile.png" alt="KAD profiling plot of U00096.1" width="600" />

#### PART III: KAD distribution 
In this step, k-mers are aligned to each original assembly. Based on alignments, distribution of each group of k-mers, particularly error k-mers on assembled sequences is determined. This analysis generates KAD landscape plots.

**Running script for U00096.1**
```
perl <path-to-kad>/KADdist.pl \
     --prefix U00096.1dist \
     --kad ./MG1655/MG1655_4_kad.txt \
     --aid U00096.1 --asm U00096.1.fasta
```
**Outputs in the output directory (U00096.1dist)**  
├── pdf  
│   └── U00096.1.kad.dist.pdf  
├── U00096.1dist_1_kmer.fasta  
├── U00096.1dist_2_kmer.aln  
├── U00096.1dist_3_kmer.kad.bed  
├── U00096.1dist_4_error.pos.bed  
├── U00096.1dist_4_error.pos.merge.bed  
├── U00096.1dist_5_asm.lengths  
├── U00096.1dist_5_kad.bigwig  
├── U00096.1dist_5_kad.wig  
├── U00096.1dist_6_KADtype.dist.txt  
└── U00096.1dist_KADdist.log

The file "U00096.1dist_6_KADtype.dist.txt" contains number of k-mers for each k-mer group (Good, Error, OverRep, LowUnderRep (lounrep), and HighUnderRep (hiunrep), and uncategorized (others)). The data in this file are used to produce the KAD landscape distrition plot "U00096.1.kad.dist.pdf". Exact positions of error k-mers can be found in two BED files "U00096.1dist_4_error.pos.bed" and "U00096.1dist_4_error.pos.merge.bed". The former BED file includes all error k-mers. The latter one has merged data if neighbor error k-mers overlap, which are referred to as error regions.

Here are the landscape plot for the assembly U00096.1:  

<img src="https://github.com/liu3zhenlab/KAD/blob/master/case_Ecoli/figures/U00096.1.kad.dist.png" alt="KAD landscape plot" width="500"/>

**KAD distribution can be generated for U00096.3 using the following script**
```
perl <path-to-kad>/KADdist.pl \
     --prefix U00096.1dist \
     --kad ./MG1655/MG1655_4_kad.txt \
     --aid U00096.3 --asm U00096.3.fasta
```

Here are the landscape plot for the assembly U00096.3:

<img src="https://github.com/liu3zhenlab/KAD/blob/master/case_Ecoli/figures/U00096.3.kad.dist.png" alt="KAD landscape plot" width="500"/>

#### PART IV. KAD comparison between two assemblies, if applied
This analysis compares KAD values for KAD comparison of two assemblies using the output from KAD profiling. K-mers showing distinct KAD values are extracted. The KAD comparison plot visulizes KAD profiles of these k-mers from the two assemblies.

**Run script**
```
perl KADcompare.pl \
     --set1 U00096.1
     --set2 U00096.3
     --prefix U00096.1vs3
     MG1655/MG1655_4_kad.txt
```

**Outputs in the output directory (U00096.1vs3)**  
├── report  
│   └── U00096.1vs3_U00096.1-U00096.3.report.html  
├── U00096.1vs3_1_U00096.1-U00096.3.KADdiff.txt  
├── U00096.1vs3_2_U00096.1-U00096.3.bincount.txt  
├── U00096.1vs3_2_U00096.1-U00096.3.Fig1.raw.pdf  
├── U00096.1vs3_2_U00096.1-U00096.3.Fig2.cuberoot.pdf  
└── U00096.1vs3_X_Rmd.render.R

In the output, "U00096.1vs3_1_U00096.1-U00096.3.KADdiff.txt" lists all k-mers with distinct KAD values of the two assemblies. "U00096.1vs3_2_U00096.1-U00096.3.bincount.txt" include counts of k-mers with distinct KADs in each small interval (bincount) for each of the assembly. Figures are visualized version of the bincount data in the raw and cuberoot scales.

Here are the KAD comparison plot for the two assemblies U00096.1 and U00096.2:  

<img src="https://github.com/liu3zhenlab/KAD/blob/master/case_Ecoli/figures/U00096.1vs3_2_U00096.1-U00096.3.Fig1.raw.png" alt="KAD comparison plot" width="500"/>

**Conclusion**  
The KAD comparison result indicates U00096.3 has less errors as compared to U000096.1. Disagreements between reads and each assembly were found at both assemblies regions that need further examination.

