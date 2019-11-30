# KAD
Assessing genome assemblies by comparing k-mer copies in assemblies and k-mer abundance in reads

### Introduction
KAD is designed for evaluating the accuracy of nucleotide base quality of genome assemblies. Briefly, abundance of k-mers are quantified for both sequencing reads and assembly sequences. Comparison of the two values results in a single value per k-mer, K-mer Abundance Difference (KAD), which indicates how well the assembly matches read data for each k-mer.


<img src="https://latex.codecogs.com/svg.latex?\Large&space;KAD=log_{2}\begin{pmatrix}\frac{c+m}{m(n+1)}\end{pmatrix}" />

where, _c_ is the count of a k-mer from reads, _m_ is the mode of counts of read k-mers, and _n_ is the copy of the k-mer in the assembly. 

### Requirements
The script was written with Perl and R is invoked. Both Perl and R are generally installed. If needed, please refer to [Perl](https://www.perl.org/) and [R](https://www.r-project.org/) for installation guides. To generate reports, [pandoc](https://pandoc.org) and R packages [knitr](https://github.com/yihui/knitr) and [rmarkdown](https://rmarkdown.rstudio.com) are needed to be installed.

[Jellyfish](https://www.cbcb.umd.edu/software/jellyfish/) was used to generate k-mers by using either FASTA or FASTQ data. The binary Jellyfish executable is included in the [bin](https://github.com/liu3zhenlab/KAD/edit/master/bin/) directory of the KAD package.

To run KADdist.pl, [bowtie](http://bowtie-bio.sourceforge.net/index.shtml) and [bedtools](https://bedtools.readthedocs.io/en/latest/) are required.

### Installation
```
git clone https://github.com/liu3zhenlab/KAD.git  
perl ./KAD/seqKADprofile.pl
```

### conda intallation
Aternatively, all required packages can be installed through conda
```
conda create -n kad
conda activate kad
conda install -c anaconda perl
conda install -c r r-base r-knitr r-rmarkdown
conda install -c bioconda pandoc bowtie bedtools
```

### Data requirements
**1. Read data**  
Illumina sequencing reads with 30x or higher sequence depth. Trimmed clean reads or error corrected reads are preferred. Raw data without trimming are not recommended.

**2. Assembly data**  
Assembly sequencing data in FASTA format. Each assembly is in a single FASTA file.

### Scripts and options
1. [seqKADprofile.pl](seqKADprofile.pl): producing KAD profiles for input assemblies.  
**Usage**: perl seqKADprofile.pl [options]  
**[Options]**  
    --**read** <file>:	\*FASTQ/A read file for k-mer generation; the parameter can be used multiple times; zip files with the suffix of .gz are allowed; required.
    --**minc** <num>:	minimal number of counts per k-mer from reads; k-mers with counts smaller than <num> are not output. default=5.  
    --**asm** <file>:	\*FASTA sequence file for k-mer generation; the parameter can be used multiple times to allow using multiple FASTA file; required.  
    --**rid** <str>:	ID used in the header of the k-mer table generated from reads.  
    --**aid** <str>:	ID used in the header of the asm k-mer table to be generated from each assembly; the parameter can be used multiple times to match --asm input. By default, a header ID is generated from the file name of each assembly by removing PATH and the suffix of .fa, .fas, or .fasta.  
    --**prefix** <str>:the output directory and the prefix for output files; default=kad.  
    --**klen** <num>:  length of k-mers; default=25.  
    --**readdepth** <num>: estimated depth of reads; not required; if specified, it will be compared to the mode of read k-mers.  
    --**kadcutoff** <str of nums>: a set of numbers to define k-mer categories; default="-2 -0.5, 0.5, 0.75, 2".  
    --**binlen** <num>:		bin length to count KAD; default=0.05.  
    --**threads** <num>:		number of cpus; default=1.  
    --**version**:		version  
    --**help**:			help information

2. [KADdist.pl](KADdist.pl): generating distributions of error and other k-mers on contigs or chromosomes.  
**Usage** KADdist.pl [options]  
**[Options]**  
    --**kad**|k <file>:      KAD output file from seqKADprofile.pl; required.  
    --**aid**|i <str>:       assembly ID in the header of KAD file; required.  
    --**asm**|a <file>:      assembly FASTA file, including path; required.  
    --**mincopy**|m <num>: k-mers  with at least --mincopy in the assembly will be aligned to the assembly; default=1.  
    --**maxcopy**|i <num>: k-mers  with at most --maxcopy in the assembly will be aligned to the assembly; default=100.  
    --**winsize**|w <num>: window size on which the number of each KAD type is counted; default=50000.  
    --**kadcutoff**|s <str>: same to --kadcutoff in the seqKADprofile.pl; default="-2 -0.5 0.5 0.75 2".  
    --**prefix**|p <str>:  the output directory and the prefix for output files; default=KADdist.  
    --**minwin4plot**|n <num>: contigs or chromosomes with minimum window number (--minwin4plot) will be plotted; default=10.  
    --**pdfoutdir**|o <str>: the subdirectory under --prefix directory for PDF outputs.  
    --**help**:            help information.  
	
3. [KADcompare.pl](KADcompare.pl): comparing unequal KADs between two assemblies.  
**Usage**: perl KADcompare.pl [options] \<kad\>  
*note: \<kad\> is the output generated by seqKADprofile.pl, which containing KAD values of k-mers.*  
**[Options]**  
    --**set1**:	\*assembly ID 1 in the \<kad\> file; required.  
    --**set2**:	\*assembly ID 2 in the \<kad\> file; required.  
    --**prefix**:	the output directory and the prefix for output filess; default=KC.  
    --**binlen**:	bin length to count KAD; default=0.05.  
    --**force**:	overwrite the existing directory if specified; if not, quit if the output directory exists.  
    --**help**: 	help information

### Examples:
**Analysis 1. KAD profiling**  
This analysis will generate KAD profiles for each input assembly. You will see numbers of total k-mers, error k-mers, and possible missing k-mers of each assembly.

_**how to run**_  
Let us say you have three assembly versions, as shown in the *data* directory:
1. asm0.fas
2. asm1.fas
3. asm2.fas

You also have a read set:
1. read1.fq.gz
2. read2.fq.gz

Assuming the Perl script was in the directory of _scriptpath_, run the following script to generate KAD profiles for all three assemblies.
```
perl scriptpath/seqKADprofile.pl --read read1.fq.gz --read read2.fq.gz \
                                 --asm asm0.fas --asm asm1.fas --asm asm2.fas
```

You might want to assign names of all three assemblies with new names different from file names, such as a0, a1, and a2.
```
perl scriptpath/seqKADprofile.pl --read read1.fq --read read2.fq \
                                 --asm asm0.fas --asm asm1.fas --asm asm2.fas \
                                 --aid a0 --aid a1 --aid a2
```
The number and order of _--aid_ inputs MUST match with _--asm_ inputs.

The parameter _--minc_ might need to change to avoid the interference from a great number of low counts (e.g. 1-3) from error sequences. By default, it is set to 5. However, if high-depth data (e.g., >100x) are generated, the number needs to be increased. Approximately 1/10 of the estimated depth might be a reasonable cutoff.

If corrected reads are used, _--minc_ can be set to a small number (e.g., 3).

```
perl scriptpath/seqKADprofile.pl --read read1.fq --read read2.fq \
                                 --asm asm0.fas --asm asm1.fas --asm asm2.fas \
                                 --aid a0 --aid a1 --aid a2 --minc 15
```
A html report in the _report_ subdirectory is generated from each run. Check this report [example](examples/result_KADprofile.report.pdf).


**Analysis 2. k-mer distribution on contigs or chromosomes of an assembly**  
Based on KAD values of k-mers from Analysis 1, problematic k-mers can be categorized into "error", "overRep", "lowUnderRep", and "highUnderRep", representing k-mers with errors, over-represention, low levels of under-representation, and high levels of under-representation in the assembly. The script [KADdist.pl](KADdist.pl) maps these k-mers to the assembly and combines the KAD value each k-mer to produce:  
1. a [bigwig](https://genome.ucsc.edu/goldenpath/help/bigWig.html) file for visualization  
2. mapping location of error k-mers  
3. plots of distributions of problematic k-mers  

_**how to run**_  
For example, from analysis 1, the assembly *a0* (asm0.fas) was KAD profiled. We now can analyze the ouput file suffixed with "kad.txt" to examine distributions of problematic k-mers. Below is the example script:

```
perl scriptpath/KADdist.pl \ 
--kad a0_4_kad.txt --prefix a0KD \
--aid a0 --asm asm0.fas
```

**Analysis 3. KAD comparison between two assemblies**  
This analysis will directly compare two assemblies based on KADs of the subset of k-mers that have unequal KADs.

After running the analysis using [seqKADprofile.pl](seqKADprofile.pl), KAD values are generated. In the [examples](https://github.com/liu3zhenlab/KAD/tree/master/examples) directory, the file **result_4_kad.txt** contains KAD values. We now can select any two assemblies in this file to compare.

_**how to run**_  
Assuming again the Perl script was in the directory of _scriptpath_, the following run compares a0 with a2. Note that the input _--set1_ and _--set2_ should match the assembly names used in the KAD file.

```
perl scriptpath/KADcompare.pl --set1 a0 --set2 a2 --prefix a0_2 result_4_kad.txt
```
A html report in the _report_ subdirectory is generated from each run. Check this report [example](examples/a0_2_a0-a2.report.pdf).

**Notes**: Here are what analysis 3 does:  
First, the script extracts k-mers with unequal copies in the two assemblies. Two KADs per k-mer of the two assemblies are therefore different. Of two KADs per k-mer, one KAD may be NA because zero count of the k-mer from both reads and the assembly. These NAs are converted to 0 due to the agreement between reads and assembly data.

Second, the script counts KAD per defined bin, which, by default, is 0.05. Separate counts per bin of two assemblies are used for visualizing the two KAD profiles of k-mers with unequal KADs.
