# KAD
Assessing genome assemblies by comparing k-mer copies in assemblies and K-mer abundance in reads

### Introduction
KAD is designed for evaluating the accuracy of nucleotide base quality of genome assemblies. Briefly, abundance of k-mers are quantified for both sequencing reads and assembly sequences. Comparison of the two values results in a single value per k-mer, K-mer Abundance Difference (KAD), which indicates how well the assembly matches read data for each k-mer.


<img src="https://latex.codecogs.com/svg.latex?\Large&space;KAD=log2\frac{c+m}{m*(n+1)}" />

where, _c_ is the count of a k-mer from reads, _m_ the mode of counts of read k-mers, _n_ is the copy of the k-mer in the assembly. 

### Citation

### Requirements
The script was written with Perl and R is invoked. Both Perl and R are generally installed. To generate reports, the R packages [knitr](https://github.com/yihui/knitr) and [rmarkdown](https://rmarkdown.rstudio.com) are needed to be installed.

### Installation
git clone https://github.com/liu3zhenlab/KAD.git  
perl ./KAD/seqKADprofile.pl

### Running guide

### Walk-through example



