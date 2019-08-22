# KAD
Assessing genome assemblies by comparing k-mer copies in assemblies and K-mer abundance in reads

### Introduction
KAD is designed for evaluating the accuracy of nucleotide base quality of genome assemblies. Briefly, abundance of k-mers are quantified for both sequencing reads and assembly sequences. Comparison of the two values results in a single value per k-mer, K-mer Abundance Difference (KAD), which indicates how well the assembly matches read data for each k-mer.


<img src="https://latex.codecogs.com/svg.latex?\Large&space;KAD=log2(\frac{c+m}{m*(n+1)}>

<img src="https://latex.codecogs.com/svg.latex?\Large&space;x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />

### Citation

### Requirements
The script was written with Perl and R is invoked. Both Perl and R are generally installed. To generate reports, the package [knitr](https://github.com/yihui/knitr) and [rmarkdown](https://rmarkdown.rstudio.com) are needed to be installed in R.

### Installation

### Running guide

### Walk-through example



