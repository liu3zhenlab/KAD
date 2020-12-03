#!/bin/bash
fas=crispr/3.gRNA.fas
blastn -query $fas -db ~/references/A188Ref1/genome/blast+/A188Ref1 -word_size 10  -outfmt 6 -evalue 0.1 > 1b2o_gRNA.blastn.A188Ref1
