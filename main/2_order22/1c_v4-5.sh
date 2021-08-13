#!/bin/bash
grep -f 1i_gene20.list /homes/liu3zhen/references/maizePanGenes/MaizeGDB_B73_pangene_2020_11.tsv | sed 's/.*B73//g' | cut -f 2- > 1i_gene20.v4-5.txt

