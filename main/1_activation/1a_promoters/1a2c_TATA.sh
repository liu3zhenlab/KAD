perl ~/scripts/fasta/pattern.search.pl -I 1a1o_B73Ref4.ensembl46.promoter.fasta -P TAT[AT]A[AT] -O F > 1a2o_B73Ref4.ensembl46.promoter.TATA.txt
cat 1a2o_B73Ref4.ensembl46.promoter.TATA.txt | sed 's/\:.*//g' | grep "^Z" | sort | uniq | wc -l
cat 1a2o_B73Ref4.ensembl46.promoter.TATA.txt | cut -f 6 | grep TAT | sort | uniq -c

