library(rmarkdown)
library(knitr)

render('/bulk/liu3zhen/research/projects/HTcrispr/crisprDesign/utils/crRNAonGene.Rmd',
  params = list(
    wd="/bulk/liu3zhen/research/projects/HTcrispr/main/1_activation/1b_crispr",
    crdesign="crispr/4.gRNA.design",
    gtf="crispr/1.gtf",
    prefix="crispr",
    nperfect="1",
    offcutoff="10"),
  knit_root_dir=getwd(),
  output_dir=getwd(),
  output_format="html_document",
  output_file="crispr.gRNAdesign.report.html")
