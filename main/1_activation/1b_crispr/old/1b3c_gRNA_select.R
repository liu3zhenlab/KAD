setwd("/bulk/liu3zhen/research/projects/HTcrispr/main/1_activation/1b_crispr")

##################################################################################################
# modules
##################################################################################################
geneinfo <- function(seqname) {
# determine gene starts
# ZeamMp004::Mt:8601-8751(+)
  genename <- gsub("\\:.*", "", seqname)
  
  is.plus <- grepl("\\+", seqname)
  promoter_range <- gsub(".*\\:", "", seqname)
  promoter_range <- gsub("\\(.*", "", promoter_range)
  promoter_range <- strsplit(promoter_range, "\\-")[[1]]
  promoter_start <- as.numeric(promoter_range[1]) + 1
  promoter_end <- as.numeric(promoter_range[2])
  if (is.plus) {
    genestart <- promoter_end + 1
    geneOrientation <- "+"
  } else {
    genestart <- promoter_start - 1
    geneOrientation <- "-"
  }
  # return: gene name, orientation, start
  c(genename, geneOrientation, genestart)
}


# Zm00001d043937-Zm00001d043937-minus-104_82-3_214207480_214207458
grna_pos <- function(sourceInfo) {
  grnaLoc <- gsub(".*\\-", "", sourceInfo)
  grnaLoc <- strsplit(grnaLoc, "_")[[1]]
  grnaChr <- grnaLoc[1]
  grnaStart <-  as.numeric(grnaLoc[2])
  grnaEnd <-  as.numeric(grnaLoc[3])
  grnaPos <- round((grnaStart + grnaEnd) / 2)
  c(grnaChr, grnaStart, grnaEnd, grnaPos)
}


##################################################################################################
# passed gRNAs
passed <- read.delim("crispr/5.crispr.BsmBI.oligos.txt", comment.char="#")
head(passed)

# all gRNAs
grnaInfo <- read.delim("crispr/4.gRNA.design", comment.char="#")
head(grnaInfo)
grnaInfo <- grnaInfo[grnaInfo$crOligo_name %in% passed$crOligo_name, ]
nrow(grnaInfo)
grnaInfo <- grnaInfo[, c("crOligo_name", "crOligo", "Source")]
grnaInfo$Gene <- gsub("\\:.*", "", grnaInfo$crOligo_name)


grnaLocations <- sapply(grnaInfo$Source, grna_pos)
head(grnaLocations)

grnaInfo$crLocation <- grnaLocations[4, ]

# A188 alignment
a188aln <- read.delim("1b2o_gRNA.blastn.A188Ref1", header=F)
a188aln_filter <- a188aln[a188aln[,3]==100 & a188aln[,4]==23, c(1, 2, 9, 10)]
head(a188aln_filter)
colnames(a188aln_filter) <- c("crOligo_name", "A188chr", "A188start", "A188end")

# merge
passed_grnaInfo <- merge(grnaInfo, a188aln_filter, by="crOligo_name")
head(passed_grnaInfo)
# TATA 
ta <- read.delim("../1a_promoters/1a2o_B73Ref4.ensembl46.promoter.TATA.txt", comment.char="#")
ta.geneinfo <- sapply(ta$Seq_name, geneinfo)
ta.geneinfo  <- data.frame(t(ta.geneinfo))
colnames(ta.geneinfo) <- c("Gene", "GeneOrientation", "GeneStart")
ta2 <- cbind(ta, ta.geneinfo)
head(ta2)
dim(ta2)
ta2$TATApos <- as.numeric(as.character(ta2$GeneStart)) - 150 + round((ta2$Start + ta2$End)/2)
ta2 <- ta2[, c("Gene", "TATApos")]
#ta.geneinfo <- ta.geneinfo[!duplicated(ta.geneinfo$Gene), ]



head(ta.geneinfo)

### final merge
final_grna <- merge(passed_grnaInfo, ta2, by="Gene", all.x=T)


