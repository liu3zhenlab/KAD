setwd("/bulk/liu3zhen/research/projects/HTcrispr/main/1_activation/1b_13N.crispr")

select_num <- 2
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

# merge gene information with A188 uniquely aligned data:
passed_grnaInfo <- merge(grnaInfo, a188aln_filter, by="crOligo_name")

# filtered good set
good_gRNA <- merge(passed, passed_grnaInfo, by="crOligo_name")
nrow(good_gRNA)
head(good_gRNA)

# final selection:
# at most two guide RNAs per gene
# crRNA oligo with A at the 1st base is preferred
# GC contain (30-80% GC): as.numeric(sapply(as.character(gsub("[AT]", "", good_gRNA$crOligo)), nchar)) / as.numeric(sapply(as.character(good_gRNA$crOligo), nchar))
# repeated bases in a row (<=4 polyhomo):  grepl("A{5}|T{5}|G{5}|C{5}", good_gRNA$crOligo)

# 
good_gRNA$isA <- 1
good_gRNA$isA[grep("^A", good_gRNA$crOligo)] <- 10 # guideRNA with A at the 1st base is preferred due to the compatible desigh with pYPQ141D2.0 using the riceU3 promoter and requiring A at 1st

# GC
gcbases <- as.character(gsub("[AT]", "", good_gRNA$crOligo))
ngc <- as.numeric(sapply(gcbases, nchar))
gcPerc <- ngc / 23

# polyhomo
is.high.polyhomo <- grepl("A{5}|T{5}|G{5}|C{5}", good_gRNA$crOligo)

good_gRNA2 <- good_gRNA[gcPerc > 0.3 & gcPerc < 0.8 & !is.high.polyhomo, ]
nrow(good_gRNA2)

# random selection
genes <- unique(good_gRNA2$Gene)

out <- NULL
for (eg in genes) {
  egdata <- good_gRNA2[good_gRNA2$Gene == eg, ]
  seldata <- egdata[sample(nrow(egdata), min(nrow(egdata), select_num), prob=egdata$isA), ]
  out <- rbind(out, seldata)
}

out <- out[order(out$Order), ]
nrow(out)

### output
write.table(out, "1b3o_regeneration.cand.13Nact.guideRNAs.txt", row.names=F, quote=F, sep="\t")
