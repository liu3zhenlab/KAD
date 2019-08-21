#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#' Function for calculate cube root value
#' cube root of input numbers (singular or vector) that could be 0, NA,
#' or any positive or negative numbers
cbrt <- function(x) {
  y <- rep(0, length(x))
  not0 <- !is.na(x) & x!=0
  y[is.na(x)] <- NA
  y[not0] <- sign(x[not0]) * exp(log(abs(x[not0])) / 3)
  y
}

#' Function for plotting two KAD curves from k-mers with unequal KADs
diffkadplot <- function(kdfile, transform = c("raw", "cuberoot"),
                        xlabel = "KAD", ylabel = "Counts of k-mers",
                        pmain = "KAD Profiles of unequal KADs",
                        cols = c("red", "darkgreen"), xrange = NULL,
                        yrange = NULL, legendtext = NULL,
                        plotout = FALSE, plotfile = NULL, ...) {
  # check parameters
  stopifnot(sum(transform %in% c("raw", "cuberoot")) >= 1)
  if (length(transform) > 1) { transform = "raw" }
  
  kd <- read.delim(kdfile, stringsAsFactors = F)
  xval <- kd$BIN
  y1 <- kd[, 2]
  y2 <- kd[, 3]
  if (transform == "cuberoot") { 
    y1 <- cbrt(y1)
    y2 <- cbrt(y2)
	ylabel <- paste(ylabel, "(cube root)")
  } else if (max(y1, y2) > 1000000000) {
  	y1 <- y1 / 1000000000
	y2 <- y2 / 1000000000
	ylabel <- paste(ylabel, "(x10^9)")
  } else if (max(y1, y2) > 1000000) {
	y1 <- y1 / 1000000
	y2 <- y2 / 1000000
	ylabel <- paste(ylabel, "(x1,000,000)")
  }
  
  
  # plot
  if (plotout) {
    if (is.null(plotfile)) {
      plotfile <- gsub(".txt$", "", kdfile)
      plotfile <- paste0(plotfile, ".pdf")
    }
    pdf(plotfile, width = 6, height = 5)
  }
  
  # xrange for plotting
  if (is.null(xrange)) {
    xrange <- range(xval)
  }
  
  # yrange for plotting
  if (is.null(yrange)) {
    yrange <- range(c(y1, y2))
  }
  
  # plot
  plot(NULL, NULL, xlim = xrange, ylim = yrange,
       xlab = xlabel, ylab = ylabel, main = pmain, ...)
  lines(xval, y1, col = cols[1], lwd = 1.5)
  lines(xval, y2, col = cols[2], lwd = 1.5)
  
  # legends:
  if (is.null(legendtext)) {
    legendtext <- colnames(kd)[2:3]
  }
  legend("topleft", legend = legendtext, col = cols,
         bty = "n", lty = 1, lwd = 1.2)
  
  if (plotout) { dev.off() }
}

# input arguments
kdfile <- args[1]
outbase <- gsub(".*\\/", "", kdfile)
outbase <- gsub("bincount\\.txt$", "", outbase)
wd <- args[2]

# plot PDFs
pdf1 <- paste0(wd, "/", outbase, "Fig1.raw.pdf")
diffkadplot(kdfile = kdfile, transform = c("raw", "cuberoot"),
            xlabel = "KAD", ylabel = "Counts of k-mers",
            pmain = "KAD Profiles of unequal KADs",
            cols = c("red", "darkgreen"), xrange = NULL,
            yrange = NULL, legendtext = NULL,
            plotout = T, plotfile = pdf1)

pdf2 <- paste0(wd, "", outbase, "Fig2.cuberoot.pdf")
diffkadplot(kdfile = kdfile, transform = "cuberoot",
            xlabel = "KAD", ylabel = "Counts of k-mers",
            pmain = "KAD Profiles of unequal KADs",
            cols = c("red", "darkgreen"), xrange = NULL,
            yrange = NULL, legendtext = NULL,
            plotout = T, plotfile = pdf2)
