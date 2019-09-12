#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#############################################################################
# axis
#############################################################################
plotaxis <- function(x, y, at.ticks=NULL, ticklen=NULL, tick.labs=NULL, at.sub=NULL, sublen=NULL) {
  # axis line
  lines(x, y, xpd = 1, col = "gray50")
  
  # ticks and subticks
  if (y[1] == y[2]) {
    # horizontal axis
    # main ticks
    if (!is.null(at.ticks)) {
      if (is.null(tick.labs)) {
        tick.labs <- at.ticks
      }
      for (i in 1:length(at.ticks)) {
        lines(rep(at.ticks[i], 2), c(y[1], y[1]-ticklen), xpd = 1, col = "gray50")
        text(x=at.ticks[i], y=y[1]-ticklen, labels=tick.labs[i], xpd=1, pos=1, cex=0.9, srt=90)
      }
      text(mean(at.ticks), y[1]-ticklen, labels="#k", xpd=1, pos=1)
    }
    # sub ticks
    if (!is.null(at.sub)) {
      for (i in 1:length(at.sub)) {
        lines(x = rep(at.sub[i], 2), y = c(y[1], y[1] - sublen), col = "gray50")
      }
    }
  } else {
    # vertical axis
    # main ticks
    if (!is.null(at.ticks)) {
      if (is.null(tick.labs)) {
        tick.labs <- at.ticks
      }
      for (i in 1:length(at.ticks)) {
        lines(c(x[1]-ticklen, x[1]), rep(at.ticks[i], 2), xpd = 1, col = "gray50")
        text(x=x[1]-ticklen, y=at.ticks[i], labels=tick.labs[i], xpd = 1, pos = 2, cex=0.9, srt=90)
      }
      text(x[1]-ticklen, mean(at.ticks), labels="#  \nkmer", xpd=1, pos=2)
    }
    
    # sub ticks
    if (!is.null(at.sub)) {
      for (i in 1:length(at.sub)) {
        lines(x=c(x[1], x[1]-sublen), y=rep(at.sub[i], 2), col = "gray50")
      }
    }
  }
}


#############################################################################
# number transformation
#############################################################################
numtr <- function(y, ymax, base, yplotmax = 5, ori = c("up", "down")) {
  #' transform number relative to plot position and orientation
  if (sum(ori %in% "up") >= 1) {
    trx <- base + y/ymax * yplotmax
  } else if (ori == "down") {
    trx <- base - y/ymax * yplotmax
  }
  trx
}

#############################################################################
# determine xaxis
#############################################################################
smart.axis <- function(chr_length) {
  numdigits <- nchar(chr_length)
  unit <- 10 ^ (numdigits - 1) / (2- round((chr_length / 10 ^ numdigits), 0)) # 1 or 5 e (numdigits - 1)
  subunit <- unit / 5
  
  numsat <- unit * (0:10)
  numsat <- numsat[numsat < chr_length]
  
  if (unit >= 6) {
    numlabels <- numsat / 1000000
    label.scale <- "Mb"
  } else if (numdigits < 6) {
    numlabels <- unit / 1000
    label.scale <- "Kb"
  }
  
  subunits <- seq(0, chr_length, by = subunit)
  subunits <- subunits[!subunits %in% c(numsat, 1)]
  # return
  list(numsat, numlabels, label.scale, subunits)
}

#############################################################################
# plot distribution
#############################################################################
dist <- function(counts.df, plotcols, ymax, base, yplotmax = 5, ori, color = "lightblue",
                 axislabel="") {
  counts <- rowSums(counts.df[plotcols])
  yvals <- numtr(counts, ori=ori, base = base, ymax = ymax, yplotmax = yplotmax)
  yvals <- c(yvals, base, base)
  xvals <- (counts.df$start + counts.df$end) / 2
  xvals <- c(xvals, max(counts.df$end), 0)
  polygon(y = yvals, x = xvals, col = color, border = NA)
  
  ymax.display.value <- ymax
}

#############################################################################
# kaddistPlot
#############################################################################
kaddistPlot <- function(kdfile, minWins=10, outdir=".") {
  colscheme <- c("tan", "tan3", "slategray2", "slategray4")
  colgroups <- c("overRep", "error", "lowUnderRep", "highUnderRep")
  
  kd <- read.delim(kdfile, stringsAsFactors = F)
  chrlen <- tapply(kd$end, kd$chr, max)
  numwins <- table(kd$chr)
  plotchrs <- names(numwins[numwins >= minWins])
  
  # each chr with at least minWins windows
  par(mar=c(1.5, 0, 3.5, 0.5))
  for (pchr in plotchrs) {
    ckd <- kd[kd$chr == pchr, ]
	# output pdf
    outpdf <- paste0(outdir, "/", pchr, ".kad.dist.pdf")
    pdf(outpdf, width = 6, height = 6)
    curchrlen <- as.numeric(chrlen[pchr]) # length
    # plot
    plot(NULL, NULL, xlim=c(1, curchrlen),
         ylim=c(-1.8, 10.5),
         xlab="", ylab="", bty="n",
         xaxt="n", yaxt="n", main=pchr)
    
    # chromosome
    base0 = 5.5
    rect(xleft = 1, ybottom = base0 - 0.1, xright = chrlen, ytop = base0 + 0.1, col = "gray50", border = NA)
    
    # xaxis
    # yaxis
    smrt.yaxis <- smart.axis(curchrlen)
    xat.coords <- smrt.yaxis[[1]]
    xat.labels <- smrt.yaxis[[2]]
    xat.scale <- smrt.yaxis[[3]]
    xat.subunits <- smrt.yaxis[[4]]
    
    # axis
    axis.base <- -0.75
    lines(c(1, curchrlen), y = rep(axis.base, 2))
    
    # main ticks
    for (i in 1:length(xat.coords)) {
      lines(x = rep(xat.coords[i], 2), y = c(axis.base-0.2, axis.base))
      lines(x = rep(xat.coords[i], 2), y = c(base0-0.2, base0+0.2))
      text(x=xat.coords[i], y = axis.base-0.2, labels = xat.labels[i], xpd=1, pos=1, cex=0.9)
    }
    
    # sub ticks
    for (i in 1:length(xat.subunits)) {
      lines(x = rep(xat.subunits[i], 2), c(axis.base-0.05, axis.base))
      lines(x = rep(xat.subunits[i], 2), c(base0-0.1, base0+0.1))
    }
    # xaxis unit:
    #text(x = 0, y = axis.base, labels = xat.scale, xpd = 1, pos = 2)
    
    xunit <- curchrlen/100
    # error and overRep
    base1 <- 6
    countmax1 <- max(rowSums(ckd[, 6:7]))
    dist(counts.df = ckd, plotcols = c(6, 7), ymax=countmax1, yplotmax=5, base=base1,
         ori = "up", color = colscheme[1])
    dist(counts.df = ckd, plotcols = 6, ymax=countmax1, yplotmax=5, base=base1,
         ori = "up", color = colscheme[2])
    label.vals1 <- c(0, countmax1)
    yaxis.vals1 <- numtr(label.vals1, ymax = countmax1, base = base1, yplotmax = 5, ori = "up")
    plotaxis(x=c(-xunit, -xunit), y=yaxis.vals1, at.ticks=yaxis.vals1, ticklen=xunit, tick.labs=label.vals1)
    
    # underRep
    base2 <- 5
    countmax2 <- max(rowSums(ckd[, 8:9]))
    dist(counts.df = ckd, plotcols = c(8,9), ymax=countmax2, yplotmax=5, base=base2,
         ori = "down", color = colscheme[3])
    dist(counts.df = ckd, plotcols = 9, ymax=countmax2, yplotmax=5, base=base2,
         ori = "down", color = colscheme[4])
    label.vals2 <- c(0, countmax2)
    yaxis.vals2 <- numtr(label.vals2, ymax = countmax2, base = base2, yplotmax = 5, ori = "down")
    plotaxis(x=c(-xunit, -xunit), y=yaxis.vals2, at.ticks=yaxis.vals2, ticklen=xunit, tick.labs=label.vals2)
    text(curchrlen / 2, -0.8, labels = paste0("position (", xat.scale, ")"), cex = 1, xpd = 1, pos=3)
    
    # legends
    legend_base0 <- -2.55
    rect(0, legend_base0-0.15, 2*xunit, legend_base0+0.15, col = colscheme[1], xpd = 1)
    text(2*xunit, legend_base0, labels = colgroups[1], pos = 4, xpd = 1)
      
    rect(0, legend_base0-1.05, 2*xunit, legend_base0-0.75, col = colscheme[2], xpd = 1)
    text(2*xunit, legend_base0-0.9, labels = colgroups[2], pos = 4, xpd = 1)
    
    rect(curchrlen*1/3, legend_base0-0.15, curchrlen*1/3+2*xunit, legend_base0+0.15, col = colscheme[3], xpd = 1)
    text(curchrlen*1/3+2*xunit, legend_base0, labels = colgroups[3], pos = 4, xpd = 1)
    
    rect(curchrlen*1/3, legend_base0-1.05, curchrlen*1/3+2*xunit, legend_base0-0.75, col = colscheme[4], xpd = 1)
    text(curchrlen*1/3+2*xunit, legend_base0-0.9, labels = colgroups[4], pos=4, xpd = 1)
    
    # window size
    winlen <- median(ckd$end - ckd$start + 1) /1000
    text(curchrlen, legend_base0-0.45, labels=paste0("win:", winlen, "kb"), pos=2, xpd = 1)
    # close device
    dev.off()
  }
}

#############################################################################
# main
#############################################################################
kdfile <- args[1]
minWins <- args[2]
outdir <- args[3]
kaddistPlot(kdfile=kdfile, minWins=minWins, outdir=outdir)

