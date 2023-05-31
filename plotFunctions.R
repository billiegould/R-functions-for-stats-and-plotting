# author: Chao Dai
# date:   Nov 25, 2020
# description: Generalized plot functions.
library(RColorBrewer)

##' generate a named list matching input classes with randomly selected colors
##' B.Gould 3/2023
##' @param classes a vector of class labels to match with colors
##' 
##' @return a named list matching input classes with random colors
##'
get_random_color_dict <- function(classes){
  n <- length(unique(classes))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  var.colors = unlist(sample(col_vector, n))
  names(var.colors) <- c(unique(classes))
  return(var.colors)
}



##' plot corrleation between two numeric vectors with density color
##' @param x a numeric vector on x axis
##' @param y a numeric vector on y axis

##' @return directly plot without return
##'
plot_colorByDensity <- function(x, y, xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), xlab="", ylab="", main=""){
  df <- data.frame(x,y)
  result <- densCols(x, y, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(result)[1,] + 1L
  cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
  df$col <- cols[df$dens]
  plot(y~x, data=df[order(df$dens),], ylim=ylim, xlim=xlim, pch=20, col=col,
       cex=1, xlab=xlab, ylab=ylab, main=main)
}

##' plot normal and cancer sample separation in MDS
##' @param samples.feature.matrix feature matrix of combined normal and cancer samples
##' @param normal.samples normal samples name
##' @param main figure title
##' @param point.cex point size in MDS
##' @param legendpos legend position
##' @param logDist apply log on sample distance or not 
##' @param addText add cancer sample name to plot or not
##' 
##' @return no value, direct plot
##'
plotMDS <- function(samples.feature.matrix, normal.samples=NULL, main=NULL, point.cex=0.8, legendpos="bottomleft", logDist=F, addText=F, text.cex=1){
  samples <- colnames(samples.feature.matrix)
  if (is.null(normal.samples)){
    normal.indexs <- grep("Normal", samples)
  }else{
    normal.indexs <- which(samples %in% normal.samples)
  }
  samples.relabel <- samples
  samples.relabel[normal.indexs] <- "healthy"
  colnames(samples.feature.matrix) <- samples.relabel
  samples.dist <- dist(t(samples.feature.matrix), method="manhattan") 
  if (logDist == T){
    samples.dist <- log2(1+samples.dist)
  }
  fit <- cmdscale(samples.dist, eig=TRUE, k=2) # k is the number of dim
  x <- fit$points[,1]
  y <- fit$points[,2]
  mixed.samples <- rownames(fit$points)
  normal.indexs <- grep("healthy", mixed.samples)
  highrisk.indexs <- grep("Highrisk", mixed.samples)
  likely.indexs <- grep("Likely", mixed.samples)
  lowrisk.indexs <- grep("Lowrisk", mixed.samples)
  samples.pch <- rep(16, length(mixed.samples))
  samples.color <- rep(NA, length(mixed.samples))
  legend.cancertext <- NULL
  legend.cancercolor <- NULL
  if (length(highrisk.indexs) > 0){
    samples.color[highrisk.indexs] <- "red"
    legend.cancercolor <- c(legend.cancercolor, "red")
    legend.cancertext <- c(legend.cancertext, "highrisk")
  }
  if (length(likely.indexs) > 0){
    samples.color[likely.indexs] <- "blue"
    legend.cancercolor <- c(legend.cancercolor, "blue")
    legend.cancertext <- c(legend.cancertext, "likelyCancer")
  }
  if (length(lowrisk.indexs) > 0){
    samples.color[lowrisk.indexs] <- "green"
    legend.cancercolor <- c(legend.cancercolor, "green")
    legend.cancertext <- c(legend.cancertext, "lowrisk")
  }
  samples.pch[normal.indexs] <- 4
  samples.color[normal.indexs] <- "black"
  plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", cex=point.cex, pch=samples.pch, 
       col=samples.color, main=main)
  legend(legendpos, legend=c("healthy", legend.cancertext), pch=c(4, rep(16, length(legend.cancertext))), bty='n', cex=1,
         col=c("black", legend.cancercolor), horiz=TRUE)
  if(addText){
    healthy.indexs <- grep("healthy", rownames(fit$points))
    indexs <- setdiff(1:nrow(fit$points), healthy.indexs)
    samples.label <- sub(".*risk_", "", rownames(fit$points)[indexs])
    samples.label <- sub("Likely_", "", samples.label)
    text(x[indexs], y[indexs], labels=samples.label, adj=0, pos=1, offset=0.25, cex=text.cex)
  }
}

##' plot bin copy number deviation heatmap, along with arm level copy number deviation 
##' @param copydev.data a data frame of (row, col, colorValue)
##' @param maintitle heatmap main title
##' @param chroms.breakpoints a numeric vector to indicate break points between chromosomes
##' @param qarm.breakpoints a numeric vector to indicate chromosome centromere
##' @param chroms.centerpoints a numeric vector to indicate the centric locations to mark chromosome name
##' @param col.limit a numeric value to indicate heatmap color range limit
##' @param sample.fontsize font size of sample names
##' @param multiCancers.label indicate whether the heatmap data contains multiple cancer types or not
##' 
##' @return no value, direct plot
##'
plot_copyDev_heatmap <- function(copydev.data, maintitle, chroms.breakpoints, qarm.breakpoints, chroms.centerpoints, 
                                 col.limit=0.5, sample.fontsize=6, multiCancers.label=F){
  if (multiCancers.label == F){
    copydev.data$SampleID <- sub(";.*", "", copydev.data$SampleID)
  }else{
    samples <- as.vector(copydev.data$SampleID)
    cancers <- as.vector(sapply(samples, function(x) strsplit(x, split=";")[[1]][2]))
    sampleIDs <- as.vector(sapply(samples, function(x) strsplit(x, split=";")[[1]][1]))
    copydev.data$SampleID <- sampleIDs
    copydev.data$cancer <- cancers
    cancers.unique <- sort(unique(cancers))
    copydev.data <- copydev.data[sort.list(copydev.data[, "cancer"]), ]
    samples.delimit <- NULL
    for (cancer in cancers.unique){
      indexs <- grep(cancer, cancers)
      last.sample <- sampleIDs[indexs[length(indexs)]]
      samples.delimit <- c(samples.delimit, last.sample)
    }
    copydev.data <- copydev.data %>% dplyr::select(-c(cancer))
    temp <- copydev.data %>% pivot_wider(names_from = SampleID, values_from = "value") %>%
      dplyr::select(-c(Regions, arm))
    temp <- t(as.data.frame(temp))
    cancers.breakpoint <- match(samples.delimit, rownames(temp))
    cancers.textpoint <- NULL
    temppoint <- c(0, cancers.breakpoint)
    for (i in 1:(length(temppoint)-1)){
      textpoint <- ceiling((temppoint[i]+temppoint[i+1])/2)
      cancers.textpoint <- c(cancers.textpoint, textpoint)
    }
  }
  copydev.data$value[copydev.data$value >= col.limit] <- col.limit
  copydev.data$value[copydev.data$value <= -col.limit] <- -col.limit
  col.range <- col.limit * c(-1, 1)
  copydev.data$SampleID <- factor(copydev.data$SampleID, levels=unique(copydev.data$SampleID))
  copydev.plot <- ggplot(data=copydev.data, aes(x=Regions, y=SampleID, fill=value)) + geom_raster(hjust=0, vjust=0.5) + 
    ggtitle(maintitle) + labs(fill = "copy number log2 deviation") + 
    scale_fill_distiller(palette="RdBu", limit=col.range) + 
    theme(axis.text.x = element_text(angle=90, hjust=1), 
          axis.text.y = element_text(size=sample.fontsize),
          plot.title = element_text(hjust = 0.5),
          legend.position = "top")
  copydev.plot <- copydev.plot + geom_vline(xintercept = chroms.breakpoints, col = "blue", lty = 3)
  copydev.plot <- copydev.plot + geom_vline(xintercept = qarm.breakpoints, col = "gray", lty = 5, size=0.25)
  chrom.labels <- unique(sub("p|q", "", unique(copydev.data$arm)))
  copydev.plot <- copydev.plot + annotate(geom="text", x=chroms.centerpoints, y=1, label=chrom.labels, color="black")
  if (multiCancers.label == T){
    copydev.plot <- copydev.plot + geom_hline(yintercept = cancers.breakpoint, col = "green", lty = 3)
    copydev.plot <- copydev.plot + annotate(geom="text", x=300, y=cancers.textpoint, label=cancers.unique, color="black", size=1.5)
  }
  copydev.plot <- ggplotly(copydev.plot)
  copydev.plot
}

##' plot arm level copy number deviation heatmap 
##' @param copydev.data a data frame of (sample, arm, colorValue)
##' @param maintitle heatmap main title
##' 
##' @return no value, direct plot
##'
plot_armCopyDev_heatmap <- function(copydev.data, maintitle){
  arms.set1 <- as.vector(sapply(1:12, function(x) c(paste0(x, "p"), paste0(x, "q"))))
  arms.set2 <- c("13q", "14q", "15q")
  arms.set3 <- as.vector(sapply(16:20, function(x) c(paste0(x, "p"), paste0(x, "q"))))
  arms.set4 <- c("21q", "22q", "Xp", "Xq")
  copydev.abs.avg <- copydev.data %>% group_by(arm) %>% 
    summarise(copydevAvg = mean(abs(value)))
  copydev.abs.avg <- as.data.frame(copydev.abs.avg)
  abs.limit <- c(0, max(copydev.abs.avg$copydevAvg))
  arms.level <- intersect(c(arms.set1, arms.set2, arms.set3, arms.set4), copydev.abs.avg$arm)
  copydev.abs.avg$arm <- factor(copydev.abs.avg$arm, levels=arms.level)
  abs.plot <- ggplot(data = copydev.abs.avg, aes(x = arm, y = copydevAvg)) + 
    geom_bar(stat = "identity", aes(fill = copydevAvg)) + theme_gray() + 
    labs(x="Chromosome arm", y="Copy number log2 deviation") +
    scale_fill_distiller(name = "copydevAvg", type = "seq", palette="Greys", limit=abs.limit, direction = 1) + 
    ggtitle(paste(maintitle, "CNV segments")) + theme(plot.title = element_text(hjust = 0.5))
  print(abs.plot)
  
  copydev.gain.avg <- copydev.data %>% group_by(arm) %>% 
    summarise(copydevAvg = mean(value)) %>% filter(copydevAvg > 0) 
  copydev.gain.avg <- as.data.frame(copydev.gain.avg)
  gain.limit <- c(0, max(copydev.gain.avg$copydevAvg))
  arms.level <- intersect(c(arms.set1, arms.set2, arms.set3, arms.set4), copydev.gain.avg$arm)
  copydev.gain.avg$arm <- factor(copydev.gain.avg$arm, levels=arms.level)
  gain.plot <- ggplot(data = copydev.gain.avg, aes(x = arm, y = copydevAvg)) + 
    geom_bar(stat = "identity", aes(fill = copydevAvg)) + theme_gray() + 
    labs(x="Chromosome arm", y="Copy number log2 deviation") +
    scale_fill_distiller(name = "copydevAvg", type = "seq", palette="Oranges", limit=gain.limit, 
                         direction = 1) + 
    ggtitle(paste(maintitle, "CNA segments")) + theme(plot.title = element_text(hjust = 0.5))
  print(gain.plot)
  
  copydev.loss.avg <- copydev.data %>% group_by(arm) %>% 
    summarise(copydevAvg = mean(value)) %>% filter(copydevAvg < 0) 
  copydev.loss.avg <- as.data.frame(copydev.loss.avg)
  loss.limit <- c(min(copydev.loss.avg$copydevAvg), 0)
  arms.level <- intersect(c(arms.set1, arms.set2, arms.set3, arms.set4), copydev.loss.avg$arm)
  copydev.loss.avg$arm <- factor(copydev.loss.avg$arm, levels=arms.level)
  loss.plot <- ggplot(data = copydev.loss.avg, aes(x = arm, y = copydevAvg)) + 
    geom_bar(stat = "identity", aes(fill = copydevAvg)) + theme_gray() + 
    labs(x="Chromosome arm", y="Copy number log2 deviation") +
    scale_fill_distiller(name = "copydevAvg", type = "seq", palette="Greens", limit=loss.limit, 
                         direction = -1) + 
    ggtitle(paste(maintitle, "CNL segments")) + theme(plot.title = element_text(hjust = 0.5))
  print(loss.plot)
}

