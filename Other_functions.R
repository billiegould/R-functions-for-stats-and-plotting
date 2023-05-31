
##' plotSampleRelation plot the sample relationship in hierarchical cluster or 2D MDS plot (from Bioconductor lumi package)
##' @param x a matrix with each column corresponding to a sample
##' @param subset the subset probes used to determine the sample relations. If it is one number, then randomly selected "number" of probes will be used. If not provide, all the probes will be used.
##' @param cv.Th the threshold of the coefficient of variance of probes used to select probes to estimate sample relations
##' @param standardize standardize the expression profiles or not 
##' @param method "MDS" or "hierarchical clustering"
##' @param dimension the principle components to visualize the MDS plot
##' @param color the color for each sample during plot. Only support the "mds" method
##' @param main the title of the plot
##' @param pch use symbols instead of text to label the samples
##' @param addLegend Whether to add legend to MDS (two-dimensional PCA) plot
##' @param \dots Other parameters used by plot function.
##' 
##' @return Invisibly return the hierarchical clustering results (if 'cluster' method used) or coordinates of the mds plot (if 'mds' method used) 
plotSampleRelation <- function(x, subset=NULL, cv.Th=0.1, standardize=TRUE, method=c('cluster', 'mds'), dimension=c(1,2), color=NULL, main=NULL, pch=NULL, addLegend=TRUE, ...) {
  if (is(x, 'ExpressionSet')) {
    dataMatrix <- exprs(x)
  } else if (is.matrix(x)) {
    dataMatrix <- x
  } else {
    stop('The class of "x" should be matrix or LumiBatch!')
  }
  
  ## Standardize each sample
  if (standardize) dataMatrix <- scale(dataMatrix)
  
  if (is.null(subset)) {
    ## Filter the genes with most of the experiments "Absent"
    probeList <- rownames(dataMatrix)
    if (is.null(probeList)) probeList <- 1:nrow(dataMatrix)
    if (cv.Th > 0) {
      cv.gene <- apply(dataMatrix, 1, function(x) sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE))
      subset <- probeList[abs(cv.gene) > cv.Th]
      if (is.null(main)) main <- paste('Sample relations based on', length(subset), 'genes with sd/mean >', cv.Th)
    } else {
      subset <- probeList
      if (is.null(main)) main <- paste('Sample relations based on', length(subset), 'genes')
    }
  } else {
    if (length(subset) == 1 && is.numeric(subset)) {
      subset <- sample(1:nrow(dataMatrix), min(subset, nrow(dataMatrix)), replace=TRUE)
    }
    if (is.null(main)) main <- paste('Sample relations based on', length(subset), 'selected genes')
  }
  
  dd <- dist(t(dataMatrix[subset,]))
  method <- match.arg(method)
  if (method == 'cluster') {
    hc = hclust(dd, 'ave')
    plot(hc, xlab='Sample', main=main, ...)
    attr(hc, 'geneNum') <- length(subset)
    attr(hc, 'threshold') <- cv.Th
    return(invisible(hc))	
  } else {
    ## Multi-Dimension Scaling
    mds.result <- cmdscale(dd, k=max(dimension), eig=TRUE)
    ppoints <- mds.result$points
    eig <- mds.result$eig
    percent <- round(eig/sum(eig) * 100, 1)
    
    colorLegend <- NULL
    if (is.null(color)) {
      color <- 1
    } else {
      if (!is.numeric(color)) {
        allColor <- colors()
        if (!all(is.element(color, allColor))) {
          colorLegend <- unique(color)
          color <- as.numeric(factor(color, levels=colorLegend))
        } 
      }
    }
    if (missing(pch)) {
      plot(ppoints[,dimension[1]], ppoints[,dimension[2]], type='n', 
           xlab=paste('Principal Component ', dimension[1], " (", percent[dimension[1]], "%)", sep=""),
           ylab=paste('Principal Component ', dimension[2], " (", percent[dimension[2]], "%)", sep=""), 
           main=main, ...)
      text(ppoints[,dimension[1]], ppoints[,dimension[2]], col=color, labels=colnames(dataMatrix), cex=1)
    } else {
      plot(ppoints[,dimension[1]], ppoints[,dimension[2]],  
           xlab=paste('Principal Component ', dimension[1], " (", percent[dimension[1]], "%)", sep=""),
           ylab=paste('Principal Component ', dimension[2], " (", percent[dimension[2]], "%)", sep=""), 
           main=main, col=color, pch=pch, ...)
    }
    attr(ppoints, 'geneNum') <- length(subset)
    attr(ppoints, 'threshold') <- cv.Th
    
    ## add legend if color is a factor
    if (!is.null(colorLegend) && addLegend) {
      if (!missing(pch)) {
        legend('topleft', legend=colorLegend, col=unique(color), pch=unique(pch))
      } else {
        legend('topleft', legend=colorLegend, text.col=1:length(colorLegend))
      }
    }
    
    return(invisible(ppoints))	
  }
}


##' plotSampleRelationSNP plot sample relationship based on SNP variants
##' @param SNPInfo a VRanges object keeps the SNP variant information
##' @param SNPInfo2 an optional VRanges object2, which keeps the SNP variant information from another group of samples
##' @param suffix the suffix attached to the end of sample labels when SNPInfo2 is provided
##' @param plotAF whether plot based on VariantFreq or binary
##' @param plotCol colum names to be plot
##' @param mdsPlot whether plot 2D MDS plot or hierarchical cluster
##' @param title the title of the plot
##' 
##' @return invisibly return the matrix used for ploting
##' 
plotSampleRelationSNP <- function(SNPInfo, SNPInfo2=NULL, suffix=NULL, plotAF=TRUE, plotCol='VariantFreq', mdsPlot=FALSE, title=NULL, ...) {
  
  if (!is(SNPInfo, 'VRanges')) SNPInfo <- asVRanges(SNPInfo)
  if (is.null(SNPInfo$VariantID)) SNPInfo$VariantID <- getVariantIDs(SNPInfo)
  if (plotCol == 'VariantFreq') {
    if (is.null(SNPInfo$VariantFreq)) SNPInfo$VariantFreq <- altFraction(SNPInfo) * 100
  } else {
    if (!(plotCol %in% names(values(SNPInfo)))) stop(paste(plotCol, 'not in SNPInfo!'))
    SNPInfo$VariantFreq <- values(SNPInfo)[, plotCol]
    plotAF <- TRUE
  }
  
  if (!is.null(SNPInfo2)) {
    if (!is(SNPInfo2, 'VRanges')) SNPInfo2 <- asVRanges(SNPInfo2)
    if (is.null(SNPInfo2$VariantID)) SNPInfo2$VariantID <- getVariantIDs(SNPInfo2)
    if (is.null(SNPInfo2$VariantFreq)) SNPInfo2$VariantFreq <- altFraction(SNPInfo2) * 100
    
    commNames <- names(values(SNPInfo))[names(values(SNPInfo)) %in% names(values(SNPInfo2))]
    values(SNPInfo) <- values(SNPInfo)[,commNames]
    values(SNPInfo2) <- values(SNPInfo2)[,commNames]
    ## add suffix to differentiate two groups and avoid overlapping sample names
    if (!is.null(suffix)) {
      if (length(suffix) == 1) suffix <- rep(suffix, 2)
    } else {
      suffix <- c('1', '2')
    }
    sampleNames(SNPInfo) <- paste0(sampleNames(SNPInfo), '_', suffix[1])
    sampleNames(SNPInfo2) <- paste0(sampleNames(SNPInfo2), '_', suffix[2])
    SNPInfo <- c(SNPInfo, SNPInfo2)
  }
  
  varBySample <- table(SNPInfo$VariantID, as.character(sampleNames(SNPInfo)))
  afBySample <- varBySample
  afBySample[cbind(SNPInfo$VariantID, as.character(sampleNames(SNPInfo)))] <- SNPInfo$VariantFreq
  ## remove all NA samples
  all.na.ind <- apply(afBySample, 2, function(x) all(is.na(x)))
  afBySample <- afBySample[,!all.na.ind,drop=FALSE]
  
  if (mdsPlot) {
    if (plotAF) {
      if (is.null(title)) {
        title <- 'Sample relations based on SNP AFs'
      }
      plotSampleRelation(afBySample, standardize = FALSE, method='mds', cv.Th=0, main=title, ...)
      #plotSampleRelation(afBySample, standardize = FALSE, method='mds', cv.Th=0, dimension=c(2,3))
    } else {
      if (is.null(title)) {
        title <- 'Sample relations based on SNP (0,1)'
      }
      plotSampleRelation(varBySample, standardize = FALSE, method='mds', cv.Th=0, main=title, ...)
    }
  } else {
    if (plotAF) {
      if (is.null(title)) {
        title <- 'Sample relations based on SNP AFs'
      }
      hc <- plotSampleRelation(afBySample, standardize = FALSE, method='cluster', cv.Th = 0, main=title, ...)
    } else {
      if (is.null(title)) {
        title <- 'Sample relations based on SNP (0,1)'
      }
      hc <- plotSampleRelation(varBySample, standardize = FALSE, method='cluster', cv.Th=0, main=title, ...)
    }
  }
  
  if (mdsPlot) {
    if (plotAF) {
      return(invisible(afBySample))
    } else {
      return(invisible(varBySample))
    }
  } else {
    if (plotAF) {
      return(invisible(afBySample[,hc$order]))
    } else {
      return(invisible(varBySample[,hc$order]))
    }
  }
}


##' plotSampleRelationSNP_BC plot the sample relationship in hierarchical cluster based on sample SNP barcodes
##' @param SNP_BC a named vector of SNP barcodes
##' @param SNP_BC2 an optional  named vector of SNP barcodes, which are from another group of samples
##' @param suffix the suffix attached to the end of sample labels when SNPInfo2 is provided
##' @param plotCol colum names to be plot
##' @param labelCol colum names to be plot as sample labels
##' @param mdsPlot whether plot 2D MDS plot or hierarchical cluster
##' @param title the title of the plot
##' 
##' @return invisibly return the matrix used for ploting
##' 
plotSampleRelationSNP_BC <- function(SNP_BC, SNP_BC2=NULL, suffix=NULL, plotCol='SNP_BC', labelCol='SampleID', mdsPlot=FALSE, title=NULL, ...) {
  
  if (is.data.frame(SNP_BC)) {
    if (all(c(plotCol, labelCol) %in% names(SNP_BC))) {
      sampleIDs <- SNP_BC[, labelCol]
      SNP_BC <- SNP_BC[, plotCol]
    } else {
      stop('plotCol and labelCol need to be provided!')
    }
  } else {
    sampleIDs <- names(SNP_BC)
  }
  sampleIDs <- sampleIDs[!is.na(SNP_BC)]
  SNP_BC <- SNP_BC[!is.na(SNP_BC)]
  if (length(unique(nchar(SNP_BC))) > 1) {
    warning('SNP_BC does not have equal lengths!')
    return(FALSE)
  }
  SNP_BC_known <- gsub('[-_?]', '', SNP_BC)
  known.ratio <- nchar(SNP_BC_known)/nchar(SNP_BC)
  SNP_BC <- SNP_BC[known.ratio > 0.5]
  sampleIDs <- sampleIDs[known.ratio > 0.5]
  
  if (length(SNP_BC) < 3) {
    warning('Minimum of three SNP_BC are needed for plotting!')
    return(FALSE)
  }
  if (is.null(sampleIDs)) stop('names of SNP_BC is required for labeling!')
  SNP_BC <- as.character(SNP_BC)
  names(SNP_BC) <- sampleIDs
  SNP_BC.value <- strsplit(SNP_BC, split='')
  SNP_BC.value <- lapply(SNP_BC.value, function(x) {
    x[x == '0'] <- '-1'
    x[x == '1'] <- '0'
    x[x == '2'] <- '1'
    suppressWarnings(x <- as.numeric(x))
    return(x)
  })
  varBySample <- do.call('cbind', SNP_BC.value)
  
  if (!is.null(SNP_BC2)) {
    if (is.data.frame(SNP_BC2)) {
      if (all(c(plotCol, labelCol) %in% names(SNP_BC2))) {
        sampleIDs.2 <- SNP_BC2[, labelCol]
        SNP_BC2 <- SNP_BC2[, plotCol]
      } else {
        stop('plotCol and labelCol need to be provided!')
      }
    } else {
      sampleIDs.2 <- names(SNP_BC2)
    }
    if (is.null(sampleIDs.2)) stop('names of SNP_BC is required for labeling!')
    SNP_BC2 <- as.character(SNP_BC2)
    names(SNP_BC2) <- sampleIDs.2
    SNP_BC2.value <- strsplit(SNP_BC2, split='')
    SNP_BC2.value <- lapply(SNP_BC2.value, function(x) {
      x[x == '0'] <- '-1'
      x[x == '1'] <- '0'
      x[x == '2'] <- '1'
      suppressWarnings(x <- as.numeric(x))
      return(x)
    })
    varBySample2 <- do.call('cbind', SNP_BC2.value)
    
    if (!is.null(suffix)) {
      if (length(suffix) == 1) suffix <- rep(suffix, 2)
    } else {
      suffix <- c('1', '2')
    }
    
    colnames(varBySample) <- paste0(colnames(varBySample), '_', suffix[1])
    colnames(varBySample2) <- paste0(colnames(varBySample2), '_', suffix[2])
    varBySample <- cbind(varBySample, varBySample2)
  }
  ## remove all NA samples
  all.na.ind <- apply(varBySample, 2, function(x) all(is.na(x)))
  varBySample <- varBySample[,!all.na.ind,drop=FALSE]
  
  # varBySample <- varBySample[!is.na(rowSums(varBySample)),,drop=FALSE]
  if (mdsPlot) {
    if (is.null(title)) {
      title <- 'Sample relations based on SNP barcodes'
    }
    plotSampleRelation(varBySample, standardize = FALSE, method='mds', cv.Th=0, main=title, ...)
  } else {
    if (is.null(title)) {
      title <- 'Sample relations based on SNP barcodes'
    }
    hc <- plotSampleRelation(varBySample, standardize = FALSE, method='cluster', cv.Th=0, main=title, ...)
  }
  
  if (mdsPlot) {
    return(invisible(varBySample))
  } else {
    return(invisible(varBySample[,hc$order]))
  }
}


##' findFileBySample summarize the bySample output files output by pipeline
##' @param sampleIDs a vector of Predicine sample IDs or a data.frame with SampleID column
##' @param file.suffix the character suffix used for filename match in OutputBySample folder, support files in consensus result folder and bamfile. Return the sample consensus result folder, if file.suffix is NULL.
##' @param mergedBam whether search files based on the results of mergedBam file
##' @param returnAllMatch whether return all matches of sampleIDs (mutiple runs may exist for the same sample)
##' @param returnRealPath whether return real path instead of symlinks (if exists)
##' @param sampleDirs outputBySample directory or individual directories for each sample
##' @param pversion the pipeline version used to mapped in outputBySample folder
##' 
##' @return a vector of file links
##' 
findFileBySample <- function(sampleIDs, file.suffix="_summary.*csv", mergedBam=FALSE, returnAllMatch=FALSE, returnRealPath=FALSE,
                             sampleDirs=NULL, pversion=NULL) {
  
  if (is.null(pversion)) {
    pversion <- ''
  } else {
    pversion <- unique(c(pversion, ''))
  }
  if (is.null(sampleDirs)) {
    sampleDirs <- options('outputBySampleDir')[[1]]
  } 
  if (length(sampleDirs) == 1) {
    # sampleDirs <- rep(sampleDirs, length(sampleIDs))
    sampleDirInfo <- findSampleDirs(sampleIDs, outputBySampleDir=sampleDirs, returnAllMatch=returnAllMatch, pversion=pversion) 
    sampleDirs <- sampleDirInfo$sampleDirs
    missingDirs <- sampleDirInfo$missingDirs
    sampleIDs <- sampleIDs[sampleIDs %in% names(sampleDirs)]
  } else if (length(sampleDirs) != length(sampleIDs)) {
    stop('sampleDirs does not match sampleIDs!')
  } else if (length(sampleDirs) != length(sampleIDs)) {
    if (any(!dir.exists(sampleDirs))) {
      missingDirs <- sampleDirs[!dir.exists(sampleDirs)]
      warning(paste("cannot find directory: ", paste(missingDirs, collapse = ',')))
      sampleIDs <- sampleIDs[dir.exists(sampleDirs)]
      sampleDirs <- sampleDirs[dir.exists(sampleDirs)]
    }
  }
  if (length(sampleDirs) == length(sampleIDs) && length(sampleIDs) > 0) names(sampleDirs) <- sampleIDs
  
  bam.suffix <- ifelse(mergedBam, '_merged', '_consensus')
  missingSamples <- mapSamples <- filelist <- NULL
  if (length(sampleIDs) == 0)
    return(list(filelist=filelist, mapSamples=mapSamples, missingSamples=missingSamples))
  
  for (j in 1:length(sampleIDs)) {
    sample.j <- sampleIDs[j]
    sampleDir.j <- sampleDirs[[sample.j]]
    
    for (k in 1:length(sampleDir.j)) {
      sampleDir.jk <- sampleDir.j[k]
      sampleName.jk <- basename(sampleDir.jk)
      if (is.null(file.suffix)) {
        file.jk <- file.path(sampleDir.jk, paste0(basename(sampleDir.jk), bam.suffix))
      } else {
        if (grepl(bam.suffix, basename(sampleDir.jk))) {
          file.jk <- dir(file.path(sampleDir.jk, 'results'), pattern=file.suffix, full.names=T)
        } else {
          file.jk <- dir(sampleDir.jk, pattern=file.suffix, full.names=T)
        }
        if (length(file.jk) == 0) {
          if (grepl('bam', file.suffix)) {
            file.jk <- dir(file.path(sampleDir.jk, 'bamfiles'), pattern=file.suffix, full.names=T)
          } else if (grepl('fastq', file.suffix)) {
            file.jk <- dir(file.path(sampleDir.jk, 'fastq'), pattern=file.suffix, full.names=T)
          } else {
            file.jk <- dir(file.path(sampleDir.jk, paste0(sampleName.jk, bam.suffix)), pattern=file.suffix, full.names=T)
          }
          if (length(file.jk) == 0) {
            file.jk <- ''
          } 
        }
        if (grepl('fastq', file.suffix)) {
          warning('fastq files are not fully supported yet!')
          if(length(file.jk) > 2) {
            warning(paste('More than one files were found, only the first one will be used for', sampleName.jk, file.suffix))
            file.jk <- file.jk[1:2]
          }  
        } else {
          if(length(file.jk) > 1) {
            warning(paste('More than one files were found, only the first one will be used for', sampleName.jk, file.suffix))
            file.jk <- file.jk[1]
          }
        }
      }
      
      if (!file.exists(file.jk[1])) {
        missingSamples <- c(missingSamples, sampleDir.jk)
        next
      } else {
        filelist <- c(filelist, file.jk)
        mapSamples <- c(mapSamples, rep(sample.j, length(file.jk)))
      }
    }
  }
  if (returnRealPath) {
    filelist <- Sys.readlink(filelist)
  }
  # names(filelist) <- mapSamples
  return(list(filelist=filelist, mapSamples=mapSamples, missingSamples=missingSamples))
}


##' getVariantBySample get variant by sample
##' @param sampleIDs a vector of Predicine sample IDs or a data.frame with SampleID column
##' @param variantIDs a list of variantIDs (e.g., chr22:42526562:42526562:G:C) or GRanges object for individual SampleID used to retrive variants. return all variants if it is NULL.
##' @param file.suffix the character suffix used for filename match in OutputBySample folder
##' @param AF.th minimum VariantFreq threshold to filter variants
##' @param altCount.th minimum supporting variant count to filter variants
##' @param returnAllMatch whether return all matches of sampleIDs (mutiple runs may exist for the same sample)
##' @param returnWildType return WildType variant coverage and SNP call in varCoverage
##' @param sampleDirs outputBySample directory or individual directories for each sample
##' @param pversion pipeline version. If provided, only results using matched pipeline version will be returned
##' 
##' @return a list with following elements: varInfo, varCoverage, filelist, missingSamples and missingVariants
##' 
getVariantBySample <- function(sampleIDs, variantIDs=NULL, file.suffix="_general.csv", AF.th=0, altCount.th=0, returnAllMatch=FALSE, returnWildType=FALSE,
                               sampleDirs=NULL, pversion=NULL) {
  
  if (is.null(sampleDirs)) {
    sampleDirs <- options('outputBySampleDir')[[1]]
  } 
  if (length(sampleDirs) == 1) {
    sampleDirs <- rep(sampleDirs, length(sampleIDs))
  } else if (length(sampleDirs) != length(sampleIDs)) {
    stop('sampleDirs does not match sampleIDs!')
  }
  if (!any(dir.exists(sampleDirs))) {
    sampleDirs[!dir.exists(sampleDirs)] <- file.path(options('outputBySampleDir')[[1]], sampleDirs[!dir.exists(sampleDirs)])
    if (!any(dir.exists(sampleDirs))) stop(paste("cannot find directory: ", paste(sampleDirs, collapse = ',')))
  }
  
  varInfo.flat <- varCoverage <- varInfo.all <- missingVariants <- allfiles <- missingSamples <- missingVarInfo <- commNames <- varCoverage <- NULL
  if (length(sampleIDs) == 0)
    return(list(varInfo=varInfo.flat, varCoverage=varCoverage, filelist=allfiles, missingSamples=missingSamples, missingVarInfo=missingVarInfo))
  
  for (j in 1:length(sampleIDs)) {
    sample.j <- sampleIDs[j]
    fileInfo.j <- findFileBySample(sample.j, file.suffix=file.suffix, returnAllMatch=returnAllMatch, pversion=pversion, sampleDirs=sampleDirs[j])
    filelist.j <- fileInfo.j$filelist
    missingSamples <- c(missingSamples, fileInfo.j$missingSamples)
    if (length(filelist.j) == 0) {
      cat(paste('Sample', sample.j, 'is missing!\n'))
      missingSamples <- c(missingSamples, sample.j)
      next
    }
    
    if (!is.null(variantIDs)) {
      variantIDs.j <- variantIDs
      if (is.list(variantIDs)) {
        variantIDs.j <- variantIDs[[sample.j]]
      }
      variantIDs.j <- checkVariantID(variantIDs.j)
    } else {
      variantIDs.j <- NULL
    }
    ## get files
    for (file.i in filelist.j) {
      if (grepl('RData$', file.i, ignore.case = TRUE)) {
        tt <- load(file.i); varInfo.i <- get(tt); rm(list=tt)
      } else if (grepl('rds$', file.i, ignore.case = TRUE)) {
        varInfo.i <- readRDS(file.i)
      } else {
        varInfo.i <- read.csv(file.i, as.is=TRUE, check.names=FALSE, fill=TRUE, comment='#', colClasses=list(ref='character', alt='character'), stringsAsFactors=FALSE)
        if (nrow(varInfo.i) == 0) {
          missingVariants <- c(missingVariants, list(variantIDs.j))
          next
        }
        varInfo.i <- varInfo.i[, colnames(varInfo.i) != '']
        varInfo.i <- asVRanges(varInfo.i)
      }
      sampleNames(varInfo.i) <- sample.j
      varInfo.i$VariantID <- getVariantIDs(varInfo.i)
      
      missingVarInfo.j <- NULL
      varCoverage.j <- NULL
      if (!is.null(variantIDs.j)) {
        if (is(variantIDs.j, 'GRanges')) {
          suppressWarnings(varInfo.i <- varInfo.i[overlapsAny(varInfo.i, variantIDs.j)])
          suppressWarnings(missingVariants.i <- variantIDs.j[!overlapsAny(variantIDs.j, varInfo.i)])
        } else {
          varInfo.i <- varInfo.i[varInfo.i$VariantID %in% variantIDs.j]
          missingVariants.i <- variantIDs.j[!(variantIDs.j %in% varInfo.i$VariantID)]
        }
        
        ## check coverage to know whether it is well covered
        if (returnWildType && !is(variantIDs.j, 'GRanges')) {
          covfile.j <- findFileBySample(sample.j, file.suffix="coverage.RData", sampleDirs=sampleDirs[j])
          if (!is.null(covfile.j$filelist)) {
            varCoverage.j <- checkVariantCoverage(variantIDs.j, varInfo.i, coverageInfo=covfile.j$filelist)
          }
        } else {
          missingVarInfo.j <- missingVariants.i
        }
      } 
      if (is.null(varInfo.i$VariantFreq) || any(is.na(varInfo.i$VariantFreq))) varInfo.i$VariantFreq <- signif(altFraction(varInfo.i), 3) * 100
      if (length(varInfo.i) > 0) {
        varInfo.i$SampleDir <- dirname(file.i)
        sampleNames(varInfo.i) <- as.character(sampleNames(varInfo.i))
        varInfo.i$SampleID <- sample.j
        if (is.null(commNames)) {
          commNames <- colnames(values(varInfo.i))
        } else {
          commNames <- intersect(commNames, colnames(values(varInfo.i)))
        }
        varInfo.i <- varInfo.i[altDepth(varInfo.i) >= altCount.th & varInfo.i$VariantFreq >= AF.th]
        varInfo.all <- c(varInfo.all, varInfo.i)
      }
      missingVarInfo <- c(missingVarInfo, list(missingVarInfo.j))
      varCoverage <- c(varCoverage, list(varCoverage.j))
    }
    allfiles <- c(allfiles, filelist.j)
  }
  tmp <- lapply(varInfo.all, function(x) {
    values(x) <- values(x)[,commNames]
    return(x)
  })
  varInfo.flat <- unlist(GRangesList(tmp))
  if (!is.null(variantIDs) && length(allfiles) > 0) {
    names(missingVarInfo) <- basename(allfiles)
  } else {
    missingVarInfo <- NULL
  }
  
  return(list(varInfo=varInfo.flat, varCoverage=varCoverage, filelist=allfiles, missingSamples=missingSamples, missingVarInfo=missingVarInfo))
}


##' findProjRunLoc find Project and Run location information by matching keyword in project names
##' @param keywords keywords used to search project names
##' @param runRootDir the run root directory used by DeepSea pipeline
##' @param RNAProject whether it is RNA project or DNA project
##' @param expand whether to separate different project results in different rows or combine them by comma
##' @return a data frame with run and project columns
##' 
findProjRunLoc <- function(keywords, runRootDir=NULL, RNAProject=FALSE, expand=TRUE) {
  
  if (is.null(runRootDir)) runRootDir <- options('runRootDir')[[1]]
  if (!dir.exists(runRootDir)) stop(paste("cannot find directory: ", runRootDir))
  
  runlist <- dir(runRootDir, pattern='^[0-9]+', include.dirs = TRUE, full.names = TRUE)
  runProjInfo <- NULL
  for (run.i in runlist) {
    
    lbwf.dir.i <- dir(run.i, pattern='^(lbwfresult)|(dsrun)', full.names = TRUE)
    if (length(lbwf.dir.i) > 0) {
      # nextflow <- TRUE
      # pversion <- gsub('(dsrun)|(lbwfresult)', '', basename(lbwf.dir.i))
      lbwf.dir.i <- dir(lbwf.dir.i, pattern='^(lbwfresult)|(dsrun)', full.names = TRUE)
    } else {
      lbwf.dir.i <- run.i
    }
    
    if (RNAProject) {
      projlist.i <- list.dirs(file.path(lbwf.dir.i, 'rna'), recursive = FALSE)
    } else {
      projlist.i <- list.dirs(lbwf.dir.i, recursive = FALSE)
    }
    keywords <- paste0('(', paste(keywords, collapse=')|(') ,')')
    proj.i <- projlist.i[grep(keywords, basename(projlist.i), ignore.case = TRUE)]
    if (length(proj.i) > 0) {
      if (length(proj.i) > 1) {
        if (!expand) {
          proj.i <- paste(proj.i, collapse = ',')
        } else {
          runProjInfo.i <- cbind(rep(run.i, length(proj.i)), proj.i)
        }
      } else {
        runProjInfo.i <- c(run.i, proj.i)
      }
      if (!is.null(runProjInfo)) {
        runProjInfo <- rbind(runProjInfo, runProjInfo.i)
      } else {
        runProjInfo <- runProjInfo.i
      }
    }
  }
  if (!is.null(runProjInfo)) {
    if (!is.matrix(runProjInfo)) {
      runProjInfo <- matrix(runProjInfo, nrow=1)
    }
    colnames(runProjInfo) <- c('run', 'project')
    rownames(runProjInfo) <- NULL
  }
  return(runProjInfo)
}


##' findSampleRunLoc find Project and Run location information by matching keyword in project names
##' @param SampleID SampleIDs used to search project names
##' @param runRootDir the run root directory used by DeepSea pipeline
##' @param outputBySampleDir the outputBySampleDir folder
##' @param pversion pipeline version, e.g., '1.7.0', '1.6.2'
##' 
##' @return a data frame of sampleInfo wiht run column
##' 
findSampleRunLoc <- function(sampleIDs, runRootDir=NULL, outputBySampleDir=NULL, pversion=NULL) {
  
  pattern.ids <- paste0('(', paste(sampleIDs, collapse=')|('), ')')
  if (is.null(runRootDir)) runRootDir <- ''
  
  runlist <- NULL
  if (file.exists(runRootDir)) {
    runlist <- dir(runRootDir, pattern='^[0-9]+', include.dirs = TRUE, full.names = TRUE)
  } else {
    if (is.null(outputBySampleDir) || !file.exists(outputBySampleDir)) {
      outputBySampleDir <- options('outputBySampleDir')[[1]]
    }
    if (!is.null(pversion)) {
      outputBySampleDir <- file.path(outputBySampleDir, pversion)
    }
    sampledirlist <- dir(outputBySampleDir, full.names=TRUE)
    sampledirlist <- sampledirlist[grep(pattern.ids, basename(sampledirlist))]
    for (sampledir.i in sampledirlist) {
      run.i <- dir(sampledir.i, pattern='^[0-9]+_', include.dirs = TRUE, full.names=TRUE)
      run.i <- Sys.readlink(run.i)
      runlist <- c(run.i, runlist)
    }
  }
  
  sampleInfo.all <- NULL
  for (run.i in runlist) {
    run.org.i <- run.i
    if (grepl('(lbwfresult)|(dsrun)', basename(run.i))) {
      run.i <- dirname(run.i)
      if (grepl('(lbwfresult)|(dsrun)', basename(run.i))) {
        run.i <- dirname(run.i)
      }
    }
    samplesheet.i <- file.path(run.i, 'sample.tab')
    if (!file.exists(samplesheet.i)) next
    sampleInfo.i <- read.table(samplesheet.i, sep='\t', head=TRUE, as.is=TRUE, fill=TRUE)
    sampleInfo.i <- sampleInfo.i[!duplicated(sampleInfo.i$sample_id),,drop=FALSE]
    rownames(sampleInfo.i) <- sampleInfo.i$sample_id
    selInd.i <- grep(pattern.ids, sampleInfo.i$sample_id)
    if (length(selInd.i) > 0) {
      sampleInfo.i$run <- run.org.i
      sampleInfo.all <- combineGRanges(sampleInfo.all, sampleInfo.i[selInd.i,,drop=FALSE])
    }
  }
  
  return(sampleInfo.all)
}


##' parseRunInfoFile parse the RunInfoFile output by DeepSea R pipeline
##' @param runInfoFile the text file output by evaluation_pipeline_bam_sh
##' 
##' @return a list of paramters: panel, dataType
##' 
parseRunInfoFile <- function(runInfoFile) {
  suppressWarnings(runInfo <- readLines(runInfoFile))
  panel <- runInfo[grep('^Panel', runInfo, ignore.case = TRUE)]
  panel <- gsub('Panel[: ]+', '', panel, ignore.case=TRUE)
  dataType <- runInfo[grep('dataType', runInfo, ignore.case = TRUE)]
  if (length(dataType) > 0)
    dataType <- gsub('.*dataType[:= ]+', '', dataType[1], ignore.case=TRUE)
  params <- list(panel=panel, dataType=dataType)
  return(params)
}


##' summarizeBySample summarize the bySample output files output by pipeline
##' @param sampleIDs a vector of Predicine sample IDs or a data.frame with SampleID column
##' @param qc.suffix the character suffix of NGS QC file used for filename match
##' @param cnv.suffix the character suffix of CNV file used for filename match
##' @param fusion.suffix the character suffix of Fusion file used for filename match
##' @param variant.suffix the character suffix of Variant file used for filename match
##' @param tmb.suffix the character suffix of TMB Variant file used for filename match
##' @param snp.suffix the character suffix of TMB Variant file used for filename match
##' @param otherVariant.suffix the character suffix of other variant file (like dist2end and general) used for filename match
##' @param addQC_Cols additional lab QC columns from samplesheet
##' @param sampleInfo a data.frame keeps the sample information, which includes SampleID, PatientID and Specimen_type/SampleType
##' @param otherSampleCols other columns in sampleInfo will be included in the combined output files
##' @param pversion pipeline version. If provided, results using matched pipeline version will be returned first
##' @param removeDuplicate whether remove duplicated variants/CNV based on variantIDs
##' @param returnAllMatch whether return all matches of sampleIDs (mutiple runs may exist for the same sample)
##' @param simpleReturn whether only return simplified variants results
##' @param sampleDirs outputBySample directory or individual directories for each sample
##' @param outputDir the location of output file, if it is NULL, not files will be saved
##' 
##' @details the function organize the pipeline output results by samples.  
##'          If sampleInfo is provided with PatientID, it will also check variant concordant.
##'          If SpecimenType/SampleType is also provided in sampleInfo, it will further check germline variants (if any SpecimenType is pbmc, normal or germline)
##' @return a list with following elements: QC, Variant, TMB, CNV, Fusion and missingFileInfo
##' 
summarizeBySample <- function(sampleIDs, qc.suffix="_summary.*csv", cnv.suffix="_CNV.*csv",
                              fusion.suffix="_fusion.*csv", variant.suffix='_variants.csv', variant.short.suffix='_variants_short.csv',
                              tmb.suffix='_TMB.*csv', snp.suffix='_annotated_SNP.csv',  otherVariant.suffix=NULL,
                              sampleInfo=NULL, otherSampleCols=NULL, pversion=NULL, removeDuplicate=FALSE, returnAllMatch=FALSE, 
                              simpleReturn=FALSE, sampleDirs=NULL, outputDir=NULL, verbose=FALSE) {
  
  if (is.null(pversion)) {
    pversion <- ''
  } else {
    pversion <- unique(c(pversion, ''))
  }
  if (is.data.frame(sampleIDs)) {
    if (is.null(sampleIDs$SampleID)) stop('SampleID is required!')
    sampleIDs <- sampleIDs$SampleID
  }
  if (length(sampleDirs) == 1) {
    sampleIDs <- unique(as.character(sampleIDs))
  }
  simpleVarCols <- c("sampleNames", "seqnames", "start", "end", "ref", "alt", "totalDepth", "refDepth", "altDepth",
                     "VariantFreq", "SYMBOL", "HGVSc", "HGVSp", "VariantID", "sampleFolder")
  
  if (!is.null(sampleInfo)) {
    if (is.null(sampleInfo$SampleID)) stop('SampleID is required!')
    if (any(duplicated(sampleInfo[,'SampleID']))) {
      warning(paste('Duplicated SampleIDs:', paste(sampleInfo[duplicated(sampleInfo[,'SampleID']), 'SampleID'], collapse=',')))
      sampleInfo <- sampleInfo[!duplicated(sampleInfo[,'SampleID']),,drop=FALSE]
    }
    rownames(sampleInfo) <- sampleInfo[,'SampleID']
    if (is.null(sampleInfo$PatientID)) {
      if (!is.null(sampleInfo$External.Patient.ID)) {
        sampleInfo$PatientID <- sampleInfo$External.Patient.ID
      } else if (!is.null(sampleInfo$"External Patient ID")) {
        sampleInfo$PatientID <- sampleInfo$"External Patient ID"
      } else if (!is.null(sampleInfo$"External_Patient_ID")) {
        sampleInfo$PatientID <- sampleInfo$"External_Patient_ID"
      }
    }
    if (is.null(sampleInfo$PatientID)) stop('PatientID is required!')
    sampleInfo$PatientID <- as.character(sampleInfo$PatientID)
    sampleInfo$SampleID <- as.character(sampleInfo$SampleID)
    if (!is.null(sampleInfo$SampleType)) {
      sampleInfo$SpecimenType <- sampleInfo$SampleType
    } else if (!is.null(sampleInfo$Specimen_Type)) {
      sampleInfo$SpecimenType <- sampleInfo$Specimen_Type
    } else if (!is.null(sampleInfo$Specimen_type)) {
      sampleInfo$SpecimenType <- sampleInfo$Specimen_type
    } else if (!is.null(sampleInfo$specimenType)) {
      sampleInfo$SpecimenType <- sampleInfo$specimenType
    } 
    if (any(otherSampleCols %in% c('Specimen_type', 'specimenType', 'Specimen_Type', 'SampleType'))) {
      otherSampleCols <- otherSampleCols[!(otherSampleCols %in% c('Specimen_type', 'specimenType', 'Specimen_Type', 'SampleType'))]
      otherSampleCols <- c(otherSampleCols, 'SpecimenType')
    }
  }
  
  if (is.null(sampleDirs)) {
    sampleDirs <- options('outputBySampleDir')[[1]]
  } 
  missingDirs <- NULL
  if (length(sampleDirs) == 1) {
    # sampleDirs <- rep(sampleDirs, length(sampleIDs))
    sampleDirInfo <- findSampleDirs(sampleIDs, outputBySampleDir=sampleDirs, returnAllMatch=returnAllMatch, pversion=pversion) 
    sampleDirs <- sampleDirInfo$sampleDirs
    missingDirs <- sampleDirInfo$missingDirs
    sampleIDs <- sampleIDs[sampleIDs %in% names(sampleDirs)]
  } else if (length(sampleDirs) != length(sampleIDs)) {
    stop('sampleDirs does not match sampleIDs!')
  } else if (length(sampleDirs) != length(sampleIDs)) {
    if (any(!dir.exists(sampleDirs))) {
      missingDirs <- sampleDirs[!dir.exists(sampleDirs)]
      warning(paste("cannot find directory: ", paste(missingDirs, collapse = ',')))
      sampleIDs <- sampleIDs[dir.exists(sampleDirs)]
      sampleDirs <- sampleDirs[dir.exists(sampleDirs)]
    }
  }
  if (length(sampleDirs) == length(sampleIDs) && length(sampleIDs) > 0) names(sampleDirs) <- sampleIDs
  if (length(sampleIDs) == 0) {
    return(c(missingFileInfo=list(missingDirs)))
  }
  
  suffixList <- c(NGSQC=qc.suffix, Variant=variant.suffix, Variant_short=variant.short.suffix, TMB=tmb.suffix, SNP=snp.suffix, 
                  otherVariant=otherVariant.suffix, CNV=cnv.suffix, Fusion=fusion.suffix)
  combinedInfoList <- NULL
  missingFileInfo <- missingDirs
  for (i in 1:length(suffixList)) {
    step.i <- names(suffixList)[i]
    if (verbose) cat(paste('Processing', step.i, '\n'))
    suffix.i <- suffixList[i]
    if (is.null(suffix.i) || suffix.i == '') {
      combinedInfoList <- c(combinedInfoList, list(NULL))
      next
    }
    combinedInfo <- data.frame()
    excludeCols <- NULL
    for (j in 1:length(sampleIDs)) {
      sample.j <- sampleIDs[j]
      if (verbose) cat(paste('Processing', sample.j, '\n'))
      sampleDir.j <- sampleDirs[[sample.j]]
      if (length(sampleDir.j) == 0) {
        missingFileInfo <- c(missingFileInfo, paste('missing sampleDir', sample.j))
      } else {
        for (k in 1:length(sampleDir.j)) {
          sampleDir.jk <- sampleDir.j[k]
          sampleName.jk <- basename(sampleDir.jk)
          file.jk <- dir(sampleDir.jk, pattern=suffix.i, full.names=T)
          if (length(file.jk) == 0) {
            Sys.sleep(1)
            file.jk <- dir(sampleDir.jk, pattern=suffix.i, full.names=T)
          }
          if (length(file.jk) == 0) {
            consensusDir.jk <- file.path(sampleDir.jk, paste0(sampleName.jk, '_consensus'))
            file.jk <- dir(consensusDir.jk, pattern=suffix.i, full.names=T)
            if (length(file.jk) == 0) {
              Sys.sleep(1)
              file.jk <- dir(consensusDir.jk, pattern=suffix.i, full.names=T)
            }
            if (length(file.jk) == 0) {
              realDir.jk <- Sys.readlink(consensusDir.jk)
              if (is.na(realDir.jk)) {
                Sys.sleep(1)
                realDir.jk <- Sys.readlink(consensusDir.jk)
                # if (is.na(realDir.jk)) stop('Sys.readlink(consensusDir.jk) is NA.')
                if (is.na(realDir.jk)) {
                  warning(paste('Empty file:', sampleName.jk, suffix.i))
                  realDir.jk <- ''
                }
              }
              if (basename(realDir.jk) == 'results') {
                projDir.jk <- dirname(dirname(realDir.jk))
                file.jk <- dir(projDir.jk, pattern=paste0(sample.j, '.*', suffix.i), full.names=T)
                if (length(file.jk) == 0) {
                  Sys.sleep(1)
                  file.jk <- dir(projDir.jk, pattern=paste0(sample.j, '.*', suffix.i), full.names=T)
                }
              } else {
                file.jk <- ''
              }
            } 
          }
          if(length(file.jk) > 1) {
            if (any(grepl('short', basename(file.jk)))) {
              file.jk <- file.jk[grepl('_short', basename(file.jk))]
            } else {
              warning(paste('More than one files were found, only the first one will be used for', sampleName.jk, suffix.i))
              file.jk <- file.jk[1]
            }
          }
          if (length(file.jk) == 0) file.jk <- ''
          if (!file.exists(file.jk)) {
            missingFileInfo.j <- paste(sampleName.jk, suffix.i, sep=': ')
            missingFileInfo <- c(missingFileInfo, missingFileInfo.j)
          } else {
            tmp <- readLines(file.jk[1], n=10)
            tmp <- tmp[!grepl('^#', tmp)]
            if (length(tmp) > 0) {
              # combinedInfo.j <- read.csv(file.j[1], head=TRUE, comment='#', as.is=TRUE, check.names = FALSE, blank.lines.skip = FALSE)
              if ((grepl('(variant)|(tmb)|(snp)', step.i, ignore.case = TRUE))) {
                combinedInfo.jk <- read.table(file.jk[1], sep=',', head=TRUE, comment='#', as.is=TRUE, colClasses=list(ref='character', alt='character'), check.names=FALSE, fill=TRUE)
              } else {
                combinedInfo.jk <- read.table(file.jk[1], sep=',', head=TRUE, comment='#', as.is=TRUE, check.names=FALSE, fill=TRUE, stringsAsFactors=FALSE)
              }
              if (step.i == 'NGSQC') {
                if (ncol(combinedInfo.jk) < nrow(combinedInfo.jk)) {
                  missingFileInfo.j <- paste(sampleName.jk, suffix.i, sep=': ')
                  missingFileInfo <- c(missingFileInfo, missingFileInfo.j)
                  next
                } else {
                  rownames(combinedInfo.jk) <- combinedInfo.jk[,1]
                  if (is.null(combinedInfo.jk$SampleID) && tolower(colnames(combinedInfo.jk)[1]) == 'sampleid') {
                    colnames(combinedInfo.jk)[1] <- 'SampleID'
                  }
                }
              }
              if (nrow(combinedInfo.jk) > 0) {
                combinedInfo.jk$sampleFolder <- sampleName.jk
                ## if simpleReturn, only output a subset of columns
                if (simpleReturn && step.i %in% c('Variant', 'Variant_short', 'TMB', 'SNP', 'otherVariant')) {
                  commCols.jk <- simpleVarCols[simpleVarCols %in% colnames(combinedInfo.jk)]
                  combinedInfo.jk <- combinedInfo.jk[, commCols.jk, drop=FALSE]
                }
                combinedInfo.jk$SampleID.org <- sample.j
              }
              if (is.null(combinedInfo) || (j == 1 && k == 1) || nrow(combinedInfo) == 0) {
                combinedInfo <- combinedInfo.jk
              } else if (nrow(combinedInfo.jk) > 0 & ncol(combinedInfo.jk) > 2) {
                combinedInfo <- combineGRanges(combinedInfo, combinedInfo.jk)
                # colNames.all <- unique(c(colnames(combinedInfo), colnames(combinedInfo.jk)))
                # colNames.all <- colNames.all[colNames.all != '']
                # missing.jk <- colNames.all[!(colNames.all %in% colnames(combinedInfo.jk))]
                # missing.comb.jk <- colNames.all[!(colNames.all %in% colnames(combinedInfo))]
                # if (length(missing.jk) > 0) {
                #     if (!is.data.frame(combinedInfo.jk))  combinedInfo.jk <- as.data.frame(combinedInfo.jk)
                #     combinedInfo.jk[, missing.jk] <- NA 
                # }
                # if (length(missing.comb.jk) > 0) {
                #     if (!is.data.frame(combinedInfo)) combinedInfo <- as.data.frame(combinedInfo)
                #     combinedInfo[,missing.comb.jk] <- NA
                # } 
                # combinedInfo <- rbind(combinedInfo[, colNames.all], combinedInfo.jk[, colNames.all]) 
              } else {
                na.jk <- c(sampleName.jk, rep(NA, ncol(combinedInfo) - 1))
                combinedInfo <- rbind(combinedInfo, na.jk)
              }
            }
          }
        }
      }
    }
    
    if (!is.null(sampleInfo) && nrow(combinedInfo) > 0) {
      combinedInfo$PatientID <- NA
      if (!is.null(combinedInfo$SampleID)) {
        sampleId <- sub('(_consensus)|(_merged)', '', combinedInfo$SampleID)
      } else if (!is.null(combinedInfo$sampleNames)) {
        sampleId <- sub('(_consensus)|(_merged)', '', combinedInfo$sampleNames)
      } else if (!is.null(combinedInfo$sampleFolder)) {
        sampleId <- sub('_[0-9]+$', '', basename(combinedInfo$sampleFolder))
        combinedInfo$SampleID <- sampleId
      } else if (!is.null(combinedInfo$sampleDir)) {
        sampleId <- sub('_[0-9]+$', '', basename(combinedInfo$sampleDir))
        combinedInfo$SampleID <- sampleId
      } else {
        stop('SampleID does not available in combinedInfo!')
      }
      ## remove suffix _.* to match sampleInfo rownames if needed
      naInd <- which(!(sampleId %in% rownames(sampleInfo)))
      if (length(naInd) > 0) {
        if (!is.null(combinedInfo$SampleID.org)) {
          sampleId[combinedInfo$SampleID.org %in% rownames(sampleInfo)] <- combinedInfo$SampleID.org[combinedInfo$SampleID.org %in% rownames(sampleInfo)]
          naInd <- which(!(sampleId %in% rownames(sampleInfo)))
        }
        if (length(naInd) > 0) {
          sampleId.new <- gsub('_.*', '', sampleId)
          if (any(sampleId.new %in% rownames(sampleInfo))) {
            naInd.update <- naInd[sampleId.new[naInd] %in% rownames(sampleInfo)]
            sampleId[naInd.update] <- gsub('_.*', '', sampleId[naInd.update])
          }
        }
      }
      
      selInd <- which(sampleId %in% rownames(sampleInfo))
      if (!is.null(combinedInfo$SpecimenType) && !is.null(sampleInfo$SpecimenType)) {
        combinedInfo[selInd, 'SpecimenType'] <- sampleInfo[sampleId[selInd], 'SpecimenType']
        excludeCols <- c(excludeCols, 'SpecimenType')
        # otherSampleCols <- otherSampleCols[!(otherSampleCols %in% 'SpecimenType')]
      }
      if (!is.null(combinedInfo$Specimen_type) && !is.null(sampleInfo$SpecimenType)) {
        combinedInfo[selInd, 'Specimen_type'] <- sampleInfo[sampleId[selInd], 'SpecimenType']
        excludeCols <- c(excludeCols, 'SpecimenType')
        # otherSampleCols <- otherSampleCols[!(otherSampleCols %in% 'SpecimenType')]
      }
      if (!is.null(combinedInfo$PatientID) && !is.null(sampleInfo$PatientID)) {
        combinedInfo[selInd, 'PatientID'] <- sampleInfo[sampleId[selInd], 'PatientID']
        # otherSampleCols <- otherSampleCols[!(otherSampleCols %in% 'PatientID')]
        excludeCols <- c(excludeCols, 'PatientID')
      }
      
      if (!is.null(otherSampleCols)) {
        otherSampleCols.i <- otherSampleCols[!(otherSampleCols %in% colnames(combinedInfo))]
        otherSampleCols.i <- otherSampleCols.i[otherSampleCols.i %in% colnames(sampleInfo)]
        otherSampleCols.i <- otherSampleCols.i[!sapply(otherSampleCols.i, function(ii) all(is.na(sampleInfo[, ii])))]
        otherSampleCols.i <- otherSampleCols.i[!(otherSampleCols.i %in% excludeCols)]
        if (length(otherSampleCols.i) > 0) {
          otherInfo <- matrix(NA, nrow=nrow(combinedInfo), ncol=length(otherSampleCols.i))
          colnames(otherInfo) <- otherSampleCols.i
          otherInfo[selInd, otherSampleCols.i] <- as.matrix(sampleInfo[sampleId[selInd], otherSampleCols.i])
          combinedInfo <- cbind(combinedInfo, otherInfo)
        }
      }
    }
    if (all(c("seqnames", "start", "end") %in% colnames(combinedInfo)) && nrow(combinedInfo) > 0) {
      combinedInfo$VariantID <- getVariantIDs(combinedInfo)
    }
    
    if (!is.null(combinedInfo$VariantID)) {
      if (!is.null(combinedInfo$SampleID)) {
        combinedInfo$ID <- paste(combinedInfo$SampleID, combinedInfo$VariantID, sep=':')
      } else if (!is.null(combinedInfo$sampleNames)) {
        combinedInfo$ID <- paste(combinedInfo$sampleNames, combinedInfo$VariantID, sep=':')
      }
      
      ## add variant count information
      varCountAcrossSample <- table(combinedInfo$VariantID[!duplicated(combinedInfo$ID)])
      combinedInfo$highFrequent.inbatch <- as.vector(varCountAcrossSample[combinedInfo$VariantID])
    }
    
    ## save combined results
    if (grepl('(variant)|(tmb)|(snp)', step.i, ignore.case = TRUE)) {
      if (is.null(combinedInfo$VariantType)) combinedInfo$VariantType <- ''
      combinedInfo$VariantType[is.na(combinedInfo$VariantType)] <- ''
      ## check whether missing any VariantType information
      if (any(combinedInfo$VariantType == '') && step.i == 'Variant_short') {
        combinedInfo <- combinedInfo[!is.na(combinedInfo$ref),,drop=FALSE]
        combinedInfo <- as.data.frame(checkVariantType(combinedInfo))
      }
      
      ## check concordant variants across multiple samples from the sample patient
      if (!is.null(combinedInfo$PatientID)) {
        PID <- paste(combinedInfo$PatientID, combinedInfo$VariantID, sep='_')
        dup.PID <- unique(PID[duplicated(PID)])
        combinedInfo$concordant <- PID %in% dup.PID
        ## check germline or not if there are germline samples together
        if (!is.null(combinedInfo$SpecimenType)) {
          sampleType <- combinedInfo$SpecimenType
          uniType <- unique(sampleType)
          if (length(uniType) > 1 && any(c('pbmc', 'germline', 'normal') %in% tolower(uniType))) {
            normInd <- tolower(sampleType) %in% c('pbmc', 'germline', 'normal')
            normPID <- PID[normInd & combinedInfo$altDepth >= 8]
            combinedInfo$germline <- PID %in% normPID
          }
        }
      }
      
    } else if (grepl('cnv', suffix.i, ignore.case = TRUE)) {
      # combinedInfo <- cnvPostFiltering(combinedInfo)
      ## check concordant cnv across multiple samples from the sample patient
      if (!is.null(combinedInfo$PatientID)) {
        PID <- paste(combinedInfo$PatientID, combinedInfo$Gene, combinedInfo$CNV_Type, sep='_')
        dup.PID <- unique(PID[duplicated(PID)])
        combinedInfo$concordant <- PID %in% dup.PID
      }
      if (!is.null(combinedInfo$ID)) {
        combinedInfo <- combinedInfo[order(combinedInfo$ID, combinedInfo$CNL_Type, decreasing=TRUE),,drop=FALSE]
      }
    }
    if (removeDuplicate && !is.null(combinedInfo$ID)) {
      combinedInfo <- combinedInfo[!duplicated(combinedInfo$ID),,drop=FALSE]
    }
    if (!is.null(combinedInfo$ID)) {
      combinedInfo <- combinedInfo[order(combinedInfo$ID, decreasing = FALSE),,drop=FALSE]
      ## check overlapping with variant
      if (tolower(step.i) %in% c('tmb', 'variant_short')) {
        ## the second element is "Variant"
        combinedInfo.variant <- combinedInfoList[[2]]
        selInd <- which(combinedInfo.variant$ID %in% combinedInfo$ID)
        if (length(selInd) > 0) {
          status <- paste(combinedInfo.variant$finalKeep[selInd], step.i, sep=',')
          status <- gsub('((TRUE)|(FALSE)),', '', status)
          combinedInfo.variant$finalKeep[selInd] <- status
          combinedInfoList[[2]] <- combinedInfo.variant
        }
      }
    }
    if (!is.null(combinedInfo$VariantID) && any(is.na(combinedInfo$VariantID))) {
      naInd <- which(is.na(combinedInfo$VariantID))
      if (length(naInd) > 0) {
        combinedInfo <- rbind(combinedInfo[-naInd,,drop=FALSE], combinedInfo[naInd,,drop=FALSE])
      }
    }
    combinedInfo$SampleID.org <-NULL
    combinedInfoList <- c(combinedInfoList, list(combinedInfo))
    
    if (!is.null(outputDir)) {
      if (file.exists(outputDir))
        write.csv(combinedInfo, file=file.path(outputDir, paste0(names(suffixList)[i], '_allCombinedInfo.csv')), row.names=FALSE)
    }
  }
  names(combinedInfoList) <- names(suffixList)
  
  if (length(missingFileInfo) > 0) {
    if (is.data.frame(missingFileInfo) && ncol(missingFileInfo) == 2) {
      colnames(missingFileInfo) <- c('SampleID', 'qc.suffix')
      rownames(missingFileInfo) <- NULL
    }
  }
  return(c(combinedInfoList, missingFileInfo=list(unique(missingFileInfo))))
}


##' getSampleSNPBC get SNP barcode by sample
##' @param sampleIDs a vector of Predicine sample IDs or a data.frame with SampleID column
##' @param snpBCInfo snp variantIDs or snp BC file
##' @param returnAllMatch whether return all matches of sampleIDs (mutiple runs may exist for the same sample)
##' @param sampleDirs outputBySample directory or individual directories for each sample
##' @param newAnalysis whether to do updated estimation of SNP BC
##' @param saveMode whether to save the result in "xxx_SNP_BC.csv" file
##' 
##' @return a list with following elements: snpBCs and missingSamples
##' 
getSampleSNPBC <- function(sampleIDs, snpBCInfo=NULL, sampleDirs=NULL, newAnalysis=FALSE, saveMode=TRUE) {
  
  if (is.null(sampleDirs)) sampleDirs <- options('outputBySampleDir')[[1]]
  if (is.null(sampleDirs)) stop('Please either provide sampleDirs or set it in options!')
  if (any(!dir.exists(sampleDirs))) stop(paste("cannot find directory: ", sampleDirs))
  
  if (length(sampleDirs) == 1) {
    sampleIDs <- unique(as.character(sampleIDs))
  }
  if (length(sampleIDs) > 1 & length(sampleDirs) == 1) {
    sampleDirs <- rep(sampleDirs, length(sampleIDs))
  }
  
  if (is.null(snpBCInfo)) {
    targetInfoDir <- options('targetInfoDir')
    snpBCInfo <- file.path(targetInfoDir, "SNP_sampleBC_5_190616.csv")
  }
  if (is.character(snpBCInfo) && length(snpBCInfo) == 1) {
    snpBCInfo <- read.csv(snpBCInfo, as.is=TRUE, colClasses='character')
    snpVariantIDs <- checkVariantID(snpBCInfo$VariantID)
  } else {
    if (grepl('.*[0-9]+:[0-9]+', snpVariantIDs[1])) {
      snpVariantIDs <- snpBCInfo
    } else {
      stop('snpBCInfo should be either the SNP VariantIDs or the SNPInfo file!')
    }
  }
  
  snpBCs <- snpVarInfo <- missingSamples <-  NULL
  for (j in 1:length(sampleIDs)) {
    sample.j <- sampleIDs[j]
    
    ## paste0(sampleName, "_variantInfo_SNP_BC.csv")
    newAnalysis.j <- newAnalysis
    if (!newAnalysis) {
      fileInfo.j <- findFileBySample(sample.j, file.suffix='_SNP_BC.csv', returnAllMatch=FALSE, sampleDirs=sampleDirs[j])
      snpBCFile.j <- fileInfo.j$filelist
      if (length(snpBCFile.j) > 0) {
        snpVarInfo.j <- read.csv(snpBCFile.j[1], as.is=TRUE, check.names=FALSE, colClasses='character')
        snpVarInfo.j <- asVRanges(snpVarInfo.j)
        newAnalysis.j <- FALSE 
      } else {
        newAnalysis.j <- TRUE 
      }
    } 
    if (newAnalysis.j) {
      fileInfo.j <- findFileBySample(sample.j, file.suffix='coverage.RData', returnAllMatch=FALSE, sampleDirs=sampleDirs[j])
      coverageFile.j <- fileInfo.j$filelist
      if (length(fileInfo.j$missingSamples) > 0) {
        missingSamples <- c(missingSamples, fileInfo.j$missingSamples)
        next
      }
      
      fileInfo.j <- findFileBySample(sample.j, file.suffix='rawVariants.RData', returnAllMatch=FALSE, sampleDirs=sampleDirs[j])
      rawVariantFile.j <- fileInfo.j$filelist
      if (length(fileInfo.j$missingSamples) > 0) {
        missingSample.i <- basename(fileInfo.j$missingSamples)
        missingSample.i <- sub("_[0-9]+$", '', missingSample.i)
        missingSamples <- c(missingSamples, missingSample.i)
        next
      }
      rawVariants <- get(load(rawVariantFile.j))
      snpVarInfo.j <- checkVariantCoverage(snpVariantIDs, rawVariants, coverageInfo=coverageFile.j, snpCall=TRUE)
      if (saveMode) {
        snpBCFile.j <- file.path(dirname(rawVariantFile.j), sub('rawVariants.RData', 'variantInfo_SNP_BC.csv', basename(rawVariantFile.j)))
        write.csv(as.data.frame(snpVarInfo.j), file=snpBCFile.j, row.names=FALSE)
      }
    }
    
    snpBC.j <- paste(snpVarInfo.j$snpCall, collapse='')
    snpBCs <- c(snpBCs, snpBC.j)
    snpVarInfo <- c(snpVarInfo, list(snpVarInfo.j))
  }
  names(snpBCs) <- names(snpVarInfo) <- sampleIDs[!(sampleIDs %in% missingSamples)]
  
  return(list(snpBCs=snpBCs, snpVarInfo=snpVarInfo, missingSamples=missingSamples))
}


##' combineProjectSummary combine the project level summary across different runs
##' @param projDirs a vector of project directories in the run folder
##' @param qc.suffix the character suffix of NGS QC file used for filename match
##' @param cnv.suffix the character suffix of CNV file used for filename match
##' @param fusion.suffix the character suffix of Fusion file used for filename match
##' @param variant.suffix the character suffix of Variant file used for filename match
##' @param tmb.suffix the character suffix of TMB Variant file used for filename match
##' @param snp.suffix the character suffix of SNP Variant file used for filename match
##' @param combinedRunPattern the runname pattern used to know whether the run is combined or not. A suffix "_C" will added in sample names for the combined samples
##' @param rmDuplicateBy remove duplicated variants (with same VariantID and SampleID) by totalDepth, VariantFreq or NA (not remove)
##' @param sampleInfo a data.frame keeps the sample information, which includes SampleID, PatientID and Specimen_type
##' @param otherSampleCols other columns in sampleInfo will be included in the combined output files
##' @param outputDir the location of output file, if it is NULL, not files will be saved
##' 
##' @return a list with following elements: QC, Variant, TMB, CNV, Fusion and missingFileInfo
##' 
combineProjectSummary <- function(projDirs, qc.suffix="_combined(_NGSQC_)?Summary_all.csv", cnv.suffix="_significant_CNV.*csv",
                                  fusion.suffix="_fusion.*short.csv", variant.suffix='_variants.*clinicOnly_short.csv',
                                  tmb.suffix='_TMB.csv', snp.suffix='_SNP_VariantInfo.csv', combinedRunPattern='__', rmDuplicateBy=c('totalDepth', 'VariantFreq', NA), 
                                  sampleInfo=NULL, otherSampleCols=NULL, outputDir=NULL) {
  
  rmDuplicateBy <- match.arg(rmDuplicateBy)
  projDirs <- unique(projDirs)
  if (any(!file.exists(projDirs))) {
    warning(paste(projDirs[!file.exists(projDirs)], collapse = ','))
    projDirs <- projDirs[file.exists(projDirs)]
    if (length(projDirs) == 0) stop('All project directories do not exist!')
  }
  specimenTypes <- c('germline', 'plasma', 'tissue', 'ffpe', 'urine')
  if (!is.null(sampleInfo)) {
    if (any(duplicated(sampleInfo[,'SampleID']))) {
      warning(paste('Duplicated SampleIDs:', paste(sampleInfo[duplicated(sampleInfo[,'SampleID']), 'SampleID'], collapse=',')))
      sampleInfo <- sampleInfo[!duplicated(sampleInfo[,'SampleID']),,drop=FALSE]
    }
    rownames(sampleInfo) <- sampleInfo[,'SampleID']
    if (is.null(sampleInfo$PatientID)) {
      if (!is.null(sampleInfo$External.Patient.ID)) {
        sampleInfo$PatientID <- sampleInfo$External.Patient.ID
      } else if (!is.null(sampleInfo$"External Patient ID")) {
        sampleInfo$PatientID <- sampleInfo$"External Patient ID"
      } else if (!is.null(sampleInfo$"External_Patient_ID")) {
        sampleInfo$PatientID <- sampleInfo$"External_Patient_ID"
      }
    }
  }
  
  sampleIDInfo <- NULL
  suffixList <- c(NGSQC=qc.suffix, Variant=variant.suffix, TMB=tmb.suffix, SNP=snp.suffix, CNV=cnv.suffix, Fusion=fusion.suffix)
  combinedInfoList <- missingProjInfo <- NULL
  for (i in 1:length(suffixList)) {
    step.i <- names(suffixList)[i]
    suffix.i <- suffixList[i]
    combinedInfo <- NULL
    for (j in 1:length(projDirs)) {
      projDir.j <- projDirs[j]
      runName.j <- basename(dirname(projDir.j))
      file.j <- dir(projDir.j, pattern=suffix.i, full.names=T)
      if (length(file.j) > 1) {
        warning(paste('More than one files were found, only the first one will be used for', projDir.j, suffix.i))
        file.j <- file.j[1]
      }
      if (length(file.j) == 0) {
        missingProjInfo.j <- c(projDir.j, suffix.i)
        missingProjInfo <- rbind(missingProjInfo, missingProjInfo.j)
      } else {
        combinedInfo.j <- read.table(file.j[1], sep=',', head=TRUE, comment='#', as.is=TRUE, check.names = FALSE, stringsAsFactors=FALSE)
        if (step.i == 'NGSQC') {
          if (any(duplicated(combinedInfo.j[,'SampleID']))) {
            warning(paste('Duplicated SampleIDs:', paste(combinedInfo.j[duplicated(combinedInfo.j[,'SampleID']), 'SampleID'], collapse=',')))
            combinedInfo.j <- combinedInfo.j[!duplicated(combinedInfo.j[,'SampleID']),,drop=FALSE]
          }
          rownames(combinedInfo.j) <- combinedInfo.j[,'SampleID']
          
          if (grepl(combinedRunPattern, runName.j)) {
            combinedInfo.j$sampleNames <- paste0(combinedInfo.j[,1], '_C')
          } else {
            combinedInfo.j$sampleNames <- combinedInfo.j[,1]
          }
          # combinedInfo.j <- as.matrix(combinedInfo.j)
          if (!is.null(sampleInfo)) {
            ## remove suffix _.* to match sampleInfo rownames if needed
            naInd.j <- which(!(combinedInfo.j[,'SampleID'] %in% rownames(sampleInfo)))
            if (length(naInd.j) > 0) {
              naId.new.j <- gsub('_.*', '', combinedInfo.j[naInd.j,'SampleID'])
              if (any(naId.new.j %in% rownames(sampleInfo))) {
                naInd.update.j <- naInd.j[naId.new.j %in% rownames(sampleInfo)]
                combinedInfo.j[naInd.update.j,'SampleID'] <- gsub('_.*', '', combinedInfo.j[naInd.update.j,'SampleID'])
              }
            }
            if (!is.null(sampleInfo$Specimen_type)) {
              specimenType <- tolower(sampleInfo[combinedInfo.j[,'SampleID'], 'Specimen_type'])
              if (any(is.na(specimenType))) {
                naInd <- is.na(specimenType)
                naSampleID <- paste(combinedInfo.j[naInd,'SampleID'], collapse=',')
                warning(paste(naSampleID, "Some samples don't have specimenType information!"))
                combinedInfo.j[!naInd, 'Specimen_type'] <- specimenType[!naInd]
              } else {
                combinedInfo.j[, 'Specimen_type'] <- specimenType
              }
            }
            combinedInfo.j$PatientID <- NA
            if (!is.null(sampleInfo$PatientID)) {
              patientID <- sampleInfo[combinedInfo.j[,'SampleID'], 'PatientID']
              if (any(is.na(patientID))) {
                naInd <- is.na(patientID)
                naSampleID <- paste(combinedInfo.j[naInd,'SampleID'], collapse=',')
                warning(paste(naSampleID, "Some samples don't have patient ID information!"))
                combinedInfo.j[!naInd, 'PatientID'] <- patientID[!naInd]
              } else {
                combinedInfo.j[, 'PatientID'] <- patientID
              }
            }
            if (!is.null(otherSampleCols)) {
              otherSampleCols <- otherSampleCols[!(otherSampleCols %in% colnames(combinedInfo.j))]
              otherSampleCols <- otherSampleCols[otherSampleCols %in% colnames(sampleInfo)]
              if (length(otherSampleCols) > 0) {
                otherInfo.j <- matrix(NA, nrow=nrow(combinedInfo.j), ncol=length(otherSampleCols))
                rownames(otherInfo.j) <- combinedInfo.j[, 'SampleID']
                colnames(otherInfo.j) <- otherSampleCols
                commIds.j <- intersect(combinedInfo.j[, 'SampleID'], sampleInfo$SampleID)
                if (length(commIds.j) > 0) {
                  otherInfo.j[commIds.j, otherSampleCols] <- as.matrix(sampleInfo[commIds.j, otherSampleCols])
                }
                combinedInfo.j <- cbind(combinedInfo.j, otherInfo.j)
              }
            }
          }
          specimen.j <- combinedInfo.j[, 'Specimen_type']
          specimen.j <- checkSpecimenType(gsub(';.*', '', specimen.j))
          combinedInfo.j[, 'Specimen_type'] <- specimen.j
          specimen.uni.j <- unique(specimen.j)
          if (any(is.na(specimen.uni.j)) && length(specimen.uni.j) == 2) {
            specimen.uni.j <- specimen.uni.j[!is.na(specimen.uni.j)]
            combinedInfo.j[is.na(specimen.j), 'Specimen_type'] <- specimen.uni.j
          }
          for (specimen.ij in specimen.uni.j) {
            combinedInfo.ij <- combinedInfo.j[which(combinedInfo.j[, 'Specimen_type'] == specimen.ij),,drop=FALSE]
            datatype.ij <- checkSpecimenType(specimen.ij)
            if (specimen.ij %in% names(combinedInfo)) {
              if (!is.null(combinedInfo[[datatype.ij]])) {
                colNames.all <- unique(c(colnames(combinedInfo[[datatype.ij]]), colnames(combinedInfo.ij)))
                missing.j <- colNames.all[!(colNames.all %in% colnames(combinedInfo.ij))]
                missing.comb.j <- colNames.all[!(colNames.all %in% colnames(combinedInfo[[datatype.ij]]))]
                if (length(missing.j) > 0) combinedInfo.ij[, missing.j] <- NA 
                if (length(missing.comb.j) > 0) combinedInfo[[specimen.ij]][,missing.comb.j] <- NA
                combinedInfo[[datatype.ij]] <- rbind(combinedInfo[[datatype.ij]][, colNames.all], combinedInfo.ij[, colNames.all]) 
              } else {
                combinedInfo[[datatype.ij]] <- combinedInfo.ij
              }
            } else {
              combinedInfo[[datatype.ij]] <- combinedInfo.ij
            }
          }
        } else {
          if (nrow(combinedInfo.j) > 0) {
            if (!is.null(combinedInfo.j$sampleNames)) {
              if (is.null(combinedInfo.j$SampleID)) combinedInfo.j$SampleID <- combinedInfo.j$sampleNames
              if (grepl(combinedRunPattern, runName.j)) {
                combinedInfo.j$sampleNames <- paste0(combinedInfo.j$sampleNames, '_C')
              } 
            } else if (!is.null(combinedInfo.j$SampleID)) {
              if (is.null(combinedInfo.j$sampleNames)) combinedInfo.j$sampleNames <- combinedInfo.j$SampleID
              if (grepl(combinedRunPattern, runName.j)) {
                combinedInfo.j$sampleNames <- paste0(combinedInfo.j$sampleNames, '_C')
              } 
            }
          } else {
            if (is.null(combinedInfo.j$sampleNames)) combinedInfo.j$sampleNames <- logical(0)
            if (is.null(combinedInfo.j$SampleID)) combinedInfo.j$SampleID <- logical(0)
          }
          if (is.null(combinedInfo) || (j == 1)) {
            combinedInfo <- combinedInfo.j
          } else {
            colNames.all <- unique(c(colnames(combinedInfo), colnames(combinedInfo.j)))
            missing.j <- colNames.all[!(colNames.all %in% colnames(combinedInfo.j))]
            missing.comb.j <- colNames.all[!(colNames.all %in% colnames(combinedInfo))]
            if (length(missing.j) > 0) combinedInfo.j[, missing.j] <- character(nrow(combinedInfo.j)) 
            if (length(missing.comb.j) > 0) combinedInfo[,missing.comb.j] <- character(nrow(combinedInfo))
            combinedInfo <- rbind(combinedInfo[, colNames.all], combinedInfo.j[, colNames.all]) 
          }
        }
      }
    }
    
    ## save combined results
    if (grepl('(variant)|(tmb)|(snp)', suffix.i, ignore.case = TRUE)) {
      combInd <- grep('_C$', combinedInfo$sampleNames)
      if (length(combInd) > 0) {
        combinedInfo <- rbind(combinedInfo[combInd,,drop=FALSE], combinedInfo[-combInd,,drop=FALSE])
      }
      ## remove duplicated variants (having same VariantID and SampleID)
      if (!is.na(rmDuplicateBy)) {
        if (rmDuplicateBy == 'totalDepth') {
          combinedInfo <- combinedInfo[order(combinedInfo$totalDepth, decreasing = TRUE),]
        } else {
          combinedInfo <- combinedInfo[order(combinedInfo$VariantFreq, decreasing = TRUE),]
        }
        if (!is.null(combinedInfo$SampleID)) {
          combinedInfo$ID <- paste(combinedInfo$SampleID, combinedInfo$VariantID, sep='_')
        } else if (!is.null(combinedInfo$sampleNames)) {
          combinedInfo$ID <- paste(combinedInfo$sampleNames, combinedInfo$VariantID, sep='_')
        } 
        if (!is.null(combinedInfo$ID))
          combinedInfo <- combinedInfo[!duplicated(combinedInfo$ID),]
      }
      ## check concordant variants across multiple samples from the sample patient
      if (!is.null(combinedInfo$PatientID)) {
        PID <- paste(combinedInfo$PatientID, combinedInfo$VariantID, sep='_')
        dup.PID <- unique(PID[duplicated(PID)])
        combinedInfo$concordant <- PID %in% dup.PID
      }
      
      ## add variant count information
      varCountAcrossSample <- table(combinedInfo$VariantID)
      combinedInfo$highFrequent.inbatch <- as.vector(varCountAcrossSample[combinedInfo$VariantID])
    } else if (grepl('cnv', suffix.i, ignore.case = TRUE)) {
      ## check concordant cnv across multiple samples from the sample patient
      if (!is.null(combinedInfo$PatientID)) {
        PID <- paste(combinedInfo$PatientID, combinedInfo$Gene, combinedInfo$CNV_Type, sep='_')
        dup.PID <- unique(PID[duplicated(PID)])
        combinedInfo$concordant <- PID %in% dup.PID
      }
    }
    
    if (step.i == 'NGSQC') {
      selCols <- c('SampleID', 'Specimen_type', 'PatientID', otherSampleCols)
      selCols <- selCols[selCols %in% colnames(combinedInfo[[1]])]
      sampleIDInfo <- do.call(rbind, lapply(combinedInfo, function(x) x[, selCols, drop=FALSE]))
      sampleIDInfo <- sampleIDInfo[!is.na(sampleIDInfo$SampleID),,drop=FALSE]
      sampleIDInfo <- sampleIDInfo[!duplicated(sampleIDInfo$SampleID),,drop=FALSE]
      rownames(sampleIDInfo) <- sampleIDInfo$SampleID
    } else if (!is.null(sampleIDInfo)) {
      combinedInfo$SpecimenType <- combinedInfo$PatientID <- NA
      sampleId <- sub('(_consensus)|(_merged)', '', combinedInfo$SampleID)
      ## remove suffix _.* to match sampleInfo rownames if needed
      naInd <- which(!(sampleId %in% rownames(sampleInfo)))
      if (length(naInd) > 0) {
        sampleId.new <- gsub('_.*', '', sampleId)
        if (any(sampleId.new %in% rownames(sampleInfo))) {
          naInd.update <- naInd[sampleId.new[naInd] %in% rownames(sampleInfo)]
          sampleId[naInd.update] <- gsub('_.*', '', sampleId[naInd.update])
        }
      }
      
      selInd <- which(sampleId %in% rownames(sampleIDInfo))
      combinedInfo[selInd, 'SpecimenType'] <- sampleIDInfo[sampleId[selInd], 'Specimen_type']
      if (!is.null(sampleIDInfo$PatientID))
        combinedInfo[selInd, 'PatientID'] <- sampleIDInfo[sampleId[selInd], 'PatientID']
      
      if (!is.null(otherSampleCols)) {
        otherSampleCols <- otherSampleCols[!(otherSampleCols %in% colnames(combinedInfo))]
        otherSampleCols <- otherSampleCols[otherSampleCols %in% colnames(sampleIDInfo)]
        if (length(otherSampleCols) > 0) {
          otherInfo <- matrix(NA, nrow=nrow(combinedInfo), ncol=length(otherSampleCols))
          colnames(otherInfo) <- otherSampleCols
          otherInfo[selInd, otherSampleCols] <- as.matrix(sampleIDInfo[sampleId[selInd], otherSampleCols])
          combinedInfo <- cbind(combinedInfo, otherInfo)
        }
      }
    }
    combinedInfoList <- c(combinedInfoList, list(combinedInfo))
    
    if (!is.null(outputDir)) {
      if (file.exists(outputDir)) {
        date.i <- Sys.Date()
        if (step.i == 'NGSQC') {
          for (k in 1:length(combinedInfo)) {
            combinedInfo.k <- combinedInfo[[k]]
            datatype.k <- names(combinedInfo)[k]
            if (nrow(combinedInfo.i) > 1) {
              write.csv(combinedInfo.k, file=file.path(outputDir, paste0(datatype.k, '_', names(suffixList)[i], '_allCombinedInfo_', date.i, '.csv')), row.names=FALSE)
            } else {
              next
            }
          }
        } else {
          write.csv(combinedInfo, file=file.path(outputDir, paste0(names(suffixList)[i], '_allCombinedInfo_', date.i, '.csv')), row.names=FALSE)
        }
      }
    }        
  }
  names(combinedInfoList) <- names(suffixList)
  
  if (length(missingFileInfo) > 0) {
    if (is.data.frame(missingFileInfo) && ncol(missingFileInfo) == 2) {
      colnames(missingFileInfo) <- c('SampleID', 'qc.suffix')
      rownames(missingFileInfo) <- NULL
    }
  }
  
  return(c(combinedInfoList, missingFileInfo=list(missingProjInfo)))
}



##' mappingBamfile.externalID symlink internal consensus bamfile to destination directory using external sample ID
##' @param ngsQCInfo NGS QC information output by NGS pipeline, or a vector of PredicineID
##' @param destDir destination directory to save the symlink files
##' @param internalDir internal project directories or column name in ngsQCInfo, which keeps the internal consensus bam files. Skipped when ngsQCInfo is a vector of sampleIDs
##' @param externalID external IDs or column name in ngsQCInfo, which keeps the external sampleIDs
##' @param includeSampleID whether include SampleID in the new filenames
##' @param otherCols other columns to be included in the mapping info if existing in ngsQCInfo
##' @param saveSuffix saveSuffix used in destination bam files.
##' 
##' @return Invisibly return mapping info.
##' 
mappingBamfile.externalID <- function(ngsQCInfo, destDir, internalDir='ProjectDir', externalID='ExternalID', includeSampleID=FALSE, otherCols=c('Panel', 'specimenType', 'PatientID', 'trialVisitNum'), saveSuffix='') {
  
  if (is.data.frame(ngsQCInfo)) {
    colnames(ngsQCInfo) <- tolower(colnames(ngsQCInfo))
    if (!file.exists(destDir)) stop('destDir does not exists!')
    
    if (is.null(ngsQCInfo$sampleid)) {
      if (!is.null(ngsQCInfo$samplenames)) {
        ngsQCInfo$sampleid <- ngsQCInfo$samplenames
      } else {
        stop('SampleID column is missing in the ngsQCInfo!')
      }
    }
    sampleID <- ngsQCInfo$sampleid
    
    ngsQCInfo.uni <- ngsQCInfo[!duplicated(as.character(ngsQCInfo$sampleid)),]
    if (is.character(internalDir)) {
      internalDir <- tolower(internalDir)
      if (internalDir %in% colnames(ngsQCInfo.uni)) {
        internalDir <- ngsQCInfo.uni[, internalDir]
      } else {
        stop('internalDir column is missing in the ngsQCInfo!')
      }
    }
    if (is.character(externalID)) {
      externalID <- tolower(externalID)
      if (externalID %in% colnames(ngsQCInfo.uni)) {
        externalID <- ngsQCInfo.uni[, externalID]
      } else {
        warning('externalID column is missing in the ngsQCInfo. Internal ID will be used!')
        externalID <- sampleID
      }
    }
    names(externalID) <- sampleID
    bamfiles.old <- file.path(internalDir, paste0(sampleID, '_consensus.bam'))
    if (any(!file.exists(bamfiles.old))) {
      miss.ind <- which(!file.exists(bamfiles.old))
      tmp <- file.path(internalDir[miss.ind], paste0(sampleID[miss.ind], '_merged_cons_sort.bam'))
      bamfiles.old[miss.ind[file.exists(tmp)]] <- tmp[file.exists(tmp)]
      miss.ind <- miss.ind[!file.exists(tmp)]
      if (length(miss.ind) > 0) {
        tmp <- file.path(internalDir[miss.ind], paste0(sampleID[miss.ind], '_cons_sort.bam'))
        bamfiles.old[miss.ind[file.exists(tmp)]] <- tmp[file.exists(tmp)]
      }
    }
    baifiles.old <- file.path(internalDir, paste0(sampleID, '_consensus.bam.bai'))
    if (any(!file.exists(baifiles.old))) {
      miss.ind <- which(!file.exists(baifiles.old))
      tmp <- file.path(internalDir[miss.ind], paste0(sampleID[miss.ind], '_merged_cons_sort.bam.bai'))
      baifiles.old[miss.ind[file.exists(tmp)]] <- tmp[file.exists(tmp)]
      miss.ind <- miss.ind[!file.exists(tmp)]
      if (length(miss.ind) > 0) {
        tmp <- file.path(internalDir[miss.ind], paste0(sampleID[miss.ind], '_cons_sort.bam.bai'))
        baifiles.old[miss.ind[file.exists(tmp)]] <- tmp[file.exists(tmp)]
      }
    }        
    names(bamfiles.old) <- sampleID
    names(baifiles.old) <- sampleID
  } else if (is.vector(ngsQCInfo)) {
    sampleID <- as.character(ngsQCInfo)
    bamfileInfo <- findFileBySample(sampleID, file.suffix='consensus.bam')
    bamfiles.old <- bamfileInfo$filelist
    missingSamples <- bamfileInfo$missingSamples
    if (length(missingSamples) > 0) {
      warning(paste(paste(bamfileInfo$missingSamples, collapse=','), 'sample folder does not exists!'))
    }
    if (length(externalID) != length(sampleID)) stop('length of externalID and sampleID should be the same!')
    names(externalID) <- sampleID
    externalID <- externalID[sampleID %in% names(bamfiles.old)]
  }
  
  externalID <- gsub(' ', '_', externalID)
  if (any(duplicated(externalID)) && !includeSampleID) {
    includeSampleID <- TRUE
    cat('SampleID will be included in the filename to make them unique!')
  }
  bamfiles.old.all <- bamfiles.new.all <- externalID.all <- NULL
  for (i in 1:length(bamfiles.old)) {
    sample.i <- names(bamfiles.old)[i]
    bamfile.old.i <- bamfiles.old[i]
    externalID.i <- externalID[sample.i]
    externalID.all <- c(externalID.all, rep(externalID.i, 2))
    if (!file.exists(bamfile.old.i)) {
      warning(paste(paste0(sample.i, '_consensus.bam'), 'file not exists!'))
      bamfiles.new <- c(bamfiles.new, NA)
      next
    }
    if (includeSampleID && !grepl(sample.i, externalID.i)) {
      externalID.i <- paste(externalID.i, sample.i, sep='_')
    }
    if (nchar(saveSuffix) > 0) {
      bamfile.new.i <- file.path(destDir, paste0(externalID.i, '_', saveSuffix, '_consensus.bam'))
    } else {
      bamfile.new.i <- file.path(destDir, paste0(externalID.i, '_consensus.bam'))
    }
    system(paste('ln -s', bamfile.old.i, bamfile.new.i))
    baifile.old.i <- baifiles.old[i]
    baifile.new.i <- paste0(bamfile.new.i, '.bai')
    system(paste('ln -s', baifile.old.i, baifile.new.i))
    bamfiles.new.all <- c(bamfiles.new.all, bamfile.new.i, baifile.new.i)
    bamfiles.old.all <- c(bamfiles.old.all, bamfile.old.i, baifile.old.i)
  }
  
  mappingFiles <- data.frame(bamfile_old=bamfiles.old.all, bamfile_new=basename(bamfiles.new.all), SampleID=rep(names(bamfiles.old), each=2), 
                             ExternalID=externalID.all)
  ## additional columns to be included in the output
  otherCols <- otherCols[tolower(otherCols) %in% tolower(colnames(ngsQCInfo))]
  if (!is.null(otherCols) && is.data.frame(ngsQCInfo)) {
    colnames(ngsQCInfo) <- tolower(colnames(ngsQCInfo)) 
    rownames(ngsQCInfo) <- ngsQCInfo$sampleid
    selInd <- colnames(ngsQCInfo) %in% tolower(otherCols)
    otherInfo <- ngsQCInfo[mappingFiles$SampleID, selInd, drop=FALSE]
    mappingFiles <- cbind(mappingFiles, otherInfo)
  }
  return(invisible(mappingFiles))
}


##' mappingFastqfile.externalID symlink internal refined fastq file to destination directory using external sample ID
##' @param ngsQCInfo NGS QC information output by NGS pipeline, or a vector of PredicineID
##' @param destDir destination directory to save the symlink files
##' @param internalDir internal project directories or column name in ngsQCInfo, which keeps the internal consensus bam files. Skipped when ngsQCInfo is a vector of sampleIDs
##' @param externalID external IDs or column name in ngsQCInfo, which keeps the external sampleIDs
##' @param includeSampleID whether include SampleID in the new filenames
##' @param otherCols other columns to be included in the mapping info if existing in ngsQCInfo
##' @param saveSuffix saveSuffix used in destination bam files.
##' 
##' @return Invisibly return mapping info.
##' 
mappingFastqfile.externalID <- function(ngsQCInfo, destDir, internalDir='ProjectDir', externalID='ExternalID', includeSampleID=FALSE, otherCols=c('Panel', 'specimenType', 'PatientID', 'trialVisitNum'), saveSuffix='') {
  
  if (is.data.frame(ngsQCInfo)) {
    colnames(ngsQCInfo) <- tolower(colnames(ngsQCInfo))
    if (!file.exists(destDir)) stop('destDir does not exists!')
    
    if (is.null(ngsQCInfo$sampleid)) {
      if (!is.null(ngsQCInfo$samplenames)) {
        ngsQCInfo$sampleid <- ngsQCInfo$samplenames
      } else {
        stop('SampleID column is missing in the ngsQCInfo!')
      }
    }
    sampleID <- ngsQCInfo$sampleid
    
    ngsQCInfo.uni <- ngsQCInfo[!duplicated(as.character(ngsQCInfo$sampleid)),]
    if (is.character(internalDir)) {
      internalDir <- tolower(internalDir)
      if (internalDir %in% colnames(ngsQCInfo.uni)) {
        internalDir <- ngsQCInfo.uni[, internalDir]
      } else {
        stop('internalDir column is missing in the ngsQCInfo!')
      }
    }
    if (!grepl('fastq', internalDir) && grepl('lbwfresult', internalDir)) {
      internalDir <- gsub('lbwfresult.*$', '', internalDir)
      internalDir <- sapply(internalDir, function(x) {
        tt.x <- dir(x, pattern='^refinedfastq.*', full.names=TRUE)
        if (length(tt.x) == 0) tt.x <- NA
        return(tt.x)
      })
      if (length(internalDir) == 0) {
        stop("internalDir can't be found!")
      } 
    }
    if (is.character(externalID[1])) {
      externalID <- tolower(externalID)
      if (externalID %in% colnames(ngsQCInfo.uni)) {
        externalID <- ngsQCInfo.uni[, externalID]
      } else {
        warning('externalID column is missing in the ngsQCInfo. Internal ID will be used!')
        externalID <- sampleID
      }
    }
    names(externalID) <- sampleID
    fastq.r1.old <- file.path(internalDir, paste0(sampleID, '_R1_clean.fastq.gz'))
    fastq.r2.old <- file.path(internalDir, paste0(sampleID, '_R2_clean.fastq.gz'))
    names(fastq.r1.old) <- sampleID
    names(fastq.r2.old) <- sampleID
  } else {
    stop('ngsQCInfo should be a data.frame with NGS QC information!')
  }
  
  externalID <- gsub(' ', '_', externalID)
  if (any(duplicated(externalID)) && !includeSampleID) {
    includeSampleID <- TRUE
    cat('SampleID will be included in the filename to make them unique!')
  }
  fastq.old.all <- fastq.new.all <- externalID.all <- NULL
  for (i in 1:length(fastq.r1.old)) {
    sample.i <- names(fastq.r1.old)[i]
    fastq.r1.old.i <- fastq.r1.old[i]
    externalID.i <- externalID[sample.i]
    externalID.all <- c(externalID.all, rep(externalID.i, 2))
    if (!file.exists(fastq.r1.old.i)) {
      fastq.r1.old.i <- dir(internalDir, pattern=paste0(sampleID, '.*_R1_clean.fastq(.gz)?$'))
      if (length(fastq.r1.old.i) > 1) fastq.r1.old.i <- fastq.r1.old.i[1]
      if (length(fastq.r1.old.i) == 0) {
        warning(paste(paste0(sample.i, '_R1_clean.fastq.gz'), 'file not exists!'))
        fastq.new.all <- c(fastq.new.all, NA)
        next
      }
    }
    if (includeSampleID && !grepl(sample.i, externalID.i)) {
      externalID.i <- paste(externalID.i, sample.i, sep='_')
    }
    if (nchar(saveSuffix) > 0) {
      fastq.r1.new.i <- file.path(destDir, paste0(externalID.i, '_', saveSuffix, '_R1.fastq.gz'))
    } else {
      fastq.r1.new.i <- file.path(destDir, paste0(externalID.i, '_R1.fastq.gz'))
    }
    system(paste('ln -s', fastq.r1.old.i, fastq.r1.new.i))
    fastq.r2.old.i <- sub('_R1_clean', '_R2_clean', fastq.r1.old.i)
    fastq.r2.new.i <- sub('_R1', '_R2', fastq.r1.new.i)
    system(paste('ln -s', fastq.r2.old.i, fastq.r2.new.i))
    fastq.new.all <- c(fastq.new.all, fastq.r1.new.i, fastq.r2.new.i)
    fastq.old.all <- c(fastq.old.all, fastq.r1.old.i, fastq.r2.old.i)
  }
  
  mappingFiles <- data.frame(fastq_old=fastq.old.all, fastq_new=basename(fastq.old.all), SampleID=rep(names(fastq.r1.old), each=2), 
                             ExternalID=externalID.all)
  ## additional columns to be included in the output
  otherCols <- otherCols[tolower(otherCols) %in% tolower(colnames(ngsQCInfo))]
  if (!is.null(otherCols) && is.data.frame(ngsQCInfo)) {
    colnames(ngsQCInfo) <- tolower(colnames(ngsQCInfo)) 
    rownames(ngsQCInfo) <- ngsQCInfo$sampleid
    selInd <- colnames(ngsQCInfo) %in% tolower(otherCols)
    otherInfo <- ngsQCInfo[mappingFiles$SampleID, selInd, drop=FALSE]
    mappingFiles <- cbind(mappingFiles, otherInfo)
  }
  return(invisible(mappingFiles))
}


##' findSampleDirs mapping pipeline output files to OutputBySample by Sample
##' @param sampleIDs a vector of Predicine IDs 
##' @param outputBySampleDir the destination OutputBySample folder
##' @param returnAllMatch whehter return all matches or not
##' @param pversion pipeline version. If provided, results using matched pipeline version will be returned first
##' @param verbose whether print out detailed progress information
##' 
##' @return a list of sampleDirs and missingDirs
##' @examples 
##'   findSampleDirs(sampleIDs='P003294')
findSampleDirs <- function(sampleIDs, outputBySampleDir=NULL, includeWGS=FALSE, returnAllMatch=FALSE, pversion=NULL, verbose=FALSE) {
  
  if (is.null(outputBySampleDir)) outputBySampleDir <- options('outputBySampleDir')[[1]]
  if (is.null(outputBySampleDir)) stop('Please either provide outputBySampleDir or set it in options!')
  if (!dir.exists(outputBySampleDir[1])) stop(paste("cannot find directory: ", outputBySampleDir[1]))
  if (length(outputBySampleDir) == 1 && length(sampleIDs) > 1) {
    sampleDirs <- rep(outputBySampleDir, length(sampleIDs))
  } else {
    sampleDirs <- outputBySampleDir
  }
  
  if (is.null(pversion)) {
    pversion <- ''
  } else {
    pversion <- unique(c(pversion, ''))
  }
  
  sampleDirs.all <- missingDirs <- NULL
  for (j in 1:length(sampleIDs)) {
    sample.j <- sampleIDs[j]
    sample.p.j <- paste0('^', sample.j)
    if (verbose) cat(paste('Processing', sample.j, '\n'))
    if (!grepl(sample.j, basename(sampleDirs[j]))) {
      sampleDir.all.j <- NULL
      for (pv.k in pversion) {
        sampleDir.check.k <- file.path(sampleDirs[j], pv.k)
        sampleDir.k <- dir(sampleDir.check.k, pattern=paste0(sample.p.j, '(_[GR])?_[_0-9]+$'), full.names = TRUE)
        if (length(sampleDir.k) == 0) {
          if (grepl('A[0-9]$', sample.p.j)) {
            sampleDir.k <- dir(sampleDir.check.k, pattern=paste0(sub('A[0-9]+$', '', sample.p.j), '(_[GR])?_[_0-9]+$'), full.names = TRUE)
          } else {
            sampleDir.k <- dir(sampleDir.check.k, pattern=paste0(sample.p.j, '_A[0-9]+(_[GR])?_[_0-9]+$'), full.names = TRUE)
          }
        }
        if (length(sampleDir.k) == 0) {
          sampleDir.k <- dir(file.path(sampleDir.check.k, 'Other'), pattern=paste0(sample.p.j, '(_[GR])?_[_0-9]+$'), full.names = TRUE)
        }
        if (length(sampleDir.k) == 0) {
          sampleDir.k <- dir(file.path(sampleDir.check.k, 'Controls'), pattern=paste0(sample.p.j, '(_[GR])?_[_0-9]+$'), full.names = TRUE)
        }
        sampleDir.all.j <- c(sampleDir.all.j, sampleDir.k)
      }
      if (length(sampleDir.all.j) > 0) {
        if (any(grepl(sample.j, basename(sampleDir.all.j)))) {
          sampleDir.all.j <- sampleDir.all.j[grepl(sample.j, basename(sampleDir.all.j))]
        }
        sampleDir.j <- sampleDir.all.j # [1]
      } else {
        sampleDir.j <- NULL
      }
      if (length(sampleDir.j) == 0) {
        missingDirs <- c(missingDirs, sample.j)
        next
      } else {
        if (length(sampleDir.j) > 1 && grepl(sample.j, sampleDir.j)) {
          sampleDir.j <- sampleDir.j[grepl(sample.j, basename(sampleDir.j))]
          if (any(grepl('wgs', basename(sampleDir.j), ignore.case=TRUE)) && includeWGS) {
            sampleDir.wgs.j <- sampleDir.j[grepl("wgs", basename(sampleDir.j), ignore.case = T)]
            pversion.check <- gsub("\\.", "\\\\.", pversion[pversion != ""])
            if (length(pversion.check) > 1) pversion.check <- paste0('(', paste0(pversion.check, collapse=')|('), ')')
            idx <- which(grepl(pversion.check, sampleDir.wgs.j))
            if (length(idx) > 0){
              sampleDir.wgs.j <- sampleDir.wgs.j[idx]
            }
            sampleDir.j <- c(sampleDir.j[!grepl('wgs', basename(sampleDir.j), ignore.case=TRUE)], sampleDir.wgs.j)
          } else {
            sampleDir.j <- sampleDir.j[!grepl('wgs', basename(sampleDir.j), ignore.case=TRUE)]
          }
        }
        if (length(sampleDir.j) > 1 && !returnAllMatch) {
          if (any(pversion != '')) {
            pversion.j <- paste0('(', paste0(pversion[pversion != ''], collapse=')|('), ')')
            sampleDir.j <- sampleDir.j[grepl(pversion.j, sampleDir.j)]
          }
          ## handle combined samples "P003294_190613_190628"
          combinedInd.j <- grepl('_[0-9]{4,}_[0-9]{4,}$', basename(sampleDir.j))
          if (length(sampleDir.j) > 1) {
            dates.j <- gsub('^.*_', '', sampleDir.j)
            sampleDir.j <- sampleDir.j[order(combinedInd.j, dates.j, -nchar(sampleDir.j), decreasing = TRUE)]
          }
          sampleDir.j <- sampleDir.j[1]
        }  
      } 
    } else {
      sampleDir.j <- sampleDirs[j]
    }
    sampleDirs.all <- c(sampleDirs.all, list(sampleDir.j))
  }
  if (length(sampleDirs.all) > 0) {
    names(sampleDirs.all) <- sampleIDs[!(sampleIDs %in% missingDirs)]
  }
  
  return(list(sampleDirs=sampleDirs.all, missingDirs=missingDirs))
}


##' checkMissingVariant check missing variants across samples from the same patients
##' @param varInfo.org variant information in VRanges or related data.frame format
##' @param qcInfo the NGS QC information, which is used to defined the full sample list for missing variant checking
##' @param trueVariantColumn if the trueVariantColumn column is defined in variantInfo.org, the related variants will be used for missing variant check
##' @param pversion pipeline version used to search results
##' @param patients a subset of patient IDs for processing 
##' @param AF.th only check variants with variant frequency higher than AF.th
##' @param file.suffix filename suffix to search related files for missing variant checking
##' @param missVarOnly whether only return missing variant information
##' @param updateVariant whether also update the variants in addition to checking missing variants
##' @param sampleDirs sampleDir locations to search sample level pipeline output results
##' 
##' @return missing variant information or combined results with missing variang information
##' 
checkMissingVariant <- function(varInfo.org, qcInfo=NULL, trueVariantColumn='TrueVariant', pversion=NULL, patients=NULL, analysisType='cfdna_analysis', AF.th=NULL, 
                                file.suffix='general.csv', missVarOnly=FALSE, updateVariant=FALSE, sampleDirs=NULL) {
  
  if (is.null(varInfo.org$PatientID)) {
    stop('PatientID field is required!')
  }
  varInfo.org <- as.data.frame(varInfo.org)
  varInfo.org$ID <- paste(varInfo.org$sampleNames, varInfo.org$VariantID, sep='_')
  varInfo <- asVRanges(varInfo.org)
  varInfo.org$Patched <- FALSE
  varInfo$PatientID <- as.character(varInfo$PatientID)
  if (!is.null(qcInfo)) {
    qcInfo <- qcInfo[!duplicated(qcInfo$SampleID),]
    patient2sample <- split(qcInfo$SampleID, qcInfo$PatientID)
  } else {
    patient2sample <- lapply(split(as.character(sampleNames(varInfo)), varInfo$PatientID), unique)
  }
  patient.count <- sapply(patient2sample, length)
  
  if (is.null(patients)) {
    uniPatient <- names(patient.count)
  } else {
    uniPatient <- patients
  }
  thresholds <- getVariantThresholds(analysisType)
  AF.th.lod <- thresholds$postTh$AF.th
  AF.th.wl <- thresholds$postTh$minAF.th * 2
  if (is.null(AF.th)) AF.th <- AF.th.lod * 2
  varInfo$ID <- paste(sampleNames(varInfo), varInfo$VariantID, sep='_')
  missVarInfo <- varInfo.updated <- NULL
  for (patient.i in uniPatient) {
    cat(paste('Processing patient', patient.i, ':\n'))
    varInfo.i <- varInfo[varInfo$PatientID == patient.i]
    if (patient.count[patient.i] > 1) {
      sampleIDs.i <- patient2sample[[patient.i]]
      if (length(sampleIDs.i) > 0) {
        if (!is.null(trueVariantColumn) && trueVariantColumn %in% colnames(varInfo.org)) {
          checkVarIDs.i <- varInfo.i$VariantID[toupper(values(varInfo.i)[, trueVariantColumn]) %in% c('YES', 'TRUE', 'T')]
        } else {
          checkVarIDs.i <- unique(varInfo.i$VariantID[varInfo.i$VariantFreq > AF.th & as.character(sampleNames(varInfo.i)) %in% sampleIDs.i])
        }
        missVarInfo.i <- getVariantBySample(sampleIDs.i, checkVarIDs.i, file.suffix = file.suffix, pversion=pversion, sampleDirs=sampleDirs, returnAllMatch=TRUE)$varInfo
      } else {
        if (!is.null(trueVariantColumn) && trueVariantColumn %in% colnames(varInfo.org)) {
          checkVarIDs.i <- varInfo.i$VariantID[toupper(values(varInfo.i)[, trueVariantColumn]) %in% c('YES', 'TRUE', 'T')]
        } else {
          checkVarIDs.i <- unique(varInfo.i$VariantID[varInfo.i$VariantFreq > AF.th])
        }
        missVarInfo.i <- NULL
      }
      if (length(missVarInfo.i) > 0) {
        sampleNames(missVarInfo.i) <- sub('_consensus', '', sampleNames(missVarInfo.i))
        missVarInfo.i$VariantID <- getVariantIDs(missVarInfo.i)
        missVarInfo.i$PatientID <- patient.i
        missVarInfo.i$ID <- paste(sampleNames(missVarInfo.i), missVarInfo.i$VariantID, sep='_')
        missVarInfo.i <- missVarInfo.i[order(totalDepth(missVarInfo.i), decreasing = TRUE)]
        missVarInfo.i <- missVarInfo.i[!duplicated(missVarInfo.i$ID)]
        varInfo.updated.i <- missVarInfo.i[(missVarInfo.i$ID %in% varInfo.i$ID)]
        missVarInfo.i <- missVarInfo.i[!(missVarInfo.i$ID %in% varInfo.i$ID)]
        ## fill in the missing annotation information if it's using rawVariants
        if (grepl('rawVariants', file.suffix) && length(missVarInfo.i) > 0) {
          checkCols <- c('white.list', 'wl.tier', 'GC.Percent', 'tandemRepeat')
          varInfo.uni.i <- varInfo.i[!duplicated(varInfo.i$VariantID)]
          names(varInfo.uni.i) <- varInfo.uni.i$VariantID
          missCols <- colnames(values(varInfo.i))[!(colnames(values(varInfo.i)) %in% colnames(values(missVarInfo.i)))]
          tmp.ind <- max(which(missCols %in% checkCols))
          missCols <- missCols[1:tmp.ind]
          otherCols <- c("VariantType")
          missCols <- unique(c(missCols, otherCols[otherCols %in% colnames(values(varInfo.i))]))
          missAnnotation <- as.data.frame(values(varInfo.uni.i[missVarInfo.i$VariantID])[,missCols])
          values(missVarInfo.i) <- cbind(values(missVarInfo.i), missAnnotation)
        }
        missVarInfo <- c(missVarInfo, missVarInfo.i)
        varInfo.updated <- combineGRanges(varInfo.updated, varInfo.updated.i)
      }
    } else {
      varInfo.updated <- combineGRanges(varInfo.updated, varInfo.i)
    }
  }
  missVarInfo <- unlist(GRangesList(missVarInfo))
  whitelist.only <- checkWhitelist(missVarInfo, includeCosmic=FALSE)
  names(missVarInfo) <- NULL
  missVarInfo <- as.data.frame(missVarInfo)
  if (!is(varInfo.updated, 'VRanges')) {
    varInfo.updated <- unlist(GRangesList(varInfo.updated))
  }
  names(varInfo.updated) <- NULL
  varInfo.updated <- as.data.frame(varInfo.updated)
  varInfo.org$Patched <- FALSE
  varInfo.updated$sampleNames <- as.character(varInfo.updated$sampleNames)
  commCols <- colnames(varInfo.org)[colnames(varInfo.org) %in% colnames(missVarInfo)]
  missVarInfo <- missVarInfo[, commCols,drop=FALSE]
  missVarInfo$Note <- varInfo.org$Note <- ''
  missVarInfo$Note[missVarInfo$VariantFreq < AF.th.lod & !whitelist.only | missVarInfo$VariantFreq < AF.th.wl] <- 'Reference only'
  missVarInfo$Patched <- TRUE
  if (missVarOnly) {
    varInfo.all <- missVarInfo
  } else {
    if (updateVariant) {
      varInfo.org <- as.data.frame(lapply(varInfo.org, as.character))
      varInfo.updated <- as.data.frame(lapply(varInfo.updated, as.character))
      rownames(varInfo.org) <- varInfo.org$ID
      rownames(varInfo.updated) <- varInfo.updated$ID
      commIDs <- intersect(rownames(varInfo.org), rownames(varInfo.updated))
      commCols <- colnames(varInfo.org)[colnames(varInfo.org) %in% colnames(varInfo.updated)]
      varInfo.org[commIDs,commCols] <- varInfo.updated[commIDs,commCols,drop=FALSE]
    }
    varInfo.all <- combineGRanges(varInfo.org[varInfo.org$PatientID %in% uniPatient,,drop=FALSE], missVarInfo)
    varInfo.all$ID <- paste(varInfo.all$PatientID, varInfo.all$VariantID, sep=':')
    dup.id <- varInfo.all$ID[duplicated(varInfo.all$ID)]
    varInfo.all$concordant <- varInfo.all$ID %in% dup.id
    if (!is.null(varInfo.all$highFrequent.inbatch)) {
      var.count <- table(varInfo.all$VariantID)
      varInfo.all$highFrequent.inbatch <- var.count[varInfo.all$VariantID]
    }
  }
  if (!is.null(qcInfo)) {
    rownames(qcInfo) <- qcInfo$SampleID
    otherCols <- c("PatientID", "externalSampleID", "Specimen_type", "trialVisitNum")
    otherCols <- otherCols[otherCols %in% colnames(qcInfo)]
    varInfo.all <- varInfo.all[, !(colnames(varInfo.all) %in% otherCols)]
    varInfo.all <- cbind(varInfo.all, qcInfo[varInfo.all$sampleNames, otherCols])
  }
  return(varInfo.all)
}

