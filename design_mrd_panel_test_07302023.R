## annotation of design.mrd.panel v. 1.5.39
## code edited by BG - 7/23/2023

design.mrd.panel.test <- function (variants, fusionInfo = NULL, sampleID = NULL, tumorFraction = NULL, 
          minVariant = 4, markColumn = "Selected", controlProbe = NULL, 
          totalProbeNum = 50, pversion = NULL, multiSample = FALSE, 
          snpBC = "", snpBCInfo = "SNP_sampleBC_5_190616.csv", txdb = "TxDb.Hsapiens.UCSC.hg19.knownGene", 
          bsgenome = "BSgenome.Hsapiens.UCSC.hg19") 
{
  if (!is(variants, "VRanges"))
    suppressWarnings(variants <- asVRanges(variants))
  if (is.null(sampleID)) {
    sampleID <- unique(sampleNames(variants))
  }
  print("running test . . ")
  # TODO
  if (length(sampleID) > 1) {
  #   if (multiSample) { 
  #     probeInfo.all <- targetRegion.all <- NULL
  #     for (i in 1:length(sampleID)) {
  #       sample.i <- sampleID[i]
  #       panelInfo.i <- design.mrd.panel(variants, fusionInfo = fusionInfo, 
  #                                       sampleID = sample.i, tumorFraction = tumorFraction, 
  #                                       minVariant = minVariant, markColumn = markColumn, 
  #                                       controlProbe = controlProbe, totalProbeNum = totalProbeNum, 
  #                                       pversion = pversion, multiSample = multiSample, 
  #                                       snpBC = snpBC, snpBCInfo = snpBCInfo, txdb = txdb, 
  #                                       bsgenome = bsgenome)
  #       if (i > 1) {
  #         probeInfo.all <- rbind(probeInfo.all, panelInfo.i$probeInfo)
  #         targetRegion.i <- panelInfo.i$targetRegion
  #       }
  #       else {
  #         probeInfo.all <- panelInfo.i$probeInfo
  #         targetRegion.all <- panelInfo.i$targetRegion
  #       }
  #     }
  #     if (any(duplicated(probeInfo.all$VariantID))) {
  #       dupVid <- unique(probeInfo$VariantID[duplicated(probeInfo$VariantID)])
  #     }
  #     return(list(probeInfo = probeInfo.all, targetRegion = targetRegion.all))
  #   }
  #   else {
  stop("More than 1 samples were found in variants. Please provide the variant ID for processing!")
  #   }#end if multiSample
  } # end if nsampleID >1
  
  if (length(sampleID) == 1 && !is.na(sampleID)) {
    variants <- variants[sampleNames(variants) == sampleID]
  }
  if (snpBC == "") {
    snpbc.Info <- getSampleSNPBC(sampleID, pversion = pversion[1])
    if (length(snpbc.Info) > 0) {
      snpBC <- snpbc.Info$snpBCs
    }
  }
  variants <- variants[order(variants$VariantFreq, decreasing = TRUE)]
  if (is.null(tumorFraction) || is.na(tumorFraction)) {
    tumorFraction <- estimateTF.maf(variants, maxAF.th = NULL)
  }
  if (is.null(variants$VariantID)){
    variants$VariantID <- getVariantIDs(variants)
  } 
  names(variants) <- variants$VariantID
  selVarIDs.other <- selVarIDs.core <- selVarIDs.add <- NULL
  if (!is.null(markColumn)) {
    selVarIDs.all <- variants$VariantID[!(toupper(values(variants)[[markColumn]]) %in% 
                                            c("FALSE", ""))]
    selVarIDs.probe <- variants$VariantID[toupper(values(variants)[[markColumn]]) %in% 
                                            c("TRUE", "YES", "ADD")]
    selVarIDs.add <- variants$VariantID[toupper(values(variants)[[markColumn]]) %in% 
                                          c("ADD")]
    selVarIDs.other <- variants$VariantID[(toupper(values(variants)[[markColumn]]) %in% 
                                             c("FALSE", ""))]
    selVarIDs.core <- variants$VariantID[(toupper(values(variants)[[markColumn]]) %in% 
                                            c("CORE"))]
    if (length(selVarIDs.probe) > totalProbeNum) {
      selVarIDs.probe <- selVarIDs.probe[1:totalProbeNum] # edit
    }
    selVarIDs.probe <- c(selVarIDs.probe, selVarIDs.core)
  }else{
    selVarIDs.all <- selVarIDs.probe <- variants$VariantID
  }
  if (length(selVarIDs.probe) < minVariant) {
    warning(paste("At least", minVariant, "variants are required for MRD panel design!"))
    return(NULL)
  }
  snpVariantIDs <- NULL
  if (!is.null(snpBCInfo) && is.character(snpBCInfo)) {
    targetInfoDir <- options("targetInfoDir")
    snpBCInfo <- file.path(targetInfoDir, snpBCInfo)
    snpBCInfo <- read.csv(snpBCInfo, as.is = TRUE, stringsAsFactors = FALSE)
    snpVariantIDs <- checkVariantID(snpBCInfo$VariantID)
    selVarIDs.other <- unique(c(selVarIDs.other, snpVariantIDs))
  }
  if (!is.null(controlProbe)) {
    selVarIDs.probe <- c(selVarIDs.probe, controlProbe)
  }
  probeInfo.gr <- variantID2probeDesign(selVarIDs.probe)
  if (!is.null(fusionInfo)) {
    fusion.bk <- fusionJunction2breakpoints(fusionInfo)
    fusion.bk$ID <- paste(sampleNames(fusion.bk), fusion.bk$VariantID_junction, 
                          sep = ":")
    fusion.bk$ID <- gsub(":[lr]", "", fusion.bk$ID)
    probe.bk <- getFusionJunctionDNA(fusionInfo)
    rownames(probe.bk) <- probe.bk$id
    fusion.bk$sequence <- probe.bk[fusion.bk$ID, "ProbeSeq"]
    fusion.bk$VariantID <- names(fusion.bk)
    fusion.bk$mutationID <- fusion.bk$VariantID_junction
    fusion.bk$VariantFreq <- fusion.bk$AF
    fusion.bk$lowRank <- !is.na(values(fusion.bk)[, "repeat"])
    fusion.bk$filteredDSCnt <- fusion.bk$ds_molecule_count
    fusion.bk$filteredCnt <- altDepth(fusion.bk)
    selVarIDs.all <- c(selVarIDs.all, fusion.bk$VariantID_junction)
    selCols <- c("VariantName", "diversity", "partner", 
                 names(values(probeInfo.gr)))
    probeInfo.gr <- combineGRanges(probeInfo.gr, fusion.bk)
    values(probeInfo.gr) <- values(probeInfo.gr)[, names(values(probeInfo.gr)) %in% 
                                                   selCols]
  }
  selVarIDs.all <- unique(selVarIDs.all)
  if (length(selVarIDs.all) < totalProbeNum) {
    probeInfo.gr.l.60 <- variantID2probeDesign(selVarIDs.probe, 
                                               probeShift = -60)
    probeInfo.gr.r.60 <- variantID2probeDesign(selVarIDs.probe, 
                                               probeShift = 60)
    probeInfo.gr.l.30 <- variantID2probeDesign(selVarIDs.probe, 
                                               probeShift = -30)
    probeInfo.gr.r.30 <- variantID2probeDesign(selVarIDs.probe, 
                                               probeShift = 30)
    repeat.time <- floor(totalProbeNum/length(selVarIDs.probe))
    if (length(selVarIDs.core) > 0) {
      remain <- (totalProbeNum + length(selVarIDs.core) * 
                   2)%%length(selVarIDs.probe)
    }
    else {
      remain <- totalProbeNum%%length(selVarIDs.probe)
    }
    if (repeat.time > 2) {
      probeInfo.gr.all <- c(probeInfo.gr, probeInfo.gr.l.60, 
                            probeInfo.gr.r.60)
    }
    else if (repeat.time == 2) {
      probeInfo.gr.all <- c(probeInfo.gr[1:remain], probeInfo.gr.l.60[1:remain], 
                            probeInfo.gr.r.60[1:remain])
      probeInfo.gr.all <- c(probeInfo.gr.all, probeInfo.gr.l.30[-(1:remain)], 
                            probeInfo.gr.r.30[-(1:remain)])
    }
    else if (repeat.time == 1) {
      probeInfo.gr.all <- c(probeInfo.gr.l.30[1:remain], 
                            probeInfo.gr.r.30[1:remain], probeInfo.gr[-(1:remain)])
    }
    else {
      probeInfo.gr.all <- probeInfo.gr
    }
    if (length(probeInfo.gr.all) < totalProbeNum) {
      restProbeNum <- totalProbeNum - length(probeInfo.gr.all)
      if (length(selVarIDs.other) > restProbeNum) 
        selVarIDs.other <- selVarIDs.other[1:restProbeNum]
      probeInfo.gr <- variantID2probeDesign(selVarIDs.other)
      if (length(selVarIDs.other) < restProbeNum & length(selVarIDs.other) > 0) {
        probeInfo.gr.l.60 <- variantID2probeDesign(selVarIDs.other, 
                                                   probeShift = -60)
        probeInfo.gr.r.60 <- variantID2probeDesign(selVarIDs.other, 
                                                   probeShift = 60)
        probeInfo.gr.l.30 <- variantID2probeDesign(selVarIDs.other, 
                                                   probeShift = -30)
        probeInfo.gr.r.30 <- variantID2probeDesign(selVarIDs.other, 
                                                   probeShift = 30)
        repeat.time <- floor(restProbeNum/length(selVarIDs.other))
        remain <- restProbeNum%%length(selVarIDs.other)
        if (repeat.time > 2) {
          probeInfo.gr.other <- c(probeInfo.gr, probeInfo.gr.l.60, 
                                  probeInfo.gr.r.60)
        }
        else if (repeat.time == 2) {
          probeInfo.gr.other <- c(probeInfo.gr[1:remain], 
                                  probeInfo.gr.l.60[1:remain], probeInfo.gr.r.60[1:remain])
          probeInfo.gr.other <- c(probeInfo.gr.other, 
                                  probeInfo.gr.l.30[-(1:remain)], probeInfo.gr.r.30[-(1:remain)])
        }
        else if (repeat.time == 1) {
          probeInfo.gr.other <- c(probeInfo.gr.l.30[1:remain], 
                                  probeInfo.gr.r.30[1:remain], probeInfo.gr[-(1:remain)])
        }
        suppressWarnings(probeInfo.gr.all <- combineGRanges(probeInfo.gr.all, 
                                                            probeInfo.gr.other))
      }
      else {
        suppressWarnings(probeInfo.gr.all <- combineGRanges(probeInfo.gr.all, 
                                                            probeInfo.gr))
      }
      selVarIDs.probe <- c(selVarIDs.probe, selVarIDs.other)
    }
  }
  else {
    probeInfo.gr.all <- probeInfo.gr
  }
  probeInfo.gr.all.old <- probeInfo.gr.all
  probeInfo.gr.all <- sort(probeInfo.gr.all)
  probeInfo.gr.all$sampleNames <- sampleID
  names(probeInfo.gr.all) <- probeInfo.gr.all$mutationID
  
  selInd = NULL
  if (any(probeInfo.gr.all$mutationID %in% names(variants))){
    selInd <- which(probeInfo.gr.all$mutationID %in% names(variants))
  }
  #print(selInd)
  mutids <- probeInfo.gr.all[selInd]$mutationID
  #print(mutids)
  probeInfo.gr.all$PatientID <- probeInfo.gr.all$pvalue <- probeInfo.gr.all$odds.ratio <- probeInfo.gr.all$AF.Baseline <- NA
  if (is.null(probeInfo.gr.all$VariantType)) 
    {probeInfo.gr.all$VariantType <- NA}
  if (is.null(probeInfo.gr.all$VariantFreq)) 
   {probeInfo.gr.all$VariantFreq <- NA}
  if (is.null(probeInfo.gr.all$lowRank)) 
    {probeInfo.gr.all$lowRank <- NA}
  if (is.null(probeInfo.gr.all$filteredDSCnt)) 
    {probeInfo.gr.all$filteredDSCnt <- NA}
  if (is.null(probeInfo.gr.all$filteredCnt)) 
    probeInfo.gr.all$filteredCnt <- NA
  if (is.null(probeInfo.gr.all$altDepth)) 
    {probeInfo.gr.all$altDepth <- NA}
  
  probeInfo.gr.all[mutids]$VariantFreq <- signif(variants[mutids]$VariantFreq, 
                                                 5)
  if (!is.null(variants$filteredCnt)) {
    probeInfo.gr.all[mutids]$filteredCnt <- signif(variants[mutids]$filteredCnt, 
                                                   5)
    probeInfo.gr.all[mutids]$filteredDSCnt <- signif(variants[mutids]$filteredDSCnt, 
                                                     5)
  }
  fill_ <- function(col){
    #print(col)
    stopifnot(!is.null(names(probeInfo.gr.all)))
    new <- values(probeInfo.gr.all)
    if (!is.null(mcols(variants)[[col]])) {
      key <- values(variants[mutids])[[col]]
      names(key) <- rownames(values(variants[mutids]))
      new[[col]] <- NA
      new[mutids,][[col]] <- key[rownames(new[mutids,])]
    }
    return(new)
  }
  values(probeInfo.gr.all) <- fill_("VariantType")
  values(probeInfo.gr.all) <- fill_("AF.Baseline")
  values(probeInfo.gr.all) <- fill_("odds.ratio")
  values(probeInfo.gr.all) <- fill_("pvalue")
  values(probeInfo.gr.all) <- fill_("lowRank")
  values(probeInfo.gr.all) <- fill_("PatientID")
  
  names(probeInfo.gr.all) <- NULL
  probeInfo.gr.all <- probeInfo.gr.all[order(probeInfo.gr.all$VariantID, 
                                             is.na(probeInfo.gr.all$VariantType))]
  suppressWarnings(targetRegion <- prepareTargetRegion(probeInfo.gr.all, 
                                                       txdb = txdb, bsgenome = bsgenome))
  selInd <- which(targetRegion$mutationID %in% names(variants))
  if (length(selInd) == 0)  {selInd <- c()}
  #print(selInd)
  targetRegion$altDepth <- targetRegion$totalDepth <- NA
  targetRegion$altDepth[selInd] <- altDepth(variants[targetRegion$mutationID[selInd]])
  targetRegion$totalDepth[selInd] <- totalDepth(variants[targetRegion$mutationID[selInd]])
  targetRegion$tumorFraction <- signif(tumorFraction, 5)
  targetRegion$SNPBC <- snpBC
  
  if (!is.null(snpVariantIDs)) {
    targetRegion$mutationID[targetRegion$mutationID %in% 
                              snpVariantIDs] <- ""
  }
  targetRegion$VariantType[targetRegion$VariantID %in% #c(#snpVariantIDs, 
      controlProbe] <- "Control"
  if (!is.null(selVarIDs.add)) {
    targetRegion$VariantType[targetRegion$VariantID %in% 
                               selVarIDs.add] <- "Add"
  }
  probeInfo <- as.data.frame(probeInfo.gr.all)
  if (any(probeInfo$VariantID %in% selVarIDs.core)) {
    probeInfo <- rbind(probeInfo[!(probeInfo$VariantID %in% 
                                     selVarIDs.core), , drop = FALSE], probeInfo[(probeInfo$VariantID %in% 
                                                                                    selVarIDs.core), , drop = FALSE])
  }
  selCols <- c("sequence", "seqnames", "start", "end", "sampleNames", 
               "mutationID")
  probeInfo <- probeInfo[, selCols, drop = FALSE]
  totalProbeNum <- totalProbeNum + length(controlProbe)
  if (nrow(probeInfo) >= totalProbeNum) {
    probeInfo <- probeInfo[1:totalProbeNum, , drop = FALSE]
  }else {
    index <- rep(1:nrow(probeInfo), round(totalProbeNum/2))
    probeInfo <- probeInfo[index[1:totalProbeNum], , drop = FALSE]
  }
  # check columns filled correctly
  stopifnot(is.character(targetRegion$VariantType)) 
  stopifnot(is.numeric(targetRegion$filteredCnt))
  stopifnot(is.numeric(targetRegion$lowRank) | is.logical(targetRegion$lowRank))
  
  return(list(probeInfo = probeInfo, targetRegion = targetRegion))
}
