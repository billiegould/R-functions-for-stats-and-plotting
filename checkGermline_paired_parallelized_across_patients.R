## a version of matched germline analysis that is parallized with one core per patient,
# also returns crititcal error messages such as missing germline samples

### 06/28/2023
## BGould

## furhter testing needed

library(NGSutilities)
library(DeepSea)
package.version("DeepSea")
library(openxlsx)
library(dplyr)
options(stringsAsFactors = FALSE)
`%!in%` = Negate(`%in%`)

workDir <- '~/DevHome/Moffitt_NMIBC/post_analysis/Urine-WES+/'
setwd(workDir)
sink(sprintf("sessionInfo_%s.txt", Sys.Date()))
sessionInfo()
sink()

######################
## a wrapper around mcmapply that will display critical warning and error messages
# that are normally hidden when running mclapply:
safe_mcmapply <- function(X, FUN, mc.cores, stop.on.error=T, ...){
  fun <- function(x){
    res_inner <- tryCatch({
      withCallingHandlers(
        expr = {
          FUN(x, ...)
        }, 
        warning = function(e) {
          message_parallel(trimws(paste0("WARNING [element ", x,"]: ", e)))
          # this line is required to continue FUN execution after the warning
          invokeRestart("muffleWarning")
        },
        error = function(e) {
          message_parallel(trimws(paste0("ERROR [element ", x,"]: ", e)))
        }
      )},
      error = function(e){
        # error is returned gracefully; other results of this core won't be affected
        return(e)
      }
    )
    return(res_inner)
  }
  
  res <- mcmapply(X, fun, mc.cores=mc.cores, ...)
  failed <- sapply(res, inherits, what = "error")
  if (any(failed == T)){
    error_indices <- paste0(which(failed == T), collapse=", ")
    error_traces <- paste0(lapply(res[which(failed == T)], function(x) x$message), collapse="\n\n")
    error_message <- sprintf("Elements with following indices failed with an error: %s. Error messages: \n\n%s", 
                             error_indices,
                             error_traces)
    if (stop.on.error)
      stop(error_message)
    else
      warning(error_message, "\n\n### Errors will be ignored ###")
  }
  return(res[!failed])
}
#' Function which prints a message using shell echo; useful for printing messages from inside mclapply when running in Rstudio
message_parallel <- function(...){
  system(sprintf('echo "\n%s\n"', paste0(..., collapse="")))
}
##########################

## germline marking for variant short file
wkdir <- "/prednet/data03/DevHome/bgould/Moffitt_NMIBC/post_analysis/Urine-WES+"
setwd(wkdir)

qcInfo = read.csv(qcInfo.file)
qcInfo$Specimen_type = gsub("Buffy_Coat", "Buffy Coat", qcInfo$Specimen_type) #[1] "Urine_Supernatant" "Buffy_Coat" "Buffy Coat"    
qcInfo$Specimen_type <- tolower(qcInfo$Specimen_type)
#table(qcInfo$Specimen_type, qcInfo$analysis_type)
uniPatient <- unique(qcInfo$PatientID)

pversion=c("1.7.0","1.7.1")
germlineSample.type <- 'buffy coat'
tumorSample.type <-  'urine_supernatant'
mc.cores <- length(uniPatient)

############
#### germline marking for variant short file, one core per patient
###########
varInfo.short = read.csv(varInfo.short.file)
res.short <- safe_mcmapply(uniPatient, 
  function(patient.i){
    #varfile.suffix.i <- 'Variant_short'
    #print(varfile.suffix.i)
  setwd(wkdir)
  qcInfo.i <- qcInfo[which(qcInfo$PatientID == patient.i),]
  analysisType.i <- qcInfo.i$analysis_type[qcInfo.i$Specimen_type == tumorSample.type][1]
  sample.p.i <- qcInfo.i$SampleID[qcInfo.i$Specimen_type == tumorSample.type]
  if (length(sample.p.i) == 0) {
    warning(paste('No tumor sample was found for patient', patient.i))
    return(GRanges())
  }
  sample.g.i <- qcInfo.i$SampleID[qcInfo.i$Specimen_type == germlineSample.type]
  if (length(sample.g.i) == 0) {
    warning(paste('No germline sample was found for patient', patient.i))
    file.g.i <- rawVariants.g.i <- coverage.g.i <- NULL
  } else {
    if (length(sample.g.i) > 1) sample.g.i <- sample.g.i[1]
    print(sample.g.i)
    file.g.i <- findFileBySample(sample.g.i, file.suffix = '.coverage.RData', pversion=pversion)$filelist ##
    if (!is.null(file.g.i)){
      coverage.g.i <- get(load(file.g.i))
    }else{
      file.g.i <- findFileBySample(sample.g.i, file.suffix = '.coverage.rds', pversion=pversion)$filelist
      coverage.g.i <- readRDS(file.g.i)
    }
    
    file.g.i <- findFileBySample(sample.g.i, file.suffix='_dist2end.csv', pversion=pversion)$filelist ##
    if (!is.null(file.g.i)){
      rawVariants.g.i <- read.csv(file.g.i)
    }else{
      file.g.i <- findFileBySample(sample.g.i, file.suffix='_rawVariants.RData', pversion=pversion)$filelist
      rawVariants.g.i <- get(load(file.g.i))
    }
  }
  varInfo.pi = GRanges()
  varInfo.tmb.pi = GRanges()
  for (sample.ij in sample.p.i) {
    cat(paste('Processing sample', sample.ij, '\n'))
    varInfo.p.ij <- varInfo.short[which(varInfo.short$sampleNames == sample.ij),,drop=FALSE]
    if (!is.null(rawVariants.g.i) && nrow(varInfo.p.ij) > 0) {
      varInfo.p.ij <- checkGermline.paired(varInfo.p.ij, rawVariants.g.i, coverage.g.i)
    }
    varInfo.pi = combineGRanges(varInfo.pi, asGRanges(varInfo.p.ij))
  } # end samples
  return("short"=varInfo.pi) 
  }, 
  varInfo.short=varInfo.short, qcInfo=qcInfo, mc.cores=mc.cores, stop.on.error=F,
  pversion=pversion, germlineSample.type=germlineSample.type, tumorSample.type=tumorSample.type)
# double check errors
sapply(res.short, inherits, what = "try-error")
#str(res.short[1]) # view errors
vars.short <- do.call(plyr::rbind.fill, res.short)
stopifnot(all(uniPatient %in% unique(vars.short$PatientID)))
stopifnot(nrow(vars.short) == nrow(varInfo.short))
run.date <- "2023-06-29" # input file date
varInfo.short.file.marked <- sub('.csv', paste0('_', tumorSample.type, '_marked.csv'), varInfo.all.file)
write.csv(vars.short, varInfo.all.file.marked)

#########
## germline marking for all variants and tmb file
#########
varInfo.all <- read.csv(varInfo.all.file)
res.all <- mcmapply(uniPatient, function(patient.i){
  #varfile.suffix.i <- 'Variant'
  #print(varfile.suffix.i)
  setwd(wkdir)
  qcInfo.i <- qcInfo[which(qcInfo$PatientID == patient.i),]
  analysisType.i <- qcInfo.i$analysis_type[qcInfo.i$Specimen_type == tumorSample.type][1]
  sample.p.i <- qcInfo.i$SampleID[qcInfo.i$Specimen_type == tumorSample.type]
  if (length(sample.p.i) == 0) {
    warning(paste('No tumor sample was found for patient', patient.i))
    return(c(GRanges(), GRanges()))
  }
  sample.g.i <- qcInfo.i$SampleID[qcInfo.i$Specimen_type == germlineSample.type]
  if (length(sample.g.i) == 0) {
    warning(paste('No germline sample was found for patient', patient.i))
    file.g.i <- rawVariants.g.i <- coverage.g.i <- NULL
  } else {
    if (length(sample.g.i) > 1) sample.g.i <- sample.g.i[1]
    print(sample.g.i)
    file.g.i <- findFileBySample(sample.g.i, file.suffix = '.coverage.RData', pversion=pversion)$filelist ##
    if (!is.null(file.g.i)){
      coverage.g.i <- get(load(file.g.i))
    }else{
      file.g.i <- findFileBySample(sample.g.i, file.suffix = '.coverage.rds', pversion=pversion)$filelist
      coverage.g.i <- readRDS(file.g.i)
    }
    
    file.g.i <- findFileBySample(sample.g.i, file.suffix='_dist2end.csv', pversion=pversion)$filelist ##
    if (!is.null(file.g.i)){
      rawVariants.g.i <- read.csv(file.g.i)
    }else{
      file.g.i <- findFileBySample(sample.g.i, file.suffix='_rawVariants.RData', pversion=pversion)$filelist
      rawVariants.g.i <- get(load(file.g.i))
    }
  }
  varInfo.pi = GRanges()
  varInfo.tmb.pi = GRanges()
  for (sample.ij in sample.p.i) {
    cat(paste('Processing sample', sample.ij, '\n'))
    varInfo.p.ij <- varInfo.all[which(varInfo.all$sampleNames == sample.ij),,drop=FALSE]
    if (!is.null(rawVariants.g.i) && nrow(varInfo.p.ij) > 0) {
      varInfo.p.ij <- checkGermline.paired(varInfo.p.ij, rawVariants.g.i, coverage.g.i)
    }
    varInfo.pi = combineGRanges(varInfo.pi, asGRanges(varInfo.p.ij))
    ## tmb filtering
    if (!is.null(varInfo.p.ij$VariantType)) {
        varInfo.p.f.ij <- varInfo.p.ij[!grepl('(germline)|(chip)', varInfo.p.ij$VariantType, ignore.case=TRUE) | is.na(varInfo.p.ij$VariantType),]
      } else {
        varInfo.p.f.ij <- varInfo.p.ij
      }
    ## turn off germlineFiltering
    #germlineFiltering <- ifelse(grepl("_rawVariants.RData", file.g.i), TRUE, FALSE) # arg no longer exists
    varInfo.tmb.p.ij <- tmbFiltering(varInfo.p.f.ij, analysisType=analysisType.i)
    varInfo.tmb.pi = combineGRanges(varInfo.tmb.pi, asGRanges(varInfo.tmb.p.ij))
  } # end samples
  return(list("all"=data.frame(varInfo.pi), "tmb"=data.frame(varInfo.tmb.pi)))
}, varInfo.all=varInfo.all, qcInfo=qcInfo, mc.cores=mc.cores, stop.on.error=F, 
pversion=pversion, germlineSample.type=germlineSample.type, tumorSample.type=tumorSample.type)

# double check errors
sapply(res.all, inherits, what = "try-warrning")
#str(res.all[1]) # view errors

vars.tmb <- sapply(res.all, function(x) x$tmb)
vars.tmb <- do.call(plyr::rbind.fill, vars.tmb)
stopifnot(all(uniPatient %in% unique(vars.tmb$PatientID)))
varInfo.tmb.file <- "./WOP00873_WOP00875_MCC-NMIBC_WES_urine_TMB_all_2023-06-29.csv"
varInfo.tmb.file.marked <- sub('.csv', paste0('_',tumorSample.type, '_marked.csv'), varInfo.tmb.file) 
write.csv(vars.tmb, varInfo.tmb.file.marked) 

vars.all <- sapply(res.all, function(x) x$all)
vars.all <- do.call(plyr::rbind.fill, vars.all)
stopifnot(all(uniPatient %in% unique(vars.all$PatientID)))
stopifnot(nrow(res.all)==nrow(varInfo.all))
varInfo.all.file.marked <- sub('.csv', paste0('_', tumorSample.type, '_marked.csv'), varInfo.all.file)
write.csv(vars.all, varInfo.all.file.marked)