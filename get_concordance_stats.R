library(eulerr)
source("~/Desktop/puffin/R/helper_functions.R")

## for an input dataframe with sample variants, calculate the number of variants per sample and the
## number of matching variants with the reference sample type in the set. 
##' @param variant.data a data frame of SNV and/or CNV variants containing PatientID, SampleIDs, StudyVisit, SampleType and VariantID
##' @param ref_sample_type a string in the form: SampleType_StudyVisit, eg "Tissue_Preop". Indicates the reference sample type to compare other samples against within patients
##' @param df_samples dataframe containing PatientIDs and SampleIDs and clinical info to include in the output table
##' @param out.file optional path to file for writing the concordance output.
##' 
##' @return a dataframe with one row per sample summarizing the number of variants and percent of variants matching the reference sample
##' 
get_concordance_stats <- function(variant.data, ref_sample_type, df_samples, out.file=NA){
  df_samples = standardize_names(df_samples, input.type="samples")
  variant.data = standardize_names(variant.data)
  stopifnot(all(c("PatientID", "SampleID.short", "StudyVisit", "SampleType", "VariantID") %in% names(variant.data)))
  print(unique(df_samples$StudyVisit))
  print(unique(df_samples$SampleType))
  df_samples = df_samples %>% mutate(sample_type=paste0(StudyVisit,"_",SampleType)) 
  sample_data = list()
  # for each patient, extract reference variants
  for (pid in unique(df_samples$PatientID)){
    #print(pid)
    df_pid = df_samples %>% filter(PatientID==pid, !(tolower(SampleType) %in% c("buffy coat","pbmc")))
    ref_sid = df_pid %>% filter(sample_type==ref_sample_type)
    if (nrow(ref_sid) >0){
      if (nrow(ref_sid)>1){
        print(paste0("WARN: more than one ref sample, using first one: ",pid))
      }
      ref_sid = ref_sid$SampleID.short[1]
      #print(ref_sid)
      df_ref_vars = variant.data %>% filter(SampleID.short==ref_sid)
      n_ref_vars = nrow(df_ref_vars)
      for (sid.short in unique(df_pid$SampleID.short)){
        df_sample_vars = variant.data %>% filter(SampleID.short==sid.short) %>% distinct()
        n_sample_vars = nrow(df_sample_vars)
        n_concord_vars = nrow(df_sample_vars %>% filter(VariantID %in% df_ref_vars$VariantID))
        sample_data[[sid.short]] <- c(n_ref_vars, n_sample_vars, n_concord_vars)
        #print(length(sample_data))
      }
    }else{
      print(sprintf("WARN: No ref. sample for: %s", pid))
      for (sid.short in unique(df_pid$SampleID.short)){
        df_sample_vars = variant.data %>% filter(SampleID.short==sid.short) %>% distinct()
        n_sample_vars = nrow(df_sample_vars)
        sample_data[[sid.short]] <- c(NA, n_sample_vars, NA)
      }
    }
  }
  max_mafs = do.call(rbind, sample_data)
  # convert to df
  df_max_mafs <- cbind(rownames(max_mafs), data.frame(max_mafs, row.names=NULL))
  names(df_max_mafs) <- c("SampleID.short","n_ref_vars", "n_sample_vars", "n_concord_vars")
  df_concord_counts = df_max_mafs %>% left_join(df_samples, by="SampleID.short") %>%
    mutate(n_concord_vars = as.numeric(n_concord_vars), n_sample_vars = as.numeric(n_sample_vars)) %>%
    mutate(n_unique_vars= n_sample_vars - n_concord_vars)
  stopifnot(all(!(duplicated(df_concord_counts$SampleID.short))))
  df_concord_counts$pct.concord = (df_concord_counts$n_concord_vars/df_concord_counts$n_sample_vars) * 100
  df_concord_counts = df_concord_counts %>% mutate(pct.concord = ifelse((n_ref_vars==0 & n_sample_vars==0), 100.0, pct.concord)) %>%
                                            mutate(pct.concord = ifelse((n_ref_vars>0 & n_sample_vars==0), 0.0, pct.concord))
  df_concord_counts = df_concord_counts %>% select(PatientID, SampleID.short, SampleType, StudyVisit,
                                                   n_ref_vars, n_sample_vars, n_concord_vars, n_unique_vars,
                                                   pct.concord, sample_type)
  ## summarize ref samples
  ref_counts = df_concord_counts %>% filter(sample_type==ref_sample_type)
  stopifnot(identical(ref_counts$n_sample_vars,ref_counts$n_concord_vars))
  print(sprintf("reference sample: %s", ref_sample_type))
  print(sprintf("median number reference variants: %s", median(ref_counts$n_sample_vars)))
  print(sprintf("min: %s", min(ref_counts$n_sample_vars)))
  print(sprintf("max: %s", max(ref_counts$n_sample_vars)))
  print(sprintf("sd: %s", sd(ref_counts$n_sample_vars)))
  ## summarize mrd samples
  nonref_counts = df_concord_counts %>% filter(sample_type!=ref_sample_type)
  print(sprintf("MRD sample: %s", unique(df_concord_counts$sample_type[df_concord_counts$sample_type!=ref_sample_type])))
  print(sprintf("median number MRD variants: %s", median(nonref_counts$n_sample_vars)))
  print(sprintf("min: %s", min(nonref_counts$n_sample_vars)))
  print(sprintf("max: %s", max(nonref_counts$n_sample_vars)))
  print(sprintf("sd: %s", sd(nonref_counts$n_sample_vars)))
  
  if (!(is.na(out.file))){
    df_out = apply(df_concord_counts, 2, as.character)
    write.csv(df_out, file=out.file, row.names=FALSE)
  }
  
  return(df_concord_counts)
}


### function for making venn diagrams of overall sample concordance
##' @param variants a dataframe of variants (snv and/or cnvs) with variant ids and sample ids
##' @param samples optional dataframe of sample ids to include 
##' @param selectors a list of two or more strings indicating "StudyVisit_SampleType" combinations to include
##' @param colors optional dictionary matching selectors to colors.
##' @param cex optinal sizing parameter for the plot text and legend
##' 
##' @return a eulerr plot object
concordance_venn <- function(variants, samples=NULL, selectors=NULL, 
                             colors=NULL, cex=1.5, fontsize = 10){
  print(table(variants$SampleType, variants$StudyVisit))
  variants = standardize_names(variants)
  print(table(variants$SampleType, variants$StudyVisit))
  if (!is.null(samples)){
    samples = standardize_names(samples, input.type="samples")
    variants = variants %>% filter(SampleID.short %in% samples$SampleID.short) 
    variants = merge.combine(variants, samples, join.type = "left",
                             join.cols = "SampleID.short", priority = "right", 
                             warn=FALSE)
    n.patients = length(unique(variants$PatientID))
  }else{
    stop("Missing sample list input.")
  }
  print(table(variants$SampleType, variants$StudyVisit))
  variants = variants %>% mutate(label=paste0(PatientID,VariantID),
                                 Selector=paste0(StudyVisit,"_",SampleType)) %>% 
    select(Selector, label) %>%
    distinct()
  if (is.null(colors)){
    colors = get_random_color_dict(selectors)
  }
  if (any(selectors %!in% names(colors))){
    check.missing(list.ref = names(colors), list.test = selectors)
    stop("Missing colors for selectors.Colors should be a named list")
  }
  plot_list = list()
  print(table(variants$Selector))
  for (sel in selectors){
    df = variants %>% filter(Selector==sel)
    if (nrow(df)==0){
      check.missing(list.ref=variants$Selector, list.test=selectors)
      stop(glue("concordance_venn:: WARN: No variants present for {sel}"))
    }
    plot_list[[sel]]=df$label
  }
  # fix this: stopifnot(length(unlist(plot_list))==nrow(variants))
  v <- euler(plot_list) # function converts list to alphabetical order by label
  eulerr_options(padding = unit(5, "mm"))
  plot(v, quantities = list(type = "counts", cex=cex),
          fills = list(fill = colors, alpha=0.6),
          legend = list(plot=TRUE, fontsize = fontsize),
          main=glue("N Patients: {n.patients}"))
  #show(plt)
  #return(c("plot"=NULL,"colors"=colors)) # func does not return plot object, ## fix.
}