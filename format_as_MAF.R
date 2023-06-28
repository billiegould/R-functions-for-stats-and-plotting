source("~/Desktop/puffin/R/helper_functions.R")

## add required MAF file format columns to Predicine variant outputs and filter out non-somatic variants. 
## Optionaly restrict output to only MAF required column for use with MAF tools analysis.
##' @param variant.data a data frame of SNV or CNV variants with sample IDs.
##' @param variant.type "snv" or "cnv" input type
##'  @param df_samples optional data frame of sample ids and clinical data to include in output
##'  @param maf.output optional whether to only output required MAF cols. default: TRUE 
##'  @param sid.format optionally format for standardizing SampleIDs. defaut:"none" 
##'  @param warn optionally turn on warnings about NAs in merged coumns. defaut: TRUE
##' 
##' @return a MAF formatted data frame containing somatic variants (germline CNVs are removed if present)
format_as_MAF <- function(variant.data, df_samples=NULL, variant.type="unknown", maf.output=TRUE,
                          sid.format="none", warn=TRUE){
  
  add_maf_cols <- function(variant.data, maf_cols){
    for (new_col in names(maf_cols)){
      if (new_col %!in% names(variant.data)){
        missing=TRUE
        for (old_col in maf_cols[[new_col]]){
          if (old_col %in% names(variant.data)){
            variant.data = variant.data %>% mutate("{new_col}" := variant.data[[old_col]])
            missing = FALSE
            break
          }
        }
        if(missing){
          stop(message(glue("Cannot create MAF column {new_col}")))
        }
      }
    }
    return(variant.data)
  }
  
  variant.data = standardize_names(variant.data, sid.format=sid.format, warn=F)
  if (!is.null(df_samples)){
    df_samples = standardize_names(df_samples, input.type="samples",sid.format=sid.format, warn=F)
    variant.data = variant.data %>% filter(SampleID.short %in% df_samples$SampleID.short)
    variant.data = merge.combine(join.type="inner", 
                                 join.cols.left="SampleID.short", 
                                 join.cols.right="SampleID.short",
                                 df1 = variant.data,
                                 df2 = df_samples %>% select(SampleID.short,PatientID,SampleType,StudyVisit) %>%
                                                      mutate(across(everything(), as.character)),
                                 priority="right",
                                 warn=warn)
  }
  
  ### for SNVs ####
  if (variant.type=="snv"){
    maf_cols = list("Tumor_Seq_Allele2"=list("alt","Tumor_Seq_Allele","Tumor_Seq_Allele1"), 
                 "Chromosome"=list("seqnames","chr"), 
                 "Start_Position"=list("start"),
                 "End_Position"=list("end"), 
                 "Reference_Allele"=list("ref"),
                 "Mutation_Status"= list("VariantType","Variant_Type"), #first replace
                 "Variant_Type"=list("VARIANT_CLASS","VariantType"), # second replace
                 "Variant_Classification"=list("Consequence"),
                 "Hugo_Symbol"=list("SYMBOL"),
                 "Tumor_Sample_UUID"=list("SampleID.short"))
    
    if ("finalKeep" %in% names(variant.data)){
      if (any(as.character(variant.data$finalKeep) %!in% c("TRUE","highConfidence","lowConfidence","FALSE"))){
        print(unique(as.character(variant.data$finalKeep)))
        print("WARN: some variants finalKeep term unrecognized. Keeping these")
        variant.data$finalKeep[is.na(variant.data$finalKeep) | (as.character(variant.data$finalKeep) %!in% c("TRUE","highConfidence","lowConfidence","FALSE"))] <- "TRUE"
      }
      variant.data = variant.data %>% filter(grepl("TRUE|highConfidence",as.character(finalKeep), ignore.case=T))
    }
    
    variant.data_maf = add_maf_cols(variant.data, maf_cols)
    
    if (any(c(is.na(variant.data_maf$Mutation_Status),
              is.null(variant.data_maf$Mutation_Status),
              variant.data_maf$Mutation_Status==""))){
      print("WARN: some varaints is.na(Variant_Type). Treating as somatic.")
      variant.data_maf$Mutation_Status[is.na(variant.data_maf$Mutation_Status)] <- "Somatic"
      variant.data_maf$Mutation_Status[is.null(variant.data_maf$Mutation_Status)] <- "Somatic"
      variant.data_maf$Mutation_Status[variant.data_maf$Mutation_Status==""] <- "Somatic"
    }
    
    print("SNV types present (selecting somatic):")
    print(unique(variant.data_maf$Mutation_Status))
    df_snvs = variant.data_maf %>% filter(grepl("somatic", Mutation_Status, ignore.case = T))
    warn_na(df_snvs[, names(maf_cols)])
    
    if (maf.output){
      
      if (!is.null(df_samples)){ 
        missing_samples = df_samples %>% filter(SampleID.short %!in% df_snvs$SampleID.short) %>% select(SampleID.short)
        df_snvs = plyr::rbind.fill(df_snvs, missing_samples)
      }else{
        print("No sample list provided, NOT including 0 variant samples in MAF output")
      }
      
      df_snvs = df_snvs %>% mutate("Variant_Classification" = replace_na(as.character(Variant_Classification), "N.A."),
                                   "Mutation_Status"="Somatic") %>%
                            unite(Tumor_Sample_Barcode, c(SampleID.short,PatientID,StudyVisit,SampleType)) %>%
                            select(all_of(c("Tumor_Sample_Barcode", names(maf_cols)))) %>%
                            mutate(across(everything(), as.character))
    }else{
      df_snvs = df_snvs %>% filter(!is.na(Hugo_Symbol)) %>% 
                            mutate("Mutation_Status"="Somatic") %>%
                            unite(Tumor_Sample_Barcode, c(SampleID.short,PatientID,StudyVisit,SampleType), remove=F)
      warn_na(df_snvs[ ,names(maf_cols)])
    }
    print(sprintf("Selected snvs: %s", nrow(df_snvs %>% filter(is.na(Variant_Classification) | Variant_Classification != "N.A."))))
    return(df_snvs)
    
### for CNVs ###
  }else if (variant.type=="cnv"){ ## format as maftools custom CNV table input: cnTable
    maf_cols=list("Hugo_Symbol"=list("Gene"), 
                  "Type"=list("Variant_Classification","CNV_Type")
                  )
    
    variant.data_maf = add_maf_cols(variant.data, maf_cols)
    # remove germ line CNVs
    for (germ in c("buffy coat", "pbmc")){
      if (germ %in% tolower(variant.data_maf$SampleType)){
        varaint.data_maf = variant.data_maf %>%
          group_by(PatientID, Hugo_Symbol, Type) %>% filter(any(SampleType==germ))
      }
    }
    df_cnv = variant.data_maf %>% filter(!is.na(Hugo_Symbol)) %>%
        unite(Tumor_Sample_Barcode, c(SampleID.short,PatientID,StudyVisit,SampleType), remove=FALSE) 
    if (maf.output){
      df_cnv = df_cnv %>% select(all_of(c("Hugo_Symbol","Tumor_Sample_Barcode", "Type"))) %>% # cols must be only these in THIS order
                          mutate(across(everything(), as.character))
      warn_na(df_cnv)
    }else{
      for (col in c("SampleType","StudyVisit","SampleID.short",names(maf_cols))){
        warn_na(df_cnv[[col]], col)
      }
    }
    print(sprintf("Selected CNVs: %s", nrow(df_cnv)))
    return(df_cnv)

  }else{
    stop("Undefined variant.type parameter for MAF conversion. (snv or cnv)")
  }
}