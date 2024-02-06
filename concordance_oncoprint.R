source("~/Desktop/puffin/R/format_as_MAF.R")
source("~/Desktop/puffin/R/helper_functions.R")
source("~/Desktop/puffin/R/get_concordance_stats.R")
library(maftools)
library(RColorBrewer)
library(ComplexHeatmap)


##' A function for reclassifying variant annotation strings, to be used with concordance_oncoprint().
##' @param annotations Consequence or Variant_Class annotations to recode
##' @param var.reduc.type string (from c("consequence_reduced", tbd . . .)) specifying a way to reclassify variant annotations, 
##' or a custom lsit of lists containing variant terms with the last element the replacement term.
##' @param cnv.strings optional terms to use for CNG and CNL. default: c("Gain","Loss") 
##' 
##' @return return variant string replacement
recode.variants <- function(annotations, var.reduc.set=NULL, cnv.strings=c("Gain","Loss")){
  annotations <- as.character(annotations)
  # TODO: MAF Standard
  #cbioportal convention, from: https://github.com/PoisonAlien/maftools/issues/686
  # vc    vc.cbio
  # 1:  Nonstop_Mutation Truncating
  # 2:   Frame_Shift_Del Truncating
  # 3: Missense_Mutation   Missense
  # 4: Nonsense_Mutation Truncating
  # 5:       Splice_Site Truncating
  # 6:   Frame_Shift_Ins Truncating
  # 7:      In_Frame_Del   In-frame
  # 8:      In_Frame_Ins   In-frame
  # var_color_dict = list(
  #     "Missense" = "darkgreen",
  #     "Promoter" = "plum1",
  #     "Truncating" = "blue",
  #     "In_Frame" = "darkred",
  #     "UTR_5p" = "cyan3",
  #     "UTR_3p" = "pink",
  #     "Gain" = "black",
  #     "Loss" = "grey52",
  #     "Multi_hit" = "black")
  
  if (is.null(var.reduc.set)){
    warning("No variant annotation reduction set specified.")
    return(annotations)
  } else if (is.list(var.reduc.set)){
    print("Using custom var.reduc.set.")
    var.reduc.code = var.reduc.set
  }else if (var.reduc.set=="consequence_reduced"){
    
    ## note, classes are in order of priority for annotation. "N.A." string needed so that samples with
    # zero variants are retained in the oncoprint.
    var.reduc.code = list(c('synonymous',"N.A."),
                          c('frameshift','Frame_Shift','stop_gained','stop_lost', 'start_gained',
                            'start_lost','splice','nonsense','Nonstop', "Truncating"),
                          c("5\'UTR","5_prime_UTR","UTR_5p"),
                          c('3_prime_UTR','3\'UTR',"UTR_3p"),
                          c("upstream","5\'Flank","Promoter"),
                          c('inframe','In_Frame',"In_Frame"),
                          c("missense","protein_altering_variant","Missense"),
                          c("downstream","retained",'IGR','Silent','Start_Codon_SNP',
                            'Translation_Start_Site','Unknown','Intron',"3\'Flank","RNA",
                            "non_coding_transcript_exon_variant","intergenic","coding_sequence_variant","N.A."),
                          c("gain","Amp","amp",cnv.strings[[1]],"Gain"),
                          c("loss","Del","deletion",cnv.strings[[2]],"Loss"))
  }else{
    print(var.reduc.set)
    stop("variant recoding not recognized. is it formatted as a list or a str? (var.reduc.set)")
  }
  # replace each variant function 
  recode_ <- function(vt, var.reduc.code){
    for (lst in var.reduc.code){
      for (term in lst[1:(length(lst)-1)]){
        if (grepl(term, vt, ignore.case=T)){
          vt.new = unlist(lst[[length(lst)]]) # substitute last term in list
          return(vt.new)
        }
        }
      }
    return(vt) # return the original variant code if not found in lists
  }
  
  variant_classes_new = unlist(lapply(annotations, function(x) recode_(x, var.reduc.code)))
  stopifnot(length(variant_classes_new) == length(annotations))
  
  #print(var.reduc.code)
  print("Original annotations:")
  print(table(annotations))
  print("New annotations:")
  print(table(variant_classes_new))
  
  return(variant_classes_new)
}


##' A function for reclassifying variant annotations, to be used with concordance_oncoprint(). See example concordance oncoprint in /man/figures/
##' @param variant.types optional vector of variant types to set colors for.
##' @param var.colors string (from c("consequence_reduced","MAF_standard", "random")) specifying a way to reclassify variant annotations, 
##' or a custom named list or vector matching variant type to color. default="random", sets random colors
##' 
##' @return return a vector of classifications same length as input vector
set_variant_colors <- function(variant.types=NA, var.colors=NULL, cnv.strings=c("Gain","Loss")){
  print("set_variant_colors . . .")
  if ((is.vector(var.colors) | is.list(var.colors)) & !is.null(names(var.colors))) {
      print("Using input custom variant colors.")
      print(var.colors)
      return(var.colors)
  }
  cnv1 = cnv.strings[[1]]
  cnv2 = cnv.strings[[2]]
  if(class(var.colors)=="character"){
  ## Important: dict keys listed in order of priority
  if(var.colors=="consequence_reduced"){
    print("Consequence reduced colors")
    var.colors = c(
      "Missense" = "darkgreen", 
      "Truncating" = "blue",
      "In_Frame" = "darkred",
      "UTR_5p" = "cyan3",
      "UTR_3p" = "maroon",
      "Promoter" = "pink")
    var.colors[cnv1] = "black"  #outline color
    var.colors[cnv2] = "grey52" #outline color
  }else if (var.colors=="MAF_standard"){
    print("MAF_standard colors")
    var.colors = c(
      "Missense_Mutation" = "darkgreen",
      "Nonsense_Mutation" = "red",
      "5pFlank" = "maroon",
      "5pUTR" = "cyan3",
      "Frame_Shift_Del" = "blue",
      "Splice_Site" = "orange",
      "Frame_Shift_Ins" = "purple",
      "In_Frame_Ins" = "darkred",
      "In_Frame_Del" = "yellow",
      #"Multi_hit" = "black",
      "Amplification" = "black", #outline color
      "Deletion" = "grey52") #outline color
  }
  }else if (is.null(var.colors)){
    #choose random colors
    print("var.colors not specified. choosing random colors.")
    if (is.na(variant.types)){
      stop("Must supply variant list for defining colors")
    }
    var.colors = get_random_color_dict(c(unique(variant.types),cnv.strings))
    ######
  }else{
    print(var.colors)
    stop("var.colors not recognized or wrong type")
  }
  var.colors[["Multi_hit"]] = "black"
  print(var.colors)
  "complete."
  return(var.colors) # a named list
}


##########
##' Helper function for defining multi_hits for concordance_oncoprint(). See example concordance oncoprint in /man/figures/
##' @param data a data frame containing all variants in one gene for one sample, annotated under col "Variant_Type"
##' @param cnv.strings an optional custom list of strings defining CNV Variant_Type in the input file
##' 
##' @return returns a list of variant types present in a particular gene.
call_variants <- function(data, cnv.strings=NA){
  var_types = list()
  df_snv = data %>% filter(!(Variant_Type %in% cnv.strings))
  df_cnv = data %>% filter(Variant_Type %in% cnv.strings)
  
  if (nrow(df_snv) > 1){
    var_types = c(var_types,"Multi_hit")}
  else if (nrow(df_snv) == 1){
    var = as.character(df_snv[1,"Variant_Classification"])
    var_types = c(var_types, var)}
  
  if (nrow(df_cnv) > 0){
    var = as.character(df_cnv[1,"Variant_Classification"])
    var_types = c(var_types, var)}
  stopifnot(length(var_types)>0)
  return(var_types)
}


##' Function for formatting variant data for input to concordance_oncocoprint(). See example concordance oncoprint in /man/figures/
##' @param variant_data a data frame containing at all variants for oncoprint. Req'd cols: Hugo_Symbol, PatientID, Variant_Classification
##' @param sample_types a list of two "StudyVisit_SampleType" strings defining which paired samples to include in the oncoprint: e.g. "Postop1_Tissue"))
##' @param df_samples a dataframe containing SampleID.short and PatientID and clinical annotations
##' @param cnv.strings default:c("Gain","Loss")
##' @return returns a data frame containing oncoplot annotaitons, rows=genes, cols=patients.
##' 
format_oncoprint_data <- function (variant_data, df_samples, sample_types=NA, cnv.strings=c("Gain","Loss")){
  print("format_oncoprint_data . . .")
  reqd_cols = c("Hugo_Symbol", "PatientID", "Variant_Classification", "SampleID.short", "SampleType", "StudyVisit")
  if (any(reqd_cols %!in% names(variant_data))){
    print("Missing data input cols in variant data:") 
    print(reqd_cols[!(reqd_cols %in% names(variant_data))])
    stop()
  }
  if (is.na(cnv.strings[1]) |is.na(cnv.strings[2]) ){
    stop("Must define cnv.strings")
  }
  print(glue("Input paired sample types: {sample_types}"))
  variant_data = variant_data %>% filter(SampleID.short %in% df_samples$SampleID.short) %>%
    mutate("sample_type"=paste0(StudyVisit,"_",SampleType)) %>%
    filter(Variant_Classification != "N.A.", Hugo_Symbol != "UnknownGene")  #TODO: move this filter earlier in pipe
  check.missing(sample_types, variant_data$sample_type)
  variant_data = variant_data %>% filter(sample_type %in% sample_types)
  
  genes = unique(variant_data$Hugo_Symbol)
  patients = unique(df_samples$PatientID)
  plot_data = c()
  for (patient in patients){
    print(patient)
    df_patient = variant_data %>% filter(PatientID==patient) 
    if (length(unique(df_patient$SampleID.short))>2){
      print(unique(df_patient$SampleID.short))
      #print(df_patient, n=1000)
      stop(sprintf("WARN: Too many samples: patient-%s", patient))
    }
    patient_data_by_gene = vector()
    if (nrow(df_patient) == 0){
      patient_data_by_gene = rep("",length(genes))
    }else{
      for (gene in genes){
        #print(gene)
        df_gene = df_patient %>% filter(Hugo_Symbol==gene)
        gene_annot = list()
        if (nrow(df_gene) > 0){
          for (st in sample_types){
            if (st %in% df_gene$sample_type){
              type_annots = call_variants(df_gene %>% filter(sample_type==st), cnv.strings)
              if (length(type_annots) > 0) { 
                gene_annot = c(gene_annot, paste0(paste0(st,"."), type_annots, collapse=";"))
              }}}
        }else{
          gene_annot = c("")
        }
        gene_annot = paste(gene_annot, collapse=";")
        if (grepl("NULL", gene_annot)){
          print(df_gene)
          stop("fix NULL variant annotation")
        }
        #print(gene_annot)
        patient_data_by_gene = c(patient_data_by_gene, gene_annot)
      } #end genes
    } # end else
    stopifnot(length(patient_data_by_gene)==length(genes))
    plot_data[[patient]] <- as.vector(patient_data_by_gene)
  }#end patients
  
  df_plot_data = as.data.frame(do.call(cbind, plot_data))
  rownames(df_plot_data) = genes
  stopifnot(length(names(df_plot_data)) == length(unique(df_samples$PatientID)))
  #print(grep("NULL", df_plot_data))
  print("complete.")
  return(df_plot_data)
}


## a function for defining shapes in an oncoprint
##' @param variant_classes a list of variant classifications to use in plot
##' @param l plot grid line width
##' @param var.colors a named List with entries "variant annotation"=color
##' @param background_col optional background color string. default: "snow2" 
##' @param ref_sample_type an underscore delim. string indicating the sample type and study visit to use as reference: eg "Preop_Tissue"
##' @param mrd_sample_type an underscore delim. string indicating the MRD sample type and study visit to use : eg "Postop_Urine"
##' @param cnv.strings an optional list of CNV variant annotations. default: c("Gain","Loss")
##' @param concord.barplot logical value whether to show concordance barplot colors in the legend
##' 
##' @return the completed alter_fun input to Oncoprint function
create_alter_fun <- function(variant_classes, l = 0.05, var.colors="random", background_col="snow2",
                             ref_sample_type = NA, mrd_sample_type = NA, cnv.strings=c("Gain","Loss"),
                             make.legend=FALSE, concord.barplot=FALSE){
  print("create_alter_fun . . .")
  stopifnot(all(!(is.na(var.colors))))
  ref_sample_type = tolower(ref_sample_type)
  mrd_sample_type = tolower(mrd_sample_type)
  variant_types = c(unique(variant_classes), "Multi_hit")
  cnv1 = cnv.strings[1]
  cnv2 = cnv.strings[2]
  half_fill_w = (1-l)/2
  #### define an oncoprint function list with variant shapes and colors
  alter_fun = list()
  alter_fun[["background"]] = local({
    half_fill_w = half_fill_w
    background_col = background_col
    function(x, y, w, h) {
      grid.rect(x, y, w*half_fill_w*2, h*half_fill_w*2, gp = gpar(fill = background_col, col = NA))
    }
  })
  
  # make upper triangles for mrd sample type
  for (vt_str in names(var.colors)){
    name_str = paste0(mrd_sample_type, ".", vt_str)
    color = var.colors[[vt_str]]
    if (!(vt_str %in% cnv.strings)){
      alter_fun[[name_str]] = local({
        half_fill_w = half_fill_w
        color = color
        function(x, y, w, h) {
          grid.polygon( x = c(x - half_fill_w*w, x + half_fill_w*w, x + half_fill_w*w), 
                        y = c(y + half_fill_w*h, y + half_fill_w*h, y - half_fill_w*h),
                        gp = gpar(fill = color, col = "white"))
        }
      })
    }else if (vt_str %in% cnv.strings){
      alter_fun[[name_str]] = local({
        half_fill_w = half_fill_w
        color = color
        function(x, y, w, h) {
          grid.polygon( x = c(x - half_fill_w*w, x + half_fill_w*w, x + half_fill_w*w), 
                        y = c(y + half_fill_w*h, y + half_fill_w*h, y - half_fill_w*h),
                        gp = gpar(fill = NA, col = color, lwd=2))
        }
      })
    }
  }
  # make lower triangles for ref (baseline) sample type
  for (vt_str in names(var.colors)){
    name_str = paste0(ref_sample_type, ".", vt_str)
    color = var.colors[[vt_str]]
    if (!(vt_str %in% cnv.strings)){
      alter_fun[[name_str]] = local({
        half_fill_w = half_fill_w
        color = color
        function(x, y, w, h) {
          grid.polygon( x = c(x - half_fill_w*w, x + half_fill_w*w, x - half_fill_w*w), 
                        y = c(y - half_fill_w*h, y - half_fill_w*h, y + half_fill_w*h),
                        gp = gpar(fill = color, col = "white"))
        }
      })
    }else if (vt_str %in% cnv.strings){
      alter_fun[[name_str]] = local({
        half_fill_w = half_fill_w
        color = color
        function(x, y, w, h) { 
          grid.polygon( x = c(x - half_fill_w*w, x + half_fill_w*w, x - half_fill_w*w), 
                        y = c(y - half_fill_w*h, y - half_fill_w*h, y + half_fill_w*h),
                        gp = gpar(fill = NA, col = color, lwd=2))
        }
      })
    }
  }
  print("complete.")
  return(alter_fun)
}

##### make an oncoprint legend
##' @param variant_classes a list of variant classifications to use in legend
##' @param l plot grid line width
##' @param var.colors a named List with entries "variant annotation"=color
##' @param background_col optional background color string. default: "snow2" 
##' @param ref_sample_type an underscore delim. string indicating the sample type and study visit to use as reference: eg "Preop_Tissue"
##' @param mrd_sample_type an underscore delim. string indicating the MRD sample type and study visit to use : eg "Postop_Urine"
##' @param cnv.strings an optional list of CNV variant annotations. default: c("Gain","Loss")
##' @param concord.barplot logical value whether to show concordance barplot colors in the legend
##' @param l plot line width
##' 
##' @return Oncoprint legend object
make_oncoprint_legend <- function(variant_classes, var.colors="random", background_col="snow2",
                                  ref_sample_type = NA, mrd_sample_type = NA, cnv.strings=c("Gain","Loss"),
                                  concord.barplot=FALSE, l=0.05, plot.file=NA){
  print("make_oncoprint_legend . . .")
  half_fill_w = (1-l)/2
  manual_lgds = NA
  graphics = list()
  variant_classes = c(unique(variant_classes), "Multi_hit")
  stopifnot(all(variant_classes %in% names(var.colors)))
  # add SNVs
  for (vt in variant_classes){
    if (vt %!in% cnv.strings){
      func <- local({
        color = var.colors[[vt]]
        function(x, y, w, h) grid.rect(x, y, w, h, gp=gpar(fill = color, col="white"))
      })
      graphics[[vt]] = func
    }
  }
  # add CNVs
  graphics[[cnv.strings[1]]] <- local({
    color = var.colors[cnv.strings[1]]
    function(x, y, w, h) grid.rect(x, y, w, h, gp=gpar(fill = background_col, col=color, lwd=3))
  })
  graphics[[cnv.strings[2]]] <- local({
    color = var.colors[cnv.strings[2]]
    function(x, y, w, h) grid.rect(x, y, w, h, gp=gpar(fill = background_col, col=color, lwd=3))
  })
  
  variant_lgd = Legend(labels = names(graphics),
                       graphics = graphics)
  sample_type_lgd = Legend(labels = c(mrd_sample_type, ref_sample_type),
                           graphics = list(
                             function(x, y, w, h) grid.polygon(x = unit.c(x - half_fill_w*w, x + half_fill_w*w, x + half_fill_w*w), 
                                                               y = unit.c(y + half_fill_w*h, y + half_fill_w*h, y - half_fill_w*h),
                                                               gp=gpar(fill = background_col, col=background_col, lwd=1)),
                             function(x, y, w, h) grid.polygon(x = unit.c(x - half_fill_w*w, x + half_fill_w*w, x - half_fill_w*w), 
                                                               y = unit.c(y - half_fill_w*h, y - half_fill_w*h, y + half_fill_w*h),
                                                               gp=gpar(fill = background_col, col=background_col, lwd=1))))
  manual_lgds = packLegend(sample_type_lgd, variant_lgd, direction = "vertical")
  bar_lgd = NA
  if (concord.barplot){
    bar_lgd = Legend(labels = c("concordant","non-concordant"),
                     graphics = list( 
                       function(x, y, w, h) grid.rect(x, y, w, h, gp=gpar(fill = "grey39", col="grey39")), #dark grey - concordant
                       function(x, y, w, h) grid.rect(x, y, w, h, gp=gpar(fill = "grey", col="grey")) #light grey - unique
                     ) 
    )
    manual_lgds = packLegend(bar_lgd, sample_type_lgd, variant_lgd, direction = "vertical")
  }
  
  jpeg(paste0(plot.file,".jpg"), width=3, height=5, units="in", res=300)
  draw(manual_lgds)
  dev.off()
  print("complete.")
  return(manual_lgds)
}

########
##' Main plotting function
##' Wrapper around ComplexHeatmap functions to create paired sample oncoprint. See example in /man/figures/concordance_oncoprint_example.jpg
##' @param snv.data required data frame containing SNV/Indels listed by SampleID
##' @param cnv.data required data frame containing CNVs listed by SampleID
##' @param clin.data required data frame containing clinical data for each PatientID
##' @param df_samples required data frame with sampleIDs to include. req'd cols: SampleID,PatientID,StudyVisit,SampleType
##' @param patients optional list of PatientIDs to include
##' @param cnv.strings variant annotations for CNGs and CNLs. default: c("Gain","Loss")
##' @param clin.data.cols optional vector of clinical data columns to annotate
##' @param ref_sample_type an underscore delim. string indicating the sample type and study visit to use: eg "Preop_Tissue"
##' @param mrd_sample_type an underscore delim. string indicating the sample type and study visit to use: eg "Preop_Tissue"
##' @param concord.barplot logical value whether to show a concordance barplot at the top of the oncoprint. default: TRUE
##' @param min.samples.mutated optional, only show genes with at least this many samples mutated
##' @param min.patients.mutated optional, only show genes with at least this many patients mutated. 
##' @param genes optional list of genes to show in the plot.
##' @param show.top.n optional number of top genes to show.
##' @param show.clin.data logical value whether to show clinical annotations at the bottom of the plot. default:c("PatientID")
##' @param show.patient.id logical value whether to annotate PatientIDs on the plot. default: TRUE
##' @param clin.annotation.colors optional list of clinical column labels in dot format, each matched with a named vector of levels to colors (e.g. list(Stage=c("T0"="blue","T1"="red")))
##' @param make.legend logical value whether to generate a custom legend for the variants in the plot. default: FLASE
##' @param var.reduc.set string ("consequence_reduced","MAF_reduced","MAF_standard","random") or a custom named 
##' List with entries "variant class"="new variant class",**listed in order of importance. default: "random"
##' @param var.colors optional custom dictionary matching variant types with oncoplot colors or string indicating a predefined set of colors
##' otherwise randomly selected colors.
##' @param alter_fun optinal external alter_fun
##' @param out.file.name
##' @param . . . additional parameters for Oncoprint plotting function: https://jokergoo.github.io/ComplexHeatmap/reference/oncoPrint.html
##' 
##' @return named list with "oncoprint"=graphing heatmap object and "legend"=legend object.
##' \if{html}{\figure{"/man/figures/concordance_oncoprint_example.jpg"}{options: width=100 alt="R logo"}}
##' 
concordance_oncoprint <- function(snv.data=NA, cnv.data=NULL, clin.data=NA, df_samples=NA, patients=NA, 
                                  clin.data.cols=c("PatientID"),
                                  ref_sample_type=NA, mrd_sample_type=NA, 
                                  var.dict=NA, cnv.strings=c("Gain","Loss"),
                                  min.samples.mutated=NA, min.patients.mutated=NA, sid.format="none",
                                  show.clin.data=FALSE, show.patient.id=TRUE, clin.annotation.colors=list(),
                                  var.reduc.set=NULL, var.colors=NULL, make.legend=FALSE, concord.barplot=TRUE,
                                  show.top.n=NA, genes=NULL, alter_fun=NULL, out.file.name=NULL){

  print("Setting params and formatting data . . .")
  ref_sample_type = tolower(ref_sample_type)
  mrd_sample_type = tolower(mrd_sample_type)
  df_samples = standardize_names(df_samples, sid.format=sid.format, input.type="samples")
  snv.data = standardize_names(snv.data, sid.format=sid.format, warn = FALSE)
  if (is.null(cnv.data)){
    cnv.data = data.frame(matrix("NA", nrow = 1, ncol = 8))
    names(cnv.data) <- c("SampleID.short","SampleType","PatientID","StudyVisit","Type",
                         "Hugo_Symbol","Variant_Classification","VariantID")
  }else{
    cnv.data = standardize_names(cnv.data, warn = FALSE)
  }
  
  # filter input data by Patient
  stopifnot("PatientID" %in% names(clin.data))
  clin.data = clin.data %>% select(unique(c("PatientID", clin.data.cols))) %>%
    filter(PatientID %in% df_samples$PatientID)
  stopifnot("PatientID" %in% names(clin.data))
  stopifnot(all(!(duplicated(clin.data$PatientID))))
  if (!(is.na(patients))){
    df_samples = df_samples %>% filter(PatientID %in% patients)
  }else{
    #patients = unique(df_samples$PatientID)
    patients=clin.data$PatientID[clin.data$PatientID %in% unique(df_samples$PatientID)]
  }
  # filter input data by sampleID
  df_samples = df_samples %>% select(PatientID, StudyVisit, SampleType, SampleID.short) %>%
    mutate(across(c("SampleType","StudyVisit"), tolower)) %>%
    unite(Tumor_Sample_Barcode, c(SampleID.short, PatientID, StudyVisit, SampleType), remove=FALSE) %>% 
    mutate("order"=match(PatientID, clin.data$PatientID)) %>% arrange(order)
  stopifnot(all(!(duplicated(df_samples$Tumor_Sample_Barcode))))

  print("selecting snvs . . .")
  all.snv_selected = snv.data %>% filter(SampleID.short %in% df_samples$SampleID.short,
                                         !grepl("synon", Variant_Classification, ignore.case=T),
                                         !is.na(Variant_Classification)) %>%
    mutate(Variant_Classification=recode.variants(Variant_Classification, var.reduc.set, cnv.strings)) %>%
    select(-PatientID)
  all.snv_selected = merge.combine(all.snv_selected,
                                   df_samples %>% select(SampleID.short, SampleType, StudyVisit, PatientID))

  print("selecting cnvs . . .")
  all.cnv_selected = cnv.data %>% filter(SampleID.short %in% df_samples$SampleID.short, !is.na(Hugo_Symbol)) %>% 
    mutate(Variant_Classification=recode.variants(Variant_Classification, var.reduc.set, cnv.strings))
  all.cnv_selected = merge.combine(all.cnv_selected,
                                   df_samples %>% select(SampleID.short, SampleType, StudyVisit, PatientID),
                                   join.type="left", join.cols.left = "SampleID.short", join.cols.right = "SampleID.short", priority = "right")
  print(nrow(all.cnv_selected))
  # set variant colors
  print("setting colors and symbols . . .")
  variant.classes = unlist(unique(c(all.cnv_selected$Variant_Classification, all.snv_selected$Variant_Classification)))
  variant.classes = variant.classes[variant.classes != "N.A."]
  if (is.null(var.colors)){
    var.colors = set_variant_colors(variant.classes, cnv.strings=cnv.strings, var.colors=var.reduc.set)
  }
  
  if (any(variant.classes %!in% names(var.colors))){
    print("WARN: these variant classes missing from var.colors.")
    print(data.frame("class"=(variant.classes), "present"=(variant.classes) %in% names(var.colors)))
  }
  # TODO: create an option to input custom variant reduction types
  # all.snv_selected_maf$Variant_Classification = factor(all.snv_selected_maf$Variant_Classification)
  # levels(all.snv_selected_maf$Variant_Classification) <- var.dict
  if (is.null(alter_fun)){
    alter_fun <- create_alter_fun(variant_classes = variant.classes,
                                  ref_sample_type = ref_sample_type,
                                  mrd_sample_type = mrd_sample_type,
                                  var.colors = var.colors,
                                  cnv.strings = cnv.strings) 
  }
  
  # optionally use colors to make a legend
  if (is.null(out.file.name)) out.file.name <- glue("./{ref_sample_type}_vs_{mrd_sample_type}_oncoprint_{Sys.Date()}")
  plot.legend=NA
  if (make.legend){
    plot.legend = make_oncoprint_legend(variant_classes = variant.classes[],
                                        ref_sample_type = ref_sample_type,
                                        mrd_sample_type = mrd_sample_type,
                                        var.colors = var.colors,
                                        concord.barplot = concord.barplot, #TRUE/FALSE
                                        cnv.strings = cnv.strings,
                                        plot.file=paste0(out.file.name,"_legend"))
    draw(plot.legend)
  }
  
  # use MAFtools to check and filter MAF data. maf format adds N.A. variants for samples w/no variants
  # TODO, maybe make this step optional
  print("Filtering MAF input . . .")
  all.snv_selected_maf = format_as_MAF(all.snv_selected, df_samples = df_samples, variant.type="snv", sid.format = sid.format)
  all.cnv_selected_maf = format_as_MAF(all.cnv_selected, df_samples = df_samples, variant.type="cnv", sid.format = sid.format)
  print(nrow(all.cnv_selected_maf))
  vc_nonsyn = unique(all.snv_selected_maf$Variant_Classification)
  vc_nonsyn = vc_nonsyn[!is.na(vc_nonsyn) & vc_nonsyn!="N.A."]
  dat <- read.maf(all.snv_selected_maf, 
                  cnTable = all.cnv_selected_maf, 
                  vc_nonSyn = vc_nonsyn, 
                  clinicalData = df_samples)
  # reformat
  df_data = dat@data %>% filter(!is.na(Hugo_Symbol)) %>% left_join(df_samples, by="Tumor_Sample_Barcode")
  df_data$Variant_Type <- as.character(df_data$Variant_Type)
  df_data$Variant_Classification <- as.character(df_data$Variant_Classification)
  df_data$Variant_Classification <- sapply(df_data$Variant_Classification, function(x) str_replace(x, "\'","p"))
  ## generate plot data
  df_plot_data = format_oncoprint_data(variant_data = df_data,
                                       df_samples = df_samples,
                                       sample_types=c(ref_sample_type, mrd_sample_type),
                                       cnv.strings=cnv.strings)
  # explicitly order samples by patient
  df_plot_data = df_plot_data %>% relocate(patients)
  
  # make plot clinical annotations - TODO: test this
  print("formatting clinical data . . . ")
  df_annot = as.data.frame(clin.data)
  rownames(df_annot) <- seq(1,nrow(df_annot),1)
  if (!(show.clin.data)){
    df_annot = df_annot %>% select(PatientID)
  }else{
    if (!(show.patient.id)){
      df_annot=df_annot %>% select(-PatientID)
    }
  }
  #names(df_annot) <- make.names(names(df_annot))
  print(df_annot)
  stopifnot(is.list(clin.annotation.colors))
  clinical_annot = HeatmapAnnotation(df=df_annot, 
                                     col=clin.annotation.colors, # factors not specified will get random colors
                                     annotation_name_side = "left")
  print(clinical_annot)
  
  # optionally get concordance stats and make barplot annotation
  concord_muts_annot = HeatmapAnnotation(foo = anno_empty(border = FALSE))
  if (concord.barplot){
    print("concordance_barplot . . .")
    cols = intersect(names(all.snv_selected),names(all.cnv_selected))
    stopifnot(all(c("SampleID.short","VariantID") %in% cols))
    all.vars = rbind(all.snv_selected[,cols], all.cnv_selected[,cols]) %>% filter(!is.na(VariantID))
    df_concord_counts = get_concordance_stats(all.vars, 
                                              ref_sample_type = ref_sample_type,
                                              df_samples = df_samples)
    df_barplot = df_concord_counts %>% filter(sample_type==mrd_sample_type) %>% left_join(clin.data, by="PatientID") %>%
      mutate("order"=match(PatientID, clin.data$PatientID)) %>% arrange(order) %>% 
      select(n_concord_vars, n_unique_vars) %>% mutate(across(everything(), as.numeric))
    if (nrow(df_barplot) != ncol(df_plot_data)){
      print(sprintf("N Concord. barplot samples: %s", nrow(df_barplot)))
      print(sprintf("N variant data samples: %s", ncol(df_plot_data)))
      stop("Some samples not paired across input types.")
    }
    concord_muts_annot = HeatmapAnnotation(N.MRD.variants = anno_barplot(df_barplot, 
                                                                         axis = TRUE, 
                                                                         # axis_param = list(gp=gpar(fontsize=10),
                                                                         #                           at = c(25,50,75,100,300), 
                                                                         #                           labels = c(25,50,75,100,300)), 
                                                                         height = unit(7, "cm")))
    print(concord_muts_annot)
  }
  
  # optionally filter the oncoprint genes
  print("Selecting data for plot . . .")
  # count_samples <- function(cell){ ## TODO filter by n.samples muated
  #   if (cell==""){
  #     return(0)
  #   }else{
  #     vars = unlist(str_split(cell, ";"))
  #     samps = unique(sapply(str_split(vars, "."), "[",1))
  #     return(length(samps))}
  # }
  # order genes by n.patients mutated
  df_plot_data = as.data.frame(df_plot_data) %>% filter(!is.na(row.names(df_plot_data))) %>%
                    mutate("n.patients"= apply(df_plot_data, 1, function(x) sum(x != ""))) %>%
                    arrange(desc(n.patients)) %>% select(-n.patients)
  if (!(is.na(show.top.n))){
    df_plot_data = df_plot_data[1:min(show.top.n,nrow(df_plot_data)),]
  }
  if (!is.null(genes)){
    print("Using custom gene list.")
    #add rows for missing genes
    missing.genes = genes[genes %!in% row.names(df_plot_data)]
    df_missing = df_plot_data[1:length(missing.genes),]
    df_missing[,] <- ""
    rownames(df_missing) <- missing.genes
    names(df_missing) <- names(df_plot_data)
    df_plot_data = rbind(df_plot_data, df_missing)
    df_plot_data = df_plot_data[match(genes, rownames(df_plot_data)),]
  }
  df_plot_data = df_plot_data %>% select(clin.data$PatientID) # order the data by patient order in clin.data
  #print(head(df_plot_data))
  #print(clin.data)
  #stop()
  df_plot_data = as.matrix(df_plot_data)
  
  ### Main plotting function ####
  print("Constructing plot . . .")

  if (is.null(out.file.name)) out.file.name <- glue("./{ref_sample_type}_vs_{mrd_sample_type}_oncoprint_{Sys.Date()}")
  jpeg(paste0(out.file.name,".jpg"), units="in", height=10, width=nrow(clin.data)*0.25 + 3, res=300)

  oncoprint <- oncoPrint(df_plot_data, alter_fun = alter_fun,
                         remove_empty_rows = FALSE,
                         row_order = c(1:nrow(df_plot_data)),
                         col = c(), 
                         bottom_annotation = clinical_annot,
                         top_annotation = concord_muts_annot,
                         pct_side = "right", row_names_side = "left",
                         column_order = c(1:length(patients)),
                         show_heatmap_legend=FALSE)
  draw(oncoprint)
  dev.off()
  print("plot written to file: {out.file.name}.jpg")
  #draw(oncoprint) doesn't work, only pdf output works currently
  #TODO:add kwargs
  # plot_w_kwargs <- function(func, ...) {
  #   func(...)
  # }
  #plot.out <- plot_w_kwargs(oncoprint_custom, ...)
  return(c("oncoprint"=oncoprint, "legend"=plot.legend))
}
