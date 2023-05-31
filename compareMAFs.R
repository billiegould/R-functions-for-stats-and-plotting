library(tidyverse)
library(glue)
source("~/Desktop/puffin/R/helper_functions.R")

## for an input dataframe with sample variants, calculate and plot MAFs for concordant and non-concordant variants for each patient
##' @param snvs.selected a data frame of SNV variants containing PatientID, SampleIDs, StudyVisit, SampleType and VariantID
##' @param samples required data frame of PatientID paired SammpleIDs to plot
##' @param selectors a list of two sample types to use for each patient, each a string in the form: StudyVisit_SampleType, eg "Tissue_Preop".
##' @param clin.data optional dataframe containing PatientIDs and clinical info 
##' @param plot.file optional path to file for writing the concordance output.
##' @param color.by optional clinical data value to color variants. default: PatientID
##' @param colors optional named list matching color.by category to colors. default: chooses random colors
##' @param log.axes whether to use log2 MAF axes. default: FALSE
##' @param plot whether to generate a plot. default=TRUE
##' @param plot.title title if generating a plot. default=<blank>
##' @param legend.position legend position if generating a plot. can be "none". default="right"
##' @param plot.file optional file path in which to save the plot.
##' @param fill.na optional MAF to use for non-concordant variants default: 0.0%
##' @param color.by.foreground.value optional variant color.by string group to plot in the foreground
##' @param pt.sz optional point size. defaul=10
##' 
##' @return a named list with: variants = dataframe of variant frequencies, plot =  ggplot object
compareMAFs <- function(all.snvs=NA, samples=NULL, selectors=NULL, clin.data=NA, fill.na=0.0, text.size=12,
                        color.by="PatientID", log.axes=FALSE, plot=TRUE, plot.title="", colors=NULL,
                        legend.position="right", plot.file=NA, pt.sz=10, color.by.foreground.value=NA){
  # WARN: patched variants are often SampleType and StudyVisit==NA
  if (!is.null(samples)){
    counts = samples %>% group_by(PatientID) %>% summarize(counts=n())
    if (any(counts$counts != 2)){
      print(counts)
      print("WARN: some Patients with >2 or < 2 samples")
    }
  }else{
    stop("Must supply sample data frame.")
  }
  all.snvs = standardize_names(all.snvs)
  samples = standardize_names(samples, input.type="samples")
  if (!(is.na(clin.data))){
    samples = samples %>% left_join(clin.data, by="PatientID")
  }
  snvs.selected = all.snvs %>% filter(SampleID.short %in% samples$SampleID.short) %>%
    select(-SampleType, -StudyVisit, -PatientID) %>%
    left_join(samples, by="SampleID.short") %>%
    mutate(Selector=paste0(StudyVisit,"_",SampleType)) 
  color.by_ = snvs.selected[,color.by]
  snvs.selected = snvs.selected %>% mutate("Color"=color.by_) %>%
    select(PatientID, SampleID.short, StudyVisit, SampleType, Hugo_Symbol, 
           VariantID, VariantFreq, VariantType, Color, Selector)
  if (!(is.null(selectors))){
    print(sprintf("input selectors: %s", selectors))
    df1 = snvs.selected %>% filter(Selector==selectors[[1]])         
    df2 = snvs.selected %>% filter(Selector==selectors[[2]])
  }else{
    stop("Sample types and study visits not specified (selectors=NA)")
  }
  print(nrow(df1))
  print(nrow(df2))
  if (nrow(df1)==0 | nrow(df2)==0){
    print("WARN: No variants selected for one or both sample types.")
    print(glue("Selectors in the data: {unique(snvs.selected$Selector)}"))
  }
  if (log.axes & fill.na==0){
    fill.na = round(0.5 * min(snvs.selected$VariantFreq), 2)
  }
  df_compare_freq = df1 %>% full_join(df2, by=c("PatientID","VariantID")) 
  df_compare_freq = df_compare_freq %>% mutate(Concordant=factor(!is.na(VariantFreq.x) & !is.na(VariantFreq.y), 
                                                                 levels=c(TRUE,FALSE)),
                                               Color=merge_info(Color.x, Color.y)) %>%
                            replace_na(list(VariantFreq.x = fill.na, VariantFreq.y=fill.na))
  n.variants = nrow(df_compare_freq) + sum(df_compare_freq$Concordant=="TRUE")
  stopifnot(n.variants == nrow(snvs.selected))
  # if (!is.na(clin.data))){
  #   df_compare_freq = df_compare_freq %>% left_join(clin.data, by="PatientID")
  # }
  gg = "meep"
  if (plot){
    if (!(is.na(plot.file))){
      pdf(plot.file, height=10, width=12)
    }
    if (is.null(colors)){
      print(unique(df_compare_freq$Color))
      colors = get_random_color_dict(df_compare_freq$Color)
    }
    if (any(df_compare_freq$Color %!in% names(colors))){
      print(colors)
      check.missing(list.ref=names(colors), list.test=df_compare_freq$Color)
      stop("Missing plot point colors.")
    }
    upper.lim = max(c(df_compare_freq$VariantFreq.x,df_compare_freq$VariantFreq.y), na.rm=T)+3
    gg <- ggplot(df_compare_freq, aes(x=VariantFreq.x, y=VariantFreq.y, color=Color)) +
      geom_point(aes(
        size=pt.sz, 
        alpha=0.8, ## something weird here
        shape=Concordant
      )) +
      scale_shape_manual(values = c("TRUE"=19,"FALSE"=17)) +
      scale_color_manual(name="PatientID", values = colors) +
      geom_abline(slope=1, intercept=c(0,0), linetype="dashed") +
      ylim(fill.na, upper.lim) + # set axes equal
      xlim(fill.na, upper.lim) +
      labs(x=paste(selectors[1], "MAF (%)"), y=paste(selectors[2], "MAF (%)"), title=plot.title) +
      theme(text=element_text(size=text.size), legend.position=legend.position, plot.title=element_text(size=20))
    if (!is.na(color.by.foreground.value)){
      df = df_compare_freq %>% filter(Color==color.by.foreground.value)
      if (nrow(df)!=0){
        gg <- gg +
          geom_point(df, mapping=aes(size=pt.sz, alpha=0.8, shape=Concordant))
      }
    }
    if (log.axes){
      gg <- gg +
        scale_y_continuous(trans='log2', breaks=c(fill.na,0.1,0.5, 1.0, 5.0, 10.0, 25.0, 50.0,100), limits=c(fill.na, upper.lim)) +
        scale_x_continuous(trans='log2', breaks=c(fill.na,0.1,0.5, 1.0, 5.0, 10.0, 25.0, 50.0,100), limits=c(fill.na, upper.lim))
    }
    if (all(df_compare_freq$VariantFreq.x == fill.na)){
      gg <- gg + xlim(fill.na, max(df_compare_freq$VariantFreq.y))
    }
    if (all(df_compare_freq$VariantFreq.y == fill.na)){
      gg <- gg + ylim(fill.na, max(df_compare_freq$VariantFreq.x))
    }
    show(gg)
    if (!(is.na(plot.file))){
      dev.off()
    }
  }
  return(list("variants"=df_compare_freq %>% select(PatientID, SampleID.short.x, SampleID.short.y,
                                                    Selector.x, Selector.y, VariantFreq.x, VariantFreq.y,
                                                    VariantType.x, VariantType.y,
                                                    Concordant, Color) %>%
                                              rename(color_by = Color), 
              "plot"=gg,
              "colors"=colors))
}