### helper functions for post-analysis

library(tidyverse)
library(glue)
library(RColorBrewer)
library(stats)

options(stringsAsFactors = FALSE)

## useful
`%!in%` = Negate(`%in%`)

## also useful
make_names <- function(df, verbose=FALSE){
  names(df) <- trimws(make.names(names(df)))
  if (verbose){
    print(paste(unlist(names(df)), collapse=", "))
  }
  return(df)
}

## reduce sample IDs to base id number, remove aliquot and library identifiers
##' @param col A dataframe column of SampleIDs
##' @param check whether to print out the id conversion for manual checking. default: FALSE
##' @param sid.format SampleID.short style c("none","strict","remove.suffix","remove.A01")
##' 
##' @return a column of shortened SampleIDs 
make_SIDshort <- function(col, sid.format="none"){
  stopifnot(is.vector(col))
  print(glue("SampleID.short format: {sid.format}"))
  if (sid.format=="strict"){
    return(str_sub(col, 1,7))
  }else if (sid.format=="none"){
    return(col)
  }else if (sid.format=="remove.suffix"){
    return(gsub("A01", "", col))
  }else if (sid.format=="remove.suffix"){
    retunr(sapply(str_split(out, "_"),"[",1))
  }else{
    stop(glue("Unknown SampleID.short format: {sid.format}"))
  }
}

## print warning if any value is.na
##'
##'
warn_na <- function(my.vector, my.vector.name=NA){
  
  warn_ <- function(vec, name){
    if (any(is.na(vec) | is.null(vec))){
      tot = sum(is.na(vec) | is.null(vec) | vec=="")
      message(glue("warn_na: WARN: NA or missing values in {name}: {tot}"))
    }
  }
  if (is.na(my.vector.name)){
    my.vector.name = deparse(substitute(my.vector))
  }
  if (is.data.frame(my.vector)){
    for (col in names(my.vector)){
      warn_(my.vector[[col]], col)
    }
  }else{
    warn_(my.vector, my.vector.name)
  }
}

### fill in col NAs with 0s at indexes where values in another column are 0
fill_na_0 <- function(col, index_col){
  col[index_col==0] <- 0.0
  return(col)
}

## bind two data frames using only the cols they have in common
rbind.common <- function(df1, df2){
  cols = intersect(names(df1), names(df2))
  if (length(cols)==0){
    stop("no columns in common")
  }else{
    return(rbind(df1[,cols], df2[,cols]))
  }
}

## summary stats for a vector
##' @param vec a numeric vector or dataframe column
##'
summarize_vector <- function(vec){
  name = deparse(substitute(vec))
  print(name)
  if (is.numeric(vec)){
    print(sprintf("n.obsv.: %s", length(vec[!is.na(vec)])))
    warn_na(vec, my.vector.name = name)
    mean.=mean(vec, na.rm=T)
    median.=median(vec, na.rm=T)
    min.=min(vec, na.rm=T)
    max.=max(vec, na.rm=T)
    print(sprintf("mean: %s", mean.))
    print(sprintf("median: %s", median.))
    print(sprintf("min: %s", min.))
    print(sprintf("max: %s", max.))
    print(sprintf("sd: %s", sd(vec, na.rm=T)))
    print(glue("{round(mean.,2)} ({round(min.,2)}-{round(max.,2)})"))
  }else{
    print(table(vec))
  }
}

## coalesce dplyr columns replacing both is.na and "" and NULL.
##' @param colA,colB Two vectors to merge into one
##' @param priority optional string for which column to prioritize when merging columns
##' 
##' @return a single column containing the info from input columns without blanks 
merge_info <- function(colA, colB, priority="left", col_name="merge", warn=T){
  if (priority=="left"){
    # use col A unless blank, then use B
    out_col = ifelse((is.na(colA) | (colA=="") | is.null(colA)), colB, colA)
  }else if (priority=="right"){
    # use col B unless blank, then use A
    out_col = ifelse((is.na(colB) | (colB=="")| is.null(colB)), colA, colB)
  }
  if (any(is.na(out_col)) & warn){
    print("merge_info::")
    warn_na(out_col, col_name)
  }
  return(out_col)
}


## merge data frames while combining info from columns with the same name
##' @param priority string indicating which column to prioritize when merging data. Default: right
##' @param join.type type of dplyr join to perform
##' @param join.by column(s) to join by
##' @param reduce boolean whether to remove original columns or not. default: TRUE
##' 
##' @return merged data frame
merge.combine <- function(df1, df2, join.type="left", join.cols.left="SampleID.short", join.cols.right="SampleID.short", 
                          priority="right", reduce=TRUE, warn=T){
  join.cols.right = as.vector(join.cols.right)
  join.cols.left = as.vector(join.cols.left)
  if (any(is.na(c(join.cols.left, join.type, join.cols.right)))){
    stop("Must specify join.type AND join.by.left AND join.by.right column")
  }
  join_ <- function(df1=df1, df2=df2, join.type=join.type, 
                    join.cols.left=join.cols.left, join.cols.right=join.cols.right){
    join.cols = join.cols.right
    names(join.cols) <- join.cols.left
    if (join.type=="left"){
      df.merge = df1 %>% left_join(df2, by=join.cols, copy=TRUE)
    }else if(join.type=="right"){
      df.merge = df1 %>% right_join(df2, by=join.cols, copy=TRUE)
    }else if(join.type=="full"){
      df.merge = df1 %>% full_join(df2, by=join.cols, copy=TRUE)
    }else if(join.type=="inner"){
      df.merge = df1 %>% inner_join(df2, by=join.cols, copy=TRUE)
    }else{
      stop("Unrecognized join.type in merge.combine")
    }
    return(df.merge)
  }
  comm_cols = intersect(names(df1), names(df2))
  if (length(comm_cols)==0){
    print("WARN: no columns in common")
    return(df.merge)
  }
  df.merge = join_(df1, df2, join.type, join.cols.left, join.cols.right)
  for (col_str in comm_cols){
    if (col_str %!in% c(join.cols.left,join.cols.right)){
      colA = paste0(col_str, ".x")
      colB = paste0(col_str, ".y")
      df.merge = df.merge %>% mutate({{col_str}} := merge_info(df.merge[[colA]],
                                                               df.merge[[colB]],
                                                               priority=priority,
                                                               col_name=col_str,
                                                               warn=warn))  
    }
  }
  if (reduce){
    df.merge = df.merge %>% select(-contains(".x"),-contains(".y"))
  }
  return(df.merge)
}

### standardize column names
##' @param data input variant or sample data frame
##' @param sid.format optional type of SampleID.short to create using make_SIDshort.  
##' @param input.type string indicating type of input data. default: "variants"
##' @param warn whether to print warnings when missing columns are filled with NAs. default=TRUE
##' 
##' @return input data frame with added renamed columns
standardize_names <- function(data, input.type="variants", sid.format="none", warn=T){
  stopifnot(input.type %in% c("variants","samples"))
  
  col_dict = list("SampleID.short"=list("SampleID","sampleNames","RequisitionID"),
                  "PatientID"=list("SubjectID","ExternalID"),
                  "StudyVisit"=list("trialVisitNum"),
                  "SampleType"=list("SpecimenType","Specimen_type","specimenType"))
  
  if (input.type=="variants"){
    col_dict_variants = list("Chromosome"=list("seqnames","chr"),
      "Variant_Classification"=list("Consequence","CNV_Type","CNL_Type"),
      "Hugo_Symbol"=list("SYMBOL","Gene"),
      "VariantID"=list("") #special function below
    )
    col_dict = as.list(c(col_dict, col_dict_variants))
  }
  
  rename_ <- function(df, col_dict, sid.format){
    for (new_col in names(col_dict)){
      if (new_col %!in% names(df)){
        
        if (new_col == "SampleID.short"){
          missing = TRUE
          for (old_col in col_dict[[new_col]]){
            if (old_col %in% names(df)){
              df = df %>% mutate({{new_col}} := as.character(make_SIDshort(df[[old_col]], sid.format=sid.format)))
              missing = FALSE
              break
            }
          }
          if (missing) {stop(glue("Missing {new_col}"))}
          
        }else if(new_col == "VariantID"){
          df$VariantID = NA
          if (any(grepl("zscore", names(df), ignore.case=T))){
            df = tryCatch(df %>% mutate("VariantID"=paste(seqnames,start,end,Variant_Classification, sep=":")),
                          error = function(e) tryCatch(df %>% mutate("VariantID"=paste(Hugo_Symbol,Variant_Classification, sep=":")),
                                                       error = function(e) tryCatch(return(df), 
                                                                                    error = function(e) {print(e)}, # do nothing
                                                                                    finally = message("WARN: Could not create VariantID"))))
          }else{
            df = tryCatch(df %>% mutate("VariantID"=paste(seqnames,start,ref,alt,sep=":")),
                          error = function(e) tryCatch(df %>% mutate("VariantID"=paste(Chromosome,
                                                                                Start_Position,
                                                                                Reference_Allele,
                                                                                Tumor_Seq_Allele, sep=":")),
                                                       error = function(e) tryCatch(return(df), 
                                                                                    error = function(e) {print(e)}, # do nothing
                                                                                    finally = message("WARN: Could not create VariantID"))
                          )
            )
          }
          df$VariantID = as.character(df$VariantID)
          
      }else{ # for all other new_col
          missing = TRUE
          for (old_col in col_dict[[new_col]]){
              if (old_col %in% names(df)){
                  df = df %>% mutate("{new_col}" := as.character(df[[old_col]]))
                  missing = FALSE
                  break
              }
          }
          if(warn & missing){
            message(glue("Standardize.names WARN:Fill all with NA for {new_col}"))
            df = df %>% mutate("{new_col}" := NA)}
      }
    }#else{ # if new_col is already in names(df)
      #print(new_col) 
    #}
  } #end new_cols to create
  return(df)  
  } # end function
  return(rename_(data, col_dict, sid.format))
}

## generate a list of N random colors using RColorBrewer
##' @param names a vector of unique factor levels for which to assign random colors
##' 
##' @return a named list matching factor level to random colors
get_random_color_dict <- function(names) { 
  n <- length(unique(names))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  colors = unlist(sample(col_vector, n))
  names(colors) <- unique(names)
  return(colors)
}


## make a quick boxplot using one categorical predictor and one continuous repsonse and output MWU test
##'@param df data frame containing all inputs
##'@param x independent variable column name
##'@param y dependent variable column name
##'@param facet variable column name by which to reproduce plots
##'@param colors an optional named list matching levels of x to colors
##'@param print.p optional whether to calculate p value and print on the boxplot. default=TRUE
##'
##'@return ggplot object
quick_boxplot <- function(df, x, y, facet=NULL, colors=NULL, 
                          print.p=TRUE, log.axes=FALSE, plot.title=""){
  if (!is.null(facet)){
    df = df %>% select({{x}},{{y}},{{facet}})
  }else{
    df = df %>% select({{x}},{{y}})
  }
  df[df=="NA"] <- NA
  df = df[complete.cases(df),]
  print(glue("Complete cases {y} by {x}: {nrow(df)}"))
  if (is.null(colors)){
    colors = get_random_color_dict(df[[x]])
  }
  stopifnot(is.factor(df[[x]]))
  df[[y]] = as.numeric(as.factor(df[[y]]))
  xlevs = levels(df[[x]])
  counts = df %>% group_by(df[x]) %>% summarize(count=n()) 
  counts = counts %>% mutate("legend" = paste0(counts[[x]]," (n=",count,")"),
                              "color"=recode(counts[[x]], !!!colors, .default = NA_character_),
                             across(everything(), as.character))
  print(counts)
  g <- ggplot(df, aes_string(x=x, y=y, fill=x)) +
    geom_boxplot(outlier.shape=NA) +
    geom_point(position=position_jitterdodge(), size=4, pch=1) +
    scale_fill_manual(labels=counts$legend, breaks=counts[[x]], values=counts$color) +
    ggtitle(plot.title) +
    theme(text = element_text(size = 20),
          #axis.text.x = element_blank(),
          plot.title = element_text(size = 14))
  print("x")
  get.p.table <- function(df, x, xlevs, y){
    filt = ifelse(df[[x]] == xlevs[[1]], TRUE, FALSE)
    lev1_vals = df[[y]][filt]
    print(glue("{xlevs[1]} mean: {mean(lev1_vals, na.rm=T)}"))
    not_lev1_vals = df[[y]][!(filt)]
    print(glue("{xlevs[2:length(xlevs)]} mean: {mean(not_lev1_vals, na.rm=T)}"))
    res <- wilcox.test(lev1_vals, 
                       not_lev1_vals, 
                       exact = FALSE, paired=FALSE)
    print(res)
    p = as.character(round(res$p.value, 3))
    df_p = data.frame(x = xlevs, y=c(NA,NA), label=c("",""))
    if (max(lev1_vals, na.rm = T) >= max(not_lev1_vals, na.rm = T)){
      df_p[2,2] <- max(lev1_vals, na.rm = T)
      df_p[2,"label"] <- paste0("p=",p)
    }else{
      df_p[1,2] <- max(not_lev1_vals, na.rm = T)
      df_p[1,"label"] <- paste0("p=",p)
    }
    names(df_p) <- c(x, y, "label") # use strings as labels
    print(df_p)
    return(df_p)
  }
    
  if (!is.null(facet)){
    print("y")
    df = df %>% mutate({{facet}} := as.factor(df[[facet]]))
    warn_na(df[[facet]], facet)
    p_dfs = list()
    facet.levs = levels(df[[facet]])
    for (t in facet.levs){
      print(glue("Level: {t}"))
      df.facet = df %>% filter(eval(parse(text = glue("{facet} == '{t}'"))))
      df_pi = get.p.table(df.facet, x, xlevs, y)
      df_pi = df_pi %>% mutate({{facet}} := t)
      p_dfs[[t]] <- df_pi
    }
    print(length(p_dfs))
    print(p_dfs[1])
    df_p = do.call(rbind, p_dfs)
    df_p = data.frame(df_p) 
    print(df_p)
    g <- g + 
      facet_wrap(paste("~",facet)) +
      geom_text(data=df_p, aes_string(x=x, y=y, group=facet), size=6, 
                label=df_p$label, fontface="italic")
  }else{
      df_p = get.p.table(df, x, xlevs, y)
      g <- g + 
        geom_text(data=df_p, aes_string(x=x, y=y), size=6, 
                  label=df_p$label, fontface="italic")
  }
  
  if (log.axes){
    g <- g + 
      scale_y_continuous(trans='log2')
  }
  show(g)
}

## quick bar plot
# bar_plot_ <- function(df, x, y, log.axis=FALSE){
#   df.mean = summary %>% group_by(CancerStage) %>% 
#     summarize("Mean.TF"=mean(as.numeric(tumorFraction.mrd)*100000, na.rm=T),
#               "Mean.mrd.score"=mean(as.numeric(mrd.score)*100, na.rm=T))
#   gg <- ggplot(df, aes_string(x=x, y=y)) + 
#     geom_col() + 
#     theme(text=element_text(size=20))
#   if (log.axis){
#     gg <- gg + scale_y_continuous(trans='log2')
#   }
#   show(gg)
# }

## contingency plot
# stacked_bar_ <- function(df, x = "Pathologic.Response", y="mrd.status", 
#                          colors=c("steelblue3","indianred"), y.percent=TRUE){
#   stopifnot(is.factor(df[[x]]))
#   stopifnot(is.factor(df[[y]]))
#   xlevs = levels(df[[x]])
#   ylevs = levels(df[[y]])
#   n.0 = sum(df$Pathologic.Response==xlevs[1])
#   n.1 = sum(df$Pathologic.Response==xlevs[2])
#   names(colors) <- as.character(ylevs)
#   gg <- ggplot(df) +
#     geom_bar(aes(x=df[[x]], fill=df[[y]]), position="fill", color="black") +
#     scale_fill_manual(name = y, values=colors, labels=paste0(ylevs, " (n=", table(df[[y]]), ")")) +
#     scale_x_discrete(labels=paste0(xlevs, "\n(n=", table(df[[x]]), ")")) +
#     theme(text = element_text(size = 16))
#   if (y.percent){
#     gg <- gg + 
#       scale_y_continuous(labels = scales::percent) + 
#       labs(y="Percent of Patients", x=x)
#   }
#   show(gg)
#   # rows response var, cols predictor
#   tab = table(df[[y]], df[[x]])
#   print(tab)
#   fisher.test(as.matrix(tab))
# }

## show which reference values are missing from a vector
check.missing <- function(list.ref, list.test){
  list.ref = unique(list.ref)
  list.test = unique(list.test)
  list.all = unique(c(list.ref, list.test))
  if (any(list.test %!in% list.ref)){
    print(glue("WARN: missing factor levels:"))
    print(list.all)
    print(data.frame("level"=list.all, "in.test.list"=list.all %in% list.test, "in.ref.list"=list.all %in% list.ref))
  }
}
