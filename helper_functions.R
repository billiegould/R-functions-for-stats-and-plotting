### helper functions for post-analysis

library(dplyr)
library(stringr)
library(glue)
library(RColorBrewer)
library(stats)
library(ggsignif)

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
  }else if (sid.format=="remove.A01"){
    return(gsub("A01", "", col))
  }else if (sid.format=="remove.suffix"){
    return(sapply(str_split(col, "_"),"[",1))
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
      warn_(as.character(my.vector[[col]]), col)
    }
  }else{
    warn_(my.vector, my.vector.name)
  }
}

### fill in col NAs with 0s at indexes where values in another column are 0
fill_na_0 <- function(col, index_col, fill.val=0.0){
  col[index_col==0] <- fill.val
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
    out_col = ifelse((is.na(colB) | (colB=="") | is.null(colB)), colA, colB)
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
merge.combine <- function(df1, df2, join.type="left", join.cols=NULL, 
                          priority="right", reduce=TRUE, warn=T){

  if (is.null(join.cols)){
    join.cols = c("SampleID.short")
    names(join.cols) <- c("SampleID.short")
    print("Joining by SampleID.short . . .")
  }else{
    stopifnot(is.vector(join.cols) | is.list(join.cols))
  }
  join_ <- function(df1=df1, df2=df2, join.type=join.type, 
                    join.cols=join.cols){
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
  }else if (warn){
    print("Merging columns:")
    print(comm_cols)
  }
  # merge dfs by specified column names
  df.merge = join_(df1, df2, join.type, join.cols)
  #print(names(df.merge))
  # for overlapping columns, merge_info together
  for (col_str in comm_cols){
    if (col_str %!in% c(join.cols)){
      print(col_str)
      colA = paste0(col_str, ".x")
      colB = paste0(col_str, ".y")
      if (any(is.numeric(df.merge[[colA]]),is.numeric(df.merge[[colB]]))){
        df.merge = df.merge %>% mutate({{col_str}} := merge_info(as.numeric(df.merge[[colA]]),
                                                               as.numeric(df.merge[[colB]]),
                                                               priority=priority,
                                                               col_name=col_str,
                                                               warn=warn))  
      }else{
        df.merge = df.merge %>% mutate({{col_str}} := merge_info(as.character(df.merge[[colA]]),
                                                                 as.character(df.merge[[colB]]),
                                                                 priority=priority,
                                                                 col_name=col_str,
                                                                 warn=warn))
    }}
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

## get p values for plottling
get.p.table.mwu <- function(df, x, xlevs, y){
  df = df %>% filter(!is.na({{y}}))
  filt = ifelse(df[[x]] == xlevs[[1]], TRUE, FALSE)
  lev1_vals = df[[y]][filt]
  print(glue("{xlevs[1]} median: {median(lev1_vals, na.rm=T)}"))
  not_lev1_vals = df[[y]][!(filt)]
  print(glue("{xlevs[2:length(xlevs)]} median: {median(not_lev1_vals, na.rm=T)}"))
  res <- wilcox.test(lev1_vals, 
                     not_lev1_vals, 
                     exact = FALSE, paired=FALSE)
  print(res)
  if (res$p.value < 0.0001){
    p = "<0.0001"
  }else{
    p = as.character(round(res$p.value, 4))
  }
  df_p = data.frame(x = xlevs, y=c(NA,NA), label=c("",""))
  if (max(lev1_vals, na.rm = T) >= max(not_lev1_vals, na.rm = T)){
    df_p[2,2] <- max(lev1_vals, na.rm = T)
    df_p[2,"label"] <- paste0("p=",p)
  }else{
    df_p[1,2] <- max(not_lev1_vals, na.rm = T)
    df_p[1,"label"] <- paste0("p=",p)
  }
  names(df_p) <- c(x, y, "label") # use strings as labels
  #print(df_p)
  return(df_p)
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
quick_boxplot <- function(df, x, y, facet=NULL, colors=NULL, log.0.adj=NULL,
                          print.p=TRUE, log.axes=FALSE, hline=NULL, plot.title="",
                          pt.label.col=NULL, fill.pts=FALSE, pt.fill.col=NULL, y.axis.dec.places=5){
  
  cols = c(x, y, facet, pt.label.col, pt.fill.col)
  stopifnot(all(cols[!is.null(cols)] %in% names(df)))
  df = df[,cols[!is.null(cols)]]
  
  if (!is.null(facet)){
    df = df %>% select({{x}},{{y}},{{facet}})
  }
  if (!is.null(pt.label.col)){
    df = df %>% mutate("label" = as.factor(df[[pt.label.col]]))
    pt.size = 0
  }else{
    df = df %>% mutate("label" = "")
  }
  df[df=="NA"] <- NA
  df = df[complete.cases(df),]
  print(glue("Complete cases {y} by {x}: {nrow(df)}"))
  
  if (is.null(colors)){
    colors = get_random_color_dict(df[[x]])
  }
  
  df[[x]] = as.factor(df[[x]])
  df[[y]] = as.numeric(df[[y]])
  xlevs = levels(df[[x]])
  print(xlevs)
  counts = df %>% group_by(df[x]) %>% summarize(count=n()) 
  counts = data.frame(counts %>% mutate("legend" = paste0(counts[[x]]," (n=",count,")"),
                              "color"=recode(counts[[x]], !!!colors, .default = NA_character_),
                             across(everything(), as.character)))
  rownames(counts) <- counts[[x]]
  counts = counts[xlevs, ]
  print(counts)
  
  ### PLOTTING
  g <- ggplot(df, aes_string(x=x, y=y, fill=x, group=x))
  
  if (!is.null(facet)){
    stopifnot(is.factor(df[[facet]]))
    #df = df %>% mutate({{facet}} := as.factor(df[[facet]]))
    warn_na(df[[facet]], facet)
    p_dfs = list()
    facet.levs = levels(df[[facet]])
    
    if (length(levels(df[[x]])) > 2) {
      #TODO: 
      # comparisons = combn(levels(df.facet[[x]]),
      #                     2, simplify = F)
      # print(comparisons)
      # g <- g +
      #   geom_signif(test="wilcox.test", comparisons = comparisons,
      #               step_increase = 0.2, map_signif_level = TRUE, textsize = 7)
      print("todo: multi comparison P-val across facets")
    }else{
      for (t in facet.levs){
          df.facet = df %>% filter(eval(parse(text = glue("{facet} == '{t}'"))))
          print(df.facet)
          df_pi = get.p.table.mwu(df.facet, x, xlevs, y)
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
          geom_text(data=df_p, aes(x=.data[[x]], y=.data[[y]], group=.data[[facet]]), size=6, 
                    label=df_p[["label"]], fontface="italic")
    }
    
  }else{ #is.na(facet)
    if ((length(levels(df[[x]])) > 2) & print.p){
      comparisons = combn(levels(df[[x]]),
                          2, simplify = F)
      sigFunc = function(x){
        if(x < 0.001){"***"} 
        else if(x < 0.01){"**"}
        else if(x < 0.05){"*"}
        else{NA}}
      g <- g +
        geom_signif(test="wilcox.test", comparisons = comparisons,#[c(3,5)],
                    step_increase = 0.2, 
                    map_signif_level=sigFunc, 
                    textsize = 7)
    }else if (print.p){
      df_p = get.p.table.mwu(df, x, xlevs, y)
      g <- g + 
        geom_text(data=df_p, aes(x=.data[[x]], y=.data[[y]]), size=6, 
                  label=df_p$label, fontface="italic") #+
      # todo: get these labels to match plotting order
        #scale_x_discrete(labels=paste0(counts[[x]], "\n(n=", counts[["count"]],")"))
    }
  }
  
  # for plotting only
  if (log.axes){
    if (is.null(log.0.adj)){
      log.0.adj = min(df[[y]][df[[y]] != 0])/10
    }
    print(glue("for plotting log.0.adj={log.0.adj}"))
    df[[y]] = df[[y]] + log.0.adj
  }
  
  g <- g +
    geom_boxplot(outlier.shape=NA) +
    scale_fill_manual(labels=counts$legend, breaks=counts[[x]], values=counts$color) +
    ggtitle(plot.title) +
    theme(text = element_text(size = 20),
          #axis.text.x = element_blank(),
          plot.title = element_text(size = 14))
  if (fill.pts){
    g <- g + geom_point(aes(shape=.data[[pt.fill.col]]), position=position_jitterdodge(), size=4) +
      scale_shape_manual(values=c(16,1)) #filled, open circle
  }else{
    g <- g + geom_point(position=position_jitterdodge(), size=4, pch=16) #filled circles
  }
  if (log.axes){
    scaleFUN <- function(y) sprintf(glue("%.{y.axis.dec.places}f"), y)
    g <- g + 
      scale_y_continuous(trans='log2', labels=scaleFUN#, 
                         #limits=c(min(df[[y]], na.rm=TRUE), 
                         #          max(df[[y]], na.rm=TRUE))
                         )
  }
  if (!is.null(hline)){
    g <- g + geom_hline(aes(yintercept=hline), linetype="dashed", color="grey", size=1)
  }
  if (!is.null(pt.label.col)){
    g <- g + geom_text(aes(label=.data[[pt.label.col]]), check_overlap = FALSE, position=position_jitter(),size=5)
  }
  #show(g)
  return(g)
}

# boxplot with lines between boxes for each group: pre-nac vs post nac samples, outline points by recur vs. non-recur
# pre_post_plot <- function(df, y, x, x.colors=NULL, line.group.var = "direction", 
#                           line.colors=c("increase"="red","decrease"="blue"), facet=NULL){
#   print(y)
#   stopifnot(is.factor(df[[x]]))
#   
#   if (is.null(x.colors)){
#     xlevs = levels(df[[x]])
#     x.colors = rep("white", length(xlevs))
#     names(x.colors) <- xlevs
#   }
#   
#   formula = as.formula(paste0(y,"~ StudyVisit"))
#   res <- wilcox.test(formula, data = df, exact = FALSE, paired=TRUE)
#   print(res)
#   
#   if (line.group.var=="direction"){
#     df_plot = df[,unique(c(facet, y, x, "PatientID"))]
#   }else{
#     df_plot = df[,unique(c(facet, y, x, line.group.var, "PatientID"))]
#   }
#   
#   df_plot = df_plot %>% mutate("stat" = df_plot[[y]], "group"=df_plot[[x]]) %>% 
#     #arrange(PatientID, {{x}}) %>% # post, pre, alphabetical
#     group_by(PatientID) %>% 
#     mutate(diff = (stat[group=="Post-NAC"]-stat[group=="Pre-NAC"]), #/mean(stat),  
#            direction=ifelse((diff < 0),"decrease","increase")
#     ) %>% ungroup()
#   print(df_plot %>% arrange(PatientID, {{x}}) %>% relocate(diff, direction, stat), n=100)
#   
#   gg <- ggplot(df_plot, aes_string(x=x, y=y)) + 
#     geom_boxplot(outlier.shape=NA) + 
#     #scale_fill_manual(x, values=x.colors) + #TODO fix this
#     #scale_x_discrete(x, colors=x.colors) + 
#     geom_point(size=4, pch=1) +
#     
#     geom_line(aes_string(group="PatientID", colour = line.group.var), linetype = 1) +
#     scale_color_manual(values=line.colors) +
#     scale_y_continuous(trans='log2') + 
#     
#     facet_wrap(paste("~", facet)) +
#     labs(y=y, title=y) + 
#     theme(text = element_text(size = 20),
#           strip.text = element_text(size=20),
#           axis.ticks = element_blank(), 
#           #axis.text.x = element_blank(),
#           legend.position="right")
#   show(gg)
# }

## contingency plot
contingency_plot <- function(df, x = NULL, y=NULL, 
                             colors=c("FALSE"="steelblue3","TRUE"="indianred"), 
                             y.percent=TRUE, facet=NULL, title=""){
  df$X = factor(df[[x]])
  df$Y = factor(df[[y]])
  warn_na(df$X)
  warn_na(df$Y)
  xlevs = levels(df$X)
  ylevs = levels(df$Y)
  n.0 = sum(df$X==xlevs[1])
  n.1 = sum(df$X==xlevs[2])
  if (length(colors) != length(ylevs)){
    colors = get_random_color_dict(names = ylevs)
  }
  gg <- ggplot(df) +
    geom_bar(aes(x=X, fill=Y), position="fill", color="black") +
    scale_fill_manual(name = y, values=colors, labels=paste0(ylevs, " (n=", table(df$Y), ")")) +
    scale_x_discrete(name=x, labels=paste0(xlevs, "\n(n=", table(df$X), ")")) +
    theme(text = element_text(size = 16)) +
    ggtitle(title)
  if (y.percent){
    gg <- gg +
      scale_y_continuous(labels = scales::percent) +
      labs(y="Percent of Patients", x=x)
  }
  
  if (!is.null(facet)){
    ## fix this code
    stop("no facet support yet")
    xlevs = levels(df$X)
    ylevs = levels(df$Y)
    n.0 = sum(df$X==xlevs[1])
    n.1 = sum(df$X==xlevs[2])
    
    gg <- gg + scale_x_discrete(labels=paste0(xlevs)) +
      facet_wrap(paste("~",facet)) +
      geom_text(data=df_p, aes_string(x=x, y=y, group=facet), size=6, 
                                      label=df_p$label, fontface="italic")
    # print p-values per facet
    df = df %>% mutate({{facet}} := as.factor(df[[facet]]))
    warn_na(df[[facet]], facet)
    p_dfs = list()
    facet.levs = levels(df[[facet]])
    for (t in facet.levs){
      print(glue("Level: {t}"))
      df.facet = df %>% filter(eval(parse(text = glue("{facet} == '{t}'"))))
      # rows response var, cols predictor
      tab = table(df.facet$Y, df.facet$X)
      print(tab)
      res = fisher.test(as.matrix(tab))
      print(res)
      if (res$p.value < 0.001){
        p = "<0.001"
      }else{
        p = as.character(round(res$p.value, 3))
      }
      df_p = data.frame(x = facet.levs, y=rep(NA, length(facet.levs)), label=rep("",length(facet.levs)))
      #if (max(lev1_vals, na.rm = T) >= max(not_lev1_vals, na.rm = T)){
      df_p[2,2] <- 0.1 #max(lev1_vals, na.rm = T)
      df_p[2,"label"] <- paste0("p=",p)
      # }else{
      #   df_p[1,2] <- max(not_lev1_vals, na.rm = T)
      #   df_p[1,"label"] <- paste0("p=",p)
      # }
      names(df_p) <- c(x, y, "label") # use strings as labels
      print(df_p)
    }
  }else{
    tab = table(df$Y, df$X)
    print(tab)
    res = fisher.test(as.matrix(tab))
    print(res)
  }
  #show(gg)
  return(gg)
}

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


##### Get optimal and target sens and spec for a data set, optionally print ROC curve.
#####
get_sens_spec <- function(df, label_col, score_col, title=NA, thresh=NA, target_sens=NA, 
                          target_spec=NA, text.cex = 2, print.thres="best"){
  library(pROC)
  df = df %>% rename("label"={{label_col}}, "score"={{score_col}}) %>% 
    mutate("label"=as.character(label))
  print(table(df$label, useNA = "always"))
  
  if (!is.na(thresh)){
    TP = nrow(df %>% filter(label=="TRUE", score >= thresh))
    FN = nrow(df %>% filter(label=="TRUE", score < thresh))
    TN = nrow(df %>% filter(label=="FALSE", score < thresh))
    FP = nrow(df %>% filter(label=="FALSE", score >= thresh))
    sens = TP/(TP + FN)
    spec = TN/(TN + FP)
    print(paste("Threshold Sens:", sens))
    print(paste("Threshold Spec:", spec))
  }
  
  pROC_obj <- roc_(data=df, response="label", predictor="score", 
                   smooth = FALSE, plot=FALSE, direction="<") # controls score lower than cases
  if (!is.na(title)){
    plot.roc(pROC_obj, auc.polygon=TRUE, max.auc.polygon=FALSE, grid=TRUE, 
             #ci = TRUE, ci.type = "bars", #ci.type=c("bars", "shape", "no")
             print.auc=TRUE, print.thres=print.thres, main=title, asp = NA, print.auc.cex=text.cex)
  }
  print(coords(pROC_obj, x=print.thres))
  
  if (!is.na(target_sens)){
    print(sprintf("target SENS %s", target_sens))
    print(coords(pROC_obj, x=target_sens, input="sensitivity", 
                 ret=c("threshold","specificity", "sensitivity")))
  }
  if (!is.na(target_spec)){
    print(sprintf("target SPEC %s", target_spec))
    print(coords(pROC_obj, x=target_spec, input="specificity", 
                 ret=c("threshold","specificity", "sensitivity")))
  }
  print(paste0("AUC: ", auc(pROC_obj)))
}