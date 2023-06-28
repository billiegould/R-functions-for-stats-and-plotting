## code for finding mutational associations with patient groups, plotting results
## B.Gould
library(tidyverse)
library(glue)
library(RColorBrewer)
library(stats)

options(stringsAsFactors = FALSE)

get_fisher_data <- function (df_all_variants, df_samples){
  print(nrow(df_samples))
  df_variants_sel = df_all_variants %>% filter(SampleID %in% df_samples$SampleID)
  genes = unique(df_variants_sel$SYMBOL)
  #print(genes)
  patients = unique(as.character(df_samples$PatientID))
  fisher_data = list()
  for (patient in patients){
    #print(patient)
    df_patient = df_variants_sel %>% filter(PatientID==patient)
    #print(nrow(df_patient))
    if (nrow(df_patient) == 0){
      patient_data_by_gene_bool = rep(FALSE,length(genes))
    }else{
      patient_data_by_gene_bool = vector()
      for (gene in genes){
        #print(gene)
        df_gene = df_patient %>% filter(SYMBOL==gene)
        #print(nrow(df_gene))
        if (nrow(df_gene) > 0){
          mutated = TRUE
        }else{
          mutated = FALSE
        }
        patient_data_by_gene_bool = c(patient_data_by_gene_bool, mutated)
      }}
    stopifnot(length(patient_data_by_gene_bool)==length(genes))
    fisher_data[[patient]] <- patient_data_by_gene_bool
  }
  
  df_fisher_data = do.call(cbind, fisher_data)
  rownames(df_fisher_data) = genes
  return(df_fisher_data)
}

assoc_test <- function(df_fisher_data, df_clinical){
  # order stage labels by column pids
  response_filter = df_clinical$Pathologic.Response[match(colnames(df_fisher_data), df_clinical$PatientID)]
  response_filter = (response_filter=="pCR")
  print(response_filter)
  
  ft <- function(row){
    mat = matrix(row, nrow=2, byrow=TRUE)
    #print(mat)
    return(fisher.test(mat)$p.value)}
  
  # STAGE
  df_fisher_test_response = NULL
  if (sum(response_filter) >= 2){
    pcr_mutated = rowSums(df_fisher_data[,response_filter])
    #print(late_mutated)
    pcr_nonmutated = rowSums(!(df_fisher_data[,response_filter]))
    #print(late_nonmutated)
    rd_mutated = rowSums(df_fisher_data[,!(response_filter)])
    #print(early_mutated)
    rd_nonmutated = rowSums(!(df_fisher_data[,!(response_filter)]))
    #print(early_nonmutated)
    df_fisher_test_response = data.frame("rd_nonmutated"=rd_nonmutated, "pcr_nonmutated"=pcr_nonmutated,
                                         "rd_mutated"=rd_mutated, "pcr_mutated"=pcr_mutated)
    #print(head(df_fisher_test_stage))
    df_fisher_test_response$p_value = apply(df_fisher_test_response, 1, ft)
    df_fisher_test_response$signif = ifelse(df_fisher_test_response$p_value <0.05, "*", "")
    print("Pathogenic.Response")
    print(df_fisher_test_response)
  }
  # # add log odds ratios
  # df_urine_stage$OR = log(df_urine_stage$late_freq / df_urine_stage$early_freq)
  # # reorder the columns
  # df_urine_stage2 = df_urine_stage %>% 
  #     filter(pct_mutated>(2/11)) %>% 
  #     arrange(OR)                     
  return(df_fisher_test_response)
}



#### current not so good plotting code for this data #####
library(grid)
library(gridExtra)

# FFPE
df_ffpe_fisher = get_fisher_data(all.var_internal_no_patch, 
                                 df_ffpe_care %>% filter(PatientID %!in% c("UTUC-015","UTUC-019")))
#head(df_ffpe_fisher)
result = assoc_test(df_ffpe_fisher, df_clinical_updated)
df_ffpe_stage = result$stage
out.dat = do.call(cbind, result)
write.csv(out.dat, "./manuscript/ffpe_gene_mut_freq_table_03172023.csv")

# format ffpe data for plotting
df_ffpe_stage$early_freq = df_ffpe_stage$early_mutated / (df_ffpe_stage$early_mutated + df_ffpe_stage$early_nonmutated)
df_ffpe_stage$late_freq = df_ffpe_stage$late_mutated / (df_ffpe_stage$late_mutated + df_ffpe_stage$late_nonmutated)
df_ffpe_stage$gene = row.names(df_ffpe_stage)
# add log odds ratios
#df_ffpe_stage$OR = log(df_ffpe_stage$late_freq / df_ffpe_stage$early_freq)

df_ffpe_stage = df_ffpe_stage %>% mutate(pct_mutated=(early_mutated+late_mutated)/(early_mutated+late_mutated+early_nonmutated+late_nonmutated))

# reorder the columns
df_ffpe_stage2 = df_ffpe_stage %>% 
  filter(pct_mutated>(3/30)) %>% 
  arrange(late_mutated)
nrow(df_ffpe_stage2) # 16 genes
df_ffpe_stage = df_ffpe_stage2 %>% mutate(gene=factor(gene, levels=gene))

# format plasma data for plotting
df_plasma_stage$early_freq = df_plasma_stage$early_mutated / (df_plasma_stage$early_mutated + df_plasma_stage$early_nonmutated)
df_plasma_stage$late_freq = df_plasma_stage$late_mutated / (df_plasma_stage$late_mutated + df_plasma_stage$late_nonmutated)
df_plasma_stage$gene = row.names(df_plasma_stage)
# add log odds ratios
#df_plasma_stage$OR = log(df_plasma_stage$late_freq / df_plasma_stage$early_freq)
df_plasma_stage = df_plasma_stage %>% mutate(pct_mutated=(early_mutated+late_mutated)/(early_mutated+late_mutated+early_nonmutated+late_nonmutated)) %>% 
  filter(pct_mutated>(1/30)) %>% 
  arrange(late_mutated) %>% mutate(gene=factor(gene, levels=gene))
nrow(df_plasma_stage) # 9 genes

# vertical back to back bar charts
# https://stackoverflow.com/questions/18265941/two-horizontal-bar-charts-with-shared-axis-in-ggplot2-similar-to-population-pyr

gene_mut_freq_plot <- function(df_stage, plot.title="gene mutation frequencies by group"){
  # label column
  g.mid<-ggplot(df_stage, aes(x=1, y=gene)) + 
    geom_text(aes(label=gene, size=11)) +
    #geom_segment(aes(x=0.94,xend=0.96,yend=gene))+
    #geom_segment(aes(x=1.04,xend=1.065,yend=gene))+
    ggtitle("") +
    ylab(NULL) +
    scale_x_continuous(expand=c(0,0),limits=c(0.98,1.02)) +
    theme(axis.title=element_blank(),
          panel.grid=element_blank(),
          axis.text.y=element_text(),
          axis.ticks.y=element_blank(),
          panel.background=element_blank(),
          axis.text.x=element_text(color=NA),
          axis.ticks.x=element_line(color=NA),
          plot.margin = unit(c(1,-1,1,-1), "mm"),
          legend.position="none")
  #left side Invasive
  g1 <- ggplot(data = df_stage, aes(x = gene, y = late_freq)) +
    geom_bar(stat = "identity", fill="coral1") + #ggtitle("Invasive UTUC\n(pTa-pT1+N0, n=16) ") +
    theme(plot.title = element_text(size = 18, hjust=1),
          axis.line.x = element_line(colour = "black", size=1.5),
          axis.title.x = element_blank(),  
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          plot.margin = unit(c(1,0.75,1,1), "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position="none",
          axis.text = element_text(size = 24)) +
    scale_y_reverse(limits=c(1,0)) + 
    #scale_y_continuous(position="right") +
    coord_flip()
  # right side Non-invasive
  g2 <- ggplot(data = df_stage, aes(x = gene, y = early_freq)) + xlab(NULL)+
    geom_bar(stat = "identity",fill="darkturquoise") + #ggtitle("Non-Invasive UTUC\n(pT2-4/N+, n=16)") +  #n=16 tissue, n=14 plasma
    theme(plot.title = element_text(size = 18),
          axis.line.x = element_line(colour = "black", size=1.5),
          axis.title.x = element_blank(), axis.title.y = element_blank(), 
          plot.margin = unit(c(1,1,1,0.75), "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          #axis.line.y = element_line(colour = "black", size=1.5),
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          legend.position="none",
          axis.text = element_text(size = 24)) +
    scale_y_continuous(limits=c(0,1)) +
    coord_flip() #+ ylim(0,1)
  
  gg1 <- ggplot_gtable(ggplot_build(g1))
  gg2 <- ggplot_gtable(ggplot_build(g2))
  gg.mid <- ggplot_gtable(ggplot_build(g.mid))
  # grid.arrange(gg1,gg.mid,gg2,ncol=3,widths=c(6/13,1/13,6/13), 
  #              bottom=textGrob("Percent of Samples Mutated", gp=gpar(fontsize=16)))
  plt <- grid.arrange(gg1,gg2,ncol=2,widths=c(6/12,6/12), 
                      bottom=textGrob("Proportion of Samples Mutated", gp=gpar(fontsize=28)),
                      top=textGrob(plot.title, gp=gpar(fontsize=28)))
  # for manual labels
  print("Labels (top to bottom)")
  print(rev(df_stage$gene))
  return(plt)    
}
