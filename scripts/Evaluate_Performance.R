#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Evaluate performance.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Call clonality and evaluate performance
#
# Author: Jurriaan Janssen (j.janssen4@amsterdamumc.nl)
#
# TODO:
# 1) 
#
# History:
#  29-08-2023: File creation, write code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(pROC)))

#-------------------------------------------------------------------------------
# 1.1 Parse snakemake objects
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_metrics <- snakemake@input[["Metrics"]]
    input_mutations <- snakemake@input[["Shared_Mutations"]]
    input_table <- snakemake@input[["Table"]]
    output <-  snakemake@output[["Evaluation_metrics"]]
    output_boxplots <- snakemake@output[['Boxplot_nMatches']]
    output_piecharts <- snakemake@output[['Piechart_SharedMutations']]
    output_piechart_evaluablepatients <- snakemake@output[['Piechart_EvaluablePatients']]
    output_performance <- snakemake@output[['Barchart_performance']]
}else{
    input_metrics <-  c('output/KappaHyperExome/Clonality_metrics_LUAD.txt','output/KappaHyperExome/Clonality_metrics_LUSC.txt','output/InhouseLungPanel/Clonality_metrics_LUAD.txt','output/InhouseLungPanel/Clonality_metrics_LUSC.txt')
    input_mutations <-  c('output/KappaHyperExome/Shared_mutations_LUAD.txt','output/KappaHyperExome/Shared_mutations_LUSC.txt','output/InhouseLungPanel/Shared_mutations_LUAD.txt','output/InhouseLungPanel/Shared_mutations_LUSC.txt')
    input_table <- 'output/Tables/TableXX_NumberOfEvaluableSamples.txt'
    output <- 'output/Tables/TableXX_Sensitivity_Specificity.txt'
}

#-------------------------------------------------------------------------------
# 1.1 Read data 
#-------------------------------------------------------------------------------
# read metrics
Clonality_metrics <-
    tibble::tibble(file = input_metrics) %>%
    mutate(panel = purrr::map_chr(file, ~strsplit(.x,'/')[[1]][2]),
           subtype = gsub('.txt','',purrr::map_chr(file, ~strsplit(.x,'_')[[1]][3])),
           # read file and get number of unique samples
           data = purrr::map(file,read.delim)) %>%
    tidyr::unnest()

# read shared_mutations
Shared_mutations <-
    tibble::tibble(file = input_mutations) %>%
    mutate(panel = purrr::map_chr(file, ~strsplit(.x,'/')[[1]][2]),
           subtype = gsub('.txt','',purrr::map_chr(file, ~strsplit(.x,'_')[[1]][3])),
           # read file and get number of unique samples
           data = purrr::map(file,~read.delim(.x , colClasses = rep('character',5)))) %>%
    tidyr::unnest()

nEvaluable_patients <- read.delim(input_table)
#-------------------------------------------------------------------------------
# 2.1 Determine clonality 
#-------------------------------------------------------------------------------
# Determine clonality using three different methods:
# 1) pval_Jaccard < 0.05
# 2) LRpvalue < 0.05
# 3) pval_Jaccard < 0.05 & LRpvalue < 0.05
Clonality_metrics <-
    Clonality_metrics %>%
    dplyr::mutate(
               Clonality_Jaccard_test = factor(ifelse(pval_Jaccard < 0.05,1,0), levels=c(0,1)),
               Clonality_LR_test = factor(ifelse(LRpvalue < 0.05,1,0), levels=c(0,1)),
               Clonality_Combined_test =  factor(ifelse(pval_Jaccard < 0.05 & LRpvalue < 0.05,1,0), levels=c(0,1)),
               Clonality = factor(ifelse(True_clonality == 'Clonal',1,0), levels=c(0,1)))

#-------------------------------------------------------------------------------
# 2.2 Calculate specificity and sensitivity
#-------------------------------------------------------------------------------
# Fetch evaluation metrics per panel and subtype
Evaluation_metrics <-
    Clonality_metrics %>%
    group_by(panel,subtype) %>%
    summarise(
        Sensitivity_LR = caret::sensitivity(table(Clonality,Clonality_LR_test)),
        Specificity_LR = caret::specificity(table(Clonality,Clonality_LR_test)),
        Sensitivity_Jaccard = caret::sensitivity(table(Clonality,Clonality_Jaccard_test)),
        Specificity_Jaccard = caret::specificity(table(Clonality,Clonality_Jaccard_test)),
        Sensitivity_Combined = caret::sensitivity(table(Clonality,Clonality_Combined_test)),
        Specificity_Combined = caret::specificity(table(Clonality,Clonality_Combined_test))) %>%
    ungroup()

#-------------------------------------------------------------------------------
# 2.3 Perform ROC analysis with number of matching mutations
#-------------------------------------------------------------------------------
# Perform ROC analysis with number of mutations and calculate AUC
ROC_analysis <-
    Clonality_metrics %>%
    group_by(panel,subtype) %>%
    summarise(
        ROC_model = list(roc(response=Clonality,predictor=n_match)),
        AUC = purrr::map_dbl(ROC_model,auc)) 

# Join ROC AUC
Evaluation_metrics <-
    Evaluation_metrics %>%
    left_join(ROC_analysis[,c('panel','subtype','AUC')]) 

#-------------------------------------------------------------------------------
# 3.1 Plot data
#-------------------------------------------------------------------------------
# 3.1.1 Plot boxplots number of matching mutations per panel for LUAD/LUSC
#       in clonal vs non-clonal
#-------------------------------------------------------------------------------
palette <- palette.colors(palette = "Okabe-Ito")
palette <- palette[!names(palette) %in% c('black','yellow','skyblue')]
palette_clonal <- palette[c('bluishgreen','vermillion')]

pdf(output_boxplots,width=10,height=8)
#pdf('/net/beegfs/cfg/tgac/jjanssen4/boxplot_test.pdf', width = 6, height = 4)
Clonality_metrics %>%
    ggplot(aes(panel, n_match,fill=True_clonality)) +
    geom_boxplot() +
    scale_fill_manual(values = as.vector(palette_clonal)) +
    facet_wrap(~subtype) +
    labs(fill='Clonality',x='Panel',y='Number of matching mutations') +
    theme_bw(base_size = 12) +
    theme(axis.text.x=element_text(size=12,angle=45,hjust=1),
          legend.position='top') +
    scale_y_continuous(trans='log10')
dev.off()
#-------------------------------------------------------------------------------
# 3.1.2 Plot Piechart of %-age of shared mutations in per panel for LUAD/LUSC
#       in clonal vs non-clonal
#-------------------------------------------------------------------------------
pdf(output_piecharts, height = 5, width = 5)
Shared_mutations %>%
    filter(panel == 'KappaHyperExome') %>%
    # Create gene categories
    mutate(Gene = case_when(
               Hugo_Symbol == 'KRAS' ~ 'KRAS',
               Hugo_Symbol == 'EGFR' ~ 'EGFR',
               Hugo_Symbol == 'TP53' ~ 'TP53',
               Hugo_Symbol == 'PIK3CA' ~ 'PIK3CA',
               Hugo_Symbol == 'NFE2L2' ~ 'NFE2L2',
               TRUE ~ 'Other')) %>%
    # Calculate number of shared mutations
    group_by(subtype,True_clonality) %>%
    mutate(Ntot = dplyr::n()) %>%
    group_by(subtype,True_clonality,Gene) %>%
    # Calculate percentage
    dplyr::summarize(Percentage = dplyr::n() / Ntot * 100) %>%
    unique() %>%
    # Plot piecharts
    ggplot(aes(x="", y=Percentage, fill=Gene)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    facet_grid(rows=vars(subtype) , cols = vars(True_clonality), switch = 'y') +
    theme_bw(base_size = 12) +
    theme(axis.ticks=element_blank(), axis.title=element_blank(), 
          axis.text.y = element_blank(), panel.grid  = element_blank(),
          axis.text.x = element_blank()) +
    scale_fill_manual(values=as.vector(palette))
dev.off()


#-------------------------------------------------------------------------------
# 3.1.2 Plot Piechart of %-age of evaluable patients per panel/subtype
#-------------------------------------------------------------------------------
pdf(output_piechart_evaluablepatients, height = 3 , width = 12)
nEvaluable_patients %>%
    filter(subtype != 'NSCLC') %>%
    # create factor ordered by percentage
    mutate(panel = factor(panel, levels = rev(unique(nEvaluable_patients$panel)))) %>%
    # calculate percentages
    group_by(subtype,panel) %>%
    summarize(Evaluable = Nsample_pct,
           `Not Evaluable` = 100 - Nsample_pct) %>%
    select(panel,subtype,Evaluable,`Not Evaluable`) %>%
    # retrieve long format
    tidyr::pivot_longer(cols = -c(panel,subtype), values_to = 'Percentage') %>%
    # Plot piecharts
    ggplot(aes(x="", y=Percentage, fill=name)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    facet_grid(rows=vars(subtype) , cols = vars(panel), switch = 'y') +
    theme_bw(base_size = 12) +
    theme(axis.ticks=element_blank(), axis.title=element_blank(), 
          axis.text.y = element_blank(), panel.grid  = element_blank(),
          axis.text.x = element_blank(),
          legend.title= element_blank(),
          legend.position='top',
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,-10,-10),
          strip.text.x = element_text(size = 5)) +
    scale_fill_manual(values=c('darkgrey','lightgrey'))
dev.off()

 
#-------------------------------------------------------------------------------
# 3.1.2 Plot Barchart performance for each panel per subtype
#-------------------------------------------------------------------------------
pdf(output_performance,width=12, height = 4)
Evaluation_metrics %>%
    select(panel,subtype,Sensitivity_Combined,Specificity_Combined) %>%
    rename(Sensitivity = Sensitivity_Combined,Specificity =Specificity_Combined) %>%
    tidyr::pivot_longer(cols = -c(panel,subtype), names_to = 'Measure') %>%
    mutate(panel = factor(panel,levels = unique(arrange(Evaluation_metrics[Evaluation_metrics$subtype == 'LUAD',],desc(Sensitivity_Combined),desc(Specificity_Combined))$panel))) %>%
    ggplot(aes(panel, value, fill = Measure))  +
    geom_bar(position="dodge", stat="identity") +
    labs(fill='Measure',x='Panel',y='Value') +
    theme_bw(base_size = 12) +
    theme(axis.text.x=element_text(size=12,angle=45,hjust=1),
          legend.position='top',
          legend.box.margin=margin(-10,-10,-10,-10)) +
    scale_fill_manual(values=as.vector(palette)) +
    facet_grid(~subtype)
dev.off()
#-------------------------------------------------------------------------------
# 4.1 Write to file
#-------------------------------------------------------------------------------
write.table(Evaluation_metrics, file = output, sep = '\t', row.names = F,quote = F  )
