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

#-------------------------------------------------------------------------------
# 1.1 Parse snakemake objects
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_metrics <- snakemake@input[["Metrics"]]
    input_metrics_MC <- snakemake@input[["Metrics_MC"]]
    output <-  snakemake@output[["Evaluation_metrics"]]
    output_raw <-  snakemake@output[["Raw_Table"]]

    output_barchart <- snakemake@output[['Barchart_performance']]
}else{
    input_metrics <- Sys.glob('output/*/Clonality_metrics_*.txt')
    input_metrics_MC <- Sys.glob('output/*/MC_Clonality_metrics_*.txt')
    output <- 'output/Tables/TableXX_Sensitivity_Specificity.txt'
    output_raw <- 'output/Tables/TableXX_Raw_Table.txt'
    output_barchart <- 'output/Figures/FigureXX_Barchart_Performance.pdf'
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
    tidyr::unnest() %>%
    mutate(comparison = paste0(Tumor1,'-',Tumor2))  %>%
    select(-file)

Clonality_metrics_MC <-
    tibble::tibble(file = input_metrics_MC) %>%
    mutate(panel = purrr::map_chr(file, ~strsplit(.x,'/')[[1]][2]),
           subtype = gsub('.txt','',purrr::map_chr(file, ~strsplit(.x,'_')[[1]][4])),
           # read file and get number of unique samples
           data = purrr::map(file,read.delim)) %>%
    tidyr::unnest() %>%
    mutate(comparison = paste0(Sample1,'-',Sample2))  %>%
    select(-file) 




#-------------------------------------------------------------------------------
# 2.1 Determine clonality 
#-------------------------------------------------------------------------------
# Determine clonality using Two metric classification:
# 1) pval_Jaccard < 0.05
# 2) LRpvalue < 0.05
# 3) pval_Jaccard < 0.05 & LRpvalue < 0.05
#
# Also add MC algorithm results
Clonality_metrics_All <-
    Clonality_metrics %>% 
    # Join molecular classification algorithm
    left_join(Clonality_metrics_MC) %>% 
    dplyr::mutate(
               Clonality_Jaccard_test = factor(ifelse(pval_Jaccard < 0.05,1,0), levels=c(0,1)),
               Clonality_LR_test = factor(ifelse(LRpvalue < 0.05,1,0), levels=c(0,1)),
               Clonality_TwoMetric =  factor(ifelse(pval_Jaccard < 0.05 & LRpvalue < 0.05,1,0), levels=c(0,1)),
               Clonality_MC_test = factor(
                   case_when(
                       Clonality_MC == 'Clonal' ~ 1,
                       grepl('non-Clonal',Clonality_MC) ~ 0,
                       ),
                   levels=c(0,1)),
               Clonality = factor(ifelse(True_clonality == 'Clonal',1,0), levels=c(0,1)))

#-------------------------------------------------------------------------------
# 2.2 Count number of 'Inconclusive' and Misclassified by MC algorithm
#-------------------------------------------------------------------------------
Inconclusives <-
    Clonality_metrics_All %>%
    group_by(panel,subtype,Clonality_MC,.drop=FALSE) %>%
    summarise(Inconclusive = dplyr::n()) %>%
    filter(Clonality_MC == 'Inconclusive' ) %>%


Correctly_classified_MC <-
    Clonality_metrics_All %>%
    mutate(Correctly_Classified = Clonality_MC_test == Clonality) %>%
    group_by(panel,subtype,Correctly_Classified,.drop=FALSE) %>%
    summarise(Correctly_Classified_MC = dplyr::n()) %>%
    filter(Correctly_Classified == TRUE) %>% select(-Correctly_Classified)
Misclassified_MC <-
    Clonality_metrics_All %>%
    mutate(MisClassification = Clonality_MC_test != Clonality) %>%
    group_by(panel,subtype,MisClassification,.drop=FALSE) %>%
    summarise(Misclassified_MC = dplyr::n()) %>%
    filter(MisClassification == TRUE) %>% select(-MisClassification)

Misclassified_FN_MC <-
    Clonality_metrics_All %>%
    filter(Clonality_MC_test == '0' & Clonality == '1') %>%
    mutate(reason = as.factor(paste0('FN: ',reason))) %>%
    group_by(panel,subtype) %>%
    mutate(FN_MC = dplyr::n())%>%
    group_by(panel,subtype,reason,FN_MC) %>%
    summarise(count = dplyr::n()) %>%
    tidyr::pivot_wider(names_from = reason, values_from = count) %>% 
    replace(is.na(.), 0)

Misclassified_FP_MC <-
    Clonality_metrics_All %>%
    filter(Clonality_MC_test == '1' & Clonality == '0') %>%
    mutate(reason = as.factor(paste0('FP: ',reason))) %>%
    group_by(panel,subtype) %>%
    mutate(FP_MC = dplyr::n()) %>%
    group_by(panel,subtype,reason,FP_MC) %>%
    summarise(count = dplyr::n()) %>%
    tidyr::pivot_wider(names_from = reason, values_from = count) %>%
    replace(is.na(.), 0)

nTP <-
    Clonality_metrics_All %>%
    filter( Clonality_MC_test == '1' &  Clonality == '1') %>%
    group_by(panel,subtype) %>%
    summarise(nTP = dplyr::n())

nNegativePred <-
    Clonality_metrics_All %>%
    filter( Clonality_MC_test == '0') %>%
    group_by(panel,subtype) %>%
    summarise(nNegativesPred = dplyr::n())

nPositivePred <-
    Clonality_metrics_All %>%
    filter( Clonality_MC_test == '1') %>%
    group_by(panel,subtype) %>%
    summarise(nPositivesPred = dplyr::n())

nTN <-
    Clonality_metrics_All %>%
    filter(Clonality_MC_test == '0' &  Clonality == '0') %>%
    group_by(panel,subtype) %>%
    summarise(nTN = dplyr::n())

nNegativePred_TwoMetric <-
    Clonality_metrics_All %>%
    filter( Clonality_TwoMetric == '0') %>%
    group_by(panel,subtype) %>%
    summarise(nNegativesPred_TwoMetric = dplyr::n())

nPositivePred_TwoMetric <-
    Clonality_metrics_All %>%
    filter( Clonality_TwoMetric == '1') %>%
    group_by(panel,subtype) %>%
    summarise(nPositivesPred_TwoMetric = dplyr::n())

#-------------------------------------------------------------------------------
# 2.2 Calculate specificity and sensitivity and add statistics
#-------------------------------------------------------------------------------
# Fetch evaluation metrics per panel and subtype
# For now, only include MC algorithm
Evaluation_metrics <-
    Clonality_metrics_All %>%
    group_by(panel,subtype) %>%
    summarise(
        Sensitivity_MC = caret::sensitivity(table(Clonality,Clonality_MC_test), negative = '0', positive = '1' , na.rm=T),
        Specificity_MC = caret::specificity(table(Clonality,Clonality_MC_test), negative = '0', positive = '1' , na.rm=T),
        NComparisons = dplyr::n()) %>%
    ungroup() %>%
    left_join(Inconclusives) %>%
    left_join(Misclassified_MC) %>%
    left_join(Correctly_classified_MC) %>%
    left_join(Misclassified_FN_MC) %>%
    left_join(Misclassified_FP_MC) %>%
    left_join(nTP) %>%
    left_join(nTN) %>%
    left_join(nNegativePred) %>%
    left_join(nPositivePred) %>%
    replace(is.na(.), 0)

Evaluation_metrics %>% View()
#-------------------------------------------------------------------------------
# 2.3 Create barplot performance
#-------------------------------------------------------------------------------
Evaluation_metrics %>%
    ggplot(mtcars, aes(am)) +
    scale_y_continuous(labels=percent_format()) + facet_grid(~subtype)

unique(Evaluation_metrics$panel)


Plot_data <- Evaluation_metrics %>%
    filter(panel != 'KappaHyperExome') %>%
    select(subtype,panel,Correctly_Classified_MC,Inconclusive,Misclassified_MC) %>%
    tidyr::pivot_longer(cols = -c(subtype,panel)) %>%
    mutate(name = factor(recode(name, 'Correctly_Classified_MC' = 'Correct','Misclassified_MC' = 'Incorrect'), levels = c('Incorrect','Inconclusive','Correct')),
           panel = factor(panel, levels = c("MayoCompleteLungCancerPanel","InhouseLungPanel","FoundationOneCDx","IlluminaTrueSightTumor500"))) %>%
    group_by(subtype,panel) %>%
    mutate(Percentage = paste0(round(100*(value/sum(value))),'%')) %>%
    ungroup()

print(Plot_data, n=24)

pdf(output_barchart,height = 5.8*(243.885/235.320), width = 8*(38.621/73.440))
Plot_data %>%
    ggplot(aes(panel,value,fill=name)) +
    geom_bar(position="fill", stat="identity",colour="black") +
    scale_y_continuous(labels=scales::percent_format()) +
    scale_fill_manual(values = c('#ff4d55ff','#ccccccff','#51be60ff')) +
    facet_grid(~subtype) +
    theme_classic(base_size = 14) +
    theme(legend.position = 'top',axis.text.x = element_text(angle = 45, hjust=1)) +
    labs(y='',x='',fill = '')
dev.off()

#-------------------------------------------------------------------------------
# 4.1 Write to file
#-------------------------------------------------------------------------------
write.table(Evaluation_metrics, file = output, sep = '\t', row.names = F,quote = F  )
write.table(Clonality_metrics_All, file = output_raw, sep = '\t', row.names = F,quote = F  )




