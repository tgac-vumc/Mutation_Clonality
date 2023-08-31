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
suppressMessages(suppressWarnings(library(pROC)))

#-------------------------------------------------------------------------------
# 1.1 Parse snakemake objects
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_metrics <- snakemake@input[["Metrics"]]
    output <-  snakemake@output[["Evaluation_metrics"]]
}else{
    input_metrics <-  c('output/InhouseLungPanel/Clonality_metrics.txt','output/IlluminaFocusPanel/Clonality_metrics.txt')
    output <- 'output/Tables/TableXX_Sensitivity_Specificity.txt'
}

#-------------------------------------------------------------------------------
# 1.1 Read data 
#-------------------------------------------------------------------------------
# read metrics
Clonality_metrics <-
    tibble::tibble(file = input_metrics) %>%
    mutate(panel = purrr::map_chr(file, ~strsplit(.x,'/')[[1]][2]),
           # read file and get number of unique samples
           data = purrr::map(file,read.delim)) %>%
    tidyr::unnest()

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
               Clonality_Jaccard_test = ifelse(pval_Jaccard < 0.05,1,0),
               Clonality_LR_test = ifelse(LRpvalue < 0.05,1,0),
               Clonality_Combined_test =  ifelse(pval_Jaccard < 0.05 & LRpvalue < 0.05,1,0),
               Clonality = ifelse(True_clonality == 'Clonal',1,0))


#-------------------------------------------------------------------------------
# 2.2 Calculate specificity and sensitivity
#-------------------------------------------------------------------------------
# Fetch evaluation metrics per panel
Evaluation_metrics <-
    Clonality_metrics %>%
    group_by(panel) %>%
    summarise(
        Sensitivity_LR = caret::sensitivity(table(Clonality,Clonality_LR_test)),
        Specificity_LR = caret::specificity(table(Clonality,Clonality_LR_test)),
        Sensitivity_Jaccard = caret::sensitivity(table(Clonality,Clonality_Jaccard_test)),
        Specificity_Jaccard = caret::specificity(table(Clonality,Clonality_Jaccard_test)),
        Sensitivity_Combined = caret::sensitivity(table(Clonality,Clonality_Combined_test)),
        Specificity_Combined = caret::specificity(table(Clonality,Clonality_Combined_test))) 


#-------------------------------------------------------------------------------
# 2.3 Perform ROC analysis with number of matching mutations
#-------------------------------------------------------------------------------
# Perform ROC analysis with number of mutations and calculate AUC
ROC_analysis <-
    Clonality_metrics %>%
    group_by(panel) %>%
    summarise(
        ROC_model = list(roc(response=Clonality,predictor=n_match)),
        AUC = purrr::map_dbl(ROC_model,auc)) 

# Join ROC AUC
Evaluation_metrics <-
    Evaluation_metrics %>%
    left_join(ROC_analysis[,c('panel','AUC')])
    
#-------------------------------------------------------------------------------
# 3.1 Write to file
#-------------------------------------------------------------------------------
write.table(Evaluation_metrics, file = output, sep = '\t', row.names = F,quote = F  )
