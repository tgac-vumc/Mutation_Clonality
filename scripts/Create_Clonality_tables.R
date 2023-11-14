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
suppressMessages(suppressWarnings(library(openxlsx)))

#-------------------------------------------------------------------------------
# 1.1 Parse snakemake objects
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_metrics <- snakemake@input[["Clonality_metrics"]]
    input_metrics_MC <- snakemake@input[["Clonality_metrics_MC"]]
    output <-  snakemake@output[["Table"]]
}else{
    input_metrics <-  c('output/IlluminaTrueSightTumor170/Clonality_metrics_LUAD.txt','output/IlluminaTrueSightTumor170/Clonality_metrics_LUSC.txt')
    input_metrics_MC <-  c('output/IlluminaTrueSightTumor170/Clonality_Calls_MC_LUAD.txt','output/IlluminaTrueSightTumor170/Clonality_Calls_MC_LUSC.txt')
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

Clonality_metrics_MC <-
    tibble::tibble(file = input_metrics_MC) %>%
    mutate(panel = purrr::map_chr(file, ~strsplit(.x,'/')[[1]][2]),
           subtype = gsub('.txt','',purrr::map_chr(file, ~strsplit(.x,'_')[[1]][4])),
           # read file and get number of unique samples
           data = purrr::map(file,read.delim)) %>%
    tidyr::unnest() %>%
    mutate(comparison = paste0(Sample1,'-',Sample2))

#-------------------------------------------------------------------------------
# 1.1 Modify and select data
#-------------------------------------------------------------------------------
Clonality_metrics <-
    Clonality_metrics %>%
    left_join(Clonality_metrics_MC, by = c('panel','subtype','True_clonality','comparison')) %>% 
    dplyr::mutate(
               Clonality_TwoMetric_test =  factor(ifelse(pval_Jaccard < 0.05 & LRpvalue < 0.05,'Clonal','Non-Clonal')),
               Clonality_MC_test = factor(
                   case_when(
                       Clonality_MC == 'Clonal' ~ 'Clonal',
                       grepl('non-Clonal',Clonality_MC) ~ 'Non-Clonal',
                       TRUE ~ NA_character_
                       ),
                   levels=c(0,1)),
               Clonality = factor(ifelse(True_clonality == 'Clonal',1,0), levels=c(0,1))) %>%
    select(panel,subtype,comparison,True_clonality,Clonality_MC, Clonality_TwoMetric_test,n1,n2,n_match,Jaccard, reason, Shared_mutations,Driver_information_sample1, Driver_information_sample2 , Other_information_sample1 , Other_information_sample2)


#-------------------------------------------------------------------------------
# 2.1 Write to file, one sheet per subtype
#-------------------------------------------------------------------------------
OUT <- createWorkbook()
for(subtype in c('LUAD','LUSC','NSCLC')){
    addWorksheet(OUT, subtype)
    print(subtype)
    if(subtype == 'NSCLC'){
        data <- Clonality_metrics
    }else{
        data <- Clonality_metrics[Clonality_metrics$subtype == subtype,]
    }
    writeData(OUT, sheet = subtype, x = data)
}
saveWorkbook(OUT, output)
