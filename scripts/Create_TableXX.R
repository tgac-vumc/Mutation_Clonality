#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create_TableXX.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Create table of number of evaluable samples per panel
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

#-------------------------------------------------------------------------------
# 1.1 Parse snakemake objects
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_Allmutations <- snakemake@input[["AllMutations"]]
    input_mutations <- snakemake@input[["Mutations"]]
    input_overview <- snakemake@input[["SampleOverview"]]
    output <-  snakemake@output[["Table"]]
}else{
    input_Allmutations <-'data/TRACERx100_supplement/nejmoa1616288_appendix_2.xlsx' #'data/TRACERx421_supplement/mutTableAll.cloneInfo.20220726.txt'
    input_mutations <-  c('data/KappaHyperExome/Selected_mutations_LUAD.txt','data/IlluminaFocusPanel/Selected_mutations_LUAD.txt','data/KappaHyperExome/Selected_mutations_LUSC.txt','data/IlluminaFocusPanel/Selected_mutations_LUSC.txt')
    input_overview <- 'data/TRACERx100_supplement/nejmoa1616288_appendix_2.xlsx'#'data/TRACERx421_supplement/20221109_TRACERx421_all_patient_df.rds'
    output <- 'output/Tables/TableXX_NumberOfEvaluableSamples.txt'

}

#-------------------------------------------------------------------------------
# 1.1 Read data and count number of unique samples
#-------------------------------------------------------------------------------
# Read sample overview
if(grepl('.rds$',input_overview)){
    sample_overview <- readRDS(input_overview)
    Allmutations <- read.delim(input_Allmutations) %>%
        left_join(sample_overview, by = c('patient_id' = 'cruk_id')) %>%
        tidyr::separate_rows(RegionSum, sep = ';') %>%
        mutate(Region = purrr::map_chr(RegionSum,~strsplit(.x,':')[[1]][1])) %>%
        select(patient_id,Region,histology_multi_full) %>%
        unique()
    
    Ntot <- data.frame(histology_multi_full = 'NSCLC',Ntot = nrow(Allmutations)) %>%
        rbind(Allmutations %>% group_by(histology_multi_full) %>%
              summarize(Ntot = n()) %>% filter(histology_multi_full %in% c('LUAD','LUSC'))) %>%
        rename(subtype = histology_multi_full)

}else{
    sample_overview <- readxl::read_xlsx(input_overview, sheet = 'TableS2', skip = 1)
    Allmutations <- readxl::read_xlsx(input_Allmutations, sheet = 'TableS3', skip = 19) %>%
        left_join(sample_overview, by = c('SampleID' = 'TRACERxID')) %>%
        tidyr::separate_rows(RegionSum, sep = ';') %>%
        mutate(Region = purrr::map_chr(RegionSum,~strsplit(.x,':')[[1]][1])) %>%
        select(SampleID,Region,Histology) %>%
        unique() %>%
        mutate(subtype = case_when(
                   Histology == 'Invasive adenocarcinoma' ~  'LUAD',
                   Histology == 'Squamous cell carcinoma' ~ 'LUSC',
                   TRUE ~  'Other'
                   ))
     Ntot <- data.frame(subtype = 'NSCLC',Ntot = nrow(Allmutations)) %>%
        rbind(Allmutations %>% group_by(subtype) %>%
              summarize(Ntot = n()) %>% filter(subtype %in% c('LUAD','LUSC')))       
}


# Fetch table with number of samples
Table_Nsamples <-
    data.frame(file = input_mutations) %>%
    mutate(panel = purrr::map_chr(file, ~strsplit(.x,'/')[[1]][3]),
           subtype = gsub('.txt','',purrr::map_chr(file, ~strsplit(.x,'_')[[1]][3])),

           # read file and get number of unique samples
           Nsample = purrr::map_int(
                                file, ~nrow(unique(select(read.delim(.x),SampleID,Region))))) %>%
    left_join(Ntot) %>%
    mutate(fraction_evaluable_tumors = Nsample/Ntot) %>%
    select(panel,subtype,Nsample,Ntot,fraction_evaluable_tumors)


#-------------------------------------------------------------------------------
# 2.1 Read data and count number of unique samples
#-------------------------------------------------------------------------------
write.table(Table_Nsamples, file = output, sep = '\t', row.names = F,quote = F  )
