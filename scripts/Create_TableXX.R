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
    input_mutations <- snakemake@input[["Mutations"]]
    input_overview <- snakemake@input[["SampleOverview"]]

    output <-  snakemake@output[["Table"]]
}else{
    input_mutations <-  c('data/KappaHyperExome/Selected_mutations_LUAD.txt','data/IlluminaFocusPanel/Selected_mutations_LUAD.txt','data/KappaHyperExome/Selected_mutations_LUSC.txt','data/IlluminaFocusPanel/Selected_mutations_LUSC.txt')
    input_overview <- 'data/TRACERx421_supplement/20221109_TRACERx421_all_patient_df.rds'
    output <- 'output/Tables/TableXX_NumberOfEvaluableSamples.txt'

}

#-------------------------------------------------------------------------------
# 1.1 Read data and count number of unique samples
#-------------------------------------------------------------------------------
# Read sample overview
sample_overview <- readRDS(input_overview)

# read mutations
Table_Nsamples <-
    data.frame(file = input_mutations) %>%
    mutate(panel = purrr::map_chr(file, ~strsplit(.x,'/')[[1]][2]),
           subtype = gsub('.txt','',purrr::map_chr(file, ~strsplit(.x,'_')[[1]][3])),

           # read file and get number of unique samples
           Nsample = purrr::map_int(
                                file, ~nrow(unique(select(read.delim(.x),SampleID,Region)))),
           Ntot = dplyr::case_when(
                             subtype == 'LUAD' ~ 363,
                             subtype == 'LUSC' ~ 208,
                             subtype == 'NSCLC' ~ 740),
           Nsample_pct = Nsample / Ntot * 100) %>%
    select(panel,subtype,Nsample,Nsample_pct) %>%
    arrange(Nsample)

#-------------------------------------------------------------------------------
# 2.1 Read data and count number of unique samples
#-------------------------------------------------------------------------------
write.table(Table_Nsamples, file = output, sep = '\t', row.names = F,quote = F  )
