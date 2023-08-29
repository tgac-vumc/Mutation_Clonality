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
    output <-  snakemake@output[["Metrics"]]
}else{
    input_mutations <-  c('data/KappaHyperExome/Selected_mutations.txt','data/IlluminaFocusPanel/Selected_mutations.txt','data/IlluminaComprehensiveCancerPanel/Selected_mutations.txt','data/IlluminaComprehensiveCancerPanelv3/Selected_mutations.txt','data/IlluminaTrueSightTumor170/Selected_mutations.txt','data/IlluminaTrueSightTumor500/Selected_mutations.txt','data/InhouseLungPanel/Selected_mutations.txt') 
    output <- 'output/Tables/TableXX_NumberOfEvaluableSamples.txt'

}

#-------------------------------------------------------------------------------
# 1.1 Read data and count number of unique samples
#-------------------------------------------------------------------------------
# read mutations
Table_Nsamples <-
    data.frame(file = input_mutations) %>%
    mutate(panel = purrr::map_chr(file, ~strsplit(.x,'/')[[1]][2]),
           # read file and get number of unique samples
           Nsample = purrr::map_int(
                                file, ~nrow(unique(select(read.delim(.x),SampleID,Region)))),
           Nsample_pct = Nsample / 327 * 100) %>%
    select(panel,Nsample,Nsample_pct) %>%
    arrange(Nsample)

#-------------------------------------------------------------------------------
# 2.1 Read data and count number of unique samples
#-------------------------------------------------------------------------------
write.table(Table_Nsamples, file = output, sep = '\t', row.names = F,quote = F  )
