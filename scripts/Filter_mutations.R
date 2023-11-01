#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Filter_mutations.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Filter TRACERx mutations by panel regions
#
# Author: Jurriaan Janssen (j.janssen4@amsterdamumc.nl)
#
# TODO:
# 1) 
#
# History:
#  22-08-2023: File creation, write code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
suppressMessages(suppressWarnings(library(dplyr)))

#-------------------------------------------------------------------------------
# 1.1 Parse snakemake objects
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_mutations <- snakemake@input[["Mutations"]]
    input_panel <- snakemake@input[["panel"]]
    panelID <- snakemake@wildcards[["panel"]]
    subtype <- snakemake@wildcards[["subtype"]]
    output <-  snakemake@output[["Mutations"]]
}else{
    input_mutations <- 'data/TRACERx_supplement/nejmoa1616288_appendix_2.xlsx'
    input_panel <- 'manifest/TP53.bed'
    panelID <- 'TP53'
    output <- 'data/InhouseLungPanel/Selected_mutations_LUAD.txt'
    subtype <- 'LUAD'
}

#-------------------------------------------------------------------------------
# 1.1 Read data 
#-------------------------------------------------------------------------------
# read mutations object
mutations <- readxl::read_xlsx(input_mutations, sheet = 'TableS3', skip = 19)
# Read sample info
Sample_info <- readxl::read_xlsx(input_mutations, sheet = 'TableS2', skip = 1)
# read panel
panel <- read.delim(input_panel, header = F)
#-------------------------------------------------------------------------------
# 2.1 Reformat mutations
#-------------------------------------------------------------------------------
# Select samples by subtype
if(subtype == 'LUAD'){
    Sample_info <- Sample_info %>% filter(Histology == 'Invasive adenocarcinoma')
}else if(subtype == 'LUSC'){
    Sample_info <- Sample_info %>% filter(Histology == 'Squamous cell carcinoma')
}

samples <- unique(Sample_info$TRACERxID)

# Transform mutation table to indicate presence of mutation in recurrence
mutations <- mutations %>%
    # filter by subtype
    filter(SampleID %in% samples) %>%
    # split rows per recurrence
    tidyr::separate_rows(RegionSum, sep = ';') %>%
    # recover number of alt reads and region number
    mutate(Region = purrr::map_chr(RegionSum,~strsplit(.x,':')[[1]][1]),
           ALT =  purrr::map_chr(RegionSum,~strsplit(strsplit(.x,':')[[1]][2],'/')[[1]][1]),
           chr= paste0('chr',chr),
           # create boolean to indicate presence of mutation
           Mutation_present = ALT != 0) %>% 
    select(SampleID,Region,Hugo_Symbol,chr,start,stop,ref,var,func,MutationID,Mutation_present)

#-------------------------------------------------------------------------------
# 2.2 Filter mutations by panel
#-------------------------------------------------------------------------------
colnames(panel)[1:3] <- c('chr','chromStart','chromEnd')

if(panelID == 'KappaHyperExome'){
    Filtered_mutations <-
        mutations %>%
        filter(func == 'exonic')
}else if(panelID %in% c('TP53','KRAS','EGFR','MET')){
    Filtered_mutations <-
        mutations %>%
        filter(Hugo_Symbol == panelID)
}else{
    # Join panel regions and filter by panel
    Filtered_mutations <-
        mutations  %>%
        inner_join(panel) %>%
        filter(start >= chromStart & start <= chromEnd ) %>%
        select(-chromStart,-chromEnd)
}

#-------------------------------------------------------------------------------
# 3.1 Write to file
#-------------------------------------------------------------------------------
write.table(Filtered_mutations, file = output, sep = '\t', row.names = F,quote = F  )
