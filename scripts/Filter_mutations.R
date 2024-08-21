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
    input_mutations_Frankell <- snakemake@input[["Mutations_Frankell"]]
    input_mutations_Bakir <- snakemake@input[["Mutations_Bakir"]]
    input_sample_overview <- snakemake@input[["SampleOverview"]]
    input_panel <- snakemake@input[["panel"]]
    panelID <- snakemake@wildcards[["panel"]]
    subtype <- snakemake@wildcards[["subtype"]]
    output <-  snakemake@output[["Mutations"]]
}else{
    #input_mutations <- 'data/TRACERx_supplement/nejmoa1616288_appendix_2.xlsx'
    input_mutations_Frankell <- 'data/TRACERx421_supplement_Frankell/20221109_TRACERx421_mutation_table.fst'
    input_mutations_Bakir <- 'data/TRACERx421_supplement_Bakir/mutTableAll.cloneInfo.20220726.txt'

    input_sample_overview <- 'data/TRACERx421_supplement_Frankell/20221109_TRACERx421_all_patient_df.rds'
    input_panel <- 'manifest/InhouseLungPanel.bed'
    panelID <- 'InhouseLungPanel'
    output <- 'data/InhouseLungPanel/Selected_mutations_NSCLC.txt'
    subtype <- ''
}

#-------------------------------------------------------------------------------
# 1.1 Read data 
#-------------------------------------------------------------------------------
# read mutations object
mutations_Frankell <- fst::read_fst(input_mutations_Frankell)
mutations_Bakir <- read.delim(input_mutations_Bakir)
# Read sample info
Sample_info <- readRDS(input_sample_overview)
# read panel
panel <- read.delim(input_panel, header = F)

#-------------------------------------------------------------------------------
# 2.1 Reformat mutations
#-------------------------------------------------------------------------------
# Select relevant column


head(mutations_Frankell)
mutations_Frankell <- mutations_Frankell %>% select(patient_id,mutation_id,Hugo_Symbol,chr,start,stop,ref,var,func,exonic.func,NucleotideChange, AAChange,RegionSum) 
mutations_Bakir <- mutations_Bakir %>% select(patient_id,mutation_id,Hugo_Symbol,chr,start,stop,ref,var,func,exonic.func,NucleotideChange, AAChange,RegionSum)


# Select samples by subtype
if(subtype == 'LUAD'){
    Sample_info <- Sample_info %>% filter(grepl('LUAD',histology_multi_full))
}else if(subtype == 'LUSC'){
    Sample_info <- Sample_info %>% filter(grepl('LUSC',histology_multi_full))
}

samples <- unique(Sample_info$cruk_id)

# Transform mutation table to indicate presence of mutation in recurrence
# Combine Frankell and Bakir
mutations <- rbind(mutations_Frankell %>%filter(patient_id %in% samples) %>% tidyr::separate_rows(RegionSum, sep = ';'),
                   mutations_Bakir %>%filter(patient_id %in% samples) %>% tidyr::separate_rows(RegionSum, sep = ';')) %>%
    # recover number of alt reads and region number
    mutate(Region = purrr::map_chr(RegionSum,~strsplit(.x,':')[[1]][1]),
           ALT =  purrr::map_chr(RegionSum,~strsplit(strsplit(.x,':')[[1]][2],'/')[[1]][1]),
           Depth = purrr::map_chr(RegionSum,~strsplit(strsplit(.x,':')[[1]][2],'/')[[1]][2]), 
           VAF = as.numeric(ALT) / as.numeric(Depth),
           # create boolean to indicate presence of mutation
           Mutation_present = ALT != 0,
           SampleID = patient_id,
           MutationID = mutation_id) %>%
    select(SampleID,Region,Hugo_Symbol,chr,start,stop,ref,var,func,exonic.func,MutationID,NucleotideChange, AAChange,VAF,Mutation_present) %>%
    # Filter out two samples from different patients wich show 33 matching mutations
    #filter(!(SampleID %in% c('CRUK0036','CRUK0296'))) %>%
    unique() %>%
    arrange(MutationID)

#-------------------------------------------------------------------------------
# 2.2 Filter mutations by panel
#-------------------------------------------------------------------------------
colnames(panel)[1:3] <- c('chr','chromStart','chromEnd')

if(panelID == 'KappaHyperExome'){
    Filtered_mutations <-
        mutations 
}else if(panelID %in% c('TP53','KRAS','EGFR','MET')){
    Filtered_mutations <-
        mutations %>%
        filter(Hugo_Symbol == panelID)
}else if(panelID == 'FoundationOneCDx'){
    # Load FoundationOne panel (ExonGenes+panel objects)
    source('scripts/Fetch_FoundationOne_regions.R')
    Filtered_mutations <-
        mutations  %>%
        inner_join(panel) %>%
        filter((start >= chromStart & start <= chromEnd) | (func == 'exonic' & Hugo_Symbol %in% ExonGenes))
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

