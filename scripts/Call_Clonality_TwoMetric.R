#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Call_Clonality_TwoMetric.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Call clonality using one comparison for clonal and non-clonal pairs
# Scenario 1: Clonal pairs are made up out of primary and metastasis (WorstCaseScenario)
# Scenario 2: Clonal pairs are made up out of random primary and metastasis (RandomScenario
#
# For TRACERx100 take two tumors of primary under the same scenarios

# Author: Jurriaan Janssen (j.janssen4@amsterdamumc.nl)
#
# TODO:
# 1) 
#
# History:
#  12-12-2023: File creation, write code
#  19-08-2024: Changes for complete TRACERx421 and Tumor pairs
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
if(!'jaccard' %in% installed.packages()){install.packages('jaccard')}
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(jaccard)))
suppressMessages(suppressWarnings(library(Clonality)))
source('scripts/Calculate_TCGA_mutation_frequencies.R')

#-------------------------------------------------------------------------------
# 1.1 Parse snakemake objects
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_mutations <- snakemake@input[["Mutations"]]
    input_TumorPairs <- snakemake@input[["Tumor_pairs"]]
    input_metrics <- snakemake@input[["Metrics"]]
    subtype <- snakemake@wildcards[["subtype"]]
    output_shared_mutations <-snakemake@output[["Shared_Mutations"]]
    output <-  snakemake@output[["Metrics"]]
}else{
    input_mutations <-  'output/InhouseLungPanel/Selected_mutations_LUAD.txt'
    input_annotations <- 'output/OncoKB_annotations_LUAD.txt'
    input_TumorPairs <- 'data/Tumor_pairs_LUAD.txt'
    subtype <- 'LUAD'
    output <- 'output/InhouseLungPanel/Clonality_metrics_LUAD.txt'
    output_shared_mutations <- 'output/InhouseLungPanel/Shared_mutations_LUAD.txt'
}

#-------------------------------------------------------------------------------
# 1.1 Read data 
#-------------------------------------------------------------------------------
# read mutations
mutations <- read.delim(input_mutations, stringsAsFactors = F)
# read tumor pairs
Tumor_pairs <- read.delim(input_TumorPairs, stringsAsFactors = F)
    
#-------------------------------------------------------------------------------
# 2.1 Reformat mutations
#-------------------------------------------------------------------------------
mutations_broad <-
    # Create fields
    mutations %>%
    mutate(MutationID = paste(gsub('chr', '', chr),start, ref,var, sep = ' '),
           SampleID = paste0(SampleID,'_',Region),
           Mutation_present = as.integer(Mutation_present)) %>%
    select(SampleID,Hugo_Symbol,MutationID,Mutation_present) %>%
    unique() %>%
    # Get broad format and fill missing values with 0 (not present)
    tidyr::pivot_wider(id_cols = c(MutationID,Hugo_Symbol),
                       values_from = Mutation_present,
                       names_from = SampleID,
                       values_fill = 0)


# fill in missing samples with 0's
Samples <- unique(c(Tumor_pairs$Tumor1,Tumor_pairs$Tumor2))
Missing_Samples <- Samples[!Samples %in% colnames(mutations_broad)]
for(sample in Missing_Samples){
    mutations_broad[sample] = 0
}

# Create numeric matrix
mutation_matrix <-
    mutations_broad %>%
    select(-Hugo_Symbol) %>%
    tibble::column_to_rownames(var= 'MutationID') %>%
    as.matrix()

#-------------------------------------------------------------------------------
# 2.2 Fetch reference data for LR model
#-------------------------------------------------------------------------------
freq <- Calculate_TCGA_frequencies(rownames(mutation_matrix), subtype)
#-------------------------------------------------------------------------------
# 3.1 Perform Clonality testing
#-------------------------------------------------------------------------------
# Calculate LR ratio and Jaccard similarity statistics
#-------------------------------------------------------------------------------
# Initialize dataframes
Clonality_metrics <- data.frame()
Matching_mutations <- data.frame()                      
# Iterate over tumor pairs
for(i in 1:nrow(Tumor_pairs)){
    sample1 <- Tumor_pairs[i,'Tumor1']
    sample2 <- Tumor_pairs[i,'Tumor2']
     # Perform Likelihood test and calculate jaccard similarity
    set.seed(123)
    LR_output <- SNVtest(mutation_matrix[,sample1],mutation_matrix[,sample2], freq)
    Jaccard_similarity <- jaccard(mutation_matrix[,sample1], mutation_matrix[,sample2])
    Jaccard_pvalue <- jaccard.test(mutation_matrix[,sample1], mutation_matrix[,sample2], method = 'mca',accuracy=1e-05 )$pvalue

    # if both regions have no mutations
    if(sum(mutation_matrix[,sample1]) + sum(mutation_matrix[,sample2]) == 0){
        LR_output <- c(0,0,0,NA,NA,NA)
        names(LR_output) <- c('n1','n2','n_match','LRstat','maxKsi','LRpvalue')
        Jaccard_similarity <- 0
    }

    if(LR_output['n_match'] == 0){Jaccard_pvalue <- 1}
    
    
    # Add to df
    Clonality_metrics <- rbind(Clonality_metrics,data.frame(t(LR_output)) %>% mutate(comparison = paste0(sample1,'-',sample2), Jaccard = Jaccard_similarity,pval_Jaccard = Jaccard_pvalue, True_clonality = Tumor_pairs$True_Clonality[i]) )

    # Add shared mutation information to df
    shared_mutations <- rownames(mutation_matrix)[which(rowSums(mutation_matrix[,c(sample1,sample2)]) == 2)]
    Matching_mutations <-
        rbind(Matching_mutations,
              data.frame(Shared_mutations = shared_mutations) %>%
              left_join(mutations_broad[,c('MutationID','Hugo_Symbol')], by = c('Shared_mutations' = 'MutationID')) %>%
              mutate(Sample1 = sample1, Sample2 = sample2, True_clonality = Tumor_pairs$True_Clonality[i])
              )
    
}

# Create output df
Clonality_metrics_out <- cbind(Tumor_pairs, Clonality_metrics[,c('n1','n2','n_match','LRstat','LRpvalue','Jaccard','pval_Jaccard')])

#-------------------------------------------------------------------------------
# 3.1 Write to file
#-------------------------------------------------------------------------------
write.table(Clonality_metrics_out, file = output, sep = '\t', row.names = F,quote = F  )
write.table(Matching_mutations,file = output_shared_mutations, sep = '\t', row.names = F,quote = F  )
