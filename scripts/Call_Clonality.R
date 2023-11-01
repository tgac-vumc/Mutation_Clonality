#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Call_Clonality.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Call clonality using clonal and non-clonal pairs
#
# Author: Jurriaan Janssen (j.janssen4@amsterdamumc.nl)
#
# TODO:
# 1) 
#
# History:
#  25-08-2023: File creation, write code
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
    subtype <- snakemake@wildcards[["subtype"]]
    output_shared_mutations <-snakemake@output[["Shared_Mutations"]]
    output <-  snakemake@output[["Metrics"]]
}else{
    input_mutations <-  'data/KappaHyperExome/Selected_mutations_NSCLC.txt'
    subtype <- 'NSCLC'
    output <- 'output/KappaHyperExome/Clonality_metrics_LUAD.txt'
    output_shared_mutations <- 'output/KappaHyperExome/Shared_mutations_LUAD.txt'

}

#-------------------------------------------------------------------------------
# 1.1 Read data 
#-------------------------------------------------------------------------------
# read panel
mutations <- read.delim(input_mutations, stringsAsFactors = F)
#-------------------------------------------------------------------------------
# 2.1 Reformat mutations
#-------------------------------------------------------------------------------
# Fetch SampleNames
SampleNames <- unique(mutations$SampleID)

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
# 3.1 Perform Clonality testing using clonal pairs
#-------------------------------------------------------------------------------
# Create clonal pairs with different regions within one patient
#-------------------------------------------------------------------------------
# Initialize dataframes
Clonality_metrics <- data.frame()
Matching_mutations <- data.frame()                                     

# Iterate over samples
for(sample in SampleNames){
    mutation_matrix_sample <- mutation_matrix[,grepl(sample,colnames(mutation_matrix))]
    # Get unique region combinations without doublets
    comparisons <- data.table::setDT(data.table::CJ(colnames(mutation_matrix_sample),colnames(mutation_matrix_sample), unique=T)) %>%
        .[!duplicated(t(apply(., 1, sort))),] %>%
        as.data.frame() %>%
        filter(V1 != V2)
    # Iterate over comparisons
    for(i in 1:nrow(comparisons)){
        sample1 <- comparisons[i,1]
        sample2 <- comparisons[i,2]
        # Perform Likelihood test and calculate jaccard similarity
        set.seed(1)
        # skip if both regions have no mutations
        if(sum(mutation_matrix_sample[,sample1]) + sum(mutation_matrix_sample[,sample2]) == 0){next}
        LR_output <- SNVtest(mutation_matrix_sample[,sample1],mutation_matrix_sample[,sample2], freq)
        Jaccard_similarity <- jaccard(mutation_matrix_sample[,sample1], mutation_matrix_sample[,sample2])
        Jaccard_pvalue <- jaccard.test(mutation_matrix_sample[,sample1], mutation_matrix_sample[,sample2], method = 'mca',accuracy=1e-05 )$pvalue
        
        # Add to df
        Clonality_metrics <- rbind(Clonality_metrics,data.frame(t(LR_output)) %>% mutate(comparison = paste0(sample1,'-',sample2), Jaccard = Jaccard_similarity,pval_Jaccard = Jaccard_pvalue, True_clonality = 'Clonal') )

        # Add shared mutation information to df
        shared_mutations <- rownames(mutation_matrix_sample)[which(rowSums(mutation_matrix_sample) == 2)]
        Matching_mutations <-
            rbind(Matching_mutations,
              data.frame(Shared_mutations = shared_mutations) %>%
              left_join(mutations_broad[,c('MutationID','Hugo_Symbol')], by = c('Shared_mutations' = 'MutationID')) %>%
              mutate(Sample1 = sample1, Sample2 = sample2, True_clonality = 'Clonal')
              )
        
    }
}
#-------------------------------------------------------------------------------
# 3.1 Perform Clonality testing using non-clonal pairs
#-------------------------------------------------------------------------------
# Create non-clonal pairs with other patients using only one region (the first)
# for each patient
#-------------------------------------------------------------------------------
# Get unique sample combinations without doublets
comparisons <- data.table::setDT(data.table::CJ(V1=SampleNames,V2=SampleNames, unique=T)) %>%
    .[!duplicated(t(apply(., 1, sort))),] %>%
    as.data.frame() %>%
    filter(V1 != V2)

for(i in 1:nrow(comparisons)){
    sample1 <- comparisons[i,1]
    sample2 <- comparisons[i,2]
    # Fetch patients andselect first column (Region)
    mutation_matrix_sample1 <- mutation_matrix[,grepl(comparisons[i,1],colnames(mutation_matrix))][,1]
    mutation_matrix_sample2 <- mutation_matrix[,grepl(comparisons[i,2],colnames(mutation_matrix))][,1]
    # skip if both regions have no mutations
    if(sum(mutation_matrix_sample1) + sum(mutation_matrix_sample2) == 0){next}
    # perform Likelihood test and calculate jaccard similarity
    set.seed(123)
    LR_output <- SNVtest(mutation_matrix_sample1,mutation_matrix_sample2, freq)
    Jaccard_similarity <- jaccard(mutation_matrix_sample1, mutation_matrix_sample2)
    Jaccard_pvalue <- jaccard.test(mutation_matrix_sample1, mutation_matrix_sample2, method = 'mca',accuracy=1e-05)$pvalue
    # convert 0 match p values in jaccard to 1 instead of 0
    if(LR_output['n_match'] == 0){Jaccard_pvalue <- 1}

    # Add to df
    Clonality_metrics <- rbind(Clonality_metrics,data.frame(t(LR_output)) %>% mutate(comparison = paste0(comparisons[i,1],'-',comparisons[i,2]),  Jaccard = Jaccard_similarity,pval_Jaccard = Jaccard_pvalue, True_clonality = 'Non-Clonal') )


    # skip samples that have no matching mutations

    mutation_matrix_sample <- cbind(mutation_matrix_sample1,mutation_matrix_sample2)
    shared_mutations <- rownames(mutation_matrix_sample)[which(rowSums(mutation_matrix_sample) == 2)]
    if(length(shared_mutations) == 0){next}
        Matching_mutations <-
            rbind(Matching_mutations,
              data.frame(Shared_mutations = shared_mutations) %>%
              left_join(mutations_broad[,c('MutationID','Hugo_Symbol')], by = c('Shared_mutations' = 'MutationID')) %>%
              mutate(Sample1 = sample1, Sample2 = sample2, True_clonality = 'Non-Clonal')
              )
    
}

#-------------------------------------------------------------------------------
# 3.1 Write to file
#-------------------------------------------------------------------------------
write.table(Clonality_metrics, file = output, sep = '\t', row.names = F,quote = F  )
write.table(Matching_mutations,file = output_shared_mutations, sep = '\t', row.names = F,quote = F  )
