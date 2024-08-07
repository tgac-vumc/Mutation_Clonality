#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Call_Clonality_OneComparison.R
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
    input_mutationsExome <- snakemake@input[["MutationsExome"]]

    input_metrics <- snakemake@input[["Metrics"]]
    sample_overview <- snakemake@input[["sampleOverview"]]
    subtype <- snakemake@wildcards[["subtype"]]
    scenario <- snakemake@wildcards[["scenario"]]
    dataset <- snakemake@wildcards[["dataset"]]
    output_shared_mutations <-snakemake@output[["Shared_Mutations"]]
    output <-  snakemake@output[["Metrics"]]
}else{
    input_mutations <-  'data/TRACERx421/InhouseLungPanel/Selected_mutations_LUAD.txt'
    input_mutationsExome <-  'data/TRACERx421/KappaHyperExome/Selected_mutations_LUAD.txt'
    subtype <- 'LUAD'
    scenario <- 'WorseCaseScenario'
    dataset <- 'TRACERx421'
    sample_overview <- 'data/TRACERx421_supplement/sampleOverview.txt'
    output <- 'output/InhouseLungPanel/Clonality_metrics_LUAD.txt'
    output_shared_mutations <- 'output/InhouseLungPanel/Shared_mutations_LUAD.txt'
}

#-------------------------------------------------------------------------------
# 1.1 Read data 
#-------------------------------------------------------------------------------
# read exome
mutations <- read.delim(input_mutations, stringsAsFactors = F)
mutationsExome <- read.delim(input_mutationsExome, stringsAsFactors = F)

# read sample Overview
sample_overview <- read.delim(sample_overview) %>%
    select(region,sampleType)

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


mutations_broadExome <-
    # Create fields
    mutationsExome %>%
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

mutation_matrixExome <-
    mutations_broadExome %>%
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
# Create clonal pairs with different region within one patient
# Either create random pairs, or create pairs with lowest Jaccard index in exome
# For TRACERx421, select comparisons with primary and metastasis
#-------------------------------------------------------------------------------
# Initialize dataframes
Clonality_metrics <- data.frame()
Matching_mutations <- data.frame()                                     
# Iterate over samples
for(sample in SampleNames){
    jaccard_values <- c()
    mutation_matrix_sample <- mutation_matrix[,grepl(sample,colnames(mutation_matrix))]
    # Get unique region combinations without doublets
    comparisons <- data.table::setDT(data.table::CJ(colnames(mutation_matrix_sample),colnames(mutation_matrix_sample), unique=T)) %>%
        .[!duplicated(t(apply(., 1, sort))),] %>%
        as.data.frame() %>%
        filter(V1 != V2)
 
    if(scenario != 'All'){
    # For TRACERx421, only select comparisons between primary and metastasis
        if(dataset == 'TRACERx421'){
            comparisons <-
                comparisons %>%
                # Join sample type information
                left_join(sample_overview, by = c('V1' = 'region')) %>%
                left_join(sample_overview, by = c('V2' = 'region')) %>%
                filter(
                (sampleType.x == 'primary' & sampleType.y ==  'metastasis') |
                (sampleType.y == 'primary' & sampleType.x ==  'metastasis'))
            # move to next sample if no primary and metastasis is present in this patient
            #if(nrow(comparisons == 0)){next}
    }
        # Worst case is minimum jaccard in panel
        if(scenario  ==  'WorstCaseScenario'){
            for(i in 1:nrow(comparisons)){
                jaccard_values <- c(jaccard_values, jaccard(mutation_matrix[,comparisons$V1[i]],mutation_matrix[,comparisons$V2[i]]))
            }
            comparisons <- comparisons[which.min(jaccard_values),]
        # Worse case is minimum jaccard in exome
        }else if(scenario == 'WorseCaseScenario'){
            for(i in 1:nrow(comparisons)){
                jaccard_values <- c(jaccard_values, jaccard(mutation_matrixExome[,comparisons$V1[i]],mutation_matrixExome[,comparisons$V2[i]]))
            }
            comparisons <- comparisons[which.min(jaccard_values),]

    }else if(scenario == 'RandomScenario'){
        # Fetch random comparison
        set.seed(123)
        comparisons <- comparisons[sample(1:nrow(comparisons),1),]
    }
    }
    for(i in 1:nrow(comparisons)){
        sample1 <- comparisons[i,1]
        sample2 <- comparisons[i,2]

        # Perform Likelihood test and calculate jaccard similarity
        set.seed(123)
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
# Create non-clonal pairs with other one patient using  one region for each patient
# Do this either randomly, or the two samples with the highest jaccard index in exome
#-------------------------------------------------------------------------------
# Get all sample comparisons
for(sample in SampleNames){
    jaccard_values <- c()
    mutation_matrix_sample <- mutation_matrix[,grepl(sample,colnames(mutation_matrix))]
    mutation_matrix_Nosample <- mutation_matrix[,!grepl(sample,colnames(mutation_matrix))]
    # Fetch all possible combinations
    comparisons <- data.table::setDT(data.table::CJ(colnames(mutation_matrix_sample),colnames(mutation_matrix_Nosample), unique=T)) %>%
        .[!duplicated(t(apply(., 1, sort))),] %>%
        as.data.frame() %>%
        filter(V1 != V2)
    if(scenario != 'All'){
        if(scenario  ==  'WorstCaseScenario'){
            for(i in 1:nrow(comparisons)){
                jaccard_values <- c(jaccard_values,jaccard(mutation_matrix[,comparisons$V1[i]],mutation_matrix[,comparisons$V2[i]]))
            }
            comparisons <- comparisons[which.max(jaccard_values),]
        
        }else if(scenario == 'WorseCaseScenario'){
            for(i in 1:nrow(comparisons)){
                jaccard_values <- c(jaccard_values,jaccard(mutation_matrixExome[,comparisons$V1[i]],mutation_matrixExome[,comparisons$V2[i]]))
            }
            comparisons <- comparisons[which.max(jaccard_values),]

        }else if(scenario == 'RandomScenario'){
            # Fetch random comparison
            set.seed(123)
            comparisons <- comparisons[sample(1:nrow(comparisons),1),]
        }
    }
    for(i in 1:nrow(comparisons)){
        sample1 <- comparisons[i,1]
        sample2 <- comparisons[i,2]
        # Fetch patients andselect first column (Region)
        mutation_matrix_sample1 <- mutation_matrix[,sample1]
        mutation_matrix_sample2 <- mutation_matrix[,sample2]
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
