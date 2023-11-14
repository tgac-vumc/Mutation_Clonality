#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Call_Clonality_MCalgorithm.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Call clonality in clonal and non-clonal pairs using molecular classification
# (MC) algorithm
#
# Author: Jurriaan Janssen (j.janssen4@amsterdamumc.nl), Barbara Andrade Barbosa
#
# TODO:
# 1) 
#
# History:
#  09-11-2023: File creation, write code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
suppressMessages(suppressWarnings(library(dplyr)))
#-------------------------------------------------------------------------------
# 1.1 Parse snakemake objects
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_mutations <- snakemake@input[["Mutations"]]
    input_mutations_inhouse <- snakemake@input[["Mutations_inhouse"]]

    subtype <- snakemake@wildcards[["subtype"]]
    output <-  snakemake@output[["Metrics"]]
}else{
    input_mutations <-  'data/InhouseLungPanel/Selected_mutations_NSCLC.txt'
    subtype <- 'NSCLC'
    input_mutations_inhouse <- 'data/InhouseLungPanel/Selected_mutations_NSCLC.txt'
    output <- 'output/KappaHyperExome/Clonality_metrics_NSCLC.txt'
}

#-------------------------------------------------------------------------------
# 1.1 Read data 
#-------------------------------------------------------------------------------
# read mutations
mutations <- read.delim(input_mutations, stringsAsFactors = F)

# read Inhouse mutation data
Inhouse_Mutations <- read.delim(input_mutations_inhouse)

#-------------------------------------------------------------------------------
# 2.1 Reformat mutations
#-------------------------------------------------------------------------------
# Fetch SampleNames
SampleNames <- unique(mutations$SampleID)

# Add fields
mutations <- mutations %>%
    mutate(MutationID = paste(gsub('chr', '', chr),start, ref,var, sep = ' '),
           SampleID = paste0(SampleID,'_',Region),
           Mutation_present = as.integer(Mutation_present),
           MutationID = paste0(Hugo_Symbol,'_',MutationID))

mutations_broad <-
    # Create fields
    mutations %>%
    select(SampleID,MutationID,Mutation_present) %>%
    unique() %>%
    # Get broad format and fill missing values with 0 (not present)
    tidyr::pivot_wider(id_cols =c(MutationID),
                       values_from = Mutation_present,
                       names_from = SampleID,
                       values_fill = 0)

# fetch numeric matrix
mutation_matrix <-
    mutations_broad %>%
    tibble::column_to_rownames(var= 'MutationID') %>%
    as.matrix()

# Fetch mutation info
INFO <-
    mutations %>% 
    mutate(INFO = paste0(MutationID,':',NucleotideChange,';', AAChange,';VAF_',VAF)) %>%
    select(SampleID,MutationID,INFO)

#-------------------------------------------------------------------------------
# 2.1 Define MC algorithm
#-------------------------------------------------------------------------------
# define drivers genes
oncogenic_driver_genes <- c("KRAS","EGFR","ERBB2", "MET")
# Fetch base changes
oncogenic_driver_mutations <- Inhouse_Mutations %>%
    mutate(
        MutationID = paste(gsub('chr', '', chr),start, ref,var, sep = ' '),
        MutationID = paste0(Hugo_Symbol,'_',MutationID)) %>%
    filter(Hugo_Symbol %in% oncogenic_driver_genes) %>%
    pull(MutationID)

MCalgorithm <- function(comparison,True_clonality) {
    # Fetch samplenames
    sample1 <- comparison[,1]
    sample2 <- comparison[,2]
    # Fetch data based on True clonality comparison
    if(True_clonality == 'Clonal'){
        data1 <- mutation_matrix[,sample1]
        data2 <- mutation_matrix[,sample2]
    }else{
        data1 <- mutation_matrix[,grepl(sample1,colnames(mutation_matrix))][,1]
        data2 <- mutation_matrix[,grepl(sample2,colnames(mutation_matrix))][,1]
    }


    #Fetch oncogenic driver mutations
    #oncogenic1 <- data1[purrr::map_chr(names(data1),~strsplit(.x,'_')[[1]][1]) %in% oncogenic_driver_mutations]
    #oncogenic2 <- data2[purrr::map_chr(names(data1),~strsplit(.x,'_')[[1]][1]) %in% oncogenic_driver_mutations]

    oncogenic1 <- data1[names(data1) %in% oncogenic_driver_mutations]
    oncogenic2 <- data2[names(data2) %in% oncogenic_driver_mutations]

    # Fetch other mutations
    data_no_drivers1 <- data1[!names(data1) %in% names(oncogenic1)]
    data_no_drivers2 <- data2[!names(data2) %in% names(oncogenic2)]

    INFO_oncogenic1 <- paste(INFO[INFO$SampleID == sample1 & INFO$MutationID %in% names(oncogenic1[oncogenic1 == 1]),'INFO'],collapse=' --- ')
    INFO_oncogenic2 <- paste(INFO[INFO$SampleID == sample2 & INFO$MutationID %in% names(oncogenic1[oncogenic2 == 1]),'INFO'],collapse=' --- ')
    INFO_no_drivers1 <- paste(INFO[INFO$SampleID == sample1 & INFO$MutationID %in% names(data_no_drivers1[data_no_drivers1 == 1]),'INFO'],collapse=' --- ')
    INFO_no_drivers2 <- paste(INFO[INFO$SampleID == sample2 & INFO$MutationID %in% names(data_no_drivers2[data_no_drivers2 == 1]),'INFO'],collapse=' --- ')

    shared_mutations <- paste(names(data1)[which(data1 == 1 & data2 == 1)],collapse=' --- ')
    
    #-------------------------------------------------------------------------------
                                        # Call clonality
    #-------------------------------------------------------------------------------
    # If any oncogenic mutation is found in any sample
    if(any(oncogenic1 == 1) | any(oncogenic2 == 1)){
        # if samples do not have identical driver mutations it is non-clonal
        if(!identical(oncogenic1,oncogenic2)){
            Clonality <- 'non-Clonal'
            return(
                data.frame(Sample1=sample1,Sample2=sample2,True_clonality = True_clonality, Clonality_MC = Clonality, reason = 'Different oncogenic driver mutations found',Shared_mutations = shared_mutations,Driver_information_sample1 =  INFO_oncogenic1, Driver_information_sample2 =  INFO_oncogenic2, Other_information_sample1 = INFO_no_drivers1, Other_information_sample2 = INFO_no_drivers2 )
            )
            
        }

    }
    
    # If same oncogenic mutation, or both wild type 
    if(identical(oncogenic1,oncogenic2)){
        # next check if there is more or equal than 1 shared mutation (excluding driver mutations
        match_bool <- rowSums(matrix(c(data_no_drivers1,data_no_drivers2),ncol=2)) == 2
        n_match <- sum(match_bool)
        # If there is any match mutation it is clonal
        if(any(match_bool)){
            Clonality <- 'Clonal'
            if(all(oncogenic1 == 0) & all(oncogenic2 == 0)){
                reason <- paste0('Wildtype oncogenic but ',n_match,' shared mutations' )
            }else{
                reason <- paste0('Same oncogenic driver and ',n_match,' shared mutations')
            }
            return(
                data.frame(Sample1=sample1,Sample2=sample2,True_clonality = True_clonality, Clonality_MC = Clonality, reason = reason,Shared_mutations = shared_mutations,Driver_information_sample1 =  INFO_oncogenic1, Driver_information_sample2 =  INFO_oncogenic2, Other_information_sample1 = INFO_no_drivers1, Other_information_sample2 = INFO_no_drivers2)
            )
        }else{
            # In the case of no matching mutations:
            TP53_1 <- data1[purrr::map_chr(names(data1),~strsplit(.x,'_')[[1]][1]) == 'TP53']
            TP53_2 <- data2[purrr::map_chr(names(data1),~strsplit(.x,'_')[[1]][1]) == 'TP53']
            if(all(TP53_1 == 0) & all(TP53_2 == 0)){
                Clonality <- 'Inconclusive'
                return(data.frame(Sample1=sample1,Sample2=sample2,True_clonality = True_clonality, Clonality_MC = Clonality, reason = 'No shared non-oncogenic mutations found and TP53 wildtype',Shared_mutations = shared_mutations,Driver_information_sample1 =  INFO_oncogenic1, Driver_information_sample2 =  INFO_oncogenic2, Other_information_sample1 = INFO_no_drivers1, Other_information_sample2 = INFO_no_drivers2))
            }else{
                 Clonality <- 'Probably non-Clonal'
                 return(data.frame(Sample1=sample1,Sample2=sample2,True_clonality = True_clonality, Clonality_MC =Clonality, reason = 'No shared non-oncogenic mutations found and different TP53 mutations',Shared_mutations = shared_mutations,Driver_information_sample1 =  INFO_oncogenic1, Driver_information_sample2 =  INFO_oncogenic2, Other_information_sample1 = INFO_no_drivers1, Other_information_sample2 = INFO_no_drivers2))
            }
        }
    }
    return(data.frame(Sample1=sample1,Sample2=sample2,True_clonality = True_clonality, Clonality_MC = '????', reason = 'Unknown',Shared_mutations = shared_mutations,Driver_information_sample1 =  INFO_oncogenic1, Driver_information_sample2 =  INFO_oncogenic2, Other_information_sample1 = INFO_no_drivers1, Other_information_sample2 = INFO_no_drivers2))
}




#-------------------------------------------------------------------------------
# 3.1 Perform Clonality testing using clonal pairs
#-------------------------------------------------------------------------------
# Create clonal pairs with different regions within one patient
#-------------------------------------------------------------------------------
# Initialize dataframes
MC_Clonalities <- data.frame()

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
        # Perform molecular classification algorithm
        MC_Clonalities <-
            rbind(
            MC_Clonalities,
            MCalgorithm(comparisons[i,], True_clonality = 'Clonal')
            )
        
    }
}


# Example of Clonal relation with different driver mutations
sample1 <- 'CRUK0062_R1'
sample2 <- 'CRUK0062_R2'
# Here one of the sample has an ALK mutation

sample1 <- 'CRUK0071_R1'
sample2 <- 'CRUK0071_R3'
# One sample has EGFR mutation and the other doesnt




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
    # Perform molecular classification algorithm
    MC_Clonalities <-
        rbind(
            MC_Clonalities,
            MCalgorithm(comparisons[i,], True_clonality = 'Non-Clonal')
        )
}

#-------------------------------------------------------------------------------
# 3.1 Write to file
#-------------------------------------------------------------------------------
write.table( MC_Clonalities, file = output, sep = '\t', row.names = F,quote = F )
