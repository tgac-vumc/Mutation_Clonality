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
#  20-08-2024: Add OncoKB annotations, edits for tumor pairs
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
suppressMessages(suppressWarnings(library(dplyr)))
#-------------------------------------------------------------------------------
# 1.1 Parse snakemake objects
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_mutations <- snakemake@input[["Mutations"]]
    input_TumorPairs <- snakemake@input[["Tumor_pairs"]]
    input_Annotations <- snakemake@input[["Annotations"]]

    subtype <- snakemake@wildcards[["subtype"]]
    output <-  snakemake@output[["Metrics"]]
}else{
    input_mutations <-  'output/InhouseLungPanel/Selected_mutations_LUAD.txt'
    input_Annotations  <- 'output/OncoKB_annotations_LUAD.txt'
    input_TumorPairs <- 'data/Tumor_pairs_LUAD.txt'
    subtype <- 'LUAD'
    output <- 'output/KappaHyperExome/Clonality_metrics_NSCLC.txt'
}

#-------------------------------------------------------------------------------
# 1.1 Read data 
#-------------------------------------------------------------------------------
# read mutations
mutations <- read.delim(input_mutations, stringsAsFactors = F)
# read OncoKB annotations
Annotations <- read.delim(input_Annotations,stringsAsFactors = F) %>% 
    mutate(MutationID = paste0(Hugo_Symbol,'_',paste(gsub('chr', '', Chromosome),Start_Position, Reference_Allele,Tumor_Seq_Allele1, sep = ' '))) %>%
    select(-Tumor_Sample_Barcode) %>%
    unique()

# read tumor pairs
Tumor_pairs <- read.delim(input_TumorPairs, stringsAsFactors = F)
#-------------------------------------------------------------------------------
# 2.1 Reformat mutations
#-------------------------------------------------------------------------------
# Add fields and join annotations
mutations <- mutations %>% 
    mutate(MutationID = paste(gsub('chr', '', chr),start, ref,var, sep = ' '),
           SampleID = paste0(SampleID,'_',Region),
           Mutation_present = as.integer(Mutation_present),
           MutationID = paste0(Hugo_Symbol,'_',MutationID)) %>% 
    select(-c(AAChange,NucleotideChange)) %>%
    left_join(unique(Annotations[,c('MutationID','HGVSg','HGVSp','ONCOGENIC')])) 

# Fetch mutation info
INFO <-
    mutations %>% 
    mutate(INFO = paste0(MutationID,':',HGVSg,';', HGVSp,';VAF_',VAF),
           Pathogenic = dplyr::case_when(
                                   grepl('EGFR',INFO) & grepl('T790M',HGVSp) ~ FALSE, # TKI resistance mutation
                                   nchar(ref) > 2 | nchar(var) > 2 ~ TRUE, # indels
                                   !is.na(HGVSp) ~ TRUE, #nonsynonomous 
                                   TRUE ~ FALSE),
           # EGFR resistance mutations are oncogenic but should not be selected as driver mutations
           ONCOGENIC = case_when(
                   grepl('T790M|A767|M766',HGVSp) & Hugo_Symbol == 'EGFR' ~ 'Resistance Mutation',
                   TRUE ~ ONCOGENIC)) %>%
    select(SampleID,Hugo_Symbol,MutationID,INFO,Pathogenic,ONCOGENIC) %>%
    unique()


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


# fill in missing samples with 0's
Samples <- unique(c(Tumor_pairs$Tumor1,Tumor_pairs$Tumor2))
Missing_Samples <- Samples[!Samples %in% colnames(mutations_broad)]
for(sample in Missing_Samples){
    mutations_broad[sample] = 0
}

# fetch numeric matrix
mutation_matrix <-
    mutations_broad %>%
    tibble::column_to_rownames(var= 'MutationID') %>%
    as.matrix()

#-------------------------------------------------------------------------------
# 2.1 Define MC algorithm
#-------------------------------------------------------------------------------
# define drivers genes
oncogenic_driver_genes <- c("KRAS","EGFR","ERBB2", 'BRAF',"MET")
# Fetch base changes flagged as driver mutations
oncogenic_driver_mutations <-
    INFO %>% 
    filter((grepl('Oncogenic',ONCOGENIC)&(Hugo_Symbol %in% oncogenic_driver_genes))) %>%
    pull(MutationID) %>%
    unique()

MCalgorithm <- function(comparison) {
    # Fetch samplenames
    sample1 <- comparison[,'Tumor1']
    sample2 <- comparison[,'Tumor2']
    # Fetch mutations
    data1 <- mutation_matrix[,sample1]
    data2 <- mutation_matrix[,sample2]
   


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
                data.frame(Sample1=sample1,Sample2=sample2,True_clonality = comparison$True_Clonality, Clonality_MC = Clonality, reason = 'Different oncogenic driver mutations found',Shared_mutations = shared_mutations,Driver_information_sample1 =  INFO_oncogenic1, Driver_information_sample2 =  INFO_oncogenic2, Other_information_sample1 = INFO_no_drivers1, Other_information_sample2 = INFO_no_drivers2 )
            )
            
        }

    }
    
    # If same oncogenic mutation, or both wild type 
    if(identical(oncogenic1,oncogenic2)){
        # next check if there is more or equal than 1 shared mutation (excluding driver mutations)
        match_bool <- rowSums(matrix(c(data_no_drivers1,data_no_drivers2),ncol=2)) == 2
        n_match <- sum(match_bool)
        # If there is any match mutation it is clonal
        if(any(match_bool)){
            Clonality <- 'Clonal'
            if(all(oncogenic1 == 0) & all(oncogenic2 == 0)){
                reason <- paste0('Wildtype oncogenic but shared mutations' )
            }else{
                reason <- paste0('Same oncogenic driver and shared mutations')
            }
            return(
                data.frame(Sample1=sample1,Sample2=sample2,True_clonality = comparison$True_Clonality, Clonality_MC = Clonality, reason = reason,Shared_mutations = shared_mutations,Driver_information_sample1 =  INFO_oncogenic1, Driver_information_sample2 =  INFO_oncogenic2, Other_information_sample1 = INFO_no_drivers1, Other_information_sample2 = INFO_no_drivers2)
            )
        }else{
            # In the case of no matching mutations:
            TP53_1 <- data1[purrr::map_chr(names(data1),~strsplit(.x,'_')[[1]][1]) == 'TP53']
            TP53_2 <- data2[purrr::map_chr(names(data1),~strsplit(.x,'_')[[1]][1]) == 'TP53']
            if(all(TP53_1 == 0) & all(TP53_2 == 0)){
                Clonality <- 'Inconclusive'
                return(data.frame(Sample1=sample1,Sample2=sample2,True_clonality = comparison$True_Clonality, Clonality_MC = Clonality, reason = 'No shared non-oncogenic mutations found and TP53 wildtype',Shared_mutations = shared_mutations,Driver_information_sample1 =  INFO_oncogenic1, Driver_information_sample2 =  INFO_oncogenic2, Other_information_sample1 = INFO_no_drivers1, Other_information_sample2 = INFO_no_drivers2))
            }else{
                 Clonality <- 'Probably non-Clonal'
                 return(data.frame(Sample1=sample1,Sample2=sample2,True_clonality = comparison$True_Clonality, Clonality_MC =Clonality, reason = 'No shared non-oncogenic mutations found and different TP53 mutations',Shared_mutations = shared_mutations,Driver_information_sample1 =  INFO_oncogenic1, Driver_information_sample2 =  INFO_oncogenic2, Other_information_sample1 = INFO_no_drivers1, Other_information_sample2 = INFO_no_drivers2))
            }
        }
    }
    return(data.frame(Sample1=sample1,Sample2=sample2,True_clonality = comparison$True_Clonality, Clonality_MC = '????', reason = 'Unknown',Shared_mutations = shared_mutations,Driver_information_sample1 =  INFO_oncogenic1, Driver_information_sample2 =  INFO_oncogenic2, Other_information_sample1 = INFO_no_drivers1, Other_information_sample2 = INFO_no_drivers2))
}




#-------------------------------------------------------------------------------
# 3.1 Perform Clonality testing using clonal pairs
#-------------------------------------------------------------------------------
# Create clonal pairs with different regions within one patient
#-------------------------------------------------------------------------------
# Initialize dataframes
MC_Clonalities <- data.frame()

# Iterate over comparisons
for(i in 1:nrow(Tumor_pairs)){
    # Perform molecular classification algorithm
    MC_Clonalities <-rbind(MC_Clonalities, MCalgorithm(Tumor_pairs[i,]))
}

#-------------------------------------------------------------------------------
# 3.1 Write to file
#-------------------------------------------------------------------------------
write.table( MC_Clonalities, file = output, sep = '\t', row.names = F,quote = F )
