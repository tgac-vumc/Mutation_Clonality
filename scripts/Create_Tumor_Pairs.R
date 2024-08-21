#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create_Tumor_Pairs.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Create clonal and non-clonal tumor pairs of TRACERx421 tumors.
#
# Author: Jurriaan Janssen (j.janssen4@amsterdamumc.nl)
#
# TODO:
# 1) 
#
# History:
#  13-08-2023: File creation, write code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
suppressMessages(suppressWarnings(library(dplyr)))
#-------------------------------------------------------------------------------
# 1.1 Parse snakemake objects
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_mutations <- snakemake@input[["Mutations"]]
    input_patient <- snakemake@input[["PatientOverview"]]
    input_sample <- snakemake@input[["SampleOverview"]]
    input_sample_Bakir <- snakemake@input[["SampleOverview_Bakir"]]

    subtype <- snakemake@wildcards[["subtype"]]
    output <-  snakemake@output[["Tumor_pairs"]]
}else{
    input_mutations <- 'output/KappaHyperExome/Selected_mutations_NSCLC.txt'
    input_patient <-  'data/TRACERx421_supplement_Frankell/20221109_TRACERx421_all_patient_df.rds'
    input_sample <- 'data/TRACERx421_supplement_Frankell/20221109_TRACERx421_all_tumour_df.rds'
    input_sample_Bakir <- 'data/TRACERx421_supplement_Bakir/sampleOverview.txt'
    subtype <- 'LUSC'
}

#-------------------------------------------------------------------------------
# 1.1 Read data 
#-------------------------------------------------------------------------------
# read sample overviews
PatientOverview <- readRDS(input_patient)
SampleOverview <- readRDS(input_sample)
SampleOverview_Bakir <- read.delim(input_sample_Bakir)


# read mutations
Mutations <- read.delim(input_mutations) %>%
    select(SampleID,Region) %>% unique() 

#-------------------------------------------------------------------------------
# 2.1 Fetch tumor information
#-------------------------------------------------------------------------------
# Fetch SampleNames
# Select samples by subtype
if(subtype == 'LUAD'){
    PatientOverview <- PatientOverview %>% filter(grepl('LUAD',histology_multi_full)) 
}else if(subtype == 'LUSC'){
    PatientOverview <- PatientOverview %>% filter(grepl('LUSC',histology_multi_full))
}
Mutations <- Mutations %>% filter(SampleID %in% PatientOverview$cruk_id)
Tumors <-
    Mutations %>%
    mutate(SampleName = paste0(SampleID,'_',Region)) %>%
    left_join(PatientOverview, by = c('SampleID' = 'cruk_id')) %>% 
    left_join(SampleOverview_Bakir, by = c('SampleName' = 'region')) %>%
    mutate_at(vars(histology_lesion1,histology_lesion2), ~recode(.x,'Invasive adenocarcinoma'='LUAD','Squamous cell carcinoma' = 'LUSC')) %>%
    mutate(Taken_at_primary_surgery = grepl('SU',Region),
           Is_Primary_lung_tumor = grepl('SU_T1',Region),
           Is_Second_primary = grepl('T2|T3',Region),
           Is_Second_primary_lung_tumor = grepl('SU_T2|SU_T3',Region),
           Sampled_recurrence = grepl('BR',Region),
           Sampled_progression = grepl('BP',Region),
           sampleType = case_when(
               grepl('SU_T1',Region) ~ 'primary',
               grepl('SU_T2',Region) ~ 'second primary',
               grepl("SU.*LN|LN.*SU",Region) ~ 'metastasis',
               grepl("SU_T2|SU_T3",Region) ~ 'second primary',
               grepl("BP",Region) ~ 'metastasis',
               grepl("BR",Region) ~ 'metastasis',
               grepl("MR",Region) ~ 'metastasis',
               TRUE ~ 'unknown'),
           sampleTypeDetail = case_when(
               grepl('SU_T1',Region) ~ 'primary',
               grepl("SU.*LN|LN.*SU",Region) ~ 'synchronous LN met',
               grepl("SU_T2|SU_T3",Region) ~ 'synchronous primary',
               grepl("LN", Region) ~ 'metachronous LN met',
               grepl('T1|T2|T3', Region) ~ 'metachronous met',
               TRUE ~ 'unknown'),
           Relapse_cat =  coalesce(Relapse_cat, Relapse_cat_new),
           Location = case_when(
               sampleType %in%  c('synchronous primary','primary') ~ 'Intrathoracic',
               grepl('LN',Region) ~ 'Intrathoracic',
               TRUE ~ Relapse_cat),
           sampleTypeDetailLocation = paste0(Location,': ',sampleTypeDetail),

           Histology = dplyr::case_when(
                                  Is_Second_primary_lung_tumor ~histology_lesion2,
                                  !Is_Second_primary_lung_tumor ~ histology_lesion1,
                                  !histology_multi_full %in% c('LUAD','LUSC') ~ 'Other'


                              )) %>%
select(SampleID,SampleName,sampleType,sampleTypeDetail,Location,sampleTypeDetailLocation,Relapse_cat,Taken_at_primary_surgery,Is_Primary_lung_tumor,Is_Second_primary,Is_Second_primary_lung_tumor,Sampled_recurrence,Sampled_progression,histology_multi_full,Histology,histology_lesion1,histology_lesion2)


#write.table(Tumors,file="Tumors.txt",sep= "\t",row.names = F,quote = F)
#-------------------------------------------------------------------------------
# 3.1 CREATE CLONAL PAIRS
# Create pairs of clonal tumors. Each pair consists of Primary tumour and a 
# randomly selected intrapulmonary metastasis
#-------------------------------------------------------------------------------
# Fetch patients with at least 1 metastasis
metPatients <- Tumors %>%filter(sampleType== 'metastasis') %>% pull(SampleID) %>% unique()
# set random seed
set.seed(1)
Selected_pairs <-
    Tumors %>%
    # select patient with a metastasis
    filter(SampleID %in% metPatients) %>% #filter(Location == 'IntraThoracic' ) %>%
    # Create primary/metastasis groups per patient
    group_by(SampleID,sampleType,Location) %>%
    # Randomly slice one within group row
    slice_sample(n=1) %>%
    ungroup()
# Fetch Intrahoracic pairs
Intrathoracic_pairs <-
    Selected_pairs %>% 
    filter(Location == 'Intrathoracic') %>%
    tidyr::pivot_wider(id_cols = c(SampleID),values_from = SampleName,names_from = sampleType) %>%
    mutate(Pair_type = 'IPM',True_Clonality = 'Clonal')


# Fetch patients without intrathoracic pair
Patients_without_IntraMet <- Intrathoracic_pairs %>% filter(is.na(metastasis)) %>% pull(SampleID)
# Fetch extrathoracic pairs
Extrathoracic_pairs <-
    Selected_pairs %>%
    filter(SampleID %in% Patients_without_IntraMet) %>%
    tidyr::pivot_wider(id_cols = SampleID,values_from = SampleName,names_from = sampleType) %>%
    mutate(Pair_type = 'Extrathoracic metastasis',True_Clonality = 'Clonal')

# Concatenate Intra- with Extrathoracic pairs
Clonal_pairs <- rbind(Intrathoracic_pairs,Extrathoracic_pairs) %>% filter(!is.na(metastasis))

#-------------------------------------------------------------------------------
# 3.2 CREATE NON-CLONAL PAIRS
# Create pairs of Non-clonal tumors. Prioritize SPLCs (from the same histology)
# Collision tumors are included as SPLCs. Create all possible combinations for
# patient with 3 collision tumors.
#
# After creating pairs of  SPLCs, create Interpatient pairs such that the 
# same number of clonal pairs is reached. Do this by sampling WITH replacement.
#-------------------------------------------------------------------------------
# Define number of pairs
nPairs <- nrow(Clonal_pairs)

# INCLUDE:
# CRUK0030: two primaries at surgery, both LUAD
# CRUK0249: two primaries at surgery, both LUAD
# CRUK0519: two primaries at surgery, both LUAD
# CRUK0620: two primaries at surgery, both LUAD
# CRUK0881: collision at surgery, both LUAD # histology multi full says LUADx2

# EXCLUDE:
# CRUK0223: two primaries at surgery, one LUAD one LUSC
# CRUK0372: two primaries at surgery, one LUAD and one pleomorphic 
# CRUK0586: two primaries at surgery, one LUAD one LUSC
#---------------------------------------------------------------------------------------------
# INVESTIGATE:
# - CRUK0039: 2 regions sampled, apparently collisions, both LUAD, check mutations of regions
#             2 REGIONS HAVE ~1100 INTERSECT AND ~100 OUTERSECT, I ASSUME CLONAL
#
# * CRUK0704: collision at surgery, both LUAD # histology multi full says LUADx3 but only 2 T's?
#             Quite some regions are sampled, check which are the collision regions
#             JUDGING BY THE MUTATIONS FOUND WE HAVE 3 SEPERATE TUMORS:
#                    SU_T2_*
#                    SU_T1_R3-R5
#                    SU_T1_R1-R2
#
# - CRUK0057: two histologies within same tumor, LUAD/minimally invasive LUAD
#             2 REGIONS  HAVE ~1100 INTERSECT AND ~100 OUTERSECT, I ASSUME CLONAL
#
# * CRUK0301:  4 regions sampled, histology_multi_full == LUADx2, one FLN at surgery
#             JUDGING BY THE MUTATIONS FOUND WE HAVE 2 SEPERATE TUMORS AND ONE LN:
#                    SU_T1_R1-R2 (and LN which is clonal)
#                    SU_T1_R3-R4
#
# - CRUK0418  5 regions sampled, histology_multi_full == LUADx2
#             Judging by the amount of shared mutations, they seem clonal
#             3367 mutations are shared amongst all samples
#
# - CRUK0441: 2 regions sampled, histology_multi_full == LUADx2
#             410 outersect, and 172 in intersect (possibly non-clonal??)
#             PatientOverview says: Lesion2 is not sampled
#
# - CRUK0460: 3 regions sampled, histology_multi_full == LUADx2
#             327 in intersect in all 3 regions
#             PatientOverview says: Lesion 2 is not sampled

# Fetch SPLCs designated primary/second primary pairs
# No LUSC patient has 2 SPLCs with LUSC subtype  
if(subtype == 'LUSC'){
    SPLCs <- data.frame()
}else if(subtype == 'LUAD'){
    SPLCs1 <-
        Tumors %>%
        filter(SampleID %in% c('CRUK0030','CRUK0249','CRUK0519','CRUK0620','CRUK0881')) %>%
        filter(sampleType %in% c('primary','second primary')) %>%
        group_by(SampleID,sampleType) %>%
        # Randomly slice one within group row
        slice_sample(n=1) %>%
        ungroup() %>%
        tidyr::pivot_wider(id_cols = c(SampleID),values_from = SampleName,names_from = sampleType) %>%
        rename(second_primary = `second primary`)

    # Fetch SPLCs through manual curation
    # CRUK0704
    set.seed(1)
    SPLCs_CRUK0704 <-
        data.frame(SampleID = 'CRUK0704',
                   primary = c(sample(paste0('CRUK0704_SU_T1.R',seq(1,2)),1), # 1
                               sample(paste0('CRUK0704_SU_T1.R',seq(1,2)),1), # 1
                               sample(paste0('CRUK0704_SU_T1.R',seq(3,5)),1)), # 2
                   second_primary = c(sample(paste0('CRUK0704_SU_T2.R',seq(1,4)),1), # 3
                                      sample(paste0('CRUK0704_SU_T1.R',seq(3,5)),1), # 2
                                      sample(paste0('CRUK0704_SU_T2.R',seq(1,4)),1))# 3
                   )
    # CRUK0301
    set.seed(1)
    SPLCs_CRUK0301 <-
        data.frame(SampleID = 'CRUK0301',
                   primary = sample(paste0('CRUK0704_SU_T1.R',seq(1,2)),1),
                   second_primary =  sample(paste0('CRUK0704_SU_T1.R',seq(3,4)),1))


    # Concatenate SPLCs
    SPLCs <- rbind(SPLCs1,SPLCs_CRUK0704,SPLCs_CRUK0301) %>%
        mutate(Pair_type = 'SPLC',True_Clonality = 'Non-Clonal')

}
# Add Interpatient tumor pairs to match the number of clonal pairs
# Only select primary tumors
#---------------------------------------------------------------------------------------------
nInterpatient_pairs <- nPairs - nrow(SPLCs)

# Select primary tumors and shuffle randomly
Primary_tumors1 <- Tumors %>% filter(sampleType == 'primary') %>% sample_frac() %>% select(SampleID,SampleName)
# Reshuffle again
Primary_tumors2 <- Primary_tumors1 %>% sample_frac() %>% rename(SampleID2 = SampleID,SampleName2=SampleName)
# concatenate and filter out interpatient pairs
Interpatient_pairs <-
    Primary_tumors1 %>%
    cbind(Primary_tumors2) %>%
    filter(SampleID != SampleID2) %>%
    # randomly select interpatient pairs
    slice_sample(n=nInterpatient_pairs) %>%
    # Select variables
    mutate(SampleID = paste0(SampleID,'-',SampleID2)) %>%
    rename(primary = SampleName,second_primary = SampleName2) %>%
    mutate(Pair_type = 'Interpatient',True_Clonality = 'Non-Clonal') %>%
    select(SampleID,primary,second_primary,Pair_type,True_Clonality)


# Concatenate results

Non_Clonal_pairs <- rbind(SPLCs,Interpatient_pairs)

#-------------------------------------------------------------------------------
# 3.2 Combine all pairs
colnames(Non_Clonal_pairs)[c(2,3)] <- c('Tumor1','Tumor2')
colnames(Clonal_pairs)[c(2,3)] <- c('Tumor1','Tumor2')
Tumor_pairs <- rbind(Clonal_pairs,Non_Clonal_pairs)

#-------------------------------------------------------------------------------
# 4.1 Write to file
#-------------------------------------------------------------------------------
write.table(Tumor_pairs, file = output, sep = '\t', row.names = F,quote = F  )



function(x){

    library(ComplexHeatmap)
    DummyData <- fastDummies::dummy_cols(Tumors %>% select(-c(SampleID,SampleName,histology_multi_full ,histology_lesion1 ,histology_lesion2)))
    rownames(DummyData) <- Tumors$SampleName

    left_annotation <- rowAnnotation(
        Is_Primary = Tumors$Is_Primary_lung_tumor,
        Is_SPLC = Tumors$Is_Second_primary_lung_tumor,
        Location = Tumors$Location,
        Sampled_at_primary_surgery = Tumors$Taken_at_primary_surgery,
        Sampled_at_recurrence = Tumors$Sampled_recurrence,
        Sampled_at_progression = Tumors$Sampled_progression,
        SampleType = Tumors$sampleTypeDetail,
        Subtype = Tumors$histology_multi_full,
        simple_anno_size = unit(1, "cm"))
    top_annotation <- columnAnnotation(
        Is_Primary = Tumors$Is_Primary_lung_tumor,
        Is_SPLC = Tumors$Is_Second_primary_lung_tumor,
        Location = Tumors$Location,
        Sampled_at_primary_surgery = Tumors$Taken_at_primary_surgery,
        Sampled_at_recurrence = Tumors$Sampled_recurrence,
        Sampled_at_progression = Tumors$Sampled_progression,
        SampleType = Tumors$sampleTypeDetail,
        Subtype = Tumors$histology_multi_full,
        simple_anno_size = unit(1, "cm"))




    ht <- Heatmap(as.matrix(dist(DummyData,method = 'binary')),
                  cluster_rows = T,
                  cluster_columns = T,
                  top_annotation = top_annotation,
                  show_column_names = F,
                  show_row_names = F)



    png('test.png',width=2480/3,height= 2480/3)
    draw(ht)
    dev.off()
}
