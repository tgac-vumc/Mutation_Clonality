#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Evaluate performance.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Call clonality and evaluate performance
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
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(pROC)))

#-------------------------------------------------------------------------------
# 1.1 Parse snakemake objects
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_metrics <- snakemake@input[["Metrics"]]
    input_metrics_MC <- snakemake@input[["Metrics_MC"]]
    input_mutations <- snakemake@input[["Shared_Mutations"]]
    input_table <- snakemake@input[["Table"]]
    output <-  snakemake@output[["Evaluation_metrics"]]
    output_boxplots <- snakemake@output[['Boxplot_nMatches']]
    output_piecharts <- snakemake@output[['Piechart_SharedMutations']]
    output_piechart_evaluablepatients <- snakemake@output[['Piechart_EvaluablePatients']]
    output_performance <- snakemake@output[['Barchart_performance']]
}else{
    input_metrics <- c('output/TRACERx421/InhouseLungPanel/Clonality_metrics_LUSC_WorseCaseScenario.txt','output/TRACERx421/IlluminaTrueSightTumor500/Clonality_metrics_LUSC_WorseCaseScenario.txt','output/TRACERx421/FoundationOneCDx/Clonality_metrics_LUSC_WorseCaseScenario.txt','output/TRACERx421/InhouseLungPanel/Clonality_metrics_LUSC_RandomScenario.txt','output/TRACERx421/IlluminaTrueSightTumor500/Clonality_metrics_LUSC_RandomScenario.txt','output/TRACERx421/FoundationOneCDx/Clonality_metrics_LUSC_RandomScenario.txt')
    input_metrics_MC <-  c('output/TRACERx421/InhouseLungPanel/Clonality_Calls_MC_LUSC_WorseCaseScenario.txt','output/TRACERx421/IlluminaTrueSightTumor500/Clonality_Calls_MC_LUSC_WorseCaseScenario.txt','output/TRACERx421/FoundationOneCDx/Clonality_Calls_MC_LUSC_WorseCaseScenario.txt','output/TRACERx421/InhouseLungPanel/Clonality_Calls_MC_LUSC_RandomScenario.txt','output/TRACERx421/IlluminaTrueSightTumor500/Clonality_Calls_MC_LUSC_RandomScenario.txt','output/TRACERx421/FoundationOneCDx/Clonality_Calls_MC_LUSC_RandomScenario.txt')
    input_mutations <-  c('output/TRACERx421/InhouseLungPanel/Shared_mutations_LUSC_WorseCaseScenario.txt','output/TRACERx421/IlluminaTrueSightTumor500/Shared_mutations_LUSC_WorseCaseScenario.txt','output/TRACERx421/FoundationOneCDx/Shared_mutations_LUSC_WorseCaseScenario.txt','output/TRACERx421/InhouseLungPanel/Shared_mutations_LUSC_RandomScenario.txt','output/TRACERx421/IlluminaTrueSightTumor500/Shared_mutations_LUSC_RandomScenario.txt','output/TRACERx421/FoundationOneCDx/Shared_mutations_LUSC_RandomScenario.txt')
    input_table <- 'output/TRACERx100/Tables/TableXX_NumberOfEvaluableSamples.txt'
}

#-------------------------------------------------------------------------------
# 1.1 Read data 
#-------------------------------------------------------------------------------
# read metrics
Clonality_metrics <-
    tibble::tibble(file = input_metrics) %>%
    mutate(panel = purrr::map_chr(file, ~strsplit(.x,'/')[[1]][3]),
           dataset = purrr::map_chr(file, ~strsplit(.x,'/')[[1]][2]),
           subtype = gsub('.txt','',purrr::map_chr(file, ~strsplit(.x,'_')[[1]][3])),
           scenario = gsub('.txt','',purrr::map_chr(file, ~strsplit(.x,'_')[[1]][4])),
           # read file and get number of unique samples
           data = purrr::map(file,read.delim)) %>%
    tidyr::unnest() %>%
    select(-file) %>%
    filter(!is.na(scenario))


tmp <- Clonality_metrics %>% filter(subtype == 'LUSC') %>% filter(panel == 'FoundationOneCDx') %>%
    filter(scenario == 'RandomScenario')


Clonality_metrics_MC <-
    tibble::tibble(file = input_metrics_MC) %>%
    mutate(panel = purrr::map_chr(file, ~strsplit(.x,'/')[[1]][3]),
           dataset = purrr::map_chr(file, ~strsplit(.x,'/')[[1]][2]),
           subtype = gsub('.txt','',purrr::map_chr(file, ~strsplit(.x,'_')[[1]][4])),
           scenario = gsub('.txt','',purrr::map_chr(file, ~strsplit(.x,'_')[[1]][5])),

           # read file and get number of unique samples
           data = purrr::map(file,read.delim)) %>%
    tidyr::unnest() %>%
    mutate(comparison = paste0(Sample1,'-',Sample2))  %>%
    select(-file) %>%
    filter(!is.na(scenario))

# read shared_mutations
Shared_mutations <-
    tibble::tibble(file = input_mutations) %>%
    mutate(panel = purrr::map_chr(file, ~strsplit(.x,'/')[[1]][3]),
           subtype = gsub('.txt','',purrr::map_chr(file, ~strsplit(.x,'_')[[1]][3])),
           dataset = purrr::map_chr(file, ~strsplit(.x,'/')[[1]][2]),
           scenario = gsub('.txt','',purrr::map_chr(file, ~strsplit(.x,'_')[[1]][4])),
           # read file and get number of unique samples
           data = purrr::map(file,~read.delim(.x , colClasses = rep('character',5)))) %>%
    tidyr::unnest()  %>%
    select(-file) %>%
    filter(!is.na(scenario))

# Read number of evaluable samples
nEvaluable_patients <-
    tibble::tibble(file = input_table) %>%
    mutate(dataset = purrr::map_chr(file, ~strsplit(.x,'/')[[1]][2]),
            data = purrr::map(file,~read.delim(.x , colClasses = rep('character',5)))) %>%
    tidyr::unnest() %>%
    select(-file)


#-------------------------------------------------------------------------------
# 2.1 Determine clonality 
#-------------------------------------------------------------------------------
# Determine clonality using three different methods:
# 1) pval_Jaccard < 0.05
# 2) LRpvalue < 0.05
# 3) pval_Jaccard < 0.05 & LRpvalue < 0.05
Clonality_metrics <-
    Clonality_metrics %>%
    left_join(Clonality_metrics_MC, by = c('scenario','dataset','panel','subtype','True_clonality','comparison')) %>% 
    dplyr::mutate(
               Clonality_Jaccard_test = factor(ifelse(pval_Jaccard < 0.05,1,0), levels=c(0,1)),
               Clonality_LR_test = factor(ifelse(LRpvalue < 0.05,1,0), levels=c(0,1)),
               Clonality_TwoMetric =  factor(ifelse(pval_Jaccard < 0.05 & LRpvalue < 0.05,1,0), levels=c(0,1)),
               Clonality_MC_test = factor(
                   case_when(
                       Clonality_MC == 'Clonal' ~ 1,
                       grepl('non-Clonal',Clonality_MC) ~ 0,
                       ),
                   levels=c(0,1)),
               Clonality = factor(ifelse(True_clonality == 'Clonal',1,0), levels=c(0,1)),
               # Remove mirrored comparions for removal of duplicates
               comparison = purrr::map_chr(comparison,~paste(sort(strsplit(.x,'-')[[1]]),collapse='-'))) %>% 
    # Filter out duplicate rows by selecting first entry of unique clonality metrics
    group_by(comparison,panel,scenario) %>%
    filter(row_number(comparison) == 1) %>%
    ungroup()



#-------------------------------------------------------------------------------
# 2.2 Count number of 'Inconclusive' and Misclassified by MC algorithm
#-------------------------------------------------------------------------------
Inconclusives <-
    Clonality_metrics %>%
    filter(Clonality_MC == 'Inconclusive' ) %>% 
    group_by(scenario,dataset,panel,subtype) %>%
    summarise(Inconclusive = dplyr::n()) %>%
    replace(is.na(.), 0)


Misclassified_TwoMetric <-
    Clonality_metrics %>%
    filter(Clonality_TwoMetric != Clonality ) %>%
    group_by(scenario,dataset,panel,subtype) %>%
    summarise(Misclassified_TwoMetric = dplyr::n()) %>%
    replace(is.na(.), 0)

Misclassified_MC <-
    Clonality_metrics %>%
    filter(Clonality_MC_test != Clonality ) %>%
    group_by(scenario,dataset,panel,subtype) %>%
    summarise(Misclassified_MC = dplyr::n()) %>%
    replace(is.na(.), 0)

Misclassified_FN_MC <-
    Clonality_metrics %>%
    filter(Clonality_MC_test == '0' & Clonality == '1') %>%
    mutate(reason = as.factor(paste0('FN: ',reason))) %>%
    group_by(scenario,dataset,panel,subtype) %>%
    mutate(FN_MC = dplyr::n())%>%
    group_by(scenario,dataset,panel,subtype,reason,FN_MC) %>%
    summarise(count = dplyr::n()) %>%
    tidyr::pivot_wider(names_from = reason, values_from = count) %>% 
    replace(is.na(.), 0)

Misclassified_FP_MC <-
    Clonality_metrics %>%
    filter(Clonality_MC_test == '1' & Clonality == '0') %>%
    mutate(reason = as.factor(paste0('FP: ',reason))) %>%
    group_by(scenario,dataset,panel,subtype) %>%
    mutate(FP_MC = dplyr::n()) %>%
    group_by(scenario,dataset,panel,subtype,reason,FP_MC) %>%
    summarise(count = dplyr::n()) %>%
    tidyr::pivot_wider(names_from = reason, values_from = count) %>%
    replace(is.na(.), 0)

# TODO FIX reasons!!
Clonality_metrics$reason %>% unique()


Misclassified_FN_TwoMetric <-
    Clonality_metrics %>%
    filter(Clonality_TwoMetric == '0' & Clonality == '1') %>%
    group_by(scenario,dataset,panel,subtype) %>%
    summarise(FN_TwoMetric = dplyr::n()) 
    
Misclassified_FP_TwoMetric <-
    Clonality_metrics %>%
    filter(Clonality_TwoMetric == '1' & Clonality == '0') %>%
    group_by(scenario,dataset,panel,subtype) %>%
    summarise(FP_TwoMetric = dplyr::n())

    
nTP <-
    Clonality_metrics %>%
    filter( Clonality_MC_test == '1' &  Clonality == '1') %>%
    group_by(scenario,dataset,panel,subtype) %>%
    summarise(nTP = dplyr::n())

nNegativePred <-
    Clonality_metrics %>%
    filter( Clonality_MC_test == '0') %>%
    group_by(scenario,dataset,panel,subtype) %>%
    summarise(nNegativesPred = dplyr::n())

nPositivePred <-
    Clonality_metrics %>%
    filter( Clonality_MC_test == '1') %>%
    group_by(scenario,dataset,panel,subtype) %>%
    summarise(nPositivesPred = dplyr::n())

nTN <-
    Clonality_metrics %>%
    filter(Clonality_MC_test == '0' &  Clonality == '0') %>%
    group_by(scenario,dataset,panel,subtype) %>%
    summarise(nTN = dplyr::n())

nTP_TwoMetric <-
    Clonality_metrics %>%
    filter( Clonality_TwoMetric == '1' &  Clonality == '1') %>%
    group_by(scenario,dataset,panel,subtype) %>%
    summarise(nTP_TwoMetric = dplyr::n())

nNegativePred_TwoMetric <-
    Clonality_metrics %>%
    filter( Clonality_TwoMetric == '0') %>%
    group_by(scenario,dataset,panel,subtype) %>%
    summarise(nNegativesPred_TwoMetric = dplyr::n())

nPositivePred_TwoMetric <-
    Clonality_metrics %>%
    filter( Clonality_TwoMetric == '1') %>%
    group_by(scenario,dataset,panel,subtype) %>%
    summarise(nPositivesPred_TwoMetric = dplyr::n())

nTN_TwoMetric <-
    Clonality_metrics %>%
    filter(Clonality_TwoMetric == '0' &  Clonality == '0') %>%
    group_by(scenario,dataset,panel,subtype) %>%
    summarise(nTN_TwoMetric = dplyr::n())







#-------------------------------------------------------------------------------
# 2.2 Calculate specificity and sensitivity
#-------------------------------------------------------------------------------
# Fetch evaluation metrics per panel and subtype
Evaluation_metrics <-
    Clonality_metrics %>%
    group_by(scenario,dataset,panel,subtype) %>%
    summarise(
        #Sensitivity_LR = caret::sensitivity(table(Clonality,Clonality_LR_test)),
        #Specificity_LR = caret::specificity(table(Clonality,Clonality_LR_test)),
        #Sensitivity_Jaccard = caret::sensitivity(table(Clonality,Clonality_Jaccard_test)),
        #Specificity_Jaccard = caret::specificity(table(Clonality,Clonality_Jaccard_test)),
        Sensitivity_TwoMetric = caret::sensitivity(table(Clonality,Clonality_TwoMetric), negative = '0', positive = '1'  ),
        Specificity_TwoMetric = caret::specificity(table(Clonality,Clonality_TwoMetric), negative = '0', positive = '1' ),
        Sensitivity_MC = caret::sensitivity(table(Clonality,Clonality_MC_test), negative = '0', positive = '1' , na.rm=T),
        Specificity_MC = caret::specificity(table(Clonality,Clonality_MC_test), negative = '0', positive = '1' , na.rm=T),
        NComparisons = dplyr::n()) %>%
    ungroup()


#-------------------------------------------------------------------------------
# 2.3 Calculate panel size and Ngenes with detected mutation in TRACERx
#-------------------------------------------------------------------------------
# intialize empty df
PanelSize <- data.frame()
# Iterate over panels
for(panel_name in unique(Evaluation_metrics$panel)){
    # Fetch panel size
    if(panel_name == 'FoundationOneCDx'){
        source('scripts/Fetch_FoundationOne_regions.R')
        exonSize <- as.numeric(system(paste0("grep -E '", paste0(ExonGenes,collapse = "|"),"' manifest/KappaHyperExome.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'"), intern = T))
        intronSize <- sum(as.numeric(panel$chromEnd) - as.numeric(panel$chromStart))
        size <- exonSize+intronSize
    }else{
        size <- system(paste0("awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' manifest/",panel_name,".bed"), intern = T)
    }
    # Fetch number of genes with mutation in panel
    Ngenes_with_mutation <- length(unique(c(
        read.delim(paste0('data/TRACERx421/',panel_name,'/Selected_mutations_NSCLC.txt'))[,3],
        read.delim(paste0('data/TRACERx100/',panel_name,'/Selected_mutations_NSCLC.txt'))[,3])))
    # add to df
    PanelSize <- rbind(PanelSize, data.frame(Size= size , panel = panel_name,Ngenes_with_mutation = Ngenes_with_mutation))
}




#-------------------------------------------------------------------------------
# 2.4 Perform ROC analysis with number of matching mutations
#-------------------------------------------------------------------------------
# Perform ROC analysis with number of mutations as predictor and calculate AUC
#ROC_analysis <-
#    Clonality_metrics %>%
#    group_by(scenario,dataset,panel,subtype) %>%
#    summarise(
#        ROC_model = list(roc(response=Clonality,predictor=n_match)),
#        AUC = purrr::map_dbl(ROC_model,auc)) %>%
#    select(scenario, dataset, panel , subtype, AUC)


#-------------------------------------------------------------------------------
# 2.4 Perform ROC analysis with number of matching mutations
#-------------------------------------------------------------------------------
# Join ROC AUC, Inconclusives,number of evaluables, misclassifications and panel size
Evaluation_metrics <-
    Evaluation_metrics %>%
    #left_join(ROC_analysis)  %>%
    left_join(nEvaluable_patients) %>%
    left_join(Inconclusives) %>%
    left_join(Misclassified_MC) %>%
    left_join(Misclassified_FN_MC) %>%
    left_join(Misclassified_FP_MC) %>%
    left_join(Misclassified_TwoMetric) %>%
    left_join(PanelSize) %>%
    left_join(nTP) %>%
    left_join(nTN) %>%
    left_join(nNegativePred) %>%
    left_join(nPositivePred) %>%
    left_join(nTP_TwoMetric) %>%
    left_join(nTN_TwoMetric) %>%
    left_join(nNegativePred_TwoMetric) %>%
    left_join(nPositivePred_TwoMetric)

#-------------------------------------------------------------------------------
# 3.1 Plot data
#-------------------------------------------------------------------------------
# 3.1.1 Plot boxplots number of matching mutations per panel for LUAD/LUSC
#       in clonal vs non-clonal
#-------------------------------------------------------------------------------
palette <- palette.colors(palette = "Okabe-Ito")
palette <- palette[!names(palette) %in% c('black','yellow','skyblue')]
palette_clonal <- palette[c('bluishgreen','vermillion')]

pdf(output_boxplots,width=10,height=8)
#pdf('/net/beegfs/cfg/tgac/jjanssen4/boxplot_test.pdf', width = 6, height = 4)
Clonality_metrics %>%
    ggplot(aes(panel, n_match,fill=True_clonality, color = True_clonality)) +
    geom_boxplot() +
    scale_fill_manual(values = as.vector(palette_clonal)) +
    scale_color_manual(values = as.vector(palette_clonal)) +

    facet_wrap(~subtype + True_clonality + dataset, switch = 'y', scales = 'free') +
    labs(fill='Clonality',x='Panel',y='Number of matching mutations') +
    theme_bw(base_size = 12) +
    theme(axis.text.x=element_text(size=12,angle=45,hjust=1),
          legend.position='top')  + geom_jitter()
dev.off()
#-------------------------------------------------------------------------------
# 3.1.2 Plot Piechart of %-age of shared mutations in per panel for LUAD/LUSC
#       in clonal vs non-clonal
#-------------------------------------------------------------------------------
pdf(output_piecharts, height = 5, width = 5)
Shared_mutations %>%
    filter(panel == 'KappaHyperExome') %>%
    # Create gene categories
    mutate(Gene = case_when(
               Hugo_Symbol == 'KRAS' ~ 'KRAS',
               Hugo_Symbol == 'EGFR' ~ 'EGFR',
               Hugo_Symbol == 'TP53' ~ 'TP53',
               Hugo_Symbol == 'PIK3CA' ~ 'PIK3CA',
               Hugo_Symbol == 'NFE2L2' ~ 'NFE2L2',
               TRUE ~ 'Other')) %>%
    # Calculate number of shared mutations
    group_by(subtype,dataset,True_clonality) %>%
    mutate(Ntot = dplyr::n()) %>%
    group_by(subtype,dataset,True_clonality,Gene) %>%
    # Calculate percentage
    dplyr::summarize(Percentage = dplyr::n() / Ntot * 100) %>%
    unique() %>%
    # Plot piecharts
    ggplot(aes(x="", y=Percentage, fill=Gene)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    facet_grid(rows=vars(subtype,dataset) , cols = vars(True_clonality), switch = 'y') +
    theme_bw(base_size = 12) +
    theme(axis.ticks=element_blank(), axis.title=element_blank(), 
          axis.text.y = element_blank(), panel.grid  = element_blank(),
          axis.text.x = element_blank()) +
    scale_fill_manual(values=as.vector(palette))
dev.off()


#-------------------------------------------------------------------------------
# 3.1.2 Plot Piechart of %-age of evaluable patients per panel/subtype
#-------------------------------------------------------------------------------
pdf(output_piechart_evaluablepatients, height = 3 , width = 12)
nEvaluable_patients %>% 
    filter(subtype != 'NSCLC') %>%
    # create factor ordered by percentage
    mutate(panel = factor(panel, levels = rev(unique(nEvaluable_patients$panel)))) %>%
    # calculate percentages
    group_by(subtype,dataset,panel) %>%
    summarize(Evaluable = as.numeric(fraction_evaluable_tumors)*100,
           `Not Evaluable` = 100 - Evaluable) %>%
    select(panel,dataset,subtype,Evaluable,`Not Evaluable`) %>%
    # retrieve long format
    tidyr::pivot_longer(cols = -c(panel,dataset,subtype), values_to = 'Percentage') %>% 
    # Plot piecharts
    ggplot(aes(x="", y=Percentage, fill=name)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    facet_grid(rows=vars(subtype) , cols = vars(panel), switch = 'y') +
    theme_bw(base_size = 12) +
    theme(axis.ticks=element_blank(), axis.title=element_blank(), 
          axis.text.y = element_blank(), panel.grid  = element_blank(),
          axis.text.x = element_blank(),
          legend.title= element_blank(),
          legend.position='top',
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,-10,-10),
          strip.text.x = element_text(size = 5)) +
    scale_fill_manual(values=c('darkgrey','lightgrey'))
dev.off()

 
#-------------------------------------------------------------------------------
# 3.1.2 Plot Barchart performance for each panel per subtype
#-------------------------------------------------------------------------------
pdf(output_performance,width=12, height = 4)
Evaluation_metrics %>%
    select(panel,dataset,subtype,Sensitivity_TwoMetric,Specificity_TwoMetric) %>%
    rename(Sensitivity = Sensitivity_TwoMetric,Specificity =Specificity_TwoMetric) %>%
    tidyr::pivot_longer(cols = -c(panel,dataset,subtype), names_to = 'Measure') %>%
    mutate(panel = factor(panel,levels = unique(arrange(Evaluation_metrics[Evaluation_metrics$subtype == 'LUAD',],desc(Sensitivity_TwoMetric),desc(Specificity_TwoMetric))$panel))) %>%
    ggplot(aes(panel, value, fill = Measure))  +
    geom_bar(position="dodge", stat="identity") +
    labs(fill='Measure',x='Panel',y='Value') +
    theme_bw(base_size = 12) +
    theme(axis.text.x=element_text(size=12,angle=45,hjust=1),
          legend.position='top',
          legend.box.margin=margin(-10,-10,-10,-10)) +
    scale_fill_manual(values=as.vector(palette)) +
    facet_wrap(~subtype + dataset, scales = 'free')
dev.off()
#-------------------------------------------------------------------------------
# 4.1 Write to file
#-------------------------------------------------------------------------------
write.table(Evaluation_metrics, file = output, sep = '\t', row.names = F,quote = F  )

# Write table to evaluate the misclassifications
Clonality_metrics %>%
    filter(Clonality_MC_test != Clonality ) %>%
    select(-c(n1,n2,n_match,LRstat,maxKsi,LRpvalue,pval_Jaccard,Jaccard,Clonality_Jaccard_test,Clonality_LR_test,Sample1,Sample2)) %>%
    write.table(file = 'Misclassifications_MC_algorithm.txt', sep = '\t', row.names = F,quote = F )


#-------------------------------------------------------------------------------
# 4.1 Code for presentation
#-------------------------------------------------------------------------------
Misclassifications <- read.delim('Misclassifications_MC_algorithm.txt')
Evaluation_metrics <- read.delim('output/Tables/TableXX_Sensitivity_Specificity.txt')

library(scales)
pdf('Mutation_panel_sizes_and_genes.pdf', height = 6, width = 8)
Evaluation_metrics %>%
    select(panel,dataset,Ngenes_with_mutation,Size) %>%
    unique() %>%
    filter(Ngenes_with_mutation != 1,dataset == 'TRACERx421') %>%
    mutate(panel = factor(panel,levels=sort(panel,decreasing = T))) %>%
    ggplot(aes(Size, Ngenes_with_mutation,color = panel)) +
    geom_point(size = 3)+
    ggrepel::geom_text_repel(aes(Size, Ngenes_with_mutation,label=panel, color = panel),size=5) +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks() +
    theme_classic(base_size = 16) +
    theme(legend.position = "none") +
    labs(y='Number of genes',x='Size (bp)')
dev.off()

pdf('NumberOfMatchingMutations.pdf', height = 4, width = 6)
palette_clonal <- c('#6478a7','#aba148')
Clonality_metrics %>%
    ggplot(aes(True_clonality, n_match,color=True_clonality,fill= True_clonality)) +
    geom_dotplot(binaxis='y', stackdir='center',stackratio=2, dotsize=0.25) +
    scale_fill_manual(values = as.vector(palette_clonal)) +
    scale_color_manual(values = as.vector(palette_clonal)) +

    facet_wrap(~panel, scales = 'free') +
    labs(fill='Clonality',x='Panel',y='Number of matching mutations') +
    theme_bw(base_size = 14) +
    theme(axis.text.x=element_text(size=12,angle=45,hjust=1),
          legend.position='top') 
dev.off()

#-------------------------------------------------------------------------------
# 4.1 Code for Visual Abstract
#-------------------------------------------------------------------------------
plot_data_LUAD <-
    data.frame(
        Panel = factor(c(rep('MayoComplete',6),rep('InhousePanel',6),rep('TSO500',6),rep('F1CDX',6)), levels = c('MayoComplete','InhousePanel','F1CDX','TSO500')),
        Percentages = c(57,3,26,71,5,38,   #64,4,32, # Mayo
                        84,4,6,90,6,10,    #87,5,8, # Inhouse
                        92,4,4,95,2,3,     #93.5,3,3.5,# TSO500
                        91,4,5,95,2,3),     # 93,3,4), # F1CDX  
        Label = factor(rep(rep(c('Correct','Incorrect','Inconclusive'),2),4), levels = c('Incorrect','Inconclusive','Correct')),
        Subtype = 'LUAD')

plot_data_LUSC <-
    data.frame(
        Panel = factor(c(rep('MayoComplete',6),rep('InhousePanel',6),rep('TSO500',6),rep('F1CDX',6)), levels = c('MayoComplete','InhousePanel','F1CDX','TSO500')),
        Percentages = c(48,0,22,41,0,59,   #64,4,32, # Mayo
                        95,1,4,97,0,3,    #87,5,8, # Inhouse
                        95,1,4,97,0,3,     #93.5,3,3.5,# TSO500
                        95,1,4,97,0,3),     # 93,3,4), # F1CDX  
        Label = factor(rep(rep(c('Correct','Incorrect','Inconclusive'),2),4), levels = c('Incorrect','Inconclusive','Correct')),
        Subtype = 'LUSC')







plot_data_LUAD
plot_data
library(ggpubr)
rbind(plot_data_LUAD,plot_data_LUSC) %>%
    ggbarplot(
        data=.,
  x = c("Panel"), y = "Percentages", add = 'median_se',
  fill = "Label", palette = c('#ff4d55ff','#ccccccff','#51be60ff')) +
  facet_wrap(~ Subtype, strip.position = "bottom")






plot_data_LUAD %>%
ggplot(plot_data,aes(x=Panel, y=c(64,4,32,87,5,8,93.5,3,3.5,93,3,4), fill=Label)) + 
  geom_bar(stat="identity")

