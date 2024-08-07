from scripts.utils import *
from glob import glob
configfile: "config.yaml"

#+++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS,TARGET AND FILES ++++++++++++++++++++++++++++++++++++++++           
# 0.1 Specify wildcards
(panels,) = glob_wildcards("manifest/{panel}.bed")
subtypes = ['NSCLC','LUAD','LUSC']
datasets = ['TRACERx100','TRACERx421']
scenarios = ['WorstCaseScenario','WorseCaseScenario','RandomScenario']#,'All']
#----------------------------------------------------------------------------------------------------------------
# 0.2 specify target rule
rule all:
    input:
        expand('output/{dataset}/{panel}/Clonality_metrics_{subtype}_{scenario}.txt',dataset=datasets, panel = panels, subtype = subtypes,scenario = scenarios),
        expand('output/{dataset}/{panel}/Clonality_Calls_MC_{subtype}_{scenario}.txt',dataset=datasets, panel = panels, subtype = subtypes,scenario = scenarios),
	expand('output/Tables/TableXX_Sensitivity_Specificity.txt'),
        expand('output/{dataset}/Tables/Clonality_metrics_{panel}.xlsx',dataset=datasets,panel=panels),
        expand('output/{dataset}/Tables/TableXX_NumberOfEvaluableSamples.txt', dataset = datasets)
#++++++++++++++++++++++++++++++++++++++++++++++ 1 PREPROCESS DATA   +++++++++++++++++++++++++++++++++++++++++++++++++++
# 1.1 Fetch TRACERx mutations in panel regions
rule Filter_mutations:
    input:
        Mutations = lambda wildcards: 'data/TRACERx421_supplement/mutTableAll.cloneInfo.20220726.txt' if wildcards.dataset == 'TRACERx421' else 'data/TRACERx100_supplement/nejmoa1616288_appendix_2.xlsx' ,
        SampleOverview = 'data/TRACERx421_supplement/20221109_TRACERx421_all_patient_df.rds',
        panel = "manifest/{panel}.bed"
    output:
        Mutations = 'data/{dataset}/{panel}/Selected_mutations_{subtype}.txt'
    conda:
        'envs/R.yaml'
    script:
        'scripts/Filter_mutations_{wildcards.dataset}.R'

        
#++++++++++++++++++++++++++++++++++++++++++++++++ 2 CLONALITY   +++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2.1 Call clonality with using two-metrics classifier
rule Call_Clonality_TwoMetric:
    input:
        Mutations = 'data/{dataset}/{panel}/Selected_mutations_{subtype}.txt',
        MutationsExome = lambda wildcards:  'data/'+ wildcards.dataset+'/'+wildcards.panel+'/Selected_mutations_'+wildcards.subtype+'.txt' if wildcards.scenario == 'WorstCaseScenario' else 'data/'+wildcards.dataset+'/KappaHyperExome/Selected_mutations_'+wildcards.subtype+'.txt',
        sampleOverview = 'data/TRACERx421_supplement/sampleOverview.txt'
    output:
        Metrics = 'output/{dataset}/{panel}/Clonality_metrics_{subtype}_{scenario}.txt',
        Shared_Mutations = 'output/{dataset}/{panel}/Shared_mutations_{subtype}_{scenario}.txt'
    conda:
        'envs/R.yaml'
    script:
        'scripts/Call_Clonality_TwoMetric.R'

#----------------------------------------------------------------------------------------------------------------
# 2.2 Call clonality: Molecular classification algorithm
rule Call_Clonality_MC:
    input:
        Mutations = 'data/{dataset}/{panel}/Selected_mutations_{subtype}.txt',
        Mutations_inhouse = 'data/{dataset}/InhouseLungPanel/Selected_mutations_NSCLC.txt',
        MutationsExome = lambda wildcards:  'data/'+ wildcards.dataset+'/'+wildcards.panel+'/Selected_mutations_'+wildcards.subtype+'.txt' if wildcards.scenario == 'WorstCaseScenario' else 'data/'+wildcards.dataset+'/KappaHyperExome/Selected_mutations_'+wildcards.subtype+'.txt',
        Annotations = 'data/TRACERx421/KappaHyperExome/Selected_mutations_NSCLC.txt',
        sampleOverview = 'data/TRACERx421_supplement/sampleOverview.txt'
    output:
        Metrics = 'output/{dataset}/{panel}/Clonality_Calls_MC_{subtype}_{scenario}.txt'
    conda:
        'envs/R.yaml'
    script:
        'scripts/Call_Clonality_MCalgorithm.R'


#++++++++++++++++++++++++++++++++++++++++++++++++ 3 EVALUATION   +++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2.1 Call clonality
rule Evaluate_Performance:
    input:
        Shared_Mutations = expand('output/{dataset}/{panel}/Shared_mutations_{subtype}_{scenario}.txt', panel = panels, subtype = ['LUAD','LUSC'],scenario = scenarios,dataset=datasets),
        Metrics = expand('output/{dataset}/{panel}/Clonality_metrics_{subtype}_{scenario}.txt', panel = panels, subtype = ['LUAD','LUSC'],scenario = scenarios,dataset=datasets),
        Metrics_MC = expand('output/{dataset}/{panel}/Clonality_Calls_MC_{subtype}_{scenario}.txt', panel = panels, subtype = ['LUAD','LUSC'],scenario = scenarios,dataset=datasets),
        Table = expand('output/{dataset}/Tables/TableXX_NumberOfEvaluableSamples.txt',dataset=datasets),
    output:
        Evaluation_metrics = 'output/Tables/TableXX_Sensitivity_Specificity.txt',
        Boxplot_nMatches = 'output/Figures/FigureXX_Boxplot_nMatches.pdf',
        Piechart_SharedMutations = 'output/Figures/FigureXX_Piechart_Shared_Mutations.pdf',
        Piechart_EvaluablePatients = 'output/Figures/FigureXX_Piechart_Evaluable_patients.pdf',
        Barchart_performance = 'output/Figures/FigureXX_Barchart_Performance.pdf',
    conda:
        'envs/R.yaml'
    script:
        'scripts/Evaluate_Performance.R'

#++++++++++++++++++++++++++++++++++++++++++++++++ 4 MISC   +++++++++++++++++++++++++++++++++++++++++++++++++++++
# 4.1 Create table with number of evaluable tumors
rule Create_TableXX:
    input:
        AllMutations = lambda wildcards: 'data/TRACERx421_supplement/mutTableAll.cloneInfo.20220726.txt' if wildcards.dataset == 'TRACERx421' else 'data/TRACERx100_supplement/nejmoa1616288_appendix_2.xlsx' ,
        Mutations = lambda wildcards: expand('data/'+wildcards.dataset+'/{panel}/Selected_mutations_{subtype}.txt', panel = panels, subtype = subtypes),
        SampleOverview = lambda wildcards: 'data/TRACERx421_supplement/20221109_TRACERx421_all_patient_df.rds'  if wildcards.dataset == 'TRACERx421' else 'data/TRACERx100_supplement/nejmoa1616288_appendix_2.xlsx',
    output:
        Table = 'output/{dataset}/Tables/TableXX_NumberOfEvaluableSamples.txt',
    conda:
        'envs/R.yaml'
    script:
        'scripts/Create_TableXX.R'

    
# 4.2 Create excel sheets per panel with all metrics
rule Create_Clonality_tables:
    input:
        Clonality_metrics = lambda wildcards: expand('output/'+wildcards.dataset+'/'+wildcards.panel+'/Clonality_metrics_{subtype}_{scenario}.txt',subtype = subtypes, scenario= scenarios),
        Clonality_metrics_MC = lambda wildcards: expand('output/'+wildcards.dataset+'/'+wildcards.panel+'/Clonality_Calls_MC_{subtype}_{scenario}.txt',subtype = subtypes,scenario= scenarios)
    output:
        Table = 'output/{dataset}/Tables/Clonality_metrics_{panel}.xlsx',
    conda:
        'envs/R.yaml'
    script:
        'scripts/Create_Clonality_tables.R'
     
