from scripts.utils import *
from glob import glob
configfile: "config.yaml"

#+++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS,TARGET AND FILES ++++++++++++++++++++++++++++++++++++++++           
# 0.1 Specify wildcards
(panels,) = glob_wildcards("manifest/{panel}.bed")
subtypes = ['LUAD','LUSC']
#----------------------------------------------------------------------------------------------------------------
# 0.2 specify target rule
rule all:
    input:
        'output/Tables/TableXX_Evaluation_Metrics.txt',

#+++++++++++++++++++++++++++++++++++++++++++++++ 0 DOWNLOAD DATA   ++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1 Download supplementary data from Bakir/Frankell et al.
rule Download_data:
    output:
        Mutations_Frankell = 'data/TRACERx421_supplement_Frankell/20221109_TRACERx421_mutation_table.fst',
        Mutations_Bakir = 'data/TRACERx421_supplement_Bakir/mutTableAll.cloneInfo.20220726.txt',
        PatientOverview = 'data/TRACERx421_supplement_Frankell/20221109_TRACERx421_all_patient_df.rds',
        SampleOverview = 'data/TRACERx421_supplement_Frankell/20221109_TRACERx421_all_tumour_df.rds',
	SampleOverview_Bakir = 'data/TRACERx421_supplement_Bakir/sampleOverview.txt'
    shell:
        'scripts/Download_data.sh'

#++++++++++++++++++++++++++++++++++++++++++++++ 1 PREPROCESS DATA   +++++++++++++++++++++++++++++++++++++++++++++++++++
# 1.1 Subsample TRACERx mutations in panel regions
rule Filter_mutations:
    input:
        Mutations_Frankell = 'data/TRACERx421_supplement_Frankell/20221109_TRACERx421_mutation_table.fst',
        Mutations_Bakir = 'data/TRACERx421_supplement_Bakir/mutTableAll.cloneInfo.20220726.txt',
        PatientOverview = 'data/TRACERx421_supplement_Frankell/20221109_TRACERx421_all_patient_df.rds',
        panel = "manifest/{panel}.bed"
    output:
        Mutations = 'output/{panel}/Selected_mutations_{subtype}.txt'
    conda:
        'envs/R.yaml'
    script:
        'scripts/Filter_mutations.R'
        
#----------------------------------------------------------------------------------------------------------------
# 1.2 Annotate mutations with Oncogenic mutation status using OncoKB

rule Run_OncoKB:
    input:
        Mutations = 'output/KappaHyperExome/Selected_mutations_{subtype}.txt',
    output:
        MAF = temp('output/KappaHyperExome/MAF_{subtype}.txt'),
        Annotations = 'output/OncoKB_annotations_{subtype}.txt',
    conda:
        'envs/OncoKB.yaml'
    params:
        OncoKB_token = config['OncoKB']['access_token']
    shell:
        """
        python3 scripts/Create_MAF.py -i {input.Mutations} -o {output.MAF} -s {wildcards.subtype} ;
        python3 scripts/oncokb-annotator/MafAnnotator.py -i {output.MAF} -o {output.Annotations} -b {params.OncoKB_token}
        """

#----------------------------------------------------------------------------------------------------------------
# 1.2 Create Tumor pairs
rule Create_Tumor_Pairs:
    input:
        Mutations = 'output/KappaHyperExome/Selected_mutations_{subtype}.txt',
        PatientOverview = 'data/TRACERx421_supplement_Frankell/20221109_TRACERx421_all_patient_df.rds',
        SampleOverview = 'data/TRACERx421_supplement_Frankell/20221109_TRACERx421_all_tumour_df.rds',
        SampleOverview_Bakir = 'data/TRACERx421_supplement_Bakir/sampleOverview.txt',
    output:
        Tumor_pairs = 'data/Tumor_pairs_{subtype}.txt'
    conda:
        'envs/R.yaml'
    script:
        'scripts/Create_Tumor_Pairs.R'

#++++++++++++++++++++++++++++++++++++++++++++++++ 2 CLONALITY   +++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2.1 Call clonality with using two-metrics classifier
rule Call_Clonality_TwoMetric:
    input:
        Mutations = 'output/{panel}/Selected_mutations_{subtype}.txt',
        Tumor_pairs = 'data/Tumor_pairs_{subtype}.txt'
    output:
        Metrics = 'output/{panel}/Clonality_metrics_{subtype}.txt',
        Shared_Mutations = 'output/{panel}/Shared_mutations_{subtype}.txt'
    conda:
        'envs/R.yaml'
    script:
        'scripts/Call_Clonality_TwoMetric.R'

#----------------------------------------------------------------------------------------------------------------
# 2.2 Call clonality: Molecular classification algorithm
rule Call_Clonality_MC:
    input:
        Mutations = 'output/{panel}/Selected_mutations_{subtype}.txt',
        Annotations = 'output/OncoKB_annotations_{subtype}.txt',
        Tumor_pairs = 'data/Tumor_pairs_{subtype}.txt'
    output:
        Metrics = 'output/{panel}/MC_Clonality_metrics_{subtype}.txt'
    conda:
        'envs/R.yaml'
    script:
        'scripts/Call_Clonality_MCalgorithm.R'


#++++++++++++++++++++++++++++++++++++++++++++++++ 3 EVALUATION   +++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2.1 Call clonality
rule Evaluate_Performance:
    input:
        Metrics = expand('output/{panel}/Clonality_metrics_{subtype}.txt', panel = panels, subtype = ['LUAD','LUSC']),
        Metrics_MC = expand('output/{panel}/MC_Clonality_metrics_{subtype}.txt', panel = panels, subtype = ['LUAD','LUSC']),
    output:
        Evaluation_metrics = 'output/Tables/TableXX_Evaluation_Metrics.txt',
        Raw_Table = 'output/Tables/TableXX_Raw_Table.txt',
        Barchart_performance = 'output/Figures/FigureXX_Barchart_Performance.pdf',
        SankeyPlot = 'output/Figures/SankeyPlot.pdf'
    conda:
        'envs/R.yaml'
    script:
        'scripts/Evaluate_Performance.R'

