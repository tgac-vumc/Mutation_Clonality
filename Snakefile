from scripts.utils import *
from glob import glob
configfile: "config.yaml"

#+++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS,TARGET AND FILES ++++++++++++++++++++++++++++++++++++++++           
# 0.1 Specify wildcards
(panels,) = glob_wildcards("manifest/{panel}.bed")
subtypes = ['NSCLC','LUAD','LUSC']
#----------------------------------------------------------------------------------------------------------------
# 0.2 specify target rule
rule all:
    input:
        expand( 'output/{panel}/Clonality_metrics_{subtype}.txt', panel = panels, subtype = subtypes),
	#'output/Tables/TableXX_Sensitivity_Specificity.txt',
        'output/Tables/TableXX_NumberOfEvaluableSamples.txt',
        'output/Figures/FigureXX_Barchart_Performance.pdf'
        
#++++++++++++++++++++++++++++++++++++++++++++++ 1 PREPROCESS DATA   +++++++++++++++++++++++++++++++++++++++++++++++++++
# 1.1 Fetch TRACERx mutations in panel regions
rule Filter_mutations:
    input:
        Mutations = 'data/TRACERx_supplement/nejmoa1616288_appendix_2.xlsx',
        panel = "manifest/{panel}.bed"
    output:
        Mutations = 'data/{panel}/Selected_mutations_{subtype}.txt'
    conda:
        'envs/R.yaml'
    script:
        'scripts/Filter_mutations.R'

        
#++++++++++++++++++++++++++++++++++++++++++++++++ 2 CLONALITY   +++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2.1 Call clonality
rule Call_Clonality:
    input:
        Mutations = 'data/{panel}/Selected_mutations_{subtype}.txt'
    output:
        Metrics = 'output/{panel}/Clonality_metrics_{subtype}.txt',
        Shared_Mutations = 'output/{panel}/Shared_mutations_{subtype}.txt'
    conda:
        'envs/R.yaml'
    script:
        'scripts/Call_Clonality.R'


#++++++++++++++++++++++++++++++++++++++++++++++++ 3 EVALUATION   +++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2.1 Call clonality
rule Evaluate_Performance:
    input:
        Shared_Mutations = expand('output/{panel}/Shared_mutations_{subtype}.txt', panel = panels, subtype = ['LUAD','LUSC']),
        Metrics = expand('output/{panel}/Clonality_metrics_{subtype}.txt', panel = panels, subtype = ['LUAD','LUSC']),
        Table = 'output/Tables/TableXX_NumberOfEvaluableSamples.txt',
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
# 3.1 Create table with number of evaluable patients
rule Create_TableXX:
    input:
        Mutations = expand('data/{panel}/Selected_mutations_{subtype}.txt', panel = panels, subtype = subtypes)
    output:
        Table = 'output/Tables/TableXX_NumberOfEvaluableSamples.txt',
    conda:
        'envs/R.yaml'
    script:
        'scripts/Create_TableXX.R'
