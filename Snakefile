from scripts.utils import *
from glob import glob
configfile: "config.yaml"

#+++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS,TARGET AND FILES ++++++++++++++++++++++++++++++++++++++++           
# 0.1 Specify wildcards
(panels,) = glob_wildcards("manifest/{panel}.bed")

#----------------------------------------------------------------------------------------------------------------
# 0.2 specify target rule
rule all:
    input:
        expand( 'output/{panel}/Clonality_metrics.txt', panel = panels)
        
#++++++++++++++++++++++++++++++++++++++++++++++ 1 PREPROCESS DATA   +++++++++++++++++++++++++++++++++++++++++++++++++++
# 1.1 Fetch TRACERx mutations in panel regions
rule Filter_mutations:
    input:
        Mutations = 'data/TRACERx_supplement/nejmoa1616288_appendix_2.xlsx',
        panel = "manifest/{panel}.bed"
    output:
        Mutations = 'data/{panel}/Selected_mutations.txt'
    conda:
        'envs/R.yaml'
    script:
        'scripts/Filter_mutations.R'

        
#++++++++++++++++++++++++++++++++++++++++++++++++ 2 CLONALITY   +++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2.1 Call clonality
rule Call_Clonality:
    input:
        Mutations = 'data/{panel}/Selected_mutations.txt'
    output:
        Metrics = 'output/{panel}/Clonality_metrics.txt',
    conda:
        'envs/R.yaml'
    script:
        'scripts/Call_Clonality.R'

#++++++++++++++++++++++++++++++++++++++++++++++++ 3 MISC   +++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3.1 Create table with number of evaluable patients
rule Create_TableXX:
    input:
        Mutations = expand('data/{panel}/Selected_mutations.txt')
    output:
        Metrics = 'output/Tables/TableXX_NumberOfEvaluableSamples.txt',
    conda:
        'envs/R.yaml'
    script:
        'scripts/Create_TableXX.R'
