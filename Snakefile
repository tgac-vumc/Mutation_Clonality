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
        expand('data/{panel}/Selected_mutations.txt', panel = panels)
        
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

        
