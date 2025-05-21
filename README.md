# Mutation_Clonality

<p align="center">
  <img width="30%" height="30%" src="https://github.com/tgac-vumc/Mutation_Clonality/blob/main/dag.svg">
</p>

[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.esmoop.2025.105072-blue.svg)](https://doi.org/10.1016/j.esmoop.2025.105072)

Code to reproduce figures in manuscript "Performance and considerations in the use of diagnostic mutation panels for clonality testing in non-small-cell lung carcinoma" published in ESMO Open: Volume 10, Issue 5105072, May 2025.
The pipeline calculates the accuracy of several diagnostic mutation panels using public data from the TRACERx421 cohort. The user can add a custom mutation panel (in .bed format) to the manifest/
 to test the accuracy.
To run this pipeline Snakemake is required.

for easy installation you need (mini)conda.

Miniconda installation from folder where you want to install miniconda:

```
cd </path/to/files/dir/>
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

follow the instructions of the installation process, give the location where you want miniconda to be installed and answer YES to add miniconda to your path.

go to the directory where the analysis need to be performed

```
cd </path/to/analysis/dir>
git clone https://github.com/tgac-vumc/Mutation_clonality.git
cd Mutation_Clonality
```

install the snakemake environment,

```

conda env create --name snakemake --file envs/snakemake.yaml

```
activate the snakemake environment
```
source activate snakemake

```

Setup your files correctly: 
- Specify OncoKB access token in config.yaml   


navigate to Mutation_Clonality directory and start snakemake.

```
snakemake --use-conda

```
Useful snakemake options

-j , --cores, --jobs : Use at most N cores in parallel (default: 1). If N is omitted, the limit is set to the number of available cores.   
-n , --dryrun : Do not execute anything. but show rules which are planned to be performed.    
-k , --keep-going : Go on with independent jobs if a job fails.    
-f , --force : Force the execution of the selected target or the first rule regardless of already created output.   
-U , --until : Runs the pipeline until it reaches the specified rules or files. Only runs jobs that are dependencies of the specified rule or files, does not run sibling DAGs.   
-T , --timestamp : Add a timestamp to all logging output   
 
for all options go to http://snakemake.readthedocs.io/en/stable/executable.html#all-options
