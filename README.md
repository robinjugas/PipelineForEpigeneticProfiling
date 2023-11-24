# PipelineForEpigeneticProfiling
This repository contains the Snakemake pipeline for processing of the 
Dataset of Epigenetic profiles of pediatric oncology tumors. 

After downloading the pipeline scripts, 
place the sequencing reads into :
/PATH/bismark_pipeline/results/raw_reads/
Also, fill the the yaml config file and modify the shell script for running the Snakemake pipeline. 


## Snakemake environment
Snakemake workflow management can run inside [conda](https://docs.conda.io/en/latest/). Create the conda snakemake environemnt easily with conda:
```
conda create --name snakemake_environment -c bioconda -c conda-forge snakemake mamba
```
or setup the environment path:
```
conda create --prefix /path/snakemake_environment -c bioconda -c conda-forge snakemake mamba
```


## How to run the Bismark methylation pipelie
When you have filled the yaml configuration file and put the files in the folder structure you can start the Snakemake pipeline with this command:
```
snakemake --cores 12 --snakefile path_to_pipeline/Snakefile --directory path_to_DATA --configfile path/config.yaml --use-conda 
```

    --cores 12 - define number of threads
    --snakefile path_to_pipeline/Snakefile - define the path to the Bismark methylation Snakefile
    --directory path_to_DATA  - define the path to the data folder structure containing the results/raw_reads
    --configfile path/config.yaml  - define the path to the configuration  Snakefile

    --conda-frontend mamba -r   -  optional arguments, but recommended to use mamba for quicker package installation

Alternatively, you can run the example bash script **run_pipeline.sh**, which you can modify. 

**If you find some bugs, please fill the github bug report. I'll likely need example of the data to recreate the issue.**