#!/bin/bash
# Run snake file wih set config file, directory and copy files neccessary

# ## WORK DIR
WORK_DIR=`pwd`


## activate conda ENV
eval "$(conda shell.bash hook)"
conda activate snakemake

## RUN SNAKEMAKE

# run
snakemake --cores 12 --snakefile /pathToPipeline/Snakefile --directory $WORK_DIR --configfile $WORK_DIR/config.yaml --use-conda -r

conda deactivate
