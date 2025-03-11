Snakefile: Contains preprocessing, deconvolution of bulk data, state discovery

config.yaml: Specifies the directories and parameters, may not need to change
* Make sure logs/slurm/ folder exists so that snakemake can run

Snakefile needs the envs folder to run.
Individual environments are automatically created for each rule when running Snakefile. 

In the scripts folder, cluster_sc.py, sc_cf.R, and survival.R are standalone files not incorporated into the snakemake pipeline.
1. cluster_sc.py: Classify the single cell data using state loadings from bulk data
2. sc_cf.R: Plotting true cell fractions of the single cell atlas
3. survival.R: Plotting the Kaplan-Meier curves of the TCGA STAD, enrichment analysis using Chi squared test

