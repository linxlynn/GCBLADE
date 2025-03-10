configfile: "config.yaml"
from scripts.utils import *            
import os
from glob import glob

#+++++++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS AND TARGET ++++++++++++++++++++++++++++++++++++++++++++
# 0.1 Prepare wildcards and variables
data_dir = config["all"]["data_dir"]
log_dir = config["all"]["log_dir"]
datasets = major_datasets = config['all']['datasets']
expectation = ['with_expectation',]
Njob = config['BLADE']['Parameters']['Njob']
subtypes = ['STAD'] #, 'CIN', 'GS', 'EBV', 'MSI']
celltype_levels = ['cell_type2'] #,'cell_type']
hierarchical_levels = ['cell_type2']
genders = ['ALL']
# specify results folder
BLADE_run = ['results_consensus']
# specify fold (runs)
folds = list(range(1,2)) # list(range(1,11))
# specify GSE for single cell data
GSE_split = ['ALL']
# Specify cell types
cell_types = ["Lymphoid.cell", "Epithelial.cell", "Mast.cell", "Endocrine.cell",
              "Endothelial.cell", "Fibroblast", "Macrophage", "Myeloid.cell"]
cell_types2 = ["CD4..T.cell", "CD8..T.cell", "DC", "Endocrine.cell",
            "Endothelial.cell", "Epithelial.cell", "Macrophage.M1",
            "Macrophage.M2", "Mast.cell", "Memory.B.cell", "Naive.B.cell", "NK.cell",
            "Neutrophil", "Plasma.cell", "Treg", "iCAF", "myCAF"]

#-------------------------------------------------------------------------------------------------------------------
# 1.0 specify target rules, uncomment outputs that you want to generate
rule all:
    input:
        # expand(data_dir + "{major_datasets}/results/adata_phenotyped_detailed.h5ad", major_datasets = major_datasets),
        # expand(data_dir + "{major_datasets}/results/10X_Counts/", major_datasets = major_datasets),
        # expand(data_dir + "{major_dataset}/results/oncoBLADE/Signature_matrices_{celltype_level}.pickle", major_datasets = major_datasets, celltype_level = celltype_levels),
        # expand(data_dir + "{major_datasets}/results/oncoBLADE/NewSignature_{celltype_level}.pickle", major_datasets = major_datasets, celltype_level = celltype_levels),
        # expand(data_dir + "TCGA_STAD/Subtype.full.corrected.csv"),
        # expand(data_dir + "{major_datasets}/results/{celltype_level}_UMAP.png", major_datasets = major_datasets, celltype_level = celltype_levels),
        # expand(data_dir + "TCGA_STAD/raw_counts/bulk/{subtype}/{gender}/readcounts.txt", subtype = subtypes, gender = genders),
        # expand(data_dir + "TCGA_STAD/raw_counts/bulk/{subtype}/{gender}/samplesheet.txt", subtype = subtypes, gender = genders),
        # expand(data_dir + "{major_datasets}/results/oncoBLADE/{subtype}/{gender}/bulk.pickle", major_datasets = major_datasets, subtype = subtypes, gender = genders),
        # expand(data_dir + "{major_datasets}/simulated/Pseudobulk.pickle", major_datasets = major_datasets),
        # expand(data_dir + "{major_datasets}/results/oncoBLADE/True_profiles_{celltype_level}.pickle", major_datasets = major_datasets, celltype_level = celltype_levels),
        # expand(data_dir + "{major_datasets}/results/oncoBLADE/True_fractions_{celltype_level}.pickle", major_datasets = major_datasets, celltype_level = celltype_levels),
        # expand(data_dir + "{major_datasets}/results/oncoBLADE/{expectation}/AutoGeneS/{BLADE_run}/{celltype_level}_{subtype}_{gender}.oncoBLADE_output.k{fold}.pickle", major_datasets = major_datasets, expectation = expectation, BLADE_run = BLADE_run, celltype_level = celltype_levels, subtype = subtypes, gender = genders, fold = folds),
        # expand(data_dir + "{major_datasets}/results/oncoBLADE/{expectation}/AutoGeneS/{BLADE_run}/{hierarchical_level}_{subtype}_{gender}.hierarchical_oncoBLADE_output.k{fold}.pickle", major_datasets = major_datasets, expectation = expectation, BLADE_run = BLADE_run, hierarchical_level = hierarchical_levels, subtype = subtypes, gender = genders, fold = folds),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/{celltype_level}_{subtype}_{gender}.estimated_cf.k{fold}.csv", major_datasets = major_datasets, BLADE_run = BLADE_run, expectation = expectation, celltype_level = celltype_levels, subtype = subtypes, gender = genders, fold = folds),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/{hierarchical_level}_{subtype}_{gender}.hierarchical_estimated_cf.k{fold}.csv", major_datasets = major_datasets, BLADE_run = BLADE_run, expectation = expectation, hierarchical_level = hierarchical_levels, subtype = subtypes, gender = genders, fold = folds),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf/{celltype_level}_cf_simple_plot.k{fold}.png", major_datasets = major_datasets, BLADE_run = BLADE_run, expectation = expectation, celltype_level = celltype_levels,fold = folds),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf/{celltype_level}_cf_detailed_plot.k{fold}.png", major_datasets = major_datasets, BLADE_run = BLADE_run, expectation = expectation, celltype_level = celltype_levels,  fold = folds),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_norm/{celltype_level}_norm_simple_plot.k{fold}.png", major_datasets = major_datasets, BLADE_run = BLADE_run, expectation = expectation, celltype_level = celltype_levels,fold = folds),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_norm/{celltype_level}_norm_detailed_plot.k{fold}.png", major_datasets = major_datasets, BLADE_run = BLADE_run, expectation = expectation, celltype_level = celltype_levels,fold = folds),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf/{hierarchical_level}_cf_hierarchical_cumulative_plot.k{fold}.png", major_datasets = major_datasets, BLADE_run = BLADE_run, expectation = expectation, hierarchical_level = hierarchical_levels, fold = folds),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf/{hierarchical_level}_cf_hierarchical_detailed_plot.k{fold}.png", major_datasets = major_datasets, BLADE_run = BLADE_run, expectation = expectation, hierarchical_level = hierarchical_levels, fold = folds),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_norm/{hierarchical_level}_norm_hierarchical_cumulative_plot.k{fold}.png", major_datasets = major_datasets, BLADE_run = BLADE_run, expectation = expectation, hierarchical_level = hierarchical_levels, fold = folds),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_norm/{hierarchical_level}_norm_hierarchical_detailed_plot.k{fold}.png", major_datasets = major_datasets, BLADE_run = BLADE_run, expectation = expectation, hierarchical_level = hierarchical_levels, fold = folds),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_combined/{celltype_level}_cf_simple_plot.k{fold}.png", major_datasets = major_datasets, BLADE_run = BLADE_run, expectation = expectation, celltype_level = celltype_levels, fold = folds), 
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_combined/{celltype_level}_cf_detailed_plot.k{fold}.png", major_datasets = major_datasets, BLADE_run = BLADE_run, expectation = expectation, celltype_level = celltype_levels, fold = folds),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_combined_norm/{celltype_level}_norm_simple_plot.k{fold}.png", major_datasets = major_datasets, BLADE_run = BLADE_run, expectation = expectation, celltype_level = celltype_levels, fold = folds),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_combined_norm/{celltype_level}_norm_detailed_plot.k{fold}.png", major_datasets = major_datasets, BLADE_run = BLADE_run, expectation = expectation, celltype_level = celltype_levels, fold = folds),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_combined/{hierarchical_level}_cf_hierarchical_cumulative_plot.k{fold}.png", major_datasets = major_datasets, BLADE_run = BLADE_run, expectation = expectation, hierarchical_level = hierarchical_levels, fold = folds),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_combined/{hierarchical_level}_cf_hierarchical_detailed_plot.k{fold}.png", major_datasets = major_datasets, BLADE_run = BLADE_run, expectation = expectation, hierarchical_level = hierarchical_levels, fold = folds),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_combined_norm/{hierarchical_level}_norm_hierarchical_cumulative_plot.k{fold}.png", major_datasets = major_datasets, BLADE_run = BLADE_run, expectation = expectation, hierarchical_level = hierarchical_levels, fold = folds),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_combined_norm/{hierarchical_level}_norm_hierarchical_detailed_plot.k{fold}.png", major_datasets = major_datasets, BLADE_run = BLADE_run, expectation = expectation, hierarchical_level = hierarchical_levels, fold = folds),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_{subtype}_ALL.combined_estimated_cf.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels, subtype=subtypes),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_{subtype}_ALL.normalized_combined_estimated_cf.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels, subtype=subtypes),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_{subtype}_ALL.combined_metrics.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels, subtype=subtypes),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_{subtype}_ALL.hier_combined_estimated_cf.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, hier_level=hier_level, subtype=subtypes),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_{subtype}_ALL.hier_normalized_combined_estimated_cf.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, hier_level=hier_level, subtype=subtypes),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_{subtype}_ALL.hier_combined_metrics.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, hier_level=hier_level, subtype=subtypes),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_{subtype}_CIN.normalized_combined_estimated_cf.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels, subtype=subtypes),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_{subtype}_EBV.normalized_combined_estimated_cf.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels, subtype=subtypes),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_{subtype}_GS.normalized_combined_estimated_cf.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels, subtype=subtypes),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_{subtype}_MSI.normalized_combined_estimated_cf.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels, subtype=subtypes),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_{subtype}_CIN.combined_estimated_cf.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels, subtype=subtypes),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_{subtype}_EBV.combined_estimated_cf.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels, subtype=subtypes),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_{subtype}_GS.combined_estimated_cf.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels, subtype=subtypes),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_{subtype}_MSI.combined_estimated_cf.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels, subtype=subtypes),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_{subtype}_CIN.hier_normalized_combined_estimated_cf.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, hier_level=hier_level, subtype=subtypes),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_{subtype}_EBV.hier_normalized_combined_estimated_cf.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, hier_level=hier_level, subtype=subtypes),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_{subtype}_GS.hier_normalized_combined_estimated_cf.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, hier_level=hier_level, subtype=subtypes),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_{subtype}_MSI.hier_normalized_combined_estimated_cf.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, hier_level=hier_level, subtype=subtypes),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_{subtype}_CIN.hier_combined_estimated_cf.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, hier_level=hier_level, subtype=subtypes),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_{subtype}_EBV.hier_combined_estimated_cf.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, hier_level=hier_level, subtype=subtypes),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_{subtype}_GS.hier_combined_estimated_cf.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, hier_level=hier_level, subtype=subtypes),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_{subtype}_MSI.hier_combined_estimated_cf.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, hier_level=hier_level, subtype=subtypes),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_subtypes_consistency", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_subtypes_consistency/results.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_subtypes_consistency_STAD", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_subtypes_consistency_STAD/results.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_hier_subtypes_consistency_STAD", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, hier_level=hier_level),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_hier_subtypes_consistency_STAD/results.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, hier_level=hier_level),
        # expand(data_dir + "{major_datasets}/results/oncoBLADE/{expectation}/AutoGeneS/{BLADE_run}/{celltype_level}_STAD.purification.k3.pickle", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels),
        # expand(data_dir + "{major_datasets}/results/oncoBLADE/{expectation}/AutoGeneS/{BLADE_run}/{celltype_level}_STAD.hier_purification.k6.pickle", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/state/{celltype_level}_StateLoadings.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/state/{celltype_level}_StateScores.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/state/{celltype_level}_cophcor.png", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/state/{celltype_level}_heatmap.pdf", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/state/{celltype_level}_barplots/", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/state/{celltype_level}_samples_assigned_states_specK2.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels),
        # expand(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/state/{celltype_level}_genes_assigned_states_specK2.csv", major_datasets=major_datasets, BLADE_run=BLADE_run, expectation=expectation, celltype_level=celltype_levels)
         


#+++++++++++++++++++++++++++++++++++++++++++++++ 2 PREPROCESS DATA  ++++++++++++++++++++++++++++++++++++++++++++++++
# 2.1 Convert Seurat to Adata
# convert seurat object to an Adata object with iCAF1/2/3 into grouped iCAF & myCAF1/2 into myCAF
rule Convert_Seurat:
   input:
       seurat_object = data_dir + "{major_datasets}/results/GC_seuratObj.rda"
   output:
       adata_phenotyped = data_dir + "{major_datasets}/results/adata_phenotyped_detailed.h5ad",
   conda:
       "envs/Convert_Seurat.yaml"
   shell:
       """
       Rscript scripts/Convert_Seurat_to_Adata.R \
       -i {input.seurat_object} \
       -o {output.adata_phenotyped} \
       """
	
# 2.2 Convert Seurat to Matrix
1# convert phenotyped Seurat object to a matrix.mtx.gz object
rule Make_matrix:
   input:
       seurat_object =  data_dir + "{major_datasets}/results/GC_seuratObj.rda"
   output:
       output = directory(data_dir + "{major_datasets}/results/10X_Counts/")
   conda:
       "envs/Make_mtx.yaml"
   shell:
       """
       Rscript scripts/Make_matrix.R \
       --input {input.seurat_object} \
       --output {output.output} \
       """


#+++++++++++++++++++++++++++++++++++++++++++++++ 3 CELL PHENOTYPING  +++++++++++++++++++++++++++++++++++++++++++++++

# Create UMAP for single cell data with both major and minor cell types
rule UMAP_phenotyping:
    input:
        seurat_object =  data_dir + "{major_datasets}/results/GC_seuratObj.rda"
    output:
        UMAP_plot = data_dir + "{major_datasets}/results/{celltype_level}_UMAP.png"
    conda:
        "envs/umap_env.yaml"
    shell:
       """
       Rscript scripts/Create_UMAP.R \
       --input_file {input.seurat_object} \
       --cell_level {wildcards.celltype_level} \
       --output_plot {output.UMAP_plot} \
       """
        

#++++++++++++++++++++++++++++++++++++++++++ 4 DIFFERENTIAL GENE EXPRESSION +++++++++++++++++++++++++++++++++++++++++
# 4.1 Perform differential gene expression analysis for different celltypes with AutoGeneS
rule AutoGeneS:
  input:
      adata_final = data_dir + "{major_dataset}/results/adata_phenotyped_detailed.h5ad"
  output:
      DEGs = data_dir + "{major_dataset}/results/AutoGeneS/AutoGeneS_standard_genes_{celltype_level}_run6.txt"
  conda:
      "envs/AutoGeneS.yaml"
  shell:
      """
      python3 scripts/AutoGeneS.py \
      -i {input.adata_final} \
      -level {wildcards.celltype_level} \
      -o {output.DEGs} 
      """
	
#+++++++++++++++++++++++++++++++++++++++++++++++ 5 CREATE SIGNATURE  +++++++++++++++++++++++++++++++++++++++++++++++
# 5.1 Create cell type specific Signature Matrix
rule Create_Signature:
   input:
       adata_final = data_dir + "{major_dataset}/results/adata_phenotyped_detailed.h5ad",
   output:
       signature_matrices =  data_dir + "{major_dataset}/results/oncoBLADE/Signature_matrices_{celltype_level}.pickle"
   container:
       "docker://gcfntnu/scanpy:1.7.0"
   shell:
       """
       python3 scripts/Create_Signature.py \
       -i {input.adata_final} \
       -celltype {wildcards.celltype_level} \
       -o {output.signature_matrices} 
       """

# 5.2 Correct variance for mean dependence of stddev in lognormal data
rule Correct_variance:
   input:
       signature_matrices = data_dir + "{major_datasets}/results/oncoBLADE/Signature_matrices_{celltype_level}.pickle"
   output:
       signature_matrices_corrected =  data_dir + "{major_datasets}/results/oncoBLADE/NewSignature_{celltype_level}.pickle"
   conda:
       "envs/scran.yaml"
   shell:
       """
       Rscript scripts/Correct_variance.R \
       -i {input.signature_matrices} \
       -o {output.signature_matrices_corrected} \
       """        	

#++++++++++++++++++++++++++++++++++++++++++++++ 6 Preprocess bulk RNAseq +++++++++++++++++++++++++++++++++++++++++++++++++++
# 6.0 Define TCGA molecular subtype per patient
rule Define_subtype:
   input:
       subtypes_sheet = data_dir + "TCGA_STAD/Subtype.full.txt",
       auxiliary = data_dir + "TCGA_STAD/nationwidechildrens.org_auxiliary_stad.txt",
       patients = data_dir + "TCGA_STAD/nationwidechildrens.org_clinical_patient_stad.txt"
   output:
       subtypes_corrected = data_dir + "TCGA_STAD/Subtype.full.corrected.csv"
   conda:
       "envs/mrna_env.yaml"
   shell:
       """
       Rscript scripts/Define_subtypes.R \
       --subtypes {input.subtypes_sheet} \
       --auxiliary {input.auxiliary} \
       --patients {input.patients} \
       --output {output.subtypes_corrected} \
       """
 	
# 6.1 Make readcounts file from RNAseq data from the TCGA_STAD
rule make_readcounts:
    input:
        patient_data = data_dir + "TCGA_STAD/nationwidechildrens.org_clinical_patient_stad.txt", 
        subtype_sheet = data_dir + "TCGA_STAD/Subtype.full.corrected.csv",
        gdc_samplesheet = data_dir + "TCGA_STAD/raw_counts/RNAseq_data/gdc_sample_sheet.2024-03-01.tsv",
        ACE_samplesheet = data_dir + "TCGA_STAD/raw_counts/bulk/ACE_samplesheet.txt"
    output:
        readcounts = data_dir + "TCGA_STAD/raw_counts/bulk/{subtype}/{gender}/readcounts.txt",
        tumor_fractions = data_dir + "TCGA_STAD/raw_counts/bulk/{subtype}/{gender}/samplesheet.txt"
    params:
        dir_name = data_dir + "TCGA_STAD/raw_counts/RNAseq_data",
    conda:
        "envs/mrna_env.yaml"
    shell:
        """
        Rscript scripts/make_readcounts.R \
        -i {input.patient_data} \
        --subtype_sheet {input.subtype_sheet} \
        --gdc_samplesheet {input.gdc_samplesheet} \
        --ACE_samplesheet {input.ACE_samplesheet} \
        --input_dir {params.dir_name} \
        --subtype {wildcards.subtype} \
        --gender {wildcards.gender} \
        -o {output.readcounts} \
        --tumor_fractions {output.tumor_fractions} \
        """


# 6.2 Preprocess bulk transcriptome data for use in BLADE
rule preprocess_bulk_RNAseq:
   input:
       readcounts = data_dir + "TCGA_STAD/raw_counts/bulk/{subtype}/{gender}/readcounts.txt"
   output:
       bulk = data_dir + "{major_datasets}/results/oncoBLADE/{subtype}/{gender}/bulk.pickle"
   params:
       cache_dir = ".cache/{major_datasets}/",
   conda:
       "envs/preprocess_env.yaml"
   shell:
       """
       Rscript scripts/Preprocess_bulk_RNAseq.R \
       -i {input.readcounts} \
       -o {output.bulk} \
       """
	
#+++++++++++++++++++++++++++++++++++ 6.3 Simulated bulk data from single cell data  +++++++++++++++++++++++++++++++++++++
# 6.3.1 Create simulated bulk transcriptome data

rule Simulate_bulk:
   input:
       adata_final = data_dir + "{major_datasets}/results/adata_phenotyped_detailed.h5ad", 
       matrix = data_dir + "{major_datasets}/results/10X_Counts/matrix.mtx.gz",
   output:
       simulated_bulk =  data_dir + "{major_datasets}/simulated/Pseudobulk.pickle"
   params:
       cache_dir = ".cache/{major_datasets}/",
       input_dir = data_dir + "{major_datasets}/results/10X_Counts/"
   container:
       "docker://gcfntnu/scanpy:1.7.0"
   shell:
       """
       python3 scripts/Simulate_bulk.py \
       -i {input.adata_final} \
       -d {params.input_dir} \
       -o {output.simulated_bulk} \
       -cache_dir {params.cache_dir}
       """
     
# 6.3.2 Collect true profiles

rule Gather_True_Values:
   input:
       adata_final = data_dir + "{major_datasets}/results/adata_phenotyped_detailed.h5ad",
   output:
       True_profiles =  data_dir + "{major_datasets}/results/oncoBLADE/True_profiles_{celltype_level}.pickle",
       True_fractions =  data_dir + "{major_datasets}/results/oncoBLADE/True_fractions_{celltype_level}.pickle"
   container:
       "docker://gcfntnu/scanpy:1.7.0"
   shell:
       """
       python3 scripts/Gather_True_values.py \
       -i {input.adata_final} \
       -celltypes {wildcards.celltype_level} \
       -o_profiles {output.True_profiles} \
       -o_fractions {output.True_fractions}
       """
	    
#++++++++++++++++++++++++++++++++++++++++++++++++ 7 DECONVOLUTION +++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++ 7.1 Deconvolution on TCGA data +++++++++++++++++++++++++++++++++++++++
# 7.1.1 Run oncoBLADE with AutoGeneS
rule oncoBLADE:
    input:
        bulk = data_dir + "{major_datasets}/results/oncoBLADE/{subtype}/{gender}/bulk.pickle",
        signature_matrices_corrected =  data_dir + "{major_datasets}/results/oncoBLADE/NewSignature_{celltype_level}.pickle",
        tumor_fractions = data_dir + "TCGA_STAD/raw_counts/bulk/{subtype}/{gender}/samplesheet.txt",
        DEGs_data  = data_dir +  "{major_datasets}/results/AutoGeneS/AutoGeneS_standard_genes_{celltype_level}_run6.txt"
    output:
        oncoBLADE_output = data_dir + "{major_datasets}/results/oncoBLADE/{expectation}/AutoGeneS/{BLADE_run}/{celltype_level}_{subtype}_{gender}.oncoBLADE_output.k{fold}.pickle"
    params:
        Parameters= "config.yaml",
        cohort = "TCGA"
    threads: Njob
    conda:
        "envs/BLADE_env.yaml"
    shell:
        """
        python3 scripts/BLADE_AutoGeneS.py  \
        -i_bulk {input.bulk} \
        -i_sign {input.signature_matrices_corrected} \
	-i_TF {input.tumor_fractions} \
        -i_autogenes {input.DEGs_data} \
	-expectation {wildcards.expectation}\
        -parameters {params.Parameters} \
        -cohort {params.cohort} \
        -o {output.oncoBLADE_output} \
        """

# 7.1.2 Run Hierarchical oncoBLADE
rule Hierarchical_oncoBLADE:
    input:
        bulk = data_dir + "{major_datasets}/results/oncoBLADE/{subtype}/{gender}/bulk.pickle",
        signature_matrices_corrected2 =  data_dir + "{major_datasets}/results/oncoBLADE/NewSignature_cell_type.pickle",
        signature_matrices_corrected =  data_dir + "{major_datasets}/results/oncoBLADE/NewSignature_cell_type2.pickle",
        DEGs_data  = data_dir + '{major_datasets}/results/AutoGeneS/AutoGeneS_standard_genes_{hierarchical_level}_run6.txt',
        oncoBLADE_output = data_dir + "{major_datasets}/results/oncoBLADE/{expectation}/AutoGeneS/{BLADE_run}/cell_type_{subtype}_{gender}.oncoBLADE_output.k{fold}.pickle"
    output:
        hierarchical_oncoBLADE_output = data_dir + "{major_datasets}/results/oncoBLADE/{expectation}/AutoGeneS/{BLADE_run}/{hierarchical_level}_{subtype}_{gender}.hierarchical_oncoBLADE_output.k{fold}.pickle"
    params:
        Parameters= "config.yaml",
        cohort = "TCGA"
    threads: Njob
    conda:
        "envs/BLADE_env.yaml"
    shell:
        """
        python3 scripts/Hierarchical_oncoBLADE.py  \
        -i_bulk {input.bulk} \
        -i_sign2 {input.signature_matrices_corrected2} \
        -i_sign {input.signature_matrices_corrected} \
        -i_DEGs {input.DEGs_data} \
        -i_oncoBLADE {input.oncoBLADE_output} \
        -parameters {params.Parameters} \
        -expectation {wildcards.expectation} \
        -hierarchical_level {wildcards.hierarchical_level} \
        -cohort {params.cohort} \
        -o {output.hierarchical_oncoBLADE_output}
        """     

#++++++++++++++++++++++++++++++++++++++++++++++++ 8 EVALUATION +++++++++++++++++++++++++++++++++++++++++++++++++
# 8.1 Evaluate results from oncoBLADE AutoGeneS
rule Evaluate_results:
    input:
        oncoBLADE_output = data_dir + "{major_datasets}/results/oncoBLADE/{expectation}/AutoGeneS/{BLADE_run}/{celltype_level}_{subtype}_{gender}.oncoBLADE_output.k{fold}.pickle", 
        samplesheet = data_dir + "TCGA_STAD/raw_counts/bulk/{subtype}/{gender}/samplesheet.txt",
        signature_matrices = data_dir + "{major_datasets}/results/oncoBLADE/NewSignature_{celltype_level}.pickle"
    output:
        df_estimated_cf = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/{celltype_level}_{subtype}_{gender}.estimated_cf.k{fold}.csv"
    shell:
        """
        python3 scripts/Evaluate_results.py \
        -i {input.oncoBLADE_output} \
        --samples {input.samplesheet} \
        -i_sign {input.signature_matrices} \
        --celltype {wildcards.celltype_level} \
        --subtype {wildcards.subtype} \
        --gender {wildcards.gender} \
        -o {output.df_estimated_cf} \
        """

# 8.2 Evaluate results from Hierarchical_oncoBLADE 
rule Evaluate_hierarchical_results:
    input:
        hierarchical_oncoBLADE_output = data_dir + "{major_datasets}/results/oncoBLADE/{expectation}/AutoGeneS/{BLADE_run}/{hierarchical_level}_{subtype}_{gender}.hierarchical_oncoBLADE_output.k{fold}.pickle",
        samplesheet = data_dir + "TCGA_STAD/raw_counts/bulk/{subtype}/{gender}/samplesheet.txt",
        signature_matrices = data_dir + "{major_datasets}/results/oncoBLADE/NewSignature_cell_type2.pickle"
    output:
        df_estimated_cf = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/{hierarchical_level}_{subtype}_{gender}.hierarchical_estimated_cf.k{fold}.csv"
    params:
        Parameters = "hierarchical"
    shell:
        """
        python3 scripts/Evaluate_results.py \
        -i {input.hierarchical_oncoBLADE_output} \
        --samples {input.samplesheet} \
        -i_sign {input.signature_matrices} \
        --celltype {params.Parameters} \
        --subtype {wildcards.subtype} \
        --gender {wildcards.gender} \
        -o {output.df_estimated_cf} \
        """
#++++++++++++++++++++++++++++++++++++++++++++++++ 9 COMPARING TCGA CELL FRACTIONS +++++++++++++++++++++++++++++++++++++	
#++++++++++++++++++++++++++++++++++++++  COMPARING CELL FRACTIONS BETWEEN SUBTYPES ++++++++++++++++++++++++++++++++
# 9.1 Compare the cell fractions between subtypes (subtypes given separately)
rule Compare_CF_subtype:
    input:
        CIN_cell_fractions = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/{celltype_level}_CIN_ALL.estimated_cf.k{fold}.csv",
        EBV_cell_fractions = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/{celltype_level}_EBV_ALL.estimated_cf.k{fold}.csv",
        GS_cell_fractions = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/{celltype_level}_GS_ALL.estimated_cf.k{fold}.csv",
        MSI_cell_fractions = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/{celltype_level}_MSI_ALL.estimated_cf.k{fold}.csv"
    output:
        cf_plot_simple = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf/{celltype_level}_cf_simple_plot.k{fold}.png",
        cf_plot_detailed = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf/{celltype_level}_cf_detailed_plot.k{fold}.png",
        norm_plot_simple = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_norm/{celltype_level}_norm_simple_plot.k{fold}.png",
        norm_plot_detailed = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_norm/{celltype_level}_norm_detailed_plot.k{fold}.png"
    conda:
        "envs/compare_cf.yaml"
    shell:
        """
        Rscript scripts/Compare_subtype_cf.R \
        --CIN_cf {input.CIN_cell_fractions} \
        --EBV_cf {input.EBV_cell_fractions} \
        --GS_cf {input.GS_cell_fractions} \
        --MSI_cf {input.MSI_cell_fractions} \
        --cell_level {wildcards.celltype_level} \
        --output_simple {output.cf_plot_simple} \
        --output_detailed {output.cf_plot_detailed} \
        --output_norm_simple {output.norm_plot_simple} \
        --output_norm_detailed {output.norm_plot_detailed} \
        """

# 9.2 Compare cell fractions between subtypes for hierarchical model (subtypes given separately)
rule Compare_CF_subtype_hierarchical:
    input:   
        CIN_cell_fractions = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/{hierarchical_level}_CIN_ALL.hierarchical_estimated_cf.k{fold}.csv",
        EBV_cell_fractions = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/{hierarchical_level}_EBV_ALL.hierarchical_estimated_cf.k{fold}.csv",
        GS_cell_fractions = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/{hierarchical_level}_GS_ALL.hierarchical_estimated_cf.k{fold}.csv",
        MSI_cell_fractions = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/{hierarchical_level}_MSI_ALL.hierarchical_estimated_cf.k{fold}.csv"
    output:
        cf_plot_simple = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf/{hierarchical_level}_cf_hierarchical_cumulative_plot.k{fold}.png",
        cf_plot_detailed = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf/{hierarchical_level}_cf_hierarchical_detailed_plot.k{fold}.png",
        norm_plot_simple = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_norm/{hierarchical_level}_norm_hierarchical_cumulative_plot.k{fold}.png",
        norm_plot_detailed = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_norm/{hierarchical_level}_norm_hierarchical_detailed_plot.k{fold}.png"
    params:
        hierarchical = "YES"
    conda:
        "envs/compare_cf.yaml"
    shell:
        """
        Rscript scripts/Compare_subtype_cf_hierarchical.R \
        --CIN_cf {input.CIN_cell_fractions} \
        --EBV_cf {input.EBV_cell_fractions} \
        --GS_cf {input.GS_cell_fractions} \
        --MSI_cf {input.MSI_cell_fractions} \
        --cell_level {wildcards.hierarchical_level} \
        --hierarchical {params.hierarchical} \
        --output_simple {output.cf_plot_simple} \
        --output_detailed {output.cf_plot_detailed} \
        --output_norm_simple {output.norm_plot_simple} \
        --output_norm_detailed {output.norm_plot_detailed} \
        """   
       
# 9.3 Compare CF b/t subtypes (data given combined)
rule Compare_CF_subtype_combined:
    input:
        STAD_cell_fractions = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/{celltype_level}_STAD_ALL.estimated_cf.k{fold}.csv",
        subtype_corrected = data_dir + "TCGA_STAD/Subtype.full.corrected.csv"
    output:
        cf_plot_simple = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_combined/{celltype_level}_cf_simple_plot.k{fold}.png",
        cf_plot_detailed =data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_combined/{celltype_level}_cf_detailed_plot.k{fold}.png",
        norm_plot_simple = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_combined_norm/{celltype_level}_norm_simple_plot.k{fold}.png",
        norm_plot_detailed = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_combined_norm/{celltype_level}_norm_detailed_plot.k{fold}.png"
    conda:
        "envs/compare_cf.yaml"
    shell:
        """
        Rscript scripts/Compare_subtype_combined.R \
        --STAD_cf {input.STAD_cell_fractions} \
        --subtype_full {input.subtype_corrected} \
        --cell_level {wildcards.celltype_level} \
        --output_simple {output.cf_plot_simple} \
        --output_detailed {output.cf_plot_detailed} \
        --output_norm_simple {output.norm_plot_simple} \
        --output_norm_detailed {output.norm_plot_detailed} \
        """

# 9.4 Compare CF b/t subtypes for hierarchical model (data given combined)
rule Compare_CF_subtype_hierarchical_combined:
    input:   
        STAD_cell_fractions = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/{hierarchical_level}_STAD_ALL.hierarchical_estimated_cf.k{fold}.csv",
        subtype_corrected = data_dir + "TCGA_STAD/Subtype.full.corrected.csv"
    output:
        cf_plot_simple = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_combined/{hierarchical_level}_cf_hierarchical_cumulative_plot.k{fold}.png",
        cf_plot_detailed = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_combined/{hierarchical_level}_cf_hierarchical_detailed_plot.k{fold}.png",
        norm_plot_simple = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_combined_norm/{hierarchical_level}_norm_hierarchical_cumulative_plot.k{fold}.png",
        norm_plot_detailed = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/compare_cf_combined_norm/{hierarchical_level}_norm_hierarchical_detailed_plot.k{fold}.png"
    params:
        hierarchical = "YES"
    conda:
        "envs/compare_cf.yaml"
    shell:
        """
        Rscript scripts/Compare_subtype_combined_hier.R \
        --STAD_cf {input.STAD_cell_fractions} \
        --subtype_full {input.subtype_corrected} \
        --cell_level {wildcards.hierarchical_level} \
        --hierarchical {params.hierarchical} \
        --output_simple {output.cf_plot_simple} \
        --output_detailed {output.cf_plot_detailed} \
        --output_norm_simple {output.norm_plot_simple} \
        --output_norm_detailed {output.norm_plot_detailed} \
        """

#++++++++++++++++++++++++++++++++++++++++ 10 CONSENSUS +++++++++++++++++++++++++++++++++++++
# 10.1 Combine multiple runs into one file
rule combine_kfolds:
    input:
        cf_files = lambda wildcards: [f"{data_dir}{wildcards.major_datasets}/results/evaluate/{wildcards.BLADE_run}/{wildcards.expectation}/{wildcards.celltype_level}_{wildcards.subtype}_ALL.estimated_cf.k{fold}.csv" for fold in folds],
        metrics = lambda wildcards: [f"{data_dir}{wildcards.major_datasets}/results/evaluate/{wildcards.BLADE_run}/{wildcards.expectation}/{wildcards.celltype_level}_{wildcards.subtype}_ALL.evaluation_metrics.k{fold}.csv" for fold in folds],
        subtype_full = data_dir + "TCGA_STAD/Subtype.full.corrected.csv"
    output:
        comb_data = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_{subtype}_ALL.combined_estimated_cf.csv",
        norm_data = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_{subtype}_ALL.normalized_combined_estimated_cf.csv",
        obj_data = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_{subtype}_ALL.combined_metrics.csv"
    conda:
        "envs/combine_kfold.yaml"
    params:
        if_STAD = "YES"
    shell:
        """
        Rscript scripts/Combine_k_folds.R \
        --cf_file {input.cf_files}\
        --metrics {input.metrics} \
        --if_STAD {params.if_STAD} \
        --subtype_file {input.subtype_full} \
        --output_comb {output.comb_data} \
        --output_norm {output.norm_data} \
        --output_obj {output.obj_data} \
        """

rule combine_kfolds_hier:
    input:
        cf_files = lambda wildcards: [f"{data_dir}{wildcards.major_datasets}/results/evaluate/{wildcards.BLADE_run}/{wildcards.expectation}/{wildcards.hier_level}_{wildcards.subtype}_ALL.hierarchical_estimated_cf.k{fold}.csv" for fold in folds],
        metrics = lambda wildcards: [f"{data_dir}{wildcards.major_datasets}/results/evaluate/{wildcards.BLADE_run}/{wildcards.expectation}/{wildcards.hier_level}_{wildcards.subtype}_ALL.hierarchical_evaluation_metrics.k{fold}.csv" for fold in folds],
        subtype_full = data_dir + "TCGA_STAD/Subtype.full.corrected.csv"
    output:
        comb_data = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_{subtype}_ALL.hier_combined_estimated_cf.csv",
        norm_data = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_{subtype}_ALL.hier_normalized_combined_estimated_cf.csv",
        obj_data = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_{subtype}_ALL.hier_combined_metrics.csv"
    conda:
        "envs/combine_kfold.yaml"
    params:
        if_STAD = "YES"
    shell:
        """
        Rscript scripts/Combine_k_folds.R \
        --cf_file {input.cf_files}\
        --metrics {input.metrics} \
        --if_STAD {params.if_STAD} \
        --subtype_file {input.subtype_full} \
        --output_comb {output.comb_data} \
        --output_norm {output.norm_data} \
        --output_obj {output.obj_data} \
        """

# 10.2 Divide STAD files (data deconvoluted together) into subtypes
rule divide_STAD_combined:
    input:
        # STAD_comb = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_STAD_ALL.normalized_combined_estimated_cf.csv" # normalized
        STAD_comb = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_STAD_ALL.combined_estimated_cf.csv" # original
    output:
        # subtype_CIN = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_{subtype}_CIN.normalized_combined_estimated_cf.csv",
        # subtype_EBV = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_{subtype}_EBV.normalized_combined_estimated_cf.csv",
        # subtype_GS = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_{subtype}_GS.normalized_combined_estimated_cf.csv",
        # subtype_MSI = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_{subtype}_MSI.normalized_combined_estimated_cf.csv"
        subtype_CIN = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_{subtype}_CIN.combined_estimated_cf.csv",
        subtype_EBV = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_{subtype}_EBV.combined_estimated_cf.csv",
        subtype_GS = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_{subtype}_GS.combined_estimated_cf.csv",
        subtype_MSI = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_{subtype}_MSI.combined_estimated_cf.csv"
    params:
        if_normalized = "NO"
    shell:
        """
        Rscript scripts/divide_STAD_combined.R \
        --STAD_cf {input.STAD_comb} \
        --if_normalized {params.if_normalized} \
        --cell_level {wildcards.celltype_level} \
        --CIN_cf {output.subtype_CIN} \
        --EBV_cf {output.subtype_EBV} \
        --GS_cf {output.subtype_GS} \
        --MSI_cf {output.subtype_MSI} \
        """

rule divide_STAD_combined_hier:
    input:
        # STAD_comb = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_STAD_ALL.hier_normalized_combined_estimated_cf.csv" # normalized
        STAD_comb = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_STAD_ALL.hier_combined_estimated_cf.csv" # original
    output:
        # subtype_CIN = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_{subtype}_CIN.hier_normalized_combined_estimated_cf.csv",
        # subtype_EBV = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_{subtype}_EBV.hier_normalized_combined_estimated_cf.csv",
        # subtype_GS = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_{subtype}_GS.hier_normalized_combined_estimated_cf.csv",
        # subtype_MSI = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_{subtype}_MSI.hier_normalized_combined_estimated_cf.csv"
        subtype_CIN = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_{subtype}_CIN.hier_combined_estimated_cf.csv",
        subtype_EBV = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_{subtype}_EBV.hier_combined_estimated_cf.csv",
        subtype_GS = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_{subtype}_GS.hier_combined_estimated_cf.csv",
        subtype_MSI = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_{subtype}_MSI.hier_combined_estimated_cf.csv"
    params:
        if_normalized = "NO"
    shell:
        """
        Rscript scripts/divide_STAD_combined.R \
        --STAD_cf {input.STAD_comb} \
        --if_normalized {params.if_normalized} \
        --cell_level {wildcards.hier_level} \
        --CIN_cf {output.subtype_CIN} \
        --EBV_cf {output.subtype_EBV} \
        --GS_cf {output.subtype_GS} \
        --MSI_cf {output.subtype_MSI} \
        """

# 10.3 Subtype difference using Cohen's d        
rule subtype_difference:
    input:
        # CIN_cf = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_CIN_ALL.normalized_combined_estimated_cf.csv", # subtypes files deconvoluted independently
        # EBV_cf = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_EBV_ALL.normalized_combined_estimated_cf.csv",
        # GS_cf = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_GS_ALL.normalized_combined_estimated_cf.csv",
        # MSI_cf = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_MSI_ALL.normalized_combined_estimated_cf.csv"
        CIN_cf = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_STAD_CIN.normalized_combined_estimated_cf.csv", # subtype files from STAD
        EBV_cf = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_STAD_EBV.normalized_combined_estimated_cf.csv",
        GS_cf = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_STAD_GS.normalized_combined_estimated_cf.csv",
        MSI_cf = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_STAD_MSI.normalized_combined_estimated_cf.csv"
    output:
        # result = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_subtypes_consistency/results.csv",
        # plots = directory(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_subtypes_consistency")
        result = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_subtypes_consistency_STAD/results.csv",
        plots = directory(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{celltype_level}_subtypes_consistency_STAD")
    conda:
        "envs/analyze_difference.yaml"
    shell:
        """
        Rscript scripts/subtype_difference.R \
        --CIN_cf {input.CIN_cf} \
        --EBV_cf {input.EBV_cf} \
        --GS_cf {input.GS_cf} \
        --MSI_cf {input.MSI_cf} \
        --result {output.result} \
        --plots_dir {output.plots} \
        """

rule subtype_difference_hier:
    input:
        CIN_cf = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_STAD_CIN.hier_normalized_combined_estimated_cf.csv", # normalized 
        EBV_cf = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_STAD_EBV.hier_normalized_combined_estimated_cf.csv",
        GS_cf = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_STAD_GS.hier_normalized_combined_estimated_cf.csv",
        MSI_cf = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_STAD_MSI.hier_normalized_combined_estimated_cf.csv"
        # CIN_cf = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_STAD_CIN.hier_combined_estimated_cf.csv", # original
        # EBV_cf = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_STAD_EBV.hier_combined_estimated_cf.csv",
        # GS_cf = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_STAD_GS.hier_combined_estimated_cf.csv",
        # MSI_cf = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_STAD_MSI.hier_combined_estimated_cf.csv"
    output:
        result = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_hier_subtypes_consistency_STAD/results.csv",
        plots = directory(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/Consensus/{hier_level}_hier_subtypes_consistency_STAD")
    conda:
        "envs/analyze_difference.yaml"
    shell:
        """
        Rscript scripts/subtype_difference.R \
        --CIN_cf {input.CIN_cf} \
        --EBV_cf {input.EBV_cf} \
        --GS_cf {input.GS_cf} \
        --MSI_cf {input.MSI_cf} \
        --result {output.result} \
        --plots_dir {output.plots} \
        """
        
#++++++++++++++++++++++++++++++++++++++++ 11 STATE DISCOVERY +++++++++++++++++++++++++++++++++++++
# 11.1 Refinement
rule purify_allgenes:
    input:
        # pickle_file = data_dir + "{major_datasets}/results/oncoBLADE/{expectation}/AutoGeneS/{BLADE_run}/{celltype_level}_STAD_ALL.oncoBLADE_output.k3.pickle", # cell_type
        pickle_file = data_dir + "{major_datasets}/results/oncoBLADE/{expectation}/AutoGeneS/{BLADE_run}/{celltype_level}_STAD_ALL.hierarchical_oncoBLADE_output.k6.pickle", # cell_type2
        bulk_file = data_dir + "{major_datasets}/results/oncoBLADE/STAD/ALL/bulk.pickle",
        signature = data_dir + "{major_datasets}/results/oncoBLADE/NewSignature_{celltype_level}.pickle",
        DEGs = data_dir +  "{major_datasets}/results/AutoGeneS/AutoGeneS_standard_genes_{celltype_level}_run6.txt"
    output:
        # obj = data_dir + "{major_datasets}/results/oncoBLADE/{expectation}/AutoGeneS/{BLADE_run}/{celltype_level}_STAD.purification.k3.pickle" # cell_type
        obj = data_dir + "{major_datasets}/results/oncoBLADE/{expectation}/AutoGeneS/{BLADE_run}/{celltype_level}_STAD.hier_purification.k6.pickle" # cell_type2
    params:
        ncores = 20  # Adjust as needed
    conda:
        "envs/purification.yaml"
    shell:
        """
        python3 scripts/purify_genes.py \
        --pickle_file {input.pickle_file} \
        --bulk {input.bulk_file} \
        --signature {input.signature} \
        --autogenes {input.DEGs} \
        --ncores {params.ncores} \
        --output_obj {output.obj} \
        """

# 11.2 State Discovery
rule run_StateDiscovery:
    input:
        # cf_file = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/{celltype_level}_STAD_ALL.estimated_cf.k3.csv", # cell type
        # purified_file = data_dir + "{major_datasets}/results/oncoBLADE/{expectation}/AutoGeneS/{BLADE_run}/{celltype_level}_STAD.purification.k3.pickle", # cell type
        cf_file = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/{celltype_level}_STAD_ALL.hierarchical_estimated_cf.k6.csv", # cell type 2
        purified_file = data_dir + "{major_datasets}/results/oncoBLADE/{expectation}/AutoGeneS/{BLADE_run}/{celltype_level}_STAD.hier_purification.k6.pickle", # cell type 2
    output:
        loadings = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/state/{celltype_level}_StateLoadings.csv",
        scores = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/state/{celltype_level}_StateScores.csv",
        plot = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/state/{celltype_level}_cophcor.png"
    conda:
        "envs/cell_state.yaml"
    shell:
        """
        python3 scripts/StateDiscovery.py \
        --cf_file {input.cf_file} \
        --purified_file {input.purified_file} \
        --loadings {output.loadings} \
        --scores {output.scores} \
        --plot {output.plot} \
        """

# 11.3 cluster heatmap that shows the top 1 gene for each cell state of each cell type
rule heatmap:
    input:
        loadings = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/state/{celltype_level}_StateLoadings.csv"
    output:
        heatmap = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/state/{celltype_level}_heatmap.pdf"
    params:
        cell_level = lambda wildcards: wildcards.celltype_level
    conda:
        "envs/heatmap.yaml"
    shell:
        """
        Rscript scripts/Heatmap_Loadings.R \
        --input {input.loadings} \
	--output {output.heatmap} \
        --cell_level {params.cell_level} \
        """    

# 11.4 bar plot that shows the top genes for each cell state of each cell type
rule barplot:
    input:
        loadings = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/state/{celltype_level}_StateLoadings.csv"
    output:
        directory = directory(data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/state/{celltype_level}_barplots/")
    params:
        cell_level = lambda wildcards: wildcards.celltype_level
    conda:
        "envs/barplot.yaml"
    shell:
        """
        Rscript scripts/state_barplot.R \
        --input {input.loadings} \
	--output {output.directory} \
        --cell_level {params.cell_level} \
        """    

rule assign_state:
    input:
        loadings = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/state/{celltype_level}_StateLoadings.csv",
        scores = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/state/{celltype_level}_StateScores.csv"
    output:
        output1 = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/state/{celltype_level}_genes_assigned_states.csv",
        output2 = data_dir + "{major_datasets}/results/evaluate/{BLADE_run}/{expectation}/state/{celltype_level}_samples_assigned_states.csv"
    conda:
        "envs/cell_state.yaml"
    shell:
        """
        python3 scripts/assign_state.py \
        --loadings {input.loadings} \
        --scores {input.scores} \
        --output1 {output.output1} \
        --output2 {output.output2} \
        """
