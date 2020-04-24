There is no required version of R, but R (>=3.5.0) is preferred.

To run simulation using run_twas_simulation.R, the following libraries are required:
* [coloc](https://cran.r-project.org/web/packages/coloc/index.html)
* [simsalapar](https://cran.r-project.org/web/packages/simsalapar/index.html)
* [optparse](https://www.rdocumentation.org/packages/optparse/versions/1.6.6)
* Installed [GCTA](https://cnsgenomics.com/software/gcta/#Overview)

For running *colocalization*, we first identify a list of genes in a chromosome. 
Next, we identify all the SNPs in the GWAS dataset that are within a specified window from the TSS and TES of the selected genes. 
Of these, we consider as “lead SNPs” all SNPs in the GWAS dataset with a p-value < a chosen threshold that are > 10 KB apart from each other. 
### GCTA-COJO protocol: 
* We use the following GCTA-COJO protocol in the GWAS and eQTL datasets for each lead SNP. 
#### GWAS dataset 
* We first identify SNPs around a chosen window of each lead SNP from the GWAS dataset, which we will refer to as our “GWAS-testable SNPs”. 
* For each lead SNP, we perform conditional analyses using GCTA-COJO on the corresponding GWAS-testable SNPs and obtain “conditionally independent” SNPs associated with the lipid phenotype. 
#### eQTL dataset
* Here, we identify the subset of genes (referred to as “eGenes”) whose boundaries are within the chosen window of the lead SNP. 
* For each selected eGene, we identify all the associated SNPs in the eQTL dataset, or the set of “eQTL-testable SNPs”. 
* Again, we perform conditional analyses using GCTA-COJO to obtain “conditionally independent” SNPs associated with gene expression for a given lead SNP-eGene pair in the selected tissue. 
#### GCTA steps
* We use the --cojo-slct option to perform model selection and get a list of independently associated testable SNPs (p-value < chosen threshold). 
* If there are at least two independent testable SNPs at a distance > 100KB (but < 1 MB) from the lead SNP, we run GCTA-COJO again to perform association analysis on all testable SNPs conditioned on the independently associated SNPs procured from the first GCTA run. 
* In our example dataset, we use 1000 genome EUR (chromosome 1) as reference dataset to calculate pairwise LD. 
### Coloc protocol:
* We collect the conditionally independent SNPs overlapping between GWAS and eQTL datasets for each lead SNP-eGene pair. 
* We subsequently estimate the coloc probability of H3 (alternative hypothesis that eQTL and GWAS associations correspond to independent signals) and H4 (alternative hypothesis that eQTL and GWAS associations correspond to the same signal) for the lead SNP-eGene pair. 
* In case a lead SNP corresponds to more than one eGene for a given tissue, we perform colocalization on both lead SNP-eGene pairs. 
* We assume a prior probability that a SNP is associated with (1) lipid phenotype (default=1E-04), (2) gene expression (default=0.001), and (3) both GWAS and gene expression (default=1E-06) for all coloc analyses. 
* We select the lead SNP with the lowest p[H3] and and highest p[H4] for a given eGene and use this as the corresponding p[H3] for the gene. 
* We can subsequently filter out all genes whose p[H3]>0.5 and p[H4]<0.5 for a given tissue. 


## Setup and Example

1. Clone the repository.
``` 
git clone https://github.com/RitchieLab/Remove-LD-contaminated-genes
```
2. Go to the cloned folder (set as working directory)
``` 
cd Remove-LD-contaminated-genes 
```
3. Download and unzip example data https://ritchielab.org/files/Lipid_Pleiotropy_project/LD_Contamination_example_data.tar.gz
and save downloaded data under the cloned folder
4. Run 'run_gcta_and_coloc.R'.
```
Rscript run_gcta_and_coloc.R \
  --chromosome 1 \
  --window_size 10000 \
  --gwas_p_threshold 1E-06 \
  --eqtl_file "LD_Contamination_example_data/Adipose_Subcutaneous.allpairs.chr1.txt" \
  --gwas_file "LD_Contamination_example_data/jointGwasMc_LDL_chr1_formatted.txt" \
  --genes_file "LD_Contamination_example_data/Genes_list.txt" \
  --maf 0.01 \
  --trait "LDL" \
  --tissue "Adipose_Subcutaneous" \
  --gwas_data_name "GLGC" \
  --cojo_p 1E-06 \
  --gene_of_interest "ENSG00000134243" \
  --output_folder "output" \
  --core 10 \
  --ld_folder "LD_Contamination_example_data" \
  --eqtl_sample_size "LD_Contamination_example_data/GTEx_v7_tissue_specific_sample_sizes.txt" \
  --coloc_p1 1E-04 \
  --coloc_p2 0.001 \
  --coloc_p12 1E-06 
```

Here is an explanation of the listed parameters

  * --*chromosome* The chromosome to which the given gene(s) correspond(s) 
  * --*window_size* The window chosen around the chosen gene from the TSS and TES of the gene (default = 1 Mb) 
  * --*gwas_p_threshold* The threshold for GWAS p-value for SNP variants corresponding to the chosen gene(s) (default = 1E-03) 
  * --*eqtl_file* File for eQTL summary statistics (from GTEx) for chosen tissue, split by chromosome, with path 
  * --*gwas_file* File for tab separated GWAS summary statistics data (with header) for chosen trait with column names SNP, BP, CHR, A1, A2, BETA, SE, P, FREQ, with path 
  * --*genes_file* File (with path) for tab separated list of chosen genes with column names =  ENSG_gene, gene_start_position, gene_stop_position, chromosome
  * --*maf* MAF threshold (default=0.01) 
  * --*trait* Name of the trait 
  * --*tissue* Name of tissue corresponding to the eQTL dataset 
  * --*gwas_data_name* Name of GWAS dataset (e.g. GLGC, GIANT)
  * --*cojo_p* P-value threshold for gcta-cojo (default = 1E-03) 
  * --*gene_of_interest* ENSG_gene(s) of interest from genes_file (ignore decimal point)
  * --*output_folder* Path of the output folder
  * --*core* Number of cores to run parallel tasks (default = 10) 
  * --*ld_folder* Path of the folder with plink files for LD calculation in gcta (should have chromosome number in the filename in .chromosome.bim/bed/fam format)
  * --*eqtl_sample_size* Filename (with path) of sample sizes for eQTL datasets across different tissues; has two columns corresponding to tissue name and sample size
  * --*coloc_p1* Prior probability a SNP is associated with GWAS trait (default = 1E-04)
  * --*coloc_p2* Prior probability a SNP is associated with gene expression (default = 0.001) 
  * --*coloc_p12* Prior probability a SNP is associated with GWAS trait and gene expression (default = 1E-06)
