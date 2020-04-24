There is no required version of R, but R (>=3.5.0) is preferred.

To run simulation using run_twas_simulation.R, the following libraries are required:
* [coloc](https://cran.r-project.org/web/packages/coloc/index.html)
* [simsalapar](https://cran.r-project.org/web/packages/simsalapar/index.html)
* Installed [GCTA](https://cnsgenomics.com/software/gcta/#Overview)

* For running *colocalization*, we first identify a list of TWAS-significant genes in a chromosome. 
* Next, we identify all the SNPs in the GWAS dataset that are within a 1 MB region from the TSS and TES of the TWAS-significant genes. 
* Of these, we consider as “lead SNPs” all SNPs in the GWAS dataset with a p-value < 0.001 that are > 10 KB apart from each other. 
### GCTA-COJO protocol: 
* We use the following GCTA-COJO protocol in the GWAS and eQTL datasets for each lead SNP. 
#### GWAS dataset 
* We first identify SNPs around 1 MB region of each lead SNP from the GWAS dataset, which we will refer to as our “GWAS-testable SNPs”. 
* For each lead SNP, we perform conditional analyses using GCTA-COJO on the corresponding GWAS-testable SNPs and obtain “conditionally independent” SNPs associated with the lipid phenotype. 
#### eQTL dataset
* Here, we identify the subset of TWAS-significant genes (referred to as “eGenes”) whose boundaries are within 1 MB of the lead SNP. 
* For each selected eGene, we identify all the associated SNPs in the eQTL dataset, or the set of “eQTL-testable SNPs”. 
* Again, we perform conditional analyses using GCTA-COJO to obtain “conditionally independent” SNPs associated with gene expression for a given lead SNP-eGene pair. 
#### GCTA steps
* We use the --cojo-slct option to perform model selection and get a list of independently associated testable SNPs (p-value < 0.001). 
* If there are at least two independent testable SNPs at a distance > 100KB (but < 1 MB) from the lead SNP, we run GCTA-COJO again to perform association analysis on all testable SNPs conditioned on the independently associated SNPs procured from the first GCTA run. 
* In our example dataset, we use 1000 genome EUR (chromosome 1) as reference dataset to calculate pairwise LD. 
### Coloc protocol:
* We collect the conditionally independent SNPs overlapping between GWAS and eQTL datasets for each lead SNP-eGene pair. 
* We subsequently estimate the coloc probability of H3 (alternative hypothesis that eQTL and GWAS associations correspond to independent signals) and H4 (alternative hypothesis that eQTL and GWAS associations correspond to the same signal) for the lead SNP-eGene pair. 
* In case a lead SNP corresponded to more than one eGene for a given tissue, we perform colocalization on both lead SNP-eGene pairs. 
* We assume a prior probability that a SNP is associated with (1) lipid phenotype (default=1E-04), (2) gene expression (default=0.001), and (3) both GWAS and gene expression (default=1E-06) for all coloc analyses. 
* We select the lead SNP with the lowest p[H3] and lowest GWAS p-value for a given eGene and use this as the corresponding p[H3] for the gene. 
* We subsequently filter out all genes whose p[H3]>0.5 and that had no SNPs with p[H4]>0.528 for a given combination of GWAS cohort, lipid phenotype and tissue. 
