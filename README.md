This code is to identify and remove genes for which the expression predictor SNPs and causal phenotype SNPs are different but in linkage disequilibrium. 

There is no required version of R, but R (>=3.5.0) is preferred.

The following libraries are required to run this sofware:
* [coloc](https://cran.r-project.org/web/packages/coloc/index.html)
* [simsalapar](https://cran.r-project.org/web/packages/simsalapar/index.html)
* [optparse](https://www.rdocumentation.org/packages/optparse/versions/1.6.6)
* [data.table](https://cran.r-project.org/web/packages/data.table/data.table.pdf)
* Installed [GCTA](https://cnsgenomics.com/software/gcta/#Overview) v1.26 or higher

For running *colocalization*, we first identify a list of genes in a chromosome. 
Next, we identify all the SNPs in the GWAS dataset that are within a specified window from the TSS and TES of the selected genes. 
Of these, we consider as “lead SNPs” all SNPs in the GWAS dataset with a p-value < a chosen threshold that are > 10 KB apart from each other. 
### GCTA-COJO protocol: 
* We use the following GCTA-COJO protocol in the GWAS and eQTL datasets for each lead SNP. 
#### GWAS dataset 
* We first identify SNPs around a chosen window of each lead SNP from the GWAS dataset, which we will refer to as our “GWAS-testable SNPs”. 
* For each lead SNP, we perform conditional analyses using GCTA-COJO on the corresponding GWAS-testable SNPs and obtain “conditionally independent” SNPs associated with the lipid phenotype. 

```
#sample GWAS file
SNP	BP	CHR	A1	A2	BETA	SE	P	FREQ	N
rs1172982	100230111	1	t	c	0.0043	0.0055	0.4689	0.3219	89888
rs1172981	100230197	1	t	c	0.0057	0.0103	0.7688	0.06069	89888
rs11166327	100230867	1	c	t	0.0053	0.0106	0.7725	0.06069	89888
rs4908018	100234367	1	a	c	0.0058	0.0103	0.7887	0.05937	89841
rs2392072	100234743	1	a	g	0.0039	0.0058	0.4052	0.3153	86877
```
#### eQTL dataset
* Here, we identify the subset of genes (referred to as “eGenes”) whose boundaries are within the chosen window of the lead SNP. 
* For each selected eGene, we identify all the associated SNPs in the eQTL dataset, or the set of “eQTL-testable SNPs”. 
* Again, we perform conditional analyses using GCTA-COJO to obtain “conditionally independent” SNPs associated with gene expression for a given lead SNP-eGene pair in the selected tissue. 

```
#sample eQTL file
            gene_id          variant_id tss_distance ma_samples ma_count
1 ENSG00000227232.4 1_13417_C_CGAGA_b37       -16136         50       50
2 ENSG00000227232.4     1_17559_G_C_b37       -11994          7        7
3 ENSG00000227232.4     1_54421_A_G_b37        24868         27       27
4 ENSG00000227232.4     1_54490_G_A_b37        24937        119      127
5 ENSG00000227232.4     1_61920_G_A_b37        32367         16       17
         maf pval_nominal      slope  slope_se
1 0.07225430   0.00908288  0.3556660 0.1354910
2 0.00923483   0.15611600  0.4947870 0.3480440
3 0.03678470   0.63643300 -0.0881399 0.1862840
4 0.17887300   0.43983700  0.0705390 0.0912018
5 0.02348070   0.52871600 -0.1293250 0.2050630
```

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
2. Go to the cloned folder (set as working directory).
``` 
cd Remove-LD-contaminated-genes 
```
3. Download and unzip example data https://ritchielab.org/files/Lipid_Pleiotropy_project/LD_Contamination_example_data.tar.gz
and save downloaded data under the cloned folder.
4. Add ```gcta64``` to the folder path.

5. Run ```run_gcta_and_coloc.R```.
```
Rscript run_gcta_and_coloc.R \
  --chromosome 1 \
  --window_size 1000000 \
  --gwas_p_threshold 1e-05 \
  --eqtl_file "LD_Contamination_example_data/Adipose_Subcutaneous.allpairs.chr1.txt" \
  --gwas_file "LD_Contamination_example_data/jointGwasMc_LDL_chr1_formatted.txt" \
  --genes_file "LD_Contamination_example_data/Genes_list.txt" \
  --maf 0.01 \
  --trait "LDL" \
  --tissue "Adipose_Subcutaneous" \
  --gwas_data_name "GLGC" \
  --cojo_p 0.001 \
  --gene_of_interest "ENSG00000134243" \
  --output_folder "output" \
  --core 10 \
  --ld_folder "LD_Contamination_example_data" \
  --eqtl_sample_size "LD_Contamination_example_data/GTEx_v7_tissue_specific_sample_sizes.txt" \
  --coloc_p1 1e-04 \
  --coloc_p2 0.001 \
  --coloc_p12 1E-06 
```

Following is an explanation of the listed parameters:

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
  
6. Run ```remove_ld_contaminated_genes.R```
```
Rscript remove_ld_contaminated_genes.R \
  --file_name chr1_GLGC_LDL_Adipose_Subcutaneous_colocProbs.txt \
  --coloc_p_h3_threshold 0.5 \
  --coloc_p_h4_threshold 0.5 \
  --output_folder "output"
 ```
 
Following is an explanation of the listed parameters:

  * --*file_name* Load coloc output file obtained after running run_gcta_and_coloc.R
  * --*coloc_p_h3_threshold* Threshold for minP[H3] per gene across all lead SNPs in the gene (default=0.5)
  * --*coloc_p_h4_threshold* Threshold for maxP[H4] per gene across all lead SNPs in the gene (default=0.5)
  * --*output_folder* Folder where list of genes after removing LD-contaminated genes is saved
  
## Reference
The manuscript "Unified framework identifies novel replicating links between plasma lipids and diseases from Electronic Health Records across large-scale cohorts" is currently under review.
