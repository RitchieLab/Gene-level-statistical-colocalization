This code can conduct statistical colocalization between GWAS-significant SNPs and the corresponding eQTLs (for a given tissue) on a gene-by-gene basis. This is ideally suited for situations in which we have a set of genes and tissues of interest (e.g. obtained from running Transcriptome wide Association Studies) and corresponding gene-expression and GWAS summary statistics.

There is no required version of R, but R (>=3.5.0) is preferred.

The following libraries are required to run this sofware:
* [coloc](https://cran.r-project.org/web/packages/coloc/index.html)
* [simsalapar](https://cran.r-project.org/web/packages/simsalapar/index.html)
* [optparse](https://www.rdocumentation.org/packages/optparse/versions/1.6.6)
* [data.table](https://cran.r-project.org/web/packages/data.table/data.table.pdf)
* Installed [GCTA](https://cnsgenomics.com/software/gcta/#Overview) v1.26 or higher

### Harmonization between GWAS and eQTL datasets
* Seee https://github.com/hakyimlab/MetaXcan/wiki/Tutorial:-GTEx-v8-MASH-models-integration-with-a-Coronary-Artery-Disease-GWAS for basic steps on harmonizing between GWAS and GTEx v8 datasets so that they are on the same genome build (run_harmonization code in the example also shows the basic protocol for harmonization between a GLGC lipid dataset and GTEx v8 eQTL dataset)

### Primary signals protocol:
* Note that the following protocol is applied to each combination of gene, trait (quantitative or case/control for a given GWAS dataset), and tissue. We run statistical colocalization between GWAS summary statistics and gene expression summary statistics obtained from GTEx v8.
* For running *colocalization*, we first identify all the SNPs in the GWAS dataset that are within a specified window (default = 1Mb) from the TSS and TES of the selected gene for a given tissue. Of these, we consider as “lead SNPs” all SNPs in the GWAS dataset with a p-value < a chosen threshold (default = 0.0001) that are at least a specified distance (default = 200 KB) apart from each other. 
* For each lead SNP, we collect all the SNPs within a specified distance (default = 200KB on either side of the SNP) in both GWAS and eQTL datasets. 
* We subsequently estimate the coloc probability of H3 (alternative hypothesis that eQTL and GWAS associations correspond to independent signals) and H4 (alternative hypothesis that eQTL and GWAS associations correspond to the same signal) for the lead SNP-eGene pair. 
* We assume a prior probability that a SNP is associated with (1) lipid phenotype (default=1E-04), (2) gene expression (default=1E=0-4), and (3) both GWAS and gene expression (default=1E-06) for all coloc analyses. 

### Secondary signals protocol: 
* OPTIONAL: We use the following GCTA-COJO/coloc protocol in the GWAS and eQTL datasets for each lead SNP to identify putative independent secondary associations at the locus. 
* For each lead SNP, we obtain p-values conditional on the top-eQTL (p-value < chosen threshold; default = 0.0001) in the eQTL dataset using --cojo-cond option to perform stepwise regression at that locus. We repeat this in the GWAS dataset as well if the top-eQTL also has GWAS p-value < chosen threshold (default = 0.0001).
* In our example dataset, we use 1000 genome EUR as reference dataset to calculate pairwise LD. In reality, it is recommended to use at least 5K individuals in the reference dataset.
* We run *coloc* between GWAS and eQTL datasets (same parameters as before) at the locus using conditional probabilities obtained from COJO.

### Collect and compile
* We select the lead SNP with P[H3]<0.5 that has the highest P[H4] for a given eGene and use this as the corresponding P[H3] and P[H4] for the gene. 
* We can subsequently filter out all genes whose P[H3]>0.5 for a given tissue (these could be LD contaminated; i.e. GWAS causal variant and eQTL are different but in LD). 

## Input files

```
#sample GWAS file (pre-harmonization)
  SNP          BP	  CHR	A1	A2	BETA	   SE	    P   FREQ        N
rs1172982   100230111	  1	 t	 c	0.0043	0.0055	0.4689	0.3219    89888
rs1172981   100230197	  1	 t	 c	0.0057	0.0103	0.7688	0.06069   89888
rs11166327  100230867	  1	 c	 t	0.0053	0.0106	0.7725	0.06069	  89888
rs4908018   100234367	  1	 a	 c	0.0058	0.0103	0.7887	0.05937	  89841
rs2392072   100234743	  1	 a	 g	0.0039	0.0058	0.4052	0.3153	  86877
```

```
#sample eQTL file
    gene_id             variant_id	    tss_distance	ma_samples	ma_count      maf	pval_nominal	  slope	        slope_se
ENSG00000261456.5	chr10_11501_C_A_b38       -62662	 124	        124	      0.116105	0.0400165	-0.187313	0.0909802
ENSG00000261456.5	chr10_11553_G_C_b38       -62610	 108	        108	      0.105058	0.0216029	-0.221475	0.0961113
ENSG00000261456.5	chr10_18924_A_C_b38       -55239         86	        86	      0.0776173	0.283069	-0.114165	0.106242
ENSG00000261456.5	chr10_44215_T_TTCTG_b38	  -29948	 13	        13	      0.0122411	0.608588	0.130235	0.254164
ENSG00000261456.5	chr10_45349_G_A_b38	  -28814	 21	        21	      0.0180723	0.579648	0.124365	0.224383
```

```
#sample Input dataset

SYPL2   TC      Liver                           GLGC    ENSG00000143028 1
ANGPTL3 LDL     Liver                           GLGC    ENSG00000132855 1
PSRC1   TC      Adipose_Visceral_Omentum        GLGC    ENSG00000134222 1
PSRC1   LDL     Liver                           GLGC    ENSG00000134222 1
PLTP    TG      Adipose_Visceral_Omentum        GLGC    ENSG00000100979 20
NBEAL1  TC      Adipose_Visceral_Omentum        GLGC    ENSG00000144426 2
```

```
## Pre-saved files needed to run this code (these are downloaded with the code):
#eQTL sample size file (GTEX v8)

Tissue          X..RNASeq.and.Genotyped.samples
Adipose_Subcutaneous                    581
Adipose_Visceral_Omentum                469
Adrenal_Gland                           233
Artery_Aorta                            387
Artery_Coronary                         213
Artery_Tibial                           584
Brain_Amygdala                          129
Brain_Anterior_cingulate_cortex_BA24    147

#Gene start stop file (GRCh38)

ENSG00000265692 82290046        82292814        LINC01970       17
ENSG00000260563 82293716        82294910        AC132872.1      17
ENSG00000173762 82314868        82317608        CD7             17
ENSG00000141574 82321024        82334074        SECTM1          17
ENSG00000182459 82359247        82363775        TEX19           17
ENSG00000278964 82362349        82363196        AC132938.4      17
ENSG00000181408 82371400        82375586        UTS2R           17
ENSG00000260011 82381110        82382690        AC132938.1      17
ENSG00000181396 82389210        82418637        OGFOD3          17
ENSG00000264812 82400703        82401382        AC132938.2      17
ENSG00000201239 82417226        82417321        Y_RNA           17
ENSG00000169660 82418318        82442645        HEXD            17
ENSG00000279066 82425498        82427310        HEXD-IT1        17

#Chromosome-wise reference files for COJO in {chr}.bim/.bed/.fam format (1000 Genomes EUR (b38) available in the example dataset)

```

## Output file

```
Trait                     TC                             TC                               TC                            TC
Tissue          Adipose_Visceral_Omentum        Adipose_Visceral_Omentum        Adipose_Visceral_Omentum        Adipose_Visceral_Omentum
GWAS.Data                GLGC                           GLGC                             GLGC                          GLGC
CHR                        2                              2                                2                            2
GWAS.Lead.SNP         rs7581601                     rs11694172                        rs3769676                     rs1174604
GWAS.BP               202633063                       202667581                       203401620                     203556561
GWAS.P               4.494987e-09                   1.417578e-11                    7.515972e-08                  3.270255e-08
GWAS.Beta               -0.0217                        0.0277                           0.0199                       0.0210
GWAS.SE                 0.0037                         0.0041                           0.0037                       0.0038
GWAS.N                  187300                         187092                           185112                       187263
GWAS.A1                   C                              G                                C                             C
GWAS.A2                   A                              A                                T                             T
GWAS.Freq             0.3338                          0.2164                            0.3443                       0.3496
GWAS.MAF              0.3338                          0.2164                            0.3443                       0.3496
ENSG_gene           ENSG00000144426                 ENSG00000144426                  ENSG00000144426              ENSG00000144426
Gene                  NBEAL1                           NBEAL1                           NBEAL1                        NBEAL1
Region            202437904:202832555             202469849:202867547               203202545:203600396         203357100:203747335
Gene.Start            203014879                       203014879                        203014879                    203014879
Gene.End              203226378                       203226378                        203226378                    203226378
Gene.Num.SNPs           1584                           1584                              1584                         1584
Gene.Num.GWAS.Sig.SNPs    7                             7                                  7                            7
Best.Causal.SNP     2:202797474                      2:202797474                      2:203206233                  2:203556561
TopSNP_eQTL         2:202791831                      2:202791831                      2:203206233                  2:203425314
eQTL.Beta           -0.13291200                      0.00157705                       0.01767270                   0.01257750
eQTL.SE              0.0327190                       0.0342433                        0.0306130                    0.0302590
eQTL.P               5.85398e-05                     9.63290e-01                      5.64066e-01                 6.77882e-01
eQTL.N                  469                             469                              469                          469
eQTL.A1                  A                               A                                T                            T
eQTL.A2                  C                               G                                C                            C
eQTL.MAF             0.260684                         0.250533                         0.357602                     0.352878
Coloc.Ratio          0.0012040                        0.0007724                        0.0179100                    0.0522000
Coloc.H0             1.476e-07                        9.218e-08                        4.173e-04                    4.350e-02
Coloc.H1             0.004493                         0.002836                         0.008311                     0.865400
Coloc.H2             3.267e-05                        3.239e-05                        4.659e-02                    4.143e-03
Coloc.H3             0.99430                          0.99640                          0.92780                      0.08239
Coloc.H4             0.0011980                        0.0007702                        0.0169200                    0.0045380
Coloc2.Ratio         0.0012910                        0.0007955                        0.2025000                    0.1898000
Coloc2.H0            4.699e-08                        2.820e-08                        4.549e-02                    4.616e-02
Coloc2.H1            0.0014300                        0.0008674                        0.9058000                    0.9185000
Coloc2.H2            3.277e-05                        3.245e-05                        1.879e-03                    1.385e-03
Coloc2.H3            0.99720                          0.99830                          0.03733                      0.02749
Coloc2.H4            0.0012890                        0.0007948                        0.0094800                    0.0064410

```
* Note that the output file has been transposed for easy viewing of all columns (here rows).

## Setup and Example

1. Clone the repository.
``` 
git clone https://github.com/RitchieLab/Gene-level-statistical-colocalization
```
2. Go to the cloned folder (set as working directory).
``` 
cd Gene-level-statistical-colocalization 
```
3. Download and unzip example data https://ritchielab.org/files/Lipid_Pleiotropy_project/gene_level_coloc_example_data.tar.gz
and save downloaded data under the cloned folder.

4. ```tar -xvzf gene_level_coloc_example_data.tar.gz```

5. run harmonization (see ```run_harmonization``` example code under ```harmonized_gwas``` folder)

6. Add ```gcta64``` to the folder path.

7. Run ```sh example```.

```
Rscript run_gene_level_coloc.R \
--gwas_data_name="GLGC" \
--trait="LDL" \
--tissue="Liver" \
--gene_of_interest="ENSG00000144426" \
--run_cojo=FALSE \
--cojo_maf=0.01 \
--chr=2 \
--coloc_p1=1e-04 \
--coloc_p2=1e-04 \
--coloc_p12=1e-06 \
--lead_snp_window_size=200000 \
--gene_boundary_window_size=1000000 \
--gwas_p_threshold=0.0001 \
--eqtl_p_threshold=0.0001 \
--gwas_response_type="quant" \
--gwas_file=GLGC_LDL_GWAS_harmonized.txt.gz \
--genes_file=${genes_file} \
--output_folder=${output_folder} \
--reference_folder=${reference_folder} \
--eqtl_folder=${path_gtex} \
--core=10 \
--eqtl_sample_size_file=${eqtl_sample_size_file}

```
* Note that it is not recommended to run a for loop over genes specified in Input file in practice; one can parallelize runs in an HPC environment.

Following is an explanation of the listed parameters:

  * --*chromosome* The chromosome to which the given gene(s) correspond(s) 
  * --*trait* Name of the trait 
  * --*tissue* Name of tissue corresponding to the eQTL dataset 
  * --*gwas_data_name* Name of GWAS dataset (e.g. GLGC, GIANT)
  * --*gene_of_interest* ENSG_gene(s) of interest from genes_file (ignore decimal point)
  * --*gene_boundary_window_size* The window around the chosen gene from the TSS and TES of the gene (default = 1 Mb) 
  * --*lead_snp_window_size* The window around the lead SNP used for colocalization analysis
  * --*gwas_p_threshold* The threshold for GWAS p-value for SNP variants corresponding to the chosen gene(s) (default = 1E-03) 
  * --*eqtl_p_threshold* The threshold for eQTL p-value for SNP variants mapping to the chosen eGene(s) (default = 1E-02)
  * --*gwas_response_type* The response variable type for GWAS (quant or cc)
  * --*coloc_p1* Prior probability a SNP is associated with GWAS trait (default = 1E-04)
  * --*coloc_p2* Prior probability a SNP is associated with gene expression (default = 0.001) 
  * --*coloc_p12* Prior probability a SNP is associated with GWAS trait and gene expression (default = 1E-06)
  * --*gwas_file* File for tab separated *harmomized* GWAS summary statistics data (with header) for chosen trait 
  * --*genes_file* File (with path) for tab separated list of chosen genes with column names =  ENSG_gene, gene_start_position, gene_stop_position, chromosome
  * --*snps_file* External Input file of GWAS varID to be used for analyses, i.e. chr:bp
  * --*run_cojo* Flag to run GCTA-COJO for identifying secondary signals; TRUE/FALSE (default=FALSE)
  * --*cojo_maf* MAF threshold used in GCTA-COJO(default=0.01) 
  * --*cojo_p* P-value threshold for gcta-cojo (default = 1E-03)
  * --*reference_folder* Path of the folder with plink files for LD calculation in gcta (should have chromosome number in the filename in .chromosome.bim/bed/fam format)
  * --*eqtl_folder* Path of the GTEX v8 datasets split by chromosome
  * --*output_folder* Path of the output folder
  * --*core* Number of cores to run parallel tasks (default = 10) 
  * --*ld_folder* Path of the folder with plink files for LD calculation in gcta (should have chromosome number in the filename in .chromosome.bim/bed/fam format)
  * --*eqtl_sample_size_file* Filename (with path) of sample sizes for eQTL datasets across different tissues; has two columns corresponding to tissue name and sample size
  
6. Run ```remove_non_colocalized_genes.R```
```
Rscript remove_non_colocalized_genes.R \
  --file_name chr1_GLGC_LDL_Adipose_Subcutaneous_colocProbs.txt \
  --coloc_p_h3_threshold 0.5 \
  --coloc_p_h4_threshold 0.01 \
  --output_folder "output"
 ```
 
Following is an explanation of the listed parameters:

  * --*file_name* Load coloc output file obtained after running run_gcta_and_coloc.R
  * --*coloc_p_h3_threshold* Threshold for minP[H3] per gene across all lead SNPs in the gene (default=0.5)
  * --*coloc_p_h4_threshold* Threshold for maxP[H4] per gene across all lead SNPs in the gene (default=0.01)
  * --*output_folder* Folder where list of genes after filtering out ones with no evidence of colocalization with GWAS SNPs is saved
  
## Reference
The manuscript "Unified framework identifies novel replicating links between plasma lipids and diseases from Electronic Health Records across large-scale cohorts" is currently under review.
