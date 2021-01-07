#!/usr/bin/env Rscript
rm(list=ls())

################## This code uses statistical colocalization to identify whether lead SNPs from GWAS summary statistics colocalize with cis-eQTL for a give gene. This gene could be TWAS-significant and have "LD contamination", i.e.expression predictor SNPs and GWAS causal SNPs are different but in LD. Here, we test whether GWAS significant variants in the gene have coloc P[H3]<0.5 and P[H4]>0.5  with corresponding eQTL SNPs (for a given tissue) implying the GWAS-significant variant and the eQTL are likely to correspond to the same signal. The software used for colocalization analysis is "coloc", which assumes a single causal variant in a locus. In case we find evidence for colocalization, we also run GCTA-COJO, which can condition variants mapping to the gene on the top eQTL-associated GWAS SNP. These conditional p-values are then fed into "coloc" to test for potential secondary signals.

#### Yogasudha Veturi 26October2020

#############################################
## 1. Load necessary libraries
#############################################

suppressMessages(library(coloc))
suppressMessages(library(simsalapar))
suppressMessages(library(optparse))

#############################################
## 2. Define parameters and prepare datasets
#############################################

opt_list <- list(
  make_option("--trait",type="character",default="Trait",help="Trait name (default = %default)"),
  make_option("--tissue",type="character",help="Tissue name"),
  make_option("--chromosome",type="integer",help="The chromosome to which the given gene corresponds (1 to 23)"),
  make_option("--lead_snp_window_size",type="integer",default=100000,help="The window chosen around a lead SNP for colocalization analysis (default = %default)"),
  make_option("--gwas_p_threshold",type="numeric",default=1e-03,help="The threshold for GWAS p-value for SNP variants corresponding to the chosen gene (default = %default)"),
  make_option("--gwas_response_type",type="character",help="The response variable type for GWAS (quant or cc)"),
  make_option("--gwas_data_name",type="character",default="GWAS",help="Name of GWAS dataset (e.g. GLGC, GIANT) (default = %default)"),
  make_option("--gwas_file",type="character",help="File for tab separated GWAS summary statistics data (with header) for chosen trait with column names derived after harmomization with GTEX v8: variant_id,panel_variant_id,chromosome,position,effect_allele,non_effect_allele,frequency,zscore,effect_size,standard_error,pvalue,sample_size,n_cases (for case-control), with path"),
  make_option("--genes_file",type="character",help="File (with path) for tab separated list of chosen genes with column names =  Gene, ENSG_gene, gene_start_position, gene_stop_position, chromosome (GRCh 37)"),
  make_option("--gene_of_interest",type="character",help="ENSG_gene of interest from genes_file (ignore decimal point)"),
  make_option("--gene_boundary_window_size",type="integer",default=1000000,help="The window chosen around the chosen gene from the TSS and TES of the gene (default = %default)"),
  make_option("--eqtl_p_threshold",type="numeric",default=1e-03,help="The threshold for eqtl nominal p-value for SNP variants corresponding to the chosen gene (default = %default)"),
  make_option("--cojo_maf",type="numeric",default=0.01,help="MAF threshold used in GCTA-COJO (default = %default)"),
  make_option("--eqtl_folder",type="character",help="Path of the folder with GTEX v8 eQTL allpairs data in the following format: Tissue.allpairs.chr'chr'.txt.gz"),
  make_option("--output_folder",type="character",help="Path of the output folder"),
  make_option("--reference_folder",type="character",help="Path of the folder with plink files for LD calculation in gcta (should have chromosome number in the filename in .chromosome.bim/bed/fam format)"),
  make_option("--eqtl_sample_size_file",type="character",help="Filename (with path) of sample sizes for eQTL datasets across different tissues; has two columns corresponding to tissue name and sample size"),
  make_option("--coloc_p1",type="numeric",default=1e-04,help="Prior probability a SNP is associated with GWAS trait (default = %default)"),
  make_option("--coloc_p2",type="numeric",default=1e-04,help="Prior probability a SNP is associated with gene expression (default = %default)"),
  make_option("--coloc_p12",type="numeric",default=1e-06,help="Prior probability a SNP is associated with GWAS trait and gene expression (default = %default)"),
  make_option("--core",type="integer",default=10,help="Number of cores to run parallel tasks (default = %default)"),
  make_option("--snps_file",type="character",default=NA,help="External Input file of GWAS varID, i.e. chr:bp (default = %default)")
)

opts <- parse_args(OptionParser(option_list=opt_list))

#2.1 parse values
chr <- opts$chromosome
if(is.numeric(chr)==FALSE) stop("Chromosome as integer (1 to 23) not specified")
lead.snp.window.size <- opts$lead_snp_window_size
gene.boundary.window.size <- opts$gene_boundary_window_size
gwas.p.threshold <- opts$gwas_p_threshold
eqtl.p.threshold <- opts$eqtl_p_threshold
gwas.file <- opts$gwas_file
gwas.response.type <- opts$gwas_response_type
if(gwas.response.type%in%c("cc","quant")==FALSE) stop("Response type should either be 'cc' or 'quant'")
genes.file <- opts$genes_file
cojo.maf <- opts$cojo_maf
gwas.data.name <- opts$gwas_data_name
trait <- opts$trait
output.folder <- opts$output_folder
reference.folder <- opts$reference_folder
tissue <- opts$tissue
gene.of.interest <- opts$gene_of_interest
eqtl.sample.size.file <- opts$eqtl_sample_size_file
if(length(eqtl.sample.size.file)==0) stop("eQTL sample size file for all tissues in GTEX v8 not specified")
coloc.p1 <- opts$coloc_p1
coloc.p2 <- opts$coloc_p2
coloc.p12 <- opts$coloc_p12
core <- opts$core
snps.file <- opts$snps_file
eqtl.folder <- opts$eqtl_folder
eqtl.file <- paste0(eqtl.folder,"/",tissue,".allpairs.chr",chr,".txt.gz")

#2.2 Read and prepare data files

gwas <- as.data.frame(data.table::fread(gwas.file,header=T,stringsAsFactors=FALSE,sep="\t"))
if (length(grep(paste0(c("variant_id","panel_variant_id","chromosome","position","effect_allele","non_effect_allele","frequency","zscore","effect_size","standard_error","pvalue","sample_size","n_cases"),collapse="|"),colnames(gwas),ignore.case=TRUE))<13) stop("Required columns are not present in the requested format in the GWAS dataset")
colnames(gwas) = c("SNP","panel_variant_id","CHR","BP","A1","A2","A1FREQ","ZSCORE","BETA","SE","P","N","NCASES")
gwas[,"CHR"] = substr(gwas[,"CHR"],4,nchar(gwas[,"CHR"]))
gwas <- gwas[which(gwas[,"CHR"]==chr),]
if(nrow(gwas)==0) stop("No SNPs in corresponding chromosome; check GWAS file and chromosome specification")
gwas <- gwas[order(gwas[,"BP"],decreasing=FALSE),]
gwas[,"P"] <- as.numeric(gwas[,"P"])
gwas[,"varID"] <- paste0(gwas[,"CHR"],":",gwas[,"BP"])
gwas.fin <- gwas[which(gwas[,"P"]<gwas.p.threshold),]
if(nrow(gwas.fin)==0) stop("No sig GWAS hits for given trait  in this chromosome")

genes.list <- as.data.frame(data.table::fread(genes.file,header=T,stringsAsFactors=FALSE,sep="\t"))
eqtl <- as.data.frame(data.table::fread(eqtl.file,header=T,stringsAsFactors=FALSE))
reference.files <- list.files(reference.folder)
reference.files <- reference.files[grep(paste0("[.]",chr,"[.]"),reference.files)]
reference.filename <- unique(substr(reference.files, 1, nchar(reference.files)-4))
reference.file.bim = reference.files[grep(".bim",reference.files)]
if(length(reference.file.bim)==0) stop("Reference filename is not in .chromosome.bim/bed/fam format")
if(length(reference.file.bim)>1) stop("More than one reference filename")
reference.snps <- as.data.frame(data.table::fread(paste0(reference.folder,"/",reference.file.bim),header=F,sep="\t"))
genes.fin <- genes.list[which(genes.list[,"Gene_stable_ID"]==gene.of.interest),]
if(nrow(genes.fin)==0) stop("Genes(s) of interest not in file")
N.eqtl <- read.table(eqtl.sample.size.file,header=T,stringsAsFactors=FALSE)
N.eqtl.tissue <- N.eqtl[which(N.eqtl[,1]==tissue),2]
if(length(which(N.eqtl[,1]==tissue))!=1) stop("The eqtl sample size file does not contain the tissue of interest in column 1")
eqtl.gene.id.split <- as.character(unlist(strsplit(eqtl[,"gene_id"],"[.]")))
eqtl$eqtl_gene_id <- eqtl.gene.id.split[grep("ENSG",eqtl.gene.id.split)]
gene_fin <- genes.fin[,"Gene_stable_ID"]
if(length(gene_fin)==0) stop("No genes specified")
gene_name <- genes.fin[,"Gene_name"]
message(tissue," / ", trait," / ",gwas.data.name," / ",gene_name," / ", gene_fin, " / chr=",chr)

#############################################
## 3. Isolate relevant SNPs for analysis
#############################################

dir.create(paste0(output.folder,"/intermediate"), recursive=T)
setwd(output.folder)

#3.1 Get SNPs for GWAS dataset (p-value < chosen threshold) that lie within (default=1MB) of chosen gene boundaries

if(is.na(snps.file)) snps_gwas <- gwas.fin[,"varID"] else {
  panel_variant_id_gwas <- read.table(snps.file,stringsAsFactors=FALSE)$V1
  panel_variant_id_gwas <- panel_variant_id_gwas[which(panel_variant_id_gwas%in%gwas.fin[,"panel_variant_id"])]
  snps_gwas <- gwas.fin[match(panel_variant_id_gwas,gwas.fin[,"panel_variant_id"]),"varID"]
}
pos_gwas <- as.numeric(gwas.fin[match(snps_gwas,gwas.fin[,"varID"]),"BP"])

snps_gwas2 <- NA
tmp <- unique(snps_gwas[which(pos_gwas>=(genes.fin[1,"Gene_start_(bp)"]-gene.boundary.window.size)&pos_gwas<=(genes.fin[1,"Gene_end_(bp)"]+gene.boundary.window.size))])
snps_gwas2 <- c(snps_gwas2,tmp)
snps_gwas2 <-snps_gwas2[-1]
lead_snps_gwas <- sort(unique(snps_gwas2))

#3.2 Get SNPs for eQTL dataset (p-value < chosen threshold)

eqtl.fin <- eqtl[which(eqtl$eqtl_gene_id==gene_fin),]
eqtl_chrpos_alleles <-  as.character(unlist(strsplit(as.character(unlist(strsplit(as.character(unlist(eqtl.fin[,"variant_id"])),"_b38"))),"_")))
eqtl_chr <- eqtl_chrpos_alleles[grep("chr",eqtl_chrpos_alleles)]
eqtl.fin[,"CHR"] <- substr(eqtl_chr,4,nchar(eqtl_chr))
if(length(unique(eqtl.fin[,"CHR"]))==0) stop("eQTL allpairs dataset does not correspond to the chromosome specified")
eqtl_pos_alleles <- eqtl_chrpos_alleles[-grep("chr",eqtl_chrpos_alleles)]
eqtl.fin[,"BP"] <- eqtl_pos_alleles[-grep("A|C|T|G|a|c|t|g",eqtl_pos_alleles)]
eqtl_alleles <- eqtl_pos_alleles[grep("A|C|T|G|a|c|t|g",eqtl_pos_alleles)]
eqtl.fin[,"A1"] <- eqtl_alleles[-which(c(1:length(eqtl_alleles))%%2==0)]
eqtl.fin[,"A2"] <- eqtl_alleles[which(c(1:length(eqtl_alleles))%%2==0)]
eqtl.fin[,"varID"] <- paste0(eqtl.fin[,"CHR"],":",eqtl.fin[,"BP"])
eqtl.fin[,"CHR"] = as.numeric(eqtl.fin[,"CHR"])
eqtl.fin[,"BP"] = as.numeric(eqtl.fin[,"BP"])
lead_snps_eqtl <- eqtl.fin[which(eqtl.fin[,"pval_nominal"]<eqtl.p.threshold),"varID"]

#3.3 Get final list of lead SNPs
common_lead_snps <- intersect(lead_snps_gwas,lead_snps_eqtl)
#if(length(lead_snps)==0) stop("No lead SNPs at the chosen p-value thresholds")
if(length(common_lead_snps)==0) message("No overlapping 'lead' SNPs between GWAS and eQTL (GTEXv8) datasets for given gene")

###3.3.1 Function to get lead SNPs at a distance of (default=100KB) from each other
lead_snps <- lead_snps_gwas
gwas.fin2 <- gwas.fin[which(gwas.fin[,"varID"]%in%lead_snps),]
getleadsnp.fn <- function(data,baseRow,jump=1) {
  success <- TRUE
  x <- data[baseRow,]
  bottomRow <- baseRow

  while (success) {
    # do 
    Distance <- lead.snp.window.size
    startDistance <- x["BP"]
    snp <- paste0(x["CHR"],":",x["BP"])
    bottomRow <- bottomRow+jump
    tmp1 <- abs(data[bottomRow,"BP"]-startDistance)
    success <- (tmp1 < Distance) &  bottomRow<nrow(data)
    # message("DistanceDifference=",tmp1," / bottomRow=", bottomRow)
  }
  return(bottomRow)
}

if(length(lead_snps)>1){
  rows <- i <- 1
    while(i<=length(lead_snps)-1) {
      i <- getleadsnp.fn(gwas.fin2,baseRow=i)
      rows <- c(rows,i)
    }
    snp <- NA
    for(i in 1:(length(rows)-1)){
      if(rows[i+1]-rows[i]>1) tmp <- gwas.fin2[which(gwas.fin2[,"P"]==min(gwas.fin2[(rows[i]:(rows[i+1]-1)),"P"])),"varID"] else {
        if(i<(length(rows)-1)) tmp <- gwas.fin2[rows[i],"varID"]; if(i==(length(rows)-1)) tmp <- c(gwas.fin2[rows[i],"varID"],gwas.fin2[rows[i+1],"varID"])}
      snp=c(snp,tmp)

  lead_snps <- unique(snp[-1])
 }
}

###3.3.2 LD clump lead SNPs

#  gwas.for.clump <- gwas.fin[which(gwas.fin[,"varID"]%in%lead_snps),]
#  lead_snps <- invisible(clump_data(gwas.for.clump, clump_r2 = clump.r2))[,"varID"]

#3.4 Count number of SNPs (significant vs total) that lie within (chosen window size) of chosen gene's boundaries

data2 <- gwas[which(gwas[,"P"]<5e-08),]
genes.fin$num.snps=genes.fin$num.sig.snps=NA

for(i in 1:nrow(genes.fin)){
  genes.fin[i,"num.sig.snps"] <- nrow(data2[which(data2[,"BP"]>=(genes.fin[i,"Gene_start_(bp)"]-gene.boundary.window.size)&data2[,"BP"]<=(genes.fin[i,"Gene_end_(bp)"]+gene.boundary.window.size)),])
  genes.fin[i,"num.snps"] <- nrow(gwas[which(gwas[,"BP"]>=(genes.fin[i,"Gene_start_(bp)"]-gene.boundary.window.size)&gwas[,"BP"]<=(genes.fin[i,"Gene_end_(bp)"]+gene.boundary.window.size)),])
  # print(i)
}

########################################################################################################################
## 4. Loop over all lead SNPs and prepare datasets for coloc (and GCTA-COJO to identify secondary signals, if necessary)
########################################################################################################################

#4.1 Function to find top and bottom rows in the GWAS data file that lie  within the boundaries of the chosen gene (+/- the window size surrounding the gene boundaries)

getWindow.fn <- function(data,baseRow,size=window.size,jump) {
  success <- FALSE
  while (!success) {
    # do 
    x <- data[baseRow,]
    extraDistance <- size    ### Create window around the gene
    startDistance <- x["BP"]
    snp <- paste0(x["CHR"],":",x["BP"])
    success2.1 = success2.2 = success2.3 = success2.4 = TRUE
    topRow <- bottomRow <- baseRow
    while (success2.1) {
      bottomRow <- bottomRow+jump
      tmp1 <- abs(data[bottomRow,"BP"]-startDistance)
      success2.1 <- (tmp1 < extraDistance) &  bottomRow<(nrow(data)-jump)
      #message("endDistanceDifference=",tmp1," / bottomRow=", bottomRow)
    }
    while(!success2.1&success2.2) {
      bottomRow <- bottomRow-1
      tmp1 <- abs(data[bottomRow,"BP"]-startDistance)
      success2.2 <- (tmp1 > extraDistance)&  bottomRow<nrow(data)
      # message("endDistanceDifference=",tmp1," / bottomRow=", bottomRow)
    }
    while (success2.3) {
      topRow <- topRow-jump
      tmp2 <- abs(startDistance-data[topRow,"BP"])
      success2.3 <- (tmp2 < extraDistance)&  topRow>1+jump
      #message("startDistanceDifference=",tmp2," / topRow=", topRow)
    }
    while(!success2.3&success2.4) {
      topRow <- topRow+1
      tmp2 <- abs(startDistance-data[topRow,"BP"])
      success2.4 <- (tmp2 > extraDistance)&  topRow>1
      # message("startDistanceDifference=",tmp2," / topRow=", topRow)
    }

    success <- (!success2.2)&(!success2.4)
  }
  out = list(topRow,bottomRow)
  return(out)
}

#4.2 Start for loop
if(gwas.response.type=="quant") coloc_probs_header = c("Trait","Tissue","GWAS.Data","CHR","GWAS.Lead.SNP","GWAS.BP","GWAS.P","GWAS.Beta","GWAS.SE","GWAS.N","GWAS.A1","GWAS.A2","GWAS.Freq","GWAS.MAF","ENSG_gene","Gene","Region","Gene.Start","Gene.End","Gene.Num.SNPs","Gene.Num.GWAS.Sig.SNPs","Best.Causal.SNP","TopSNP_eQTL","eQTL.Beta","eQTL.SE","eQTL.P","eQTL.N","eQTL.A1","eQTL.A2","eQTL.MAF","Coloc.Ratio","Coloc.H0","Coloc.H1","Coloc.H2","Coloc.H3","Coloc.H4","Coloc2.Ratio","Coloc2.H0","Coloc2.H1","Coloc2.H2","Coloc2.H3","Coloc2.H4")
if(gwas.response.type=="cc") coloc_probs_header = c("Trait","Tissue","GWAS.Data","CHR","GWAS.Lead.SNP","GWAS.BP","GWAS.P","GWAS.Beta","GWAS.SE","GWAS.N","GWAS.NCASES","GWAS.A1","GWAS.A2","GWAS.Freq","GWAS.MAF","ENSG_gene","Gene","Region","Gene.Start","Gene.End","Gene.Num.SNPs","Gene.Num.GWAS.Sig.SNPs","Best.Causal.SNP","TopSNP_eQTL","eQTL","eQTL.Beta","eQTL.SE","eQTL.P","eQTL.N","eQTL.A1","eQTL.A2","eQTL.MAF","Coloc.Ratio","Coloc.H0","Coloc.H1","Coloc.H2","Coloc.H3","Coloc.H4","Coloc2.Ratio","Coloc2.H0","Coloc2.H1","Coloc2.H2","Coloc2.H3","Coloc2.H4")

write.table(t(coloc_probs_header), file=paste0("chr",chr,"_",gwas.data.name,"_",trait,"_",tissue,"_",gene_name,"_colocProbs.txt"),row.names=F,col.names=F,quote=F,sep="\t")
for(j in 1:length(lead_snps)){ # Loop starts for lead snps
  varid <- lead_snps[j]
  pos <- as.numeric(as.character(unlist(strsplit(varid,":"))))[2]
  coloc_probs <- matrix(,1,length(coloc_probs_header))
  colnames(coloc_probs) <- coloc_probs_header
  coloc_probs <- data.frame(coloc_probs)
  for(i in 1:ncol(coloc_probs)) coloc_probs[,i] <- as.character(unlist(coloc_probs[,i]))

  base <- which(paste0(gwas[,"CHR"],":",gwas[,"BP"])==varid)[1]
  boundaryRows <- data.frame(getWindow.fn(baseRow=base,data=gwas,size=lead.snp.window.size,jump=100))  ## Get boundary rows for GWAS dataset; these will yield all the GWAS SNPs within chosen window size on either side of lead SNP
  x_gwas <- gwas[boundaryRows[[1]]:boundaryRows[[2]],]
  if(length(which(x_gwas[,"A1FREQ"]%in%c(0,1,NA)))>0) x_gwas <- x_gwas[-which(x_gwas[,"A1FREQ"]%in%c(0,1,NA)),]
  snp <- x_gwas[which(paste0(x_gwas[,"CHR"],":",x_gwas[,"BP"])==varid),"SNP"]
  x_gwas$varID <- paste0(chr,":",x_gwas[,"BP"])
  eqtl.fin2 <- eqtl.fin[which(eqtl.fin[,"BP"]>gwas[boundaryRows[[1]],"BP"]&eqtl.fin[,"BP"]<=gwas[boundaryRows[[2]],"BP"]),]

  x_eqtl <- cbind(eqtl.fin[which(eqtl.fin[,"varID"]%in%x_gwas[,"varID"]),],"N"=N.eqtl.tissue)
  if(length(which(x_eqtl[,"maf"]==0))>0) x_eqtl <- x_eqtl[-which(x_eqtl[,"maf"]==0),]
  if(gwas.response.type=="quant")  gwas.coloc <- x_gwas[,c("CHR","varID","BP","A1","A1FREQ","BETA","SE","P","N")];
  if(gwas.response.type=="cc")  {
    x_gwas[,"PROPCASES"] <- x_gwas[,"NCASES"]/x_gwas[,"N"]
    gwas.coloc <- x_gwas[,c("CHR","varID","BP","A1","A1FREQ","BETA","SE","P","N","PROPCASES")];
  }
  eqtl.coloc <- x_eqtl[,c("CHR","varID","BP","A1","maf","slope","slope_se","pval_nominal","N")]
  colnames(eqtl.coloc) = c("Chr","SNP","bp","refA","freq","b","se","p","n")
  if(gwas.response.type=="quant")  colnames(gwas.coloc) = colnames(eqtl.coloc)
  if(gwas.response.type=="cc") colnames(gwas.coloc) = c("Chr","SNP","bp","refA","freq","b","se","p","n","s")

  eqtl.coloc <- eqtl.coloc[which(eqtl.coloc[,"SNP"]%in%gwas.coloc[,"SNP"]),]
  gwas.coloc <- gwas.coloc[which(gwas.coloc[,"SNP"]%in%eqtl.coloc[,"SNP"]),]
  gwas.coloc <- gwas.coloc[match(eqtl.coloc[,"SNP"],gwas.coloc[,"SNP"]),]
  if(length(which(is.na(gwas.coloc[,"SNP"])))>0) gwas.coloc = gwas.coloc[-which(is.na(gwas.coloc[,"SNP"])),]
  eqtl.coloc <- eqtl.coloc[match(gwas.coloc[,"SNP"],eqtl.coloc[,"SNP"]),]
  if(length(which(is.na(eqtl.coloc[,"SNP"])))>0) eqtl.coloc = eqtl.coloc[-which(is.na(eqtl.coloc[,"SNP"])),]
  stopifnot(all.equal(as.character(unlist(eqtl.coloc[,"SNP"])),as.character(unlist(gwas.coloc[,"SNP"]))))

  eqtl.coloc$maf <- ifelse(eqtl.coloc[,"freq"]<0.5,eqtl.coloc[,"freq"],1-eqtl.coloc[,"freq"])
  gwas.coloc$maf <- ifelse(gwas.coloc[,"freq"]<0.5,gwas.coloc[,"freq"],1-gwas.coloc[,"freq"])

  ##########################################################################################
  ## 5. Run colocalization analysis between GWAS and eQTL datasets (identify primary signal)
  ##########################################################################################

  #5.1 Prepare eQTL and GWAS datasets for colocalization analysis 
  if(nrow(eqtl.coloc)>0){
    if(gwas.response.type=="quant") {
      coloc_data <- data.frame("gwas"=gwas.coloc[,c("SNP","maf","p","n","b","se")],"eqtl"=eqtl.coloc[,c("SNP","maf","p","n","b","se")])
      colnames(coloc_data) = c("gwas.SNP","gwas.maf","gwas.p","gwas.n","gwas.b","gwas.se","eqtl.SNP","eqtl.maf","eqtl.p","eqtl.n","eqtl.b","eqtl.se")
      if(length(which(is.na(coloc_data[,"eqtl.p"])|is.na(coloc_data[,"gwas.p"])))>0) coloc_data <- coloc_data[-which(is.na(coloc_data[,"eqtl.p"])|is.na(coloc_data[,"gwas.p"])),]
      coloc <- tryCatch.W.E(coloc.abf(dataset1=list(snp=coloc_data[,"gwas.SNP"],pvalues=coloc_data[,"gwas.p"], N=round(median(coloc_data[,"gwas.n"]),digits=0), MAF=coloc_data[,"gwas.maf"],type=gwas.response.type), dataset2=list(snp=coloc_data[,"eqtl.SNP"],pvalues=coloc_data[,"eqtl.p"], N=median(coloc_data[,"eqtl.n"]), type="quant",MAF=coloc_data[,"eqtl.maf"]),p1=coloc.p1,p2=coloc.p2,p12=coloc.p12))
    }
    if(gwas.response.type=="cc") {
      coloc_data <- data.frame("gwas"=gwas.coloc[,c("SNP","maf","p","n","b","se","s")],"eqtl"=eqtl.coloc[,c("SNP","maf","p","n","b","se")])
      colnames(coloc_data) = c("gwas.SNP","gwas.maf","gwas.p","gwas.n","gwas.b","gwas.se","gwas.s","eqtl.SNP","eqtl.maf","eqtl.p","eqtl.n","eqtl.b","eqtl.se")
      if(length(which(is.na(coloc_data[,"eqtl.p"])|is.na(coloc_data[,"gwas.p"])))>0) coloc_data <- coloc_data[-which(is.na(coloc_data[,"eqtl.p"])|is.na(coloc_data[,"gwas.p"])),]
      coloc <- tryCatch.W.E(coloc.abf(dataset1=list(snp=coloc_data[,"gwas.SNP"],pvalues=coloc_data[,"gwas.p"], N=round(median(coloc_data[,"gwas.n"]),digits=0), MAF=coloc_data[,"gwas.maf"],s=coloc_data[,"gwas.s"],type=gwas.response.type), dataset2=list(snp=coloc_data[,"eqtl.SNP"],pvalues=coloc_data[,"eqtl.p"], N=median(coloc_data[,"eqtl.n"]), type="quant",MAF=coloc_data[,"eqtl.maf"]),p1=coloc.p1,p2=coloc.p2,p12=coloc.p12))
    }
    warn <- inherits(coloc$warning,"warning")
    error <- inherits(coloc$value,"error")
  }
  if(coloc$value[[1]][[5]]<0.5&coloc$value[[1]][[6]]>0.75) message("\nStrong support (primary signal) for colocalization in this region!")
  best.causal.snp <- as.character(unlist(coloc$value[[2]][which(coloc$value[[2]][,"SNP.PP.H4"]==max(coloc$value[[2]][,"SNP.PP.H4"])),"snp"]))
  if(coloc$value[[1]][[5]]<0.5&coloc$value[[1]][[6]]>0.5&coloc$value[[1]][[6]]<=0.75) message("\nGood support (primary signal) for colocalization in this region")
  if(coloc$value[[1]][[6]]<0.2) message("\nNo support (primary signal) for colocalization in this region")
  if(coloc$value[[1]][[5]]<0.5&coloc$value[[1]][[6]]>0.2&coloc$value[[1]][[6]]<=0.5) message("\nModerate support (primary signal) for colocalization in this region.")
  if(coloc$value[[1]][[5]]>0.5) message("\nEvidence of LD contamination found in gene.")

  ###################################################################################################
  ## 6. Condition SNPs on the top eQTL-associated SNP using GCTA-COJO and check for secondary signals
  ###################################################################################################

  #6.1 Get topsnp for eQTL data
  eqtl_topsnp <- eqtl.fin2[which(eqtl.fin2[,"pval_nominal"]==min(eqtl.fin2[,"pval_nominal"])),"varID"][1]
  if(eqtl.fin2[which(eqtl.fin2[,"varID"]==eqtl_topsnp),"pval_nominal"]<=eqtl.p.threshold) {use.conditional.probs <- TRUE;message(paste0("eQTL top SNP has a p-value<=",eqtl.p.threshold,"; proceeding with search for secondary signals\n"))} else {message(paste0("eQTL top SNP has a p-value>",eqtl.p.threshold,"; halting search for secondary signals\n"));use.conditional.probs <- FALSE}

  if(use.conditional.probs){

    gwas.matched.ref <- match(x_gwas[,"varID"],paste0(reference.snps[,1],":",reference.snps[,4]))
    eqtl.matched.ref <- match(x_eqtl[,"varID"],paste0(reference.snps[,1],":",reference.snps[,4]))

    x_gwas.bp <- reference.snps[gwas.matched.ref,4]
    x_gwas["varID"] <- paste0(chr,":",x_gwas.bp)
    x_gwas["SNP"] <- reference.snps[gwas.matched.ref,2]
    x_gwas_gcta <- x_gwas[,c("SNP","A1","A2","A1FREQ","BETA","SE","P","N")]
    if(length(which(is.na(x_gwas_gcta[,"SNP"])))>0) x_gwas_gcta <- x_gwas_gcta[-which(is.na(x_gwas_gcta[,"SNP"])),]
    x_eqtl <- cbind("SNP"=reference.snps[eqtl.matched.ref,2],x_eqtl[,-1])
    x_eqtl_gcta <- cbind(x_eqtl[match(x_gwas_gcta[,"SNP"],x_eqtl[,"SNP"]),c("SNP","A1","A2","maf","slope","slope_se","pval_nominal")],"N"=N.eqtl.tissue)
    if(nrow(x_eqtl_gcta)==0) {message("tissue=",tissue," / j=",j," / ",trait," / ",gwas.data.name,"  / num.lead.snps=",length(lead_snps)) ; next}
    eqtl_topsnp2 <- reference.snps[which(paste0(reference.snps[,1],":",reference.snps[,4])==eqtl_topsnp),2]
    write.table(eqtl_topsnp2,file=paste0("intermediate/",trait,"_",gwas.data.name,"_chr",chr,"_pos_",pos,"_",gene_name,"_",tissue,"_topsnp.txt"),row.names=F,col.names=T,quote=F,sep="\t")
    if(length(eqtl_topsnp2)==0) {
      print("eQTL top SNP does not have a match in the reference dataset: Are reference and eQTL data on the same build?"); next}

    eqtl_topsnp_full = eqtl.fin2[which(eqtl.fin2[,"varID"]==eqtl_topsnp),]
    x_eqtl_gcta = unique(rbind(x_eqtl_gcta,cbind("SNP"=eqtl_topsnp2, eqtl_topsnp_full[,c("A1","A2","maf","slope","slope_se","pval_nominal")],"N"=N.eqtl.tissue)))
    write.table(x_eqtl_gcta,file=paste0("intermediate/eQTL_",tissue,"_chr",chr,"_pos_",pos,"_",gene_name,".txt"),row.names=F,col.names=T,quote=F,sep="\t")

    gwas_topsnp_full = gwas.fin[which(gwas.fin[,"varID"]==eqtl_topsnp),]
    if(nrow(gwas_topsnp_full)>0)
    {
      x_gwas_gcta = unique(rbind(x_gwas_gcta,cbind("SNP"=eqtl_topsnp2, gwas_topsnp_full[,c("A1","A2","A1FREQ","BETA","SE","P","N")])))
      if(x_gwas_gcta[which(x_gwas_gcta[,"SNP"]==eqtl_topsnp2),"P"]<gwas.p.threshold){
        write.table(x_gwas_gcta,file=paste0("intermediate/GWAS_",trait,"_",gwas.data.name,"_chr",chr,"_pos_",pos,"_",gene_name,".txt"),row.names=F,col.names=T,quote=F,sep="\t")

        #6.1 Run GCTA-COJO to perform association analysis of all the included SNPs conditional on the given top SNP in the GWAS dataset. Results are saved in the .cma.cojo file

        invisible(system(paste0("module load gcta ;gcta64 --bfile ",reference.folder,"/",reference.filename," --chr ",chr," --maf ",cojo.maf," --cojo-file intermediate/GWAS_",trait,"_",gwas.data.name,"_chr",chr,"_pos_",pos,"_",gene_name,".txt --cojo-cond intermediate/",trait,"_",gwas.data.name,"_chr",chr,"_pos_",pos,"_",gene_name,"_",tissue,"_topsnp.txt --out intermediate/GWAS_",trait,"_",gwas.data.name,"_chr_",chr,"_pos_",pos,"_",gene_name,"_conditionalP"),intern=TRUE))
        tmp.1 <- tryCatch.W.E(read.table(paste0("intermediate/GWAS_",trait,"_",gwas.data.name,"_chr_",chr,"_pos_",pos,"_",gene_name,"_conditionalP.cma.cojo"),header=T))
        warn.1 <- inherits(tmp.1$warning,"warning")
        error.1 <- inherits(tmp.1$value,"error")
        if(error.1==FALSE){
          gwas.coloc <- tmp.1$value
          gwas.coloc[,"SNP"] <- paste0(gwas.coloc[,"Chr"],":",gwas.coloc[,"bp"])
          gwas.coloc <- gwas.coloc[,c("Chr","SNP","bp","refA","freq","bC","bC_se","pC","n")]
          colnames(gwas.coloc) <- c("Chr","SNP","bp","refA","freq","b","se","p","n")
        }
      }
    }
    #6.2 Run GCTA-COJO to perform association analysis of all the included SNPs conditional on the given top SNP in the eQTL dataset. Results are saved in the .cma.cojo file
    invisible(system(paste0("module load gcta ;gcta64 --bfile ",reference.folder,"/",reference.filename," --chr ",chr," --maf ",cojo.maf," --cojo-file intermediate/eQTL_",tissue,"_chr",chr,"_pos_",pos,"_",gene_name,".txt --cojo-cond intermediate/",trait,"_",gwas.data.name,"_chr",chr,"_pos_",pos,"_",gene_name,"_",tissue,"_topsnp.txt --out intermediate/",tissue,"_chr_",chr,"_pos_",pos,"_",gene_name,"_conditionalP"),intern=TRUE))
    tmp.2 <-  tryCatch.W.E(read.table(paste0("intermediate/",tissue,"_chr_",chr,"_pos_",pos,"_",gene_name,"_conditionalP.cma.cojo"),header=T))
    warn.2 <- inherits(tmp.2$warning,"warning")
    error.2 <- inherits(tmp.2$value,"error")
    if(error.2==FALSE){
      eqtl.coloc <- tmp.2$value
      eqtl.coloc[,"SNP"] <- paste0(eqtl.coloc[,"Chr"],":",eqtl.coloc[,"bp"])
      eqtl.coloc <- eqtl.coloc[,c("Chr","SNP","bp","refA","freq","bC","bC_se","pC","n")]
      colnames(eqtl.coloc) <- c("Chr","SNP","bp","refA","freq","b","se","p","n")
   }

    if(exists("error.1")) {if(error.1==FALSE&error.2==FALSE) proceed.for.coloc="TRUE"}
    if(exists("error.1")==FALSE&error.2==FALSE) proceed.for.coloc="TRUE"
    if(proceed.for.coloc) {
      eqtl.coloc <- eqtl.coloc[which(eqtl.coloc[,"SNP"]%in%gwas.coloc[,"SNP"]),]
      gwas.coloc <- gwas.coloc[which(gwas.coloc[,"SNP"]%in%eqtl.coloc[,"SNP"]),]
      gwas.coloc <- gwas.coloc[match(eqtl.coloc[,"SNP"],gwas.coloc[,"SNP"]),]
      if(length(which(is.na(gwas.coloc[,"SNP"])))>0) gwas.coloc = gwas.coloc[-which(is.na(gwas.coloc[,"SNP"])),]
      eqtl.coloc <- eqtl.coloc[match(gwas.coloc[,"SNP"],eqtl.coloc[,"SNP"]),]
      if(length(which(is.na(eqtl.coloc[,"SNP"])))>0) eqtl.coloc = eqtl.coloc[-which(is.na(eqtl.coloc[,"SNP"])),]
      stopifnot(all.equal(as.character(unlist(eqtl.coloc[,"SNP"])),as.character(unlist(gwas.coloc[,"SNP"]))))

      eqtl.coloc$maf <- ifelse(eqtl.coloc[,"freq"]<0.5,eqtl.coloc[,"freq"],1-eqtl.coloc[,"freq"])
      gwas.coloc$maf <- ifelse(gwas.coloc[,"freq"]<0.5,gwas.coloc[,"freq"],1-gwas.coloc[,"freq"])
      if(gwas.response.type=="cc") {
        ncases <- gwas[match(gwas.coloc[,"SNP"],gwas[,"SNP"]),"NCASES"]
        gwas.coloc$s <- ncases/gwas.coloc[,"n"]
      }

    ##########################################################################################################
    ## 7. Run colocalization analysis between GWAS SNPs and eQTL using conditional probabilites from GCTA-COJO
    ##########################################################################################################
    #7.1 Prepare eQTL and GWAS datasets for colocalization analysis 
      if(gwas.response.type=="quant"){
        coloc_data2 <- data.frame("gwas"=gwas.coloc[,c("SNP","maf","p","n","b","se")],"eqtl"=eqtl.coloc[,c("SNP","maf","p","n","b","se")])
        colnames(coloc_data2) = c("gwas.SNP","gwas.maf","gwas.p","gwas.n","gwas.b","gwas.se","eqtl.SNP","eqtl.maf","eqtl.p","eqtl.n","eqtl.b","eqtl.se")
        if(length(which(is.na(coloc_data2[,"eqtl.p"])|is.na(coloc_data2[,"gwas.p"])))>0) coloc_data2 <- coloc_data2[-which(is.na(coloc_data2[,"eqtl.p"])|is.na(coloc_data2[,"gwas.p"])),]
        coloc2 <- tryCatch.W.E(coloc.abf(dataset1=list(pvalues=coloc_data2[,"gwas.p"], N=round(median(coloc_data2[,"gwas.n"]),digits=0), MAF=coloc_data2[,"gwas.maf"],type=gwas.response.type), dataset2=list(pvalues=coloc_data2[,"eqtl.p"], N=median(coloc_data2[,"eqtl.n"]), type="quant",MAF=coloc_data2[,"gwas.maf"]),p1=coloc.p1,p2=coloc.p2,p12=coloc.p12)$summary)
      }
      if(gwas.response.type=="cc") {
        coloc_data2 <- data.frame("gwas"=gwas.coloc[,c("SNP","maf","p","n","b","se","s")],"eqtl"=eqtl.coloc[,c("SNP","maf","p","n","b","se")])
        colnames(coloc_data2) = c("gwas.SNP","gwas.maf","gwas.p","gwas.n","gwas.b","gwas.se","gwas.s","eqtl.SNP","eqtl.maf","eqtl.p","eqtl.n","eqtl.b","eqtl.se")
        if(length(which(is.na(coloc_data2[,"eqtl.p"])|is.na(coloc_data2[,"gwas.p"])))>0) coloc_data2 <- coloc_data2[-which(is.na(coloc_data2[,"eqtl.p"])|is.na(coloc_data2[,"gwas.p"])),]
        coloc2 <- tryCatch.W.E(coloc.abf(dataset1=list(pvalues=coloc_data2[,"gwas.p"], N=round(median(coloc_data2[,"gwas.n"]),digits=0), MAF=coloc_data2[,"gwas.maf"],s=coloc_data2[,"gwas.s"],type=gwas.response.type), dataset2=list(pvalues=coloc_data2[,"eqtl.p"], N=median(coloc_data2[,"eqtl.n"]), type="quant",MAF=coloc_data2[,"gwas.maf"]),p1=coloc.p1,p2=coloc.p2,p12=coloc.p12)$summary)
      }
      warn2 <- inherits(coloc2$warning,"warning")
      error2 <- inherits(coloc2$value,"error")
      if(coloc2$value[[6]]<0.2) message("\nNo support for secondary signals at this locus\n")
      if(coloc2$value[[5]]<0.5&coloc2$value[[6]]>0.5) {message("\nStrong support for secondary signals at this locus!\n")}
     if(coloc2$value[[5]]<0.5&coloc2$value[[6]]>0.2&coloc2$value[[6]]<=0.5) {message("\nModerate support for secondary signals at this locus!\n")}
      if(coloc2$value[[5]]>0.5) message("\nNo support for secondary signals at this locus\n")
    }
  }
  ############################
  ## 8. Prepare output dataset
  ############################

  coloc_probs[,"GWAS.Lead.SNP"] <- snp
  coloc_probs[,"Trait"] <- trait
  coloc_probs[,"Tissue"] <- tissue
  coloc_probs[,"CHR"] <- chr
  coloc_probs[,"GWAS.BP"] <- pos
  coloc_probs[,"GWAS.P"] <- x_gwas[which(x_gwas[,"varID"]==varid),"P"][1]
  coloc_probs[,"GWAS.Beta"] <- x_gwas[which(x_gwas[,"varID"]==varid),"BETA"][1]
  coloc_probs[,"GWAS.SE"] <- x_gwas[which(x_gwas[,"varID"]==varid),"SE"][1]
  coloc_probs[,"GWAS.Data"] <- gwas.data.name
  coloc_probs[,"GWAS.A1"] <- x_gwas[which(x_gwas[,"varID"]==varid),"A1"][1]
  coloc_probs[,"GWAS.A2"] <- x_gwas[which(x_gwas[,"varID"]==varid),"A2"][1]
  coloc_probs[,"GWAS.Freq"] <- x_gwas[which(x_gwas[,"varID"]==varid),"A1FREQ"][1]
  coloc_probs[,"GWAS.MAF"] <- coloc_data[which(coloc_data[,"gwas.SNP"]==varid),"gwas.maf"]
  coloc_probs[,"ENSG_gene"] <- gene_fin
  coloc_probs[,"Gene"] <- gene_name
  coloc_probs[,"TopSNP_eQTL"] <- eqtl_topsnp
  coloc_probs[,"Best.Causal.SNP"] <- best.causal.snp
  coloc_probs[,"eQTL.Beta"] <- x_eqtl[which(x_eqtl[,"varID"]==varid),"slope"][1]
  coloc_probs[,"eQTL.SE"] <- x_eqtl[which(x_eqtl[,"varID"]==varid),"slope_se"][1]
  coloc_probs[,"eQTL.P"] <- x_eqtl[which(x_eqtl[,"varID"]==varid),"pval_nominal"][1]
  coloc_probs[,"eQTL.A1"] <- x_eqtl[which(x_eqtl[,"varID"]==varid),"A1"][1]
  coloc_probs[,"eQTL.A2"] <- x_eqtl[which(x_eqtl[,"varID"]==varid),"A2"][1]
  coloc_probs[,"eQTL.MAF"] <- coloc_data[which(coloc_data[,"eqtl.SNP"]==varid),"eqtl.maf"]
  coloc_probs[,"Region"] <- paste0(gwas[boundaryRows[[1]],"BP"],":",gwas[boundaryRows[[2]],"BP"])
  coloc_probs[,"Gene.Num.SNPs"] <- genes.fin[which(genes.fin[,"Gene_stable_ID"]==gene_fin),"num.snps"]
  coloc_probs[,"Gene.Num.GWAS.Sig.SNPs"] <- genes.fin[which(genes.fin[,"Gene_stable_ID"]==gene_fin),"num.sig.snps"]
  coloc_probs[,"GWAS.N"] <-  x_gwas[which(x_gwas[,"varID"]==varid),"N"][1]
  coloc_probs[,"eQTL.N"] <-  x_eqtl[which(x_eqtl[,"varID"]==varid),"N"][1]
  if(gwas.response.type=="cc")   coloc_probs[,"GWAS.NCASES"] <-  x_gwas[which(x_gwas[,"varID"]==varid),"NCASES"][1]
  coloc_probs[,"Gene.Start"] <- genes.fin[which(genes.fin[,"Gene_stable_ID"]==gene_fin),"Gene_start_(bp)"]
  coloc_probs[,"Gene.End"] <- genes.fin[which(genes.fin[,"Gene_stable_ID"]==gene_fin),"Gene_end_(bp)"]

  if(error==FALSE){
    coloc_probs[,"Coloc.Ratio"] = signif(coloc$value[[1]]["PP.H4.abf"]/(coloc$value[[1]]["PP.H3.abf"]+coloc$value[[1]]["PP.H4.abf"]),digits=4)
    coloc_probs[,c("Coloc.H0","Coloc.H1","Coloc.H2","Coloc.H3","Coloc.H4")] = signif(coloc$value[[1]][-1],digits=4)
  }

  if(use.conditional.probs)
    if(proceed.for.coloc){
      if(error2==FALSE){
        coloc_probs[,"Coloc2.Ratio"] = signif(coloc2$value["PP.H4.abf"]/(coloc2$value["PP.H3.abf"]+coloc2$value["PP.H4.abf"]),digits=4)
        coloc_probs[,c("Coloc2.H0","Coloc2.H1","Coloc2.H2","Coloc2.H3","Coloc2.H4")] = signif(coloc2$value[-1],digits=4)
        rm(error2, gwas.coloc,eqtl.coloc)
      }
  }
  message("tissue=",tissue," / lead.snp.num=",j," / lead.snp=",varid," / eqtl_topsnp=", eqtl_topsnp," / " ,trait," / ",gwas.data.name," / num.lead.snps=",length(lead_snps)," / gene=",gene_name, "\n")

  write.table(coloc_probs, file=paste0("chr",chr,"_",gwas.data.name,"_",trait,"_",tissue,"_",gene_name,"_colocProbs.txt"),row.names=F,col.names=F,quote=F,sep="\t",append=T)
  rm(error,coloc_probs,coloc_data,x_eqtl,x_gwas)
} ### end for loop


