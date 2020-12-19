rm(list=ls())

################## This code uses statistical colocalization to identify genes whose eQTLs colocalize with a chosen set of GWAS-significant SNPs. 
#Here, we find (i)conditionally independent variants in each gene in GWAS and eQTL datasets using "GCTA-COJO" and 
#(ii) estimate the colocalization probabilities P[H3] and P[H4] per gene between GWAS and eQTL datasets using conditional p-values from GCTA-COJO in "coloc".
#### Yogasudha Veturi 15Dec2020

#############################################
## 1. Load necessary libraries
#############################################

suppressMessages(library(coloc))
suppressMessages(library(simsalapar))
suppressMessages(library(optparse))
## NOTE: gcta64 must be added to the path 

#############################################
## 2. Define parameters and prepare datasets
#############################################

opt_list <- list(
  make_option("--trait",type="character",help="Trait name"),
  make_option("--tissue",type="character",help="Tissue name"),
  make_option("--chromosome",type="integer",help="The chromosome to which the given gene corresponds"),
  make_option("--lead_snp_window_size",type="integer",default=100000,help="The window chosen around a lead SNP for colocalization analysis (default = %default)"),
  make_option("--gwas_p_threshold",type="numeric",default=1e-03,help="The threshold for GWAS p-value for SNP variants corresponding to the chosen gene (default = %default)"),
  make_option("--gwas_response_type",type="character",help="The response variable type for GWAS (quant or cc)"),
  make_option("--gwas_data_name",type="character",help="Name of GWAS dataset (e.g. GLGC, GIANT)"),
  make_option("--gwas_file",type="character",help="File for tab separated GWAS summary statistics data (with header) for chosen trait with column names SNP, BP, CHR, A1, A2, BETA/OR (for case-control), SE, P, A1FREQ, NCASES (for case-control) with path"),
  make_option("--genes_file",type="character",help="File (with path) for tab separated list of chosen genes with column names =  Gene, ENSG_gene, gene_start_position, gene_stop_position, chromosome (GRCh 37)"),
  make_option("--gene_of_interest",type="character",help="ENSG_gene of interest from genes_file (ignore decimal point)"),
  make_option("--gene_boundary_window_size",type="integer",default=1000000,help="The window chosen around the chosen gene from the TSS and TES of the gene (default = %default)"),
  make_option("--eqtl_p_threshold",type="numeric",default=1e-03,help="The threshold for eqtl nominal p-value for SNP variants corresponding to the chosen gene (default = %default)"),
  make_option("--eqtl_file",type="character",help="File for eQTL summary statistics (from GTEx v8) for chosen tissue, split by chromosome, with path"),
  make_option("--cojo_maf",type="numeric",default=0.01,help="MAF threshold used in GCTA-COJO (default = %default)"),
  make_option("--output_folder",type="character",help="Path of the output folder"),
  make_option("--reference_folder",type="character",help="Path of the folder with plink files for LD calculation in gcta (should have chromosome number in the filename in .chromosome.bim/bed/fam format)"),
  make_option("--eqtl_sample_size",type="character",help="Filename (with path) of sample sizes for eQTL datasets across different tissues; has two columns corresponding to tissue name and sample size"),
  make_option("--coloc_p1",type="numeric",default=1e-04,help="Prior probability a SNP is associated with GWAS trait (default = %default)"),
  make_option("--coloc_p2",type="numeric",default=1e-04,help="Prior probability a SNP is associated with gene expression (default = %default)"),
  make_option("--coloc_p12",type="numeric",default=1e-06,help="Prior probability a SNP is associated with GWAS trait and gene expression (default = %default)"),
  make_option("--core",type="integer",default=10,help="Number of cores to run parallel tasks (default = %default)"),
  make_option("--lead_snps_file",type="character",default=NA,help="External Input file of GWAS lead varID, i.e. chr:bp (default = %default)"),
  make_option("--liftover_filename",type="character",help="Filename (with path) to UCSC chainfiles for liftover from GRCh 38 in GTEx v8 to GRCh 37")
)

opts <- parse_args(OptionParser(option_list=opt_list))

#2.1 parse values
chr <- opts$chromosome
lead.snp.window.size <- opts$lead_snp_window_size
gene.boundary.window.size <- opts$gene_boundary_window_size
gwas.p.threshold <- opts$gwas_p_threshold
eqtl.p.threshold <- opts$eqtl_p_threshold
gwas.file <- opts$gwas_file
gwas.response.type <- opts$gwas_response_type
eqtl.file <- opts$eqtl_file
genes.file <- opts$genes_file
cojo.maf <- opts$cojo_maf
gwas.data.name <- opts$gwas_data_name
trait <- opts$trait
output.folder <- opts$output_folder
reference.folder <- opts$reference_folder
tissue <- opts$tissue
gene.of.interest <- opts$gene_of_interest
N.eqtl.file <- opts$eqtl_sample_size
coloc.p1 <- opts$coloc_p1
coloc.p2 <- opts$coloc_p2
coloc.p12 <- opts$coloc_p12
core <- opts$core
lead.snps.file <- opts$lead_snps_file
liftover.filename <- opts$liftover_filename

#2.2 Read and prepare data files

gwas <- as.data.frame(data.table::fread(gwas.file,header=T,stringsAsFactors=FALSE,sep="\t"))
gwas <- gwas[which(gwas[,"CHR"]==chr),]
if(length(grep("OR",colnames(gwas),ignore.case=TRUE))==1&length(grep("BETA",colnames(gwas),ignore.case=TRUE))==0&gwas.response.type=="cc")
  gwas$BETA <- log(gwas$OR)
if(length(grep("zstat",colnames(gwas),ignore.case=TRUE))==1&length(grep("BETA",colnames(gwas),ignore.case=TRUE))==0&gwas.response.type=="cc")
    gwas$BETA <- gwas$ZSTAT/gwas$SE
gwas <- gwas[order(gwas[,"BP"],decreasing=FALSE),]
gwas[,"P"] <- as.numeric(gwas[,"P"])
gwas[,"varID"] <- paste0(gwas[,"CHR"],":",gwas[,"BP"])
gwas.fin <- gwas[which(gwas[,"P"]<gwas.p.threshold),]
if(nrow(gwas.fin)==0) stop("No sig GWAS hits for given trait  in this chromosome")

genes.list <- as.data.frame(data.table::fread(genes.file,header=T,stringsAsFactors=FALSE,sep="\t"))
genes.list <- genes.list[which(genes.list[,"chromosome"]==chr),]
eqtl <- as.data.frame(data.table::fread(eqtl.file,header=T,stringsAsFactors=FALSE))
reference.files <- list.files(reference.folder)
reference.files <- reference.files[grep(paste0("[.]",chr,"[.]"),reference.files)]
reference.filename <- unique(substr(reference.files, 1, nchar(reference.files)-4))
reference.file.bim = reference.files[grep(".bim",reference.files)]
reference.snps <- as.data.frame(data.table::fread(paste0(reference.folder,"/",reference.file.bim),header=F,sep="\t"))
genes.fin <- genes.list[which(genes.list[,"ENSG_gene"]==gene.of.interest),]
if(nrow(genes.fin)==0) stop("Genes(s) of interest not in file")
N.eqtl <- read.table(N.eqtl.file,header=T,stringsAsFactors=FALSE)
N.eqtl.tissue <- N.eqtl[which(N.eqtl[,1]==tissue),2]
eqtl.gene.id.split <- as.character(unlist(strsplit(eqtl[,"gene_id"],"[.]")))
eqtl$eqtl_gene_id <- eqtl.gene.id.split[grep("ENSG",eqtl.gene.id.split)]
gene_fin <- genes.fin[,"ENSG_gene"]
gene_name <- genes.fin[,"Gene"]
message(tissue," / ", trait," / ",gwas.data.name," / ",gene_name," / ", gene_fin, " / chr=",chr)


############################################
## 3. Isolate relevant SNPs for analysis
#############################################

#3.1 Map eQTL signif_variant_gene_pairs and all_pairs SNPs for Hg38 to Hg19 using liftover

dir.create(paste0(output.folder,"/intermediate"), recursive=T)
setwd(output.folder)
eqtl2 <- eqtl[which(eqtl$eqtl_gene_id==gene_fin),]
eqtl_chrpos_alleles <-  as.character(unlist(strsplit(as.character(unlist(strsplit(as.character(unlist(eqtl2[,"variant_id"])),"_b38"))),"_")))
eqtl_chr <- eqtl_chrpos_alleles[grep("chr",eqtl_chrpos_alleles)]
eqtl_chr2 <- substr(eqtl_chr,4,nchar(eqtl_chr))
if(length(unique(eqtl_chr2))==0) stop("eQTL allpairs dataset does not correspond to the chromosome specified")
eqtl_pos_alleles <- eqtl_chrpos_alleles[-grep("chr",eqtl_chrpos_alleles)]
eqtl_pos <- eqtl_pos_alleles[-grep("A|C|T|G|a|c|t|g",eqtl_pos_alleles)]
eqtl_alleles <- eqtl_pos_alleles[grep("A|C|T|G|a|c|t|g",eqtl_pos_alleles)]
eqtl_A1 <- eqtl_alleles[-which(c(1:length(eqtl_alleles))%%2==0)]
eqtl_A2 <- eqtl_alleles[which(c(1:length(eqtl_alleles))%%2==0)]
eqtl_for_liftover <- cbind(eqtl_chr,eqtl_pos,as.numeric(eqtl_pos)+1,eqtl_pos)
write.table(eqtl_for_liftover,file=paste0("intermediate/GTEx_v8_eQTL_",tissue,"_chr",chr,"_",gene_fin,"_allpairs_for_liftover_from_Hg38_to_Hg19.bed"),row.names=F,col.names=F,quote=F,sep="\t")
invisible(system(paste0("module load liftOver; liftOver intermediate/GTEx_v8_eQTL_",tissue,"_chr",chr,"_",gene_fin,"_allpairs_for_liftover_from_Hg38_to_Hg19.bed ",liftover.filename," intermediate/GTEx_v8_eQTL_",tissue,"_chr",chr,"_",gene_fin,"_allpairs_mapped_to_HG19_liftover_output.bed intermediate/GTEx_v8_eQTL_",tissue,"_chr",chr,"_",gene_fin,"_allpairs_mapped_to_HG19_liftover_unlifted.bed"),intern=TRUE))

#3.2 Get lead SNPs for GWAS dataset (p-value < chosen threshold) that lie within (default=1MB) of chosen gene boundaries

if(is.na(lead.snps.file)) lead_snps_gwas <- gwas.fin[,"varID"] else {
  lead_snps_gwas <- read.table(lead.snps.file,stringsAsFactors=FALSE)$V1
  lead_chrpos <- as.character(unlist(strsplit(lead_snps_gwas,":")))
  lead_chr <- as.numeric(lead_chrpos[-which(c(1:length(lead_chrpos))%%2==0)])
  lead_snps_gwas <- lead_snps_gwas[which(lead_chr==chr)]
  lead_snps_gwas <- lead_snps_gwas[which(lead_snps_gwas%in%paste0(gwas.fin[,1],":",gwas.fin[,2]))]
}
lead_chrpos <- as.character(unlist(strsplit(lead_snps_gwas,":")))
lead_chr <- as.numeric(lead_chrpos[-which(c(1:length(lead_chrpos))%%2==0)])
lead_pos <- as.numeric(lead_chrpos[which(c(1:length(lead_chrpos))%%2==0)])

lead_snps2 <- NA
for(i in 1:nrow(genes.fin)){
  tmp <- unique(lead_snps_gwas[which(lead_pos>=(genes.fin[i,"gene_start_position"]-gene.boundary.window.size)&lead_pos<=(genes.fin[i,"gene_stop_position"]+gene.boundary.window.size))])
  lead_snps2 <- c(lead_snps2,tmp)
  #print(i)
}
lead_snps2 <-lead_snps2[-1]
lead_snps_gwas <- sort(unique(lead_snps2))

#3.3 Get lead SNPs for eQTL dataset (p-value < chosen threshold)

eqtl.full.file <- read.table(paste0("intermediate/GTEx_v8_eQTL_",tissue,"_chr",chr,"_",gene_fin,"_allpairs_mapped_to_HG19_liftover_output.bed"),stringsAsFactors=FALSE)
eqtl.varid <- paste0(substr(eqtl.full.file[,1],4,nchar(eqtl.full.file[,1])),":",eqtl.full.file[,2])
eqtl.match <- match(paste0(eqtl.full.file[,1],":",eqtl.full.file[,4]),paste0(eqtl_chr,":",eqtl_pos))
eqtl.fin <- cbind(eqtl.full.file[eqtl.match,],"CHR"=chr,"BP"=eqtl.full.file[,2],"BP_b38"=eqtl.full.file[,4],"varID"=eqtl.varid[eqtl.match],"A1"=eqtl_A1[eqtl.match],"A2"=eqtl_A2[eqtl.match],"maf"=eqtl2[eqtl.match,"maf"],"slope"=eqtl2[eqtl.match,"slope"],"slope_se"=eqtl2[eqtl.match,"slope_se"],"pval_nominal"=eqtl2[eqtl.match,"pval_nominal"])
lead_snps_eqtl <- eqtl.fin[which(eqtl.fin[,"pval_nominal"]<eqtl.p.threshold),"varID"]

#3.4 Get final list of lead SNPs

lead_snps <- intersect(lead_snps_gwas,lead_snps_eqtl)
#if(length(lead_snps)==0) stop("No lead SNPs at the chosen p-value thresholds")
if(length(lead_snps)==0) lead_snps <- lead_snps_gwas

###3.4.1 Function to get lead SNPs at a distance of (default=100KB) from each other

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

#3.5 Count number of SNPs (significant vs total) that lie within (chosen window size) of chosen gene's boundaries

data2 <- gwas[which(gwas[,"P"]<5e-08),]
genes.fin$num.snps=genes.fin$num.sig.snps=NA

for(i in 1:nrow(genes.fin)){
  genes.fin[i,"num.sig.snps"] <- nrow(data2[which(data2[,"BP"]>=(genes.fin[i,"gene_start_position"]-gene.boundary.window.size)&data2[,"BP"]<=(genes.fin[i,"gene_stop_position"]+gene.boundary.window.size)),])
  genes.fin[i,"num.snps"] <- nrow(gwas[which(gwas[,"BP"]>=(genes.fin[i,"gene_start_position"]-gene.boundary.window.size)&gwas[,"BP"]<=(genes.fin[i,"gene_stop_position"]+gene.boundary.window.size)),])
  # print(i)
}

if(length(gene_fin)==0) stop("No genes specified")

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

if(gwas.response.type=="quant") coloc_probs_header = c("Trait","Tissue","GWAS.Data","CHR","Lead.SNP","Lead.BP","GWAS.P","GWAS.Beta","GWAS.SE","GWAS.N","GWAS.A1","GWAS.A2","GWAS.Freq","GWAS.MAF","ENSG_gene","Gene","TopSNP_eQTL","Gene.BPRange","Gene.Start","Gene.End","Gene.Num.SNPs","Gene.Num.GWAS.Sig.SNPs","eQTL","eQTL.Beta","eQTL.SE","eQTL.P","eQTL.N","eQTL.A1","eQTL.A2","eQTL.MAF","Coloc.Ratio","Coloc.H0","Coloc.H1","Coloc.H2","Coloc.H3","Coloc.H4","Coloc2.Ratio","Coloc2.H0","Coloc2.H1","Coloc2.H2","Coloc2.H3","Coloc2.H4")
if(gwas.response.type=="cc") coloc_probs_header = c("Trait","Tissue","GWAS.Data","CHR","Lead.SNP","Lead.BP","GWAS.P","GWAS.Beta","GWAS.SE","GWAS.N","GWAS.NCASES","GWAS.A1","GWAS.A2","GWAS.Freq","GWAS.MAF","ENSG_gene","Gene","TopSNP_eQTL","Gene.BPRange","Gene.Start","Gene.End","Gene.Num.SNPs","Gene.Num.GWAS.Sig.SNPs","eQTL","eQTL.Beta","eQTL.SE","eQTL.P","eQTL.N","eQTL.A1","eQTL.A2","eQTL.MAF","Coloc.Ratio","Coloc.H0","Coloc.H1","Coloc.H2","Coloc.H3","Coloc.H4","Coloc2.Ratio","Coloc2.H0","Coloc2.H1","Coloc2.H2","Coloc2.H3","Coloc2.H4")

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
      coloc <- tryCatch.W.E(coloc.abf(dataset1=list(pvalues=coloc_data[,"gwas.p"], N=round(median(coloc_data[,"gwas.n"]),digits=0), MAF=coloc_data[,"gwas.maf"],type=gwas.response.type), dataset2=list(pvalues=coloc_data[,"eqtl.p"], N=median(coloc_data[,"eqtl.n"]), type="quant",MAF=coloc_data[,"eqtl.maf"]),p1=coloc.p1,p2=coloc.p2,p12=coloc.p12)$summary)
    }
    if(gwas.response.type=="cc") {
      coloc_data <- data.frame("gwas"=gwas.coloc[,c("SNP","maf","p","n","b","se","s")],"eqtl"=eqtl.coloc[,c("SNP","maf","p","n","b","se")])
      colnames(coloc_data) = c("gwas.SNP","gwas.maf","gwas.p","gwas.n","gwas.b","gwas.se","gwas.s","eqtl.SNP","eqtl.maf","eqtl.p","eqtl.n","eqtl.b","eqtl.se")
      if(length(which(is.na(coloc_data[,"eqtl.p"])|is.na(coloc_data[,"gwas.p"])))>0) coloc_data <- coloc_data[-which(is.na(coloc_data[,"eqtl.p"])|is.na(coloc_data[,"gwas.p"])),]
      coloc <- tryCatch.W.E(coloc.abf(dataset1=list(pvalues=coloc_data[,"gwas.p"], N=round(median(coloc_data[,"gwas.n"]),digits=0), MAF=coloc_data[,"gwas.maf"],s=coloc_data[,"gwas.s"],type=gwas.response.type), dataset2=list(pvalues=coloc_data[,"eqtl.p"], N=median(coloc_data[,"eqtl.n"]), type="quant",MAF=coloc_data[,"eqtl.maf"]),p1=coloc.p1,p2=coloc.p2,p12=coloc.p12)$summary)
    }
    warn <- inherits(coloc$warning,"warning")
    error <- inherits(coloc$value,"error")
  }
  if(coloc$value[[5]]<0.5&coloc$value[[6]]>0.5) use.conditional.probs <- TRUE else use.conditional.probs <- FALSE

  ########################################################################################################
  ## 6. Condition SNPs on the top eQTL-associated GWAS SNP using GCTA-COJO and check for secondary signals
  ########################################################################################################

  #6.1 Get topsnp for eQTL data

  coloc_topsnps <- coloc_data[which(coloc_data[,"gwas.p"]<gwas.p.threshold&coloc_data[,"eqtl.p"]<eqtl.p.threshold),]
  if(nrow(coloc_topsnps)>0) eqtl_topsnp <- coloc_topsnps[which(coloc_topsnps[,"eqtl.p"]==min(coloc_topsnps[,"eqtl.p"])),"eqtl.SNP"][1]  else eqtl_topsnp=NA

  if(use.conditional.probs){
    rm(gwas.coloc,eqtl.coloc)
    print("check for secondary signals")

    gwas.matched.ref <- match(x_gwas[,"varID"],paste0(reference.snps[,1],":",reference.snps[,4]))
    eqtl.matched.ref <- match(x_eqtl[,"varID"],paste0(reference.snps[,1],":",reference.snps[,4]))

    x_gwas.bp <- reference.snps[gwas.matched.ref,4]
    x_gwas["varID"] <- paste0(chr,":",x_gwas.bp)
    x_gwas["SNP"] <- reference.snps[gwas.matched.ref,2]
    x_gwas_gcta <- x_gwas[,c("SNP","A1","A2","A1FREQ","BETA","SE","P","N")]
    if(length(which(is.na(x_gwas_gcta[,"SNP"])))>0) x_gwas_gcta <- x_gwas_gcta[-which(is.na(x_gwas_gcta[,"SNP"])),]
    x_eqtl <- cbind("SNP"=reference.snps[eqtl.matched.ref,2],x_eqtl[,-1])
    x_eqtl_gcta <- cbind(x_eqtl[match(x_gwas_gcta[,"SNP"],x_eqtl[,"SNP"]),c("SNP","A1","A2","maf","slope","slope_se","pval_nominal")],"N"=N.eqtl.tissue)
    write.table(x_gwas_gcta,file=paste0("intermediate/GWAS_",trait,"_",gwas.data.name,"_chr",chr,"_pos_",pos,"_",gene_name,".txt"),row.names=F,col.names=T,quote=F,sep="\t")
    if(nrow(x_eqtl_gcta)==0) {message("tissue=",tissue," / j=",j," / ",trait," / ",gwas.data.name,"  / num.lead.snps=",length(lead_snps)) ; next}

    eqtl_topsnp2 <- x_gwas[which(x_gwas[,"varID"]==eqtl_topsnp),"SNP"]
    write.table(x_eqtl_gcta,file=paste0("intermediate/eQTL_",tissue,"_chr",chr,"_pos_",pos,"_",gene_name,".txt"),row.names=F,col.names=T,quote=F,sep="\t")
    write.table(eqtl_topsnp2,file=paste0("intermediate/",trait,"_",gwas.data.name,"_chr",chr,"_pos_",pos,"_",gene_name,"_",tissue,"_topsnp.txt"),row.names=F,col.names=T,quote=F,sep="\t")

    #6.1 Run GCTA-COJO to perform association analysis of all the included SNPs conditional on the given top SNP in the GWAS dataset. Results are saved in the .cma.cojo file

    invisible(system(paste0("module load gcta ;gcta64 --bfile ",reference.folder,"/",reference.filename," --chr ",chr," --maf ",cojo.maf," --cojo-file intermediate/GWAS_",trait,"_",gwas.data.name,"_chr",chr,"_pos_",pos,"_",gene_name,".txt --cojo-cond intermediate/",trait,"_",gwas.data.name,"_chr",chr,"_pos_",pos,"_",gene_name,"_",tissue,"_topsnp.txt --out intermediate/GWAS_",trait,"_",gwas.data.name,"_chr_",chr,"_pos_",pos,"_",gene_name,"_conditionalP"),intern=TRUE))
    tmp.1 <- tryCatch.W.E(read.table(paste0("intermediate/GWAS_",trait,"_",gwas.data.name,"_chr_",chr,"_pos_",pos,"_",gene_name,"_conditionalP.cma.cojo"),header=T))
    warn.1 <- inherits(tmp.1$warning,"warning")
    error.1 <- inherits(tmp.1$value,"error")
    if(error.1==FALSE){
      gwas.coloc <- tmp.1$value
    }

    #6.2 Run GCTA-COJO to perform association analysis of all the included SNPs conditional on the given top SNP in the eQTL dataset. Results are saved in the .cma.cojo file

    invisible(system(paste0("module load gcta ;gcta64 --bfile ",reference.folder,"/",reference.filename," --chr ",chr," --maf ",cojo.maf," --cojo-file intermediate/eQTL_",tissue,"_chr",chr,"_pos_",pos,"_",gene_name,".txt --cojo-cond intermediate/",trait,"_",gwas.data.name,"_chr",chr,"_pos_",pos,"_",gene_name,"_",tissue,"_topsnp.txt --out intermediate/",tissue,"_chr_",chr,"_pos_",pos,"_",gene_name,"_conditionalP"),intern=TRUE))
    tmp.2 <-  tryCatch.W.E(read.table(paste0("intermediate/",tissue,"_chr_",chr,"_pos_",pos,"_",gene_name,"_conditionalP.cma.cojo"),header=T))
    warn.2 <- inherits(tmp.2$warning,"warning")
    error.2 <- inherits(tmp.2$value,"error")
    if(error.2==FALSE){
      eqtl.coloc <- tmp.2$value
    }

    if(error.1==FALSE&error.2==FALSE){
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
        coloc_data2 <- data.frame("gwas"=gwas.coloc[,c("SNP","maf","pC","n","bC","bC_se")],"eqtl"=eqtl.coloc[,c("SNP","maf","pC","n","bC","bC_se")])
        colnames(coloc_data2) = c("gwas.SNP","gwas.maf","gwas.p","gwas.n","gwas.b","gwas.se","eqtl.SNP","eqtl.maf","eqtl.p","eqtl.n","eqtl.b","eqtl.se")
        if(length(which(is.na(coloc_data2[,"eqtl.p"])|is.na(coloc_data2[,"gwas.p"])))>0) coloc_data2 <- coloc_data2[-which(is.na(coloc_data2[,"eqtl.p"])|is.na(coloc_data2[,"gwas.p"])),]
        coloc2 <- tryCatch.W.E(coloc.abf(dataset1=list(pvalues=coloc_data2[,"gwas.p"], N=round(median(coloc_data2[,"gwas.n"]),digits=0), MAF=coloc_data2[,"gwas.maf"],type=gwas.response.type), dataset2=list(pvalues=coloc_data2[,"eqtl.p"], N=median(coloc_data2[,"eqtl.n"]), type="quant",MAF=coloc_data2[,"gwas.maf"]),p1=coloc.p1,p2=coloc.p2,p12=coloc.p12)$summary)
      }
      if(gwas.response.type=="cc") {
        coloc_data2 <- data.frame("gwas"=gwas.coloc[,c("SNP","maf","pC","n","bC","bC_se","s")],"eqtl"=eqtl.coloc[,c("SNP","maf","pC","n","bC","bC_se")])
        colnames(coloc_data2) = c("gwas.SNP","gwas.maf","gwas.p","gwas.n","gwas.b","gwas.se","gwas.s","eqtl.SNP","eqtl.maf","eqtl.p","eqtl.n","eqtl.b","eqtl.se")
        if(length(which(is.na(coloc_data2[,"eqtl.p"])|is.na(coloc_data2[,"gwas.p"])))>0) coloc_data2 <- coloc_data2[-which(is.na(coloc_data2[,"eqtl.p"])|is.na(coloc_data2[,"gwas.p"])),]
        coloc2 <- tryCatch.W.E(coloc.abf(dataset1=list(pvalues=coloc_data2[,"gwas.p"], N=round(median(coloc_data2[,"gwas.n"]),digits=0), MAF=coloc_data2[,"gwas.maf"],s=coloc_data2[,"gwas.s"],type=gwas.response.type), dataset2=list(pvalues=coloc_data2[,"eqtl.p"], N=median(coloc_data2[,"eqtl.n"]), type="quant",MAF=coloc_data2[,"gwas.maf"]),p1=coloc.p1,p2=coloc.p2,p12=coloc.p12)$summary)
      }
      warn2 <- inherits(coloc2$warning,"warning")
      error2 <- inherits(coloc2$value,"error")
      if(coloc2$value[[6]]<0.2) print("no secondary signals")
      if(coloc2$value[[5]]<0.5&coloc2$value[[6]]>0.2) print("possible secondary signals")
      if(coloc2$value[[5]]>0.5&coloc2$value[[6]]>0.2) print("no secondary signals")
    }
  }
  ############################
  ## 8. Prepare output dataset
  ############################

  coloc_probs[,"Lead.SNP"] <- snp
  coloc_probs[,"Trait"] <- trait
  coloc_probs[,"Tissue"] <- tissue
  coloc_probs[,"CHR"] <- chr
  coloc_probs[,"Lead.BP"] <- pos
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
  coloc_probs[,"eQTL"] <- paste0("chr",chr,"_",x_eqtl[which(x_eqtl[,"varID"]==varid),"BP_b38"],"_",x_eqtl[which(x_eqtl[,"varID"]==varid),"A1"],"_",x_eqtl[which(x_eqtl[,"varID"]==varid),"A2"],"_b38")
  coloc_probs[,"eQTL.Beta"] <- x_eqtl[which(x_eqtl[,"varID"]==varid),"slope"][1]
  coloc_probs[,"eQTL.SE"] <- x_eqtl[which(x_eqtl[,"varID"]==varid),"slope_se"][1]
  coloc_probs[,"eQTL.P"] <- x_eqtl[which(x_eqtl[,"varID"]==varid),"pval_nominal"][1]
  coloc_probs[,"eQTL.A1"] <- x_eqtl[which(x_eqtl[,"varID"]==varid),"A1"][1]
  coloc_probs[,"eQTL.A2"] <- x_eqtl[which(x_eqtl[,"varID"]==varid),"A2"][1]
  coloc_probs[,"eQTL.MAF"] <- coloc_data[which(coloc_data[,"eqtl.SNP"]==varid),"eqtl.maf"]
  coloc_probs[,"Gene.BPRange"] <- genes.fin[which(genes.fin[,"ENSG_gene"]==gene_fin),"gene_stop_position"]-genes.fin[which(genes.fin[,"ENSG_gene"]==gene_fin),"gene_start_position"]
  coloc_probs[,"Gene.Num.SNPs"] <- genes.fin[which(genes.fin[,"ENSG_gene"]==gene_fin),"num.snps"]
  coloc_probs[,"Gene.Num.GWAS.Sig.SNPs"] <- genes.fin[which(genes.fin[,"ENSG_gene"]==gene_fin),"num.sig.snps"]
  coloc_probs[,"GWAS.N"] <-  x_gwas[which(x_gwas[,"varID"]==varid),"N"][1]
  coloc_probs[,"eQTL.N"] <-  x_eqtl[which(x_eqtl[,"varID"]==varid),"N"][1]
  if(gwas.response.type=="cc")   coloc_probs[,"GWAS.NCASES"] <-  x_gwas[which(x_gwas[,"varID"]==varid),"NCASES"][1]
  coloc_probs[,"Gene.Start"] <- genes.fin[which(genes.fin[,"ENSG_gene"]==gene_fin),"gene_start_position"]
  coloc_probs[,"Gene.End"] <- genes.fin[which(genes.fin[,"ENSG_gene"]==gene_fin),"gene_stop_position"]

    if(error==FALSE){
    coloc_probs[,"Coloc.Ratio"] = signif(coloc$value["PP.H4.abf"]/(coloc$value["PP.H3.abf"]+coloc$value["PP.H4.abf"]),digits=4)
    coloc_probs[,c("Coloc.H0","Coloc.H1","Coloc.H2","Coloc.H3","Coloc.H4")] = signif(coloc$value[-1],digits=4)
  }

  if(use.conditional.probs)
    if(error.1==FALSE&error.2==FALSE){
      if(error2==FALSE){
        coloc_probs[,"Coloc2.Ratio"] = signif(coloc2$value["PP.H4.abf"]/(coloc2$value["PP.H3.abf"]+coloc2$value["PP.H4.abf"]),digits=4)
        coloc_probs[,c("Coloc2.H0","Coloc2.H1","Coloc2.H2","Coloc2.H3","Coloc2.H4")] = signif(coloc2$value[-1],digits=4)
        rm(error2, gwas.coloc,eqtl.coloc)
      }
  }
  message("tissue=",tissue," / lead.snp.num=",j," / lead.snp=",varid," / eqtl_topsnp=", eqtl_topsnp," / " ,trait," / ",gwas.data.name," / num.lead.snps=",length(lead_snps)," / gene=",gene_name)

  write.table(coloc_probs, file=paste0("chr",chr,"_",gwas.data.name,"_",trait,"_",tissue,"_",gene_name,"_colocProbs.txt"),row.names=F,col.names=F,quote=F,sep="\t",append=T)
  rm(error,coloc_probs,coloc_data,x_eqtl,x_gwas)
} ### end for loop
