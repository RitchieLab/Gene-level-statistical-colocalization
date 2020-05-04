rm(list=ls())

################## This code uses statistical colocalization to identify genes that may have "LD contamination", i.e. 
#when expression predictor SNPs and GWAS causal SNPs are different but in LD. Here, we find (i)conditionally independent 
#variants in each gene in GWAS and eQTL datasets using "GCTA-COJO" and (ii) estimate the colocalization probabilities P[H3] 
#and P[H4] per gene between GWAS and eQTL datasets using conditional p-values from GCTA-COJO in "coloc".
#### Yogasudha Veturi 15April2020

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
  make_option("--chromosome",type="integer",help="The chromosome to which the given gene corresponds"),
  make_option("--window_size",type="integer",default=1000000,help="The window chosen around the chosen gene from the TSS and TES of the gene (default = %default)"),
  make_option("--gwas_p_threshold",type="numeric",default=1e-03,help="The threshold for GWAS p-value for SNP variants corresponding to the chosen gene (default = %default)"),
  make_option("--eqtl_file",type="character",help="File for eQTL summary statistics (from GTEx) for chosen tissue, split by chromosome, with path"),
  make_option("--gwas_file",type="character",help="File for tab separated GWAS summary statistics data (with header) for chosen trait with column names SNP, BP, CHR, A1, A2, BETA, SE, P, FREQ, with path"),
  make_option("--maf",type="numeric",default=0.01,help="MAF threshold (default = %default)"),
  make_option("--trait",type="character",help="Trait name"),
  make_option("--tissue",type="character",help="Tissue name"),
  make_option("--gwas_data_name",type="character",help="Name of GWAS dataset (e.g. GLGC, GIANT)"),
  make_option("--cojo_p",type="numeric",default=0.001,help="MAF threshold (default = %default)"),
  make_option("--genes_file",type="character",help="File (with path) for tab separated list of chosen genes with column names =  ENSG_gene, gene_start_position, gene_stop_position, chromosome"),
  make_option("--gene_of_interest",type="character",help="ENSG_gene(s) of interest from genes_file (ignore decimal point)"),
  make_option("--output_folder",type="character",help="Path of the output folder"),
  make_option("--ld_folder",type="character",help="Path of the folder with plink files for LD calculation in gcta (should have chromosome number in the filename in .chromosome.bim/bed/fam format)"),
  make_option("--eqtl_sample_size",type="character",help="Filename (with path) of sample sizes for eQTL datasets across different tissues; has two columns corresponding to tissue name and sample size"),
  make_option("--coloc_p1",type="numeric",default=1e-04,help="Prior probability a SNP is associated with GWAS trait (default = %default)"),
  make_option("--coloc_p2",type="numeric",default=0.001,help="Prior probability a SNP is associated with gene expression (default = %default)"),
  make_option("--coloc_p12",type="numeric",default=1e-06,help="Prior probability a SNP is associated with GWAS trait and gene expression (default = %default)"),
  make_option("--core",type="integer",default=10,help="Number of cores to run parallel tasks (default = %default)")
)

opts <- parse_args(OptionParser(option_list=opt_list))

chr <- opts$chromosome
window.size <- opts$window_size
gwas.p.threshold <- opts$gwas_p_threshold
gwas.file <- opts$gwas_file
eqtl.file <- opts$eqtl_file
genes.file <- opts$genes_file
maf <- opts$maf
gwas.data.name <- opts$gwas_data_name
trait <- opts$trait
cojo.p <- opts$cojo_p
output.folder <- opts$output_folder
ld.folder <- opts$ld_folder
tissue <- opts$tissue
gene.of.interest <- opts$gene_of_interest
N.eqtl.file <- opts$eqtl_sample_size
coloc.p1 <- opts$coloc_p1
coloc.p2 <- opts$coloc_p2
coloc.p12 <- opts$coloc_p12
core <- opts$core

#2.2 Read and prepare data files
gwas <- as.data.frame(data.table::fread(gwas.file,header=T,stringsAsFactors=FALSE,sep="\t"))
gwas <- gwas[which(gwas[,"CHR"]==chr),]
gwas <- gwas[order(gwas[,"BP"],decreasing=FALSE),]
gwas[,"P"] <- as.numeric(gwas[,"P"])
gwas[,"SNP2"] <- paste0(gwas[,"CHR"],":",gwas[,"BP"])
gwas.fin <- gwas[which(gwas[,"P"]<gwas.p.threshold),]
if(nrow(gwas.fin)==0) stop("No sig GWAS hits for given trait  in this chromosome")

lead_snps <- gwas.fin[,"SNP2"]
genes.list <- as.data.frame(data.table::fread(genes.file,header=T,stringsAsFactors=FALSE,sep="\t"))
genes.list <- genes.list[which(genes.list[,"chromosome"]==chr),]
eqtl <- as.data.frame(data.table::fread(eqtl.file,header=T,stringsAsFactors=FALSE))
ld.files <- list.files(ld.folder)
ld.files <- ld.files[grep(paste0("[.]",chr,"[.]"),ld.files)]
ld.filename <- unique(substr(ld.files, 1, nchar(ld.files)-4))
ld.file.bim = ld.files[grep(".bim",ld.files)]
ld.snps <- as.data.frame(data.table::fread(paste0(ld.folder,"/",ld.file.bim),header=F,sep="\t"))
genes.fin <- genes.list[grep(gene.of.interest,genes.list[,1]),]
if(nrow(genes.fin)==0) stop("Genes(s) of interest not in file")
N.eqtl <- read.table(N.eqtl.file,header=T,sep="\t",stringsAsFactors=FALSE)
N.eqtl.tissue <- N.eqtl[which(N.eqtl[,1]==tissue),2]


#############################################
## 3. Isolate relevant SNPs for analysis
#############################################

#3.1 Function to get lead SNPs at a distance of 10KB from each other
getleadsnp.fn <- function(data,baseRow,jump=1) {
  success <- TRUE
  x <- data[baseRow,]
  bottomRow <- baseRow

  while (success) {
    # do 
    Distance <- 10000
    startDistance <- x["BP"]
    snp <- paste0(x["CHR"],":",x["BP"])
    bottomRow <- bottomRow+jump
    tmp1 <- abs(data[bottomRow,"BP"]-startDistance)
    success <- (tmp1 < Distance) &  bottomRow<nrow(data)
   # message("DistanceDifference=",tmp1," / bottomRow=", bottomRow)
   }
  return(bottomRow)
}

rows <- i <- 1
while(i<=length(lead_snps)-1) {
  i <- getleadsnp.fn(gwas.fin,baseRow=i)
  rows <- c(rows,i)
}
snp <- NA
 for(i in 1:(length(rows)-1)){
  if(rows[i+1]-rows[i]>1) tmp <- gwas.fin[which(gwas.fin[,"P"]==min(gwas.fin[(rows[i]:(rows[i+1]-1)),"P"])),"SNP2"] else {
     if(i<(length(rows)-1)) tmp <- gwas.fin[rows[i],"SNP2"]; if(i==(length(rows)-1)) tmp <- c(gwas.fin[rows[i],"SNP2"],gwas.fin[rows[i+1],"SNP2"])}
  snp=c(snp,tmp)
 #print(i)
}

lead_snps <- unique(snp[-1])
lead_chrpos <- as.character(unlist(strsplit(lead_snps,":")))
lead_chr <- as.numeric(lead_chrpos[-which(c(1:length(lead_chrpos))%%2==0)])
lead_pos <- as.numeric(lead_chrpos[which(c(1:length(lead_chrpos))%%2==0)])
lead_snps2 <- NA
genes.fin$num.snps=genes.fin$num.sig.snps=NA

#3.2 Count number of SNPs (significant vs total) that lie within (chosen window size) of chosen gene's boundaries

data2 <- gwas[which(gwas[,"P"]<5e-08),]

for(i in 1:nrow(genes.fin)){
    genes.fin[i,"num.sig.snps"] <- nrow(data2[which(data2[,"BP"]>=(genes.list[i,"gene_start_position"]-window.size)&data2[,"BP"]<=(genes.list[i,"gene_end_position"]+window.size)),])
    genes.fin[i,"num.snps"] <- nrow(gwas[which(gwas[,"BP"]>=(genes.list[i,"gene_start_position"]-window.size)&gwas[,"BP"]<=(genes.list[i,"gene_end_position"]+window.size)),])
   # print(i)
}

#3.3 Count number of lead SNPs that lie within (chosen window size) of chosen gene boundaries

for(i in 1:nrow(genes.fin)){
  tmp <- unique(lead_snps[which(lead_pos>=(genes.fin[i,"gene_start_position"]-window.size)&lead_pos<=(genes.fin[i,"gene_end_position"]+window.size))])
  lead_snps2 <- c(lead_snps2,tmp)
  #print(i)
}
lead_snps2 <-lead_snps2[-1]
genes_fin <- genes.fin[,1]
lead_snps <- sort(unique(lead_snps2))
lead_rsid <- gwas.fin[match(lead_snps,paste0(gwas.fin[,"CHR"],":",gwas.fin[,"BP"])),"SNP"]
if(length(genes_fin)==0) stop("No genes specified")

#3.4 Function to find top and bottom rows in the GWAS data file that lie  within the boundaries of the chosen gene (+/- the window size surrounding the gene boundaries)

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

##############################################################################################################################################################
## 4. Run GCTA over lead SNPs to find "conditionally independent" SNPs in GWAS and eQTL datasets; prepare datasets for coloc
##############################################################################################################################################################

dir.create(output.folder, recursive=T)
setwd(output.folder)
coloc_probs_header = c("Trait","Tissue","GWAS.Data","CHR","Lead.SNP","Lead.BP","GWAS.P","GWAS.Beta","GWAS.SE","GWAS.N","GWAS.A1","GWAS.A2","GWAS.Freq","GWAS.MAF","Gene","Gene.BPRange","Gene.Start","Gene.End","Gene.Num.SNPs","Gene.Num.GWAS.Sig.SNPs","eQTL","eQTL.Beta","eQTL.SE","eQTL.P","eQTL.N","eQTL.A1","eQTL.A2","eQTL.MAF","Coloc.Ratio","Coloc.H0","Coloc.H1","Coloc.H2","Coloc.H3","Coloc.H4")
write.table(t(coloc_probs_header), file=paste0("chr",chr,"_",gwas.data.name,"_",trait,"_",tissue,"_colocProbs.txt"),row.names=F,col.names=F,quote=F,sep="\t")
for(j in 1:length(lead_snps)){ # Loop starts for lead snps
     snp <- lead_snps[j]
     rsid <- lead_rsid[j]
     pos <- as.numeric(as.character(unlist(strsplit(snp,":"))))[2]
     genes_fin = unique(genes.fin[which((genes.fin[,"gene_start_position"]-window.size)<=pos&(genes.fin[,"gene_end_position"]+window.size)>=pos),"ENSG_gene"])
     coloc_probs <- matrix(,1,length(coloc_probs_header))
     colnames(coloc_probs) <- coloc_probs_header
     coloc_probs <- data.frame(coloc_probs)
     for(i in 1:ncol(coloc_probs)) coloc_probs[,i] <- as.character(unlist(coloc_probs[,i]))

     base <- which(paste0(gwas[,"CHR"],":",gwas[,"BP"])==snp)[1]
     boundaryRows <- data.frame(getWindow.fn(baseRow=base,data=gwas,size=window.size,jump=100))  ## Get boundary rows for GWAS dataset; these will yield all the GWAS SNPs within chosen window size on either side of lead SNP
     x <- gwas[boundaryRows[[1]]:boundaryRows[[2]],]

     x_gcta <- x
     rsid2 <- ld.snps[match(x_gcta[,"SNP"],ld.snps[,2]),2]
     x_gcta[,"SNP"] <- rsid2
     x_gcta = x_gcta[,c("SNP","A1","A2","FREQ","BETA","SE","P","N")]

     write.table(x_gcta,file=paste0("GWAS_",trait,"_",gwas.data.name,"_chr",chr,"_pos_",pos,".txt"),row.names=F,col.names=T,quote=F,sep="\t")

#4.1 Run GCTA using system commands to get list of independently associated SNPs with GWAS p-value < cojo.p. Results after applying model selection procedure (--cojo-slct) are saved in the .jma.cojo file

     invisible(system(paste0("gcta64 --bfile ../",ld.folder,"/",ld.filename,"  --chr ",chr," --maf ",maf," --cojo-file GWAS_",trait,"_",gwas.data.name,"_chr",chr,"_pos_",pos,".txt --cojo-slct --cojo-p ",cojo.p,"  --out ",gwas.data.name,"_",trait,"_chr",chr,"_pos_",pos,"_pval_",cojo.p),intern=TRUE))
     gwas.topsnps <- read.table(paste0(gwas.data.name,"_",trait,"_chr",chr,"_pos_",pos,"_pval_",cojo.p,".jma.cojo"),header=T,stringsAsFactors=FALSE)
     x_topsnps = gwas.topsnps[which(gwas.topsnps[,"SNP"]%in%x[,"SNP"]),]

#4.2 If there are at least 2 independent SNPs  (other than the lead SNP) that are greater than 100KB from the lead SNP, run GCTA to perform association analysis of the included SNPs conditional on the given list of independent SNPs. Results are saved in the .cma.cojo file

     if(nrow(x_topsnps)>=2) {
       x_topsnps <- x_topsnps[which(x_topsnps[,"SNP"]!=rsid),]; if(nrow(x_topsnps)>=2) x_topsnps <- x_topsnps[which(abs(x_topsnps[,"bp"]-pos)>100000),]
       if(nrow(x_topsnps)>=2) {
         write.table(x_topsnps[,"SNP"],file=paste0("GWAS_",trait,"_",gwas.data.name,"_chr",chr,"_pos_",pos,"_topsnpslist.txt"),row.names=F,col.names=T,quote=F,sep="\t")
         invisible(system(paste0("gcta64 --bfile ../",ld.folder,"/",ld.filename," --chr ",chr," --maf ",maf," --cojo-file GWAS_",trait,"_",gwas.data.name,"_chr",chr,"_pos_",pos,".txt --cojo-cond GWAS_",trait,"_",gwas.data.name,"_chr",chr,"_pos_",pos,"_topsnpslist.txt --out GWAS_",trait,"_",gwas.data.name,"_chr_",chr,"_pos_",pos,"_conditionalP"),intern=TRUE))
         tmp3=tryCatch.W.E(read.table(paste0("GWAS_",trait,"_",gwas.data.name,"_chr_",chr,"_pos_",pos,"_conditionalP.cma.cojo"),header=T))
         warn3 <- inherits(tmp3$warning,"warning")
         error3 <- inherits(tmp3$value,"error")
         if(error3==FALSE){
         gwas.cond <- tmp3$value
         gwas.cond[,"p"] <- gwas.cond[,"pC"]
         }} else {gwas.cond <- x[,-which(colnames(x)=="SNP2")]; colnames(gwas.cond) <- c("SNP","bp","Chr","refA","A2","b","se","p","freq","n")}
      }

#4.3 Repeat steps 4.1 and 4.2 for all "eGenes" to get conditionally independent SNPs for the eQTL dataset instead of GWAS dataset
     for(k in 1:length(genes_fin)){
        x1 <- eqtl[which(eqtl[,"gene_id"]%in%genes_fin[k]),]
        if(nrow(x1)==0) {message("tissue=",tissue," / j=",j," / ",trait," / ",gwas.data.name," / k=", k," / num.lead.snps=",length(lead_snps)) ; next}
        x1_chrpos_alleles <-  as.character(unlist(strsplit(as.character(unlist(strsplit(as.character(unlist(x1[,2])),"_b37"))),"_")))
        x1_chrpos <- x1_chrpos_alleles[-grep("A|C|T|G|a|c|t|g",x1_chrpos_alleles)]
        x1_chr <- as.numeric(x1_chrpos[-which(c(1:length(x1_chrpos))%%2==0)])
        x1_pos <- as.numeric(x1_chrpos[which(c(1:length(x1_chrpos))%%2==0)])
        which.x1 <- which(paste0(x1_chr,":",x1_pos)%in%paste0(x[,"CHR"],":",x[,"BP"]))
        x1_alleles <- x1_chrpos_alleles[grep("A|C|T|G|a|c|t|g",x1_chrpos_alleles)]
        x1_A1 <- x1_alleles[-which(c(1:length(x1_alleles))%%2==0)]
        x1_A2 <- x1_alleles[which(c(1:length(x1_alleles))%%2==0)]

        x2 <- cbind(x1[which.x1,],x1_chr[which.x1],x1_pos[which.x1],x1_A1[which.x1],x1_A2[which.x1],paste0(x1_chr[which.x1],":",x1_pos[which.x1]))
        colnames(x2) <- c(colnames(eqtl),c("CHR","BP","A1","A2","SNP"))
        x2$N <- N.eqtl.tissue
        x2_gcta <- x2[,c("SNP","A1","A2","maf","slope","slope_se","pval_nominal","N")]
        rsids.matched <- ld.snps[match(x2_gcta[,"SNP"],paste0(ld.snps[,1],":",ld.snps[,4])),2]
        x2_gcta[,"SNP"] <- rsids.matched

        write.table(x2_gcta,file=paste0("eQTL_",tissue,"_chr",chr,"_pos_",pos,"_",genes_fin[k],".txt"),row.names=F,col.names=T,quote=F,sep="\t")
        invisible(system(paste0("gcta64 --bfile ../",ld.folder,"/",ld.filename," --chr ",chr," --maf ",maf," --cojo-file eQTL_",tissue,"_chr",chr,"_pos_",pos,"_",genes_fin[k],".txt --cojo-top-SNPs 2 --out ",tissue,"_chr_",chr,"_pos_",pos,"_",genes_fin[k]),intern=TRUE))
        eqtl.topsnps <- read.table(paste0(tissue,"_chr_",chr,"_pos_",pos,"_",genes_fin[k],".jma.cojo"),header=T)
        if(nrow(eqtl.topsnps)>=2) {
          eqtl.topsnps = eqtl.topsnps[which(eqtl.topsnps[,"SNP"]!=rsid),];
          eqtl.topsnps = eqtl.topsnps[which(abs(eqtl.topsnps[,"bp"]-pos)>100000),]
          if(nrow(eqtl.topsnps)>=2) {
            write.table(eqtl.topsnps[,2],file=paste0("eQTL_",tissue,"_chr",chr,"_pos_",pos,"_",genes_fin[k],"_topsnplist.txt"),row.names=F,col.names=T,quote=F,sep="\t")
            invisible(system(paste0("gcta64 --bfile ../",ld.folder,"/",ld.filename," --chr ",chr," --maf ",maf," --cojo-file eQTL_",tissue,"_chr",chr,"_pos_",pos,"_",genes_fin[k],".txt --cojo-cond eQTL_",tissue,"_chr",chr,"_pos_",pos,"_",genes_fin[k],"_topsnplist.txt --out ",tissue,"_chr_",chr,"_pos_",pos,"_",genes_fin[k],"_conditionalP"),intern=TRUE))
            tmp1 <-  tryCatch.W.E(read.table(paste0(tissue,"_chr_",chr,"_pos_",pos,"_",genes_fin[k],"_conditionalP.cma.cojo"),header=T))
            warn1 <- inherits(tmp1$warning,"warning")
            error1 <- inherits(tmp1$value,"error")
            if(error1==FALSE){
              eqtl.cond <- tmp1$value
              eqtl.cond[,"p"] <- eqtl.cond[,"pC"]
           }} else {eqtl.cond <- cbind(x2_gcta[,"SNP"],x2[,c("BP","CHR")],x2_gcta[,c("A1","A2","maf","slope","slope_se","pval_nominal","N")]); colnames(eqtl.cond) = c("SNP","bp","Chr","refA","A2","freq","b","se","p","n") }
        }
        eqtl.topsnps <- read.table(paste0(tissue,"_chr_",chr,"_pos_",pos,"_",genes_fin[k],".jma.cojo"),header=T)
        eqtl.cond2 <- eqtl.cond[which(eqtl.cond[,"SNP"]%in%gwas.cond[,"SNP"]),]
        gwas.cond2 <- gwas.cond[which(gwas.cond[,"SNP"]%in%eqtl.cond2[,"SNP"]),]
        gwas.cond2 <- gwas.cond2[match(eqtl.cond2[,"SNP"],gwas.cond2[,"SNP"]),]
        if(length(which(is.na(gwas.cond2[,"SNP"])))>0) gwas.cond2 = gwas.cond2[-which(is.na(gwas.cond2[,"SNP"])),]
        eqtl.cond2 <- eqtl.cond2[match(gwas.cond2[,"SNP"],eqtl.cond2[,"SNP"]),]
        if(length(which(is.na(eqtl.cond2[,"SNP"])))>0) eqtl.cond2 = eqtl.cond2[-which(is.na(eqtl.cond2[,"SNP"])),]
        stopifnot(all.equal(as.character(unlist(eqtl.cond2[,"SNP"])),as.character(unlist(gwas.cond2[,"SNP"]))))

        eqtl.maf <- ifelse(eqtl.cond2[,"freq"]<0.5,eqtl.cond2[,"freq"],1-eqtl.cond2[,"freq"])
        gwas.maf <- ifelse(gwas.cond2[,"freq"]<0.5,gwas.cond2[,"freq"],1-gwas.cond2[,"freq"])

        if(nrow(eqtl.cond2)>0){

############################################################
## 5. Run colocalization analysis between GWAS SNPs and eQTL
############################################################

#5.1 Prepare eQTL and GWAS datasets for colocalization analysis 

           coloc_data <- data.frame("gwas"=gwas.cond2[,c("freq","p","n","b","se")],"eqtl"=eqtl.cond2[,c("freq","p","n","b","se")])
           coloc_data[,"SNP"] <- gwas.cond2[,"SNP"]
           coloc_data$gwas.maf <- gwas.maf
           coloc_data$eqtl.maf <- eqtl.maf

           coloc <- tryCatch.W.E(coloc.abf(dataset1=list(pvalues=coloc_data[,"gwas.p"], N=round(median(coloc_data[,"gwas.n"]),digits=0),beta = coloc_data[,"gwas.b"], MAF=coloc_data[,"gwas.maf"],type="quant"), dataset2=list(pvalues=coloc_data[,"eqtl.p"], beta=coloc_data[,"eqtl.b"],N=median(coloc_data[,"eqtl.n"]), type="quant",MAF=coloc_data[,"gwas.maf"]),p1=coloc.p1,p2=coloc.p2,p12=coloc.p12)$summary)
           warn <- inherits(coloc$warning,"warning")
           error <- inherits(coloc$value,"error")

#5.2 Prepare output dataset 

           coloc_probs[,"Lead.SNP"] <- rsid
           coloc_probs[,"Trait"] <- trait
           coloc_probs[,"Tissue"] <- tissue
           coloc_probs[,"CHR"] <- chr
           coloc_probs[,"Lead.BP"] <- pos
           coloc_probs[,"GWAS.P"] <- x_gcta[which(x_gcta[,"SNP"]==rsid),"P"][1]
           coloc_probs[,"GWAS.Beta"] <- x_gcta[which(x_gcta[,"SNP"]==rsid),"BETA"][1]
           coloc_probs[,"GWAS.SE"] <- x_gcta[which(x_gcta[,"SNP"]==rsid),"SE"][1]
           coloc_probs[,"GWAS.Data"] <- gwas.data.name
           coloc_probs[,"GWAS.A1"] <- x_gcta[which(x_gcta[,"SNP"]==rsid),"A1"][1]
           coloc_probs[,"GWAS.A2"] <- x_gcta[which(x_gcta[,"SNP"]==rsid),"A2"][1]
           coloc_probs[,"GWAS.Freq"] <- x_gcta[which(x_gcta[,"SNP"]==rsid),"FREQ"][1]
           coloc_probs[,"GWAS.MAF"] <- coloc_data[which(coloc_data[,"SNP"]==rsid),"gwas.maf"]
           coloc_probs[,"Gene"] <- genes_fin[k]
           coloc_probs[,"eQTL"] <- eqtl.topsnps[which(eqtl.topsnps[,"p"]==min(eqtl.topsnps[,"p"])),"SNP"][1]
           coloc_probs[,"eQTL.Beta"] <- x2_gcta[which(x2_gcta[,"SNP"]==rsid),"slope"][1]
           coloc_probs[,"eQTL.SE"] <- x2_gcta[which(x2_gcta[,"SNP"]==rsid),"slope_se"][1]
           coloc_probs[,"eQTL.P"] <- x2_gcta[which(x2_gcta[,"SNP"]==rsid),"pval_nominal"][1]
           coloc_probs[,"eQTL.A1"] <- x2_gcta[which(x2_gcta[,"SNP"]==rsid),"A1"][1]
           coloc_probs[,"eQTL.A2"] <- x2_gcta[which(x2_gcta[,"SNP"]==rsid),"A2"][1]
           coloc_probs[,"eQTL.MAF"] <- coloc_data[which(coloc_data[,"SNP"]==rsid),"eqtl.maf"]
           coloc_probs[,"Gene.BPRange"] <- genes.fin[which(genes.fin[,"ENSG_gene"]==genes_fin[k]),"gene_end_position"]-genes.fin[which(genes.fin[,"ENSG_gene"]==genes_fin[k]),"gene_start_position"]
           coloc_probs[,"Gene.Num.SNPs"] <- genes.fin[which(genes.fin[,"ENSG_gene"]==genes_fin[k]),"num.snps"]
           coloc_probs[,"Gene.Num.GWAS.Sig.SNPs"] <- genes.fin[which(genes.fin[,"ENSG_gene"]==genes_fin[k]),"num.sig.snps"]
           coloc_probs[,"GWAS.N"] <-  x_gcta[which(x_gcta[,"SNP"]==rsid),"N"][1]
           coloc_probs[,"eQTL.N"] <-  x2_gcta[which(x2_gcta[,"SNP"]==rsid),"N"][1]
           coloc_probs[,"Gene.Start"] <- genes.fin[which(genes.fin[,"ENSG_gene"]==genes_fin[k]),"gene_start_position"]
           coloc_probs[,"Gene.End"] <- genes.fin[which(genes.fin[,"ENSG_gene"]==genes_fin[k]),"gene_end_position"]

           if(error==FALSE){
             coloc_probs[,"Coloc.Ratio"] = signif(coloc$value["PP.H4.abf"]/(coloc$value["PP.H3.abf"]+coloc$value["PP.H4.abf"]),digits=4)
             coloc_probs[,c("Coloc.H0","Coloc.H1","Coloc.H2","Coloc.H3","Coloc.H4")] = signif(coloc$value[-1],digits=4)
           }

           message("tissue=",tissue," / lead.snp.num=",j," / gene.num=",k," / num.genes=", length(genes_fin), " / ", trait," / ",gwas.data.name," / num.lead.snps=",length(lead_snps))
           write.table(coloc_probs, file=paste0("chr",chr,"_",gwas.data.name,"_",trait,"_",tissue,"_colocProbs.txt"),row.names=F,col.names=F,quote=F,sep="\t",append=T)
           rm(coloc_data,x1,x2,gwas.cond2,eqtl.cond2,gwas.cond,eqtl.cond)
         }
    }
  rm(x,coloc_probs)
}
                                                                                                                                                                                                                          356,1         Bot
