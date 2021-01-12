################## This code compiles results from run_gene_level_coloc.R across genes to filter out genes with no evidence of colocalization

#### Yogasudha Veturi 12January2021

#############################################
## 1. Load necessary libraries
#############################################

suppressMessages(library(optparse))
suppressMessages(library(data.table))

#############################################
## 2. Define parameters and prepare datasets
#############################################

opt_list <- list(
  make_option("--file_name",type="character",help="Load coloc output file obtained after running run_gene_level_coloc.R"),
  make_option("--coloc_p_h3_threshold",type="numeric",default=0,5,help="Threshold for minP[H3] per gene across all lead SNPs in the gene (default=%default)"),
  make_option("--coloc_p_h4_threshold",type="numeric",default=0.5,help="Threshold for maxP[H4] per gene across all lead SNPs in the gene (default=%default)"),
  make_option("--output_folder",type="character",help="Folder where list of genes is saved")
)

opts <- parse_args(OptionParser(option_list=opt_list))

#2.1 parse values
file.name <- opts$file_name
coloc.p.h3.threshold <- opts$coloc_p_h3_threshold
coloc.p.h4.threshold <- opts$coloc_p_h4_threshold
output.folder <- opts$output_folder

#Compile results
data = as.data.table(read.table(paste0(output.folder,"/",file.name),header=T,sep="\t",stringsAsFactors=FALSE))
h3 = data[data[, .I[Coloc.H3 == min(Coloc.H3)], by=Gene]$V1]
#h3 = h3[-which(duplicated(h3[,"Gene"])),]
h4 = data[data[, .I[Coloc.H4 == max(Coloc.H4)], by=Gene]$V1]
#h4 = h4[-which(duplicated(h4[,"Gene"])),]
gene.list.h3 = h3[which(h3[,"Coloc.H3"]<=coloc.p.h3.threshold),"Gene"]
gene.list.h4 = h4[which(h4[,"Coloc.H4"]>=coloc.p.h4.threshold),"Gene"]
gene.list = intersect(as.data.frame(gene.list.h3),as.data.frame(gene.list.h4))
if(length(gene.list)==0) stop("No genes that pass the LD contamination filter") 
else write.table(paste0(output.folder,"/Final_list_of_genes_no_ld_contamination.txt"),row.names=F,col.names=F,sep="\t",quote=F)
