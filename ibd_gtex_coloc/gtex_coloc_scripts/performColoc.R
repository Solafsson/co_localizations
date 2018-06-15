
# Usage: This script carries out co-localization between cd/uc and gene expression data from GTex.
#        Rscript performColoc.R <IBD subtype> <Tissue> 
#
# Input: IBD subtype should be either cd or uc. Tissue should be a GTex tissue name like "Colon_Transverse"
#        There are some hard-coded paths in this script, ibd_gwas, ibd_sumstats_path, gtex_sumstats_path and output_dir
#        These could be read in as arguments, but that is not implemented yet. 

# Output: A file titled <IBD_subtype>_<Tissue>_colocalization.txt is written to the output directory. This file contains
#        the columns "Trait1", "Trait2", "Top_SNP_Trait1","n_snps","H0","H1","H2","H3","H4" where
#        Trait1 is the IBD_subtype, Trait2 is a geneID, Top_SNP_Trait1 is an rs_id from the ibd_gwas file, n_snps is the 
#        number of SNPs found in the region for both traits and that are used for the calculations. H0-4 are the posterior
#        Probabilities as described in http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004383

## Run this with the latest R version (R-3.3.0). So instead of just typing "R" in UNIX, type "R-3.3.0"

## Algorithm to be used:
# For each chromosome: 
#   Read in IBD sumstats and eqtl sumstats for that chromosome.
#   Find the IBD variants belonging to the chromosome. For each IBD variant:
#     find the borders of the LD block and extract the variants in the ibd sumstats overlapping this region
#     Join with the eqtl sumstats and split the resulting file into genes. 
#     Test the co-localization with each gene separately

# Note that this means that if the cis-eQTLs for a gene only partially overlaps the LD block surrounding the candidate SNP, 
# then only those SNPs that overlap the block are used in the colocalization testing. What if they co-localize but there are so
# few SNPs that the evidence for colocalization is low? Might want to consider extending the region determining which SNPs are included
# to the left- and rightmost SNP tested for any gene overlapping the LD-block. Not sure how this would affect running times. 
# 



args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  trait="cd"
  Tissue="Colon_Transverse"
  print("No argument supplied, co-localizing with CD as default", call.=FALSE)
} else {
  trait=args[1] #here this is UC or CD
  tissue=args[2]
  #functional_assay=args[2]
}


#load coloc library and its dependencies
#library(colorspace,lib="/software/team152/Rpackages/")
#library(leaps,lib="/software/team152/Rpackages/")
#library(robustbase,lib="/software/team152/Rpackages/")
#library(inline,lib="/software/team152/Rpackages/")
#library(rrcov,lib="/software/team152/Rpackages/")
#library(BMA,lib="/software/team152/Rpackages/")
#library(coloc,lib="/software/team152/Rpackages/")


library(colorspace)
library(leaps)
library(robustbase)
library(inline)
library(rrcov)
library(BMA)
library(coloc)

case_prop=NULL
if( tolower(trait)=="uc"){case_prop=0.354}else
  if(tolower(trait)=="cd"){
    case_prop=0.349 
  }


##firstly, let's get in the list of GWAS loci
ibd_gwas = read.table("/lustre/scratch115/projects/coloc_and_fm/IBD_regions_from_our_2017_paper.txt",head=T,sep="\t") #here we save this data into a variable we decided to call ibd_gwas


## File paths:
ibd_sumstats_path="/lustre/scratch115/projects/coloc_and_fm/splitfiles/ibd_summary_stats/"
gtex_sumstats_path="/lustre/scratch119/humgen/teams/anderson/users/so11/co_localizations/ibd_gtex_coloc/"
output_dir="/lustre/scratch119/humgen/teams/anderson/users/so11/co_localizations/ibd_gtex_coloc/results_long_window/"
hard_window_size=500000
use_hard_windw=TRUE
#tissue="Colon_Transverse"


coloc_results=matrix(ncol=9,nrow=1,0)
coloc_results <- as.data.frame(coloc_results)
colnames(coloc_results)=c("Trait1", "Trait2", "Top_SNP_Trait1","n_snps","H0","H1","H2","H3","H4")

for(chrom in 1:22) {

  print(paste("Testing Chr", chrom))
  #Read in ibd GWAS summary statistics
  ibd_gwas_data=NULL
  if( tolower(trait)=="uc"){
    ibd_gwas_data=read.table(paste(ibd_sumstats_path, "uc_b37_45975.Chr_",chrom,".ready",sep=""),head=T,sep="\t")
  } else if(tolower(trait)=="cd"){
      ibd_gwas_data=read.table(paste(ibd_sumstats_path, "cd_b37_40266.Chr_",chrom,".ready",sep=""),head=T,sep="\t")
  } else {
    stop("Error: Did not recognize the specified trait.")
  }

  
  gtex_gwas_data <- read.table(paste(gtex_sumstats_path, "/", tissue, "/", trait,  "_chr", chrom, "_", tissue, "_forColoc.txt", sep=""), h=T)
  ## Add the variance term
  gtex_gwas_data$expr_var <- gtex_gwas_data$expr_se^2
  ibd_gwas_data$Var <- ibd_gwas_data$StdErr^2
  # Remove missing entries
  gtex_gwas_data <- subset(gtex_gwas_data, !is.na(expr_var) & !is.na(expr_P) & !is.na(expr_beta))
  ibd_gwas_data <- subset(ibd_gwas_data, !is.na(Var) & !is.na(Pvalue) & !is.na(Effect))
  
  snpList <- subset(ibd_gwas, Chr==chrom)
  
  for(i in 1:nrow(snpList)) {
    if(use_hard_windw) {
      pos <- snpList[i, 2]
      region_start <- pos - hard_window_size
      if(region_start < 0) region_start <- 0
      region_end <- pos + hard_window_size
    } else {
      region_start <- snpList[i, 4]
      region_end <- snpList[i, 5]
    }

    rsID <- as.character(snpList[i,3])
    
    ibd_subset <- subset(ibd_gwas_data, Chr==chrom & Position >= region_start & Position <= region_end)
    gtex_subset <- subset(gtex_gwas_data, CHR==chrom & BP >= region_start & BP <= region_end)
    geneList <- unique(as.character(gtex_subset$gene))
    

    ## Do the co-localization 
    for(gene in geneList) {
      gene_subset <- subset(gtex_subset, geneID==gene & as.character(gtex_subset$SNP) %in% as.character(ibd_subset$Marker))
      
      if(nrow(gene_subset) > 0 ) {
        results = coloc.abf(dataset1=list(snp=as.character(ibd_subset$Marker), beta=ibd_subset$Effect, varbeta=ibd_subset$Var, 
                                          type="cc", s=case_prop,MAF=ibd_subset$controls_maf), 
                            dataset2=list(snp=as.character(gene_subset$SNP), beta=gene_subset$expr_beta, 
                                          varbeta=gene_subset$expr_var, pvalues=as.numeric(gene_subset$expr_P), 
                                          N=gene_subset$expr_N, MAF=gene_subset$expr_maf, type="quant"), p12=1e-6)
        
        coloc_results <- rbind(coloc_results, c(trait, gene, rsID, as.numeric(results$summary)))
        
      } 
    }
  }
}

coloc_results <- coloc_results[-1, ]
write.table(coloc_results,paste(output_dir, trait, "_", tissue, "_colocalization.txt", sep=""),quote=F,col.names=T,row.names=F,sep="\t")

