
## This is an R file, so these commands should be run in an R environment (just type "R" in your unix terminal)
## In this file we will load the GWAS results for PSC, then load the summary statistics for lymphocyte counts from Astle Cell 2017 paper, and then
## perform colocalisation analysis with coloc 
## Feb 2018

## Run this with the latest R version (R-3.3.0). So instead of just typing "R" in UNIX, type "R-3.3.0"


args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("An argument must be supplied (trait: UC or CD) for the dataset to be analysed with lymphocyte counts", call.=FALSE)
}
trait=args[1] #here this is UC or CD
#functional_assay=args[2]

#load coloc library and its dependencies
library(colorspace,lib="/software/team152/Rpackages/")
library(leaps,lib="/software/team152/Rpackages/")
library(robustbase,lib="/software/team152/Rpackages/")
library(inline,lib="/software/team152/Rpackages/")
library(rrcov,lib="/software/team152/Rpackages/")
library(BMA,lib="/software/team152/Rpackages/")
library(coloc,lib="/software/team152/Rpackages/")

case_prop=NULL
if( tolower(trait)=="uc"){case_prop=0.354}else
if(tolower(trait)=="cd"){
case_prop=0.349 
}
case_prop

##first we need to read in our data from GWAS
## we do this using function read.table

##firstly, let's get in the list of GWAS loci
ibd_gwas = read.table("/lustre/scratch115/projects/psc2020/IBD_regions_from_our_2017_paper.txt",head=T,sep="\t") #here we save this data into a variable we decided to call ibd_gwas
#note that sometimes loading data can take time (for instance if our file is 200Mb or more)

#important variables
#snp_of_interest="rs1748195" #in the end, we will do this for all entries but for now let's focus on one SNP


for(j in 1:dim(ibd_gwas)[1]){
#for(j in 1:dim(ibd_gwas)[1]){snp_of_interest=as.character(ibd_gwas[j,1]);print(snp_of_interest)}
snp_of_interest=as.character(ibd_gwas[j,3])

#where to save our results:
output_directory="/lustre/scratch115/projects/psc2020/COLOC_RESULTS/IBD2017/monocyte_counts/"
output_filename=paste(output_directory,snp_of_interest,".",toupper(trait),".monocyte_counts.coloc.results",sep="")

#now, find the entry (row) where this SNP is on our table
index=which(ibd_gwas[,3]==snp_of_interest)

#have a look at this entry:
ibd_gwas[index,]
chrom=ibd_gwas[index,1]
print(paste("chrom: ",chrom,sep=""))
#####################################################################################################################################################
 #now load the data we want from the other GWAS dataset (for the other trait)
second_dataset=read.table(paste("/lustre/scratch115/projects/psc2020/splitfiles/monocyte_counts_summary_stats/mono_N170721_narrow.Chr_",chrom,".sorted.txt",sep=""),head=T,sep="\t")

#this is big, but we only care about our region around the snp we are interested in
region_start = ibd_gwas[index,4] 
region_end=ibd_gwas[index,5] 
print(paste0("region: ",region_start,"-",region_end))

#we will use the fact that our data is sorted, to find the starting_position and ending_position
idx_start=which(as.numeric(second_dataset[,2])>=region_start)[1]

tmp=which(as.numeric(second_dataset[,2])<=region_end)
idx_end=tmp[length(tmp)]
print(paste0("idx_start:",idx_start," idx_end: ",idx_end))

#now that we have that, this is all we care about (from the whole second_dataset data)
second_dataset_subset = second_dataset[idx_start:idx_end,]
#also remove entries where the pvalues are NA
toremove=which(is.na(second_dataset_subset[,4]))
toremove_missingMAF=which(is.na(second_dataset_subset[,9]))
toremove_all=unique(c(toremove,toremove_missingMAF))
if(length(toremove_all)>0){
second_dataset_subset=second_dataset_subset[-toremove_all,]
}
#we can now even remove the big data we loaded earlier
rm(second_dataset)
######################################################################################################################################################

#now load the data we want from the GWAS data.
ibd_gwas_data=NULL
if( tolower(trait)=="uc"){
ibd_gwas_data=read.table(paste("/lustre/scratch115/projects/psc2020/splitfiles/ibd_summary_stats/uc_b37_45975.Chr_",chrom,".ready",sep=""),head=T,sep="\t")
}else
if(tolower(trait)=="cd"){
ibd_gwas_data=read.table(paste("/lustre/scratch115/projects/psc2020/splitfiles/ibd_summary_stats/cd_b37_40266.Chr_",chrom,".ready",sep=""),head=T,sep="\t")
}


#we will use the fact that our data is sorted, to find the starting_position and ending_position
idx_start_gwas=which(as.numeric(ibd_gwas_data[,2])>=region_start)[1]
tmp=which(as.numeric(ibd_gwas_data[,2])<=region_end)
idx_end_gwas=tmp[length(tmp)]

#now that we have that, this is all we care about (from the whole ibd_gwas data)
ibd_gwas_data_subset = ibd_gwas_data[idx_start_gwas:idx_end_gwas,]

#we can now even remove the big data we loaded earlier
rm(ibd_gwas_data)
######################################################################################################################################################

varbeta_gwas=ibd_gwas_data_subset[,6]^2


#addition March 1st: Sometimes there is no overlap at all in the set of SNPs (!). In that case just move on.
snp_gwas=ibd_gwas_data_subset$Marker
snp_func=second_dataset_subset[,3]
tmp=which(as.character(snp_func)%in%as.character(snp_gwas))
if(length(tmp)>0){
##this is how SunGou had run coloc:
#results <- coloc.abf(dataset1=list(N=data$N.x, s=data$s.x, pvalues=data$pvalues.x, type="cc"), dataset2=list(N=data$N.y, s=data$s.y, pvalues=data$pvalues.y, type="cc"), MAF=data$MAF, p12=1e-06) # p12 is set as what Fortune et al. 2015 recommended

#this is the one which ran fine:
results = coloc.abf(dataset1=list(snp=as.character(ibd_gwas_data_subset$Marker),beta=ibd_gwas_data_subset[,5],varbeta=varbeta_gwas,type="cc",s=case_prop,
 MAF=ibd_gwas_data_subset$controls_MAF),
 dataset2=list(snp=as.character(second_dataset_subset[,3]),pvalues=as.numeric(second_dataset_subset[,4]),N=second_dataset_subset$SampleSize,
 MAF=second_dataset_subset$FreqAltAllele, type="quant")
,p12=1e-06)

}
coloc_results=matrix(ncol=6,nrow=1,0)
coloc_results[1,]=results$summary
colnames(coloc_results)=c("n_snps","H0","H1","H2","H3","H4")
write.table(coloc_results,output_filename,quote=F,col.names=T,row.names=F,sep="\t")
#now we know from them the README file from Louella that the beta is column 5 and the SE column 9.

}#here we close the for loop for the snp_of_interest

#########################################



#END
