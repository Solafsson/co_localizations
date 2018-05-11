
## This is an R file, so these commands should be run in an R environment (just type "R" in your unix terminal)
## In this file we will load the GWAS results for PSC, then load the functional analyses results from Blueprint and then
## perform colocalisation analysis with coloc (at first)
## Feb 2017

## Run this with the latest R version (R-3.3.0). So instead of just typing "R" in UNIX, type "R-3.3.0"

args = commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  stop("Three arguments must be supplied. First is cell type of interest, second is Blueprint assay, third is UC/CD", call.=FALSE)
}
celltype_of_interest=args[1]
functional_assay=args[2]
trait=args[3]
print(paste("Colocalisation of ",functional_assay," in ",celltype_of_interest," and ",trait,sep=""))
case_prop=NULL
if( tolower(trait)=="uc"){case_prop=0.354}else
if(tolower(trait)=="cd"){
case_prop=0.349
}
case_prop


#load coloc library and its dependencies
library(colorspace,lib="/software/team152/Rpackages/")
library(leaps,lib="/software/team152/Rpackages/")
library(robustbase,lib="/software/team152/Rpackages/")
library(inline,lib="/software/team152/Rpackages/")
library(rrcov,lib="/software/team152/Rpackages/")
library(BMA,lib="/software/team152/Rpackages/")
library(coloc,lib="/software/team152/Rpackages/")

##first we need to read in our data from GWAS
## we do this using function read.table

##firstly, let's get in the list of GWAS loci
ibd_gwas = read.table("/lustre/scratch115/projects/psc2020/IBD_regions_from_our_2017_paper.txt",head=T,sep="\t") #here we save this data into a variable we decided to call ibd_gwas
#note that sometimes loading data can take time (for instance if our file is 200Mb or more)


#important variables
#snp_of_interest="rs2836883" #in the end, we will do this for all entries but for now let's focus on one SNP


#for(j in 1:dim(ibd_gwas)[1]){snp_of_interest=as.character(ibd_gwas[j,1]);print(snp_of_interest)}
for(j in 1:dim(ibd_gwas)[1]){
snp_of_interest=as.character(ibd_gwas[j,3])
#celltype_of_interest="tcel" #this could also be "tcel" or "mono"
#functional_assay="gene_nor_combat" #this could be "K4ME1_log2rpm" for instance

#where to save our results:
output_directory=paste("/lustre/scratch115/projects/psc2020/COLOC_RESULTS/IBD2017/",celltype_of_interest,"/",functional_assay,"/",sep="")
output_filename=paste(output_directory,snp_of_interest,".",toupper(trait),".",celltype_of_interest,".",functional_assay,".coloc.results",sep="")

#depending on the celltype and assay, we want to get the number of counts correct:
func_counts=NULL
if(celltype_of_interest=="mono"){
	if(functional_assay=="gene_nor_combat"){
		func_counts=194
	}else if(functional_assay=="meth_M"){
		func_counts=196
	}else if(functional_assay=="K4ME1_log2rpm"){
		func_counts=172
	}else if(functional_assay=="K27AC_log2rpm"){
		func_counts=162
	}else{
		func_counts=194
	}		
}else if(celltype_of_interest=="neut"){
	if(functional_assay=="gene_nor_combat"){
		func_counts=192
	}else if(functional_assay=="meth_M"){
		func_counts=197
	}else if(functional_assay=="K4ME1_log2rpm"){
		func_counts=173
	}else if(functional_assay=="K27AC_log2rpm"){
		func_counts=174
	}else{
		func_counts=192
	}
}else if(celltype_of_interest=="tcel"){
	if(functional_assay=="gene_nor_combat"){
		func_counts=171
	}else if(functional_assay=="meth_M"){
		func_counts=133
	}else if(functional_assay=="K4ME1_log2rpm"){
		func_counts=104       
	}else if(functional_assay=="K27AC_log2rpm"){
                func_counts=142
	}else{
		func_counts=171
	}
}
else{
	print("not recognized celltype");
	func_counts=NA
}	

print(paste("func counts: ",func_counts,sep=""))
#now, find the entry (row) where this SNP is on our table
index=which(ibd_gwas[,3]==snp_of_interest)

#have a look at this entry:
ibd_gwas[index,]
chrom=ibd_gwas[index,1]
print(paste("chrom: ",chrom,sep=""))
#####################################################################################################################################################
#now load the data we want from the functional dataset
blueprint_func=NULL
if(functional_assay=="meth_M"){
	fullfile1=paste("/lustre/scratch115/projects/psc2020/splitfiles/",celltype_of_interest,"/full_meth_files/",celltype_of_interest,"_",functional_assay,"_peer_10_all_summary.Chr_",chrom,".sorted.txt",sep="")
	fullfile2=paste("/lustre/scratch115/projects/psc2020/splitfiles/",celltype_of_interest,"/",celltype_of_interest,"_",functional_assay,"_peer_10_all_summary.Chr_",chrom,".sorted.txt",sep="")
	if(file.exists(fullfile1)){
		blueprint_func=read.table(fullfile1,head=F)
	}
	else{
                blueprint_func=read.table(fullfile2,head=F)
	}
}else if(functional_assay=="K27AC_log2rpm"){
	fullfile1=paste("/lustre/scratch115/projects/psc2020/splitfiles/",celltype_of_interest,"/full_K27AC_log2rpm_files/",celltype_of_interest,"_",functional_assay,"_peer_10_all_summary.Chr_",chrom,".sorted.txt",sep="")
	fullfile2=paste("/lustre/scratch115/projects/psc2020/splitfiles/",celltype_of_interest,"/",celltype_of_interest,"_",functional_assay,"_peer_10_all_summary.Chr_",chrom,".sorted.txt",sep="")
	if(file.exists(fullfile1)){
                blueprint_func=read.table(fullfile1,head=F)
        }
        else{
                blueprint_func=read.table(fullfile2,head=F)
        }
}else{
blueprint_func=read.table(paste("/lustre/scratch115/projects/psc2020/splitfiles/",celltype_of_interest,"/",celltype_of_interest,"_",functional_assay,"_peer_10_all_summary.Chr_",chrom,".sorted.txt",sep=""),head=F)
}
#this is big, but we only care about our region around the snp we are interested in
region_start = ibd_gwas[index,4] 
region_end=ibd_gwas[index,5] #you could also get this using the name of the column, for instance here: ibd_gwas$end_position[index]

#we will use the fact that our data is sorted, to find the starting_position and ending_position
idx_start=which(as.numeric(blueprint_func[,2])>=region_start)[1]
tmp=which(as.numeric(blueprint_func[,2])<=region_end)
idx_end=tmp[length(tmp)]

#now that we have that, this is all we care about (from the whole blueprint_func data)
blueprint_func_subset = blueprint_func[idx_start:idx_end,]

#we can now even remove the big data we loaded earlier
rm(blueprint_func)
######################################################################################################################################################

#now load the data we want from the GWAS data.
ibd_gwas_data=NULL
if( tolower(trait)=="uc"){
ibd_gwas_data=read.table(paste("/lustre/scratch115/projects/psc2020/splitfiles/ibd_summary_stats/uc_b37_45975.Chr_",chrom,".ready",sep=""),head=T,sep="\t")
}else
if(tolower(trait)=="cd"){
ibd_gwas_data=read.table(paste("/lustre/scratch115/projects/psc2020/splitfiles/ibd_summary_stats/cd_b37_40266.Chr_",chrom,".ready",sep=""),head=T,sep="\t")
}
print(paste("region start: ",region_start,", and region end: ",region_end,sep=""))
#we will use the fact that our data is sorted, to find the starting_position and ending_position
idx_start_gwas=which(as.numeric(ibd_gwas_data[,2])>=region_start)[1]
tmp=which(as.numeric(ibd_gwas_data[,2])<=region_end)
idx_end_gwas=tmp[length(tmp)]

print(paste("idx start:",idx_start,", idx end: ",idx_end,sep=""))
#now that we have that, this is all we care about (from the whole ibd_gwas data)
ibd_gwas_data_subset = ibd_gwas_data[idx_start_gwas:idx_end_gwas,]

#we can now even remove the big data we loaded earlier
rm(ibd_gwas_data)
######################################################################################################################################################

varbeta_gwas=ibd_gwas_data_subset[,6]^2


#we need to split the functional file by transcript
all_transcripts=unique(as.character(blueprint_func_subset[,5]))
coloc_results=matrix(nrow=length(all_transcripts),ncol=7,NA)
coloc_results[,1]=all_transcripts
##iterate by transcript, and run coloc for each of them.
for(i in 1:length(all_transcripts)){
transcript_index=which(as.character(blueprint_func_subset[,5])==all_transcripts[i])
if(length(transcript_index)>1){ #otherwise we get errors

#addition March 1st: Sometimes there is no overlap at all in the set of SNPs (!). In that case just move on.
snp_gwas=ibd_gwas_data_subset$Marker
snp_func=blueprint_func_subset[transcript_index,4]
tmp=which(as.character(snp_func)%in%as.character(snp_gwas))
if(length(tmp)>0){
##prepare varbetas:
varbeta_func=blueprint_func_subset[transcript_index,11]^2


##this is how SunGou had run coloc:
#results <- coloc.abf(dataset1=list(N=data$N.x, s=data$s.x, pvalues=data$pvalues.x, type="cc"), dataset2=list(N=data$N.y, s=data$s.y, pvalues=data$pvalues.y, type="cc"), MAF=data$MAF, p12=1e-06) # p12 is set as what Fortune et al. 2015 recommended
#this is the one which ran fine:
results = coloc.abf(dataset1=list(snp=as.character(ibd_gwas_data_subset$Marker),beta=ibd_gwas_data_subset[,5], varbeta=varbeta_gwas,type="cc",s=case_prop, 
MAF=ibd_gwas_data_subset$controls_maf), 
dataset2=list(snp=blueprint_func_subset[transcript_index,4],beta=blueprint_func_subset[transcript_index,7],N=func_counts, 
MAF=blueprint_func_subset[transcript_index,10],varbeta=varbeta_func, type="quant"),p12=1e-06)

coloc_results[i,2:7]=results$summary
}}
}
colnames(coloc_results)=c("Gene/Feature","n_snps","H0","H1","H2","H3","H4")
write.table(coloc_results,output_filename,quote=F,col.names=T,row.names=F,sep="\t")
#now we know from them the README file from Louella that the beta is column 5 and the SE column 9.

}#here we close the for loop for the snp_of_interest

#########################################
#########################################
# So what is remaining? To get a paper first and a nobel prize down the line.



#END
