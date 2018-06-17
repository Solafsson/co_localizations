#!/bin/bash

## Usage: This is a wrapper script for the processing of GTEX data that needs to be done before I can do colocalization
##        analysis of IBD SNPs and epithelial eQTLs. Done as part of my rotation with Carl.
##        processGTexData.sh <GtexTrait_list> <diseaseList>

## Input: GtexTrait_list is a two column text file without a header that gives the name of a GTex file (found under
##        the gtex_path) in the first column, and the sample size in the second column.
##        diseaseList is a two column text file without a header that gives the name of a disease in the first column
##        and the case-fraction in the second column. It's important that the name of the disease matches the prefix
##        of the summary statistics files.
##        Expect the summary statistics files to match the regular expression: <disease>_b37*Chr_${CHR}.ready

##        example: GtexTrait_list:
#           Colon_Sigmoid.allpairs.txt.gz 203
#           Colon_Transverse.allpairs.txt.gz 246
#         example: diseaseList:
#            uc 0.354
#            cd 0.349

GtexTrait_list=$1
diseaseList=$2

gtex_path=/lustre/scratch115/resources/GTEx/AnalysisV7/GTEx_Analysis_v7_eQTL_all_associations/
ibd_assoc_path=/lustre/scratch115/projects/coloc_and_fm/splitfiles/ibd_summary_stats/
#working_dir=/lustre/scratch113/teams/anderson/users/so11/gtex_coloc/
working_dir=/lustre/scratch119/humgen/teams/anderson/users/so11/co_localizations/ibd_gtex_coloc/testing_forLiz/
script_dir=/nfs/users/nfs_s/so11/phd/co_localizations/ibd_gtex_coloc/gtex_coloc_scripts/


mkdir -p ${working_dir}cluster_output_files
clusterErrorOutput=${working_dir}cluster_output_files/
nrTissues=$( wc -l ${GtexTrait_list} | awk '{print $1}' )

bsub  -J"copy_and_unzip[1-${nrTissues}]" -W 60 -M100 -R'span[hosts=1] select[mem>100] rusage[mem=100]' \
-e ${clusterErrorOutput}copy_and_unzip_errors.%J.%I -o ${clusterErrorOutput}copy_and_unzip_output.%J.%I \
bash ${script_dir}processGTexData_wrapper.sh ${script_dir}processGTexData_worker.sh copy_and_unzip ${GtexTrait_list} \
${gtex_path} ${ibd_assoc_path} ${working_dir} ${GtexTrait_list}



# Remove .gz from tissue name
sed 's/\.gz//g' ${GtexTrait_list} > ${working_dir}traitList2.txt


## Divide files up by chromosomes
bsub -w 'done("copy_and_unzip")' -J"break_into_chromosomes[1-${nrTissues}]" -W 240 -M1000 -R'span[hosts=1] select[mem>1000] rusage[mem=1000]' \
-e ${clusterErrorOutput}break_into_chromosomes_errors.%J.%I -o ${clusterErrorOutput}break_into_chromosomes_output.%J.%I \
bash ${script_dir}processGTexData_wrapper.sh ${script_dir}processGTexData_worker.sh break_into_chromosomes ${working_dir}traitList2.txt \
${gtex_path} ${ibd_assoc_path} ${working_dir} ${working_dir}traitList2.txt




## Join files based on position
## Use decode script to match effect alleles
rm -f ${output_dir}join_match_jobFile.txt
while read tissue sampleSize; do

    tissueName=$( echo ${tissue} | cut -f 1 -d "." )
    for i in {1..22}; do
        while read disease caseFrac; do
            echo "${disease} ${i} ${tissueName} ${sampleSize}" >> ${working_dir}join_match_jobFile.txt
        done < ${diseaseList}
    done

    # No IBD associations on the X chromosome?
    #echo "cd X ${tissueName}" >> ${working_dir}join_match_jobFile.txt
    #echo "uc X ${tissueName}" >> ${working_dir}join_match_jobFile.txt
done < ${GtexTrait_list}

nrJoinMatchJobs=$( wc -l ${working_dir}join_match_jobFile.txt | awk '{print $1}' )


bsub -w 'done("break_into_chromosomes")' -J"join_and_match[1-${nrJoinMatchJobs}]" -M2000 -R'span[hosts=1] select[mem>2000] rusage[mem=2000]' \
-e ${clusterErrorOutput}join_and_match_errors.%J.%I -o ${clusterErrorOutput}join_and_match_output.%J.%I \
bash ${script_dir}processGTexData_wrapper.sh ${script_dir}processGTexData_worker.sh join_and_match ${working_dir}traitList2.txt \
${gtex_path} ${ibd_assoc_path} ${working_dir} ${working_dir}join_match_jobFile.txt


rm -f ${working_dir}coloc_jobFile.txt
while read tissue sampleSize; do
    tissueName=$( echo ${tissue} | cut -f 1 -d "." )
    while read disease caseFrac; do
        echo "${tissueName} ${sampleSize} ${disease} ${caseFrac}" >> ${working_dir}coloc_jobFile.txt
    done < ${diseaseList}
done < ${GtexTrait_list}
nrColocJobs=$( wc -l ${working_dir}coloc_jobFile.txt | awk '{print $1}' )

## Run the colocalization
bsub -w 'done("join_and_match")' -J"perform_coloc[1-${nrColocJobs}]" -W 240 -M5000 -R'span[hosts=1] select[mem>5000] rusage[mem=5000]' \
-e ${clusterErrorOutput}perform_coloc_errors.%J.%I -o ${clusterErrorOutput}perform_coloc_output.%J.%I \
bash ${script_dir}processGTexData_wrapper.sh ${script_dir}processGTexData_worker.sh coloc ${working_dir}traitList2.txt \
${gtex_path} ${ibd_assoc_path} ${working_dir} ${working_dir}coloc_jobFile.txt






