#!/bin/bash

## Usage: This script takes in a directory of co-localization results. For each SNP with H4>0.8, it fetches the
##        summary statistics for both the ibd_subtype and the gene and writes them to the output_dir

## Input: coloc_dir should be a directory containing files ending with *_colocalization.txt
#         Should have this format:
#Trait1	Trait2	Top_SNP_Trait1	n_snps	H0	H1	H2	H3	H4
#cd	ENSG00000008128.18	rs12103	454	0.628091413876393	0.329490264733835	0.0256528916955941	0.0134238258361399	0.00334160385803789
#cd	ENSG00000008130.11	rs12103	454	0.628617561336236	0.329766276253818	0.0256890013223619	0.0134514272736094	0.00247573381397493
#cd	ENSG00000067606.11	rs12103	454	0.628947852039626	0.329939543343441	0.0252666161427951	0.0132284301282121	0.00261755834592539

coloc_dir=/lustre/scratch113/teams/anderson/users/so11/gtex_coloc/results/
ibd_gwas=/lustre/scratch115/projects/coloc_and_fm/IBD_regions_from_our_2017_paper.txt
output_dir=/lustre/scratch113/teams/anderson/users/so11/gtex_coloc/for_locusZoom/
working_dir=/lustre/scratch113/teams/anderson/users/so11/gtex_coloc/
hard_H4_threshold=0.8
soft_H4_threshold=0.5
soft_H4_H4plusH3_threshold=0.9

ls ${coloc_dir}*_colocalization.txt > ${output_dir}coloc_fileList.tmp


## Gather all the results in one file
rm -f ${output_dir}hard_thres_${hard_H4_threshold}_coloc_hits.txt
while read colocFile; do
    ibd_subtype=$( echo ${colocFile} | awk 'BEGIN {FS="/"} {print $NF}' | cut -f 1 -d "_" )
    tissue=$( echo ${colocFile} | awk 'BEGIN {FS="/"} {print substr($NF, 4, length($NF)-22)}' )
    awk -v H4="$hard_H4_threshold" -v Tissue="$tissue" '$9>H4 && $9~/^[0-9][0-9]*/ {print $1, $2, $3, Tissue}' < \
    ${colocFile} >> ${output_dir}hard_thres_${hard_H4_threshold}_coloc_hits.txt

done < ${output_dir}coloc_fileList.tmp


rm -f ${output_dir}soft_thres_${soft_H4_threshold}_${soft_H4_H4plusH3_threshold}_coloc_hits.txt
while read colocFile; do
    ibd_subtype=$( echo ${colocFile} | awk 'BEGIN {FS="/"} {print $NF}' | cut -f 1 -d "_" )
    tissue=$( echo ${colocFile} | awk 'BEGIN {FS="/"} {print substr($NF, 4, length($NF)-22)}' )
    awk -v H4="$soft_H4_threshold" -v ratio="$soft_H4_H4plusH3_threshold" -v Tissue="$tissue" '$9>H4 && $9~/^[0-9][0-9]*/ && ($9/($9+$8))>ratio {print $1, $2, $3, Tissue}' < \
    ${colocFile} >> ${output_dir}soft_thres_${soft_H4_threshold}_${soft_H4_H4plusH3_threshold}_coloc_hits.txt

done < ${output_dir}coloc_fileList.tmp



## For comparison purposes we want to extract the data for the cartesian product of ibd_subtype x tissue x gene?
## Loop over all the hits and get the summary statistics
while read ibd_subtype gene topSNP tissue; do

    CHR=$( grep -w -m 1 ${topSNP} ${ibd_gwas} | awk 'BEGIN {FS="\t"} {print $1}' )
    LD_left=$( grep -w -m 1 ${topSNP} ${ibd_gwas} | awk 'BEGIN {FS="\t"} {print $4}' )
    LD_right=$( grep -w -m 1 ${topSNP} ${ibd_gwas} | awk 'BEGIN {FS="\t"} {print $5}' )

    for dis in uc cd; do
        for Tiss in Colon_Sigmoid Colon_Transverse Small_Intestine_Terminal_Ileum; do

             if [ ! -s ${output_dir}${dis}_chr${CHR}_${Tiss}_${gene}_${topSNP}.txt ]; then
                awk -v GENE="$gene" 'NR==1 || $7==GENE {print}' < ${working_dir}${Tiss}/${dis}_chr${CHR}_${Tiss}_forColoc.txt > \
                ${output_dir}${dis}_chr${CHR}_${Tiss}_${gene}_${topSNP}.txt
             fi

             if [ ! -s ${output_dir}${dis}_chr${CHR}_${Tiss}_${dis}_${topSNP}.txt ]; then
                 awk -v LD_r="$LD_right" -v LD_l="$LD_left" 'NR==1 || ($9>=LD_l && $9<=LD_r) {print $1, $12, $8, $9}' < \
                 ${working_dir}${Tiss}/${dis}_chr${CHR}_${Tiss}.effect_alleleMatched | sort -u > \
                 ${output_dir}${dis}_chr${CHR}_${Tiss}_${dis}_${topSNP}.txt
             fi
        done
    done


done < ${output_dir}all_coloc_hits.txt



