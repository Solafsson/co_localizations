#!/bin/bash

## Note: The geneKeyName has to be created manually. You would run extract_coloc_sumstats, look at your genes in the
##       all_coloc_hits.txt file and create this list by googling the gene IDs before running this script.

results_dir=/lustre/scratch119/humgen/teams/anderson/users/so11/co_localizations/ibd_gtex_coloc/results_long_window/
geneList=/lustre/scratch119/humgen/teams/anderson/users/so11/co_localizations/ibd_gtex_coloc/for_locusZoom_long_window/all_coloc_hits.txt
geneKeyName=/lustre/scratch119/humgen/teams/anderson/users/so11/co_localizations/ibd_gtex_coloc/for_locusZoom_long_window/geneNameKey.txt


## Extract the posterior probabilities for each gene
rm -rf ${results_dir}posterior_probs/
mkdir -p ${results_dir}posterior_probs
while read gene; do
    rm -f ${results_dir}posterior_probs/${gene}_postProb.txt
    for Tiss in Colon_Sigmoid Colon_Transverse Small_Intestine_Terminal_Ileum; do
        grep ${gene} ${results_dir}cd_${Tiss}_colocalization.txt | awk -v T="$Tiss" '{print $0, T}' >> ${results_dir}posterior_probs/${gene}_postProb.txt
        grep ${gene} ${results_dir}uc_${Tiss}_colocalization.txt | awk -v T="$Tiss" '{print $0, T}' >> ${results_dir}posterior_probs/${gene}_postProb.txt
    done
done < <(awk '{print $2}' < ${geneList} | sort -u )

while read id gene; do
    rename "s/${id}/${gene}/" ${results_dir}posterior_probs/*
done < ${geneKeyName}