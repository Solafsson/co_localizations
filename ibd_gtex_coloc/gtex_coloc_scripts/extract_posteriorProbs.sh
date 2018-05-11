#!/bin/bash



results_dir=/lustre/scratch113/teams/anderson/users/so11/gtex_coloc/results/
geneList=/lustre/scratch113/teams/anderson/users/so11/gtex_coloc/for_locusZoom/all_coloc_hits.txt
geneKeyName=/lustre/scratch113/teams/anderson/users/so11/gtex_coloc/for_locusZoom_backup/geneNameKey.txt


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