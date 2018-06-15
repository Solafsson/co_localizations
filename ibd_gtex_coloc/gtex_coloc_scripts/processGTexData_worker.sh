#!/bin/bash

run_function=$1
traitList=$2
gtex_path=$3
ibd_assoc_path=$4
working_dir=$5
jobFile=$6
jobIndex=$7


indexTissue=$( awk -v jI="$jobIndex" 'NR==jI {print $1; exit}' < ${traitList} )
indexTissue_N=$( awk -v jI="$jobIndex" 'NR==jI {print $2; exit}' < ${traitList} )
tissueName=$(echo ${indexTissue} | cut -f 1 -d "." )

## Usage: Pretty self explanatory.
## Input: tissue is the name of a file found under the gtex_path (e.g Colon_Sigmoid.allpairs.txt.gz)
function copy_and_unzip {

    tissue=$1
    working_dir=$2
    gtex_path=$3
    tissueName=$4

    mkdir -p ${working_dir}${tissueName}
    cp ${gtex_path}${tissue} ${working_dir}${tissueName}
    gunzip ${working_dir}${tissueName}/${tissue}
}

## Usage: The purpose of this function is to split GTex summary statistics (~120M lines) into more manageable chromosome-specific files
function break_into_chromosomes {

    gtex_file=$1
    output_dir=$2
    tissueName=$3
    sample_size=$4
    fileName=$(echo ${gtex_file} | awk 'BEGIN {FS="/"} {print $NF}' | cut -f 1 -d "." )

    awk -v N="$sample_size" '{print $2"_"N}' < ${gtex_file} | awk 'BEGIN {FS="_"; print "Chr bp A1 A2 sample_size"} NR>1 {print $1, $2, $3, $4, $5}' | paste - ${gtex_file} | \
    awk -v out="$output_dir" -v file="$fileName" '{for(i=1; i<=22; i++) \
    {if ($1==i) {print > out "chr" i "_" file}} {if($1=="X") {print > out "chrX_" file}}}'
}


# Usage: To get mafs and rsIDs for the variants in the GTex summary data, I join
#
# Output: Give an output on a format that is compatible with the match_effect_allele script written in decode.
function join_with_ibd_sumstats {

    ibd_assoc_path=$1
    working_dir=$2
    jobFile=$3
    jobIndex=$4

    ibd_subtype=$( awk -v jI="$jobIndex" 'NR==jI {print $1; exit}' < ${jobFile} )
    CHR=$( awk -v jI="$jobIndex" 'NR==jI {print $2; exit}' < ${jobFile} )
    tissue=$( awk -v jI="$jobIndex" 'NR==jI {print $3; exit}' < ${jobFile} | cut -f 1 -d "." )

    by_chr_gtex_path=${working_dir}${tissue}/

    join -1 2 -2 2 <(sort -k2,2 ${ibd_assoc_path}${ibd_subtype}_*Chr_${CHR}* ) \
    <( sort -k2,2 ${by_chr_gtex_path}chr${CHR}_${tissue} ) | \
    awk 'BEGIN {print "SNP ibd_eff_All ibd_other_All ibd_beta expr_eff_All expr_other_All expr_beta \
    ibd_maf expr_maf CHR BP ibd_se expr_se ibd_P expr_P ibd_N expr_N geneID"} {print $3, toupper($8), toupper($7), $5, \
    $13, $12,$22, $10, $20, $2, $1, $6, $23, $4, $21, $9, $14, $15}' > ${by_chr_gtex_path}${ibd_subtype}_chr${CHR}_${tissue}_formatching.txt

}

function match_effect_allele {

    ibd_assoc_path=$1
    working_dir=$2
    jobFile=$3
    jobIndex=$4
    effectType=$5
    frq_incl=$6

    ibd_subtype=$( awk -v jI="$jobIndex" 'NR==jI {print $1; exit}' < ${jobFile} )
    CHR=$( awk -v jI="$jobIndex" 'NR==jI {print $2; exit}' < ${jobFile} )
    tissue=$( awk -v jI="$jobIndex" 'NR==jI {print $3; exit}' < ${jobFile} | cut -f 1 -d "." )
    sampleSize=$( awk -v jI="$jobIndex" 'NR==jI {print $4; exit}' < ${jobFile} | cut -f 1 -d "." )

    by_chr_gtex_path=${working_dir}${tissue}/

    /nfs/users/nfs_s/so11/rotations/andersonLab/scripts/gtex_coloc_scripts/match_effect_allele.sh ${by_chr_gtex_path}${ibd_subtype}_chr${CHR}_${tissue}_formatching.txt \
    ${effectType} ${frq_incl} ${by_chr_gtex_path}${ibd_subtype}_chr${CHR}_${tissue}

    awk '{print $13, $15, $7, $5, $11, $1, $16, $8, $9}' < ${by_chr_gtex_path}${ibd_subtype}_chr${CHR}_${tissue}.effect_alleleMatched > \
    ${by_chr_gtex_path}${ibd_subtype}_chr${CHR}_${tissue}_forColoc.txt


    sed -i "s/b37/${sampleSize}/g" ${by_chr_gtex_path}${ibd_subtype}_chr${CHR}_${tissue}_forColoc.txt

}

function coloc {

    jobFile=$1
    jobIndex=$2

    ibd_subtype=$( awk -v jI="$jobIndex" 'NR==jI {print $3; exit}' < ${jobFile} )
    tissue=$( awk -v jI="$jobIndex" 'NR==jI {print $1; exit}' < ${jobFile} | cut -f 1 -d "." )


    /software/R-3.3.0/bin/Rscript /nfs/users/nfs_s/so11/rotations/andersonLab/scripts/gtex_coloc_scripts/performColoc.R ${ibd_subtype} ${tissue}

}
if [ "$run_function" == "copy_and_unzip" ]; then
    copy_and_unzip ${indexTissue} ${working_dir} ${gtex_path} ${tissueName}
fi

if [ "$run_function" == "break_into_chromosomes" ];  then
    break_into_chromosomes ${working_dir}${tissueName}/${tissueName}.allpairs.txt ${working_dir}${tissueName}/ ${tissueName} ${indexTissue_N}
fi

## Combine to one job joining and matching
if [ "$run_function" == "join_and_match" ];  then
    join_with_ibd_sumstats ${ibd_assoc_path} ${working_dir} ${working_dir}join_match_jobFile.txt ${jobIndex}
    match_effect_allele ${ibd_assoc_path} ${working_dir} ${working_dir}join_match_jobFile.txt ${jobIndex} BETA true
fi

if [ "$run_function" == "coloc" ]; then
    coloc ${working_dir}coloc_jobFile.txt ${jobIndex}

fi


exit $?


