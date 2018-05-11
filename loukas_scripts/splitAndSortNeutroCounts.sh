#!/bin/bash
t1dcc="/lustre/scratch115/projects/psc2020/splitfiles/neutrophil_counts_summary_stats/neut_N170702_narrow_form.tsv"

outputcc="/lustre/scratch115/projects/psc2020/splitfiles/neutrophil_counts_summary_stats/neut_N170702_narrow"
for((i=1;i<=22;i++));
do
# If the results are very different go back to it.
cat ${t1dcc} | awk -v chr=$i '{split($8,arr,">");if($3==chr){ print $3"\t"$4"\t"$2"\t"$11"\t"exp($9)"\t"$5"\t"$6"\t"170702"\t"$14}}' > ${outputcc}.Chr_${i}.tmp
	cat ${outputcc}.Chr_${i}.tmp |  sort -nk2 >  ${outputcc}.Chr_${i}.sorted.txt
        sed -i '1iChr\tPosition\tMarker\tPvalue\tOR(AltAllele)\tRefAllele\tAltAllele\tSampleSize\tFreqAltAllele' ${outputcc}.Chr_${i}.sorted.txt #add header
done
rm ${outputcc}*tmp

#I also need to add frequencies for the minor allele...
