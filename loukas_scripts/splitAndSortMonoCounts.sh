#!/bin/bash
t1dcc="/lustre/scratch115/projects/psc2020/splitfiles/monocyte_counts_summary_stats/mono_N170721_narrow_form.tsv"

outputcc="/lustre/scratch115/projects/psc2020/splitfiles/monocyte_counts_summary_stats/mono_N170721_narrow"
for((i=1;i<=22;i++));
do
# I've put the sample sizes based on the 2015 Nat Genet paper but it was not clear what they are exactly.
# If the results are very different go back to it.
cat ${t1dcc} | awk -v chr=$i '{split($8,arr,">");if($3==chr){ print $3"\t"$4"\t"$2"\t"$11"\t"exp($9)"\t"$5"\t"$6"\t"170721"\t"$14}}' > ${outputcc}.Chr_${i}.tmp
	cat ${outputcc}.Chr_${i}.tmp |  sort -nk2 >  ${outputcc}.Chr_${i}.sorted.txt
        sed -i '1iChr\tPosition\tMarker\tPvalue\tOR(AltAllele)\tRefAllele\tAltAllele\tSampleSize\tFreqAltAllele' ${outputcc}.Chr_${i}.sorted.txt #add header
done
rm ${outputcc}*tmp

#I also need to add frequencies for the minor allele...
