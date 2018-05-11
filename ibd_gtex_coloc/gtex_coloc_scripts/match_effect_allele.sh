#!/bin/bash

# Usage: match_effect_allele.sh <comb_file> <Effect_type> <FRQ_incl> <output_file_prefix>
# Before: comb_file contains summary statistics from two associtation studies of the same trait.(Already joined on variant name).
#	  The first columns of this file should be
#	  <SNP> <Effect_allele_study1> <Other_allele_study1> <Effect_study1> <Effect_allele_study2> <Other_allele_study2> <Effect_study2>
#	  and optionally the next two columns can be <FRQ_study1> <FRQ_study2> if <FRQ_incl> equals true.
#	  The file can include other columns but the first columns must be as specified above.
#
#	  Effect_type is either OR or BETA
#	  FRQ_incl is either true or false.
# After: ${output_file_prefix}.effect_alleleMatched contains <SNP> <Ref_allele> <Other_allele> <Effect_study1> <Effect_study2> ... <any other columns in the file>
#	 ${output_file_prefix}.badSNPs contains all variants that could not be matched between the two studies.


# Notes about the usage: (Don't depend on these facts, they might change later)
# SNPs for which the alleles don't match after accounting for strand differences are omitted.
# Triallelic SNPs that are in view_assoc data coded as !G for example are omitted.
# You must be carefull that alleles are coded in capital letters and note that if the data contains indels or deletions
# these are only handled when coded exactly the same and alleles are reported on the same strand.

# FRQ_incl states whether the frequency of the allele is supposed to be taken into account. If FRQ_incl is false then AT/GC SNPs are omitted.
# If the alleles of any AT/GC SNPs in the data have more than 25% frequency difference then we expect matching these based on frequency will be
# alright


input_file=$1
Effect_type=$2
FRQ_incl=$3
output_file_prefix=$4


if [ $Effect_type != "OR" -a $Effect_type != "or" -a $Effect_type != "BETA" -a $Effect_type != "beta" ]; then
	echo "Error: Unrecognised Effect type ${Effect_type}. Use OR or BETA"
  	echo "Quitting"
  	exit 1
fi


if [ "$FRQ_incl" != "true" -a "$FRQ_incl" != "false" ]
	then
	echo "Error: Unrecognised parameter FRQ_incl ${$FRQ_incl}. Use true or false to state if AT/GC snps should be matched based on frequency."
  	echo "Quitting"
  	exit 1
fi

head -1 ${input_file} | awk '{ out=$1" "$2" "$3" "$4" "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out }' > matced.header.tmp

#	  <SNP> <Effect_allele_study1> <Other_allele_study1> <Effect_study1> <Effect_allele_study2> <Other_allele_study2> <Effect_study2>
#	  and optionally the next two columns can be <FRQ_study1> <FRQ_study2> if <FRQ_incl> equals true

if [ \( $Effect_type == "OR" -o $Effect_type == "or" \) -a "$FRQ_incl" == "true" ]
then
	awk -v out_file="$output_file_prefix" '{if(!(($2~/[AT]/ && $3~/[AT]/) || ($2~/[GC]/ && $3~/[GC]/)) && $2==$5 && $3==$6){ out=$1" "$2" "$3" "$4" "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched" > out_file".effect_alleleMatched"} \
	else if($2=="A" && $3!="T" && $5=="T" && $6!="A" && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4" "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2=="T" && $3!="A" && $5=="A" && $6!="T" && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4" "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2=="C" && $3!="G" && $5=="G" && $6!="C" && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4" "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2=="G" && $3!="C" && $5=="C" && $6!="G" && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4" "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2~/[AT]/ && $3~/[CG]/ && $5~/[CG]/ && $6~/[AT]/ && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4" "1/$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2~/[GC]/ && $3~/[AT]/ && $5~/[AT]/ && $6~/[GC]/ && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4" "1/$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2~/[AT]/ && $3~/[AT]/ && $5~/[AT]/ && $6~/[AT]/ && length($6)==1  && length($5)==1 && length($3)==1 && length($2)==1 && $8!="." && $9!="." && ($8>=0.7 || $8<=0.3) && \
	($8-$9>= -0.25 && $8-$9<=0.25)){ out=$1" "$2" "$3" "$4" "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2~/[CG]/ && $3~/[CG]/ && $5~/[CG]/ && $6~/[CG]/  && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1 && $8!="." && $9!="." && ($8>=0.7 || $8<=0.3) && \
	($8-$9>= -0.25 && $8-$9<=0.25)){ out=$1" "$2" "$3" "$4" "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2~/[AT]/ && $3~/[AT]/ && $5~/[AT]/ && $6~/[AT]/ && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1 && $8!="." && $9!="." && ($8>=0.7 || $8<=0.3) && \
	($8-$9<= -0.25 || $8-$9>=0.25)){ out=$1" "$2" "$3" "$4" "1/$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2~/[CG]/ && $3~/[CG]/ && $5~/[CG]/ && $6~/[CG]/ && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1 && $8!="." && $9!="." && ($8>=0.7 || $8<=0.3) && \
	($8-$9<= -0.25 || $8-$9>=0.25)){ out=$1" "$2" "$3" "$4" "1/$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if((length($2)!=1 || length($3)!=1) && $2==$5 && $3==$6) { out=$1" "$2" "$3" "$4" "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if((length($2)!=1 || length($3)!=1) && $2==$6 && $3==$5) { out=$1" "$2" "$3" "$4" "1/$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else {print > out_file".badSNPs"}}' <  ${input_file}



elif [ \( $Effect_type == "BETA" -o $Effect_type == "beta" \) -a "$FRQ_incl" == "true" ]
then
	awk -v out_file="$output_file_prefix" '{if(!(($2~/[AT]/ && $3~/[AT]/) || ($2~/[GC]/ && $3~/[GC]/)) && $2==$5 && $3==$6){ out=$1" "$2" "$3" "$4"  "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2=="A" && $3!="T" && $5=="T" && $6!="A" && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4"  "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2=="T" && $3!="A" && $5=="A" && $6!="T" && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4"  "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2=="C" && $3!="G" && $5=="G" && $6!="C" && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4"  "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2=="G" && $3!="C" && $5=="C" && $6!="G" && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4"  "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2~/[AT]/ && $3~/[CG]/ && $5~/[CG]/ && $6~/[AT]/ && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4" "$7*(-1); for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2~/[GC]/ && $3~/[AT]/ && $5~/[AT]/ && $6~/[GC]/ && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4" "$7*(-1); for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2~/[AT]/ && $3~/[AT]/ && $5~/[AT]/ && $6~/[AT]/ && length($6)==1  && length($5)==1 && length($3)==1 && length($2)==1 && $8!="." && $9!="." && ($8>=0.7 || $8<=0.3) && \
	($8-$9>= -0.25 && $8-$9<=0.25)){ out=$1" "$2" "$3" "$4"  "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2~/[CG]/ && $3~/[CG]/ && $5~/[CG]/ && $6~/[CG]/  && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1 && $8!="." && $9!="." && ($8>=0.7 || $8<=0.3) && \
	($8-$9>= -0.25 && $8-$9<=0.25)){ out=$1" "$2" "$3" "$4"  "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2~/[AT]/ && $3~/[AT]/ && $5~/[AT]/ && $6~/[AT]/ && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1 && $8!="." && $9!="." && ($8>=0.7 || $8<=0.3) && \
	($8-$9<= -0.25 || $8-$9>=0.25)){ out=$1" "$2" "$3" "$4"  "$7*(-1); for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2~/[CG]/ && $3~/[CG]/ && $5~/[CG]/ && $6~/[CG]/ && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1 && $8!="." && $9!="." && ($8>=0.7 || $8<=0.3) && \
	($8-$9<= -0.25 || $8-$9>=0.25)){ out=$1" "$2" "$3" "$4"  "$7*(-1); for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if((length($2)!=1 || length($3)!=1) && $2==$5 && $3==$6) { out=$1" "$2" "$3" "$4"  "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if((length($2)!=1 || length($3)!=1) && $2==$6 && $3==$5) { out=$1" "$2" "$3" "$4" "$7*(-1); for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else {print > out_file".badSNPs"}}' <  ${input_file}


elif [ \( $Effect_type == "BETA" -o $Effect_type == "beta" \) -a "$FRQ_incl" != "true" ]
then
	echo "Frequency not included. All AT/GC SNPs will be discarded"
	echo "To match AT/GC snps by frequency, set third parameter as true and include frequency of the SNPs in columns 8 and 9 of the input file"


	awk -v out_file="$output_file_prefix" '{if(!(($2~/[AT]/ && $3~/[AT]/) || ($2~/[GC]/ && $3~/[GC]/)) && $2==$5 && $3==$6){ out=$1" "$2" "$3" "$4" "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2=="A" && $3!="T" && $5=="T" && $6!="A" && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4" "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2=="T" && $3!="A" && $5=="A" && $6!="T" && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4" "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2=="C" && $3!="G" && $5=="G" && $6!="C" && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4" "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2=="G" && $3!="C" && $5=="C" && $6!="G" && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4" "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2~/[AT]/ && $3~/[CG]/ && $5~/[CG]/ && $6~/[AT]/ && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4" "$7*(-1); for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2~/[GC]/ && $3~/[AT]/ && $5~/[AT]/ && $6~/[GC]/ && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4" "$7*(-1); for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if((length($2)!=1 || length($3)!=1) && $2==$5 && $3==$6) { out=$1" "$2" "$3" "$4" "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if((length($2)!=1 || length($3)!=1) && $2==$6 && $3==$5) { out=$1" "$2" "$3" "$4" "$7*(-1); for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else {print > out_file".badSNPs"}}' <  ${input_file}


else
	echo "Frequency not included. All AT/GC SNPs will be discarded"
	echo "To match AT/GC snps by frequency, set third parameter as true and include frequency of the SNPs in columns 8 and 9 of the input file"

	awk  -v out_file="$output_file_prefix" '{if(!(($2~/[AT]/ && $3~/[AT]/) || ($2~/[GC]/ && $3~/[GC]/)) && $2==$5 && $3==$6) { out=$1" "$2" "$3" "$4" "$7; for(i=8;i<=NF;i++) {out=out" "$i}; print out > out_file".effect_alleleMatched"}
else if($2=="A" && $3!="T" && $5=="T" && $6!="A" && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4" "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2=="T" && $3!="A" && $5=="A" && $6!="T" && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4" "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2=="C" && $3!="G" && $5=="G" && $6!="C" && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4" "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2=="G" && $3!="C" && $5=="C" && $6!="G" && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4" "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2~/[AT]/ && $3~/[CG]/ && $5~/[CG]/ && $6~/[AT]/ && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4" "1/$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if($2~/[GC]/ && $3~/[AT]/ && $5~/[AT]/ && $6~/[GC]/ && length($6)==1 && length($5)==1 && length($3)==1 && length($2)==1){ out=$1" "$2" "$3" "$4" "1/$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if((length($2)!=1 || length($3)!=1) && $2==$5 && $3==$6) { out=$1" "$2" "$3" "$4" "$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else if((length($2)!=1 || length($3)!=1) && $2==$6 && $3==$5) { out=$1" "$2" "$3" "$4" "1/$7; for(i=8;i<=NF;i++){out=out" "$i}; print out > out_file".effect_alleleMatched"} \
	else {print > out_file".badSNPs"}}' <  ${input_file}

fi


cat matced.header.tmp ${output_file_prefix}.effect_alleleMatched > ${output_file_prefix}.effect_alleleMatched.tmp
mv ${output_file_prefix}.effect_alleleMatched.tmp ${output_file_prefix}.effect_alleleMatched
#rm matced.header.tmp
