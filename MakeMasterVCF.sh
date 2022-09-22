#!/bin/bash
#go to where all the good VCFs are 
cd /net/dunham/vol2/Renee/yEvo_seq/220526_RB_caffeine/WorkDirectory/

#grab all the RBcaff folder names and store it into 
for directory in ./*/; do basename $directory; done | awk '$0 ~ "RBcaff" {print $0}' > /net/dunham/vol2/Leah/hackathon/MasterVCF/RBCaff.txt

File="/net/dunham/vol2/Leah/hackathon/MasterVCF/RBCaff.txt"
Lines=$(cat $File)
for Line in $Lines
do
	cd /net/dunham/vol2/Renee/yEvo_seq/220526_RB_caffeine/WorkDirectory/$Line
	cp *samtools_filtered.vcf.annotated /net/dunham/vol2/Leah/hackathon/MasterVCF/
done

cd /net/dunham/vol2/Leah/hackathon/MasterVCF/

for Line in $Lines
do
	awk -vOFS="\t" '$0 !~ "##" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$11,$12,$13,$14="SAMPLE"}' ${Line}_samtools_filtered.vcf.annotated | awk -v awkvar=$Line 'NR>1 {$13=awkvar}1' OFS='\t' > ${Line}_VCFclean.txt 
done

awk FNR!=1 *_VCFclean.txt > temp.txt

awk -vOFS="\t" '$0 !~ "##" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$11,$12,$13,$14="SAMPLE"}' RBcaff_1blue1_S9_samtools_filtered.vcf.annotated | head -n 1 > header.txt 

sed -i 's/#//' header.txt

cat header.txt temp.txt > MasterVCF.txt
 



