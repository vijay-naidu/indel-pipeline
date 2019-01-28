#!/bin/bash

names='ERR502366'
first='_1_val_1'
second='_2_val_2'
for name in $names
do
echo $name
bwa mem AL123456.fa $name$first'.fq' $name$second'.fq' > $name'.sam'
samtools view -bS $name'.sam' > $name'.bam'
samtools sort $name'.bam' -o $name'.sorted'
samtools index $name'.sorted'
java -jar pilon-1.22.jar --genome AL123456.fa --frags $name'.sorted' --output $name --vcf
vcftools --vcf $name'.vcf' --keep-only-indels --out 'indels'$name --recode --recode-INFO-all
vcftools --vcf $name'.vcf' --out 'snps'$name --recode --recode-INFO-all --non-ref-ac-any 1
done
echo All done
