#!/bin/bash
set -e
set -x 
set -o pipefail

##################################

min_dp=8
max_dp=60

#monomorphic site
plink2 --vcf /disk3/call_snp/bcftools.raw.out.vcf.gz -aec --geno 0.1 --vcf-min-dp  8 --vcf-max-dp 60 --min-alleles 1 --max-alleles 1 --var-min-qual 30  --export vcf-4.2  --out /disk3/call_snp/plink/hkqm.${min_dp}.${max_dp}.mono
#polymorphic site
plink2 --vcf /disk3/call_snp/bcftools.raw.out.vcf.gz -aec --geno 0.1 --vcf-min-dp  8 --vcf-max-dp 60 --min-alleles 2  --max-alleles 2 --var-min-qual 30 --export vcf-4.2  --out /disk3/call_snp/plink/hkqm.${min_dp}.${max_dp}

bcftools filter /disk3/call_snp/plink/hkqm.${min_dp}.${max_dp}.vcf --SnpGap 5 --threads 70 -Oz -o /disk3/call_snp/plink/hkqm.${min_dp}.${max_dp}.gap.vcf.gz

#remove indels
plink2 --vcf /disk3/call_snp/plink/hkqm.${min_dp}.${max_dp}.gap.vcf.gz --allow-extra-chr --snps-only --export vcf-4.2  --out /disk3/call_snp/plink/hkqm.${min_dp}.${max_dp}.gap.snponly
plink2  --vcf /disk3/call_snp/plink/hkqm.${min_dp}.${max_dp}.gap.snponly.vcf --aec --maf 0.001 --keep /disk3/call_snp/71.sample.txt  --set-all-var-ids @:#  --rm-dup force-first   --export vcf-4.2  --out  /disk3/call_snp/plink/vcf71.${min_dp}.${max_dp}

####################
#filtering for population strucutre analysis
#maf0.05
plink2 --vcf /disk3/call_snp/plink/vcf71.${min_dp}.${max_dp}.vcf --aec --maf 0.05  --export vcf-4.2 --out /disk3/call_snp/plink/vcf71.${min_dp}.${max_dp}.maf.0.05

#LD
plink2 --vcf /disk3/call_snp/plink/vcf71.${min_dp}.${max_dp}.maf.0.05.vcf --allow-extra-chr --set-all-var-ids @:#  --indep-pairwise 50 10 0.1  --out /disk3/call_snp/plink/vcf71.8.60.maf.0.05.ld
plink2 --vcf /disk3/call_snp/plink/vcf71.${min_dp}.${max_dp}.maf.0.05.vcf --allow-extra-chr --set-all-var-ids @:# --rm-dup force-first --extract /disk3/call_snp/plink/vcf71.8.60.maf.0.05.ld.prune.in --recode vcf-4.2 --out /disk3/call_snp/plink/vcf71.8.60.maf.0.05.ld

#####################
#filtering for genotype-environment association analysis
#no missing
plink2 --vcf /disk3/call_snp/plink/vcf71.${min_dp}.${max_dp}.vcf --aec --geno 0  --export vcf-4.2  --out /disk3/call_snp/plink/vcf71.nomiss.${min_dp}.${max_dp}

#maf 0.05
plink2 --vcf /disk3/call_snp/plink/vcf71.nomiss.${min_dp}.${max_dp}.vcf --aec --maf 0.05 --export vcf-4.2  --out /disk3/call_snp/plink/vcf71.nomiss.${min_dp}.${max_dp}.maf.0.05
bgzip /disk3/call_snp/plink/vcf71.nomiss.${min_dp}.${max_dp}.maf.0.05.vcf
tabix /disk3/call_snp/plink/vcf71.nomiss.${min_dp}.${max_dp}.maf.0.05.vcf.gz

###include outgroup
plink2  --vcf /disk3/call_snp/plink/hkqm.${min_dp}.${max_dp}.gap.snponly.vcf --aec --maf 0.001 --keep /disk3/call_snp/81.sample.txt  --set-all-var-ids @:#  --rm-dup force-first --export vcf-4.2  --out /disk3/call_snp/plink/vcf81.${min_dp}.${max_dp}
plink2 --vcf  /disk3/call_snp/plink/vcf81.${min_dp}.${max_dp}.vcf --allow-extra-chr --set-all-var-ids @:#  --indep-pairwise 50 10 0.1  --out /disk3/call_snp/plink/vcf81.${min_dp}.${max_dp}.ld
plink2 --vcf  /disk3/call_snp/plink/vcf81.${min_dp}.${max_dp}.vcf --allow-extra-chr --set-all-var-ids @:# --rm-dup force-first --extract /disk3/call_snp/plink/vcf81.${min_dp}.${max_dp}.ld.prune.in --recode vcf-4.2 --out /disk3/call_snp/plink/vcf81.${min_dp}.${max_dp}.ld

#maf 0.05 
#plink2 --vcf  /disk3/call_snp/plink/vcf81.${min_dp}.${max_dp}.vcf --aec --maf 0.05 --export vcf-4.2  --out/disk3/call_snp/plink/vcf81.${min_dp}.${max_dp}.maf.0.05
plink2 --vcf  /disk3/call_snp/plink/vcf81.8.60.maf.0.05.vcf --allow-extra-chr --set-all-var-ids @:#  --indep-pairwise 50 10 0.1  --out /disk3/call_snp/plink/vcf81.8.60.maf.0.05.ld
plink2 --vcf /disk3/call_snp/plink/vcf81.8.60.maf.0.05.vcf --allow-extra-chr --set-all-var-ids @:# --rm-dup force-first --extract /disk3/call_snp/plink/vcf81.8.60.maf.0.05.ld.prune.in --recode vcf-4.2 --out /disk3/call_snp/plink/vcf81.8.60.maf.0.05.ld





