#!/bin/bash
set -x
set -e
set -o pipefail

###### PCA analysis #####
plink2 --vcf /disk3/call_snp/plink/vcf71.${min_dp}.${max_dp}.maf.0.05.ld.vcf --aec --pca --out ./pca

###### Admixture analysis for K=1 to 5 #####
plink2 --vcf /disk3/call_snp/plink/vcf71.${min_dp}.${max_dp}.maf.0.05.ld.vcf  --aec 0 --make-bed -out /disk3/call_snp/plink/vcf71.${min_dp}.${max_dp}.maf.0.05.ld.str
#
for k in {1..5}
do
  /home/share/admixture_linux-1.3.0/admixture --cv /disk3/call_snp/plink/vcf71.${min_dp}.${max_dp}.maf.0.05.ld.str.bed $k -j10 | tee ./admix.log${k}.out &
done
wait

###### Extract cross-validation (CV) errors #####
grep -h CV ./admix.log*.out > ./admix.CV.txt
