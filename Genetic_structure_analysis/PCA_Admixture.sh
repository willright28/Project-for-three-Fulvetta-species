#!/bin/bash
set -x
set -e
set -o pipefail

###### PCA analysis #####
plink2 --vcf ./maf.0.05.ld.vcf --aec --pca --out ./pca

###### Admixture analysis for K=1 to 5 #####
for k in {1..5}
do
  /home/share/admixture_linux-1.3.0/admixture --cv ./str.bed $k -j10 | tee ./admix.log${k}.out &
done

wait

###### Extract cross-validation (CV) errors #####
grep -h CV ./admix.log*.out > ./admix.CV.txt
