#!/bin/bash

python VCF_processing/parseVCF.py -i /mnt/disk1/martin/vcf81.8.60.maf.0.05.vcf -o /mnt/disk1/martin//output.geno.gz

python3 /mnt/disk1/martin/genomics_general-master/ABBABABAwindows.py \
-g /mnt/disk1/martin/output.geno.gz -f diplo \
-o /mnt/disk1/martin/1k_set1.csv \
-P1 Hainan -P2 East -P3 Central -O htsm \
--popsFile /mnt/disk1/martin/species_sets.txt -w 1000 -m 10 --T 80

python3 /mnt/disk1/martin/genomics_general-master/ABBABABAwindows.py \
-g /mnt/disk1/martin/output.geno.gz -f diplo  \
-o  /mnt/disk1/martin/1k_set2.csv \
-P1 East -P2 Central -P3 West -O htsm \
--popsFile /mnt/disk1/martin/species_sets.txt -w 1000 -m 10 --T 80
