#!/bin/bash
set -x
set -e
set -o pipefail

### Generate a pairwise genetic distance matrix using VCF2Dis ###
/home/share/VCF2Dis-1.42/bin/VCF2Dis \
  -InPut /disk3/call_snp/plink/vcf71.${min_dp}.${max_dp}.maf.0.05.ld.vcf \
  -OutPut ./vcf2dis.dist

### Construct the neighbor-joining (NJ) tree using the FastME web server ###
# Upload the output file `vcf2dis.dist` to the FastME online platform:
# http://www.atgc-montpellier.fr/fastme/
# Choose the BioNJ algorithm

