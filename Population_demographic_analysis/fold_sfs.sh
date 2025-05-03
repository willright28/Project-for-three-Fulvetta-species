#!/bin/bash
set -e
set -x
set -o pipefail

# Step 1: Identify the values for projecting down each species
nohup ./easySFS/easySFS.py \
  -i ./fsc28/sfs/neutral/region/filter.vcf \
  -p ./fsc28/sfs/species_sets.1.txt \
  --preview -a -v 1> ./fsc28/sfs/neutral/region/proj.log.txt 2>&1 &

# Step 2: Generate the fold SFS using selected projection values
nohup ./easySFS/easySFS.py \
  -i ./fsc28/sfs/neutral/region/filter.vcf \
  -v -a \
  -p ./fsc28/sfs/species_sets.1.txt \
  --proj 34,58,36 \
  -o ./fsc28/sfs/neutral/region \
  --prefix neutral.proj & 