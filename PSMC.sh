#!/bin/bash
ref=$1
id=$2
bcftools mpileup -C 50 -q 20 -Ou -f ${ref} ./realig_bam/${id}.mem_map.sort.mkd.realign.bam | bcftools call -c - |\
/home/share/bcftools-1.9/misc/vcfutils.pl vcf2fq -d 4 -D 50 > ./psmc_data/${id}.psmc.fq

fq2psmcfa -q20 ./psmc_data/${id}.psmc.fq > ./psmc_data/${id}.psmcfa
psmc -N30 -t5 -r5 -p "4+30*2+4+6+10" -o ./psmc_data/${id}.psmc ./psmc_data/${id}.psmcfa
