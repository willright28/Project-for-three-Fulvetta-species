#!/bin/bash
###################
###################

set -e
set -x 
set -o pipefail

workdir=/disk3/hkqm/call_snp/vcf96
ref=/disk3/zhangshuai3/hkqm/ref_hkqm/CL100120897_L02.gapcloser.FINAL.fasta

#awk '{print $1, 1 ,$2}' /disk3/zhangshuai3/hkqm/ref_hkqm/CL100120897_L02.gapcloser.FINAL.fasta.fai > /disk3/hkqm/call_snp/vcf96/split.txt

#awk '{print $1}' /disk3/hkqm/call_snp/vcf96/callsnp.list.txt | sed "s:^:/disk3/zhangshuai3/hkqm/realig_bam/:" |sed "s:$:.mem_map.sort.mkd.realign.bam:" > /disk3/hkqm/call_snp/vcf96/bam.out.list.txt

#while read line ; do { if [ -f $line ];then echo "yes" ; else echo "no" ; fi }; done < /disk3/hkqm/call_snp/vcf96/bam.out.list.txt

thread=60
temfifo=$$.fifo
mkfifo $temfifo
exec 9<>$temfifo
rm $temfifo

for ((i=1;i<=$thread;i++))
do
    echo >&9
done

while read chr sta end
do
if [ $chr ]
then
read -u9
{

#分开callsnp
#mpileup:producing genotype likelihoods;#http://samtools.github.io/bcftools/bcftools.html

time  bcftools mpileup -f $ref -r ${chr}:${sta}-${end} -b ${workdir}/bam.out.list.txt \
--annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR | bcftools call -m -Oz -f GQ -o ${workdir}/raw_vcf_out/${chr}.bcftools.chr.vcf.gz
    echo >&9
}&
fi
done< /disk3/hkqm/call_snp/vcf96/split.more.txt
wait

awk '{print $4}' /disk3/hkqm/call_snp/vcf96/split.more.txt | sed "s:^:/disk3/hkqm/call_snp/vcf96/raw_vcf_out/:" | sed "s:$:.bcftools.chr.vcf.gz:" > /disk3/hkqm/call_snp/vcf96/chr.raw.list.txt 

time bcftools concat --threads 70 -Oz -o /disk3/hkqm/call_snp/vcf96/htcw.bcftools.raw.out.vcf.gz -f  /disk3/hkqm/call_snp/vcf96/chr.raw.list.txt

exec 9>$-

