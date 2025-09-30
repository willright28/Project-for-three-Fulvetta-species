#!/bin/bash
###################
workdir=/disk3/call_snp/vcf96
ref=/disk3/ref/CL100120897_L02.gapcloser.FINAL.fasta

#SNP calling was run in parallel on split chromosome files, and the per-chromosome VCFs were concatenated to  the final whole-genome VCF
awk '{print $1, 1 ,$2}' /disk3/ref/CL100120897_L02.gapcloser.FINAL.fasta.fai > ${workdir}/split.txt

awk '{print $1}' ${workdir}/id.list.txt | sed "s:^:/disk3/realig_bam/:" |sed "s:$:.mem_map.sort.mkd.realign.bam:" > ${workdir}/bam.out.list.txt
while read line ; do { if [ -f $line ];then echo "yes" ; else echo "no" ; fi }; done < ${workdir}/bam.out.list.txt


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
#call snp
    bcftools mpileup -f $ref -r ${chr}:${sta}-${end} -b ${workdir}/bam.out.list.txt \
    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR | bcftools call -m -Oz -f GQ -o ${workdir}/raw_vcf_out/${chr}.bcftools.chr.vcf.gz
    echo >&9
}&
fi
done< ${workdir}/split.txt
wait

awk '{print $4}' ${workdir}/split.txt | sed "s:^:${workdir}/raw_vcf_out/:" | sed "s:$:.bcftools.chr.vcf.gz:" > ${workdir}/chr.raw.list.txt 
time bcftools concat --threads 70 -Oz -o  ${workdir}/all.bcftools.raw.out.vcf.gz -f   ${workdir}/chr.raw.list.txt

exec 9>$-

