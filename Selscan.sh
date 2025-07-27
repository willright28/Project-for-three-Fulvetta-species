#!/bin/bash
#######################ihs########################

cd /disk3/hkqm/select_sweep
while read line
do
 		while (true)
 			do
		number=`ps -ef | grep /disk3/hkqm/select_sweep/ihs | wc -l`
 		if [ $number -lt 10 ];
 		then

 /home/share/selscan/bin/linux/selscan --ihs \
 --vcf /disk3/hkqm/phase/phase_output/central/central.${line}.recode.vcf  \
 --map /disk3/hkqm/phase/phase_output/central/central.${line}.recode.map \
 --threads 50 \
 --pmap  \
 --out /disk3/hkqm/select_sweep/ihs_central/${line}.central.ihs &

 /home/share/selscan/bin/linux/selscan --ihs \
 --vcf /disk3/hkqm/phase/phase_output/western/western.${line}.recode.vcf  \
 --map /disk3/hkqm/phase/phase_output/western/western.${line}.recode.map \
 --threads 50 \
 --pmap  \
 --out /disk3/hkqm/select_sweep/ihs_western/${line}.western.ihs &

 /home/share/selscan/bin/linux/selscan --ihs \
 --vcf /disk3/hkqm/phase/phase_output/eastern/eastern.${line}.recode.vcf  \
 --map /disk3/hkqm/phase/phase_output/eastern/eastern.${line}.recode.map \
 --threads 50 \
 --pmap  \
 --out /disk3/hkqm/select_sweep/ihs_eastern/${line}.eastern.ihs &

 		break
 		fi
 		sleep 30s
 done
 done < /disk3/hkqm/phase/chr_filter.txt

 wait

 /home/share/selscan/bin/linux/norm --ihs --files /disk3/hkqm/select_sweep/ihs_eastern/*.eastern.ihs.ihs.out --bp-win --winsize 100000 
 /home/share/selscan/bin/linux/norm --ihs --files /disk3/hkqm/select_sweep/ihs_western/*.western.ihs.ihs.out --bp-win --winsize 100000 
 /home/share/selscan/bin/linux/norm --ihs --files /disk3/hkqm/select_sweep/ihs_central/*.central.ihs.ihs.out --bp-win --winsize 100000 

#########################nsl##############################
!/bin/bash
 cd /disk3/hkqm/select_sweep

  while read line
  do
  		while (true)
 			do
		number=`ps -ef | grep /disk3/hkqm/select_sweep/nsl | wc -l`
 		if [ $number -lt 10 ];
  		then

 /home/share/selscan/bin/linux/selscan --nsl \
 --vcf /disk3/hkqm/phase/phase_output/central/central.${line}.recode.vcf  \
 --threads 50 \
 --out /disk3/hkqm/select_sweep/nsl_central/${line}.central.nsl &

 /home/share/selscan/bin/linux/selscan --nsl \
 --vcf /disk3/hkqm/phase/phase_output/western/western.${line}.recode.vcf  \
 --threads 50 \
 --out /disk3/hkqm/select_sweep/nsl_western/${line}.western.nsl &

 /home/share/selscan/bin/linux/selscan --nsl \
 --vcf /disk3/hkqm/phase/phase_output/eastern/eastern.${line}.recode.vcf  \
 --threads 50 \
 --out /disk3/hkqm/select_sweep/nsl_eastern/${line}.eastern.nsl &

  		break
  		fi
 		sleep 30s
 done
 done < /disk3/hkqm/phase/chr_filter.txt
