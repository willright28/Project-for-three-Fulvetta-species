#!/bin/bash
set -e
set -x 
set -o pipefail

model_pre="model_prefix"

#run Fastsimcoal
for i in {1..50}
do
{
mkdir -p run$i
   cp ${model_pre}.tpl ${model_pre}.est ${model_pre}_jointMAFpop1_0.obs ${model_pre}_jointMAFpop2_0.obs  ${model_pre}_jointMAFpop2_1.obs run$i"/"
   cd run$i
./fsc28/fsc28_linux64/fsc28 -t ${model_pre}.tpl -e ${model_pre}.est -m -0 -C 10 -n 50000 -L 40 -s 0 -M -y4 -c3
cd ..
}
done

wait

echo "finish"






