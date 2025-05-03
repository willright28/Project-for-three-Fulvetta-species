#!/bin/bash
##Using Fratercula and davidi as an example, the same code was applied to the other species pair
#1. Calculate summary statistics for each input file in our training simulations
cat ./FILET/simdata/noMig.w2c.msOut | ./FILET/FILET-master/twoPopnStats_forML 18 34 | python2 ./FILET/FILET-master/normalizeTwoPopnStats.py None 10000 > ./FILET/test/trainingSimsStats/noMig.msOut 
cat ./FILET/simdata/mig12.w2c.msOut | ./FILET/FILET-master/twoPopnStats_forML 18 34 | python2 ./FILET/FILET-master/normalizeTwoPopnStats.py None 10000 > ./FILET/test/trainingSimsStats/mig12.msOut 
cat ./FILET/simdata/mig21.w2c.msOut | ./FILET/FILET-master/twoPopnStats_forML 18 34 | python2 ./FILET/FILET-master/normalizeTwoPopnStats.py None 10000 > ./FILET/test/trainingSimsStats/mig21.msOut 

#2. Aggregate our training data into one labeled data matrix to be used for training.
python2  ./FILET/FILET-master/buildThreeClassTrainingSet.py ./FILET/test/trainingSimsStats/ ./FILET/test/trainingSets/threeClass.fvec
#3. Train the classifier.
python2 ./FILET/FILET-master/trainFiletClassifier.py ./FILET/test/trainingSets/threeClass.fvec ./FILET/test/classifier/threeClass.p pi1 pi2 tajd1 tajd2 hetVar1 hetVar2 ss1 ss2 ddRank1 ddRank2 zx Fst gmin dxy_mean dxy_min snn 

#4. Input  real phased sequence data set
#Conducted by each chromosome
while read line
do
./FILET/FILET-master/pgStatsBedSubpop_forML  ./phase/fa.split/eastern_fa/eastern.${line}.fa  ./phase/fa.split/central_fa/central.${line}.fa  ./phase/fa.split/hkqm.${line}.fasta  ./phase/fa.split/${line}.win.bed 0.5 | python2 ./FILET/FILET-master/normalizePgStats.py 10000 > ./FILET/test/featureVectorsToClassify/e2c.${line}.ss
done <  /disk3/hkqm/phase/chr_.txt

#5. Use classifier to produce classifications on the real data set
python2 /disk3/hkqm/FILET/FILET-master/classifyChromosome.py ./FILET/test/classifier/threeClass.p ./FILET/test/featureVectorsToClassify/  0.01 /disk3/hkqm/FILET/test/results
