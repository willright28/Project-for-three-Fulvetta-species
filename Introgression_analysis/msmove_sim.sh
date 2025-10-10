#!/bin/bash
Training data were generated for each of the
454 three classes through coalescent simulations using demographic parameters inferred by Fastsimcoal
under the best-fit model using Msmove (https://github.com/geneva/msmove).

#Use msmove to generate training data for FILET. See https://github.com/genevalab/msmove for more details on how to use msmove.
#fratercula and davidi
msmove 104 10000 -t 7.805 -r 7.678 1000 -I 2 36 68 -n 1 0.736 -n 2 1.913  -en 0.241 2 2.672 -ej 0.243 1 2 -en 0.243 2 1 > /disk3/hkqm/FILET/simdata/1kb/noMig.w2c.msOut &
msmove 104 10000 -t 7.805 -r 7.678 1000 -I 2 36 68 -n 1 0.736 -n 2 1.913 -m 1 2 8  -m 2 1 0 -em 0.0511 1 2 0  -em 0.0511 2 1 0 -en 0.241 2 2.672 -ej 0.243 1 2 -en 0.243 2 1 >  /disk3/hkqm/FILET/simdata/1kb/mig12.w2c.msOut  &
msmove 104 10000 -t 7.805 -r 7.678 1000 -I 2 36 68 -n 1 0.736 -n 2 1.913 -m 2 1 18.2 -m 1 2 0 -em 0.0511 1 2 0  -em 0.0511 2 1 0 -en 0.241 2 2.672 -ej 0.243 1 2 -en 0.243 2 1 >  /disk3/hkqm/FILET/simdata/1kb/mig21.w2c.msOut &

#hueti and davidi
msmove 106 10000 -t 20.8632 -r 21.0446 1000 -I 2 38 68 -n 1 0.525 -n 2 0.716 -ej 0.09 2 1 -en 0.09 1 1 > /disk3/hkqm/FILET/simdata/1kb/noMig.e2c.model5.msOut &
msmove 106 10000 -t 20.8632 -r 21.0446 1000 -I 2 38 68 -n 1 0.525 -n 2 0.716 -m 1 2 7.668  -m 2 1 0 -em 0.019 1 2 0  -em 0.019 2 1 0 -ej 0.09 2 1 -en 0.09 1 1 > /disk3/hkqm/FILET/simdata/1kb/mig12.e2c.model5.msOut &
msmove 106 10000 -t 20.8632 -r 21.0446 1000 -I 2 38 68 -n 1 0.525 -n 2 0.716 -m 2 1 11.663 -m 1 2 0 -em 0.019 1 2 0  -em 0.019 2 1 0 -ej 0.09 2 1 -en 0.09 1 1 > /disk3/hkqm/FILET/simdata/1kb/mig21.e2c.model5.msOut &



