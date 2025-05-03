#!/bin/bash

#fratercula and davidi
msmove 52 2000 -t 78.052 -r 76.781 10000 -I 2 18 34 -n 1 0.736 -n 2 1.913  -en 0.241 2 2.672 -ej 0.243 1 2 -en 0.243 2 1 > ./FILET/simdata/noMig.w2c.msOut &
msmove 52 2000 -t 78.052 -r 76.781 10000 -I 2 18 34 -n 1 0.736 -n 2 1.913 -m 1 2 0.8  -m 2 1 0 -em 0.0511 1 2 0  -em 0.0511 2 1 0 -en 0.241 2 2.672 -ej 0.243 1 2 -en 0.243 2 1 > ./FILET/simdata/mig12.w2c.msOut &
msmove 52 2000 -t 78.052 -r 76.781 10000 -I 2 18 34 -n 1 0.736 -n 2 1.913 -m 2 1 1.824 -m 1 2 0 -em 0.0511 1 2 0  -em 0.0511 2 1 0 -en 0.241 2 2.672 -ej 0.243 1 2 -en 0.243 2 1 > ./FILET/simdata/mig21.w2c.msOut &

#hueti and davidi
msmove 53 2000 -t 208.632 -r 210.446 10000 -I 2 19 34 -n 1 0.525 -n 2 0.716 -ej 0.09 2 1 -en 0.09 1 1 > ./FILET/simdata/noMig.model5.msOut &
msmove 53 2000 -t 208.632 -r 210.446 10000 -I 2 19 34 -n 1 0.525 -n 2 0.716 -m 1 2 7.668  -m 2 1 0 -em 0.019 1 2 0  -em 0.019 2 1 0 -ej 0.09 2 1 -en 0.09 1 1 > ./FILET/simdata/mig12.e2c.msOut &
msmove 53 2000 -t 208.632 -r 210.446 10000 -I 2 19 34 -n 1 0.525 -n 2 0.716 -m 2 1 11.663 -m 1 2 0 -em 0.019 1 2 0  -em 0.019 2 1 0 -ej 0.09 2 1 -en 0.09 1 1 > ./FILET/simdata/mig21.e2c.msOut &
