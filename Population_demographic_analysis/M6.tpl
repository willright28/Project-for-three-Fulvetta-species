//Parameters for the coalescence simulation program : fastsimcoal.exe
3 samples to simulate :
//Population effective sizes (number of genes)
Nhueti$
Ndavidi$
Nfratercula$
//Haploid samples sizes 
34
58
36
//Growth rates: negative growth implies population expansion
0
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
0	   MIG01$	   0  	
MIG10$ 0     MIG12$
0    MIG21$    0
//Migration matrix 1
0	0	0  	
0 0 0
0 0 0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
3 historical event
Tm1$ 0 0 0 1 0 1
Tdiv_EC$ 0 1 1 RES_EC$ 0 1
Tdiv_WC$ 1 2 1 RES_WC$ 0 1
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 2.3e-8 OUTEXP
