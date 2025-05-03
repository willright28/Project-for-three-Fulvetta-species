#!/bin/bash
set -x

wkdir=${1} ## Working directory
run=${2}  ## Run ID
nsites=${3} ## Number of climate-adaptive SNPs (non-linked)
effect_size=${4} ## Effect size of mutations
esd=${5} # SD of environmental effect on fitness (non-heritable component)
sigma=${6} ## Breadth of tolerance (plasticity)
env_data=${7} ## Environmental data

cat > ${wkdir}/run${run}/slim.run${run}.txt << EOM

initialize() {
    // Define constants
    defineConstant("startTime", clock()); 
    defineConstant("migMatrix", "./mig.txt");   // Migration matrix
    defineConstant("Environment", "${env_data}");   // Environmental variable data file
    initializeSLiMOptions(preventIncidentalSelfing=T); 
    defineConstant("Mrate", 2.3e-8); // Mutation rate
    initializeMutationRate(Mrate*100); 
    defineConstant("popscale", 1); // Scaling factor for testing
    defineConstant("extinctionLimit", 0.01); // Threshold of mean fitness below which population is considered extinct

    defineConstant("esd", ${esd}); // SD of environmental noise in phenotype
    defineConstant("pl", ${sigma}); // Phenotypic plasticity (SD of stabilizing selection curve)
    defineConstant("mutEffect", ${effect_size}); // SD of mutation effect size at Climate-adaptive mutations (mean = 0)
    defineConstant("C", ${nsites}); // Number of climate-adaptive mutations
    defineConstant("N", 1); // Number of neutral unlinked sites

    defineConstant("maxfit", dnorm(0.0, 0.0, pl)); 
    defineConstant("fitcushion", 0.001); // Minimum possible fitness to prevent zero
    
    // Define mutation and genomic element types
    initializeMutationType("m1", 0.5, "f", 0.0); // Neutral mutations
    initializeMutationType("m2", 0.5, "n", 0, mutEffect); // Climate-adaptive mutations
    initializeGenomicElementType("g1", m1, 1.0);
    initializeGenomicElementType("g2", m2, 1.0);
    initializeGenomicElement(g2, 0, C - 1);
    initializeGenomicElement(g1, C, C + N - 1);  
    initializeRecombinationRate(0.5);        
    m2.mutationStackPolicy = "s"; 	
    m2.convertToSubstitution = F;	
}

// Read migration matrix with population sizes on the diagonal
1 early() { 
    migrations = readFile(migMatrix);
    defineConstant("npops", size(migrations)); 
    
    // Read population sizes
    popsizes = NULL;
    for (p in 0:(npops - 1)) {
        ms = strsplit(migrations[p], sep="\t");
        popsizes = c(popsizes, asInteger(ms[p]));
    }
    sim.setValue("popsizes", popsizes); 
    
    // Create subpopulations
    for (p in 0:(npops - 1)) {
        sim.addSubpop(p, asInteger(popscale * popsizes[p] / 100));
    }
    sim.setValue("deadpops", rep(0, npops)); 

    // Set migration rates
    for (p in 0:(npops - 1)) {
        ms = asFloat(strsplit(migrations[p], sep="\t"));
        sim.subpopulations[p].setValue("popmig", ms); 
        for (m in 0:(npops - 1)) {
            if (m != p) {
                sim.subpopulations[m].setMigrationRates(sim.subpopulations[p], ms[m]);
            }
        }
    }
}

// Load environmental data (pops as rows, generations as columns)
1 early() { 
    env = readFile(Environment);
    for (i in 0:(npops - 1)) {
        sim.setValue(paste(c(i, ".env"), sep=""), asFloat(strsplit(env[i + 1], sep="\t"))); 
    }

    // Set fitness callbacks using stabilizing selection model
    for (i in 0:(npops - 1)) {
        envValues = sim.getValue(paste(c(i, ".env"), sep=""));
        for (gen in 1:20200) {
            env = envValues[(gen - 1)];
            src = "{ return " + fitcushion + " + " + "dnorm(" + env + " - individual.tagF, 0.0, " + pl + ") / " + (maxfit / (1 - fitcushion)) + "; }";
            sim.registerFitnessEffectCallback(NULL, src, i, gen, gen); 
        }
    }
}

// Mutation effect
mutationEffect(m2) { return 1.0; }

// Calculate phenotype (QTLs + noise)
1: late() {
    inds = sim.subpopulations.individuals;
    inds.tagF = inds.sumOfMutationsOfType(m2) + rnorm(size(inds), 0, esd);
    writeFile("${wkdir}/run${run}/cycle.txt", paste0(sim.cycle));
}

// Adjust population size scaling at generation 10000
10000 late() {
    for (p in 0:(npops - 1)) {
        sim.subpopulations[p].setSubpopulationSize(asInteger(popscale * sim.getValue("popsizes")[p] / 10));
    }
    sim.chromosome.setMutationRate(Mrate * 10);
}

// Restore population size scaling at generation 15001
15001 late() {
    for (p in 0:(npops - 1)) {
        sim.subpopulations[p].setSubpopulationSize(asInteger(popscale * sim.getValue("popsizes")[p]));
    }
    sim.chromosome.setMutationRate(Mrate);
}

// Output data
19951: early() {
    m2_obj = unique(sim.mutationsOfType(m2));
    all_m2_number = size(m2_obj);
    new_m2_number = size(m2_obj[m2_obj.originTick >= 20000]);

    all_ai_site = c();
    for (p in 0:(size(sim.subpopulations) - 1)) {
        pop = sim.subpopulations[p];
        mut_obj = unique(pop.individuals.genomes.mutationsOfType(m2));
        ai_site = mut_obj[mut_obj.subpopID != p];
        all_ai_site = unique(c(ai_site, all_ai_site));

        cat(sim.cycle + "\t" + p + "\t" + mean(pop.cachedFitness(NULL)) + "\t" + pop.individualCount + "\t" + mean(pop.individuals.tagF) + "\t" + size(mut_obj) + "\t" + sd(pop.individuals.sumOfMutationsOfType(m2)) + "\n");
        fitness_out_string = paste0(sim.cycle + "\t" + p + "\t" + mean(pop.cachedFitness(NULL)) + "\t" + pop.individualCount + "\t" + mean(pop.individuals.tagF) + "\t" + size(mut_obj) + "\t" + sd(pop.individuals.sumOfMutationsOfType(m2)));
        writeFile(paste0("${wkdir}/run${run}/pop_fitness_", ${run}, ".txt"), fitness_out_string, append=T);
    }

    all_ai_site_number = size(unique(all_ai_site));
    cat(sim.cycle + '\t' + all_m2_number + '\t' + all_ai_site_number + '\t' + all_ai_site_number / all_m2_number + '\t' + new_m2_number + '\t' + mean(sim.subpopulations.cachedFitness(NULL)) + '\n');
    site_out_string = paste0(sim.cycle + '\t' + all_m2_number + '\t' + all_ai_site_number + '\t' + all_ai_site_number / all_m2_number + '\t' + new_m2_number + '\t' + mean(sim.subpopulations.cachedFitness(NULL)));
    writeFile(paste0("${wkdir}/run${run}/mutation_sites_", ${run}, ".txt"), site_out_string, append=T);
}

// Record average fitness for adjusting population sizes
19989 early() {
    sim.setValue("ps", rep(0, size(sim.subpopulations)));
}
19989 : 19999 early() {
    means = sapply(sim.subpopulations, "mean(applyValue.cachedFitness(NULL));");
    sim.setValue("ps", sim.getValue("ps") + means);
    sim.setValue("fc", rep(1.0, size(sim.subpopulations)));
    if (sim.cycle == 19990) {
        fc = sim.getValue("fc");
        for (i in 0:(size(sim.subpopulations) - 1)) {
            fc[i] = sim.getValue("ps")[i] / 10;
        }
        sim.setValue("fc", fc);
    }
}

// Adjust population sizes based on fitness above
20000 : early() {
    subpops = sim.subpopulations;
    popsizes = sim.getValue("popsizes");
    deadpops = sim.getValue("deadpops");
    for (i in 0:(size(subpops) - 1)) {
        pop = subpops[i];
        if (deadpops[i] == 1) {
            next;
        }
        popsize = popsizes[i];
        meanFitness = mean(pop.cachedFitness(NULL));
        newSize = asInteger(popscale * popsize * meanFitness / sim.getValue("fc")[i]);
        if (meanFitness < extinctionLimit & deadpops[i] == 0) {
            deadpops[i] = 1;
            sim.setValue("deadpops", deadpops);
            cat("#" + sim.cycle + "\t" + i + "\textinct.\n");
            newSize = 10;
        }
        if (sum(deadpops) == size(subpops)) {
            cat("#" + sim.cycle + "\t" + "ALL\textinct.\n");
            sim.simulationFinished();
        }
        pop.setSubpopulationSize(newSize);

        // set migration rates for extinct populations
        for (j in 0:(size(subpops) - 1)) {
            if (i != j) {
                if (deadpops[i] == 1 | deadpops[j] == 1) {
                    subpops[j].setMigrationRates(subpops[i], 0); 
                    subpops[i].setMigrationRates(subpops[j], 0); 
                }
            }
        }
    }
}

20200 late() {}

EOM
