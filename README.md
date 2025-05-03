# Project for Three Fulvetta Species

Scripts for Zhang et al. (2025): *Interspecific introgression mitigates climate change risk of the mountainous birds*

---

## 1. Genetic Structure Analysis

- **NJ_tree.sh**: Generates a pairwise genetic distance matrix using VCF2Dis and constructs a neighbor-joining (NJ) tree using FastME.
- **PCA_Admixture.sh**: Performs principal component analysis (PCA) and admixture analysis.

---

## 2. Population Demographic Analysis

- **fold_sfs.sh**: Generates the folded site frequency spectrum (SFS).
- **M\*.tpl**: Fastsimcoal2 model template files.
- **M\*.est**: Fastsimcoal2 model parameter files.
- **fsc28.sh**: Runs demographic models using Fastsimcoal2.

---

## 3. Introgression Analysis

- **msmove_sim.sh**: Generates training data through coalescent simulations using parameters inferred from Fastsimcoal2.
- **FILET.sh**: Identifies introgressed genomic regions using FILET.

---

## 4. Ecological Analysis

- **Niche_divergence_test_analysis.R**: Evaluates ecological niche similarity between species.
- **Ecological_niche_models_analysis.R**: Predicts species distribution areas using ensemble ecological niche models.

---

## 5. Genotype–Climate Association Analysis

- **LFMM2.r**: Identifies SNPs associated with climatic variables using LFMM2.
- **RDA.r**: Identifies SNP–environment associations using redundancy analysis (RDA).

---

## 6. Allele Frequency Shift Analysis

- **Allele_frequency_shifts_analysis.R**: Models the rate of allele frequency change under future climate scenarios.

---

## 7. Genetic Offset Analysis

- **Genetic_offset_analysis.R**: Calculates genetic offset using two methods: gradient forest and geometric offset approaches.

---

## 8. SLiM4 Simulation Analysis

- **mig.txt**: Migration matrix derived from Fastsimcoal2 demographic inference.
- **env_example_data.txt**: Linear climate change input used to simulate future warming.
- **slim_geneflow_model.sh**: SLiM simulation script allowing interspecific gene flow.
- **slim_none_geneflow_model.sh**: SLiM simulation script prohibiting interspecific gene flow.

---

## 9. Landscape connectivity analysis

- **Landscape_connectivity_analysis.R**: Tests a set of alternative hypotheses regarding landscape 
connectivity

