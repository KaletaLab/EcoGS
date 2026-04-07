# EcoGS  
**by Georgios Marinos and Samer Kadib Alban**

A suite to predict ecological relationships between microbial genome-scale metabolic models by comparing growth when models are simulated alone versus when simulated within a community model.

**EcoGS** operates on existing **sybil** [1] metabolic models originated from the **gapseq** [2] pipeline. Further model reconstruction is not required to run the package. Please use [**sybilSBML**](https://www.cs.hhu.de/en/research-groups/computational-cell-biology/software-contributions/sybil) to import **gapseq** models from SBML files, as EcoGS does not operate on **gapseq** models in [**cobrar**](https://github.com/Waschina/cobrar) format.

Growth rates are calculated using flux balance analysis via the **R** [3] package **sybil** [1]. Community models are assembled and simulated using **MicrobiomeGS2** [4].

---

## Ecological interaction classification

Ecological relationships are inferred from growth changes observed when two models are simulated together compared to when simulated alone.

| Species 1 | Species 2 | Relation      |
|-----------|-----------|---------------|
| 0         | 0         | Neutralism    |
| 0         | 1         | Commensalism  |
| 0         | -1        | Amensalism    |
| 1         | 1         | Mutualism     |
| 1         | -1        | Parasitism    |
| -1        | -1        | Competition   |

0: no growth change  
1: increase of growth  
-1: decrease of growth  

---

## Workflow overview

To run EcoGS on user-provided models:

1. `metabolic_interactions_with_MicrobiomeGS2`  
   Simulates growth of individual models and all pairwise combinations.

2. `make_eco_mat`  
   Calculates the ecological interaction matrix and the matrix of growth changes.

3. `plot_relations`  
   Plots the number of model pairs in each ecological category.

4. `relation_per_sample`  
   Calculates ecological interaction frequencies per sample (requires user-provided OTU table).

5. `relation_ratios`  
   Transforms ecological relation frequencies into log10 ratios.

---

## Installation in a conda environment
Create a conda environment and install EcoGS from GitHub:
```r
# create the conda environment
conda env update -f requirements.yml
conda activate EcoGS
# install sybil
Rscript -e "devtools::install_github('SysBioChalmers/sybil')"
# install sybil
Rscript -e "devtools::install_github('SysBioChalmers/sybil-SBML')"
# install MicrobiomeGS2
Rscript -e "devtools::install_github('Waschina/MicrobiomeGS2')"

# install EcoGS without vignette
Rscript -e "devtools::install_github('KaletaLab/EcoGS')"

#OR

# install EcoGS including vignette (longer build time)
Rscript -e "devtools::install_github('KaletaLab/EcoGS', build_vignettes = TRUE)"
```

### Requirements

1. R (>= 4.0.0)  
2. data.table (>= 1.13.4)  
3. sybil (>= 2.1.5)  
4. stringr (>= 1.4.0)  
5. utils  
6. MicrobiomeGS2  
7. glpkAPI  
8. parallel  
9. doParallel  
10. foreach  
11. graphics  

The default solver is **GLPK**. Choosing **CPLEX** is also possible and simulations are faster, but it requires a working installation of **IBM ILOG CPLEX Optimization Studio**.

For installation of **MicrobiomeGS2** and **cplexAPI**:  
https://github.com/Waschina/MicrobiomeGS2/blob/main/README.md  

For installing **sybil** and **sybilSBML**:  
https://www.cs.hhu.de/en/research-groups/computational-cell-biology/software-contributions/sybil  

---

## Technical details

**EcoGS** builds upon **MicrobiomeGS2**. For each pair of bacterial models, two models are merged into a community model and treated as distinct interconnected compartments.

Each model retains its own biomass objective function. Both objectives are optimised simultaneously. No predefined abundances are assumed; relative growth is inferred directly from optimisation of each biomass reaction.

Although such optimisation problems may admit non-unique solutions, **sybil** in combination with the **CPLEX** solver returned stable and reproducible results in practice.

---

## Manual and tutorial

A detailed step-by-step tutorial, including data initialisation and downstream statistical analysis, is provided as a vignette.

To access the manual from within **R**:

```r
library(EcoGS)
vignette("vignette", package = "EcoGS")
```

The SIHUMIx [5] metabolic models utilised in the tutorial are bundled within the package:

```r
path <- system.file("extdata", "SIHUMIx_gapseq_EcoGS.RDS", package = "EcoGS")
models <- readRDS(path)
```

The tutorial is also available as a PDF:  
https://github.com/KaletaLab/EcoGS/blob/main/vignettes/vignette.pdf  

---
## Quick start

The following example runs immediately after installation. Please make sure that you choose the installation method that includes the vignette. It uses existing SIHUMIx **sybil** metabolic models and does **not** require **gapseq** or any external reconstruction pipeline.

```r
# Load EcoGS
library(EcoGS)

# Load the pre-defined SIHUMIx metabolic models
path <- system.file("extdata", "SIHUMIx_gapseq_EcoGS.RDS", package = "EcoGS")
SIHUMIx_gapseq_EcoGS <- readRDS(path)

set.seed(42)
# Create a random abundance table for 8 species across 10 mice
abundance <- as.data.frame(matrix(0,8,10))
row.names(abundance) <- names(SIHUMIx_gapseq_EcoGS)
names(abundance) <- paste0("mouse",1:10)
for (i in 1:8) {
abundance[i,] <- sample(1:100, 10)
}

# Normalise each column by its sum to obtain relative abundance
abundance <- as.data.frame(apply(abundance, 2, function(x) x / sum(x)))

# Step 1: Simulate pairwise metabolic interactions
step1 <- metabolic_interactions_with_MicrobiomeGS2(list_of_models = SIHUMIx_gapseq_EcoGS,
cores = 1, save_pair = T, solver = "glpkAPI")

# Step 2: Construct the ecological matrix
step2 <- make_eco_mat(growth_file = step1, solver = "glpkAPI")

# Step 3: Visualise the predicted relations
step3 <- plot_relations(eco_mat = step2$eco_mat)

# Step 4: Weigh interactions
step4_min <- relation_per_sample(OTU_table = abundance, eco_mat = step2$eco_mat,
weighing_method = "min")
step4_multi <- relation_per_sample(OTU_table = abundance, eco_mat = step2$eco_mat,
weighing_method = "multi" )

# Step 5: Calculate log10 interaction ratios
step5_min <- relation_ratios(relations_table = step4_min$weighed_relations)
step5_multi <- relation_ratios(relations_table = step4_multi$weighed_relations)

```
## References and acknowledgements

We thank the iTREAT consortium and the German Research Foundation (DFG) under Germany`s Excellence Strategy for their support within the Collaborative Research Centre “Origin and Function of Metaorganisms” (CRC 1182), the Research Group miTarget (RU 5042) , the Excellence Cluster for Precision Medicine in Chronic Inflammation (PMI – EXC 2167/2 – 390884018), and the Project ExoMod (518920252).

This research was supported in part through high-performance computing resources available at the Kiel University Computing Centre (DFG Project Number 40395346). We thank Dr Johannes Zimmermann and Dr Robin Koch for helpful discussions and support.

1. Gelius-Dietrich G, Desouki AA, Fritzemeier CJ, Lercher MJ. sybil – Efficient constraint-based modelling in R. BMC Syst Biol. 2013;7:125.  
2. Zimmermann J, Kaleta C, Waschina S. gapseq: informed prediction of bacterial metabolic pathways and reconstruction of accurate metabolic models. Genome Biol. 2021;22:81. 
3. R Core Team. R: A Language and Environment for Statistical Computing. Vienna, Austria: R Foundation for Statistical Computing.  
4. Waschina S. Analysis and simulation of genome- & ecosystem-scale microbial metabolism. (https://github.com/Waschina/MicrobiomeGS2/)  
5. Becker N, Kunath J, Loh G, Blaut M. Human intestinal microbiota: characterisation of a simplified and stable gnotobiotic rat model. Gut Microbes. 2011;2:25–33.

---

## Citation

Georgios Marinos, Karlis Arturs Moors, Kristina Schlicht, Malte Rühlemann, Silvio Waschina, Wolfgang Lieb, Andre Franke, Matthias Laudes, Mathieu Groussin, Mathilde Poyet, Christoph Kaleta*, and A. Samer Kadibalban*, 2025.  
Genome-scale metabolic models predict diet- and lifestyle-driven shifts of ecological interactions in the gut microbiome.  
https://doi.org/10.1101/2025.09.23.678088  

*Shared corresponding authors; enquiries to c.kaleta@iem.uni-kiel.de and/or s.kadibalban@iem.uni-kiel.de*

---

GNU General Public License version 3.0 (GPLv3) is applied to all copyrightable parts of this software.

Contact details:  
https://www.iem.uni-kiel.de/de/medizinische-systembiologie/medizinische-systembiologie
