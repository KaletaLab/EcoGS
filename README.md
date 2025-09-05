# EcoGS **by Georgios Marinos and Samer Kadib Alban**

A suite to predict the ecological relationships between gapseq models by comparing the achieved growth when the models are alone and the achieved growth when the models are part of a community model. The software supports genome-scale metabolic models that are reconstructed using gapseq [1]. The achieved growth of each model was calculated models based on flux balance analysis using the R [2] package sybil [3]. The community model was assembled and its growth was assessed using the R package MicrobiomeGS2 [4].


|Species 1|Species 2|Relation|
|-----|-----|-----|
|0|0|Neutralism|
|0|1|Commensalism|
|0|-1|Amensalism|
|1|1|Mutualism|
|1|-1|Parasitism|
|-1|-1|Competition|

0: no growth change, 1: increase of growth, -1: decrease of growth

## Quick start
To run the suite, run the following functions:

1. `metabolic_interactions_with_MicrobiomeGS2` to simulate growth of the models
2. `make_eco_mat` to calculate the ecological matrix and the matrix of growth change for all models
3. `plot_relations` to plot the number of pairs in each type of ecological relations
4. `relation_per_sample` to calculate the prediction of the frequency of ecological relationship between pairs of co-grown bacteria in each sample. It requires a user imported OTU table.
5. `relation_ratios` to transform the ecological relation frequencies into log10 ratios.

## Installation

```r
library(devtools)
install_github("KaletaLab/EcoGS")
```

It requires the following software / R packages:
  
  1. R (>= 4.0.0)
  2. data.table (>= 1.13.4)
  3. sybil (>= 2.1.5)
  4. stringr (>= 1.4.0)
  5. utils
  6. MicrobiomeGS2
  7. cplexAPI
  8. parallel
  9. doParallel
  10. foreach
  11. graphics

For further instructions regarding the installation of MicrobiomeGS2 and cplexAPI, please consult the following page: https://github.com/Waschina/MicrobiomeGS2/blob/main/README.md


## Technical details
We used MicrobiomeGS2 as the basis for EcoGS, in each case, two bacterial models were merged into one and regarded as interconnected distinct compartments of the merged model.
Importantly, we kept the objective function of each model independent from the other and optimised each objective simultaneously.
The objective function was the biomass reaction of each model. Therefore, we did not assume any abundance but we inferred it from the optimization of each biomass reaction. Depending on both resulting growth rates, we estimated the ecological interactions. Although there are non-unique solutions to such optimization problems, the sybil software in R and the CPLEX solver return stable results. Therefore, in practice, we did not face this issue.



## References and acknowlegments
This research was supported in part through high-performance computing resources available at the Kiel University Computing Centre.

1. Zimmermann J, Kaleta C, Waschina S. gapseq: informed prediction of bacterial metabolic pathways and reconstruction of accurate metabolic models. Genome Biol. 2021 Mar 10;22(1):81.
2. R Core Team. R: A Language and Environment for Statistical Computing [Internet]. Vienna, Austria: R Foundation for Statistical Computing; Available from: https://www.R-project.org/
3.	Gelius-Dietrich G, Desouki AA, Fritzemeier CJ, Lercher MJ. sybil – Efficient constraint-based modelling in R. BMC Syst Biol. 2013 Dec;7(1):125.
4.	Waschina S. Analysis and simulation of genome- & ecosystem-scale microbial metabolism [Internet]. Available from: www.github.com/Waschina/MicrobiomeGS2


## Citation
Georgios Marinos, Karlis Arturs Moors, Malte Rühlemann, Silvio Waschina, Wolfgang Lieb, Andre Franke, Matthias Laudes, Mathieu Groussin, Mathilde Poyet, Christoph Kaleta*, and A. Samer Kadibalban*, 2025, Genome-scale metabolic models predict diet- and lifestyle-driven shifts of ecological interactions in the gut microbiome, in preparation

*Shared corresponding authors; inquires to c.kaleta@iem.uni-kiel.de and/or s.kadibalban@iem.uni-kiel.de

GNU General Public License version 3.0 (GPLv3) is applied to all copyrightable parts of this software.

Contact details: https://www.iem.uni-kiel.de/de/medizinische-systembiologie/medizinische-systembiologie


