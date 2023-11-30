# EcoGS
_**by Samer Alban Kadib, Georgios Marinos, Johannes Zimmermann, Silvio Waschina, and Christoph Kaleta**_

A suite to predict the ecological relationships between gapseq models by comparing the achieved growth when the models are alone and the achieved growth when the models are part of a community model.

|Species 1|Species 2|Relation|
|-----|-----|-----|
|0|0|Neutralism|
|0|1|Commensalism|
|0|-1|Amensalism|
|1|1|Mutualism|
|1|-1|Parasitism|
|-1|-1|Competition|

0: no growth change, 1: increase of growth, 2: decrease of growth

## Quick start
To run the suite, run the following functions:

1. `metabolic_interactions_with_MicrobiomeGS2()` to simulate growth of the models
2. `make_eco_mat` to calculate the ecological matrix and the matrix of growth change for all models
3. `plot_relations` to plot the number of pairs in each type of ecological relations
4. `relation_per_sample` to calculate the prediction of the frequency of ecological relationship between pairs of co-grown bacteria in each sample. It requires a user imported OTU table.

## Installation
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

## References and acknowlegments
This research was supported in part through high-performance computing resources available at the Kiel University Computing Centre.

## Citation
Alban Kadib, Marinos _et al._, **Genome–scale metabolic models predict associations among human dietary compounds and microbial ecological interactions**, in preparation
