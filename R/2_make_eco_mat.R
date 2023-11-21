#' Calculation of ecological relationships
#'
#' Ecological relationship between pairs of co-grown bacteria
#' based on the list that is exported by the function metabolic_interactions_with_MicrobiomeGS2
#'
#' @param growth_file a list, as exported by the function metabolic_interactions_with_MicrobiomeGS2 and contains:
#'
#'  a. The result of the biomass optimization of each model alone
#'
#'  b. The result of the biomass optimization and the metabolic interactions for each pair of models
#' @param alpha the margin in the relative growth (co_growth - alone_growth)/alone_growth that would be counted as a considerable change, default is 0.05
#' @return
#' (a list object in the R environment)
#'
#' a. matrix of all_bacteria X all_bacteria with each cell represents the ecological relation between the two bacteria
#'
#' b. matrix of the change in growth for all models in columns in the presence of all other models in rows
#'
#' @details
#' A reference table for developers:
#'
#' | Species 1| Species 2| Sum  | Relation  |
#' | ------------- |:-------------:| :-----:|-----:|
#' | 0      | 0 | 0 |Neutralism|
#' | 0      | 1      |   1 |Commensalism|
#' | 0      | -1      |   -1 |Amensalism|
#' | 1      | 1      |   2 |Mutualism|
#' | 1      | -1      |   0 |Parasitism|
#' | -1      | -1      |   2 |Competition|
#' @md
#' @export

make_eco_mat <- function(growth_file, alpha = 0.05){
  growth_list <- growth_file
  pair_growth <- growth_list[["growth_in_pairs"]]
  pair_growth$Model_1_Growth[pair_growth$Model_1_Growth < 10^-6] = 0
  pair_growth$Model_2_Growth[pair_growth$Model_2_Growth < 10^-6] = 0
  single_growth <- growth_list[["single_growth"]]
  single_remove_ind <- single_growth$Biomass_rate_alone < 10^-6
  single_remove_names <- single_growth$Model.name[single_remove_ind]
  pair_remove_ind <- pair_growth$Model_1 %in% single_remove_names | pair_growth$Model_2 %in% single_remove_names
  single_growth <- single_growth[!single_remove_ind,]
  pair_growth <- pair_growth[!pair_remove_ind,]
  untrusted_pairs <- pair_growth[pair_growth$Optimisation_status != 1,]
  growth_mat <- as.data.frame(matrix(NA, nrow(single_growth), nrow(single_growth)))
  rownames(growth_mat) <- single_growth$Model.name
  colnames(growth_mat) <- single_growth$Model.name
  for (i in 1:nrow(pair_growth)) {
    if (i %% 1000 == 0) {
      print(paste0("Reading, pair number: ", i))
    }
    sp1 <- pair_growth$Model_1[i]
    sp2 <- pair_growth$Model_2[i]
    growth_mat[sp2, sp1] <- pair_growth$Model_1_Growth[i]
    growth_mat[sp1, sp2] <- pair_growth$Model_2_Growth[i]
  }
  change_mat <- as.data.frame(t(t(growth_mat) - single_growth$Biomass_rate_alone))
  change_mat[is.na(change_mat)] <- 0
  relative_change_mat <- as.data.frame(t(t(abs(change_mat))/single_growth$Biomass_rate_alone))
  change_mat[relative_change_mat < alpha] <- 0
  sign_mat <- as.matrix(sign(change_mat))
  if (nrow(untrusted_pairs)) {
    for (l in 1:nrow(untrusted_pairs)) {
      mod1 <- untrusted_pairs$Model_1[l]
      mod2 <- untrusted_pairs$Model_2[l]
      sign_mat[mod1,mod2] = NA
    }
  }
  eco_mat <- as.data.frame(sign_mat)
  i <- 0
  for (sp1 in single_growth$Model.name) {
    i <- i + 1
    print(paste0("creating eco matrix, species Number: ", i))
    for (sp2 in single_growth$Model.name) {
      sign_sum <- sum(sign_mat[sp2, sp1], sign_mat[sp1, sp2])
      if (is.na(sign_sum)) {
        eco_mat[sp2, sp1] <- "Uncertain"
      } else if (sign_sum == 1) {
        eco_mat[sp2, sp1] <- "Commensalism"
      } else if (sign_sum == -1) {
        eco_mat[sp2, sp1] <- "Amensalism"
      } else if (sign_sum == 2) {
        eco_mat[sp2, sp1] <- "Mutualism"
      } else if (sign_sum == -2) {
        eco_mat[sp2, sp1] <- "Competition"
      } else if (sign_mat[sp2, sp1] == 0 & sign_mat[sp1, sp2] == 0) {
        eco_mat[sp2, sp1] <- "Neutralism"
      } else if (sign_mat[sp2, sp1] == 1 & sign_mat[sp1, sp2] == -1) {
        eco_mat[sp2, sp1] <- "Antagonism"
      } else if (sign_mat[sp2, sp1] == -1 & sign_mat[sp1, sp2] == 1) {
        eco_mat[sp2, sp1] <- "Antagonism"
      }
      else{
        eco_mat[sp2, sp1] <- "Wrong"
      }
    }
    eco_mat[sp1, sp1] <- "Alone"
  }
  out_mats <- list(eco_mat, change_mat)
  names(out_mats) <- c("eco_mat", "change_mat")
  return(out_mats)
}
