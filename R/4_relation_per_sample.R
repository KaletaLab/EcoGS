#' Frequency prediction of ecological relationships
#'
#' The prediction of the frequency of ecological relationship between pairs of co-grown bacteria in each sample (community)
#' as exported by the function make_eco_mat and a user imported OTU table
#'
#' @param OTU_table a matrix of the abundances of each species in each sample (columns are samples and rows are bacterial species)
#' @param eco_mat a list, as exported by the function make_eco_mat and which contains: matrix of all_models X all_models with each cell represents the ecological relation between the two strains
#' @return
#' (a list object in the R environment)
#'
#' a. matrix of the relative frequency for each ecological relation in each sample
#'
#' b. matrix of the relative frequency for each ecological relation in each sample weighted by the abundances from the OTU table
#'
#' @export
#'
relation_per_sample <- function(OTU_table, eco_mat){
  relations <- unique(unlist(eco_mat))
  relations <- relations[relations != "Alone"]
  relations_per_samp <- as.data.frame(matrix(0, ncol(OTU_table), length(relations)))
  names(relations_per_samp) <- relations
  rownames(relations_per_samp) <- names(OTU_table)
  weighed_relations <- relations_per_samp
  i <- 0
  id <- names(OTU_table)[1]
  for (id in names(OTU_table)) {
    i <- i + 1
    print(paste0("Sample Number: ", i))
    bac_ind <- OTU_table[,id] > 0
    if (sum(bac_ind) < 2) {
      next
    }
    bac_in_samp <- rownames(OTU_table)[bac_ind]
    samp_abund <- data.frame(species = bac_in_samp, abund = OTU_table[bac_in_samp,id])
    #  samp_abund <- samp_abund[order(samp_abund$abund, decreasing = T),]
    samp_abund <- samp_abund[order(samp_abund$abund),]
    rownames(samp_abund) <- samp_abund$species
    bac_in_samp[!bac_in_samp %in% names(eco_mat)]
    samp_eco_mat <- eco_mat[bac_in_samp, bac_in_samp]
    samp_abund_mat <- samp_eco_mat
    samp_abund_mat[] <- 0
    for (sp in rownames(samp_abund)) {
      samp_abund_mat[sp,] <- samp_abund[sp, "abund"]
      samp_abund_mat[,sp] <- samp_abund[sp, "abund"]
    }
    for (relat in relations) {
      ind <- samp_eco_mat == relat
      relations_per_samp[id, relat] = sum(ind)/2
      weighed_relations[id, relat] <- sum(samp_abund_mat[ind])/2
    }
  }
  out_mats <- list(relations_per_samp, weighed_relations)
  names(out_mats) <- c("relations_per_samp", "weighed_relations")
  return(out_mats)
}
