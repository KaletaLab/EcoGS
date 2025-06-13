#' Metabolic interactions with MicrobiomeGS2
#'
#' Metabolic interactions and growth capacity of metabolic models
#' based on MicrobiomeGS2 and sybil packages
#'
#' @param list_of_models a NAMED list of models in sybil format on the same diet
#' @param cores the number of available cores
#' @param save_pair set T to save each pair of models on the working directory
#'
#' @return
#' (to be returned as an object in the R environment)
#'
#' a. the result of the biomass optimisation of each model alone
#'
#' b. the result of the biomass optimisation and the metabolic interactions for each pair of models
#'
#'=======================================
#'
#' (to be saved as a RDS file inside the working directory, if save_pair set to T)
#'
#' For each pair, a list with:
#'
#'    a.the joined model
#'
#'    b.the common exchange reactions
#'
#'    c.the fluxes of the meta-model
#'
#'    d.the metabolic interactions among the sub-models
#' @export



metabolic_interactions_with_MicrobiomeGS2 <- function(list_of_models,
                                             cores,save_pair) {
  #definitions
  `%dopar%` <- foreach::`%dopar%`
  i <- 1:length(list_of_models)
  #PART A  Growth of Each Model Alone
  #create clusters parallel computation
  cl <- parallel::makeCluster(cores, type = "PSOCK")
  parallel::clusterExport(cl, c("list_of_models"), envir = environment())
  doParallel::registerDoParallel(cl)

  # call sybil inside each cluster
  data_growth_alone <- foreach::foreach(i = 1:length(list_of_models),
                                        .combine = rbind) %dopar% {
                                 sybil::SYBIL_SETTINGS("SOLVER", "cplexAPI") #MicrobiomeGS Default
                                 sybil::SYBIL_SETTINGS("METHOD", "hybbaropt") #MicrobiomeGS Default

                                 optimum <- sybil::optimizeProb(list_of_models[[i]], algorithm = "fba")
                                 #Data frame to save the growth information of each model alone
                                 data.frame("Model name" = names(list_of_models)[i],
                                            "Biomass_rate_alone" = sybil::mod_obj(optimum))
                               }
  parallel::stopCluster(cl)

  #PART B  Growth of Pairs of  Models

  #get all possible pairs of two models accross all models
  total_combinations <- t(utils::combn(length(list_of_models),2))

  #create clusters parallel computation
  cl <- parallel::makeCluster(cores, type = "PSOCK")
  parallel::clusterExport(cl, c("list_of_models","total_combinations","save_pair"),
                          envir = environment())
  doParallel::registerDoParallel(cl)

  # call MicrobiomeGS2 and sybil inside each cluster
  data_growth_pairs <- foreach::foreach(i = 1:length(total_combinations[,1]),
                                        .combine = rbind) %dopar% {

                                sybil::SYBIL_SETTINGS("SOLVER", "cplexAPI") #MicrobiomeGS Default
                                sybil::SYBIL_SETTINGS("METHOD", "hybbaropt") #MicrobiomeGS Default

                                 # select combinations
                                 mycomb <- total_combinations[i,]

                                 # join models
                                 mod.joined <- MicrobiomeGS2::join_mult_models(list_of_models[mycomb], scale.boundaries = 1)

                                 # focus on the model
                                 mod_probio <- mod.joined$modj

                                 # remove objective coefficients if any
                                 mod_probio@obj_coef <- rep(0, length(mod_probio@react_id))

                                 # add objective coefficient to each biomass reaction
                                 mod_probio@obj_coef[grep("M[0-9]+_EX_cpd11416_c0",mod_probio@react_id)] <- 1

                                 # coupling
                                 cpl_probio <- MicrobiomeGS2::get_coupling_constraints_mult(mod_probio,
                                                                             cpl_c = 400,
                                                                             cpl_u = 1e-6)

                                 # warm FBA #MicrobiomeGS Default
                                 wrm_probio <- sybil::sysBiolAlg(mod_probio,
                                                          algorithm = "mtfEasyConstraint2",
                                                          easyConstraint = cpl_probio,
                                                          pFBAcoeff = 1e-6,
                                                          scaling = 1)
                                 # solution FBA
                                 sol_probio <- sybil::optimizeProb(wrm_probio)

                                 # biomass fluxes
                                 ub_lb <- sol_probio$fluxes[grep("cpd11416",mod_probio@react_id,value = F)]

                                 # did both model grow?
                                 both_grew <- F
                                 if (ub_lb[[1]] >  1E-6 & ub_lb[[2]] > 1E-6 ) {
                                   both_grew <- T
                                 }

                                 # extract reaction data
                                 # A data frame with all fluxes
                                 data_inside <- data.frame("REACTION_ID" = mod_probio@react_id,
                                                           "REACTION_NAME" = mod_probio@react_name,
                                                           "REACTION_MTF" = sol_probio$fluxes[1:length(mod_probio@react_id)])

                                 #a data frame with all fluxes of M1 model
                                 data_inside_M1 <- data_inside[grep("^M1_EX",data_inside$REACTION_ID),]
                                 # extracting and adding cpd ids to the dataframe
                                 cutm1 <- gsub(pattern = "M[0-9]_EX_cpd",replacement = "cpd",x = data_inside_M1$REACTION_ID)
                                 cutm3 <- unique(gsub(pattern = "_e0",replacement = "",x = cutm1,fixed = T))
                                 data_inside_M1$cpd <- cutm3
                                 # extracting positive and negative fluxes, if their absolute value > 1e-6
                                 data_inside_M1p <- data_inside_M1[which(data_inside_M1$REACTION_MTF > 0 & abs(data_inside_M1$REACTION_MTF) > 1e-6),]
                                 data_inside_M1n <- data_inside_M1[which(data_inside_M1$REACTION_MTF < 0 & abs(data_inside_M1$REACTION_MTF) > 1e-6),]

                                 #a data frame with all fluxes of M2 model
                                 data_inside_M2 <- data_inside[grep("^M2_EX",data_inside$REACTION_ID),]
                                 # extracting and adding cpd ids to the dataframe
                                 cutM2 <- gsub(pattern = "M[0-9]_EX_cpd",replacement = "cpd",x = data_inside_M2$REACTION_ID)
                                 cutm3 <- unique(gsub(pattern = "_e0",replacement = "",x = cutM2,fixed = T))
                                 data_inside_M2$cpd <- cutm3
                                 # extracting positive and negative fluxes, if their absolute value > 1e-6
                                 data_inside_M2p <- data_inside_M2[which(data_inside_M2$REACTION_MTF > 0 & abs(data_inside_M2$REACTION_MTF) > 1e-6),]
                                 data_inside_M2n <- data_inside_M2[which(data_inside_M2$REACTION_MTF < 0 & abs(data_inside_M2$REACTION_MTF) > 1e-6),]

                                 #combine M1 and M2 reactions that are of opposite direction (i.e., metabolic exchange=
                                 data_inside_M1p_M2n <- merge.data.frame(data_inside_M1p,data_inside_M2n,
                                                                         all = F, by = "cpd")
                                 data_inside_M2p_M1n <- merge.data.frame(data_inside_M2p,data_inside_M1n,
                                                                         all = F, by = "cpd")

                                 if (save_pair == T) {
                                   #a list to save the growth information of each pair locally
                                   list_internal <- vector(mode = "list", length = 5)
                                   list_internal[[1]] <- mod.joined$model.IDs
                                   names(list_internal)[1] <- "joined_model"
                                   list_internal[[2]] <- mod.joined$ex.rxns
                                   names(list_internal)[2] <- "common_exchange_reactions"
                                   list_internal[[3]] <- data_inside
                                   names(list_internal)[3] <- "all_fluxes"
                                   list_internal[[4]] <- data_inside_M1p_M2n
                                   names(list_internal)[4] <- "fluxes_from_M1_to_M2"
                                   list_internal[[5]] <- data_inside_M2p_M1n
                                   names(list_internal)[5] <- "fluxes_from_M2_to_M1"
                                   saveRDS(list_internal, paste0(paste0(names(list_of_models[mycomb]),
                                                                        collapse = "_+_"), ".RDS"))

                                 }

                                 #a data frame to save the growth information of each pair
                                 data.frame("Model_1" = names(list_of_models[mycomb][1]),
                                            "Model_2" = names(list_of_models[mycomb][2]),
                                            "Model_1_Growth" = ub_lb[1],
                                            "Model_2_Growth" = ub_lb[2],
                                            "Interaction_from_M1_to_M2_IDs" = paste0(data_inside_M1p_M2n$cpd, collapse = ", "),
                                            "Interaction_from_M1_to_M2_Names" = paste0(data_inside_M1p_M2n$REACTION_NAME.x, collapse = ", "),
                                            "Interaction_from_M2_to_M1_IDs" = paste0(data_inside_M2p_M1n$cpd, collapse = ", "),
                                            "Interaction_from_M2_to_M1_Names" = paste0(data_inside_M2p_M1n$REACTION_NAME.x, collapse = ", "),
                                            "Optimisation_status" = sol_probio$stat,
                                            "Both_Grew" = both_grew)
                               }
  parallel::stopCluster(cl)
  #finally, the list to be returned....
  list_export <- vector("list",length = 2)
  list_export[[1]] <- data_growth_pairs
  names(list_export)[[1]] <- "growth_in_pairs"
  list_export[[2]] <- data_growth_alone
  names(list_export)[[2]] <- "single_growth"
  return(list_export)
}
