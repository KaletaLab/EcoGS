#' Transformation for the ecological relationships into log10 of the pairwise ratios
#'
#'
#' @param relations_table the ecological relation frequencies per sample, exported by the function "relation_per_sample"
#' @return
#' 
#' @export
#'
relation_ratios <- function(relations_table){
  relations_table$Uncertain <- NULL
  column_combinations <- combn(colnames(relations_table), 2, simplify = TRUE)
  # Create a new dataframe with the result of the division for each pair
  relative_relations <- data.frame((apply(column_combinations, 2, function(cols) {
    relations_table[[cols[1]]] / relations_table[[cols[2]]]
  })))
  colnames(relative_relations) <- apply(column_combinations, 2, function(cols) {
    paste(cols, collapse = '/')
  })
  relative_relations <- log10(relative_relations)
  relative_relations[relative_relations=="NaN" | relative_relations =="Inf"] <- 0
  rownames(relative_relations) <- rownames(relations_table)
  return(relative_relations)
}