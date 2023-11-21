#' Plotting ecological relations
#'
#' To plot the number of pairs in each type of ecological relations as exported by the function make_eco_mat
#'
#' @param eco_mat A list, as exported by the function make_eco_mat. The matrix of all_models X all_models with each cell represents the ecological relation between the two strains.
#' @return A plot of the number of pairs in each type of ecological relations from the eco_mat
#'
#' @export
#'
plot_relations <- function(eco_mat){
  freq <- table(unlist(eco_mat))/2
  freq <- freq[names(freq) != "Alone"]
  freq <- sort(freq)
  colors <- data.frame(relation = c("Uncertain", "Mutualism", "Neutralism",  "Competition",
                                    "Commensalism",   "Antagonism",   "Amensalism"),
                       col = c("black", "orange", "yellow", "pink", "blue", "red", "darkgreen"))
  rownames(colors) <- colors$relation
  colors <- colors[names(freq),]
  myplot <- graphics::barplot(freq,
          main = "Eco relations",
          xlab = "Categories",
          ylab = "Frequency",
          col = colors$col)
  return(freq)
}
