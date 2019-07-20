# May 21, 2019
#' Make matrix of contrast 
#' 
#' @description 
#' This function creates a matrix of contrasts between the individuals not in
#' the training set (i.e. calibration set) and the mean of the whole population.
#' 
#' @param not_phenotyped index of unphenotyped matrix in the relationship matrix.
#' 
#' @export
#'


make.contrast <- function(not_phenotyped) {
  mat=matrix(-1/n_total,n_total,n_total-n_selected)
  for (i in 1:ncol(mat)) {
    mat[not_phenotyped[i],i]=1-1/n_total
  }
  return(mat)
}# Close the function