# May 21, 2019
#' Estimate additive relationship matrix
#' 
#' @description 
#' Uses A.mat function in RR-BLUP to estimate additive relationship matrix of  
#' genotypic. This is the first step for the CDmean algorithm. 
#' 
#' @param pheno A dataframe of phenotypes for the whole population.
#' @param geno An incidence matrix of genotypes for the whole population.
#' 
#' 
#' @import rrBLUP
#' @import MASS
#' 
#' @export
#' 
make.Amatrix <- function(pheno, geno){
                            
  
  # Deal with input
  linename <- as.data.frame(pheno[,1])
  colnames(linename)[1] <- "ID"
  
  rownames(geno) = geno[,1]
  geno = geno[,-1]
  geno = as.matrix(geno)
  
  A_matrix <- A.mat(geno)
  row.names(A_matrix) <- NULL
  A_matrix <- as.matrix(A_matrix)
  
  # Return the data
  return(A_matrix)
  
  
} # Close the function