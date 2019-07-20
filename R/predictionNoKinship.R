# May 21, 2019
#' Make genomic predictions with no kinship model
#'
#' @description
#' Uses a RR-BLUP mixed model to predict marker effects given phenotypic and
#' genotypic data from a training population. Then uses the marker effects
#' and genotypic data from selection candidates to predict genotypic values for
#' those candidates. The REML method of estimating variances is used.
#'
#' @references Endelman, J. B. 2011. Ridge Regression and Other Kernels for Genomic Selection
#' with R Package rrBLUP. Plant Genome 4:250-255. doi:10.3835/plantgenome2011.08.0024
#'
#' @param pheno_train_trait a matrix of phenotypes for trait of interest in training population.
#' @param geno_train An incidence matrix of genotypes for the training population.
#' @param geno_validate An incidence matrix of genotypes for the prediction population
#' (i.e. selection candidates).
#' @param pheno_validate A dataframe of observed phenotypes for the validation population
#' to be used for calculating predicting accuracy of the trained model
#'
#' @return
#' The prediction accuracy of the trained model
#'
#' @import dplyr
#' @import rrBLUP
#'
#'
#' @export
#'
make.predictions.no.kinship <- function(pheno_train_trait, geno_train, pheno_validate, geno_validate, model = "RRBLUP") {
  
  # Deal with input
  pheno_train_trait <- as.matrix(pheno_train_trait)
  geno_train <- as.matrix(geno_train)
  geno_validate <- as.matrix(geno_validate)
  
  # Solve the mixed model
  solve_out <- mixed.solve(y = pheno_train_trait, Z = geno_train, method = "REML")
  marker_effects <- solve_out$u
  
  # Calculate GEBVs
  GEBV <- geno_validate %*% marker_effects
  
  row.names(GEBV) <- row.names(geno_validate)
  
  # find correlations between the predicted and observed values in the validation population
  pred_accuracy <- cor(GEBV, pheno_validate, use = "complete.obs")
  
  
  # Return the data
  return(pred_accuracy)
  
} # Close the function
