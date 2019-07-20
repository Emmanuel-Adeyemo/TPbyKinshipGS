#' May 21, 2019
#' Edited June 7, 2019
#' Make genomic predictions
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
#' @param pheno_train A matrix of phenotypes for the training population.
#' @param pheno_train_trait a matrix of phenotypes for trait of interest in training population.
#' @param geno_train An incidence matrix of genotypes for the training population.
#' @param geno_validate An incidence matrix of genotypes for the prediction population
#' (i.e. selection candidates).
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
make.predictions <- function(pheno_train, pheno_train_trait, geno_train, geno_validate, model = "RRBLUP") {

  # need a TP size of at least 5 individuals
  if (nrow(pheno_train) < 5){
    GEBV <- "NA"
    GEBV <- as.numeric(GEBV)
  }else{
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
    #pred_accuracy <- cor(GEBV, pheno_validate, use = "complete.obs")
  }

  # Return the data
  return(GEBV)

} # Close the function
