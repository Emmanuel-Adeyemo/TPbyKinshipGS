# July 15, 2019
#' Calculate lambda from REML
#'
#' @description
#' Uses a RR-BLUP mixed model to estimate the marker variance (Vm) and the residual variance
#' (Ve) to be used to calculate lambda
#'
#' @references Endelman, J. B. 2011. Ridge Regression and Other Kernels for Genomic Selection
#' with R Package rrBLUP. Plant Genome 4:250-255. doi:10.3835/plantgenome2011.08.0024
#'
#' @param pheno_train_trait a matrix of phenotypes for trait of interest in training population.
#' @param geno_train An incidence matrix of genotypes for the training population.
#'
#' @return
#' lambda
#'
#' @import rrBLUP
#'
#' @export
#'
calculate.lambda <- function(pheno_train_trait, geno_train, model = "RRBLUP") {

  # Deal with input
  rownames(geno_train) <- geno_train[,1]
  geno_train <- geno_train[,-1]

  pheno_train_trait <- as.matrix(pheno_train_trait)
  geno_train <- as.matrix(geno_train)

  # Solve the mixed model
  solve_out <- mixed.solve(y = pheno_train_trait, Z = geno_train, method = "REML")

  V_m <- solve_out$Vu  # marker variance
  V_e <- solve_out$Ve  # residual varaince
  N_m <- ncol(geno_train)  # number of markers

  V_a <- V_m * N_m  # additive variance
  lambda <- V_e/V_a

  # Return the data
  return(lambda)

} # Close the function



# July 15, 2019
#' Calculate heritability from REML
#'
#' @description
#' Uses a RR-BLUP mixed model to estimate the marker variance (Vm) and the residual variance
#' (Ve) to be used to calculate heritability
#'
#' @references Endelman, J. B. 2011. Ridge Regression and Other Kernels for Genomic Selection
#' with R Package rrBLUP. Plant Genome 4:250-255. doi:10.3835/plantgenome2011.08.0024
#'
#' @param pheno_train_trait a matrix of phenotypes for trait of interest in training population.
#' @param geno_train An incidence matrix of genotypes for the training population.
#'
#' @return
#' lambda
#'
#' @import rrBLUP
#'
#' @export
#'
calculate.heritability <- function(pheno_train_trait, geno_train, model = "RRBLUP") {

  # Deal with input
  rownames(geno_train) <- geno_train[,1]
  geno_train <- geno_train[,-1]

  pheno_train_trait <- as.matrix(pheno_train_trait)
  geno_train <- as.matrix(geno_train)

  # Solve the mixed model
  solve_out <- mixed.solve(y = pheno_train_trait, Z = geno_train, method = "REML")

  V_m <- solve_out$Vu  # marker variance
  V_e <- solve_out$Ve  # residual varaince
  N_m <- ncol(geno_train)  # number of markers

  V_a <- V_m * N_m  # additive variance
  h2 <- (2*V_a)/ ((2*V_a) + V_e)

  # Return the data
  return(h2)

} # Close the function
