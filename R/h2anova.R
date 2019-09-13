# July 19, 2019
#' Calculate narrow sense heritability
#'
#' @description
#' Calculates narroesense heritability from Analysis of Variance
#'
#'
#' @param pheno A dataframe of lines with their phenotypic data
#' @param r Number of reps in the study
#' @param e Number of environments in the study
#' @param trait The trait of interest. Should be written as a character in quotation marks.
#'
#' @return
#'
#' @export
#'
calculate.heritability.anova <- function(trait, pheno, r, e){

  out <- anova(lm(pheno[,trait] ~ Rep + ID + ID * Rep, data=pheno))
  msp <- out[2,3]
  mspe <- out[3,3]

  V_g <- (msp - mspe)/(r*e)
  V_y <- mspe/(r*e)
  V_e <- out[4,3]
  h2 <- V_g/(V_g + V_y)

  #lambda <- V_e/V_g

  return(h2)

}
