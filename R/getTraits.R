#' May 29, 2019
#'
#' Get train traits
#'
#' @description
#' Gets the each trait for training population
#'
#'
#' @param train_pheno A dataframe of phenotypic data for the training set
#'
#' @export
#'
get.train.traits <- function(train_pheno){

  StP_INC <- train_pheno$StP_INC
  StP_SEV <- train_pheno$StP_SEV
  StP_DIS <- train_pheno$StP_DIS
  StP_Microtwt <- train_pheno$StP_Microtwt
  StP_VSK <- train_pheno$StP_VSK

  Crk_INC <- train_pheno$Crk_INC
  Crk_SEV <- train_pheno$Crk_SEV
  Crk_DIS <- train_pheno$Crk_DIS
  Crk_Microtwt <- train_pheno$Crk_Microtwt
  Crk_VSK <- train_pheno$Crk_VSK

} # Close the function



#' May 29, 2019
#'
#' Get validation traits
#'
#' @description
#' Gets the each trait for the individual to be validated
#'
#'
#' @param validate_pheno A dataframe of phenotypic data for the training set
#'
#' @export
#'
get.validation.traits <- function(validate_pheno){

  StP_INC_vp = validate_pheno$StP_INC
  StP_SEV_vp = validate_pheno$StP_SEV
  StP_DIS_vp = validate_pheno$StP_DIS
  StP_Microtwt_vp = validate_pheno$StP_Microtwt
  StP_VSK_vp = validate_pheno$StP_VSK

  Crk_INC_vp = validate_pheno$Crk_INC
  Crk_SEV_vp = validate_pheno$Crk_SEV
  Crk_DIS_vp = validate_pheno$Crk_DIS
  Crk_Microtwt_vp = validate_pheno$Crk_Microtwt
  Crk_VSK_vp = validate_pheno$Crk_VSK

  return(StP_INC_vp, StP_SEV_vp, StP_DIS_vp, StP_Microtwt_vp, StP_VSK_vp,
             Crk_INC_vp, Crk_SEV_vp, Crk_DIS_vp, Crk_Microtwt_vp, Crk_VSK_vp)

} # Close the function
