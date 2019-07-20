# July 18, 2019
#' Calculate predictive ability
#'
#' @description
#' Uses a RR-BLUP mixed model to estimate GEBV of lines
#'
#' @references Endelman, J. B. 2011. Ridge Regression and Other Kernels for Genomic Selection
#' with R Package rrBLUP. Plant Genome 4:250-255. doi:10.3835/plantgenome2011.08.0024
#'
#' @param validate_id a list of line names thet were not selected by CD algorithm (VP)
#' @param train_kinship A matrix of kinship relationship among lines in the training population.
#' @param selected The number of lines (TP) selected from the CD algorithm
#' @param pheno A dataframe of lines with their phenotypic data
#' @param geno a matrix of genotypic data for lines in training population.
#' @param trait The trait of interest. Should be written as a character in quotation marks.
#' @param threshold The least amount of similarity the lines used for TP should have.
#'
#' @return
#' outer.list
#'
#' @import rrBLUP
#'
#' @export
#'
calculate.prediction.ability <- function(validate_id, 
                                         train_kinship, 
                                         selected, 
                                         pheno, 
                                         geno,
                                         trait,
                                         threshold) {
  
  # Prediction with Training population selected by kinship
  
  # make an empty dataframe
  predictions_kinship <- data.frame()
  
  predictions_random <- data.frame()
  
  training_size <- data.frame()
  
  for (val_line in validate_id){
    
    train_kinship_selected = get.train.by.kinship(val_line, train_kinship, threshold) # selects training pop
    # # by kinship used to train model # 0.00 is the threshold provided
    
    train_random_selected = get.train.by.random.sample(selected, train_kinship_selected)
    
    train_random_pheno = do.innerjoin(train_random_selected, pheno) # gets the phenotypic data for training set selected by kinship
    train_random_geno = do.innerjoin(train_random_selected, geno)
    
    train_pheno = do.innerjoin(train_kinship_selected, pheno) # gets the phenotypic data for training set selected by kinship
    train_geno = do.innerjoin(train_kinship_selected, geno) # gets genotypic data for training set selected by kinship
    
    validate_pheno = get.validation.data(pheno, val_line) # gets phenotypic data for line being predicted
    validate_geno = get.validation.data(geno, val_line) # gets genotypic data for line being predicted
    
    predicted_value = make.predictions(train_pheno,train_pheno[,trait], train_geno, validate_geno)
    observed_value = validate_pheno[,trait]
    
    combined <- data.frame(predicted = predicted_value, 
                           observed = observed_value)
    
    predicted_random_value = make.predictions(train_random_pheno,train_random_pheno[,trait], train_random_geno, validate_geno)
    
    
    combined_random <- data.frame(predicted_random = predicted_random_value,
                                  observed = observed_value)
    
    predictions_random <- rbind(predictions_random, combined_random)
    
    predictions_kinship <- rbind(predictions_kinship, combined)
    
    kinship_random_training_size <- get.training.size(train_kinship_selected, train_random_selected)
    
    training_size <- rbind(training_size, kinship_random_training_size)
    
    
    ##########################################################################################################
    # Prediction with Training population selected by CDmean alone i.e. no kinship model
    
    train_pheno = do.innerjoin(selected, pheno) # gets the phenotypic data for training set selected by kinship
    train_geno = do.innerjoin(selected, geno) # gets genotypic data for training set selected by kinship
    
    validate_pheno = do.antijoin(selected, pheno) # gets phenotypic data for line being predicted
    validate_geno = do.antijoin(selected, geno) # gets genotypic data for line being predicted
    
    no_kinship = make.predictions.no.kinship(train_pheno$DIS, train_geno, validate_pheno$DIS, validate_geno)
    
    out.list <- list(predictions_random,predictions_kinship, kinship_random_training_size, training_size, no_kinship )
    
  }
  
  random <- as.data.frame(out.list[[1]])
  kinship <- as.data.frame(out.list[[2]])
  training_size <- as.data.frame(out.list[[3]])
  no_kinship <- as.data.frame(out.list[[5]])

  random_pred_ability <- cor(random$predicted_random, random$observed, use = "complete.obs") # finds correlation i.e. prediction accuracy
  kinship_pred_ability <- cor(kinship$predicted, kinship$observed, use = "complete.obs") # finds correlation i.e. prediction accuracy
  training_size <- training_size$kinship_size
  kinship_pred_ability <- cor(kinship$predicted, kinship$observed, use = "complete.obs")

  outer.list <- list(kinship_pred_ability, random_pred_ability, training_size, no_kinship)
  
  return(outer.list)
} # Close the function



