#' June 10, 2019
#'
#' Get a random training set of same size as that selected by kinship threshold
#'
#' @description
#' Gets customized training set for each line in the validation set
#' selected randomly with its size equal to the size of training set
#' selected by kinship threshold.
#'
#'
#' @param train_by_kinship The dataframe of individuals selected by kinship threshold to be used for predicting a line
#' @param selected The dataframe of individuals selected by CDmean, serving as the training set.
#'
#' @import dplyr
#' @import data.table
#'
#' @return
#' A dataframe with the list of individuals selected randomly.
#'
#' @export
#'
get.train.by.random.sample <- function(selected, train_by_kinship){

  train_by_random_sample <- selected[sample(nrow(selected), nrow(train_by_kinship)),]
  train_by_random_sample <- as.data.frame(train_by_random_sample)
  colnames(train_by_random_sample) <- "ID"

  return(train_by_random_sample)
  #}

} # Close the function

