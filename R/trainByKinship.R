# May 23, 2019
#'
#' Get training set based on kinship threshold set
#'
#' @description
#' Gets customized training set for each line in the validation set
#' based on the kinship threshold set i.e. a line in the validation
#' set would only be predicted with a model trained with similar lines
#' only.
#'
#'
#' @param line_id A vector list of line IDs of individuals not selected by CDmean algorithm (i.e. validation population)
#' @param ki_train A dataframe of kinship values of individuals selected by the CDmean algorithm
#' @param threshold The kinship threshold used to filter lines to include in the training model. Should be provided by user
#'
#' @import dplyr
#' @import data.table
#'
#' @return
#' A dataframe with the list of individuals selected by the threshold
#'
#' @export
#'
get.train.by.kinship <- function(line_id, ki_train, threshold){


  for(i in line_id){
    train_by_kinship <- ki_train %>% dplyr::select(ID = ID, i) %>% filter(ki_train[i] >= threshold) %>% dplyr::select(-i)  # selects only column
    # that == i(line ID) and ID(contains rownames)
    # filters out individuals less than the threshold provided
    # lastly, takes out the column i, only id is left
    return(train_by_kinship)
  }

} # Close the function


