# May 21, 2019
#' Hidden functions
#'
#' Makes inverse additive matrix
#'
#' @description
#' Makes an inverse matrix of additive relationship in the whole population
#'
#' @param A_matrix A matrix of additive relationship among individuals in the whole population
#'
#' @return
#' A inverse matrix of the additive relationship
#'
#' @import MASS
#' @import rrBLUP
#'
#' @export
#'
make.Ainverse <- function(A_matrix){

  invA1=ginv(A_matrix) # Inverse of the covariance matrix
  # Return the data
  return(invA1)

} # Close the function




#'May 21, 2019
#'
#'   Assigns total number of lines
#'
#' @description
#' Assigns the total number of lines in the population to a variable - n_total
#'
#' @param A_matrix A matrix of additive relationship among individuals in the whole population
#'
#' @export
#'
assign.popsize <- function(A_matrix){

  n_total=nrow(A_matrix) # Inverse of the covariance matrix
  # Return the data
  return(n_total)

} # Close the function




#' May 21, 2019
#'
#' Do inner join
#'
#' @description
#' Does an inner join to extract selected lines with CDmean algorithm would serve as
#' TP for that cluster
#'
#' @param selected A dataframe of individuals selected by CDmean algorithm
#' @param pop A dataframe of phenotype or genotype of all individuals in the population
#'
#' @import dplyr
#' @import data.table
#'
#' @export
#'
do.innerjoin <- function(selected, pop){

  if(ncol(pop) < 20){ # checks to make sure pop is phenotypic data of the entire population
    k_train <- inner_join(selected, pop, by = "ID")
    # Return the data
    return(k_train)
  }
  else if(ncol(pop) > 20 && ncol(pop) < 1000){ #checks to make sure pop is kinship of indiviuals
                                                #(kin data column is not more than 500)
    k_train <- inner_join(selected, pop, by = "ID")
    return(k_train) # returns kinship matrix of individuals in the training pop
  }
  else if(ncol(pop) > 1000){ # checks to make sure pop is genotypic data of the entire population
    k_train <- inner_join(selected, pop, by = "ID")
    rownames(k_train) <- k_train[,1]
    k_train <- k_train[,-1]
    k_train <- as.matrix(k_train)
    return(k_train)
  }

} # Close the function




#' May 21, 2019
#'
#' Do anti join
#'
#' @description
#' Does an anti join to extract selected lines with CDmean algorithm would serve as
#' VP for that cluster
#'
#' @param selected A dataframe of individuals selected by CDmean algorithm
#' @param pop A dataframe of phenotype or genotype of all individuals in the
#' population
#'
#' @import dplyr
#' @import data.table
#'
#' @export
#'
do.antijoin <- function(selected, pop){

  if(ncol(pop) < 20){ # checks to make sure pop is phenotypic data of the entire population (pheno data has columns < 20)
    k_validate <- anti_join(pop, selected, by = "ID")
    # Return the data
    return(k_validate)
  }
  else if(ncol(pop) > 20 && ncol(pop) < 1000){ #checks to make sure pop is kinship of indiviuals
                                              #(kin data column is not more than 500)
    k_validate <- anti_join(pop, selected, by = "ID")
    return(k_validate) # returns kinship matrix of individuals in the validation pop
  }

  else if(ncol(pop) > 1000){ # checks to make sure pop is genotypic data of the entire population (geno data has columns > 1000)
    k_validate <- anti_join(pop, selected, by = "ID")
    rownames(k_validate) <- k_validate[,1]
    k_validate <- k_validate[,-1]
    k_validate <- as.matrix(k_validate)
  }


} # Close the function




#' May 22, 2019
#'
#' Get line id
#'
#' @description
#' Gets the line id of the individuals in the validation population which
#' would later be used to select the individuals to include in the training
#' model (i.e. training model would be built by kinship)
#'
#'
#' @param selected A dataframe of kinship values of individuals not selected by the CDmean algorithm (nvalidationpop x npop dataframe)
#'
#' @import dplyr
#' @import data.table
#'
#' @export
#'
get.lineID <- function(not_selected_kinship){

  transpose_not_selected_kinship <- t(not_selected_kinship) # this is to make the notselected lines as column names (i.e. npop x nvalidationpop)
  transpose_not_selected_kinship <- as.data.frame(transpose_not_selected_kinship)
  colnames(transpose_not_selected_kinship) <- (unlist(transpose_not_selected_kinship[1,])) # the first row will be the header
  transpose_not_selected_kinship <- transpose_not_selected_kinship[-1, ]  #removes first row

  line_id <- colnames(transpose_not_selected_kinship)

  # for (line in line_id){
  #   line_df <- data.frame(ID = line)
  # }
  return(line_id) # returns each line in the validation population as a vector

} # Close the function



#' May 23, 2019
#' Edit: May 24, 2019
#'
#' Get the validation data
#'
#' @description
#' Gets the phenotypic or genotypic data for each validation line
#'
#'
#' @param data A dataframe of individuals with phenotypic or genotypic data
#' @param line_id A vector list of line IDs of individuals not selected by
#' CDmean algorithm
#'
#' @export
#'
get.validation.data <- function(data_field, line_id){

  if(ncol(data_field) < 20){
    for(line in line_id){
      k_validate_data <- subset(data_field, data_field$ID == line)

    }

    return(k_validate_data)
  }
  else if(ncol(data_field) > 20){
    for(line in line_id){
      k_validate_data <- subset(data_field, data_field$ID == line)
      rownames(k_validate_data) <- k_validate_data[,1]
      k_validate_data <- k_validate_data[,-1]
      k_validate_data <- as.matrix(k_validate_data)
    }

    return(k_validate_data)
  }


} # Close the function



#' May 24, 2019
#'
#' Custom anti join to get validation pop
#'
#' @description
#' Does an anti join to get validation population and returns only the
#' line id
#'
#' @param selected A dataframe of individuals selected by CDmean algorithm
#' @param cluster_pop A dataframe of phenotype or genotype of all individuals
#' a specific cluster and not the whole population
#'
#' @import dplyr
#' @import data.table
#'
#' @export
#'
get.not.selected <- function(selected, cluster_pop){

  k_not_selected <- anti_join(cluster_pop, selected, by = "ID") %>% dplyr::select(ID)

  # Return the data
  return(k_not_selected)

} # Close the function



# May 29, 2019
#'
#' Removes extra rownames autogenerated in R
#'
#' @description
#' Removes the "annoying" extra rownames autogenerated in R when importing data.
#' This is especially useful when creating a matrix from the dataframe.
#'
#' @param X A dataframe of individuals with extra column to the left
#'
#' @return
#' A matrix without the "annoying" autogenerated column numbers
#'
#'
#' @export
#'
rm.col <- function(X){

  rownames(X) <- X[,1]
  X <- X[,-1]
  X <- as.matrix(X)

  return(X)

} # Close the function



#' May 29, 2019
#'
#' Custom-custom anti join to get validation pop
#'
#' @description
#' Does an anti join to get validation population and returns the
#' line id and phenotypic data of traits
#'
#' @param selected A dataframe of individuals selected by CDmean algorithm
#' @param cluster_pop A dataframe of phenotype or genotype of all individuals
#' a specific cluster and not the whole population
#'
#' @import dplyr
#' @import data.table
#'
#' @export
#'
get.not.selected.traits <- function(selected, cluster_pop){

  k_not_selected_traits <- anti_join(cluster_pop, selected, by = "ID")

  # Return the data
  return(k_not_selected_traits)

} # Close the function



#' May 31, 2019
#'
#' Calculates the mean prediction accuracy
#'
#' @description
#' Calculates the mean prediction accuracy from the fifty iterations
#'
#'
#' @import dplyr
#' @import data.table
#'
#' @export
#'

find.mean.prediction <- function(){

  files <- dir(pattern = "*.txt")
  combined_predictions <- lapply(files, read.delim)
  combined_predictions <- as.data.frame(combined_predictions)
  mean_predictions <- rowMeans(combined_predictions)
  mean_predictions <- as.data.frame(mean_predictions)
  return(mean_predictions)

}


#' June 9, 2019
#'
#' makes plot for predictions
#'
#' @description
#' makes plot for predictions with ggplot
#'
#' @param file The file containing information to be ploted
#'
#' @import ggplot
#' @import reshape2
#' @import plyr
#' @import ggfortify
#' @import tidyr
#' @import scales
#' @import gridExtra
#' @import extrafont
#'
#' @export
#'

make.plot <- function(file){

  file_long <- gather(file, "Trait", "Prediction_Ability",
                            Inc, Sev, Dis,  Microtwt, VSK)

  file_long$Trait <- factor(file_long$Trait,levels = c("Inc", "Sev", "Dis",
                                                                   "Microtwt", "VSK"))
  k_plot <- ggplot(file_long, aes(x=Trait, y=Prediction_Ability, colour=threshold, group=threshold)) +
    geom_line(size=1.3,linetype = "dashed") + xlab("") + ylab(" ") + geom_point(size=2.3)

  k_plot_bw <- k_plot + theme_bw() +
    theme(axis.text.x = element_text(size=18),
          axis.text.y = element_text(size=18),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text=element_text(family="Candara", size=15,color = "grey20"))+
    theme(legend.justification = c(1, 0), legend.position = c(1, 0),
          legend.box.margin=margin(c(50,5,20,50)),legend.background=element_blank())


  return(k_plot_bw)

}



#' June 11, 2019
#'
#' Gets the training size
#'
#' @description
#' Gets the training size for each line in the validation set
#'
#' @param train_kinship_selected A dataframe of lines selected by kinship threshold to be used as training set
#' @param train_random_selected A dataframe of lines selected by random with size equal to that selected by kinship threshold
#'
#'
#' @export
#'

get.training.size <- function(train_kinship_selected, train_random_selected){

  kinship_size <- nrow(train_kinship_selected)
  random_size <- nrow(train_random_selected)

  combined_size <- data.frame(kinship_size, random_size)
  return(combined_size)

}


#' June 11, 2019
#'
#' Calculates the mean training size
#'
#' @description
#' Calculates the mean training from the fifty iterations
#'
#'
#' @import dplyr
#' @import data.table
#'
#' @export
#'

find.mean.training_size <- function(){

  files <- dir(pattern = "*.txt")
  combined_training_size <- lapply(files, read.delim)
  combined_training_size <- as.data.frame(combined_training_size)
  mean_training_size <- colMeans(combined_training_size)
  mean_training_size <- as.data.frame(mean_training_size)

  mean_training_size_all <- colMeans(mean_training_size)
  mean_training_size_all <- as.data.frame(mean_training_size_all)
  return(mean_training_size_all)

}



#' June 13, 2019
#'
#' Takes out rows with NA
#'
#' @description
#' Takes out rows having NA cells
#'
#' @param dataFile A dataframe with missing values denoted as NAs
#'
#' @import tidyr
#' @import dplyr
#'
#' @export
#'

remove.na <- function(dataFile){

  clean_data <- dataFile %>% drop_na()

  rownames(clean_data) <- clean_data[,1]
  clean_data <- clean_data[,-1]

  return(clean_data)

}


#' June 12, 2019
#' Creates genetic distance
#'
#' @description
#' Creates genetic distance based on euclidian distance
#'
#' @param geno SNPs file
#'
#'
#' @export
#'

get.genetic.distance <- function(geno){

  geno <- as.matrix(geno)
  geno_distance <- dist(geno)

  rownames(geno_distance) <- geno_distance[,1]
  geno_distance <- geno_distance[,-1]

  return(geno_distance)
}


#' June 14, 2019
#'
#' Inner join without futher processing
#'
#' @description
#' Does an inner join without further processing in previous
#'
#' @param selected A dataframe of individuals selected by CDmean algorithm
#' @param cluster_pop A dataframe of phenotype or genotype of all individuals
#' a specific cluster and not the whole population
#'
#' @import dplyr
#' @import data.table
#'
#' @export
#'
do.innerjoin.only <- function(selected, cluster_pop){

  innerjoin_only <- inner_join(selected, cluster_pop, by = "ID")

  # Return the data
  return(innerjoin_only)

} # Close the function



#' June 20, 2019
#'
#' Test for significance between two correlations
#'
#' @description
#' Two tail test for significance between two correlations
#'
#' @param r1 Correlation coefficient for group one
#' @param r2 Correlation coefficient for group two
#' @param n1 Sample size for group one
#' @param n2 Sample size for group two
#'
#'
#' @export
#'
diff.corr <- function( r1, n1, r2, n2 ){

  Z1 <- 0.5 * log( (1+r1)/(1-r1) )
  Z2 <- 0.5 * log( (1+r2)/(1-r2) )

  diff   <- Z1 - Z2
  SEdiff <- sqrt( 1/(n1 - 3) + 1/(n2 - 3) )
  diff.Z  <- diff/SEdiff

  p <- 2*pnorm( abs(diff.Z), lower=F)

  return(p)

} # close the function


#' June 9, 2019
#'
#' makes plot for predictions
#'
#' @description
#' makes plot for predictions with ggplot
#'
#' @param file The file containing information to be ploted
#'
#' @import ggplot
#' @import reshape2
#' @import plyr
#' @import ggfortify
#' @import tidyr
#' @import scales
#' @import gridExtra
#' @import extrafont
#'
#' @export
#'

make.plot.pa <- function(file){

  file_long <- gather(file, "Trait", "Prediction_Ability",
                      Dis,  Microtwt, VSK)

  file_long$Trait <- factor(file_long$Trait,levels = c( "Dis","VSK",
                                                        "Microtwt"))
  # k_plot <- ggplot(file_long, aes(x=Trait, y=Prediction_Ability, colour=threshold, group=threshold)) +
  #   geom_line(size=1.3,linetype = "dashed") + xlab("") + ylab(" ") + geom_point(size=2.3)
  #
  # k_plot_bw <- k_plot + theme_bw() +
  #   theme(axis.text.x = element_text(size=18),
  #         axis.text.y = element_text(size=18),
  #         axis.title.x = element_text(size=16),
  #         axis.title.y = element_text(size=16),
  #         text=element_text(family="Candara", size=15,color = "grey20"))#+
    #theme(legend.justification = c(1, 0), legend.position = c(1, 0),
          #legend.box.margin=margin(c(50,5,20,50)),legend.background=element_blank())



  k_plot_bw <- ggline(file_long, x=Trait, y=Prediction_Ability, add = "mean_se",
         color = "threshold", palette = "jco")+
    stat_compare_means(aes(group = threshold), label = "p.signif"
                       )

  return(k_plot_bw)

}

