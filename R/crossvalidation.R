#' July 26, 2019
#' Make five-folds cross validation with RRBLUP, BB and RHKS
#'
#' @description
#' Uses both parametric and semi-parametric models to make predictions
#'
#' @references Endelman, J. B. 2011. Ridge Regression and Other Kernels for Genomic Selection
#' with R Package rrBLUP. Plant Genome 4:250-255. doi:10.3835/plantgenome2011.08.0024
#' @references P?rez P, de los Campos G. Genome-wide regression and prediction with the BGLR
#' statistical package. Genetics. 2014;198(2):483-495. doi:10.1534/genetics.114.164442
#'
#' @param pheno_train_trait a matrix of phenotypes for trait of interest in training population.
#' @param geno_train An incidence matrix of genotypes for the training population.
#' @param model Default is all models which uses RRBLUP, BB and RKHS. You can also impute
#' individual models
#'
#' @return
#' The prediction accuracy of the trained model
#'
#' @import dplyr
#' @import rrBLUP
#' @import BGLR
#'
#' @export
#'


cross.validation.models <- function(pheno_train_trait, geno_train, model = "All_models") {


  # Deal with input
  pheno_train_trait <- as.matrix(pheno_train_trait)
  geno_train <- as.matrix(geno_train)
  phenoTrait <- pheno_train_trait
  nTST <- nTST; n <- nrow(geno_train); p <-ncol(geno_train)
  nRep <- nRep; nIter<-nIter; burnIn<-burnIn
  G <- tcrossprod(geno)/p

  if (model == "RRBLUP"){

    rownames(geno_train) <- 1:nrow(geno_train) #give "gid" to marker data
    genoD <- as.matrix(dist(geno_train)) #calculate genetic distance

    phenoTrait <- pheno_train_trait

    traits <- traits
    cycles <- cycles                                                                                    #----cycle number of CV
    pred_ability <- matrix(nrow=cycles, ncol=traits)
    validNum <- nrow(geno_train) %*% test_portion
    for (r in 1:cycles){

      valid <- as.matrix(sample(1:nrow(geno_train),validNum, replace=FALSE))
      phenoTraitNA <- phenoTrait
      phenoTraitNA[valid] <- NA
      data1 <- data.frame(y=phenoTraitNA, gid=1:nrow(geno_train))
      solve.model <- kin.blup(data1, geno="gid", pheno="y", K=genoD, GAUSS=TRUE)
      pred_ability[r,1] <- cor(as.vector(solve.model$g[valid]), as.numeric(phenoTrait[valid]), use="complete")          #accuracy

    }
    return(pred_ability)
  }


  else if (model == "BB") {

    phenoTrait <- pheno_train_trait                    ##is vector
    y <- scale(phenoTrait, center=TRUE, scale=TRUE)  ## is matrix

    ETA<-list(MRK=list(X=geno_train,model="BayesB", probIn=probIn))

    pred_ability <- matrix(nrow=nRep, ncol=1, NA)
    for(i in 1:nRep){
      tst <- sample(1:n, size=nTST, replace=FALSE)
      yNA <- y; yNA[tst] <- NA
      fmBB <- BGLR(y=yNA,ETA=ETA, nIter=nIter, burnIn=burnIn,saveAt="BB_")
      pred_ability[i,1] <- cor(y[tst],fmBB$yHat[tst], use="complete")
    }

    return(pred_ability)
  }

  else if (model == "RHKS") {

    phenoTrait <- pheno_train_trait# is vector
    y <- scale(phenoTrait, center=TRUE, scale=TRUE) # is matrix
   
    ETA<-list(MRK=list(K=G,model="RKHS"))

    pred_ability <- matrix(nrow=nRep, ncol=1, NA)
    for(i in 1:nRep){
      tst <- sample(1:n, size=nTST, replace=FALSE)
      yNA <- y; yNA[tst]=NA
      fmRKHS <-BGLR(y=yNA,ETA=ETA, nIter=nIter, burnIn=burnIn,saveAt="RKHS_")
      pred_ability[i,1] <- cor(y[tst],fmRKHS$yHat[tst], use="complete")
    }

    return(pred_ability)
  }

  else if (model == "All_models") {

    rownames(geno_train) <- 1:nrow(geno_train) #give "gid" to marker data
    genoD <- as.matrix(dist(geno_train)) #calculate D # runs forever!

    phenoTrait <- pheno_train_trait

    traits <- traits
    cycles <- cycles                                                                                    #----cycle number of CV
    pred_ability_rrblup <- matrix(nrow=cycles, ncol=traits)
    validNum <- nrow(geno_train) %*% test_portion
    for (r in 1:cycles){

      valid <- as.matrix(sample(1:nrow(geno_train),validNum, replace=FALSE))
      phenoTraitNA <- phenoTrait
      phenoTraitNA[valid] <- NA
      data1 <- data.frame(y=phenoTraitNA, gid=1:nrow(geno_train))
      ans1 <- kin.blup(data1, geno="gid", pheno="y", K=genoD, GAUSS=TRUE)
      pred_ability_rrblup[r,1] <- cor(as.vector(ans1$g[valid]), as.numeric(phenoTrait[valid]), use="complete")          #accuracy

    }

    # BB
    y <- scale(phenoTrait, center=TRUE, scale=TRUE)  ## is matrix

    ETA<-list(MRK=list(X=geno_train,model="BayesB", probIn=probIn))

    pred_ability_bb <- matrix(nrow=nRep, ncol=1, NA)
    for(i in 1:nRep){
      tst <- sample(1:n, size=nTST, replace=FALSE)
      yNA <- y; yNA[tst] <- NA
      fmBB <- BGLR(y=yNA,ETA=ETA, nIter=nIter, burnIn=burnIn,saveAt="BB_")
      pred_ability_bb[i,1] <- cor(y[tst],fmBB$yHat[tst], use="complete")
    }

    # RKHS
    ETA<-list(MRK=list(K=G,model="RKHS"))

    pred_ability_rkhs <- matrix(nrow=nRep, ncol=1, NA)
    for(i in 1:nRep){
      tst <- sample(1:n, size=nTST, replace=FALSE)
      yNA <- y; yNA[tst]=NA
      fmRKHS <-BGLR(y=yNA,ETA=ETA, nIter=nIter, burnIn=burnIn,saveAt="RKHS_")
      pred_ability_rkhs[i,1] <- cor(y[tst],fmRKHS$yHat[tst], use="complete")
    }

    pred_ability <- data.frame(RRBLUP = pred_ability_rrblup,
                           BB = pred_ability_bb,
                           RHKS = pred_ability_rkhs)
    return(pred_ability)
  }

}
