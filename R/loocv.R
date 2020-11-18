
#' @description Leave-one-out cross validation for genomic prediction
#'
#' @param geno dataframe of genotypes
#' @param trait column with trait name
#' @param trait_id trait name in quotation marks
#' @param pheno dataframe of phenotypes


loocv <- function(geno, trait, trait_id, pheno){

  require(rrBLUP)

  nokinship_prediction <- vector(length=nrow(geno))

  for(i in 1:nrow(geno)){

    train_pheno <- pheno
    trait[i] <- NA
    train_pheno_remaining <- train_pheno[-i,]
    train_lines_remaining <- train_pheno_remaining[1]
    colnames(train_pheno_remaining)[1] <- "ID"
    train_geno <- geno[-i, ]
    train_geno_remaining <- train_geno
    valid_geno <- geno[i,]
    valid_line_id <- train_pheno[i,1]
    valid_line_id <- as.vector(valid_line_id)


    fit_remaining <- mixed.solve(y = train_pheno_remaining[,trait_id], Z = as.matrix(train_geno_remaining), method = "REML")


    nokinship_u <- setNames(fit_remaining$u, colnames(train_geno_remaining))
    valid_geno <- as.matrix(valid_geno)

    nokinship_g <- valid_geno %*% nokinship_u
    nokinship_prediction[i] <- nokinship_g


  }

  accuracy_nokinship = cor(nokinship_prediction, pheno[,trait_id], use = "complete.obs")
  return(accuracy_nokinship)
}


loocv_fixed <- function(geno, trait, trait_id, pheno){

  require(rrBLUP)

  nokinship_prediction <- vector(length=nrow(geno))

  for(i in 1:nrow(geno)){

    train_pheno <- pheno
    trait[i] <- NA
    train_pheno_remaining <- train_pheno[-i,]
    train_lines_remaining <- train_pheno_remaining[1]
    colnames(train_pheno_remaining)[1] <- "ID"
    train_geno <- geno[-i, ]
    train_geno_remaining <- train_geno
    valid_geno <- geno[i,]
    valid_line_id <- train_pheno[i,1]
    valid_line_id <- as.vector(valid_line_id)
    mf <- model.frame(train_pheno_remaining[,trait_id] ~ ID + Year, data = train_pheno_remaining)
    colnames(mf)[1] <- trait_id
    y <- model.response(mf)
    mf$Year <- as.factor(mf$Year)
    X <- model.matrix(ID ~ Year, data = mf)

    fit_remaining <- mixed.solve(y = y, Z = as.matrix(train_geno_remaining), X = X,  method = "REML")


    nokinship_u <- setNames(fit_remaining$u, colnames(train_geno_remaining))
    valid_geno <- as.matrix(valid_geno)

    nokinship_g <- valid_geno %*% nokinship_u
    nokinship_prediction[i] <- nokinship_g


  }

  accuracy_nokinship = cor(nokinship_prediction, pheno[,trait_id], use = "complete.obs")
  return(accuracy_nokinship)
}

predict.breedingPop <- function(train_geno, trait, trait_id, train_pheno, breed_geno){

  require(rrBLUP)

  breed_geno <- as.matrix(breed_geno)

  fit_remaining <- mixed.solve(y = train_pheno[,trait_id], Z = as.matrix(train_geno), method = "REML")

  marker_effects <- as.data.frame(fit_remaining$u)
  marker_effects <- as.matrix(marker_effects)

  breed_pred <- breed_geno %*% marker_effects

  breed_pred_final <- (breed_pred[,1]) + fit_remaining$beta

  breed_pred_final <- as.data.frame(breed_pred_final)
  breed_pred_final <- cbind(ID = rownames(breed_pred_final), breed_pred_final)
  colnames(breed_pred_final)[2] <- trait_id

  return(breed_pred_final)
}


loocv.gblup <- function(geno, trait, trait_id, pheno){

  require(rrBLUP)

  nokinship_prediction <- vector(length=nrow(geno))

  for(i in 1:nrow(geno)){

    train_pheno <- pheno
    trait[i] <- NA
    train_pheno_remaining <- train_pheno[-i,]
    train_lines_remaining <- train_pheno_remaining[1]
    colnames(train_pheno_remaining)[1] <- "ID"
    train_geno <- geno[-i, ]
    train_geno_remaining <- train_geno
    valid_geno <- geno[i,]
    valid_line_id <- train_pheno[i,1]
    valid_line_id <- as.vector(valid_line_id)

    train_geno_remaining <- as.matrix(train_geno_remaining)

    crossP <- tcrossprod(train_geno_remaining)/ncol(train_geno_remaining)
    idenMat <- diag(nrow(train_geno_remaining))

    fit_remaining <- mixed.solve(y = train_pheno_remaining[,trait_id], K = A.mat(train_geno_remaining), method = "REML")


    nokinship_u <- setNames(fit_remaining$u, rownames(train_geno_remaining))
    valid_geno <- as.matrix(valid_geno)

    nokinship_g <- valid_geno %*% nokinship_u
    nokinship_prediction[i] <- nokinship_g


  }

  accuracy_nokinship = cor(nokinship_prediction, pheno[,trait_id], use = "complete.obs")
  return(accuracy_nokinship)
}
