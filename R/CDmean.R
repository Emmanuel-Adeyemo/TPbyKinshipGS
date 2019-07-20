# May 21, 2019
#' CDmean algorithm
#'
#' @description
#' Calculates CDmean from contrasts matrix
#'
#' @references
#' #' Rincent, R., Laloe, D., Nicolas, S., Altmann, T., Brunel, D., Revilla, P.,
#' Moreau, L. (2012). Maximizing the Reliability of Genomic Selection by
#' Optimizing the Calibration Set of Reference Individuals: Comparison of
#' Methods in Two Diverse Groups of Maize Inbreds (Zea mays L.). Genetics,
#' 192(2), 715-728. http://doi.org/10.1534/genetics.112.141473
#'
#' @param n_selected Choose a size for your calibration set. You should provide a number
#' @param lambda Lambda is calculated as ve/va, it is needed to estimate the CDmean.
#' @param n_total Total number of lines in the population. Use the variable name "n_total"
#' @param A_matrix A matrix of additive relationship among individuals in the whole population. Use the variable name "A_matrix"
#' @param invA1 A inverse matrix of the additive relationship. Use the variable name "invA1"
#' @param pheno A dataframe of all lines
#' @param max_iter The number of iterations for the exchange algorithm. Provide a number.
#'
#'
#' @return
#' A dataframe of selected lines
#'
#' @import dplyr
#' @import rrBLUP
#' @import data.table
#'
#' @export
#'
solve.CDmean <- function(lambda, n_selected, n_total, A_matrix, invA1,pheno, max_iter){


  # Design matrices
  # M matrix or the orthogonal projector
  Ident<-diag(n_selected) # Identity matrix of TP size

  X<-rep(1,n_selected)

  M <- Ident- (X%*%solve(t(X)%*%X) %*% t(X)) ## Inverse of X is simply 1/n_total, but this is the generalize inverse

  tp_phenotype<-sample(n_total,n_selected) #Calibration set initialization

  save_tp=tp_phenotype # tp_phenotype = Sample1, save_tp = SaveSample1, NotSampled1 = all_phenotype, NotSampled = not_phenotyped

  all_phenotype<-seq(1:n_total)

  not_phenotyped<-all_phenotype[-tp_phenotype] # Initial validation set

  Z=matrix(0,n_selected,n_total)

  for (i in 1:length(tp_phenotype)) {

    Z[i,tp_phenotype[i]]=1

    }

  T<-make.contrast(not_phenotyped)   # T matrix of contrasts

  # Calculate of CDmean of the initial set
  matCD<-(t(T)%*%(A_matrix-lambda*ginv(t(Z)%*%M%*%Z + lambda*invA1))%*%T)/(t(T)%*%A_matrix%*%T)

  CD=diag(matCD)
  CDmeanSave=mean(CD)

  CDmeanMax1=rep(NA,max_iter)

  # Exchange algorithm (maximize CDmean)
  cpt2=1
  cpt=0
  while (cpt2<max_iter) {  # Make sure that 1000 is enough in your case (that you reached a plateau), for this look at CDmeanMax1.

    not_phenotyped=all_phenotype[-tp_phenotype]
    cpt2=cpt2+1

    # Remove one individual (randomly choosen) from the sample selected for phenotyping (calibration set) :
    Sample2=sample(tp_phenotype,1)

    # Select one individual (randomly choosen) from the individuals that are not in the Calibration set :
    Sample3=sample(not_phenotyped,1)

    # New calibration set : sample 3 and evrything in the calibration set minus sample 2
    Sample4=c(Sample3,tp_phenotype[tp_phenotype!=Sample2])

    # Calculate the mean CD of the new calibration set :
    Z=matrix(0,n_selected,n_total)

    for (i in 1:length(Sample4)) {

      Z[i,Sample4[i]]=1

      }
    not_phenotyped=all_phenotype[-Sample4]
    T<-make.contrast(not_phenotyped)

    matCD<-(t(T)%*%(A_matrix-lambda*ginv(t(Z)%*%M%*%Z + lambda*invA1))%*%T)/(t(T)%*%A_matrix%*%T)
    CD=diag(matCD)

    if (mean(CD)>CDmeanSave ) { tp_phenotype=Sample4 # Accept the new Calibration set if CDmean is increased, reject otherwise.

    CDmeanSave=mean(CD)
    cpt=0

    } else {

      cpt=cpt+1

    }
    CDmeanMax1[cpt2-1]=CDmeanSave
  }  #End of Loop



  SampleOptimized <- tp_phenotype # SampleOptimized is the optimized calibration set

  # End
  linename <- as.data.frame(pheno[,1])
  colnames(linename)[1] <- "ID"

  geno.lines2 <- linename$ID
  optim.lines <- geno.lines2[SampleOptimized]
  optim.lines <- as.data.frame(optim.lines)
  colnames(optim.lines) <- "ID"

  k_selected = optim.lines
  return(k_selected)



} # Close the function





