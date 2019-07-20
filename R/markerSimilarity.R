#' July 1, 2019
#'
#' Calculates marker similarity
#'
#' @description
#' Calculates marker similarity between two lines and do a pairwise comparison for
#' all lines present in the dataset
#' 
#'
#' @param geno A dataframe of lines with marker data. Cols = Line names, Rows = SNPs
#' 
#'
#' @export
#'
calculate.marker.similarity <- function(geno){
  
  geno <- geno[-1]
  
  geno2 <- geno
  colnames(geno2) <- NULL
  
  pairs <- combn(names(geno), 2, simplify=FALSE) # does name pairiing
  pairs_mar <- combn(geno2, 2,simplify = FALSE) # does marker pairing
  
  #Walking through the parents-Molecular Markers for each pop, obtaining CS
  coefPopAll <- NULL
  
  for(h in 1:length(pairs_mar)){ #for each pop
    
    p1 <- pairs_mar[[h]][, 1] #parents from each pop
    p2 <- pairs_mar[[h]][, 2]
    
    pTog <- rbind(p1, p2) #parents together
    
    coef.sim.all <- NULL
    for(m in 1:ncol(pTog)){
      if(!is.na(pTog[1,m]) & !is.na(pTog[2,m])){ #to avoid NA
        if(pTog[1,m] == pTog[2,m]){ #equal homozygous state 
          coef.sim <- 1 #coefficient of 1
        }
        if(pTog[1,m] != pTog[2,m]){ #different homozygous state (MM vs mm)
          coef.sim <- 0 #coefficient of 0
        }
        
      }else{
        coef.sim <- "NA"
      }
      coef.sim.all <- c(coef.sim.all, coef.sim) #saving coefs
    }
    coef.sim.all <- as.numeric(coef.sim.all) #CS between parents for the pop
    coef.all.pairs.av <- round(mean(coef.sim.all, na.rm = T), 2) #average
    coefPop <- c(coef.all.pairs.av, h)
    coefPopAll <- rbind(coefPopAll, coefPop)
    
  }   
  
  
  coef.df <- data.frame(line_1=rep(0,length(pairs)), line_2=rep(0,length(pairs)), Pop=c(1:length(pairs)))
  
  for(i in 1:length(pairs)){
    coef.df[i, 1] <- pairs[[i]][1]
    coef.df[i, 2] <- pairs[[i]][2]
    
  }
  
  
  coefPopAll <- as.data.frame(coefPopAll)
  colnames(coefPopAll) <- c("CS", "Pop")
  
  marker_sim <- inner_join(coef.df, coefPopAll, by = "Pop")
  
  return(marker_sim)
  
  
}
