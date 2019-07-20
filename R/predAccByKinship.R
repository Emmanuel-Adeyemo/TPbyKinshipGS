# May 31, 2019
#' CD mean algorithm
#'
#' @description
#' This is the heart of the package.
#'
#'
#' @param lambda Lambda is calculated as ve/va, it is needed to estimate the CDmean. You should provide a number
#' @param n_selected The number of lines to be selected as training population. Plese provide a number
#'
#'
#' @return
#' A dataframe of prediction accuracy for traits
#'
#' @import dplyr
#' @import rrBLUP
#' @import data.table
#'
#' @export
#'
get.pred.accuracy.by.kinship <- function(lambda, n_selected){


  for (m in 1:2){
    # CDmean starts
    ki_k = rm.col(ki)
    A_matrix = ki_k
    invA1 = make.Ainverse(A_matrix)
    n_total = assign.popsize(A_matrix)
    n_selected = 30
    #selected = solve.CDmean(1.3, 30, n_total, A_matrix, invA1,pheno) # returns training population for cluster
    selected = solve.CDmean(1.3, 30, assign.popsize(A_matrix), rm.col(ki), make.Ainverse(A_matrix),pheno) # returns training population for cluster

    # CDmean ends
    not_selected = get.not.selected(selected, pheno) # gets validation population within cluster
    train_kinship = do.innerjoin(selected, ki) #gets kinship for training pop
    validate_kinship = do.innerjoin(not_selected, ki) # gets kinship for validation pop
    validate_id = get.lineID(validate_kinship) # gets line_id for validation pop

    ##########################################################################################################
    # Prediction with Training population selected by kinship

    # make an empty dataframe
    predictions_all <- data.frame()

    for (val_line in validate_id){

      train_kinship_selected = get.train.by.kinship(val_line, train_kinship, 0.25) # selects training pop
      # by kinship used to train model # 0.25 is the threshold provided

      train_pheno = do.innerjoin(train_kinship_selected, pheno) # gets the phenotypic data for training set selected by kinship
      train_geno = do.innerjoin(train_kinship_selected, geno) # gets genotypic data for training set selected by kinship

      validate_pheno = get.validation.data(pheno, val_line) # gets phenotypic data for line being predicted
      validate_geno = get.validation.data(geno, val_line) # gets genotypic data for line being predicted

      u_stp_inc = make.predictions(train_pheno, train_pheno$StP_INC, train_geno, StP_INC_vp, validate_geno)
      v_stp_inc = validate_pheno$StP_INC #print(u_stp_inc,v_stp_inc)

      u_stp_sev = make.predictions(train_pheno,train_pheno$StP_SEV, train_geno, StP_SEV_vp, validate_geno)
      v_stp_sev = validate_pheno$StP_SEV#; put(u_stp_sev, v_stp_sev)

      u_stp_dis = make.predictions(train_pheno,train_pheno$StP_DIS, train_geno, StP_DIS_vp, validate_geno)
      v_stp_dis = validate_pheno$StP_DIS#; put(u_stp_dis, v_stp_dis)

      u_stp_microtwt = make.predictions(train_pheno,train_pheno$StP_Microtwt, train_geno, StP_Microtwt_vp, validate_geno)
      v_stp_microtwt = validate_pheno$StP_Microtwt#; put(u_stp_microtwt, v_stp_microtwt)

      u_stp_vsk = make.predictions(train_pheno,train_pheno$StP_VSK, train_geno, StP_VSK_vp, validate_geno)
      v_stp_vsk = validate_pheno$StP_VSK#; put(u_stp_vsk, v_stp_vsk)

      u_crk_inc = make.predictions(train_pheno,train_pheno$Crk_INC, train_geno, Crk_INC_vp, validate_geno)
      v_crk_inc = validate_pheno$Crk_INC#; put(u_crk_inc,v_crk_inc)

      u_crk_sev = make.predictions(train_pheno,train_pheno$Crk_SEV, train_geno, Crk_SEV_vp, validate_geno)
      v_crk_sev = validate_pheno$Crk_SEV#; put(u_crk_sev,v_crk_sev)

      u_crk_dis = make.predictions(train_pheno,train_pheno$Crk_DIS, train_geno, Crk_DIS_vp, validate_geno)
      v_crk_dis = validate_pheno$Crk_DIS#; put(u_crk_dis,v_crk_dis)

      u_crk_microtwt = make.predictions(train_pheno,train_pheno$Crk_Microtwt, train_geno, Crk_Microtwt_vp, validate_geno)
      v_crk_microtwt = validate_pheno$Crk_Microtwt#; put(u_crk_microtwt,v_crk_microtwt)

      u_crk_vsk = make.predictions(train_pheno,train_pheno$Crk_VSK, train_geno, Crk_VSK_vp, validate_geno)
      v_crk_vsk = validate_pheno$Crk_VSK#; put(u_crk_vsk,v_crk_vsk)

      combined <- data.frame(STP_INC = u_stp_inc, STP_SEV = u_stp_sev, STP_DIS = u_stp_dis, STP_Microtwt = u_stp_microtwt, STP_VSK = u_stp_vsk,
                             STPV_INC = v_stp_inc, STPV_SEV = v_stp_sev, STPV_DIS = v_stp_dis, STPV_Microtwt = v_stp_microtwt, STPV_VSK = v_stp_vsk,
                             Crk_INC = u_crk_inc, Crk_SEV = u_crk_sev, Crk_DIS = u_crk_dis, Crk_Microtwt = u_crk_microtwt, Crk_VSK = u_crk_vsk,
                             CrkV_INC = v_crk_inc, CrkV_SEV = v_crk_sev, CrkV_DIS = v_crk_dis, CrkV_Microtwt = v_crk_microtwt, CrkV_VSK = v_crk_vsk)


      predictions_all <- rbind(predictions_all, combined)

    }

    stp_inc_r <- cor(predictions_all$STP_INC, predictions_all$STPV_INC, use = "complete.obs") # finds correlation i.e. prediction accuracy
    stp_sev_r <- cor(predictions_all$STP_SEV, predictions_all$STPV_SEV, use = "complete.obs") # finds correlation i.e. prediction accuracy
    stp_dis_r <- cor(predictions_all$STP_DIS, predictions_all$STPV_DIS, use = "complete.obs") # finds correlation i.e. prediction accuracy
    stp_microtwt_r <- cor(predictions_all$STP_Microtwt, predictions_all$STPV_Microtwt, use = "complete.obs") # finds correlation i.e. prediction accuracy
    stp_vsk_r <- cor(predictions_all$STP_VSK, predictions_all$STPV_VSK, use = "complete.obs") # finds correlation i.e. prediction accuracy
    crk_inc_r <- cor(predictions_all$Crk_INC, predictions_all$CrkV_INC, use = "complete.obs") # finds correlation i.e. prediction accuracy
    crk_sev_r <- cor(predictions_all$Crk_SEV, predictions_all$CrkV_SEV, use = "complete.obs") # finds correlation i.e. prediction accuracy
    crk_dis_r <- cor(predictions_all$Crk_DIS, predictions_all$CrkV_DIS, use = "complete.obs") # finds correlation i.e. prediction accuracy
    crk_microtwt_r <- cor(predictions_all$Crk_Microtwt, predictions_all$CrkV_Microtwt, use = "complete.obs") # finds correlation i.e. prediction accuracy
    crk_vsk_r <- cor(predictions_all$Crk_VSK, predictions_all$CrkV_VSK, use = "complete.obs") # finds correlation i.e. prediction accuracy

    all_predictions = rbind(stp_inc_r, stp_sev_r, stp_dis_r, stp_microtwt_r, stp_vsk_r,
                            crk_inc_r, crk_sev_r, crk_dis_r, crk_microtwt_r, crk_vsk_r)

    vpred = file.path(paste0("data/results/clusterThree/2017_F5_VP_predictions_accuaracies_from_kinship_TP_iter", m, ".txt"))
    write.table(all_predictions, file = vpred, sep = "\t")

  }

}
