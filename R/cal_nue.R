cal_nue <- function(Tg_pred,PPFD_pred,vpd_pred,fAPAR_pred,age_pred,CNrt_pred,LMA_pred,vcmax25_pred){
  npp_f <- summary(bp_model)$coefficients[1,1] +  
    summary(bp_model)$coefficients[2,1] * Tg_pred +
    summary(bp_model)$coefficients[3,1] * fAPAR_pred +
    summary(bp_model)$coefficients[4,1] * log(PPFD_pred) +
    summary(bp_model)$coefficients[5,1] * log(CNrt_pred)+
    summary(bp_model)$coefficients[6,1] * log(age_pred)
  
  npp_f[npp_f<=0] <- NA
  
  anpp_f <- npp_f * 
    (1/(1 + exp(-(summary(anpp_tnpp_model)$coefficients[1,1]+
                    summary(anpp_tnpp_model)$coefficients[2,1] * log(CNrt_pred) +
                    summary(anpp_tnpp_model)$coefficients[3,1] * log(PPFD_pred) + 
                    summary(anpp_tnpp_model)$coefficients[4,1] * Tg_pred+
                    summary(anpp_tnpp_model)$coefficients[5,1] * log(age_pred)))))
  
  bnpp_f <- npp_f - anpp_f
  
  lnpp_f <- anpp_f * (1/(1 + exp(-(summary(anpp_leafnpp_model)$coefficients[1,1]+
                                     summary(anpp_leafnpp_model)$coefficients[2,1]* fAPAR_pred +
                                     summary(anpp_leafnpp_model)$coefficients[3,1] * log(vpd_pred) + 
                                     summary(anpp_leafnpp_model)$coefficients[4,1] * log(PPFD_pred)))))
  
  wnpp_f <- anpp_f - lnpp_f
  
  leafnc_f <- (summary(n1)$coefficients[1,1]/leaf_c_forest) + 
    (summary(n1)$coefficients[2,1]/leaf_c_forest) * vcmax25_pred/LMA_pred
  
  nre_f <- (1/(1+exp(-(summary(nre_model)$coefficients[1,1] +
                         summary(nre_model)$coefficients[2,1] *Tg_pred +
                         summary(nre_model)$coefficients[3,1] * log(vpd_pred)))))
  
  lnf_f <- (1-nre_f)* leafnc_f * lnpp_f
  wnf_f <- wood_tissue_percentage*wnpp_f/wood_cn_forest

  bnf_f <- bnpp_f/root_cn_forest

  nuptake_f <- lnf_f + wnf_f + bnf_f

  new_nue <-available_grid2*npp_f/nuptake_f
  
  nue_pft_ratio <- log(nue_final/new_nue)
  
  return(nue_pft_ratio)
}