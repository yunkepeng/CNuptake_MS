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
  
  leafnc_f <- (summary(n1)$coefficients[1,1]/0.47) + 
    (summary(n1)$coefficients[2,1]/0.47) * vcmax25_pred/LMA_pred
  
  nre_f <- (1/(1+exp(-(summary(nre_model)$coefficients[1,1] +
                         summary(nre_model)$coefficients[2,1] *Tg_pred +
                         summary(nre_model)$coefficients[3,1] * log(vpd_pred)))))
  
  lnf_f <- (1-nre_f)* leafnc_f * lnpp_f
  wnf_f <- wnpp_f/100
  #100 is constant wood c/n
  bnf_f <- bnpp_f/94
  #94 is constant root c/n
  nuptake_f <- lnf_f + wnf_f + bnf_f
  
  #grass
  npp_g <- summary(bp_grass_model)$coefficients[1,1] +  
    summary(bp_grass_model)$coefficients[2,1] * log(PPFD_pred) +
    summary(bp_grass_model)$coefficients[3,1] * Tg_pred 
  
  anpp_g <- npp_g* 0.49
  bnpp_g <- npp_g-anpp_g
  
  leafnc_g <- 1/18
  
  nre_g <- 0.69
  
  lnf_g <- anpp_g*leafnc_g*(1-nre_g)
  
  bnf_g <- bnpp_g *(1/41)
  #41 is constant root c/n
  nuptake_g <- lnf_g + bnf_g
  
  #combine into 2 pfts
  npp_pft <- available_grid2* (npp_f*forest_percent +npp_g*grass_percent)
  npp_forest <- available_grid2* (npp_f*forest_percent)
  npp_grass <- available_grid2* (npp_g*grass_percent)
  
  anpp_pft <- available_grid2*(anpp_f*forest_percent +anpp_g*grass_percent)
  anpp_forest <- available_grid2* (anpp_f*forest_percent)
  anpp_grass <- available_grid2* (anpp_g*grass_percent)
  
  lnpp_forest <- available_grid2*lnpp_f*forest_percent
  
  wnpp_forest <- available_grid2*wnpp_f*forest_percent
  
  bnpp_pft <- available_grid2*(bnpp_f*forest_percent +bnpp_g*grass_percent)
  bnpp_forest <- available_grid2* (bnpp_f*forest_percent)
  bnpp_grass <- available_grid2* (bnpp_g*grass_percent)
  
  leafcn_pft <- 1/(available_grid2*(leafnc_f*forest_percent +leafnc_g*grass_percent))
  leafcn_forest <- 1/available_grid2*leafnc_f*forest_percent
  leafcn_grassland <- 1/available_grid2*leafnc_g*grass_percent
  summary(leafcn_pft)
  
  nre_pft <- available_grid2*(nre_f*forest_percent +nre_g*grass_percent)
  nre_forest <-  available_grid2*nre_f*forest_percent
  nre_grassland <- available_grid2*nre_g*grass_percent
  summary(nre_pft)
  
  lnf_pft <- available_grid2*(lnf_f*forest_percent +lnf_g*grass_percent)
  lnf_forest <- available_grid2* (lnf_f*forest_percent)
  lnf_grass <- available_grid2* (lnf_g*grass_percent)
  
  wnf_forest <- available_grid2*wnf_f*forest_percent
  
  bnf_pft <- available_grid2*(bnf_f*forest_percent +bnf_g*grass_percent)
  bnf_forest <- available_grid2* (bnf_f*forest_percent)
  bnf_grass <- available_grid2* (bnf_g*grass_percent)
  
  nuptake_pft <- available_grid2*(nuptake_f*forest_percent +nuptake_g*grass_percent)
  nuptake_forest <- available_grid2* (nuptake_f*forest_percent)
  nuptake_grass <- available_grid2* (nuptake_g*grass_percent)
  
  new_nue <-npp_pft/nuptake_pft
  
  nue_pft_ratio <- log(nue_final/new_nue)
  
  return(nue_pft_ratio)
}