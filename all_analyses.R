rm(list=ls())
library(MLmetrics)
library(tidyverse) 
library(ncmeta)
library(viridis)
library(ggthemes)
library(LSD)
library(yardstick)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(gplots)
library(tidyselect)
library(extrafont)
library(caret)
library(recipes)
devtools::load_all("/Users/yunpeng/yunkepeng/rbeni/")
#library(rbeni)
library(raster)
library(maps)
library(rworldmap)
library(cowplot)
library(ncdf4)
library(scales)
library(lme4)
library(lmerTest)
library("PerformanceAnalytics")
library(MuMIn)
library(tidyverse)
library(ggplot2)
library(lme4)
library(visreg)
library(ggpubr)
library(car)
library("ggplotify")
library(remotes)
library(tune)
library(relaimpo)

#reset validation metrics info

analyse_modobs2 <- function(
    df,
    mod,
    obs,
    type       = "points",
    filnam     = NA,
    relative   = FALSE,
    xlim       = NULL,
    ylim       = NULL,
    use_factor = NULL,
    shortsubtitle = FALSE,
    plot_subtitle = TRUE,
    plot_linmod = TRUE,
    ...
){
  
  require(ggplot2)
  require(dplyr)
  require(LSD)
  require(ggthemes)
  require(RColorBrewer)
  
  #if (identical(filnam, NA)) filnam <- "analyse_modobs.pdf"
  
  ## rename to 'mod' and 'obs' and remove rows with NA in mod or obs
  df <- df %>%
    as_tibble() %>%
    ungroup() %>%
    dplyr::select(mod=mod, obs=obs) %>%
    tidyr::drop_na(mod, obs)
  
  ## get linear regression (coefficients)
  linmod <- lm( obs ~ mod, data=df )
  
  ## construct metrics table using the 'yardstick' library
  df_metrics <- df %>%
    yardstick::metrics(obs, mod) %>%
    dplyr::bind_rows( tibble( .metric = "n",        .estimator = "standard", .estimate = summarise(df, numb=n()) %>% unlist() ) ) %>%
    dplyr::bind_rows( tibble( .metric = "slope",    .estimator = "standard", .estimate = coef(linmod)[2]) ) %>%
    # dplyr::bind_rows( tibble( .metric = "nse",      .estimator = "standard", .estimate = hydroGOF::NSE( obs, mod, na.rm=TRUE ) ) ) %>%
    dplyr::bind_rows( tibble( .metric = "mean_obs", .estimator = "standard", .estimate = summarise(df, mean=mean(obs, na.rm=TRUE)) %>% unlist() ) ) %>%
    dplyr::bind_rows( tibble( .metric = "prmse",    .estimator = "standard",
                              .estimate = dplyr::filter(., .metric=="rmse") %>% dplyr::select(.estimate) %>% unlist() /
                                dplyr::filter(., .metric=="mean_obs") %>% dplyr::select(.estimate) %>% unlist() ) ) %>%
    dplyr::bind_rows( tibble( .metric = "pmae",    .estimator = "standard",
                              .estimate = dplyr::filter(., .metric=="mae") %>% dplyr::select(.estimate) %>% unlist() /
                                dplyr::filter(., .metric=="mean_obs") %>% dplyr::select(.estimate) %>% unlist() ) ) %>%
    dplyr::bind_rows( tibble( .metric = "bias",        .estimator = "standard", .estimate = dplyr::summarise(df, mean((mod-obs), na.rm=TRUE    )) %>% unlist() ) ) %>%
    dplyr::bind_rows( tibble( .metric = "pbias",       .estimator = "standard", .estimate = dplyr::summarise(df, mean((mod-obs)/obs, na.rm=TRUE)) %>% unlist() ) )
  
  rsq_val <- df_metrics %>% dplyr::filter(.metric=="rsq") %>% dplyr::select(.estimate) %>% unlist() %>% unname()
  rmse_val <- df_metrics %>% dplyr::filter(.metric=="rmse") %>% dplyr::select(.estimate) %>% unlist() %>% unname()
  mae_val <- df_metrics %>% dplyr::filter(.metric=="mae") %>% dplyr::select(.estimate) %>% unlist() %>% unname()
  bias_val <- df_metrics %>% dplyr::filter(.metric=="bias") %>% dplyr::select(.estimate) %>% unlist() %>% unname()
  slope_val <- df_metrics %>% dplyr::filter(.metric=="slope") %>% dplyr::select(.estimate) %>% unlist() %>% unname()
  n_val <- df_metrics %>% dplyr::filter(.metric=="n") %>% dplyr::select(.estimate) %>% unlist() %>% unname()
  
  if (relative){
    rmse_val <- rmse_val / mean(df$obs, na.rm = TRUE)
    bias_val <- bias_val / mean(df$obs, na.rm = TRUE)
  }
  
  rsq_lab <- format( rsq_val, digits = 2 )
  rmse_lab <- format( rmse_val, digits = 3 )
  mae_lab <- format( mae_val, digits = 3 )
  bias_lab <- format( bias_val, digits = 3 )
  slope_lab <- format( slope_val, digits = 3 )
  n_lab <- format( n_val, digits = 3 )
  
  results <- tibble( rsq = rsq_val, rmse = rmse_val, mae = mae_val, bias = bias_val, slope = slope_val, n = n_val )
  
  if (shortsubtitle){
    subtitle <- bquote( italic(R)^2 == .(rsq_lab) ~~
                          RRMSE == .(rmse_lab) )
  } else {
    subtitle <- bquote( italic(R)^2 == .(rsq_lab) ~~
                          RRMSE == .(rmse_lab) ~~
                          bias == .(bias_lab) ~~
                          slope == .(slope_lab) ~~
                          italic(N) == .(n_lab) )
  }
  
  if (type=="heat"){
    
    # if (!identical(filnam, NA)) dev.off()
    # source("~/LSD/R/LSD.heatscatter.R")
    
    gg <- heatscatter(
      df$mod,
      df$obs,
      xlim=xlim,
      ylim=ylim,
      main="",
      ggplot=TRUE )
    
    gg <- gg +
      geom_abline(intercept=0, slope=1, linetype="dotted") +
      theme_classic() +
      labs(x = mod, y = obs)
    
    if (plot_linmod) gg <- gg + geom_smooth(method='lm', color="red", size=0.5, se=FALSE)
    if (plot_subtitle) gg <- gg + labs(subtitle = subtitle)
    
    if (!identical(filnam, NA)) {
      ggsave(filnam, width=5, height=5)
    }
    
  } else if (type=="hex"){
    
    ## ggplot hexbin
    gg <- df %>%
      ggplot2::ggplot(aes(x=mod, y=obs)) +
      geom_hex() +
      scale_fill_gradientn(
        colours = colorRampPalette( c("gray65", "navy", "red", "yellow"))(5)) +
      geom_abline(intercept=0, slope=1, linetype="dotted") +
      # coord_fixed() +
      # xlim(0,NA) +
      # ylim(0,NA) +
      theme_classic() +
      labs(x = mod, y = obs)
    
    if (plot_subtitle) gg <- gg + labs(subtitle = subtitle)
    if (plot_linmod) gg <- gg + geom_smooth(method='lm', color="red", size=0.5, se=FALSE)
    
    if (!identical(filnam, NA)) {
      ggsave(filnam, width=5, height=5)
    }
    
  } else if (type=="points") {
    
    ## points
    gg <- df %>%
      ggplot(aes(x=mod, y=obs)) +
      geom_point() +
      geom_abline(intercept=0, slope=1, linetype="dotted") +
      # coord_fixed() +
      # xlim(0,NA) +
      # ylim(0,NA) +
      theme_classic() +
      labs(x = mod, y = obs)
    
    if (plot_subtitle) gg <- gg + labs(subtitle = subtitle)
    if (plot_linmod) gg <- gg + geom_smooth(method='lm', color="red", size=0.5, se=FALSE)
    
    if (!identical(filnam, NA)) {
      ggsave(filnam, width=5, height=5)
    }
    
  } else if (type=="density") {
    
    ## points
    gg <- df %>%
      ggplot(aes(x=mod, y=obs)) +
      
      stat_density_2d(aes(fill = after_stat(nlevel)), geom = "polygon") +
      scale_fill_gradientn(colours = colorRampPalette( c("gray65", "navy", "red", "yellow"))(5),
                           guide = "legend") +
      
      geom_abline(intercept=0, slope=1, linetype="dotted") +
      # coord_fixed() +
      # xlim(0,NA) +
      # ylim(0,NA) +
      theme_classic() +
      labs(x = mod, y = obs)
    
    if (plot_subtitle) gg <- gg + labs(subtitle = subtitle)
    if (plot_linmod) gg <- gg + geom_smooth(method='lm', color="red", size=0.5, se=FALSE)
    
    if (!identical(filnam, NA)) {
      ggsave(filnam, width=5, height=5)
    }
    
  }
  
  return(list(df_metrics=df_metrics, gg=gg, linmod=linmod, results = results))
}

white <- theme(plot.background=element_rect(fill="white", color="white"))

larger_size <- theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
                     plot.subtitle=element_text(size=15))

#1. input complete dataset and do some re-processing
stepwise <- function(df_input,target_var){
  #-----------------------------------------------------------------------
  # Input:  whole dataframe and target variable
  #assume that site_a is the only random factor
  #-----------------------------------------------------------------------
  target <- target_var
  df <- df_input
  
  preds <- df %>% dplyr::select(-c(target,site_a)) %>% 
    names()
  
  r_list <- c()
  #For loop functions, include all predictor's r2 at the end
  for (var in preds){
    forml <- paste( 'lmer(', target, '~', var, '+(1|site_a), data = df)')
    fit_lin <- eval(parse(text = forml)) 
    rsq <- r.squaredGLMM(fit_lin)[1]
    r_list <- c(r_list,rsq)
  }
  
  #convert to a dataframe, including all r2
  All_rsquare <- data.frame (
    preds = factor(preds,levels=preds), 
    rsq = r_list)
  
  #select max r2 in all predictors
  max(r_list)
  
  new_All_rsquare <- All_rsquare %>% 
    # desc orders from largest to smallest
    arrange(desc(rsq))
  
  #2. stepwise regression selection
  
  ## list
  list_aic <- list()
  list_bic <- list()
  list_R <- list()
  list_variable <- list()
  
  # predictors retained in the model firstly
  preds_retained <- as.character(new_All_rsquare[1,1])
  preds_candidate <- preds[-which(preds == preds_retained)] 
  
  
  for (a in 1:(length(preds)-1)){
    rsq_candidates <- c()
    linmod_candidates <- list()
    for (i in 1:length(preds_candidate)){
      pred_add <- c(preds_retained, preds_candidate[i])
      forml  <- paste( 'lmer(', target, '~', paste(pred_add, collapse = '+'), '+(1|site_a), data = df)')
      # create a function and make its format available to output in for loop
      fit_lin <- eval(parse(text = forml))
      linmod_candidates[[ i ]] <- fit_lin
      # obtain multiple r2 at each selection, and find the best one at the end
      rsq <- r.squaredGLMM(fit_lin)[1]
      rsq_candidates[i] <- rsq
    }
    pred_max <- preds_candidate[ which.max(rsq_candidates) ]
    # include best factors in retained factor
    preds_retained <- c(preds_retained, pred_max)
    list_variable[[a]] <- pred_max 
    # include AIC, BIC, adjusted R2, R2, cross-validated R2 and RMSE at each k 
    list_aic[[  a ]] <- AIC(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained, collapse = '+'),  '+(1|site_a), data = df)'))))
    
    list_bic[[ a ]] <- BIC(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained, collapse = '+'),  '+(1|site_a), data = df)'))))
    
    list_R[[ a ]] <- r.squaredGLMM(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained, collapse = '+'),  '+(1|site_a), data = df)'))))[1]
    preds_candidate <- preds_candidate[-which(preds_candidate == pred_max)]
  }
  
  
  R_null <- r.squaredGLMM(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained[1], collapse = '+'),  '+(1|site_a), data = df)'))))[1]
  AIC_null <- AIC(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained[1], collapse = '+'),  '+(1|site_a), data = df)'))))
  BIC_null <- BIC(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained[1], collapse = '+'),  '+(1|site_a), data = df)'))))
  variable_null <- preds_retained[1]
  
  R_all <- round(as.numeric(c(R_null,list_R)),2)
  AIC_all <- round(as.numeric(c(AIC_null,list_aic)),2)
  BIC_all <- round(as.numeric(c(BIC_null,list_bic)),2)
  variable_all <- (as.character(c(variable_null,list_variable)))
  
  df1 <- as.data.frame(cbind(variable_all,R_all,AIC_all,BIC_all))
  
  #Adjusted-R
  p1 <- ggplot() + 
    geom_point(data = df1, aes(x = factor(variable_all,level = variable_all), y = R_all)) 
  #AIC
  p2 <- ggplot() + 
    geom_point(data = df1, aes(x = factor(variable_all,level = variable_all), y = AIC_all)) 
  #BIC
  p3 <- ggplot() + 
    geom_point(data = df1, aes(x = factor(variable_all,level = variable_all), y = BIC_all))
  
  output_list <- list(p1,p2,p3)
  
  return(output_list)
  #-----------------------------------------------------------------------
  # Output: four figures 
  #-----------------------------------------------------------------------
}

#stepwise lm
stepwise_lm <- function(df_input,target_var){
  #-----------------------------------------------------------------------
  # Input:  whole dataframe and target variable
  #-----------------------------------------------------------------------
  target <- target_var
  df <- df_input
  
  preds <- df %>% dplyr::select(-c(target)) %>% 
    names()
  
  r_list <- c()
  #For loop functions, include all predictor's r2 at the end
  for (var in preds){
    forml <- paste( 'lm(', target, '~', var, ', data = df)')
    fit_lin <- eval(parse(text = forml)) 
    rsq <- r.squaredGLMM(fit_lin)[1]
    r_list <- c(r_list,rsq)
  }
  
  #convert to a dataframe, including all r2
  All_rsquare <- data.frame (
    preds = factor(preds,levels=preds), 
    rsq = r_list)
  
  #select max r2 in all predictors
  max(r_list)
  
  new_All_rsquare <- All_rsquare %>% 
    # desc orders from largest to smallest
    arrange(desc(rsq))
  
  #2. stepwise regression selection
  
  ## list
  list_aic <- list()
  list_bic <- list()
  list_R <- list()
  list_variable <- list()
  
  # predictors retained in the model firstly
  preds_retained <- as.character(new_All_rsquare[1,1])
  preds_candidate <- preds[-which(preds == preds_retained)] 
  
  
  for (a in 1:(length(preds)-1)){
    rsq_candidates <- c()
    linmod_candidates <- list()
    for (i in 1:length(preds_candidate)){
      pred_add <- c(preds_retained, preds_candidate[i])
      forml  <- paste( 'lm(', target, '~', paste(pred_add, collapse = '+'), ', data = df)')
      # create a function and make its format available to output in for loop
      fit_lin <- eval(parse(text = forml))
      linmod_candidates[[ i ]] <- fit_lin
      # obtain multiple r2 at each selection, and find the best one at the end
      rsq <- r.squaredGLMM(fit_lin)[1]
      rsq_candidates[i] <- rsq
    }
    pred_max <- preds_candidate[ which.max(rsq_candidates) ]
    # include best factors in retained factor
    preds_retained <- c(preds_retained, pred_max)
    list_variable[[a]] <- pred_max 
    # include AIC, BIC, adjusted R2, R2, cross-validated R2 and RMSE at each k 
    list_aic[[  a ]] <- AIC(eval(parse(text = paste( 'lm(', target, '~', paste(preds_retained, collapse = '+'),  ', data = df)'))))
    
    list_bic[[ a ]] <- BIC(eval(parse(text = paste( 'lm(', target, '~', paste(preds_retained, collapse = '+'),  ', data = df)'))))
    
    list_R[[ a ]] <- r.squaredGLMM(eval(parse(text = paste( 'lm(', target, '~', paste(preds_retained, collapse = '+'),  ', data = df)'))))[1]
    preds_candidate <- preds_candidate[-which(preds_candidate == pred_max)]
  }
  
  
  R_null <- r.squaredGLMM(eval(parse(text = paste( 'lm(', target, '~', paste(preds_retained[1], collapse = '+'),  ', data = df)'))))[1]
  AIC_null <- AIC(eval(parse(text = paste( 'lm(', target, '~', paste(preds_retained[1], collapse = '+'),  ', data = df)'))))
  BIC_null <- BIC(eval(parse(text = paste( 'lm(', target, '~', paste(preds_retained[1], collapse = '+'),  ', data = df)'))))
  variable_null <- preds_retained[1]
  
  R_all <- round(as.numeric(c(R_null,list_R)),2)
  AIC_all <- round(as.numeric(c(AIC_null,list_aic)),2)
  BIC_all <- round(as.numeric(c(BIC_null,list_bic)),2)
  variable_all <- (as.character(c(variable_null,list_variable)))
  
  df1 <- as.data.frame(cbind(variable_all,R_all,AIC_all,BIC_all))
  
  #Adjusted-R
  p1 <- ggplot() + 
    geom_point(data = df1, aes(x = factor(variable_all,level = variable_all), y = R_all)) 
  #AIC
  p2 <- ggplot() + 
    geom_point(data = df1, aes(x = factor(variable_all,level = variable_all), y = AIC_all)) 
  #BIC
  p3 <- ggplot() + 
    geom_point(data = df1, aes(x = factor(variable_all,level = variable_all), y = BIC_all))
  
  output_list <- list(p1,p2,p3)
  
  return(output_list)
  #-----------------------------------------------------------------------
  # Output: four figures 
  #-----------------------------------------------------------------------
}

#remove alpha throughout all analyses because 
#(1) D and alpha shows contridictory results for NUE 
#(2) D and alpha repeated in model selection

NPP_all <- read.csv("~/data/NPP_Yunke/NPP_Nmin_dataset_with_predictors.csv")

#summarise number of sites
dim(subset(NPP_all,is.na(Nmin)==TRUE) %>% group_by(site)  %>% summarise(mean = mean(lon)))

NPP_all$NPP.foliage[NPP_all$NPP.foliage==0] <-NA
NPP_all$NPP.wood[NPP_all$NPP.wood==0] <-NA
NPP_all$ANPP_2[NPP_all$ANPP_2==0] <-NA
NPP_all$age[NPP_all$age==0] <-NA

NPP_all$site_a <- NPP_all$site
NPP_all$tnpp_a <- NPP_all$TNPP_1
NPP_all$anpp_a <- NPP_all$ANPP_2

NPP_all$anpp_tnpp_a <-  log((NPP_all$ANPP_2/NPP_all$tnpp_a)/(1-(NPP_all$ANPP_2/NPP_all$tnpp_a)))
NPP_all$anpp_leafnpp_a <-  log((NPP_all$NPP.foliage/NPP_all$ANPP_2)/(1-(NPP_all$NPP.foliage/NPP_all$ANPP_2)))

NPP_all$soilCN_a <- log(NPP_all$soilCN)
NPP_all$observedfAPAR_a <- NPP_all$observedfAPAR
NPP_all$obs_age_a <- log(NPP_all$age)

NPP_all$age_a <- log(NPP_all$mapped_age)
NPP_all$Tg_a <- NPP_all$Tg
NPP_all$PPFD_a <- log(NPP_all$PPFD)
NPP_all$vpd_a <- log(NPP_all$vpd)
NPP_all$fAPAR_a <- NPP_all$fAPAR
NPP_all$CNrt_a <- log(NPP_all$CNrt)
NPP_all$LMA_a <- log(NPP_all$LMA)
NPP_all$vcmax25_a <- log(NPP_all$vcmax25)
NPP_all$ndep_a <- log(NPP_all$ndep)

NPP_forest <- subset(NPP_all,pft=="Forest")

#check why some grassland BP and ANPP is so high
outliers <- subset(NPP_all,pft=="Grassland" & ANPP_2>500)
newmap <- getMap(resolution = "low")
sp::plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)
points(outliers$lon,outliers$lat, col="green", pch=16,cex=2)

BP_dataset <- na.omit(NPP_forest[,c("tnpp_a","obs_age_a","observedfAPAR_a","soilCN_a","Tg_a","PPFD_a","vpd_a","site_a")])
#model1 <- stepwise(BP_dataset,"tnpp_a")
#model1[[1]]
#model1[[2]]
#bp_model <- (lmer(tnpp_a~Tg_a+observedfAPAR_a+obs_age_a+PPFD_a+alpha_a+(1|site_a),data=BP_dataset))
#summary(bp_model)

#ndep check - BP adding ndep - significant
BP_dataset_ndep <- na.omit(NPP_forest[,c("tnpp_a","ndep_a","age_a","fAPAR_a","CNrt_a","Tg_a","PPFD_a","vpd_a","site_a")])
dim(BP_dataset_ndep)
stepwise(BP_dataset_ndep,"tnpp_a")[[1]]
stepwise(BP_dataset_ndep,"tnpp_a")[[3]]
bp_model_ndep <- (lmer(tnpp_a~Tg_a+fAPAR_a+ndep_a+CNrt_a+age_a+(1|site_a),data=BP_dataset_ndep))
summary(bp_model_ndep)
r.squaredGLMM(bp_model_ndep)
#ndep check - ANPP/BP adding ndep - works!
anpp_tnpp_dataset_ndep <- na.omit(NPP_forest[,c("anpp_tnpp_a","ndep_a","age_a","fAPAR_a","CNrt_a","Tg_a","PPFD_a","vpd_a","site_a")])
stepwise(anpp_tnpp_dataset_ndep,"anpp_tnpp_a")[[1]]
anpp_tnpp_model_ndep <- (lmer(anpp_tnpp_a~ndep_a+CNrt_a+PPFD_a+Tg_a+(1|site_a),data=anpp_tnpp_dataset_ndep))
summary(anpp_tnpp_model_ndep)
r.squaredGLMM(anpp_tnpp_model_ndep)
#ndep check - leaf.npp/ANPP adding ndep - not improved! Ndep is non-significant
anpp_leafnpp_dataset_noage_ndep <- na.omit(NPP_forest[,c("anpp_leafnpp_a","ndep_a","Tg_a","PPFD_a","vpd_a","site_a")])
stepwise(anpp_leafnpp_dataset_noage_ndep,"anpp_leafnpp_a")[[1]]
r.squaredGLMM(lmer(anpp_leafnpp_a~Tg_a+PPFD_a+vpd_a+(1|site_a),data=anpp_leafnpp_dataset_noage_ndep))
summary(lmer(anpp_leafnpp_a~ndep_a+vpd_a+Tg_a+(1|site_a),data=anpp_leafnpp_dataset_noage_ndep))
r.squaredGLMM(lmer(anpp_leafnpp_a~ndep_a+vpd_a+Tg_a+(1|site_a),data=anpp_leafnpp_dataset_noage_ndep))

#now, start works
BP_dataset2 <- na.omit(NPP_forest[,c("tnpp_a","age_a","fAPAR_a","CNrt_a","Tg_a","PPFD_a","vpd_a","site_a")])
dim(BP_dataset2)
a2 <- stepwise(BP_dataset2,"tnpp_a")
a2[[1]]
a2[[2]]
a2[[3]]
bp_model <- (lmer(tnpp_a~Tg_a+fAPAR_a+PPFD_a+CNrt_a+age_a+(1|site_a),data=BP_dataset2))
summary(bp_model)
r.squaredGLMM(bp_model)

#check how many data were removed
nrow(BP_dataset)/nrow(BP_dataset2)

vif_bp <- vif((lmer(tnpp_a~Tg_a+fAPAR_a+PPFD_a+CNrt_a+age_a+vpd_a+(1|site_a),data=BP_dataset2)))

#anpp_tnpp_dataset <- na.omit(NPP_forest[,c("anpp_tnpp_a","obs_age_a","observedfAPAR_a","soilCN_a","Tg_a","PPFD_a","vpd_a","site_a")])
#dim(anpp_tnpp_dataset)
#model2 <- stepwise(anpp_tnpp_dataset,"anpp_tnpp_a")
#model2[[1]]
#model2[[2]]
#anpp_tnpp_model <- (lmer(anpp_tnpp_a~soilCN_a+obs_age_a+observedfAPAR_a+(1|site_a),data=anpp_tnpp_dataset))
#summary(anpp_tnpp_model)

#mapped
anpp_tnpp_dataset2 <- na.omit(NPP_forest[,c("anpp_tnpp_a","age_a","fAPAR_a","CNrt_a","Tg_a","PPFD_a","vpd_a","site_a")])
dim(subset(NPP_forest,ANPP_2>0))
model2a <- stepwise(anpp_tnpp_dataset2,"anpp_tnpp_a")
model2a[[1]]
model2a[[2]]
model2a[[3]]
anpp_tnpp_model <- (lmer(anpp_tnpp_a~CNrt_a+PPFD_a+Tg_a+age_a+(1|site_a),data=anpp_tnpp_dataset2))
summary(anpp_tnpp_model)
r.squaredGLMM(anpp_tnpp_model)

vif_anpp_tnpp <- vif((lmer(anpp_tnpp_a~Tg_a+fAPAR_a+PPFD_a+CNrt_a+age_a+vpd_a+(1|site_a),data=anpp_tnpp_dataset2)))

#with age, but age should be removed since it shows higher AIC and lower R2
anpp_leafnpp_dataset_age <- na.omit(NPP_forest[,c("anpp_leafnpp_a","age_a","fAPAR_a","CNrt_a","Tg_a","PPFD_a","vpd_a","site_a")])
model3a <- stepwise(anpp_leafnpp_dataset_age,"anpp_leafnpp_a")
model3a[[1]]
model3a[[2]]
model3a[[3]]
test <- (lmer(anpp_leafnpp_a~age_a+PPFD_a+vpd_a+(1|site_a),data=anpp_leafnpp_dataset_age)) 
summary(test)

#AIC: 1359
r.squaredGLMM(test)
AIC(test)
BIC(test)

#AIC: 1351
r.squaredGLMM(lmer(anpp_leafnpp_a~Tg_a+PPFD_a+vpd_a+(1|site_a),data=anpp_leafnpp_dataset_age)) 
AIC(lmer(anpp_leafnpp_a~Tg_a+PPFD_a+vpd_a+(1|site_a),data=anpp_leafnpp_dataset_age)) 
BIC(lmer(anpp_leafnpp_a~Tg_a+PPFD_a+vpd_a+(1|site_a),data=anpp_leafnpp_dataset_age)) 

#AIC: 1343
r.squaredGLMM(lmer(anpp_leafnpp_a~fAPAR_a+PPFD_a+vpd_a+(1|site_a),data=anpp_leafnpp_dataset_age)) 
AIC(lmer(anpp_leafnpp_a~fAPAR_a+PPFD_a+vpd_a+(1|site_a),data=anpp_leafnpp_dataset_age)) 
BIC(lmer(anpp_leafnpp_a~fAPAR_a+PPFD_a+vpd_a+(1|site_a),data=anpp_leafnpp_dataset_age)) 

#without age - re-selection - this is best
anpp_leafnpp_dataset <- na.omit(NPP_forest[,c("anpp_leafnpp_a","fAPAR_a","CNrt_a","Tg_a","PPFD_a","vpd_a","site_a")])
model3 <- stepwise(anpp_leafnpp_dataset,"anpp_leafnpp_a")
model3[[1]]
model3[[3]]
anpp_leafnpp_model <- (lmer(anpp_leafnpp_a~fAPAR_a+vpd_a+PPFD_a+(1|site_a),data=anpp_leafnpp_dataset)) 
r.squaredGLMM(anpp_leafnpp_model)
AIC(anpp_leafnpp_model)
summary(anpp_leafnpp_model)
r.squaredGLMM(anpp_leafnpp_model)
vif_anpp_leafnpp <- vif((lmer(anpp_leafnpp_a~Tg_a+fAPAR_a+PPFD_a+CNrt_a+age_a+vpd_a+(1|site_a),data=anpp_leafnpp_dataset_age)))

#check tnpp grassland
#not filtering any management/non-management! while previous do so
#removing tiandi's grassland
NPP_grassland <- subset(NPP_all,pft=="Grassland" & is.na(Nmin)==TRUE)
grassland_sitemean <- aggregate(NPP_grassland,by=list(NPP_grassland$site), FUN=mean, na.rm=TRUE) 

BP_dataset_grass <- na.omit(grassland_sitemean[,c("tnpp_a","Tg_a","PPFD_a","vpd_a","CNrt_a","fAPAR_a")])
model_g1 <- stepwise_lm(BP_dataset_grass,"tnpp_a")
model_g1[[1]]
model_g1[[2]]
model_g1[[3]]

bp_grass_model <- (lm(tnpp_a~PPFD_a+Tg_a,data=BP_dataset_grass))
summary(bp_grass_model)
r.squaredGLMM(bp_grass_model)

vif_bp_grass <- vif((lm(tnpp_a~Tg_a+PPFD_a+vpd_a+CNrt_a+fAPAR_a,data=BP_dataset_grass)))

#anpp/tnpp
dim(subset(grassland_sitemean,TNPP_1>0))
dim(subset(grassland_sitemean,TNPP_1>0 & ANPP_2>0))

anpp_tnpp_dataset_grass <- na.omit(grassland_sitemean[,c("anpp_tnpp_a","Tg_a","PPFD_a","vpd_a","CNrt_a","fAPAR_a")])
model_g2 <- stepwise_lm(anpp_tnpp_dataset_grass,"anpp_tnpp_a")
model_g2[[1]]
summary(lm(anpp_tnpp_a~Tg_a,data=anpp_tnpp_dataset_grass))
#non-significant! so alternatively using constant ratio
summary((lm(anpp_a~-1+tnpp_a,data=grassland_sitemean))) # 0.49 for anpp, so 0.51 for bnpp

#leaf Nmass
###2. leaf Nmass basing on a site-species model
SP_input <- read.csv(file="~/data/leaf_traits/combined_leaf_traits_updated.csv") #new one
summary(SP_input$C_percent)
#median of value of leaf C = 0.47 or 47%

#remove bahar, as repeated to atkin - and filter a few points with vcmax25<0
SP_input <- subset(SP_input,source!="Bahar et al 2017 New Phytologist")
SP_input <- subset(SP_input,Vcmax25>0)
SP_input$Vcmax.25 <- SP_input$Vcmax25
SP_input$Elevation <- SP_input$z

SP_input2 <- SP_input[,c("lat","lon","Elevation","Vcmax.25","narea","lma")]
sitemean <- aggregate(SP_input2,by=list(SP_input2$lon,SP_input2$lat), FUN=mean, na.rm=TRUE) #site-mean

sitemean$sitename <- paste0("s", 1:nrow(sitemean),sep="") # define sitename (s1,s2..s276)

SP_site1 <- sitemean[,c("lon","lat","sitename")]
SP_final1 <- merge(SP_input,SP_site1,by=c("lat","lon"),all.x=TRUE) #merged sitename to SP data

SP_Vcmax.25 <- aggregate(Vcmax.25~sitename+species,SP_final1,mean) #umol/m2/s
SP_Elevation <- aggregate(Elevation~sitename+species,SP_final1,mean)
SP_narea<- aggregate(narea~sitename+species,SP_final1,mean) # g/m2
SP_lma<- aggregate(lma~sitename+species,SP_final1,mean) # g/m2
SP_lat<- aggregate(lat~sitename+species,SP_final1,mean)
SP_lon<- aggregate(lon~sitename+species,SP_final1,mean)

#merging all observed traits in a site-species dataset.
sitespecies_final <-Reduce(function(x,y) merge(x = x, y = y, by = c("sitename","species"),all.x=TRUE), 
                           list(SP_lon,SP_lat,SP_Elevation,SP_Vcmax.25,
                                SP_narea,SP_lma))

#obtain Nrubisco and Nstructural from this large dataset
#firstly, for site-species data
nmass_a <- sitespecies_final$narea/sitespecies_final$lma
vcmax25_lma_a <- sitespecies_final$Vcmax.25/sitespecies_final$lma
sitename_a <- sitespecies_final$sitename
species_a <- sitespecies_final$species

hist(sitespecies_final$narea) # g/m2
hist(sitespecies_final$lma) # g/m2
hist(sitespecies_final$Vcmax.25) # umol/m2/s

#Fit (Nmass) ~ Ns + Nr * (Vcmax25/LMA) - for site-species data
n1 <- lmer(nmass_a~vcmax25_lma_a + (1|sitename_a)+(1|species_a))
summary(n1)
r.squaredGLMM(n1)

#validation directly

sitemean$pred_nmass <- (summary(n1)$coefficients[1,1]) + (summary(n1)$coefficients[2,1])* sitemean$Vcmax.25/sitemean$lma
sitemean$obs_nmass <- sitemean$narea/sitemean$lma

p11 <- analyse_modobs2(sitemean,"pred_nmass","obs_nmass", type = "points",relative=TRUE)$gg + larger_size+coord_obs_pred()+
  labs(y = ~paste("Forest leaf ", N[obs.], " (g g"^-1,")")) +labs(x = ~paste("Forest leaf ", N[pred.], " (g g"^-1,")"))

length(na.omit(sitemean$obs_nmass))

###3. NRE model basing site-mean (lm)
NRE_climate <- read.csv("~/data/NRE_various/NRE_dataset.csv")
NRE_climate$nre_a <- log(NRE_climate$nre/(1-NRE_climate$nre))
NRE_climate$Tg_a <- NRE_climate$Tg
NRE_climate$vpd_a <- log(NRE_climate$vpd)
NRE_climate$PPFD_a <- log(NRE_climate$PPFD)
NRE_climate$ndep_a <- log(NRE_climate$ndep)
NRE_climate2 <- na.omit(NRE_climate[,c("nre_a","vpd_a","Tg_a","PPFD_a","ndep_a")])
stepwise_lm(NRE_climate2,"nre_a")[[1]]
nre_model <- lm(nre_a~Tg_a+vpd_a,data=NRE_climate2)
summary(nre_model)
r.squaredGLMM(nre_model)
length(na.omit(NRE_climate$nre_a))
#validation directly
NRE_climate$pred_nre <- (1/(1+exp(-(summary(nre_model)$coefficients[1,1] + summary(nre_model)$coefficients[2,1] *NRE_climate$Tg_a + 
                                      summary(nre_model)$coefficients[3,1] * NRE_climate$vpd_a))))

p12 <- analyse_modobs2(NRE_climate,"pred_nre","nre", type = "points",relative=TRUE)$gg + larger_size+coord_obs_pred()+
  labs(y = ~paste("Forest ", NRE[obs.])) +labs(x = ~paste("Forest ", NRE[pred.]))


# Cmass constant = 47%
SP_input <- read.csv(file="~/data/leaf_traits/combined_leaf_traits_updated.csv") #new one
SP_input_mean <- aggregate(SP_input,by=list(SP_input$lon,SP_input$lat), FUN=mean, na.rm=TRUE)
summary(SP_input_mean$C_percent,na.rm=TRUE)
dim(subset(SP_input_mean,C_percent>0))

# Root C/N = 94 - as derived from median values of collected samples
NPP_forest_sitemean <- aggregate(NPP_forest,by=list(NPP_forest$site), FUN=mean, na.rm=TRUE) 
summary(NPP_forest_sitemean$CN_root_final) # using median = 94
dim(subset(NPP_forest,CN_root_final>0))

# Wood C/N = 100 - as derived median values of TRY database
wood_cn <- read.csv("~/data/CN_wood/wood_cn.csv")
wood_cn_sitemean <- aggregate(wood_cn,by=list(wood_cn$lon,wood_cn$lat), FUN=mean, na.rm=TRUE) 
summary(wood_cn_sitemean$OrigValueStr)
dim(wood_cn_sitemean)

##2. For Grassland
#leaf c/n model. median = 18. 
#root c/n model median = 41.
#see line 887-888 in Forest_site_org, all derived from "TianDi Grassland", but its BP, ANPP, BNPP data not used in our study
#summary(aggregate(subset(dataset6,pft=="Grassland"),by=list(subset(dataset6,pft=="Grassland")$site), FUN=mean, na.rm=TRUE)$CN_leaf_final)
#summary(aggregate(subset(dataset6,pft=="Grassland"),by=list(subset(dataset6,pft=="Grassland")$site), FUN=mean, na.rm=TRUE)$CN_root_final)

#NRE = 69%
NRE_Du <- read.csv(file="~/data/NRE_various/NRE_Du/NRE_Du.csv")
NRE_Dong <- read.csv(file="~/data/NRE_various/NRE_Deng/NRE_Deng.csv")

#first - make forest model only
NRE_Du_df <- NRE_Du[,c("lon","lat","NRE","MAT","MAP","VegeType")]
#vegetype =1 is woody ecosystem (all assumed as forest here)
#vegetype = 2 is grassland ecosystem
NRE_Du_df <- subset(NRE_Du_df,VegeType==2)
NRE_Du_df <- aggregate(NRE_Du_df,by=list(NRE_Du_df$lon,NRE_Du_df$lat), FUN=mean, na.rm=TRUE) #site-mean
NRE_Du_df <- NRE_Du_df[,c(3:7)]
head(NRE_Du_df)
dim(NRE_Du_df)

NRE_Dong_df <- NRE_Dong[,c("Longitude","Latitude","NRE.nitrogen.resorption.efficiency.","MAT","MAP","Biome.abbreviation...")]
names(NRE_Dong_df) <- c("lon","lat","NRE","MAT","MAP","biome")
NRE_Dong_df <- subset(NRE_Dong_df,biome=="Grs")
#Forest is TRF,STF,TF,BF
#Desert is Des; Tundra is TUN
#Grassland is Grs
NRE_Dong_df <- aggregate(NRE_Dong_df,by=list(NRE_Dong_df$lon,NRE_Dong_df$lat), FUN=mean, na.rm=TRUE) #site-mean
head(NRE_Dong_df)
NRE_Dong_df <- NRE_Dong_df[,c(3:7)]
dim(NRE_Dong_df)

NRE_Dong_df$source <- "Dong"
NRE_Du_df$source <- "Du"
NRE_df <- rbind(NRE_Du_df,NRE_Dong_df)
summary(NRE_df$NRE)
length(NRE_df$NRE)
#median of NRE in grassland is 0.69, site  =26

#final model look
#bp_model,anpp_tnpp_model,anpp_leafnpp_model,bp_grass_model,0.49/0.51,n1,nre_model

#check model vif
#vif_bp,vif_anpp_tnpp,vif_anpp_leafnpp,vif_bp_grass

b1 <- as.ggplot(~barplot(vif_bp, main = "VIF of Forest BP model", horiz = TRUE, col = "steelblue",
                         names.arg = c("Tg", "fAPAR", "ln PPFD", "ln soil C/N", "ln age", "ln vpd")))
b2 <- as.ggplot(~barplot(vif_anpp_tnpp, main = "VIF of Forest ANPP/BP model", horiz = TRUE, col = "steelblue",
                         names.arg = c("Tg", "fAPAR", "ln PPFD", "ln soil C/N", "ln age", "ln vpd")))
b3 <- as.ggplot(~barplot(vif_anpp_leafnpp, main = "VIF of Forest leaf-NPP/ANPP model", horiz = TRUE, col = "steelblue",
                         names.arg = c("Tg", "fAPAR", "ln PPFD", "ln soil C/N", "ln age", "ln vpd")))
b4 <- as.ggplot(~barplot(vif_bp_grass, main = "VIF of Grassland BP model", horiz = TRUE, col = "steelblue",
                         names.arg = c("Tg", "ln PPFD", "ln vpd", "ln soil C/N", "fAPAR")))
plot_grid(b1,b2,b3,b4,
          labels = c('(a)','(b)','(c)','(d)'),
          ncol=2,label_x = 0.9,label_y=0.92)+white
ggsave(paste("~/data/output/newphy_vif_figs.jpg",sep=""),width = 10, height = 13)

#forest validation
NPP_forest$pred_npp <- summary(bp_model)$coefficients[1,1] +  
  summary(bp_model)$coefficients[2,1] * NPP_forest$Tg_a +
  summary(bp_model)$coefficients[3,1] * NPP_forest$fAPAR_a +
  summary(bp_model)$coefficients[4,1] * NPP_forest$PPFD_a +
  summary(bp_model)$coefficients[5,1] * NPP_forest$CNrt_a+
  summary(bp_model)$coefficients[6,1] * NPP_forest$age_a

NPP_forest$pred_anpp <- NPP_forest$pred_npp * 
  (1/(1 + exp(-(summary(anpp_tnpp_model)$coefficients[1,1]+
                  summary(anpp_tnpp_model)$coefficients[2,1] * NPP_forest$CNrt_a +
                  summary(anpp_tnpp_model)$coefficients[3,1] * NPP_forest$PPFD_a + 
                  summary(anpp_tnpp_model)$coefficients[4,1] * NPP_forest$Tg_a+
                  summary(anpp_tnpp_model)$coefficients[5,1] * NPP_forest$age_a))))

NPP_forest$pred_bnpp <- NPP_forest$pred_npp - NPP_forest$pred_anpp

NPP_forest$pred_lnpp <- NPP_forest$pred_anpp *
  (1/(1 + exp(-(summary(anpp_leafnpp_model)$coefficients[1,1]+
                  summary(anpp_leafnpp_model)$coefficients[2,1]* NPP_forest$fAPAR_a +
                  summary(anpp_leafnpp_model)$coefficients[3,1] * NPP_forest$vpd_a + 
                  summary(anpp_leafnpp_model)$coefficients[4,1] * NPP_forest$PPFD_a))))

NPP_forest$pred_wnpp <- NPP_forest$pred_anpp - NPP_forest$pred_lnpp

NPP_forest$pred_leafnc <- (summary(n1)$coefficients[1,1]/0.47) + 
  (summary(n1)$coefficients[2,1]/0.47) * NPP_forest$vcmax25/NPP_forest$LMA

NPP_forest$pred_lnf <- NPP_forest$pred_lnpp*NPP_forest$pred_leafnc
NPP_forest$pred_nre <- (1/(1+exp(-(summary(nre_model)$coefficients[1,1] +
                                     summary(nre_model)$coefficients[2,1] *NPP_forest$Tg_a +
                                     summary(nre_model)$coefficients[3,1] * NPP_forest$vpd_a))))

NPP_forest$pred_wnf <- NPP_forest$pred_wnpp/100
NPP_forest$pred_bnf <- NPP_forest$pred_bnpp/94

NPP_forest$pred_nuptake <- NPP_forest$pred_lnf*(1-NPP_forest$pred_nre)+NPP_forest$pred_wnf+NPP_forest$pred_bnf

NPP_forest$lnf_obs_final <-NPP_forest$NPP.foliage/NPP_forest$CN_leaf_final
NPP_forest$bnf_obs_final  <- NPP_forest$BNPP_1/NPP_forest$CN_root_final
NPP_forest$wnf_obs_final  <- NPP_forest$NPP.wood/NPP_forest$CN_wood_final

#aggregate

p1 <- analyse_modobs2(NPP_forest, "pred_npp","TNPP_1",type = "points",relative=TRUE)$gg + larger_size+coord_obs_pred()+
  labs(y = ~paste("Forest ", BP[obs.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("Forest ", BP[pred.], " (gC m"^-2,"yr"^-1,")"))

p2 <- analyse_modobs2(NPP_forest,"pred_anpp", "ANPP_2",type = "points",relative=TRUE)$gg+ larger_size+coord_obs_pred()+
  labs(y = ~paste("Forest ", ANPP[obs.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("Forest ", ANPP[pred.], " (gC m"^-2,"yr"^-1,")"))

p3 <- analyse_modobs2(NPP_forest,"pred_lnpp","NPP.foliage", type = "points",relative=TRUE)$gg + larger_size+coord_obs_pred()+
  labs(y = ~paste("Forest leaf ", NPP[obs.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("Forest leaf ", NPP[pred.], " (gC m"^-2,"yr"^-1,")"))

p4 <- analyse_modobs2(NPP_forest, "pred_wnpp","NPP.wood",type = "points",relative=TRUE)$gg+ larger_size+coord_obs_pred()+
  labs(y = ~paste("Forest wood ", NPP[obs.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("Forest wood ", NPP[pred.], " (gC m"^-2,"yr"^-1,")"))

p5 <- analyse_modobs2(NPP_forest,"pred_bnpp","BNPP_1", type = "points",relative=TRUE)$gg+ larger_size+coord_obs_pred()+
  labs(y = ~paste("Forest ", BNPP[obs.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("Forest ", BNPP[pred.], " (gC m"^-2,"yr"^-1,")"))

p6 <- analyse_modobs2(NPP_forest,"pred_lnf","lnf_obs_final", type = "points",relative=TRUE)$gg+ larger_size+coord_obs_pred()+
  labs(y = ~paste("Forest leaf N ", flux[obs.], " (gN m"^-2,"yr"^-1,")")) +labs(x = ~paste("Forest leaf N ", flux[pred.], " (gN m"^-2,"yr"^-1,")"))

p7 <- analyse_modobs2(NPP_forest,"pred_nuptake","Nmin", type = "points",relative=TRUE)$gg+ larger_size+coord_obs_pred()+
  labs(y = ~paste("Net ", minerlization[obs.], " (gN m"^-2,"yr"^-1,")")) +labs(x = ~paste("Forest N ", uptake[pred.], " (gN m"^-2,"yr"^-1,")"))

#check how many samples = 225, sites = 84
dim(subset(NPP_forest,is.na(pred_nuptake)==F & is.na(Nmin)==F))
dim(unique(subset(NPP_forest,pred_nuptake>0 & Nmin>0)[c("lon","lat")]))

#grass validation
NPP_grassland$grassland_pred_npp <- summary(bp_grass_model)$coefficients[1,1] +  
  summary(bp_grass_model)$coefficients[2,1] * NPP_grassland$PPFD_a +
  summary(bp_grass_model)$coefficients[3,1] * NPP_grassland$Tg_a 

NPP_grassland$grassland_pred_anpp <- NPP_grassland$grassland_pred_npp * 0.49
NPP_grassland$grassland_pred_lnf <- NPP_grassland$grassland_pred_anpp/18
NPP_grassland$grassland_pred_bnpp <- NPP_grassland$grassland_pred_npp - NPP_grassland$grassland_pred_anpp
NPP_grassland$grassland_pred_bnf <- NPP_grassland$grassland_pred_bnpp/41
NPP_grassland$grassland_pred_nuptake <- NPP_grassland$grassland_pred_lnf*(1-0.69)+NPP_grassland$grassland_pred_bnf

p8 <- analyse_modobs2(NPP_grassland,"grassland_pred_npp","TNPP_1", type = "points",relative=TRUE)$gg + larger_size+coord_obs_pred()+
  labs(y = ~paste("Grassland ", BP[obs.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("Grassland ", BP[pred.], " (gC m"^-2,"yr"^-1,")"))

p9 <- analyse_modobs2(NPP_grassland,"grassland_pred_anpp","ANPP_2", type = "points",relative=TRUE)$gg + larger_size+coord_obs_pred()+
  labs(y = ~paste("Grassland ", ANPP[obs.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("Grassland ", ANPP[pred.], " (gC m"^-2,"yr"^-1,")"))

p10 <- analyse_modobs2(NPP_grassland,"grassland_pred_bnpp","BNPP_1", type = "points",relative=TRUE)$gg + larger_size+coord_obs_pred()+
  labs(y = ~paste("Grassland ", BNPP[obs.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("Grassland ", BNPP[pred.], " (gC m"^-2,"yr"^-1,")"))

white <- theme(plot.background=element_rect(fill="white", color="white"))

#check mabe
mean(abs((NPP_forest$TNPP_1-NPP_forest$pred_npp)/NPP_forest$TNPP_1),na.rm=TRUE) * 100
mean(abs((NPP_forest$ANPP_2-NPP_forest$pred_anpp)/NPP_forest$ANPP_2),na.rm=TRUE) * 100
mean(abs((NPP_forest$NPP.foliage-NPP_forest$pred_lnpp)/NPP_forest$NPP.foliage),na.rm=TRUE) * 100
mean(abs((NPP_forest$NPP.wood-NPP_forest$pred_wnpp)/NPP_forest$NPP.wood),na.rm=TRUE) * 100
mean(abs((NPP_forest$BNPP_1-NPP_forest$pred_bnpp)/NPP_forest$BNPP_1),na.rm=TRUE) * 100
mean(abs((NPP_forest$lnf_obs_final-NPP_forest$pred_lnf)/NPP_forest$lnf_obs_final),na.rm=TRUE) * 100
mean(abs((NPP_forest$Nmin-NPP_forest$pred_nuptake)/NPP_forest$Nmin),na.rm=TRUE) * 100
mean(abs((NPP_grassland$TNPP_1-NPP_grassland$grassland_pred_npp)/NPP_grassland$TNPP_1),na.rm=TRUE) * 100
mean(abs((NPP_grassland$ANPP_2-NPP_grassland$grassland_pred_anpp)/NPP_grassland$ANPP_2),na.rm=TRUE) * 100
mean(abs((NPP_grassland$BNPP_1-NPP_grassland$grassland_pred_bnpp)/NPP_grassland$BNPP_1),na.rm=TRUE) * 100
mean(abs((NRE_climate$nre-NRE_climate$pred_nre)/NRE_climate$nre),na.rm=TRUE) * 100
mean(abs((sitemean$obs_nmass-sitemean$pred_nmass)/sitemean$obs_nmass),na.rm=TRUE) * 100

#fig.2 validation
plot_grid(p1,p2,p5,
          p8,p9,p10,
          p3,p4,p11,
          p6,p12,p7, 
          labels = c('(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)'),
          ncol=3,label_x = 0.9,label_y=0.92,label_size = 25)+white
ggsave(paste("~/data/output/newphy_fig2.jpg",sep=""),width = 24, height = 25)

#now, inputting all predictors
vcmax25_df <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/vcmax25.nc"),
  varnam = "vcmax25"))

Tg <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/Tg.nc"),
  varnam = "Tg"))

PPFD <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/PPFD.nc"),
  varnam = "PPFD"))

vpd <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/vpd.nc"),
  varnam = "vpd"))

fAPAR <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/fAPAR.nc"),
  varnam = "fAPAR"))

age <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/age.nc"),
  varnam = "age"))

CNrt <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/CNrt.nc"),
  varnam = "CNrt"))

LMA <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/LMA.nc"),
  varnam = "LMA"))

###input land cover
ncin <- nc_open("~/data/landcover/modis_landcover_halfdeg_2010_FILLED.nc")
lon <- ncvar_get(ncin,"lon")
nlon <- dim(lon) 
lat<-ncvar_get(ncin,"lat")
nlat<-dim(lat)
pftcover <-ncvar_get(ncin,"pftcover")
nc_close(ncin)
pftcover_long <- as.vector(pftcover)
pftcover <- as.data.frame(matrix(pftcover_long, nrow = nlon * nlat, ncol = 10))
#see get_fpc_grid function: https://github.com/stineb/sofun/blob/db7a9e8e486f576fd7b9f1f74edb1df7a8d2c4f7/src/forcing_global_wmodel.mod.f90 
#it clarified that: 1-6 is forest, 8 is grassland
forest_percent <- rowSums(pftcover[,1:6],na.rm=TRUE)/rowSums(pftcover[,c(1:6,8)],na.rm = TRUE)
grass_percent <- pftcover[,8]/rowSums(pftcover[,c(1:6,8)],na.rm = TRUE)
summary(grass_percent + forest_percent) # check - their sum = 1, perfect!

###3. calculate weighted-sum
#firstly - filter na points - so that all output map has same numbers of NA.
all_predictors <- as.data.frame(cbind(Tg$myvar,PPFD$myvar,vpd$myvar,
                                      fAPAR$myvar,age$myvar,
                                      CNrt$myvar,LMA$myvar,vcmax25_df$myvar))
all_predictors$available_grid = rowMeans(all_predictors)
#just to find all na columns
all_predictors$available_grid[is.na(all_predictors$available_grid)==FALSE] <- 1
summary(all_predictors$available_grid)
available_grid2 <- all_predictors$available_grid

#represent grids when stand-age is especially in NA, but others are fine
names(all_predictors) <- c("Tg","PPFD","vpd","fAPAR","age","CNrt","LMA","vcmax25","available_grid")
all_predictors$lon <- vcmax25_df$lon
all_predictors$lat <- vcmax25_df$lat
summary(all_predictors)

#final calculation - now divide into forest, grassland and pft
#available_grid2 here was used as a list of data to identify if a grid is available (=1) or any prediction fields shown as NA 

#bp_model,anpp_tnpp_model,anpp_leafnpp_model,bp_grass_model,0.49/0.51,n1,nre_model

npp_f <- summary(bp_model)$coefficients[1,1] +  
  summary(bp_model)$coefficients[2,1] * Tg$myvar +
  summary(bp_model)$coefficients[3,1] * fAPAR$myvar +
  summary(bp_model)$coefficients[4,1] * log(PPFD$myvar) +
  summary(bp_model)$coefficients[5,1] * log(CNrt$myvar)+
  summary(bp_model)$coefficients[6,1] * log(age$myvar)

#show npp<=0 points - high lat!
#newmap <- getMap(resolution = "low")
#plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)
#points(subset(as.data.frame(cbind(vcmax25_df,npp_f)),npp_f<=0)$lon,subset(as.data.frame(cbind(vcmax25_df,npp_f)),npp_f<=0)$lat, col="red", pch=16,cex=1)

#see how many girds are NA - 0.3%
(length(npp_f[npp_f<=0])-length(npp_f[is.na(npp_f)==TRUE]))/(length(npp_f[npp_f>0])-length(npp_f[is.na(npp_f)==TRUE]))

#convert to NA!
npp_f[npp_f<=0] <- NA

anpp_f <- npp_f * 
  (1/(1 + exp(-(summary(anpp_tnpp_model)$coefficients[1,1]+
                  summary(anpp_tnpp_model)$coefficients[2,1] * log(CNrt$myvar) +
                  summary(anpp_tnpp_model)$coefficients[3,1] * log(PPFD$myvar) + 
                  summary(anpp_tnpp_model)$coefficients[4,1] * Tg$myvar+
                  summary(anpp_tnpp_model)$coefficients[5,1] * log(age$myvar)))))

bnpp_f <- npp_f - anpp_f

lnpp_f <- anpp_f * (1/(1 + exp(-(summary(anpp_leafnpp_model)$coefficients[1,1]+
                                   summary(anpp_leafnpp_model)$coefficients[2,1]* fAPAR$myvar +
                                   summary(anpp_leafnpp_model)$coefficients[3,1] * log(vpd$myvar) + 
                                   summary(anpp_leafnpp_model)$coefficients[4,1] * log(PPFD$myvar)))))

wnpp_f <- anpp_f - lnpp_f

leafnc_f <- (summary(n1)$coefficients[1,1]/0.47) + 
  (summary(n1)$coefficients[2,1]/0.47) * vcmax25_df$myvar/LMA$myvar

nre_f <- (1/(1+exp(-(summary(nre_model)$coefficients[1,1] +
                       summary(nre_model)$coefficients[2,1] *Tg$myvar +
                       summary(nre_model)$coefficients[3,1] * log(vpd$myvar)))))

lnf_f <- (1-nre_f)* leafnc_f * lnpp_f
wnf_f <- wnpp_f/100
#100 is constant wood c/n
bnf_f <- bnpp_f/94
#94 is constant root c/n
nuptake_f <- lnf_f + wnf_f + bnf_f

#LMG test in forest
df_all <- (as.data.frame(cbind(nuptake_f,npp_f,all_predictors)))
df_all$nue <- df_all$npp_f/df_all$nuptake_f

#key conponents
df_all$anpp_bp <- anpp_f/npp_f
df_all$leafnpp_anpp <- lnpp_f/anpp_f
df_all$leafcn <- leafnc_f
df_all$leafnre <- nre_f
head(df_all)
df_all <- na.omit(df_all)

gpa.model<-lm(nuptake_f~npp_f+anpp_bp+leafnpp_anpp+leafcn+leafnre, data=df_all)
calc.relimp(gpa.model, rela=TRUE)

gpa.model<-lm(nuptake_f~Tg+PPFD+vpd+fAPAR+age+CNrt+LMA+vcmax25, data=df_all)
calc.relimp(gpa.model, rela=TRUE)

gpa.model<-lm(nue~npp_f+anpp_bp+leafnpp_anpp+leafcn+leafnre, data=df_all)
calc.relimp(gpa.model, rela=TRUE)

gpa.model<-lm(nue~Tg+PPFD+vpd+fAPAR+age+CNrt+LMA+vcmax25, data=df_all)
calc.relimp(gpa.model, rela=TRUE)



#grass
npp_g <- summary(bp_grass_model)$coefficients[1,1] +  
  summary(bp_grass_model)$coefficients[2,1] * log(PPFD$myvar) +
  summary(bp_grass_model)$coefficients[3,1] * Tg$myvar 

#show npp<=0 points - high lat!
#newmap <- getMap(resolution = "low")
#plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)
#points(subset(as.data.frame(cbind(vcmax25_df,npp_g)),npp_g<=0)$lon,subset(as.data.frame(cbind(vcmax25_df,npp_g)),npp_g<=0)$lat, col="red", pch=16,cex=1)

#test how many are NA
(length(npp_g[npp_g<=0])-length(npp_g[is.na(npp_g)==TRUE]))/(length(npp_g[npp_g>0])-length(npp_g[is.na(npp_g)==TRUE]))

#convert to NA!
npp_g[npp_g<=0] <- NA

summary(npp_g)

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
nuptake_pft_final <- nuptake_pft
nuptake_forest <- available_grid2* (nuptake_f*forest_percent)
nuptake_grass <- available_grid2* (nuptake_g*grass_percent)

all_maps <- as.data.frame(cbind(vcmax25_df,npp_pft,npp_forest,npp_grass,
                                anpp_pft,anpp_forest,anpp_grass,
                                bnpp_pft,bnpp_forest,bnpp_grass,
                                lnpp_forest,wnpp_forest,wnf_forest,
                                leafcn_pft,leafcn_forest,leafcn_grassland,
                                nre_pft,nre_forest,nre_grassland,
                                lnf_pft,lnf_forest,lnf_grass,
                                bnf_pft,bnf_forest,bnf_grass,
                                nuptake_pft,nuptake_forest,nuptake_grass))

summary(all_maps)

#####area_m2 to show each grid's area in m2
calc_area <- function( lat, dx=1, dy=1 ){
  r_earth <- 6370499.317638  # to be consistent with how Ferret calculates areas of spheres (https://www.pmel.noaa.gov/maillists/tmap/ferret_users/fu_2016/msg00155.html)
  area <- 4 * r_earth^2 * 0.5 * dx * pi/180 * cos( abs(lat) * pi/180 ) * sin( 0.5 * dy * pi/180 )
  return(area)
}
lonlat <- vcmax25_df[,c("lon","lat")]
area_m2 <- calc_area(lonlat$lat,0.5,0.5)
#fland - to show each grid's land cover percentage
nc <- read_nc_onefile("~/data/fland/global.fland.nc") #Input nc
output_fland <- nc_to_df(nc, varnam = "fland")
fland <- output_fland$myvar
#include conversion factor (from g to Pg)
conversion <- area_m2 * fland /1e+15

#Fig.3 global maps
white <- theme(plot.background=element_rect(fill="white", color="white"))
npp_f1 <- (NPP_forest %>% filter(TNPP_1>0) %>% filter(pred_npp>0))[,c("lon","lat")]
npp_g1 <- (NPP_grassland %>% filter(TNPP_1>0) %>% filter(grassland_pred_npp>0))[,c("lon","lat")]
total_value <- round(sum(all_maps[,"npp_pft"]*conversion,na.rm=TRUE),2)

gg <- plot_map3(na.omit(all_maps[,c("lon","lat","npp_pft")]),
                varnam = "npp_pft",latmin = -65, latmax = 85,combine=FALSE)

a3 <- gg$ggmap +
  geom_point(data=npp_f1,aes(lon,lat),col="red",size=1.5)+
  geom_point(data=npp_g1,aes(lon,lat),col="blue",size=1.5)+
  labs(title = paste("BP:", total_value, "PgC/yr", sep=" " ))+
  theme_grey(base_size = 12)+ white

a4 <- gg$gglegend+labs(title = ~paste("gC m"^-2,"yr"^-1))+ white

#3.anpp
anpp_f1 <- (NPP_forest %>% filter(ANPP_2>0) %>% filter(pred_anpp>0))[,c("lon","lat")]
anpp_g1 <- (NPP_grassland %>% filter(ANPP_2>0) %>% filter(grassland_pred_anpp>0))[,c("lon","lat")]
total_value <- round(sum(all_maps[,"anpp_pft"]*conversion,na.rm=TRUE),2)

gg <- plot_map3(na.omit(all_maps[,c("lon","lat","anpp_pft")]),
                varnam = "anpp_pft",latmin = -65, latmax = 85,combine=FALSE)

a5 <- gg$ggmap +
  geom_point(data=anpp_f1,aes(lon,lat),col="red",size=1.5)+
  geom_point(data=anpp_g1,aes(lon,lat),col="blue",size=1.5)+
  labs(title = paste("ANPP:", total_value, "PgC/yr", sep=" " ))+
  theme_grey(base_size = 12)+ white

a6 <- gg$gglegend+labs(title = ~paste("gC m"^-2,"yr"^-1))+ white

#4. leaf c/n
SP_input <- read.csv(file="~/data/leaf_traits/combined_leaf_traits_updated.csv") #new one 
SP_input2 <- SP_input[,c("lat","lon","z","Vcmax25","narea","lma")]
sitemean <- aggregate(SP_input2,by=list(SP_input2$lon,SP_input2$lat), FUN=mean, na.rm=TRUE) 
sitemean$pred_leafn <- (summary(n1)$coefficients[1,1]) + (summary(n1)$coefficients[2,1])* sitemean$Vcmax25/sitemean$lma
sitemean$obs_leafn <- sitemean$narea/sitemean$lma

laefcn <- (sitemean %>% filter(pred_leafn>0) %>% filter(obs_leafn>0))[,c("lon","lat")]

gg <- plot_map3(na.omit(all_maps[,c("lon","lat","leafcn_pft")]),
                varnam = "leafcn_pft",latmin = -65, latmax = 85,combine=FALSE)

total_value <- round(mean(leafcn_pft,na.rm=TRUE),2)
a7 <- gg$ggmap +
  geom_point(data=laefcn,aes(lon,lat),col="red",size=1.5)+
  labs(title = paste("Leaf C/N: ", total_value))+
  theme_grey(base_size = 12)+ white

a8 <- gg$gglegend+ white

#5. NRE
nre_site <- (NRE_climate %>% filter(pred_nre>0) %>% filter(NRE>0))[,c("lon","lat")]

gg <- plot_map3(na.omit(all_maps[,c("lon","lat","nre_pft")]),
                varnam = "nre_pft",latmin = -65, latmax = 85,combine=FALSE)

total_value <- round(mean(nre_pft,na.rm=TRUE),2)

a9 <- gg$ggmap +
  geom_point(data=nre_site,aes(lon,lat),col="red",size=1.5)+
  labs(title = paste("NRE: ", total_value))+
  theme_grey(base_size = 12)+ white

a10 <- gg$gglegend+ white

#6. nuptake
nuptake_f1 <- (NPP_forest %>% filter(pred_nuptake>0) %>% filter(Nmin>0))[,c("lon","lat")]

total_value <- 1000*round(sum(all_maps[,"nuptake_pft"]*conversion,na.rm=TRUE),2) #unit convert from PgN/yr to TgN/yr

gg <- plot_map3(na.omit(all_maps[,c("lon","lat","nuptake_pft")]),
                varnam = "nuptake_pft",latmin = -65, latmax = 85,combine=FALSE)

a11 <- gg$ggmap +
  geom_point(data=nuptake_f1,aes(lon,lat),col="red",size=1.5)+
  labs(title = paste("N uptake: ", total_value, "TgN/yr", sep=" " ))+
  theme_grey(base_size = 12)+ white

a12 <- gg$gglegend+labs(title = ~paste("gN m"^-2,"yr"^-1))+ white

#NUE
all_maps$NUE <- all_maps$npp_pft/all_maps$nuptake_pft
summary(all_maps$NUE)
gg <- plot_map3(na.omit(all_maps[,c("lon","lat","NUE")]),
                varnam = "NUE",latmin = -65, latmax = 85,combine=FALSE)

total_value<- round(sum(all_maps[,"npp_pft"]*conversion,na.rm=TRUE)/sum(all_maps[,"nuptake_pft"]*conversion,na.rm=TRUE),2)

a13 <- gg$ggmap +
  labs(title = paste("NUE: ", total_value, " ", sep=" " ))+
  theme_grey(base_size = 12)+ white

a14 <- gg$gglegend+labs(title = ~paste("gN m"^-2,"yr"^-1))+ white

plot_grid(a3,a4,a5,a6,a7,a8,
          a9,a10,a11,a12,a13,a14,
          nrow=2,
          rel_widths = c(3/12, 1/12,3/12,1/12,3/12,1/12),
          labels = c('(a)',' ','(b)',' ','(c)',' ',
                     '(d)',' ','(e)',' ','(f)',' '))

#this is for showing points
ggsave(paste("~/data/output/newphy_fig3_SI.jpg",sep=""),width = 20, height = 10*(2/3))


#preparing one without showing points - in MS
total_value <- round(sum(all_maps[,"npp_pft"]*conversion,na.rm=TRUE),2)
gg <- plot_map3(na.omit(all_maps[,c("lon","lat","npp_pft")]),
                varnam = "npp_pft",latmin = -65, latmax = 85,combine=FALSE)
a3 <- gg$ggmap +
  labs(title = paste("BP:", total_value, "PgC/yr", sep=" " ))+
  theme_grey(base_size = 12)+ white
a4 <- gg$gglegend+labs(title = ~paste("gC m"^-2,"yr"^-1))+ white

#3.anpp
total_value <- round(sum(all_maps[,"anpp_pft"]*conversion,na.rm=TRUE),2)
gg <- plot_map3(na.omit(all_maps[,c("lon","lat","anpp_pft")]),
                varnam = "anpp_pft",latmin = -65, latmax = 85,combine=FALSE)
a5 <- gg$ggmap +
  labs(title = paste("ANPP:", total_value, "PgC/yr", sep=" " ))+
  theme_grey(base_size = 12)+ white
a6 <- gg$gglegend+labs(title = ~paste("gC m"^-2,"yr"^-1))+ white

#4. leaf c/n
gg <- plot_map3(na.omit(all_maps[,c("lon","lat","leafcn_pft")]),
                varnam = "leafcn_pft",latmin = -65, latmax = 85,combine=FALSE)

total_value <- round(mean(leafcn_pft,na.rm=TRUE),2)
a7 <- gg$ggmap +
  labs(title = paste("Leaf C/N: ", total_value))+
  theme_grey(base_size = 12)+ white
a8 <- gg$gglegend+ white

#5. NRE
gg <- plot_map3(na.omit(all_maps[,c("lon","lat","nre_pft")]),
                varnam = "nre_pft",latmin = -65, latmax = 85,combine=FALSE)
total_value <- round(mean(nre_pft,na.rm=TRUE),2)

a9 <- gg$ggmap +
  labs(title = paste("NRE: ", total_value))+
  theme_grey(base_size = 12)+ white
a10 <- gg$gglegend+ white

#6. nuptake
total_value <- 1000*round(sum(all_maps[,"nuptake_pft"]*conversion,na.rm=TRUE),2) #unit convert from PgN/yr to TgN/yr
gg <- plot_map3(na.omit(all_maps[,c("lon","lat","nuptake_pft")]),
                varnam = "nuptake_pft",latmin = -65, latmax = 85,combine=FALSE)
a11 <- gg$ggmap +
  labs(title = paste("N uptake: ", total_value, "TgN/yr", sep=" " ))+
  theme_grey(base_size = 12)+ white

a12 <- gg$gglegend+labs(title = ~paste("gN m"^-2,"yr"^-1))+ white

#NUE
all_maps$NUE <- all_maps$npp_pft/all_maps$nuptake_pft
summary(all_maps$NUE)
gg <- plot_map3(na.omit(all_maps[,c("lon","lat","NUE")]),
                varnam = "NUE",latmin = -65, latmax = 85,combine=FALSE)

total_value<- round(sum(all_maps[,"npp_pft"]*conversion,na.rm=TRUE)/sum(all_maps[,"nuptake_pft"]*conversion,na.rm=TRUE),2)

a13 <- gg$ggmap +
  labs(title = paste("NUE: ", total_value, " ", sep=" " ))+
  theme_grey(base_size = 12)+ white

a14 <- gg$gglegend+labs(title = ~paste("gN m"^-2,"yr"^-1))+ white

plot_grid(a3,a4,a5,a6,a7,a8,
          a9,a10,a11,a12,a13,a14,
          nrow=2,
          rel_widths = c(3/12, 1/12,3/12,1/12,3/12,1/12),
          labels = c('(a)',' ','(b)',' ','(c)',' ',
                     '(d)',' ','(e)',' ','(f)',' '))
ggsave(paste("~/data/output/newphy_fig3.jpg",sep=""),width = 20, height = 10*(2/3))

#work on effect of each factor on NUE
#now, work on NUE
nue_final <- all_maps$NUE # by default
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

nue_standard <- cal_nue(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,
                        age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$myvar)
summary(nue_standard)

mean_Tg <- mean(Tg$myvar[is.na(available_grid2)==FALSE]);mean_Tg
mean_PPFD <- mean(PPFD$myvar[is.na(available_grid2)==FALSE]);mean_PPFD
mean_vpd <- mean(vpd$myvar[is.na(available_grid2)==FALSE]);mean_vpd
mean_fAPAR <- mean(fAPAR$myvar[is.na(available_grid2)==FALSE]);mean_fAPAR
mean_age <- mean(age$myvar[is.na(available_grid2)==FALSE]);mean_age
mean_CNrt <- mean(CNrt$myvar[is.na(available_grid2)==FALSE]);mean_CNrt
mean_LMA <- mean(LMA$myvar[is.na(available_grid2)==FALSE]);mean_LMA
mean_vcmax25 <- mean(vcmax25_df$myvar[is.na(available_grid2)==FALSE]);mean_vcmax25

nue_Tg <- cal_nue(rep(mean_Tg,259200),PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$myvar)
nue_PPFD <- cal_nue(Tg$myvar,rep(mean_PPFD,259200),vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$myvar)
nue_vpd <- cal_nue(Tg$myvar,PPFD$myvar,rep(mean_vpd,259200),fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$myvar)
nue_fAPAR <- cal_nue(Tg$myvar,PPFD$myvar,vpd$myvar,rep(mean_fAPAR,259200),age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$myvar)
nue_age <- cal_nue(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,rep(mean_age,259200),CNrt$myvar,LMA$myvar,vcmax25_df$myvar)
nue_CNrt <- cal_nue(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,rep(mean_CNrt,259200),LMA$myvar,vcmax25_df$myvar)
nue_LMA <- cal_nue(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,rep(mean_LMA,259200),vcmax25_df$myvar)
nue_vcmax25 <- cal_nue(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,rep(mean_vcmax25,259200))

nue_all <- as.data.frame(cbind(all_predictors,nue_Tg,nue_PPFD,nue_vpd,nue_fAPAR,nue_age,nue_CNrt,nue_LMA,nue_vcmax25))
summary(nue_all)
nue_all$NUE <- all_maps$NUE
#show soil C/N effects on NUE
#aggregate
nue_all$name[nue_all$lat< -60] <- "90°S ~ 60°S"
nue_all$name[nue_all$lat >= -60 & nue_all$lat < -30] <- "60°S ~ 30°S"
nue_all$name[nue_all$lat >= -30 & nue_all$lat < 0] <- "30°S ~ 0"
nue_all$name[nue_all$lat >= 0 & nue_all$lat < 30]<- "0 ~ 30°N"
nue_all$name[nue_all$lat >= 30 & nue_all$lat < 60]<- "30°N ~ 60°N"
nue_all$name[nue_all$lat >= 60]<- "60°N ~ 90°N"

nue_all$name <- factor(nue_all$name,levels = c("90°S ~ 60°S","60°S ~ 30°S","30°S ~ 0",
                                               "0 ~ 30°N","30°N ~ 60°N","60°N ~ 90°N"))
nue_all <- subset(nue_all,name!="90°S ~ 60°S")
a1 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = nue_CNrt),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="blue")+ geom_vline(xintercept=0, linetype="dashed")+labs(y= " ", x = "Soil C/N effect") +larger_size
a2 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = nue_age),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="blue")+ geom_vline(xintercept=0, linetype="dashed")+labs(y= " ", x = "Age effect") +larger_size
a3 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = nue_fAPAR),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="blue")+ geom_vline(xintercept=0, linetype="dashed")+labs(y= " ", x = "fAPAR effect") +larger_size
a4 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = nue_Tg),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="blue")+ geom_vline(xintercept=0, linetype="dashed")+labs(y= " ", x = "Tg effect") +larger_size
a5 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = nue_vpd),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="blue")+ geom_vline(xintercept=0, linetype="dashed")+labs(y= " ", x = "D effect") +larger_size
a6 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = nue_PPFD),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="blue")+ geom_vline(xintercept=0, linetype="dashed")+labs(y= " ", x = "PPFD effect") +larger_size
a7 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = nue_vcmax25),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="blue")+ geom_vline(xintercept=0, linetype="dashed")+labs(y= " ", x = "Vcmax25 effect") +larger_size
a8 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = nue_LMA),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="blue")+ geom_vline(xintercept=0, linetype="dashed")+labs(y= " ", x = "LMA effect") +larger_size
a9 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = NUE),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="blue")+labs(y= " ", x = "NUE") +larger_size

aa1 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = CNrt),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="black")+labs(y= " ", x = "Soil C/N") +larger_size
aa2 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = age),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="black")+labs(y= " ", x = "Age") +larger_size
aa3 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = fAPAR),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="black")+labs(y= " ", x = "fAPAR") +larger_size
aa4 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = Tg),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="black")+labs(y= " ", x = "Tg") +larger_size
aa5 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = vpd),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="black")+labs(y= " ", x = "D") +larger_size
aa6 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = PPFD),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="black")+labs(y= " ", x = "PPFD")+larger_size 
aa7 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = vcmax25),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="black")+labs(y= " ", x = "Vcmax25") +larger_size
aa8 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = LMA),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="black")+labs(y= " ", x = "LMA") +larger_size

plot_grid(a1,a2,a3,a4,a5,a6,a7,a8,a9,
          labels = c('(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'),label_size = 15)+white
ggsave(paste("~/data/output/newphy_fig4.jpg",sep=""),width = 22, height = 10)

plot_grid(aa1,aa2,aa3,aa4,aa5,aa6,aa7,aa8,
          labels = c('(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'),label_size = 15)+white
ggsave(paste("~/data/output/newphy_fig4_parallel.jpg",sep=""),width = 22, height = 10)

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_CNrt")]),varnam = "nue_CNrt",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3"),
                breaks = seq(-0.05,0.05,0.01))
g1 <- gg$ggmap +labs(title = "Soil C/N")+theme_grey(base_size = 15);g2 <- gg$gglegend

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_age")]),varnam = "nue_age",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4","royalblue3","wheat","tomato3"),
                breaks = seq(-0.05,0.05,0.01))
g3 <- gg$ggmap +labs(title = "Age")+theme_grey(base_size = 15);g4 <- gg$gglegend+labs(title = ~paste("years"))

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_fAPAR")]),varnam = "nue_fAPAR",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4","royalblue3", "wheat","tomato3"),
                breaks = seq(-0.04,0.04,0.01))
g5 <- gg$ggmap +labs(title = "fAPAR")+theme_grey(base_size = 15);g6 <- gg$gglegend

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_PPFD")]),varnam = "nue_PPFD",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3","tomato4"),
                breaks = seq(-0.08,0.12,0.02))
g7 <- gg$ggmap +labs(title = "PPFD")+theme_grey(base_size = 15);g8 <- gg$gglegend + labs(title =~paste(mu, "mol m"^-2,"s"^-1))

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_Tg")]),varnam = "nue_Tg",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3"),
                breaks = seq(-0.30,0.30,0.10))
g9 <- gg$ggmap +labs(title =~paste(T[g]))+theme_grey(base_size = 15);g10 <- gg$gglegend + labs(title =~paste("\u00B0C"))

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_vpd")]),varnam = "nue_vpd",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3"),
                breaks = seq(-0.20,0.30,0.05))
g11 <- gg$ggmap +labs(title = "D")+theme_grey(base_size = 15);g12 <- gg$gglegend+ labs(title =~paste("kPa"))

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_vcmax25")]), varnam = "nue_vcmax25",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3"),
                breaks = seq(-0.2,0.2,0.04))
g15 <- gg$ggmap +labs(title =~paste(V[cmax25]))+theme_grey(base_size = 15);g16 <- gg$gglegend+ labs(title =~paste(mu, "mol m"^-2,"s"^-1))

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_LMA")]),varnam = "nue_LMA",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3"),
                breaks = seq(-0.05,0.05,0.01))
g17 <- gg$ggmap +labs(title = "LMA")+theme_grey(base_size = 15);g18 <- gg$gglegend+ labs(title =~paste("g ","m"^-2))

plot_grid(g1,g2,g3,g4,g5,g6,
          g7,g8,g9,g10,g11,g12,
          g15,g16,g17,g18,
          nrow=3,rel_widths = c(3/12, 1/12,3/12,1/12,3/12,1/12),
          labels = c('(a)',' ','(b)',' ','(c)',' ',
                     '(d)',' ','(e)',' ','(f)',' ',
                     '(g)',' ','(h)',' '),label_size = 15)+white

ggsave(paste("~/data/output/newphy_fig4_not_used.jpg",sep=""),width = 20, height = 10)


#figs2 representing all predictors
gg <- plot_map3(na.omit(CNrt[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85,combine=FALSE)
d1 <- gg$ggmap +labs(title = "Soil C/N")+theme_grey(base_size = 15)+ white;d2 <- gg$gglegend+ white

gg <- plot_map3(na.omit(age[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85,combine=FALSE)
d3 <- gg$ggmap +labs(title = "Age")+theme_grey(base_size = 15)+ white;d4 <- gg$gglegend+labs(title = ~paste("years"))+ white

gg <- plot_map3(na.omit(fAPAR[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85,combine=FALSE)
d5 <- gg$ggmap +labs(title = "fAPAR")+theme_grey(base_size = 15)+ white;d6 <- gg$gglegend+ white

gg <- plot_map3(na.omit(PPFD[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85,combine=FALSE,breaks = seq(100,800,100))
d7 <- gg$ggmap +labs(title = "PPFD")+theme_grey(base_size = 15)+ white;d8 <- gg$gglegend + labs(title =~paste(mu, "mol m"^-2,"s"^-1))+ white

gg <- plot_map3(na.omit(Tg[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85,combine=FALSE)
d9 <- gg$ggmap +labs(title =~paste(T[g]))+theme_grey(base_size = 15)+ white;d10 <- gg$gglegend + labs(title =~paste("\u00B0C"))+ white

gg <- plot_map3(na.omit(vpd[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85,combine=FALSE)
d11 <- gg$ggmap +labs(title = "D")+theme_grey(base_size = 15)+ white;d12 <- gg$gglegend+ labs(title =~paste("kPa"))+ white

gg <- plot_map3(na.omit(vcmax25_df[,c("lon","lat","myvar")]), varnam = "myvar",latmin = -65, latmax = 85,combine=FALSE)
d13 <- gg$ggmap +labs(title =~paste(V[cmax25]))+theme_grey(base_size = 15)+ white;d14 <- gg$gglegend+ labs(title =~paste(mu, "mol m"^-2,"s"^-1))+ white

gg <- plot_map3(na.omit(LMA[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85,combine=FALSE)
d15 <- gg$ggmap +labs(title = "LMA")+theme_grey(base_size = 15)+ white;d16 <- gg$gglegend+ labs(title =~paste("g ","m"^-2))+ white

plot_grid(d1,d2,d3,d4,d5,d6,
          d7,d8,d9,d10,d11,d12,
          d13,d14,d15,d16,
          nrow=3,rel_widths = c(3/12, 1/12,3/12,1/12,3/12,1/12),
          labels = c('(a)',' ','(b)',' ','(c)',' ',
                     '(d)',' ','(e)',' ','(f)',' ',
                     '(g)',' ','(h)',' '),label_size = 15)+ white

ggsave(paste("~/data/output/newphy_figS2.jpg",sep=""),width = 20, height = 10)

#trendy
#model output
#CABLE-POP
CABLE_GPP <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/CABLE-POP_S2_gpp_ANN_mean.nc"), varnam = "gpp"))
CABLE_NPP <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/CABLE-POP_S2_npp_ANN_mean.nc"), varnam = "npp"))

#CABLE-POP
CLASS_GPP <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/CLASS-CTEM_S2_gpp_ANN_mean.nc"), varnam = "gpp"))
CLASS_NPP <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/CLASS-CTEM_S2_npp_ANN_mean.nc"), varnam = "npp"))

#CLM
#CLM_fNup <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/CLM5.0_S2_fNup_ANN_mean.nc"), varnam = "fNup")) 
#this product is confusing (1) unit needs /1000 to get gn/m2/year? (2)still many values very high
CLM_GPP <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/CLM5.0_S2_gpp_ANN_mean.nc"), varnam = "gpp"))
CLM_GPP$lon[CLM_GPP$lon>180] <- CLM_GPP$lon[CLM_GPP$lon>180]-360
CLM_NPP <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/CLM5.0_S2_npp_ANN_mean.nc"), varnam = "npp"))
CLM_NPP$lon[CLM_NPP$lon>180] <- CLM_NPP$lon[CLM_NPP$lon>180]-360

#ISAM
#(given in kgC m-2 month-1 - they might not need to multiply with 12
#cdo seltimestep,127/156 ISAM_S2_fNup.nc a1.nc
#cdo -O timmean a1.nc ISAM_S2_fNup_ANN_mean.nc

ISAM_fNup <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/ISAM_S2_fNup_ANN_mean.nc"), varnam = "fNup"))
ISAM_fNup$myvar<- ISAM_fNup$myvar*1000*31556952

ISAM_gpp <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/ISAM_S2_gpp_ANN_mean.nc"), varnam = "gpp"))
ISAM_gpp$myvar<- ISAM_gpp$myvar*1000*31556952

ISAM_npp <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/ISAM_S2_npp_ANN_mean.nc"), varnam = "npp"))
ISAM_npp$myvar<- ISAM_npp$myvar*1000*31556952

#ISBA
ncin <- nc_open("~/data/trendy/v8/ISBA-CTRIP_S2_gpp_ANN_mean.nc")
lon <- ncvar_get(ncin,"lon_FULL");nlon <- dim(lon) 
lat <- ncvar_get(ncin,"lat_FULL");nlat <- dim(lat) 
ISBA_GPP <- ncvar_get(ncin,"gpp")
nc_close(ncin)
ISBA_GPP <- as.vector(ISBA_GPP)
lonlat <- expand.grid(lon,lat)
ISBA_GPP <- as.data.frame(cbind(lonlat,ISBA_GPP))
names(ISBA_GPP) <- c("lon","lat","GPP")

ncin <- nc_open("~/data/trendy/v8/ISBA-CTRIP_S2_npp_ANN_mean.nc")
ISBA_NPP <- ncvar_get(ncin,"npp")
nc_close(ncin)
ISBA_NPP <- as.vector(ISBA_NPP)
ISBA_NPP <- as.data.frame(cbind(lonlat,ISBA_NPP))
names(ISBA_NPP) <- c("lon","lat","NPP")

#JSBACH
JSBACH_fNup <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/JSBACH_S2_fNup_ANN_mean.nc"), varnam = "fNup"))
JSBACH_GPP <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/JSBACH_S2_gpp_ANN_mean.nc"), varnam = "gpp"))
JSBACH_NPP <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/JSBACH_S2_npp_ANN_mean.nc"), varnam = "npp"))
plot_map3(na.omit(JSBACH_fNup[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85)

#JULES
# data from output/JULES-ES-1.0/JULES-ES.1p0.vn5.4.50.CRUJRA2.TRENDYv8.365.S2_Monthly_npp.nc
JULES_GPP <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/jule_gpp.nc"), varnam = "gpp"))
JULES_GPP$lon[JULES_GPP$lon>180] <- JULES_GPP$lon[JULES_GPP$lon>180]-360
JULES_NPP <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/jule_npp.nc"), varnam = "npp"))
JULES_NPP$lon[JULES_NPP$lon>180] <- JULES_NPP$lon[JULES_NPP$lon>180]-360
#looks ok 
#plot_map3(na.omit(JULES_GPP[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85)

#LPJ
LPJ_GPP <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/LPJ-GUESS_S2_gpp_ANN_mean.nc"), varnam = "gpp"))
LPJ_NPP <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/LPJ-GUESS_S2_npp_ANN_mean.nc"), varnam = "npp"))

#ORCHIDEE
ORCHIDEE_GPP <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/ORCHIDEE_S2_gpp_ANN_mean.nc"), varnam = "gpp"))
ORCHIDEE_NPP <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/ORCHIDEE_S2_npp_ANN_mean.nc"), varnam = "npp"))

#ORCHIDEE-CNP
nc_flip_lat <- function(nc){
  
  nc$lat <- rev(nc$lat)
  
  # nlat <- length(nc$lat)
  # nc$vars[[1]] <- nc$vars[[1]][,nlat:1]
  
  arr_flip_lat <- function(arr){
    nlat <- dim(arr)[2]
    arr <- arr[,nlat:1]
    return(arr)
  }
  nc$vars <- purrr::map(nc$vars[1], ~arr_flip_lat(.))
  
  return(nc)
}
ORCHICNP_fNup <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/ORCHIDEE-CNP_S2_fNup_ANN_mean.nc"), varnam = "fNup"))
ORCHICNP_GPP <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/ORCHIDEE-CNP_S2_gpp_ANN_mean.nc"), varnam = "gpp"))
ORCHICNP_NPP <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/ORCHIDEE-CNP_S2_npp_ANN_mean.nc"), varnam = "npp"))

#SDGVM
SDGVM_GPP <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/SDGVM_S2_gpp_ANN_mean.nc"), varnam = "gpp"))
SDGVM_NPP <- as.data.frame(nc_to_df(read_nc_onefile("~/data/trendy/v8/SDGVM_S2_npp_ANN_mean.nc"), varnam = "npp"))

#all_maps
not_used <- list(CLASS_GPP,CLASS_NPP,CLM_GPP,CLM_NPP,JSBACH_fNup,JSBACH_GPP,JSBACH_NPP)

raster_CLASS_npp <- raster("~/data/trendy/v8/CLASS-CTEM_S2_npp_ANN_mean.nc")
raster_CLM_npp <- raster("~/data/trendy/v8/CLM5.0_S2_npp_ANN_mean.nc")
raster_JSBACH_npp <- raster("~/data/trendy/v8/JSBACH_S2_npp_ANN_mean.nc")
raster_JSBACH_fNup <- raster("~/data/trendy/v8/JSBACH_S2_fNup_ANN_mean.nc")
raster_LPX_npp <- raster("~/data/trendy/v8/LPX-Bern_S2_npp_ANN_mean.nc")
raster_LPX_fNup <- raster("~/data/trendy/v8/LPX-Bern_S2_fNup_ANN_mean.nc")

allmaps<- list(CABLE_GPP,CABLE_NPP,ISAM_fNup,ISAM_gpp,ISAM_npp,ISBA_GPP,ISBA_NPP,
               JULES_GPP,JULES_NPP,LPJ_GPP,LPJ_NPP,ORCHIDEE_GPP,ORCHIDEE_NPP,
               ORCHICNP_fNup,ORCHICNP_GPP,ORCHICNP_NPP,SDGVM_GPP,SDGVM_NPP)
obj_name <- c("CABLE_GPP","CABLE_NPP","ISAM_fNup","ISAM_gpp","ISAM_npp","ISBA_GPP","ISBA_NPP",
              "JULES_GPP","JULES_NPP","LPJ_GPP","LPJ_NPP","ORCHIDEE_GPP","ORCHIDEE_NPP",
              "ORCHICNP_fNup","ORCHICNP_GPP","ORCHICNP_NPP","SDGVM_GPP","SDGVM_NPP")
#calculate nue 
summary( ISAM_npp$myvar/ISAM_fNup$myvar)
summary(na.omit(ORCHICNP_NPP$myvar/ORCHICNP_fNup$myvar))

#aggregate based on lon and lat firstly
sitemean <- unique(NPP_forest[,c("lon","lat")])
sp_sites <- SpatialPoints(sitemean) # only select lon and lat

sitemean_final <- unique(NPP_forest[,c("lon","lat")])

for (i in c(1:length(allmaps))){
  df1 <- allmaps[[i]]
  names(df1) <- c("lon","lat",obj_name[i])
  coordinates(df1) <- ~lon+lat 
  gridded(df1) <- TRUE
  df1_global <- raster(df1, obj_name[i]) 
  
  sp_sites_new <- raster::extract(df1_global, sp_sites, sp = TRUE) %>% as_tibble() %>% 
    right_join(sitemean, by = c("lon", "lat"))
  sitemean_final[,i+2] <- sp_sites_new[,1]
}

#CLM has some problem here for coordinates - need convert # CLM_NPP$lon[CLM_NPP$lon>180] <- CLM_NPP$lon[CLM_NPP$lon>180]-360
sitemean2 <- sitemean
sitemean2$lon[sitemean2$lon<0] <- sitemean2$lon[sitemean2$lon<0]+360
sp_sites2 <- SpatialPoints(sitemean2) # only select lon and lat
CLM_NPP <- (raster::extract(raster_CLM_npp, sp_sites2, sp = TRUE) %>% as_tibble() %>% right_join(sitemean2, by = c("lon", "lat")))[,1]


CLASS_NPP <- (raster::extract(raster_CLASS_npp, sp_sites, sp = TRUE) %>% as_tibble() %>% right_join(sitemean, by = c("lon", "lat")))[,1]
JSBACH_NPP <- (raster::extract(raster_JSBACH_npp, sp_sites, sp = TRUE) %>% as_tibble() %>% right_join(sitemean, by = c("lon", "lat")))[,1]
JSBACH_fNup <- (raster::extract(raster_JSBACH_fNup, sp_sites, sp = TRUE) %>% as_tibble() %>% right_join(sitemean, by = c("lon", "lat")))[,1]
LPX_NPP <- (raster::extract(raster_LPX_npp, sp_sites, sp = TRUE) %>% as_tibble() %>% right_join(sitemean, by = c("lon", "lat")))[,1]
LPX_fNup <- (raster::extract(raster_LPX_fNup, sp_sites, sp = TRUE) %>% as_tibble() %>% right_join(sitemean, by = c("lon", "lat")))[,1]
new_added <- as.data.frame(cbind(CLASS_NPP,CLM_NPP,JSBACH_NPP,JSBACH_fNup,LPX_NPP,LPX_fNup))
names(new_added) <- c("CLASS_NPP","CLM_NPP","JSBACH_NPP","JSBACH_fNup","LPX_NPP","LPX_fNup")

sitemean_final <- as.data.frame(cbind(sitemean_final,new_added))

summary(sitemean_final)

sitemean_final$CABLE_NPP[sitemean_final$CABLE_NPP==0] <- NA
sitemean_final$ISAM_fNup[sitemean_final$ISAM_fNup==0] <- NA
sitemean_final$ISAM_gpp[sitemean_final$ISAM_gpp==0] <- NA
sitemean_final$ISAM_npp[sitemean_final$ISAM_npp==0] <- NA
sitemean_final$JULES_GPP[sitemean_final$JULES_GPP==0] <- NA
sitemean_final$JULES_NPP[sitemean_final$JULES_NPP==0] <- NA
sitemean_final$LPX_NPP[sitemean_final$LPX_NPP==0] <- NA
sitemean_final$LPX_fNup[sitemean_final$LPX_fNup==0] <- NA

#first, our model - using predicted data
BP_dataset3 <- na.omit(NPP_forest[,c("pred_npp","age_a","fAPAR_a","CNrt_a","Tg_a","PPFD_a","vpd_a","site_a","lon","lat")])
BP_dataset3<- aggregate(BP_dataset3,by=list(BP_dataset3$lon,BP_dataset3$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
bp_model_ours <- (lm(pred_npp~Tg_a+PPFD_a+vpd_a,data=BP_dataset3))
summary(bp_model_ours)

bp1a <- visreg(bp_model_ours,"Tg_a",type="contrast")
bp1b <- visreg(bp_model_ours,"PPFD_a",type="contrast")
bp1c <- visreg(bp_model_ours,"vpd_a",type="contrast")

#trendy model
NPP_statistical <- merge(NPP_forest,sitemean_final,by=c("lon","lat"),all.x=TRUE)
NPP_statistical <- subset(NPP_statistical,TNPP_1>0)

#not fitting model basing on stepwise, just directly include all
d1 <- na.omit(NPP_statistical[,c("lon","lat","CABLE_NPP","Tg_a","PPFD_a","vpd_a","site_a")])
d1<- aggregate(d1,by=list(d1$lon,d1$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t1 <- stepwise_lm(d1,"CABLE_NPP")
t1[[2]]
#mod1 <- (lm(CABLE_NPP~fAPAR_a+Tg_a+vpd_a,data=d1))
mod1 <- (lm(CABLE_NPP~Tg_a+PPFD_a+vpd_a,data=d1))
summary(mod1)

d2 <- na.omit(NPP_statistical[,c("lon","lat","ISAM_npp","Tg_a","PPFD_a","vpd_a","site_a")])
d2<- aggregate(d2,by=list(d2$lon,d2$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t2 <- stepwise_lm(d2,"ISAM_npp")
t2[[1]]
t2[[2]]
mod2 <- (lm(ISAM_npp~Tg_a+PPFD_a+vpd_a,data=d2))
summary(mod2)

d3 <- na.omit(NPP_statistical[,c("lon","lat","ISBA_NPP","Tg_a","PPFD_a","vpd_a","site_a")])
d3<- aggregate(d3,by=list(d3$lon,d3$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t3 <- stepwise_lm(d3,"ISBA_NPP")
t3[[2]]
mod3 <- (lm(ISBA_NPP~Tg_a+PPFD_a+vpd_a,data=d3))
summary(mod3)

d4 <- na.omit(NPP_statistical[,c("lon","lat","JULES_NPP","Tg_a","PPFD_a","vpd_a","site_a")])
d4<- aggregate(d4,by=list(d4$lon,d4$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t4 <- stepwise_lm(d4,"JULES_NPP")
t4[[2]]
mod4 <- (lm(JULES_NPP~Tg_a+PPFD_a+vpd_a,data=d4))
summary(mod4)

d5 <- na.omit(NPP_statistical[,c("lon","lat","LPJ_NPP","Tg_a","PPFD_a","vpd_a","site_a")])
d5<- aggregate(d5,by=list(d5$lon,d5$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t5 <- stepwise_lm(d5,"LPJ_NPP")
t5[[2]]
mod5<- (lm(LPJ_NPP~Tg_a+PPFD_a+vpd_a,data=d5))
summary(mod5)

d6 <- na.omit(NPP_statistical[,c("lon","lat","ORCHIDEE_NPP","Tg_a","PPFD_a","vpd_a","site_a")])
d6<- aggregate(d6,by=list(d6$lon,d6$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t6 <- stepwise_lm(d6,"ORCHIDEE_NPP")
t6[[2]]
mod6<- (lm(ORCHIDEE_NPP~Tg_a+PPFD_a+vpd_a,data=d6))
summary(mod6)

d7 <- na.omit(NPP_statistical[,c("lon","lat","ORCHICNP_NPP","Tg_a","PPFD_a","vpd_a","site_a")])
d7<- aggregate(d7,by=list(d7$lon,d7$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t7 <- stepwise_lm(d7,"ORCHICNP_NPP")
t7[[2]]
mod7 <- (lm(ORCHICNP_NPP~Tg_a+PPFD_a+vpd_a,data=d7))
summary(mod7)

d8 <- na.omit(NPP_statistical[,c("lon","lat","SDGVM_NPP","Tg_a","PPFD_a","vpd_a","site_a")])
d8<- aggregate(d8,by=list(d8$lon,d8$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t8 <- stepwise_lm(d8,"SDGVM_NPP")
t8[[2]]
mod8<- (lm(SDGVM_NPP~Tg_a+PPFD_a+vpd_a,data=d8))
summary(mod8)

d9 <- na.omit(NPP_statistical[,c("lon","lat","CLASS_NPP","Tg_a","PPFD_a","vpd_a","site_a")])
d9<- aggregate(d9,by=list(d9$lon,d9$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t9 <- stepwise_lm(d9,"CLASS_NPP")
t9[[2]]
mod9<- (lm(CLASS_NPP~Tg_a+PPFD_a+vpd_a,data=d9))
summary(mod9)

d10 <- na.omit(NPP_statistical[,c("lon","lat","CLM_NPP","Tg_a","PPFD_a","vpd_a","site_a")])
d10<- aggregate(d10,by=list(d10$lon,d10$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t10 <- stepwise_lm(d10,"CLM_NPP")
t10[[2]]
mod10<- (lm(CLM_NPP~Tg_a+PPFD_a+vpd_a,data=d10))
summary(mod10)

d11 <- na.omit(NPP_statistical[,c("lon","lat","JSBACH_NPP","Tg_a","PPFD_a","vpd_a","site_a")])
d11<- aggregate(d11,by=list(d11$lon,d11$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t11 <- stepwise_lm(d11,"JSBACH_NPP")
t11[[2]]
mod11<- (lm(JSBACH_NPP~Tg_a+PPFD_a+vpd_a,data=d11))
summary(mod11)

d12 <- na.omit(NPP_statistical[,c("lon","lat","LPX_NPP","Tg_a","PPFD_a","vpd_a","site_a")])
d12<- aggregate(d12,by=list(d12$lon,d12$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t12 <- stepwise_lm(d12,"LPX_NPP")
t12[[2]]
mod12<- (lm(LPX_NPP~Tg_a+PPFD_a+vpd_a,data=d12))
summary(mod12)

mod1$coefficients;mod2$coefficients;mod3$coefficients;mod4$coefficients;
mod5$coefficients;mod6$coefficients;mod7$coefficients;mod8$coefficients;
mod9$coefficients;mod10$coefficients;mod11$coefficients;mod12$coefficients

t1a <- visreg(mod1,"Tg_a",type="contrast");t2a <-visreg(mod2,"Tg_a",type="contrast");t3a <-visreg(mod3,"Tg_a",type="contrast");t4a <- visreg(mod4,"Tg_a",type="contrast");t5a <- visreg(mod5,"Tg_a",type="contrast");t6a <-visreg(mod6,"Tg_a",type="contrast");t7a <- visreg(mod7,"Tg_a",type="contrast");t8a <-visreg(mod8,"Tg_a",type="contrast");t9a <-visreg(mod9,"Tg_a",type="contrast");t10a <-visreg(mod10,"Tg_a",type="contrast");t11a <-visreg(mod11,"Tg_a",type="contrast");t12a <-visreg(mod12,"Tg_a",type="contrast")
p1a <- visreg(mod1,"PPFD_a",type="contrast");p2a <- visreg(mod2,"PPFD_a",type="contrast");p3a <- visreg(mod3,"PPFD_a",type="contrast");p4a <- visreg(mod4,"PPFD_a",type="contrast");p5a <- visreg(mod5,"PPFD_a",type="contrast");p6a <- visreg(mod6,"PPFD_a",type="contrast");p7a <- visreg(mod7,"PPFD_a",type="contrast");p8a <- visreg(mod8,"PPFD_a",type="contrast");p9a <- visreg(mod9,"PPFD_a",type="contrast");p10a <- visreg(mod10,"PPFD_a",type="contrast");p11a <- visreg(mod11,"PPFD_a",type="contrast");p12a <- visreg(mod12,"PPFD_a",type="contrast")
v1a <- visreg(mod1,"vpd_a",type="contrast");v2a <- visreg(mod2,"vpd_a",type="contrast");v3a <- visreg(mod3,"vpd_a",type="contrast");v4a <- visreg(mod4,"vpd_a",type="contrast");v5a <- visreg(mod5,"vpd_a",type="contrast");v6a <- visreg(mod6,"vpd_a",type="contrast");v7a <- visreg(mod7,"vpd_a",type="contrast");v8a <- visreg(mod8,"vpd_a",type="contrast");v9a <- visreg(mod9,"vpd_a",type="contrast");v10a <- visreg(mod10,"vpd_a",type="contrast");v11a <- visreg(mod11,"vpd_a",type="contrast");v12a <- visreg(mod12,"vpd_a",type="contrast")

fits_tg <- dplyr::bind_rows(mutate(bp1a$fit, plt = "Measurement"),mutate(t1a$fit, plt = "CABLE"),mutate(t2a$fit, plt = "ISAM"),mutate(t3a$fit, plt = "ISBA"),mutate(t4a$fit, plt = "JULES"),mutate(t5a$fit, plt = "LPJ"),mutate(t6a$fit, plt = "ORCHIDEE"),mutate(t7a$fit, plt = "ORCHICNP"),mutate(t8a$fit, plt = "SDGVM"),mutate(t9a$fit, plt = "CLASS"),mutate(t10a$fit, plt = "CLM"),mutate(t11a$fit, plt = "JSBACH"),mutate(t12a$fit, plt = "LPX"))

fits_PPFD <- dplyr::bind_rows(mutate(bp1b$fit, plt = "Measurement"),mutate(p1a$fit, plt = "CABLE"),mutate(p2a$fit, plt = "ISAM"),mutate(p3a$fit, plt = "ISBA"),mutate(p4a$fit, plt = "JULES"),mutate(p5a$fit, plt = "LPJ"),mutate(p6a$fit, plt = "ORCHIDEE"),mutate(p7a$fit, plt = "ORCHICNP"),mutate(p8a$fit, plt = "SDGVM"),mutate(p9a$fit, plt = "CLASS"),mutate(p10a$fit, plt = "CLM"),mutate(p11a$fit, plt = "JSBACH"),mutate(p12a$fit, plt = "LPX"))

fits_vpd <- dplyr::bind_rows(mutate(bp1c$fit, plt = "Measurement"),mutate(v1a$fit, plt = "CABLE"),mutate(v2a$fit, plt = "ISAM"),mutate(v3a$fit, plt = "ISBA"),mutate(v4a$fit, plt = "JULES"),mutate(v5a$fit, plt = "LPJ"),mutate(v6a$fit, plt = "ORCHIDEE"),mutate(v7a$fit, plt = "ORCHICNP"),mutate(v8a$fit, plt = "SDGVM"),mutate(v9a$fit, plt = "CLASS"),mutate(v10a$fit, plt = "CLM"),mutate(v11a$fit, plt = "JSBACH"),mutate(v12a$fit, plt = "LPX"))

final1 <- ggplot() +xlim(8,30)+geom_line(data = subset(fits_tg,plt=="Measurement"), aes(Tg_a, visregFit, group=plt, color=plt),size=4)+geom_line(data = fits_tg, aes(Tg_a, visregFit, group=plt, color=plt),size=1) +theme_classic()+theme(text = element_text(size=20),legend.position="none")+  geom_ribbon(data = bp1a$fit,aes(Tg_a, ymin=visregLwr, ymax=visregUpr),fill="gray",alpha=0.5)+scale_colour_manual(values=c(Measurement="black",CABLE="pink",ISAM="red",ISBA="#339999",JULES="#663399",LPJ="#0066CC",ORCHIDEE="#FF9933",ORCHICNP="cyan",SDGVM="yellow",CLASS="brown",CLM="darkgoldenrod",JSBACH="burlywood1",LPX="darkgreen"))+labs(y = ~paste(BP, " (gC m"^-2,"yr"^-1,")"))+xlab("Tg (°C)")
final2 <- ggplot()+xlim(5.5,6.3)+geom_line(data = subset(fits_PPFD,plt=="Measurement"), aes(PPFD_a, visregFit, group=plt, color=plt),size=4) +geom_line(data = fits_PPFD, aes(PPFD_a, visregFit, group=plt, color=plt),size=1) + ylab(" ")+theme_classic()+theme(text = element_text(size=20),legend.position="none")+geom_ribbon(data = bp1b$fit,aes(PPFD_a, ymin=visregLwr, ymax=visregUpr),fill="gray",alpha=0.5)+ scale_colour_manual(values=c(Measurement="black",CABLE="pink",ISAM="red",ISBA="#339999",JULES="#663399",LPJ="#0066CC",ORCHIDEE="#FF9933",ORCHICNP="cyan",SDGVM="yellow",CLASS="brown",CLM="darkgoldenrod",JSBACH="burlywood1",LPX="darkgreen"))+xlab(~paste("ln PPFD (", mu, "mol m"^-2,"s"^-1, ")"))
final3 <- ggplot()+xlim(-0.8,1) +geom_line(data = subset(fits_vpd,plt=="Measurement"), aes(vpd_a, visregFit, group=plt, color=plt),size=4)+geom_line(data = fits_vpd, aes(vpd_a, visregFit, group=plt, color=plt),size=1) + ylab(" ") +theme_classic()+theme(text = element_text(size=20),legend.position="none")+ geom_ribbon(data = bp1c$fit,aes(vpd_a, ymin=visregLwr, ymax=visregUpr),fill="gray",alpha=0.5)+scale_colour_manual(values=c(Measurement="black",CABLE="pink",ISAM="red",ISBA="#339999",JULES="#663399",LPJ="#0066CC",ORCHIDEE="#FF9933",ORCHICNP="cyan",SDGVM="yellow",CLASS="brown",CLM="darkgoldenrod",JSBACH="burlywood1",LPX="darkgreen"))+xlab("ln D (kPa)")

#show legend
final1_legend <- ggplot() +geom_line(data = fits_PPFD, aes(PPFD_a, visregFit, group=plt, color=plt),size=2) + xlab("ln PPFD") + ylab(" ")+theme_classic()+theme(text = element_text(size=20))+
  scale_colour_manual(" ",values=c(Measurement="black",CABLE="pink",ISAM="red",ISBA="#339999",JULES="#663399",LPJ="#0066CC",ORCHIDEE="#FF9933",ORCHICNP="cyan",SDGVM="yellow",CLASS="brown",CLM="darkgoldenrod",JSBACH="burlywood1",LPX="darkgreen"), 
                      labels = c("Our model","CABLE","ISAM","ISBA","JULES","LPJ","ORCHIDEE","ORCHICNP","SDGVM","CLASS","CLM","JSBACH","LPX"))

legend_info <- as_ggplot(get_legend(final1_legend))

#nuptake figure
#N minerlization
Nmin_statistical <- subset(NPP_forest,Nmin>0)

sitemean2 <- unique(Nmin_statistical[,c("lon","lat")])
sp_sites2 <- SpatialPoints(sitemean2) # only select lon and lat

sitemean2_final <- sitemean2

allmaps2 <- list(ORCHICNP_fNup,ISAM_fNup)
obs_name2 <- c("ORCHICNP_fNup","ISAM_fNup")

for (i in c(1:length(allmaps2))){
  df1 <- allmaps2[[i]]
  names(df1) <- c("lon","lat",obs_name2[i])
  coordinates(df1) <- ~lon+lat 
  gridded(df1) <- TRUE
  df1_global <- raster(df1, obs_name2[i]) 
  
  sp_sites_new <- raster::extract(df1_global, sp_sites2, sp = TRUE) %>% as_tibble() %>% 
    right_join(sitemean2, by = c("lon", "lat"))
  sitemean2_final[,i+2] <- sp_sites_new[,1]
}

JSBACH_fNup <- (raster::extract(raster_JSBACH_fNup, sp_sites2, sp = TRUE) %>% as_tibble() %>% right_join(sitemean2, by = c("lon", "lat")))[,1]
LPX_fNup <- (raster::extract(raster_LPX_fNup, sp_sites2, sp = TRUE) %>% as_tibble() %>% right_join(sitemean2, by = c("lon", "lat")))[,1]
new_added2 <- as.data.frame(cbind(JSBACH_fNup,LPX_fNup))
names(new_added2) <- c("JSBACH_fNup","LPX_fNup")

sitemean2_final <- as.data.frame(cbind(sitemean2_final,new_added2))

sitemean2_final$ORCHICNP_fNup[sitemean2_final$ORCHICNP_fNup==0] <- NA
sitemean2_final$ISAM_fNup[sitemean2_final$ISAM_fNup==0] <- NA

Nmin_statistical <- merge(Nmin_statistical,sitemean2_final,by=c("lon","lat"),all.x=TRUE)
Nmin_statistical$Nmin_a <- Nmin_statistical$Nmin
Nmin_statistical_final <- na.omit(Nmin_statistical[,c("pred_nuptake","Tg_a","PPFD_a","vpd_a","site_a","lon","lat")])
Nmin_statistical_final<- aggregate(Nmin_statistical_final,by=list(Nmin_statistical_final$lon,Nmin_statistical_final$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))

#TRENDY1
m1 <- na.omit(Nmin_statistical[,c("lon","lat","ORCHICNP_fNup","Tg_a","PPFD_a","vpd_a","site_a")])
m1<- aggregate(m1,by=list(m1$lon,m1$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))

m2 <- na.omit(Nmin_statistical[,c("lon","lat","ISAM_fNup","Tg_a","PPFD_a","vpd_a","site_a")])
m2<- aggregate(m2,by=list(m2$lon,m2$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))

m3 <- na.omit(Nmin_statistical[,c("lon","lat","JSBACH_fNup","Tg_a","PPFD_a","vpd_a","site_a")])
m3<- aggregate(m3,by=list(m3$lon,m3$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))

m4 <- na.omit(Nmin_statistical[,c("lon","lat","LPX_fNup","Tg_a","PPFD_a","vpd_a","site_a")])
m4<- aggregate(m4,by=list(m4$lon,m4$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))


#final figure for Nup
mod_n1 <- lm(pred_nuptake~Tg_a+PPFD_a+vpd_a,data=Nmin_statistical_final)
r.squaredGLMM(mod_n1)
summary(mod_n1)
nn1b <- visreg(mod_n1,"Tg_a",type="contrast")
nn1c <- visreg(mod_n1,"PPFD_a",type="contrast")
nn1d <- visreg(mod_n1,"vpd_a",type="contrast")

mod_n2 <- lm(ORCHICNP_fNup~Tg_a+PPFD_a+vpd_a,data=m1)
nn2b <- visreg(mod_n2,"Tg_a",type="contrast")
nn2c <- visreg(mod_n2,"PPFD_a",type="contrast")
nn2d <- visreg(mod_n2,"vpd_a",type="contrast")

mod_n3 <- lm(ISAM_fNup~Tg_a+PPFD_a+vpd_a,data=m2)
nn3b <- visreg(mod_n3,"Tg_a",type="contrast")
nn3c <- visreg(mod_n3,"PPFD_a",type="contrast")
nn3d <- visreg(mod_n3,"vpd_a",type="contrast")

mod_n4 <- lm(JSBACH_fNup~Tg_a+PPFD_a+vpd_a,data=m3)
nn4b <- visreg(mod_n4,"Tg_a",type="contrast")
nn4c <- visreg(mod_n4,"PPFD_a",type="contrast")
nn4d <- visreg(mod_n4,"vpd_a",type="contrast")

mod_n5 <- lm(LPX_fNup~Tg_a+PPFD_a+vpd_a,data=m4)
nn5b <- visreg(mod_n5,"Tg_a",type="contrast")
nn5c <- visreg(mod_n5,"PPFD_a",type="contrast")
nn5d <- visreg(mod_n5,"vpd_a",type="contrast")

fits_tg <- dplyr::bind_rows(mutate(nn1b$fit, plt = "Measurement"),
                            mutate(nn2b$fit, plt = "ORCHICNP"),
                            mutate(nn3b$fit, plt = "ISAM"),
                            mutate(nn4b$fit, plt = "JSBACH"),
                            mutate(nn5b$fit, plt = "LPX"))

fits_ppfd <- dplyr::bind_rows(mutate(nn1c$fit, plt = "Measurement"),
                              mutate(nn2c$fit, plt = "ORCHICNP"),
                              mutate(nn3c$fit, plt = "ISAM"),
                              mutate(nn4c$fit, plt = "JSBACH"),
                              mutate(nn5c$fit, plt = "LPX"))

fits_vpd1 <- dplyr::bind_rows(mutate(nn1d$fit, plt = "Measurement"),
                              mutate(nn2d$fit, plt = "ORCHICNP"),
                              mutate(nn3d$fit, plt = "ISAM"),
                              mutate(nn4d$fit, plt = "JSBACH"),
                              mutate(nn5d$fit, plt = "LPX"))

final1b <- ggplot()+xlim(8,30) + geom_line(data = subset(fits_tg,plt=="Measurement"),aes(Tg_a, visregFit, group=plt, color=plt),size=4) +geom_line(data = fits_tg, aes(Tg_a, visregFit, group=plt, color=plt),size=1) + 
  geom_ribbon(data=nn1b$fit,aes(Tg_a, ymin=visregLwr, ymax=visregUpr),fill="gray",alpha=0.5)+
  xlab("Tg (°C)") +labs(y = ~paste("N uptake", " (gN m"^-2,"yr"^-1,")")) +
  theme_classic()+theme(text = element_text(size=20),legend.position="none")+
  scale_colour_manual(values=c(Measurement="black",ISAM="red",ORCHICNP="cyan",JSBACH="burlywood1",LPX="darkgreen"))

final1c <- ggplot()+xlim(5.5,6.3) +geom_line(data = subset(fits_ppfd,plt=="Measurement"),aes(PPFD_a, visregFit, group=plt, color=plt),size=4) +geom_line(data = fits_ppfd, aes(PPFD_a, visregFit, group=plt, color=plt),size=1) + 
  geom_ribbon(data=nn1c$fit,aes(PPFD_a, ymin=visregLwr, ymax=visregUpr),fill="gray",alpha=0.5)+
  xlab(~paste("ln PPFD (", mu, "mol m"^-2,"s"^-1, ")")) + ylab(" ")+
  theme_classic()+theme(text = element_text(size=20),legend.position="none")+
  scale_colour_manual(values=c(Measurement="black",ISAM="red",ORCHICNP="cyan",JSBACH="burlywood1",LPX="darkgreen"))

final1d <- ggplot()+xlim(-0.8,1) +geom_line(data = subset(fits_vpd1,plt=="Measurement"),aes(vpd_a, visregFit, group=plt, color=plt),size=4) +geom_line(data = fits_vpd1, aes(vpd_a, visregFit, group=plt, color=plt),size=1) + 
  geom_ribbon(data=nn1d$fit,aes(vpd_a, ymin=visregLwr, ymax=visregUpr),fill="gray",alpha=0.5)+
  xlab("ln D (kPa)") + ylab(" ")+
  theme_classic()+theme(text = element_text(size=20),legend.position="none")+
  scale_colour_manual(values=c(Measurement="black",ISAM="red",ORCHICNP="cyan",JSBACH="burlywood1",LPX="darkgreen"))

#nue model for TRENDY
#nuptake figure
#N minerlization
NUE_test <- NPP_all

sitemean3 <- unique(NUE_test[,c("lon","lat")])
sp_sites3 <- SpatialPoints(sitemean3) # only select lon and lat

sitemean3_final <- sitemean3

allmaps3 <- list(ORCHICNP_fNup,ISAM_fNup,ORCHICNP_NPP,ISAM_npp)
obs_name3 <- c("ORCHICNP_fNup","ISAM_fNup","ORCHICNP_NPP","ISAM_npp")

for (i in c(1:length(allmaps3))){
  df1 <- allmaps3[[i]]
  names(df1) <- c("lon","lat",obs_name3[i])
  coordinates(df1) <- ~lon+lat 
  gridded(df1) <- TRUE
  df1_global <- raster(df1, obs_name3[i]) 
  
  sp_sites_new <- raster::extract(df1_global, sp_sites3, sp = TRUE) %>% as_tibble() %>% 
    right_join(sitemean3, by = c("lon", "lat"))
  sitemean3_final[,i+2] <- sp_sites_new[,1]
}

JSBACH_NPP <- (raster::extract(raster_JSBACH_npp, sp_sites3, sp = TRUE) %>% as_tibble() %>% right_join(sitemean3, by = c("lon", "lat")))[,1]
JSBACH_fNup <- (raster::extract(raster_JSBACH_fNup, sp_sites3, sp = TRUE) %>% as_tibble() %>% right_join(sitemean3, by = c("lon", "lat")))[,1]
LPX_NPP <- (raster::extract(raster_LPX_npp, sp_sites3, sp = TRUE) %>% as_tibble() %>% right_join(sitemean3, by = c("lon", "lat")))[,1]
LPX_fNup <- (raster::extract(raster_LPX_fNup, sp_sites3, sp = TRUE) %>% as_tibble() %>% right_join(sitemean3, by = c("lon", "lat")))[,1]
new_added3 <- as.data.frame(cbind(JSBACH_NPP,JSBACH_fNup,LPX_NPP,LPX_fNup))
names(new_added3) <- c("JSBACH_NPP","JSBACH_fNup","LPX_NPP","LPX_fNup")

sitemean3_final <- as.data.frame(cbind(sitemean3_final,new_added3))

sitemean3_final$ORCHICNP_fNup[sitemean3_final$ORCHICNP_fNup==0] <- NA
sitemean3_final$ISAM_npp[sitemean3_final$ISAM_npp==0] <- NA
sitemean3_final$LPX_NPP[sitemean3_final$LPX_NPP==0] <- NA
sitemean3_final$LPX_fNup[sitemean3_final$LPX_fNup==0] <- NA

NUE_test <- merge(NUE_test,sitemean3_final,by=c("lon","lat"),all.x=TRUE)
NUE_test$NUE_ORCHICNP <- NUE_test$ORCHICNP_NPP/NUE_test$ORCHICNP_fNup
NUE_test$NUE_ISAM <- NUE_test$ISAM_npp/NUE_test$ISAM_fNup
NUE_test$NUE_JSBACH <- NUE_test$JSBACH_NPP/NUE_test$JSBACH_fNup
NUE_test$NUE_LPX <- NUE_test$LPX_NPP/NUE_test$LPX_fNup

NUE_test_sitemean <- aggregate(NUE_test,by=list(NUE_test$lon,NUE_test$lat), FUN=mean, na.rm=TRUE)

NUE_test_final <- na.omit(NUE_test_sitemean[,c("NUE_ORCHICNP","Tg_a","PPFD_a","vpd_a")])

NUE_test_final2 <- na.omit(NUE_test_sitemean[,c("NUE_ISAM","Tg_a","PPFD_a","vpd_a")])

NUE_test_final3 <- na.omit(NUE_test_sitemean[,c("NUE_JSBACH","Tg_a","PPFD_a","vpd_a")])

NUE_test_final4 <- na.omit(NUE_test_sitemean[,c("NUE_LPX","Tg_a","PPFD_a","vpd_a")])

#NUE_ORCHICNP
mod_nue1 <- (lm(NUE_ORCHICNP~vpd_a+Tg_a+PPFD_a,data=NUE_test_final))
nue4b <- visreg(mod_nue1,"vpd_a",type="contrast")
nue6b <- visreg(mod_nue1,"Tg_a",type="contrast")
nue8b <- visreg(mod_nue1,"PPFD_a",type="contrast")

#NUE_ISAM
mod_nue2 <- (lm(NUE_ISAM~vpd_a+Tg_a+PPFD_a,data=NUE_test_final2))
nue4a <- visreg(mod_nue2,"vpd_a",type="contrast")
nue6a <- visreg(mod_nue2,"Tg_a",type="contrast")
nue8a <- visreg(mod_nue2,"PPFD_a",type="contrast")

#NUE_JSBACH
mod_nue3 <- (lm(NUE_JSBACH~vpd_a+Tg_a+PPFD_a,data=NUE_test_final3))
nue4d <- visreg(mod_nue3,"vpd_a",type="contrast")
nue6d <- visreg(mod_nue3,"Tg_a",type="contrast")
nue8d <- visreg(mod_nue3,"PPFD_a",type="contrast")

#NUE_LPX
mod_nue4 <- (lm(NUE_LPX~vpd_a+Tg_a+PPFD_a,data=NUE_test_final4))
nue4e <- visreg(mod_nue4,"vpd_a",type="contrast")
nue6e <- visreg(mod_nue4,"Tg_a",type="contrast")
nue8e <- visreg(mod_nue4,"PPFD_a",type="contrast")

#also, our simulations
NPP_grassland_remove <- NPP_grassland
forest_NUE <- NPP_forest[,c("lon","lat","pred_nuptake","pred_npp","Tg_a","PPFD_a","vpd_a")]
grassland_NUE <- NPP_grassland_remove[,c("lon","lat","grassland_pred_nuptake","grassland_pred_npp","Tg_a","PPFD_a","vpd_a")]
names(grassland_NUE) <- c("lon","lat","pred_nuptake","pred_npp","Tg_a","PPFD_a","vpd_a")
nue_allplots <- as.data.frame(rbind(forest_NUE,grassland_NUE))
nue_allplots$nue <-nue_allplots$pred_npp/nue_allplots$pred_nuptake
nue_allplots <- na.omit(nue_allplots[,c("nue","Tg_a","vpd_a","PPFD_a","lon","lat")])
nue_allplots <- aggregate(nue_allplots,by=list(nue_allplots$lon,nue_allplots$lat), FUN=mean, na.rm=TRUE)

mod_nue3 <- (lm(nue~vpd_a+Tg_a+PPFD_a,data=nue_allplots))
summary(mod_nue3)
nue4c <- visreg(mod_nue3,"vpd_a",type="contrast")
nue6c <- visreg(mod_nue3,"Tg_a",type="contrast")
nue8c <- visreg(mod_nue3,"PPFD_a",type="contrast")

fits_vpd <- dplyr::bind_rows(mutate(nue4c$fit, plt = "Measurement"),mutate(nue4b$fit, plt = "ORCHICNP"),mutate(nue4a$fit, plt = "ISAM"),mutate(nue4d$fit, plt = "JSBACH"),mutate(nue4e$fit, plt = "LPX"))
fits_Tg <- dplyr::bind_rows(mutate(nue6c$fit, plt = "Measurement"),mutate(nue6b$fit, plt = "ORCHICNP"),mutate(nue6a$fit, plt = "ISAM"),mutate(nue6d$fit, plt = "JSBACH"),mutate(nue6e$fit, plt = "LPX"))
fits_PPFD <- dplyr::bind_rows(mutate(nue8c$fit, plt = "Measurement"),mutate(nue8b$fit, plt = "ORCHICNP"),mutate(nue8a$fit, plt = "ISAM"),mutate(nue8d$fit, plt = "JSBACH"),mutate(nue8e$fit, plt = "LPX"))


final1_nue4 <- ggplot() +xlim(-0.8,1)+geom_line(data = nue4c$fit, aes(vpd_a, visregFit),size=4) +
  geom_line(data = fits_vpd, aes(vpd_a, visregFit, group=plt, color=plt),size=1) +
  xlab("ln D (kPa)") + ylab(" ")+
  geom_ribbon(data=nue4c$fit,aes(vpd_a, ymin=visregLwr, ymax=visregUpr),fill="gray",alpha=0.5)+
  theme_classic()+theme(text = element_text(size=20),legend.position="none")+scale_colour_manual(values=c(Measurement="black",ISAM="red",ORCHICNP="cyan",JSBACH="burlywood1",LPX="darkgreen"))

final1_nue6 <- ggplot() +xlim(8,30)+geom_line(data = nue6c$fit, aes(Tg_a, visregFit),size=4) +
  geom_line(data = fits_Tg, aes(Tg_a, visregFit, group=plt, color=plt),size=1) +
  xlab("Tg (°C)") + ylab("NUE")+
  geom_ribbon(data=nue6c$fit,aes(Tg_a, ymin=visregLwr, ymax=visregUpr),fill="gray",alpha=0.5)+
  theme_classic()+theme(text = element_text(size=20),legend.position="none")+scale_colour_manual(values=c(Measurement="black",ISAM="red",ORCHICNP="cyan",JSBACH="burlywood1",LPX="darkgreen"))
 
final1_nue8 <- ggplot()+xlim(5.5,6.3) +geom_line(data = nue8c$fit, aes(PPFD_a, visregFit),size=4) +
  geom_line(data = fits_PPFD, aes(PPFD_a, visregFit, group=plt, color=plt),size=1) +
  xlab(~paste("ln PPFD (", mu, "mol m"^-2,"s"^-1, ")")) + ylab(" ")+
  geom_ribbon(data=nue8c$fit,aes(PPFD_a, ymin=visregLwr, ymax=visregUpr),fill="gray",alpha=0.5)+
  theme_classic()+theme(text = element_text(size=20),legend.position="none")+scale_colour_manual(values=c(Measurement="black",ISAM="red",ORCHICNP="cyan",JSBACH="burlywood1",LPX="darkgreen"))

plot_grid(final1,final2,final3,legend_info,
          final1b,final1c,final1d,white,
          final1_nue6,final1_nue8,final1_nue4,white,
          nrow=3,label_x = 0.8, label_y = 0.8)+white

ggsave(paste("~/data/output/newphy_fig5.jpg",sep=""), width = 20, height = 15)

#nue and nuptake relation
#also, our simulations
NPP_grassland_remove <- NPP_grassland
forest_NUE <- NPP_forest[,c("lon","lat","pred_nuptake","pred_npp","age_a","Tg_a","PPFD_a","vpd_a","fAPAR_a",
                            "CNrt_a","LMA_a","vcmax25_a")]
grassland_NUE <- NPP_grassland_remove[,c("lon","lat","grassland_pred_nuptake","grassland_pred_npp","age_a","Tg_a","PPFD_a","vpd_a","fAPAR_a",
                                         "CNrt_a","LMA_a","vcmax25_a")]
names(grassland_NUE) <- c("lon","lat","pred_nuptake","pred_npp","age_a","Tg_a","PPFD_a","vpd_a","fAPAR_a",
                          "CNrt_a","LMA_a","vcmax25_a")
nue_allplots <- as.data.frame(rbind(forest_NUE,grassland_NUE))
nue_allplots$nue <-nue_allplots$pred_npp/nue_allplots$pred_nuptake

mod_nue3 <- (lm(nue~LMA_a+fAPAR_a+age_a+vpd_a+CNrt_a+Tg_a+vcmax25_a+PPFD_a,data=nue_allplots))
summary(mod_nue3)
nue1c <- visreg(mod_nue3,"LMA_a",type="contrast")
nue2c <- visreg(mod_nue3,"fAPAR_a",type="contrast")
nue3c <- visreg(mod_nue3,"age_a",type="contrast")
nue4c <- visreg(mod_nue3,"vpd_a",type="contrast")
nue5c <- visreg(mod_nue3,"CNrt_a",type="contrast")
nue6c <- visreg(mod_nue3,"Tg_a",type="contrast")
nue7c <- visreg(mod_nue3,"vcmax25_a",type="contrast")
nue8c <- visreg(mod_nue3,"PPFD_a",type="contrast")

final1_nue1 <- ggplot() +geom_line(data = nue1c$fit, aes(LMA_a, visregFit),size=2) +ylab("NUE") +  xlab("ln LMA") +theme_classic()+theme(text = element_text(size=20),legend.position="none")
final1_nue2 <- ggplot() +geom_line(data = nue2c$fit, aes(fAPAR_a, visregFit),size=2) + xlab("fAPAR") + ylab(" ")+theme_classic()+theme(text = element_text(size=20),legend.position="none")
final1_nue3 <- ggplot() +geom_line(data = nue3c$fit, aes(age_a, visregFit),size=2) + xlab("ln age") + ylab(" ")+theme_classic()+theme(text = element_text(size=20),legend.position="none")
final1_nue4 <- ggplot() +geom_line(data = nue4c$fit, aes(vpd_a, visregFit),size=2) + xlab("ln D") + ylab(" ")+theme_classic()+theme(text = element_text(size=20),legend.position="none")
final1_nue5 <- ggplot() +geom_line(data = nue5c$fit, aes(CNrt_a, visregFit),size=2)+ylab("NUE")  + xlab("ln soil C/N") +theme_classic()+theme(text = element_text(size=20),legend.position="none")
final1_nue6 <- ggplot() +geom_line(data = nue6c$fit, aes(Tg_a, visregFit),size=2) + xlab("Tg") + ylab(" ")+theme_classic()+theme(text = element_text(size=20),legend.position="none")
final1_nue7 <- ggplot() +geom_line(data = nue7c$fit, aes(vcmax25_a, visregFit),size=2) + xlab("ln vcmax25") + ylab(" ")+theme_classic()+theme(text = element_text(size=20),legend.position="none")
final1_nue8 <- ggplot() +geom_line(data = nue8c$fit, aes(PPFD_a, visregFit),size=2) + xlab("ln PPFD") + ylab(" ")+theme_classic()+theme(text = element_text(size=20),legend.position="none")

plot_grid(final1_nue1,final1_nue2,final1_nue3,final1_nue4,
          final1_nue5,final1_nue6,final1_nue7,final1_nue8,
          nrow=2,label_x = 0.8, label_y = 0.8)+white

ggsave(paste("~/data/output/newphy_add1.jpg",sep=""), width = 20, height = 10)

# n uptake
mod_nue3 <- (lm(pred_nuptake~LMA_a+fAPAR_a+age_a+vpd_a+CNrt_a+Tg_a+vcmax25_a+PPFD_a,data=nue_allplots))
summary(mod_nue3)
nue1c <- visreg(mod_nue3,"LMA_a",type="contrast")
nue2c <- visreg(mod_nue3,"fAPAR_a",type="contrast")
nue3c <- visreg(mod_nue3,"age_a",type="contrast")
nue4c <- visreg(mod_nue3,"vpd_a",type="contrast")
nue5c <- visreg(mod_nue3,"CNrt_a",type="contrast")
nue6c <- visreg(mod_nue3,"Tg_a",type="contrast")
nue7c <- visreg(mod_nue3,"vcmax25_a",type="contrast")
nue8c <- visreg(mod_nue3,"PPFD_a",type="contrast")

final1_nue1 <- ggplot() +geom_line(data = nue1c$fit, aes(LMA_a, visregFit),size=2) +ylab("N uptake") +  xlab("ln LMA") +theme_classic()+theme(text = element_text(size=20),legend.position="none")
final1_nue2 <- ggplot() +geom_line(data = nue2c$fit, aes(fAPAR_a, visregFit),size=2) + xlab("fAPAR") + ylab(" ")+theme_classic()+theme(text = element_text(size=20),legend.position="none")
final1_nue3 <- ggplot() +geom_line(data = nue3c$fit, aes(age_a, visregFit),size=2) + xlab("ln age") + ylab(" ")+theme_classic()+theme(text = element_text(size=20),legend.position="none")
final1_nue4 <- ggplot() +geom_line(data = nue4c$fit, aes(vpd_a, visregFit),size=2) + xlab("ln D") + ylab(" ")+theme_classic()+theme(text = element_text(size=20),legend.position="none")
final1_nue5 <- ggplot() +geom_line(data = nue5c$fit, aes(CNrt_a, visregFit),size=2)+ylab("N uptake")  + xlab("ln soil C/N") +theme_classic()+theme(text = element_text(size=20),legend.position="none")
final1_nue6 <- ggplot() +geom_line(data = nue6c$fit, aes(Tg_a, visregFit),size=2) + xlab("Tg") + ylab(" ")+theme_classic()+theme(text = element_text(size=20),legend.position="none")
final1_nue7 <- ggplot() +geom_line(data = nue7c$fit, aes(vcmax25_a, visregFit),size=2) + xlab("ln vcmax25") + ylab(" ")+theme_classic()+theme(text = element_text(size=20),legend.position="none")
final1_nue8 <- ggplot() +geom_line(data = nue8c$fit, aes(PPFD_a, visregFit),size=2) + xlab("ln PPFD") + ylab(" ")+theme_classic()+theme(text = element_text(size=20),legend.position="none")

plot_grid(final1_nue1,final1_nue2,final1_nue3,final1_nue4,
          final1_nue5,final1_nue6,final1_nue7,final1_nue8,
          nrow=2,label_x = 0.8, label_y = 0.8)+white

ggsave(paste("~/data/output/newphy_add2.jpg",sep=""), width = 20, height = 10)

#validation - BP
NPP_statistical$Measured_BP <- NPP_statistical$TNPP_1

pp2 <- analyse_modobs2(NPP_statistical,"pred_npp","Measured_BP", type = "points",relative=TRUE)$gg + larger_size+
  labs(y = ~paste("Forest ", BP[obs.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("Forest ", BP[pred.], " (gC m"^-2,"yr"^-1,")"))

pp4 <- analyse_modobs2(NPP_statistical,"CABLE_NPP","Measured_BP", type = "points",relative=TRUE)$gg + larger_size+
  labs(y = ~paste("Forest ", BP[obs.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("CABLE ", BP[pred.], " (gC m"^-2,"yr"^-1,")"))

pp5 <- analyse_modobs2(NPP_statistical,"ISAM_npp","Measured_BP", type = "points",relative=TRUE)$gg + larger_size+
  labs(y = ~paste("Forest ", BP[obs.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("ISAM ", BP[pred.], " (gC m"^-2,"yr"^-1,")"))

pp6 <- analyse_modobs2(NPP_statistical,"ISBA_NPP","Measured_BP", type = "points",relative=TRUE)$gg + larger_size+
  labs(y = ~paste("Forest ", BP[obs.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("ISBA ", BP[pred.], " (gC m"^-2,"yr"^-1,")"))

pp7 <- analyse_modobs2(NPP_statistical,"JULES_NPP","Measured_BP", type = "points",relative=TRUE)$gg + larger_size+
  labs(y = ~paste("Forest ", BP[obs.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("JULES ", BP[pred.], " (gC m"^-2,"yr"^-1,")"))

pp8 <- analyse_modobs2(NPP_statistical,"LPJ_NPP","Measured_BP", type = "points",relative=TRUE)$gg + larger_size+
  labs(y = ~paste("Forest ", BP[obs.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("LPJ ", BP[pred.], " (gC m"^-2,"yr"^-1,")"))

pp9 <- analyse_modobs2(NPP_statistical,"ORCHIDEE_NPP","Measured_BP", type = "points",relative=TRUE)$gg + larger_size+
  labs(y = ~paste("Forest ", BP[obs.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("ORCHIDEE ", BP[pred.], " (gC m"^-2,"yr"^-1,")"))

pp10 <- analyse_modobs2(NPP_statistical,"ORCHICNP_NPP","Measured_BP", type = "points",relative=TRUE)$gg + larger_size+
  labs(y = ~paste("Forest ", BP[obs.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("ORCHICNP ", BP[pred.], " (gC m"^-2,"yr"^-1,")"))

pp11 <- analyse_modobs2(NPP_statistical,"SDGVM_NPP","Measured_BP", type = "points",relative=TRUE)$gg + larger_size+
  labs(y = ~paste("Forest ", BP[obs.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("SDGVM ", BP[pred.], " (gC m"^-2,"yr"^-1,")"))

pp12 <- analyse_modobs2(NPP_statistical,"CLASS_NPP","Measured_BP", type = "points",relative=TRUE)$gg + larger_size+
  labs(y = ~paste("Forest ", BP[obs.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("CLASS ", BP[pred.], " (gC m"^-2,"yr"^-1,")"))

pp13 <- analyse_modobs2(NPP_statistical,"CLM_NPP","Measured_BP", type = "points",relative=TRUE)$gg + larger_size+
  labs(y = ~paste("Forest ", BP[obs.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("CLM ", BP[pred.], " (gC m"^-2,"yr"^-1,")"))

pp14 <- analyse_modobs2(NPP_statistical,"JSBACH_NPP","Measured_BP", type = "points",relative=TRUE)$gg + larger_size+
  labs(y = ~paste("Forest ", BP[obs.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("JSBACH ", BP[pred.], " (gC m"^-2,"yr"^-1,")"))

pp15 <- analyse_modobs2(NPP_statistical,"LPX_NPP","Measured_BP", type = "points",relative=TRUE)$gg + larger_size+
  labs(y = ~paste("Forest ", BP[obs.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("LPX ", BP[pred.], " (gC m"^-2,"yr"^-1,")"))

plot_grid(pp2,pp4,pp5,pp6,pp7,pp8,pp9,pp10,pp11,pp12,pp13,pp14,pp15,
          labels = c('(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)','(m)'),
          nrow=5,label_x = 0.9,label_y=0.92,label_size = 20)+white

ggsave(paste("~/data/output/newphy_figs5.jpg",sep=""), width = 23, height = 25)

#validation - N uptake
Nmin_statistical$pred_nuptake

Nmin_statistical$Nuptake_from_Nup_model <- summary(mod_n1)$coef[1,1]+
  summary(mod_n1)$coef[2,1]*Nmin_statistical$fAPAR_a+
  summary(mod_n1)$coef[3,1]*Nmin_statistical$Tg_a+
  summary(mod_n1)$coef[4,1]*Nmin_statistical$LMA_a

ppp1 <- analyse_modobs2(Nmin_statistical,"pred_nuptake","Nmin", type = "points",relative=TRUE)$gg +larger_size+
  labs(y = ~paste("Net ", minerlization[obs.], " (gN m"^-2,"yr"^-1,")")) +labs(x = ~paste("Forest N ", uptake[pred.], " (gN m"^-2,"yr"^-1,")"))

#remove direct nuptake model since it is irrelevant
#ppp2 <- analyse_modobs2(Nmin_statistical,"Nuptake_from_Nup_model","Nmin", type = "points")

ppp3 <- analyse_modobs2(Nmin_statistical,"ORCHICNP_fNup","Nmin", type = "points",relative=TRUE)$gg +larger_size+
  labs(y = ~paste("Net ", minerlization[obs.], " (gN m"^-2,"yr"^-1,")")) +labs(x = ~paste("ORCHICNP N ", uptake[pred.], " (gN m"^-2,"yr"^-1,")"))

ppp4 <- analyse_modobs2(Nmin_statistical,"ISAM_fNup","Nmin", type = "points",relative=TRUE)$gg +larger_size+
  labs(y = ~paste("Net ", minerlization[obs.], " (gN m"^-2,"yr"^-1,")")) +labs(x = ~paste("ISAM N ", uptake[pred.], " (gN m"^-2,"yr"^-1,")"))

ppp5 <- analyse_modobs2(Nmin_statistical,"JSBACH_fNup","Nmin", type = "points",relative=TRUE)$gg +larger_size+
  labs(y = ~paste("Net ", minerlization[obs.], " (gN m"^-2,"yr"^-1,")")) +labs(x = ~paste("JSBACH N ", uptake[pred.], " (gN m"^-2,"yr"^-1,")"))

ppp6 <- analyse_modobs2(Nmin_statistical,"LPX_fNup","Nmin", type = "points",relative=TRUE)$gg +larger_size+
  labs(y = ~paste("Net ", minerlization[obs.], " (gN m"^-2,"yr"^-1,")")) +labs(x = ~paste("LPX N ", uptake[pred.], " (gN m"^-2,"yr"^-1,")"))


plot_grid(ppp1,ppp3,ppp4,ppp5,ppp6,
          labels = c('(a)','(b)','(c)','(d)','(e)'),
          nrow=2,label_x = 0.9,label_y=0.92,label_size = 20)+white

ggsave(paste("~/data/output/newphy_figs6.jpg",sep=""), width = 20, height = 10)

#MAPE
mean(abs((NPP_statistical$TNPP_1-NPP_statistical$CLM_NPP)/NPP_statistical$TNPP_1),na.rm=TRUE) * 100
mean(abs((NPP_statistical$TNPP_1-NPP_statistical$ISAM_npp)/NPP_statistical$TNPP_1),na.rm=TRUE) * 100
mean(abs((NPP_statistical$TNPP_1-NPP_statistical$LPJ_NPP)/NPP_statistical$TNPP_1),na.rm=TRUE) * 100
mean(abs((NPP_statistical$TNPP_1-NPP_statistical$LPX_NPP)/NPP_statistical$TNPP_1),na.rm=TRUE) * 100
mean(abs((NPP_statistical$TNPP_1-NPP_statistical$ISBA_NPP)/NPP_statistical$TNPP_1),na.rm=TRUE) * 100
mean(abs((NPP_statistical$TNPP_1-NPP_statistical$CABLE_NPP)/NPP_statistical$TNPP_1),na.rm=TRUE) * 100
mean(abs((NPP_statistical$TNPP_1-NPP_statistical$CLASS_NPP)/NPP_statistical$TNPP_1),na.rm=TRUE) * 100
mean(abs((NPP_statistical$TNPP_1-NPP_statistical$JULES_NPP)/NPP_statistical$TNPP_1),na.rm=TRUE) * 100
mean(abs((NPP_statistical$TNPP_1-NPP_statistical$SDGVM_NPP)/NPP_statistical$TNPP_1),na.rm=TRUE) * 100
mean(abs((NPP_statistical$TNPP_1-NPP_statistical$JSBACH_NPP)/NPP_statistical$TNPP_1),na.rm=TRUE) * 100
mean(abs((NPP_statistical$TNPP_1-NPP_statistical$ORCHICNP_NPP)/NPP_statistical$TNPP_1),na.rm=TRUE) * 100
mean(abs((NPP_statistical$TNPP_1-NPP_statistical$ORCHIDEE_NPP)/NPP_statistical$TNPP_1),na.rm=TRUE) * 100
mean(abs((NPP_statistical$TNPP_1-NPP_statistical$pred_npp)/NPP_statistical$TNPP_1),na.rm=TRUE) * 100


mean(abs((Nmin_statistical$Nmin-Nmin_statistical$ORCHICNP_fNup)/Nmin_statistical$Nmin),na.rm=TRUE) * 100
mean(abs((Nmin_statistical$Nmin-Nmin_statistical$ISAM_fNup)/Nmin_statistical$Nmin),na.rm=TRUE) * 100
mean(abs((Nmin_statistical$Nmin-Nmin_statistical$JSBACH_fNup)/Nmin_statistical$Nmin),na.rm=TRUE) * 100
mean(abs((Nmin_statistical$Nmin-Nmin_statistical$LPX_fNup)/Nmin_statistical$Nmin),na.rm=TRUE) * 100
mean(abs((Nmin_statistical$Nmin-Nmin_statistical$pred_nuptake)/Nmin_statistical$Nmin),na.rm=TRUE) * 100


#visreg

a1 <- ~{
  p1a <- visreg(bp_model,"Tg_a",type="contrast")
  plot(p1a,ylab="Forest BP",xlab="Tg",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a2 <- ~{
  p1a <- visreg(bp_model,"PPFD_a",type="contrast")
  plot(p1a,ylab="Forest BP",xlab="ln PPFD",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a3 <- ~{
  p1a <- visreg(bp_model,"CNrt_a",type="contrast")
  plot(p1a,ylab="Forest BP",xlab="ln soil C/N",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a4 <- ~{
  p1a <- visreg(bp_model,"age_a",type="contrast")
  plot(p1a,ylab="Forest BP",xlab="ln age",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a5 <- ~{
  p1a <- visreg(bp_model,"fAPAR_a",type="contrast")
  plot(p1a,ylab="Forest BP",xlab="fAPAR",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a6 <- ~{
  p1a <- visreg(anpp_tnpp_model,"Tg_a",type="contrast")
  plot(p1a,ylab="logit ANPP/BP",xlab="Tg",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a7 <- ~{
  p1a <- visreg(anpp_tnpp_model,"PPFD_a",type="contrast")
  plot(p1a,ylab="logit ANPP/BP",xlab="ln PPFD",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a8 <- ~{
  p1a <- visreg(anpp_tnpp_model,"CNrt_a",type="contrast")
  plot(p1a,ylab="logit ANPP/BP",xlab="ln soil C/N",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a9 <- ~{
  p1a <- visreg(anpp_tnpp_model,"age_a",type="contrast")
  plot(p1a,ylab="logit ANPP/BP",xlab="ln age",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a10 <- ~{
  p1a <- visreg(anpp_leafnpp_model,"fAPAR_a",type="contrast")
  plot(p1a,ylab="logit leaf-NPP/ANPP",xlab="fAPAR",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a11 <- ~{
  p1a <- visreg(anpp_leafnpp_model,"PPFD_a",type="contrast")
  plot(p1a,ylab="logit leaf-NPP/ANPP",xlab="ln PPFD",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a12 <- ~{
  p1a <- visreg(anpp_leafnpp_model,"vpd_a",type="contrast")
  plot(p1a,ylab="logit leaf-NPP/ANPP",xlab="ln D",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a13 <- ~{
  p1a <- visreg(nre_model,"Tg_a",type="contrast")
  plot(p1a,ylab="logit NRE",xlab="Tg",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a14 <- ~{
  p1a <- visreg(nre_model,"vpd_a",type="contrast")
  plot(p1a,ylab="logit NRE",xlab="ln D",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a15 <- ~{
  p1a <- visreg(bp_grass_model,"Tg_a",type="contrast")
  plot(p1a,ylab="Grassland BP",xlab="Tg",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a16 <- ~{
  p1a <- visreg(bp_grass_model,"PPFD_a",type="contrast")
  plot(p1a,ylab="Grassland BP",xlab="ln PPFD",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}


plot_grid(a1,a2,a3,a4,a5,
          a6,a7,a8,a9,white,
          a10,a11,a12,white,white,
          a13,a14,white,white,white,
          a15,a16,white,white,white,
          nrow=5)+white

ggsave(paste("~/data/output/newphy_figs1.jpg",sep=""), width = 20, height = 20)

#table s1
names(all_maps)
for (i in 4:30){
  varname <- names(all_maps)[i]
  total_value <- round(sum(all_maps[,i]*conversion,na.rm=TRUE),2)
  print(varname)
  print(total_value)
}

#it seems that forest + grassland doesn't completely equal to pft (with decimal degree uncertainty 0.01)
#in this way, we re-calculate grassland = pft - forest 

#error propagation - needs updates fapar

# forest,grassland and 2-pft npp uncertainty
uncertainty_npp <- summary(bp_model)$sigma 

uncertainty_npp_grassland <- summary(bp_grass_model)$sigma 

# for all logit function model (mod_tnpp, mod_anpp, mod_lnpp, nre_model): a = 1/(1+exp(-b)), where a is the ratio (e.g. npp/gpp) and b is the regression (e.g. mod =tnpp)
#we calculate uncertainty of b firstly, based on each regression
#now we need to calculate deriative a / deriative b, based on a = 1/(1+exp(-b)). After calculation it is: exp(-b) / ( (1 + exp(-b)) ^2)
#forest anpp
mod_anpp_uncertainty <- summary(anpp_tnpp_model)$sigma
anpp_b <- summary(anpp_tnpp_model)$coefficients[1,1]+
  summary(anpp_tnpp_model)$coefficients[2,1] * log(CNrt$myvar) +
  summary(anpp_tnpp_model)$coefficients[3,1] * log(PPFD$myvar) + 
  summary(anpp_tnpp_model)$coefficients[4,1] * Tg$myvar+
  summary(anpp_tnpp_model)$coefficients[5,1] * log(age$myvar)

#deriative a (the ratio) = deriative b (the regression) * deriative a / deriative b
mod_anpp_uncertainty <- mod_anpp_uncertainty *  exp(-anpp_b) / ( (1 + exp(-anpp_b)) ^2)
summary(mod_anpp_uncertainty) #uncertainty of anpp/npp

uncertainty_anpp <- anpp_f * sqrt( (uncertainty_npp/npp_f)^2 +
                                     (mod_anpp_uncertainty/(anpp_f/npp_f))^2)

#uncertainty of bnpp. BNPP = NPP - ANPP
uncertainty_bnpp <- sqrt((uncertainty_npp)^2 + (uncertainty_anpp)^2)
summary(uncertainty_bnpp/bnpp_f)

## Uncertainty of LNPP/ANPP
mod_lnpp_uncertainty <- summary(anpp_leafnpp_model)$sigma
lnpp_b <- summary(anpp_leafnpp_model)$coefficients[1,1]+
  summary(anpp_leafnpp_model)$coefficients[2,1]* fAPAR$myvar +
  summary(anpp_leafnpp_model)$coefficients[3,1] * log(vpd$myvar) + 
  summary(anpp_leafnpp_model)$coefficients[4,1] * log(PPFD$myvar)

mod_lnpp_uncertainty <- mod_lnpp_uncertainty *  exp(-lnpp_b) / ( (1 + exp(-lnpp_b)) ^2)
summary(mod_lnpp_uncertainty) #uncerainty of leafnpp to anpp

#uncertainty of leaf npp = npp * (anpp/npp) * (leaf/anpp)
uncertainty_lnpp <- lnpp_f * sqrt( (uncertainty_npp/npp_f)^2 +
                                     (mod_anpp_uncertainty/(anpp_f/npp_f))^2 +
                                     (mod_lnpp_uncertainty/(lnpp_f/anpp_f))^2)
summary(uncertainty_lnpp/lnpp_f)

#uncertainty of wood npp = npp * (anpp/npp) * (1- leaf/anpp)
uncertainty_wnpp <- wnpp_f * sqrt( (uncertainty_npp/npp_f)^2 +
                                     (mod_anpp_uncertainty/(anpp_f/npp_f))^2 +
                                     (mod_lnpp_uncertainty/(1 - (lnpp_f/anpp_f)))^2)
summary(uncertainty_wnpp/wnpp_f)

#leaf n/c - ALREADY checked that inputted vcmax25_df + LMA calculated from (1) fortran and (2) R for predicting leaf n/c, could output the same prediction for leaf n/c Please note! using a.vcmax25.nc rather than annualvcmax25.nc. See difference in: https://www.notion.so/computationales/annualvcmax25-vs-a-vcmax25-282bddd63bba4d7b9cfcc75805b45964
n1 <- lmer(nmass_a~vcmax25_lma_a + (1|sitename_a)+(1|species_a))
mod_leafnc_uncertainty <-summary(n1)$sigma/0.47
summary(mod_leafnc_uncertainty)

#NRE
mod_nre_uncertainty <- summary(nre_model)$sigma
nre_b <- summary(nre_model)$coefficients[1,1]  +
  summary(nre_model)$coefficients[2,1] * Tg$myvar + summary(nre_model)$coefficients[3,1] * log(vpd$myvar)
mod_nre_uncertainty <- mod_nre_uncertainty *  exp(-nre_b) / ( (1 + exp(-nre_b)) ^2)
summary(mod_nre_uncertainty)

#now, leaf N flux = NPP * (ANPP/NPP) * (leafNPP/ANPP) * (leaf n/c) * (1-NRE)
uncertainty_lnf <- lnf_f * sqrt( (uncertainty_npp/npp_f)^2 +
                                   (mod_anpp_uncertainty/(anpp_f/npp_f))^2 +
                                   (mod_lnpp_uncertainty/(lnpp_f/anpp_f))^2 +
                                   (mod_leafnc_uncertainty/leafnc_f)^2 +
                                   (mod_nre_uncertainty/(1-nre_f))^2)
summary(uncertainty_lnf/lnf_f)

#now, uncertainty of wood n uptake flux (assuming wood c/n is a constant = 100)
uncertainty_wnf <- uncertainty_wnpp/100
summary(uncertainty_wnf/wnf_f)

#root N flux = uncertainty of bnpp * root N/C
uncertainty_bnf <- uncertainty_bnpp/94

#nuptake
uncertainty_nuptake <- sqrt(uncertainty_lnf^2 + uncertainty_wnf^2 + uncertainty_bnf^2) #N uptake

#grassland
uncertainty_npp_grassland
uncertainty_anpp_grassland <- uncertainty_npp_grassland*0.49
uncertainty_bnpp_grassland <- uncertainty_npp_grassland*0.51
uncertainty_lnf_grassland <- uncertainty_npp_grassland*0.49*(1-0.69)/18
uncertainty_bnf_grassland <- uncertainty_npp_grassland*0.51/41
uncertainty_nuptake_grassland <- sqrt(uncertainty_lnf_grassland^2 + uncertainty_bnf_grassland^2)

#all statistics
#forest
sum(uncertainty_npp*(forest_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_anpp*(forest_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_bnpp*(forest_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_lnpp*(forest_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_wnpp*(forest_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_lnf*(forest_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_wnf*(forest_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_bnf*(forest_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_nuptake*(forest_percent *conversion)*available_grid2,na.rm=TRUE)

#grassland
sum(uncertainty_npp_grassland*(grass_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_anpp_grassland*(grass_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_bnpp_grassland*(grass_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_lnf_grassland*(grass_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_bnf_grassland*(grass_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_nuptake_grassland*(grass_percent *conversion)*available_grid2,na.rm=TRUE)

#2-pft
sqrt(sum(uncertainty_npp*(forest_percent *conversion)*available_grid2,na.rm=TRUE)^2 + sum(uncertainty_npp_grassland*(grass_percent *conversion)*available_grid2,na.rm=TRUE)^2)
sqrt(sum(uncertainty_anpp*(forest_percent *conversion)*available_grid2,na.rm=TRUE)^2 + sum(uncertainty_anpp_grassland*(grass_percent *conversion)*available_grid2,na.rm=TRUE)^2)
sqrt(sum(uncertainty_bnpp*(forest_percent *conversion)*available_grid2,na.rm=TRUE)^2 + sum(uncertainty_bnpp_grassland*(grass_percent *conversion)*available_grid2,na.rm=TRUE)^2)
sqrt(sum(uncertainty_lnf*(forest_percent *conversion)*available_grid2,na.rm=TRUE)^2 + sum(uncertainty_lnf_grassland*(grass_percent *conversion)*available_grid2,na.rm=TRUE)^2)
sqrt(sum(uncertainty_bnf*(forest_percent *conversion)*available_grid2,na.rm=TRUE)^2 + sum(uncertainty_bnf_grassland*(grass_percent *conversion)*available_grid2,na.rm=TRUE)^2)
sqrt(sum(uncertainty_nuptake*(forest_percent *conversion)*available_grid2,na.rm=TRUE)^2 + sum(uncertainty_nuptake_grassland*(grass_percent *conversion)*available_grid2,na.rm=TRUE)^2)

#final nue - calculated from final results
#forest nue:
(55.74/0.79) 
#grassland nue:
(16.37/0.34)
#total nue:
(72.06/1.13) 
#forest nue uncertainty:
(55.74/0.79) * sqrt( (10.62/55.74)^2 +(0.23/0.79)^2)
#grassland nue uncertainty:
(16.37/0.34) * sqrt( (9.65/16.37)^2 +(0.15/0.34)^2)
#pft
(72.06/1.13) * sqrt( (14.35/72.06)^2 +(0.27/1.13)^2)


#partial residual figure
BP_dataset <- na.omit(NPP_forest[,c("tnpp_a","obs_age_a","observedfAPAR_a","soilCN_a","Tg_a","PPFD_a","vpd_a","site_a")])
#model1 <- stepwise(BP_dataset,"tnpp_a")
#model1[[1]]
#model1[[2]]
bp_model <- (lmer(tnpp_a~Tg_a+PPFD_a+soilCN_a+obs_age_a+observedfAPAR_a+(1|site_a),data=BP_dataset))
summary(bp_model)

#check obs. and pred. soil C/N
plot(NPP_forest$soilCN_a~NPP_forest$CNrt_a)

anpp_tnpp_dataset <- na.omit(NPP_forest[,c("anpp_tnpp_a","obs_age_a","observedfAPAR_a","soilCN_a","Tg_a","PPFD_a","vpd_a","site_a")])
dim(anpp_tnpp_dataset)
#model2 <- stepwise(anpp_tnpp_dataset,"anpp_tnpp_a")
#model2[[1]]
#model2[[2]]
anpp_tnpp_model <- (lmer(anpp_tnpp_a~Tg_a+PPFD_a+soilCN_a+obs_age_a++(1|site_a),data=anpp_tnpp_dataset))
summary(anpp_tnpp_model)

a1 <- ~{
  p1a <- visreg(bp_model,"Tg_a",type="contrast")
  plot(p1a,ylab="Forest BP",xlab="Tg",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a2 <- ~{
  p1a <- visreg(bp_model,"PPFD_a",type="contrast")
  plot(p1a,ylab="Forest BP",xlab="ln PPFD",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a3 <- ~{
  p1a <- visreg(bp_model,"soilCN_a",type="contrast")
  plot(p1a,ylab="Forest BP",xlab="ln soil C/N",line=list(col="white",lwd=0.01,lty=9), band=FALSE,
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a4 <- ~{
  p1a <- visreg(bp_model,"obs_age_a",type="contrast")
  plot(p1a,ylab="Forest BP",xlab="ln age",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a5 <- ~{
  p1a <- visreg(bp_model,"observedfAPAR_a",type="contrast")
  plot(p1a,ylab="Forest BP",xlab="fAPAR",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a6 <- ~{
  p1a <- visreg(anpp_tnpp_model,"Tg_a",type="contrast")
  plot(p1a,ylab="logit ANPP/BP",xlab="Tg",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a7 <- ~{
  p1a <- visreg(anpp_tnpp_model,"PPFD_a",type="contrast")
  plot(p1a,ylab="logit ANPP/BP",xlab="ln PPFD",line=list(col="white",lwd=0.01,lty=9), band=FALSE,
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a8 <- ~{
  p1a <- visreg(anpp_tnpp_model,"soilCN_a",type="contrast")
  plot(p1a,ylab="logit ANPP/BP",xlab="ln soil C/N",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a9 <- ~{
  p1a <- visreg(anpp_tnpp_model,"obs_age_a",type="contrast")
  plot(p1a,ylab="logit ANPP/BP",xlab="ln age",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}



plot_grid(a1,a2,a3,a4,a5,
          a6,a7,a8,a9,white,
          nrow=2)+white

ggsave(paste("~/data/output/newphy_figs1b.jpg",sep=""), width = 20, height = 10)