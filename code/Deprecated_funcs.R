#### These are mostly deprecated functions, some are slight version of currently used functions, but some are experimental and rabbit holes. 
### There are possible duplicates (or altered versions). 
#### These functions probably do not work on the latest version of things. 

### DEPRECATED AND OLD FUNCTIONS #####
rabbit_hole <- fuction(){reduct_five %>% filter(Accepted_Name%in% all_of(target_spp)) %>% 
    distinct(Accepted_Name,CellID,year,.keep_all = T) %>% 
    mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(10),"Agriculture",.x))) %>%
    mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(20),"Forest",.x))) %>%
    mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(30,40),"Grassland_Shrub",.x))) %>%
    mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(70,100,90),"Other",.x))) %>%
    mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(0),NA,.x))) %>% 
    # group_by(Accepted_Name) %>% select(CellID,year,Altitude,starts_with("LUC_")) %>% 
    mutate(year = as.numeric(sub("Year_","",year))) -> prueba
  
  target_year <- as.character(1982:2011)
  name_lul <- paste0("LUC_",target_year)
  name_meantas <- paste0("mean_tas_",target_year)
  name_maxtas <- paste0("max_tas_",target_year)
  name_mintas <- paste0("min_tas_",target_year)
  
  name_meanpr <- paste0("mean_pr_",target_year)
  name_maxpr <- paste0("max_pr_",target_year)
  name_minpr <- paste0("min_pr_",target_year)
  
  name_meancmi <- paste0("mean_cmi_",target_year)
  name_maxcmi <- paste0("max_cmi_",target_year)
  name_mincmi <- paste0("min_cmi_",target_year)
  
  prueba %<>% mutate(LUC_year = NA,.after=year) %>% mutate(minPR_year = NA,.after=year) %>% 
    mutate(maxPR_year = NA,.after=year) %>% mutate(meanPR_year = NA,.after=year) %>% mutate(minTAS_year = NA,.after=year) %>% 
    mutate(maxTAS_year = NA,.after=year) %>% mutate(meanTAS_year = NA,.after=year) %>% 
    mutate(minCMI_year = NA,.after=year) %>% 
    mutate(maxCMI_year = NA,.after=year) %>% mutate(meanCMI_year = NA,.after=year)
  
  for(i in 1:length(target_year)){
    cat(i,"\r")
    prueba %<>% rowwise() %>% mutate(LUC_year =  ifelse(!is.na(LUC_year),LUC_year,ifelse(year==target_year[i],!!rlang::sym(name_lul[i]),NA)),.after=year) %>% 
      mutate(LUC_year =  ifelse(!is.na(LUC_year),LUC_year,ifelse(year==target_year[i],!!rlang::sym(name_lul[i]),NA)),.after=year) %>%
      
      mutate(minPR_year =  ifelse(!is.na(minPR_year),minPR_year,ifelse(year==target_year[i],!!rlang::sym(name_minpr[i]),NA)),.after=year) %>% 
      mutate(maxPR_year =  ifelse(!is.na(maxPR_year),maxPR_year,ifelse(year==target_year[i],!!rlang::sym(name_maxpr[i]),NA)),.after=year) %>% 
      mutate(meanPR_year =  ifelse(!is.na(meanPR_year),meanPR_year,ifelse(year==target_year[i],!!rlang::sym(name_meanpr[i]),NA)),.after=year) %>% 
      
      mutate(minTAS_year =  ifelse(!is.na(minTAS_year),minTAS_year,ifelse(year==target_year[i],!!rlang::sym(name_mintas[i]),NA)),.after=year) %>% 
      mutate(maxTAS_year =  ifelse(!is.na(maxTAS_year),maxTAS_year,ifelse(year==target_year[i],!!rlang::sym(name_maxtas[i]),NA)),.after=year) %>% 
      mutate(meanTAS_year =  ifelse(!is.na(meanTAS_year),meanTAS_year,ifelse(year==target_year[i],!!rlang::sym(name_meantas[i]),NA)),.after=year) %>% 
      
      mutate(minCMI_year =  ifelse(!is.na(minCMI_year),minCMI_year,ifelse(year==target_year[i],!!rlang::sym(name_mincmi[i]),NA)),.after=year) %>% 
      mutate(maxCMI_year =  ifelse(!is.na(maxCMI_year),maxCMI_year,ifelse(year==target_year[i],!!rlang::sym(name_maxcmi[i]),NA)),.after=year) %>% 
      mutate(meanCMI_year =  ifelse(!is.na(meanCMI_year),meanCMI_year,ifelse(year==target_year[i],!!rlang::sym(name_meancmi[i]),NA)),.after=year)
    
  }
  
  prueba %<>% filter(year >= 1982 & year <=2011) %>% select(CellID:Accepted_Name, decimalLatitude:LUC_year) %>% rename(ind=Accepted_Name) %>% separate(ind,into = c("G","S","extra"),extra = "merge",sep=" ") %>% select(-extra) %>% unite("ind",G:S,sep="_")
  prueba %<>% mutate(tme = year - min(year))
  
  prueba %>% distinct(tme)
  prueba %>% distinct(Altitude)
  prueba %>% select(meanTAS_year) %>% mutate(meanTAS_year = meanTAS_year*.1 - 273.15)
  
  lmer(Altitude ~ 1 + tme + (0 + tme|ind), data = prueba %>% mutate(meanTAS_year = meanTAS_year*.1 - 273.15) %>%  filter(LUC_year %in% c("Altered_Forest","Agriculture","Forest")) %>%
         mutate(LUC_year = ifelse(LUC_year == "Agriculture","Altered_Forest",LUC_year)) %>% 
         mutate(LUC_year=factor(LUC_year,levels=c("Altered_Forest","Forest"))), control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4))) -> fit_raw
  broom.mixed::tidy(fit_raw)
}

## This one is not even from this paper, not sure what is doing here.
bi_plot <- function(){ 
  
  library(biscale)
  colors <- bi_class(esto %>% filter(!is.na(Completeness),!is.na(SR)), x = Completeness, y = SR, style = "fisher", dim = 4)
  colors %>% count(bi_class)
  legend <- bi_legend(pal = "PurpleGrn",dim = 4,
                      xlab = "Wallacean shortfall",ylab = "Species richness",
                      size = 10,rotate_pal = FALSE,flip_axes = FALSE,
                      breaks = bi_class_breaks(esto %>% filter(!is.na(Completeness),!is.na(SR)), x = Completeness, y = SR, style = "fisher", dim = 4))
  colors %>% 
    ggplot(aes(x=Longitude.x,y=Latitude.x,fill=bi_class)) +
    #geom_point(shape=15) +
    xlim(xlimits) + ylim(ylimits) + 
    geom_tile(colour="white") +
    #scale_fill_manual(values=MetBrewer::met.brewer("Cassatt2",n=16))
    bi_scale_fill(pal = "PurpleGrn", dim = 4,rotate_pal=FALSE,flip_axes = FALSE) +
    #bi_scale_color(pal = "GrPink", dim = 4) +
    theme(panel.background = element_blank(),panel.grid = element_blank(),
          legend.position = "none",
          legend.background = element_rect(fill = NA),
          legend.spacing.y = unit(c(0.2),"cm"),legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond"),
          axis.title = element_blank(),
          axis.text =  element_blank(),
          axis.line=element_blank(),
          axis.ticks = element_blank(),
          legend.margin=margin(t=-25),
          legend.key.size = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal") +
    geom_sf(data=roads,inherit.aes = FALSE,colour="black",fill=NA,size=0.2) +
    #cowplot::draw_plot(legend,x=-120,y=8,width = 12,height = 12) +
    inset_element(legend,left=0.0,right=0.6,bottom=-0.2,top=0.5,on_top=FALSE) +
    NULL
  
  colors <- bi_class(esto %>% filter(!is.na(SR)), x = age, y = SR, style = "jenks", dim = 4)
  colors %>% count(bi_class)
  legend <- bi_legend(pal = "PurpleGrn",
                      dim = 4,
                      xlab = "Temporal shortfall",
                      ylab = "Species richness",
                      size = 10,rotate_pal = FALSE,flip_axes = FALSE,
                      breaks = bi_class_breaks(esto %>% filter(!is.na(SR)), x = age, y = SR, style = "jenks", dim = 4))
  
  colors %>% 
    ggplot(aes(x=Longitude.x,y=Latitude.x,fill=bi_class)) +
    #geom_point(shape=15) +
    xlim(xlimits) + ylim(ylimits) + 
    geom_tile(colour="white") +
    #scale_fill_manual(values=MetBrewer::met.brewer("Cassatt2",n=16))
    bi_scale_fill(pal = "PurpleGrn", dim = 4,rotate_pal=FALSE,flip_axes = FALSE) +
    #bi_scale_color(pal = "GrPink", dim = 4) +
    theme(panel.background = element_blank(),panel.grid = element_blank(),
          legend.position = "none",
          legend.background = element_rect(fill = NA),
          legend.spacing.y = unit(c(0.2),"cm"),legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond"),
          axis.title = element_blank(),
          axis.text =  element_blank(),
          axis.line=element_blank(),
          axis.ticks = element_blank(),
          legend.margin=margin(t=-25),
          legend.key.size = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal") +
    geom_sf(data=roads,inherit.aes = FALSE,colour="black",fill=NA,size=0.2) +
    #cowplot::draw_plot(legend,x=-120,y=8,width = 12,height = 12) +
    inset_element(legend,left=0.0,right=0.6,bottom=-0.2,top=0.5,on_top=FALSE) +
    NULL
  
  colors <- bi_class(esto %>% filter(!is.na(SR)), x = seqs, y = SR, style = "jenks", dim = 4)
  colors %>% count(bi_class)
  legend <- bi_legend(pal = "PurpleGrn",
                      dim = 4,
                      ylab = "Species richness",
                      xlab = "Darwinian shortfall",
                      size = 10,rotate_pal = FALSE,flip_axes = FALSE,
                      breaks = bi_class_breaks(esto %>% filter(!is.na(SR)) , x = seqs, y = SR, style = "jenks", dim = 4))
  
  colors %>% 
    ggplot(aes(x=Longitude.x,y=Latitude.x,fill=bi_class)) +
    #geom_point(shape=15) +
    xlim(xlimits) + ylim(ylimits) + 
    geom_tile(colour="white") +
    #scale_fill_manual(values=MetBrewer::met.brewer("Cassatt2",n=16))
    bi_scale_fill(pal = "PurpleGrn", dim = 4,rotate_pal=FALSE,flip_axes = FALSE) +
    #bi_scale_color(pal = "GrPink", dim = 4) +
    theme(panel.background = element_blank(),panel.grid = element_blank(),
          legend.position = "none",
          legend.background = element_rect(fill = NA),
          legend.spacing.y = unit(c(0.2),"cm"),legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond"),
          axis.title = element_blank(),
          axis.text =  element_blank(),
          axis.line=element_blank(),
          axis.ticks = element_blank(),
          legend.margin=margin(t=-25),
          legend.key.size = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),legend.direction="horizontal") +
    geom_sf(data=roads,inherit.aes = FALSE,colour="black",fill=NA,size=0.2) +
    #cowplot::draw_plot(legend,x=-120,y=8,width = 12,height = 12) +
    inset_element(legend,left=0.0,right=0.6,bottom=-0.2,top=0.5,on_top=FALSE) +
    NULL
  
  
}


predicted_alts_deprecated <- function(data = "output/Splines_v.2.Rdata", 
                                      path_resample = "interim/splines_partial", 
                                      path_downsample = "interim", start = 1979, end = 2019, year_bins = 1,var_name = "SumSlopes_Min",output="output/Predicted_SlopesMin_v.2.RData") { 
  ### READ NULL MODELS AND OBSERVED, PREDICT AND MERGE ####
  ##  CHANGE HERE: SumSlopes SumSlopes_Max SumSlopes_Min SumSlopes_Range
  ##  CHANGE HERE: Slopes Slopes_Max Slopes_Min  Slopes_Range
  # This is deprecated: it relates to spline modelling.....
  load(data)
  ## Estimate models for the range: max - min 
  mean_alts <- mean_alts %>% mutate(Type="OBS",.before=data) #%>% rowwise(Accepted_Name) %>% mutate(Slopes_Range = list(splines_mods(data$YearCat,data$RollAlt_Max - data$RollAlt_Min,jit = NULL,fit = NULL)),.after=Slopes_Min) %>% mutate(SumSlopes_Range = list(splines_sums(Slopes_Range))) 
  
  list_arch <- list.files(path=here(path_resample),pattern = "_Partial_",full.names = T)
  list_arch_ds <- list.files(path=here(path_downsample),pattern = "_Partial_",full.names = T)
  
  all_models <- list()
  for (i in 1:length(list_arch_ds)){
    cat("Start processing file",i,"\n")
    cat(" ........ Reading...... ","\n")
    null_models <- mget(load(list_arch_ds[i]))
    test_null <- null_models %>% .[[1]] %>% ungroup() %>% mutate(Accepted_Name=sub("_sim_.*","",Accepted_Name)) %>% mutate(Type="NULL_downsamp",.before=data) #%>% rowwise(Accepted_Name) %>% mutate(Slopes_Range = list(splines_mods(data$YearCat,data$RollAlt_Max - data$RollAlt_Min,jit = NULL,fit = NULL)),.after=Slopes_Min) %>% mutate(SumSlopes_Range = list(splines_sums(Slopes_Range))) %>% ungroup()
    rm(null_models)
    if(i!=length(list_arch_ds)){test <- mean_alts}
    if(i==length(list_arch_ds)){
      test_null %>% distinct(Accepted_Name) %>% pull() -> whichas
      mean_alts %>% filter(Accepted_Name %in% all_of(whichas)) -> test ## filter to fit nulls
      test %>% bind_rows(.,test_null) -> test
    }
    # get significant trends and slopes
    lapply(test[[var_name]],"[",2,4) %>% unlist() -> sigs
    lapply(test[[var_name]],"[",2,1) %>% unlist() -> slopes
    test %>% ungroup() %>% mutate(Pvalue=sigs,Trend=slopes,.before=data) %>% rowwise() -> test
    lapply(test$data, min_data_NA,year_bins = year_bins) %>% do.call(rbind,.) %>% na_if(.,0) -> mininois_years ## get recorded years per species
    colnames(mininois_years) <- paste("Year",seq(start,end,year_bins),sep="_") 
    cat("...... calculating derivatives .......","\n")
    
    test %>% rowwise() %>% mutate(Predicted_Rates = list(splines_predict_bis(Slopes_Min,der=1,pred=seq(start,end,year_bins)))) -> preds ##  CHANGE HERE: Slopes Slopes_Max Slopes_Min  Slopes_Range
    preds %>% rowwise() %>% mutate(Predicted_Vals = list(splines_predict_bis(Slopes_Min,der=0,pred=seq(start,end,year_bins)))) -> preds ##  CHANGE HERE: Slopes Slopes_Max Slopes_Min  Slopes_Range
    #### bind_rows results
    
    lapply(preds$Predicted_Rates,"[",,2) -> fer 
    do.call(rbind,fer) -> fer
    colnames(fer) <- paste("Year",seq(start,end,year_bins),sep="_") 
    (fer*mininois_years) %>% as_tibble(.) -> fer
    preds %>% select(Accepted_Name,Type) %>% bind_cols(.,fer) -> preds_rates
    preds_rates %>% bind_cols(.,test %>% ungroup() %>% select(c(3:4))) -> preds_rates
    preds_rates %>% mutate(Derivative="First",.after=Type) -> preds_rates
    
    lapply(preds$Predicted_Vals,"[",,2) -> fer 
    do.call(rbind,fer) -> fer
    colnames(fer) <- paste("Year",seq(start,end,year_bins),sep="_") 
    (fer*mininois_years) %>% as_tibble(.) -> fer
    preds %>% select(Accepted_Name,Type) %>% bind_cols(.,fer) -> preds_vals
    preds_vals %>% bind_cols(.,test %>% ungroup() %>% select(c(3:4))) -> preds_vals
    preds_vals %>% mutate(Derivative="Raw",.after=Type)-> preds_vals
    preds <- rbind(preds_rates,preds_vals)
    cat("       DONE!     ",i,"\n")
    all_models[[i]] <- preds
    
  }
  all_models %>% do.call(rbind,.) -> down_models
  all_models <- list()
  for (i in 1:length(list_arch)){
    cat("Start processing file",i,"\n")
    cat(" ........ Reading...... ","\n")
    null_models <- mget(load(list_arch[i]))
    test_null <- null_models %>% .[[1]] %>% ungroup() %>% mutate(Accepted_Name=sub("_sim_.*","",Accepted_Name)) %>% mutate(Type="NULL_resamp",.before=data) #%>% rowwise(Accepted_Name) %>% mutate(Slopes_Range = list(splines_mods(data$YearCat,data$RollAlt_Max - data$RollAlt_Min,jit = NULL,fit = NULL)),.after=Slopes_Min) %>% mutate(SumSlopes_Range = list(splines_sums(Slopes_Range))) %>% ungroup()
    rm(null_models)
    test_null %>% distinct(Accepted_Name) %>% pull() -> whichas
    mean_alts %>% filter(Accepted_Name %in% all_of(whichas)) -> test ## filter to fit nulls
    test %>% bind_rows(.,test_null) -> test
    # get significant trends and slopes
    lapply(test[[var_name]],"[",2,4) %>% unlist() -> sigs
    lapply(test[[var_name]],"[",2,1) %>% unlist() -> slopes
    test %>% ungroup() %>% mutate(Pvalue=sigs,Trend=slopes,.before=data) %>% rowwise() -> test
    lapply(test$data, min_data_NA,year_bins = year_bins) %>% do.call(rbind,.) %>% na_if(.,0) -> mininois_years ## get recorded years per species
    colnames(mininois_years) <- paste("Year",seq(start,end,year_bins),sep="_") 
    cat("...... calculating derivatives .......","\n")
    
    test %>% rowwise() %>% mutate(Predicted_Rates = list(splines_predict_bis(Slopes,der=1,pred=seq(start,end,year_bins)))) -> preds ##  CHANGE HERE: Slopes Slopes_Max Slopes_Min  Slopes_Range
    preds %>% rowwise() %>% mutate(Predicted_Vals = list(splines_predict_bis(Slopes,der=0,pred=seq(start,end,year_bins)))) -> preds ##  CHANGE HERE: Slopes Slopes_Max Slopes_Min  Slopes_Range
    #### bind_rows results
    
    lapply(preds$Predicted_Rates,"[",,2) -> fer 
    do.call(rbind,fer) -> fer
    colnames(fer) <- paste("Year",seq(start,end,year_bins),sep="_") 
    (fer*mininois_years) %>% as_tibble(.) -> fer
    preds %>% select(Accepted_Name,Type) %>% bind_cols(.,fer) -> preds_rates
    preds_rates %>% bind_cols(.,test %>% ungroup() %>% select(c(3:4))) -> preds_rates
    preds_rates %>% mutate(Derivative="First",.after=Type) -> preds_rates
    
    lapply(preds$Predicted_Vals,"[",,2) -> fer 
    do.call(rbind,fer) -> fer
    colnames(fer) <- paste("Year",seq(start,end,year_bins),sep="_") 
    (fer*mininois_years) %>% as_tibble(.) -> fer
    preds %>% select(Accepted_Name,Type) %>% bind_cols(.,fer) -> preds_vals
    preds_vals %>% bind_cols(.,test %>% ungroup() %>% select(c(3:4))) -> preds_vals
    preds_vals %>% mutate(Derivative="Raw",.after=Type)-> preds_vals
    preds <- rbind(preds_rates,preds_vals)
    cat("       DONE!     ",i,"\n")
    all_models[[i]] <- preds
    
  }
  all_models %>% do.call(rbind,.) -> resamp_models
  resamp_models %>% count(Type)
  down_models %>% count(Type)
  all_models <- list(down_models,resamp_models)
  names(all_models) = c("down_models","resamp_models")
  save(all_models,file=here(output))
  
}
####

## Extract LUL data
lulc_extract_deprecated <- function(data="output/Historical_data_FULL.v.1.R",lulc_data = "LULC_850_2015") { 
  load(here(data))
  lista_luc <- list.files(path = here("data",lulc_data),pattern = ".tif$",full.names = T) %>% .[grep("1979_",.):grep("2019_",.)]
  lucs <- raster::stack(lista_luc)
  names(lucs) <- lista_luc %>% strsplit(.,"LULC_") %>% lapply(.,"[",3) %>% sub("_.*","",.) %>% paste("LUC",.,sep="_")
  xy <- reduct_five %>% select(decimalLongitude,decimalLatitude)
  raster::extract(lucs,xy) -> lucs_ext
  reduct_five <- as_tibble(lucs_ext) %>% bind_cols(reduct_five,.)
  save(reduct_five,file = here("output/Historical_data_FULL_LUL.v.1.R"))
}

## Create weight for sampling and do null models
## If trend = 0, then there is no shift in the distribution, and the only difference with the resampling would be that occurrences are sampled with weights according to baseline elevations.
## If trend != 0, then the normal distribution is shifted every year.
## If distorion = TRUE, the normal distribution is multiplied by a right-skewing beta distribution.
create_weight_deprecated <- function(x, y, shape = 3, trend = 3.37 * 1, skewed = FALSE) {
  get <- x - start_year - 1
  # generate a random beta distribution
  if(!skewed) {
    kappa = 2
    w = 1
    alpha = w*(kappa - 2) + 1
    beta =  (1-w)*(kappa - 2) + 1
    dens <- rbeta(10000,alpha + shape, beta + shape)
  }
  # the 0.011 skew-parameter effectively raise the probability from 0.5 to 0.52 by the end of the period
  if(skewed){
    kappa = 2 + (0.03*get)
    w = 1
    alpha = w*(kappa - 2) + 1
    beta =  (1-w)*(kappa - 2) + 1
    dens <- rbeta(10000,alpha + shape, beta + shape)}
  
  # transform beta distribution into elevation
  rr <- range_focal[,-1] + (trend * get)
  dens <- rr$Lower + (dens * (rr$Upper - rr$Lower))
  
  # get probability distribution and estimate probability for each elevation class
  dens_beta <- diff(ecdf(dens)(seq(0,5000,100)),lag=1,difference = 1)
  # The first estimate corresponds to 0-100 meters, and so on....
  
  # Set sampling weights by defining elevation classes.
  y %>% mutate(Weight = cut(Altitude,breaks=seq(0,5000,100),labels = dens_beta),.before=Altitude) %>% mutate(Weight = as.numeric(as.character(Weight)))
  
} 

# Generate null models with wieghted sampling. The weights are defined with base distributions defined from initial empirircal values; base parameters vary temporarily according to the observed average trends for cloud forests (coefficient from the linear mixed model).
alt_weighted_nulls <- function(data="output/Historical_data_FULL_LUL.v.1.R",
                               model_prefix = "Splines_Nulls_Shift",
                               n_sim=500,start_year=1979,records=100,workers=8) {
  cat("Loading data","\n")
  load(here(data))
  reduct_five %>% select(CellID:Accepted_Name,ID,class:family.x,decimalLatitude:BIOME) %>% 
    mutate(YearCat = as.numeric(sub("Year_","",year)),.after=1) %>% 
    filter(!is.na(Altitude),YearCat >= start_year) %>% 
    distinct(Accepted_Name,CellID,year,.keep_all = T) %>% group_by(Accepted_Name) -> esta
  esta %>% summarise(N=n()) -> tata
  target <- tata %>% filter(N >= records) %>% pull(Accepted_Name) %>% as.character(.)
  
  # Estimate elevation range of species: base parameters
  esta %>% filter(Accepted_Name %in% all_of(target)) %>% droplevels() %>% 
    group_by(Accepted_Name) %>%
    filter(YearCat <= 1990) %>% summarise(Lower=quantile(Altitude,c(0.025),na.rm=T),Upper=quantile(Altitude,c(0.975),na.rm=T)) -> range_r1
  
  # Estimate number of samples per year for all species
  esta %>% filter(Accepted_Name %in% all_of(target)) %>% group_by(YearCat,Accepted_Name) %>% summarise(N=n()) %>% arrange(Accepted_Name) -> n_years
  
  nulls_m <- list()
  count=0
  cat("Setting up multisessions","\n")
  plan(multisession,workers=workers)
  for (i in 1:length(target)){
    # Identify records (global) within altitudinal range of focal species
    cat(i,toupper(target[i]),"\n")
    count <- count + 1
    range_r1 %>% filter(Accepted_Name==target[i]) -> range_focal
    n_years %>% filter(Accepted_Name==target[i]) -> years_focal
    
    reduct_five %>% select(CellID:Accepted_Name,ID,class:family.x,decimalLatitude:BIOME) %>% 
      mutate(YearCat = as.numeric(sub("Year_","",year)),.after=1) %>% 
      filter(!is.na(Altitude),YearCat >= start_year) %>% 
      distinct(Accepted_Name,CellID,year,.keep_all = TRUE) %>%
      filter(Altitude >= range_focal[[2]],Altitude <= range_focal[[3]]) %>% 
      droplevels() %>% group_by(YearCat) %>% nest() %>% ungroup() %>% 
      filter(YearCat %in% years_focal[[1]]) -> r_null
    
    # Create vector of sampling weights for records as a function of altitude.
    r_null <- r_null %>% rowwise() %>% mutate(data = list(create_weight(YearCat,data))) %>% ungroup()
    
    cat("     Performing",n_sim,"simulations","\n")
    ## simulate data by resampling
    simulations <- 1:n_sim %>% future_map_dfr(function(j){
      
      # Sample all records with the same per-year intensity as the observed for the focal species
      # Use sampling weights per year.
      r_null %>%  ungroup() %>% mutate(N = years_focal$N, samp = map2(data, N, replace=FALSE,weight=Weight,sample_n)) %>% 
        select(-2) %>% unnest(samp) %>% droplevels() %>% group_by(YearCat) %>% 
        summarise(Variable=median(Altitude,na.rm=T),MaxVariable = quantile(Altitude,probs=0.975,na.rm=T),MinVariable=quantile(Altitude,probs=0.025,na.rm=T),n=n(),class=first(class)) %>% 
        mutate_if(is.double, ~na_if(., -Inf)) %>% 
        mutate_if(is.double, ~na_if(., Inf)) %>% 
        mutate(Accepted_Name=paste(target[i],"_sim_",j,sep="")) %>% arrange(Accepted_Name,YearCat) %>%
        mutate(RollAlt=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=mean,na.rm=TRUE,partial=TRUE,align="center")}),use.names = FALSE)) %>% 
        mutate(RollAlt_Max=unlist(tapply(MaxVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
        mutate(RollAlt_Min=unlist(tapply(MinVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
        
        mutate(RollAlt_sd=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=sd,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
        mutate(RollAlt_Max_sd=unlist(tapply(MaxVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = sd,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
        mutate(RollAlt_Min_sd=unlist(tapply(MinVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = sd,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
        
        mutate(Rolln=unlist(tapply(n,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=sum,na.rm=T,partial=T,align="center")}),use.names = F)) %>% filter(Rolln >= min_by_year) -> r1
    },.progress = TRUE,.options = furrr_options(seed = 123)) ## set seed to avoid warning!
    cat("     Summarizing","simulations","\n")
    simulations %>% ungroup() %>% nest_by(Accepted_Name) -> nulls_m[[i]] 
    names(nulls_m)[[i]] <- as.character(target[i])
    cat("------","\n")
    if(count == 100) { 
      cat("Saving models (partial)","\n")
      nulls_m %>% do.call(bind_rows,.) -> null_models
      nombre_arch <- paste(here("interim/splines_partial_v3",model_prefix),i,".Rdata",sep="")
      save(null_models,file=nombre_arch) 
      count=0
      rm(nulls_m,null_models)
      nulls_m <- list()
    }
    if(i==length(target)){
      cat("Saving models (final)","\n")
      nulls_m %>% do.call(bind_rows,.) -> null_models
      nombre_arch <- paste(here("interim/splines_partial_v3",model_prefix),i,".Rdata",sep="")
      save(null_models,file=nombre_arch) 
      rm(nulls_m,null_models)
    }
    
  }
  plan(sequential)
  ## revert to unparalleled session
}

nullWeighted_alts <- function(data = "output/Splines_v.3.Rdata", 
                              path_Wsample = "interim/splines_partial_v3",
                              start = 1979, end = 2019, year_bins = 1, 
                              output="output/SplinesNULLS_DistortedNoShift_v.3.Rdata") { 
  load(data)
  mean_alts <- mean_alts %>% mutate(Type="OBS",.before=data)
  list_arch <- list.files(path=here(path_Wsample),pattern = "Nulls_DistortedNoShift_",full.names = TRUE)
  
  all_models <- list()
  for (i in 1:length(list_arch)){
    cat("Start processing file",i,"\n")
    cat(" ........ Reading...... ","\n")
    null_models <- mget(load(list_arch[i]))
    test_null <- null_models %>% .[[1]] %>% ungroup() %>% select(Accepted_Name,data) %>% mutate(Accepted_Name=sub("_sim_.*","",Accepted_Name)) %>% mutate(Type="NULL_Weightedresamp",.before=data)
    rm(null_models)
    all_models[[i]] <- test_null %>% select(Accepted_Name,data)
  }
  all_models %>% do.call(rbind,.) -> Weighted_models
  all_models <- list(Weighted_models)
  names(all_models) = c("DistortedNoShift_models")
  save(all_models,file=here(output))
  
}

do_null <- function(){ nulls_m <- list()
count=0
reduct_five %>% select(CellID:Accepted_Name,ID,class:family.x,decimalLatitude:BIOME) %>% 
  mutate(YearCat = as.numeric(sub("Year_","",year)),.after=1) %>% 
  filter(!is.na(Altitude),YearCat >= start_year) %>% 
  distinct(Accepted_Name,CellID,year,.keep_all = T) %>% group_by(Accepted_Name) -> esta
esta %>% summarise(N=n()) -> tata
target <- tata %>% filter(N >= records) %>% pull(Accepted_Name) %>% as.character(.)
esta %>% filter(Accepted_Name %in% all_of(target)) %>% droplevels() %>% group_by(Accepted_Name) %>% summarise(Lower=quantile(Altitude,c(0.025),na.rm=T),Upper=quantile(Altitude,c(0.975),na.rm=T)) -> range_r1
# Estimate number of samples per year for all species
esta %>% filter(Accepted_Name %in% all_of(target)) %>% group_by(YearCat,Accepted_Name) %>% summarise(N=n()) %>% arrange(Accepted_Name) -> n_years
plan(multisession,workers=8)
for (i in 1:length(target)){ 
  # Identify records (global) within altitudinal range of focal species
  cat(i,toupper(target[i]),"\n")
  count <- count + 1
  range_r1 %>% filter(Accepted_Name==target[i]) -> range_focal
  n_years %>% filter(Accepted_Name==target[i]) -> years_focal
  reduct_five %>% select(CellID:Accepted_Name,ID,class:family.x,decimalLatitude:BIOME) %>% 
    mutate(YearCat = as.numeric(sub("Year_","",year)),.after=1) %>% 
    filter(!is.na(Altitude),YearCat >= start_year) %>% 
    distinct(Accepted_Name,CellID,year,.keep_all = T) %>%
    filter(Altitude >= range_focal[[2]],Altitude <= range_focal[[3]]) %>% 
    droplevels() %>% group_by(YearCat) %>% nest() %>% ungroup() %>% 
    filter(YearCat %in% years_focal[[1]]) -> r_null
  
  cat("     Performing",n_sim,"simulations","\n")
  ## simulate data by resampling
  simulations <- 1:n_sim %>% future_map_dfr(function(j){
    # Sample all records with the same per-year intensity as the observed for the focal species
    r_null %>%  mutate(N=years_focal$N,samp = map2(data, N, replace=F,sample_n)) %>% select(-2) %>% unnest(samp) %>% droplevels() %>% group_by(YearCat) %>% mutate(Accepted_Name=paste(target[i],"_sim_",j,sep=""))-> r1
  },.progress = TRUE,.options = furrr_options(seed = 123)) ## set seed to avoid warning!
  cat("     Summarizing","simulations","\n")
  simulations %>% distinct(CellID,Accepted_Name,year) -> nulls_m[[i]]
  
  
  names(nulls_m)[[i]] <- as.character(target[i])
  cat("------","\n")
  if(i==length(target)){
    cat("Saving models (final)","\n")
    nulls_m %>% do.call(bind_rows,.) -> null_occurrences
    nombre_arch <- paste(here("interim","Occurrences_ResampNulls"),".Rdata",sep="")
    save(nulls_m,file=nombre_arch) 
  }
  
}
plan(sequential)

nulls_md <- list()
esta %>% filter(Accepted_Name %in% all_of(target)) %>% count(year,Accepted_Name) %>%
  mutate(n = ifelse(n <= median(n), n, median(n)))-> N
esta %>% filter(Accepted_Name %in% all_of(target)) %>% group_by(year,Accepted_Name) %>% 
  nest() %>% left_join(.,N,by=c("year","Accepted_Name")) -> esta
plan(multisession,workers=8)
#system.time({
cat("     Performing",n_sim,"simulations","\n")
## simulate data by resampling
simulations <- 1:n_sim %>% future_map_dfr(function(j){
  esta %>% mutate(data=map2(data,n,replace=FALSE,sample_n)) %>% unnest(data) %>% mutate(Accepted_Name=paste(Accepted_Name,"_sim_",j,sep="")) %>% arrange(Accepted_Name,YearCat) -> r1
},.progress = TRUE,.options = furrr_options(seed = 123)) ## set seed to avoid warning!
cat("     Summarizing","simulations","\n")
simulations %>% distinct(CellID,Accepted_Name,year) -> nulls_md[[i]]
cat("Saving models (final)","\n")
nulls_md %>% do.call(bind_rows,.) -> null_md_occurrences
nombre_arch <- paste(here("interim","Occurrences_DownNulls"),".Rdata",sep="")
save(nulls_md,file=nombre_arch) 
#})
plan(sequential)
## revert ot unparallelized session

nulls_m %>% do.call(bind_rows,.) -> null_occurrences

null_md_occurrences %>% group_by(Accepted_Name) %>% select(CellID:Accepted_Name,year) %>% 
  #distinct(CellID,Accepted_Name,year) %>% 
  filter(grepl("sim_3$",Accepted_Name)) %>% 
  
  mutate(year=sub("Year_","",year)) %>% 
  left_join(., aver,by=c("CellID","year")) %>% 
  #filter(!is.na(Corrected_values)) %>% 
  group_by(Accepted_Name) %>%  
  mutate(year = as.numeric(sub("Year_","",year))) -> full_siono

full_siono %>%
  rename(Land_use = Corrected_values,Elevation=Altitude) %>% 
  select(Accepted_Name,year,Land_use,Elevation) %>% 
  filter(!is.na(Land_use)) %>%
  filter(Land_use %in% c("Altered_Forest","Agriculture","Forest")) %>%
  mutate(Land_use = ifelse(Land_use == "Agriculture","Altered_Forest",Land_use)) %>% 
  mutate(Land_use=factor(Land_use,levels=c("Altered_Forest","Forest"))) %>%
  rename(ind=Accepted_Name) %>% mutate(tme = year - min(year)) %>% 
  mutate(Elevation_Cat = cut(Elevation,breaks=seq(0,4000,1000),labels=seq(0,4000,1000)[-1], include.lowest=TRUE)) %>%
  separate(ind,into = c("G","S","extra"),extra = "merge",sep=" ") %>% select(-extra) %>% unite("ind",G:S,sep="_") %>% 
  #glm(Land_use ~ 1 + tme,family =  binomial, data = .) -> fit_null
  glmer(Land_use ~ 1 + tme + (0 + tme|ind),family =  binomial, data = ., nAGQ = 0,control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4)))-> fit_null
broom.mixed::tidy(fit_null)

} 

summary_LMMnulls_elevation <- function(){
  nulls="output/SplinesNULLS_v.2.Rdata"
  output="interim/resampSummary_LMM.Rdata"
  output_fig = "plots_tables/resampSummary_LMM.pdf"
  which=1
  n_sim=500
  end_year=2019
  start_year=1979
  load(here(nulls)) ### null models
  aver <- all_models[[which]] %>% group_by(Accepted_Name) %>% mutate(Type=paste("SIMS",1:n_sim,sep="_"),.before=data) %>% group_by(Accepted_Name,Type) %>% select(2:3) %>% 
    unnest(data)
  aver %<>% ungroup() %>% 
    select(RollAlt,YearCat,Accepted_Name,Type,RollAlt_error) %>% 
    mutate(YearCat=as.numeric(YearCat)) %>% rename_with(~c("y","tme","ind","sim","error"))
  aver %<>% filter(tme <= end_year) %>% mutate(year=tme) %>% mutate(tme = tme - start_year)
  otra <- aver %>% mutate(tme1 = ifelse(tme > 11,11,tme)) %>% mutate(tme2 = tme - 10) %>% mutate(tme2=ifelse(tme2 < 0,0,tme2)) %>% mutate(tme2 = ifelse(tme2 > 11,11,tme2)) %>% mutate(tme3 = tme - 20) %>% mutate(tme3=ifelse(tme3<0,0,tme3))
  ### rescale variables
  otra <- transform(otra,tme_sc=scale(tme),error_sc=scale(error),tme1_sc=scale(tme1),tme2_sc=scale(tme2),tme3_sc=scale(tme3)) %>% as_tibble()
  otra <- otra %>% filter(year<=2010)
  otra %<>% group_by(sim) %>% group_split()
  e_nulls <- list()
  sum_nulls <- list()
  for(i in 1:length(otra)){
    poruna <- otra[[i]]
    cat("Simulation -- ",i,"\n")
    cat("Fitting model","\r")
    fit1 <- lmer(y ~ 1 + tme + error + (1|ind) + (0 + tme|ind),data = poruna, control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4)))
    sum_nulls[[i]] <- broom.mixed::tidy(fit1)
    ef <-  coef(fit1)$ind %>% as_tibble(rownames = "ind") %>% select(tme,ind) %>% pivot_longer(cols = 1,names_to = "trend",values_to = "estimate" ) %>% mutate(trend = sub("tme","1979-2019",trend)) %>% mutate(sim=paste0("SIMS_",i))
    e_nulls[[i]] <- ef
  }
  do.call(rbind,e_nulls) %>% as_tibble() -> e_nulls
  trends <- lapply(sum_nulls,"[",2,4) %>% unlist()
  mean(trends)
  pp <- trends %>% as_tibble() %>% ggplot(aes(y=value)) +
    ylim(c(-2,3)) +
    ggdist::stat_dots(scale=0.9,fill = MetBrewer::met.brewer("Hiroshige",n=9,direction=1)[8], color= colorspace::darken(MetBrewer::met.brewer("Hiroshige",n=9,direction=1)[8],amount=0.2)) +
    geom_point(x=0.01,y=2.72,shape=21,size=2.5,color=colorspace::darken(MetBrewer::met.brewer("Hiroshige",n=9,direction=1)[1],amount=0.2),fill=MetBrewer::met.brewer("Hiroshige",n=9,direction=1)[1]) +
    theme(legend.position = "none",panel.background = element_blank(),axis.line = element_line(),legend.background = element_blank(),legend.key.height = unit(0.01,"cm"),legend.key.width = unit(0.01,"cm"),
          axis.title = element_text(family="EB Garamond"),
          axis.text = element_text(family="EB Garamond"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          legend.spacing.y = unit(c(1.2),"cm"),
          legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond")) +
    labs(x="",y=bquote("General elevation shift in meters"^-year)) + 
    NULL
  ggsave(filename = here(output_fig),plot = pp,height = 7,width = 8.5)
  save(e_nulls,file=here(output))
  
}

velocity_extract <- function(data="output/Historical_data_FULL.v.1.R",lulc_data = "GLASS-GLC") { 
  load(here(data))
  lista_luc <- list.files(path = here("data/climate_data_CHELSA/Velocidad"),pattern = ".tif$",full.names = T)
  lucs <- raster::stack(lista_luc)
  xy <- reduct_five %>% select(decimalLongitude,decimalLatitude)
  raster::extract(lucs,xy) -> lucs_ext
  reduct_five <- as_tibble(lucs_ext) %>% bind_cols(reduct_five,.)
  save(reduct_five,file = here("output/Historical_data_FULL_GLASS_VEL.v.1.R"))
  
}

## process and plot land-cover anc climate
do_LandCover_deprecated <- function() { 
  
  reduct_five %>% select(CellID:ind,year,starts_with("mean_cmi_") ) %>%
    mutate(year=as.numeric(sub("Year_","",year))) %>% 
    select(-Altitude,-year,-ind) %>% distinct(CellID,.keep_all = TRUE) %>% 
    pivot_longer(-c(1), names_to = c("Variable","year"),names_sep = "_cmi_",values_to = "Values") %>% 
    group_by(CellID,year) %>% 
    #mutate(TSI = (abs(Values[2] - min(Values))) / (abs(Values[2] - max((Values))))) %>% 
    rename(Climate = Values) %>% 
    select(-Variable) %>% 
    pivot_wider(names_from = year,names_prefix = "Climate",values_from = Climate) %>% 
    right_join(reduct_five %>% select(CellID:ind),.,by="CellID") -> aqui 
  
  reduct_five %>% 
    select(CellID:ind,year,starts_with("LUC_") ) %>%
    mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(10),"Agriculture",.x))) %>%
    mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(20),"Forest",.x))) %>%
    mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(30,40),"Grassland_Shrub",.x))) %>%
    mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(70,100,90),"Other",.x))) %>%
    mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(0),NA,.x))) %>% 
    select(-Altitude,-year,-ind) %>% distinct(CellID,.keep_all = TRUE) %>% 
    pivot_longer(cols = -c(1),names_to = "year",values_to = "Values") %>% 
    group_by(CellID,year) %>% 
    separate(year,sep="_",into=c("Var","year")) -> aver
  
  aver %>% mutate(New_values = ifelse(Values == "Forest" & lag(Values,n=1) != "Forest","Altered_Forest",Values)) %>% mutate(New_values = ifelse(is.na(New_values),Values,New_values))-> uno
  
  uno %<>% mutate(Corrected_values = New_values)
  uno %<>% nest()
  for(j in 1:nrow(uno)){
    cat(j,"\r")
    x <- uno$data[[j]]
    for(i in 2:nrow(x)){ 
      if(is.na(x$New_values[i])) next
      val <- ifelse(x$Corrected_values[i] == "Forest" & x$Corrected_values[i-1] != "Forest","Altered_Forest",x$Corrected_values[i])
      x$Corrected_values[i] <- val
    }
    x -> uno$data[[j]]
  }
  aver <- uno %>% unnest(data)
  rm(uno)
  aver %<>% mutate(Corrected_values = ifelse(is.na(Corrected_values),New_values,Corrected_values))
  
  
  aver %>% select(-Values,-New_values) %>% rename(LandCover = Corrected_values) %>% 
    pivot_wider(names_from = year,names_prefix = "LandCover",values_from = LandCover) %>% 
    right_join(reduct_five %>% select(CellID:ind),.,by="CellID") %>% 
    group_by(ind) %>% 
    select(-Altitude,-Var) %>% 
    pivot_longer(-c(ind,CellID),names_to = "year",values_to = "LandCover",names_prefix="LandCover") %>% 
    mutate(tme = as.numeric(year) - 1982) -> aqui
  
  left_join(aqui,nulas %>% select(ind,RSamp_SES,RSamp_SES_cat,DSamp_SES_cat,CAT,RSamp_SES,DSamp_SES,RSamp_ES,DSamp_ES,DSamp_Inter_Cat),by="ind") %>% 
    mutate(RSamp_SES_cat = cut(RSamp_SES,breaks=c(-10,-1.96,-1.28,0,1.28,1.96,10),labels=c("Negative","Negative","Weak","Weak","Positive","Positive"),include.lowest=TRUE)) %>% 
    mutate(DSamp_SES_cat = cut(DSamp_SES,breaks=c(-10,-1.96,-1.28,0,1.28,1.96,10),labels=c("Sensitive","Sensitive","Robust","Robust","Sensitive","Sensitive"),include.lowest=TRUE)) %>% 
    filter(DSamp_SES_cat %in% c("Robust")) -> fita
  
  reduct_five %>% select(CellID:Accepted_Name,year) %>% 
    #distinct(CellID,Accepted_Name,year) %>% 
    filter(Accepted_Name %in% all_of(target_spp)) %>% 
    mutate(year=sub("Year_","",year)) %>% 
    left_join(., aver,by=c("CellID","year")) %>% 
    rename(Altitude = Altitude.x) %>% 
    mutate(year = as.numeric(sub("Year_","",year))) %>%
    rename(Land_use = Corrected_values,Elevation=Altitude) %>% 
    select(Accepted_Name,year,Land_use,Elevation) %>% 
    filter(!is.na(Land_use)) %>%
    filter(Land_use %in% c("Altered_Forest","Agriculture","Forest")) %>%
    mutate(Land_use = ifelse(Land_use == "Agriculture","Altered_Forest",Land_use)) %>% 
    mutate(Land_use=factor(Land_use,levels=c("Altered_Forest","Forest"))) %>%
    rename(ind=Accepted_Name) %>% mutate(tme = year - min(year)) %>% 
    mutate(Elevation_Cat = cut(Elevation,breaks=seq(0,4000,1000),labels=seq(0,4000,1000)[-1], include.lowest=TRUE)) %>%
    separate(ind,into = c("G","S","extra"),extra = "merge",sep=" ") %>% select(-extra) %>% unite("ind",G:S,sep="_") %>% 
    left_join(., nulas,by="ind") %>%
    mutate(RSamp_SES_cat = cut(RSamp_SES,breaks=c(-10,-1.96,-1.28,0,1.28,1.96,10),labels=c("Negative","Negative","Weak","Weak","Positive","Positive"),include.lowest=TRUE)) %>% 
    mutate(DSamp_SES_cat = cut(DSamp_SES,breaks=c(-10,-1.96,-1.28,0,1.28,1.96,10),labels=c("Sensitive","Sensitive","Robust","Robust","Sensitive","Sensitive"),include.lowest=TRUE)) %>% filter(DSamp_SES_cat %in% c("Robust")) -> aqui
  
  fita %<>% filter(!is.na(LandCover)) %>%
    filter(LandCover %in% c("Altered_Forest","Agriculture","Forest")) %>%
    mutate(LandCover = ifelse(LandCover == "Agriculture","Altered_Forest",LandCover)) %>% 
    mutate(LandCover=factor(LandCover,levels=c("Altered_Forest","Forest")))
  
  
  fita %>% glmer(LandCover ~ 1 + tme + (1|ind) + (0 + tme|ind),family =  binomial, data = ., nAGQ = 0,control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4))) -> fit_lul_sum
  fita %>% glmer(LandCover ~ 1 + tme + tme:RSamp_SES_cat + (1|ind) + (0 + tme|ind),family =  binomial, data = ., nAGQ = 0,control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4))) -> fit_lul
  broom.mixed::tidy(fit_lul_sum)
  broom.mixed::tidy(fit_lul) %>% flextable::flextable() %>% flextable::save_as_docx(path = here("plots_tables/land_use_glmm.docx"))
  fit_lul_predict <- ggpredict(fit_lul, type = "random",terms=c("tme [all]","RSamp_SES_cat"),interval = "confidence") %>% tibble()
  fit_lul_predict %>% 
    ggplot(aes(x=x+1982,y=predicted,color=group,fill=group)) +
    geom_line() +
    geom_ribbon(aes(ymin = conf.low,ymax=conf.high),alpha=0.1,color="NA") +
    scale_color_met_d(name="Hiroshige",direction = -1)+
    scale_fill_met_d(name="Hiroshige",direction = -1) +
    geom_line(aes(x=x+1982,y=predicted),data = ggpredict(fit_lul_sum, type = "random",terms=c("tme [all]"),interval = "confidence"),color="black") +
    geom_ribbon(aes(ymin = conf.low,ymax=conf.high),alpha=0.1,color="NA",fill="black",data = ggpredict(fit_lul_sum, type = "random",terms=c("tme [all]"),interval = "confidence")) + 
    theme(legend.position = "right",panel.background = element_blank(),axis.line = element_line(),legend.background = element_blank(),
          axis.title = element_text(family="EB Garamond"),
          axis.text = element_text(family="EB Garamond"),
          #axis.text.x = element_blank(),
          legend.spacing.y = unit(c(1.2),"cm"),
          legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond")) +
    labs(x="Time",y="Predicted probabilities (Forest land-cover class)") +
    NULL
  
  ggsave(filename = "plots_tables/land_use_glmm.pdf",plot=last_plot())
  
  
}
do_climate_deprecated <- function() {
  data="output/Historical_data_FULL_GLASS.v.1.R"
  load(here(data))
  reduct_five %>% separate(Accepted_Name,into = c("G","S","extra"),extra = "merge",sep=" ") %>% 
    select(-extra) %>% unite("ind",G:S,sep="_") %>% 
    right_join(.,nulas %>% filter(SES == "SES_mid")%>% distinct(ind,.keep_all = TRUE),by="ind")
  
  reduct_five %>% select(CellID:ind,year,starts_with("max_tas_") ) %>%
    mutate(year=as.numeric(sub("Year_","",year))) %>% select(-Altitude,-year,-ind) %>% distinct(CellID,.keep_all = TRUE) %>% 
    pivot_longer(-c(1), names_to = c("Variable","year"),names_sep = "_pr_",values_to = "Values") %>% 
    group_by(CellID,year) %>% rename(Climate = Values) %>% select(-Variable) %>% 
    pivot_wider(names_from = year,names_prefix = "Climate",values_from = Climate) %>% 
    right_join(reduct_five %>% select(CellID:ind),.,by="CellID") -> aqui 
  
  aqui %>% group_by(ind) %>% 
    select(-Altitude) %>% 
    pivot_longer(-c(ind,CellID),names_to = "year",values_to = "Climate",names_prefix="Climate") %>% 
    mutate(tme = as.numeric(year) - 1979) -> aqui
  
  left_join(aqui,nulas %>% select(ind,RSamp_SES,RSamp_SES_cat,DSamp_SES_cat,CAT,RSamp_SES,DSamp_SES,RSamp_ES,DSamp_ES,DSamp_Inter_Cat),by="ind") %>% 
    mutate(RSamp_SES_cat = cut(RSamp_SES,breaks=c(-10,-1.96,-1.28,0,1.28,1.96,10),labels=c("Negative","Negative","Weak","Weak","Positive","Positive"),include.lowest=TRUE)) %>% 
    mutate(DSamp_SES_cat = cut(DSamp_SES,breaks=c(-10,-1.96,-1.28,0,1.28,1.96,10),labels=c("Sensitive","Sensitive","Robust","Robust","Sensitive","Sensitive"),include.lowest=TRUE)) %>% 
    filter(DSamp_SES_cat %in% c("Robust")) -> fita
  lmer(Climate ~ 1 + tme + tme:RSamp_SES_cat + (1|ind), data=fita, control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5))) -> fit_aver
  lmer(Climate ~ 1 + tme + (1|ind), data=fita, control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5))) -> fit_sum
  broom.mixed::tidy(fit_aver)
  broom.mixed::tidy(fit_aver) %>% flextable::flextable() %>% flextable::save_as_docx(path = here("plots_tables/cmi_glmm.docx"))
  ggpredict(fit_aver, type = "random",terms=c("tme [0,5,10,15,20,30,40]","RSamp_SES_cat"),interval = "confidence") %>% tibble() -> predict_fit
  predict_fit %>% 
    ggplot(aes(x=x,y=predicted,color=group,fill=group)) +
    geom_line() +
    geom_ribbon(aes(ymin = conf.low,ymax=conf.high),alpha=0.1,color="NA") +
    scale_color_met_d(name="Hiroshige",direction = -1)+
    scale_fill_met_d(name="Hiroshige",direction = -1) +
    geom_line(data=ggpredict(fit_sum, type = "random",terms=c("tme [all]"),interval = "confidence") %>% tibble(),aes(x=x,y=predicted),color="black",inherit.aes = FALSE) +
    theme(legend.position = "bottom",panel.background = element_blank(),axis.line = element_line(),legend.background = element_blank(),
          axis.title = element_text(family="EB Garamond"),
          axis.text = element_text(family="EB Garamond"),
          #axis.text.x = element_blank(),
          legend.spacing.y = unit(c(1.2),"cm"),
          legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond")) +
    labs(x="Time",y="Predicted mean tas") +
    NULL
  
  ggsave(filename = "plots_tables/cmi_glmm.pdf",plot=last_plot())
  
  
}

### SUMMARIZE ALTITUDE PER SPECIES AND MERGE WITH TRY #####
## First load and process try data
sum_species_deprecated <- function(TRY="output/growth_try.Rdata",pattern="Predicted_Slope",
                        quien = "Raw",output="output/Species_summaryRaw.RData"){
  load(here(TRY))
  list_archs <- list.files(path=here("output"),pattern = pattern,full.names = T)
  list_archs <- list_archs[grep("v.2",list_archs)]
  list_archs %>% sub(".*Predicted_","",.) %>% sub("_v.2","",.) %>% sub(".RData","",.) -> nombres
  species_means <- list()
  for (k in 1:length(list_archs)){ 
    load(list_archs[k])
    cat(list_archs[k],"\n")
    all_models$resamp_models %>% ungroup() %>% filter(Derivative==quien) %>% group_by(Accepted_Name) %>% summarise(across(starts_with("Year"), list(mean=mean,sd=sd))) -> sum_sims_resamp
    all_models$down_models %>% ungroup() %>% filter(Derivative==quien) %>% group_by(Accepted_Name) %>% summarise(across(starts_with("Year"), list(mean=mean,sd=sd))) -> sum_sims_down
    
    all_models$resamp_models %>% ungroup() %>% filter(Type=="OBS",Derivative==quien) -> observed
    
    sum_sims_resamp %>% pivot_longer(cols=-1,names_to = c("YearCat","Variable"),names_pattern = "Year_?(.*)_(.*)",values_to = "Estimate") %>% pivot_wider(names_from = Variable,values_from = Estimate) %>% rename(mean_resamp = mean,sd_resamp = sd) %>% mutate(YearCat=as.numeric(sub("Year_","",YearCat)))  -> sims_resamp_vals
    sum_sims_down %>% pivot_longer(cols=-1,names_to = c("YearCat","Variable"),names_pattern = "Year_?(.*)_(.*)",values_to = "Estimate") %>% pivot_wider(names_from = Variable,values_from = Estimate) %>% 
      rename(mean_down = mean,sd_down = sd)  %>% mutate(YearCat=as.numeric(sub("Year_","",YearCat))) -> sims_down_vals
    
    observed %>% pivot_longer(cols=-c(Accepted_Name,Type,Pvalue,Trend,Derivative),names_to = c("YearCat"),values_to = "Observed") %>% mutate(YearCat=as.numeric(sub("Year_","",YearCat))) %>% arrange(Accepted_Name,YearCat) %>% 
      left_join(.,sims_resamp_vals,by=c("Accepted_Name","YearCat")) %>%  
      left_join(.,sims_down_vals,by=c("Accepted_Name","YearCat")) %>% 
      mutate(SES_resamp = (Observed - mean_resamp) / sd_resamp) %>% 
      mutate(SES_down = (Observed - mean_down) / sd_down) -> to_plot
    
    all_models$resamp_models %>% ungroup() %>% filter(Type == "OBS",Derivative==quien) %>% pull(Accepted_Name) -> target_spp
    #### edit......
    to_plot %>% pull(Accepted_Name) -> especies
    especies %>% strsplit(.," ") %>% lapply(.,"[",1:2) %>% lapply(.,paste,collapse="_") %>% unlist() -> especies
    to_plot <- to_plot %>% mutate(binomial = especies,.after="Accepted_Name")
    especies <- unique(especies)
    growth_try %>% filter(Accepted_Name %in% all_of(especies)) %>% nest_by(Accepted_Name) -> target_try
    target_try %>% rowwise() %>% mutate(LifeForm=(Mode(data$OrigValueStr))) -> target_try
    names(target_try)[1] <- "binomial"
    left_join(to_plot,target_try,"binomial") %>% select(-data) -> to_plot
    to_plot %>% mutate(LifeForm=as.factor(LifeForm)) -> to_plot
    levels(to_plot$LifeForm) <- c("Climbers","Epiphytes","Herbs","Shrubs","Shrubs","Trees","Trees")
    ### NEED TO ESTIMATE BASELINE TREND AND THEN DO THE DIFFERENCE AND PLOT..... PER SPECIES
    species_means[[k]] <- to_plot %>% filter(Accepted_Name %in% all_of(target_spp)) %>% 
      mutate(SES_resamp_sig = SES_resamp >= 1.96 | SES_resamp <= -1.96) %>% 
      mutate(SES_down_sig = SES_down >= 1.96 | SES_down <= -1.96)
  }
  names(species_means) <- nombres
  save(species_means,file = here(output))
}

## Estimate granger causality
do_Granger_deprecated <- function(data="output/GrangerInput_Mid.Rdata",output="output/GrangerFinal_Mid.Rdata") {
  load(here(data))
  to_plot %>% ungroup() %>% pivot_wider(names_from = "Variable",values_from = "Series") %>% group_by(Resolved_ACCEPTED) %>% nest() %>% rowwise() %>% 
    ## Land use
    mutate(gci_luc = list(granger_mult(x.1 = data$LUC,y = data$Observed,n_x = 1))) %>% 
    ## Temperature
    mutate(gci_temps = list(granger_mult(x.1 = data$AnnTmean,x.2 = data$SeasTmean,x.3=data$MaxTmax,x.4 = data$MinTmin, y = data$Observed,n_x = 4))) %>% 
    ## Precipitation
    mutate(gci_precs = list(granger_mult(x.1 = data$AnnPrec,x.2 = data$SeasPrec,x.3=data$MaxPrec,x.4 = data$MinPrec, y = data$Observed,n_x = 4))) %>% 
    ## Overall climate
    mutate(gci_climate = list(granger_mult(data$AnnTmean,x.2 = data$SeasTmean,x.3=data$MaxTmax,x.4 = data$MinTmin, x.5 = data$AnnPrec,x.6 = data$SeasPrec,x.7=data$MaxPrec,x.8 = data$MinPrec,y = data$Observed,n_x = 8))) %>% 
    ## Overall climate | land use
    mutate(pgci_climate = list(granger_mult(data$AnnTmean,x.2 = data$SeasTmean,x.3=data$MaxTmax,x.4 = data$MinTmin, x.5 = data$AnnPrec,x.6 = data$SeasPrec,x.7=data$MaxPrec,x.8 = data$MinPrec,z = data$LUC, y = data$Observed,n_x = 8))) %>% 
    ## Climate plus land use
    mutate(gci_all = list(granger_mult(data$AnnTmean,x.2 = data$SeasTmean,x.3=data$MaxTmax,x.4 = data$MinTmin, x.5 = data$AnnPrec,x.6 = data$SeasPrec,x.7=data$MaxPrec,x.8 = data$MinPrec,x.9 = data$LUC, y = data$Observed,n_x = 5))) -> final
  save(final,file=here(output))
}

process_climate_series <- function(data_lul="output/Historical_data_FULL_LUL.v.1.R",
                                   output="ClimateSeries_v.2.RData",cual_n =1){
  list_archs <- list.files(path=here("output"),pattern="Predicted_Slope",full.names = T) %>% 
    .[grep("v.2",.)] 
  list_archs %>% sub(".*Predicted_","",.) %>% sub(".RData","",.) %>% sub("_v.2","",.) -> nombres
  load(here(data_lul))
  load(list_archs[cual_n])
  cat(nombres[cual_n],"\n")
  reduct_five %>% select(CellID,decimalLongitude,decimalLatitude,Altitude,starts_with("mean_"),starts_with("max_"),starts_with("min_")) %>% 
    distinct(CellID,.keep_all = T) %>% droplevels() %>% 
    pivot_longer(cols = -c(1:4),names_to = "year",values_to = "Values") -> aqui
  aqui %>% mutate(year=sub("^mean_","mean.",year)) %>% mutate(year=sub("^max_","max.",year)) %>% 
    mutate(year=sub("^min_","min.",year)) %>%
    separate(year,sep="_",into=c("Var","year")) %>% arrange(CellID,year) -> aqui
  
  aqui %>% group_by(Var) %>% arrange(CellID,year,.by_group = TRUE) %>% 
    mutate(Roll_est = unlist(tapply(Values,CellID,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean, na.rm=T, partial=T, align="center")}),use.names = F)) %>% 
    rename(YearCat=year) %>% group_by(CellID) %>% nest() %>% rowwise() -> aqui
  
  aqui %>% ungroup() %>%
    mutate(data=dif_clim_2(data,index = 6:7,flex_pt = 1990)) %>% unnest(data) -> aqui
  save(aqui,file = here("interim",output))
}

process_lul_series <- function(data_lul="output/Historical_data_FULL_LUL.v.1.R",
                               output="LULSeries_v.2.RData",cual_n =1){
  load(here(data_lul))
  list_archs <- list.files(path=here("output"),pattern="Predicted_Slope",full.names = T) %>% .[grep("v.2",.)] 
  list_archs %>% sub(".*Predicted_","",.) %>% sub(".RData","",.) %>% sub("_v.2","",.) -> nombres
  load(list_archs[cual_n])
  cat(nombres[cual_n],"\n")
  
  reduct_five %>% distinct(Accepted_Name,CellID,year,.keep_all = T) %>% group_by(Accepted_Name) %>% select(CellID,year,starts_with("LUC_")) %>% mutate(year=sub("Year_","",year)) %>% rename(YearCat=year) %>% 
    #filter(Accepted_Name=="Asplenium abscissum Willd.") %>% 
    ungroup() %>% select(-Accepted_Name) %>% 
    pivot_longer(cols = -c(1:2),names_to = "year",values_to = "Values") %>%  separate(year,sep="_",into=c("Var","year")) %>% arrange(CellID,year) %>% group_by(CellID) -> aver
  
  aver %<>% mutate(Values=case_when(Values %in% c(1:6, 9, 12) ~"AU",Values %in% c(7)~"PF",Values %in% c(10,11)~"SV",Values %in% c(8)~"NF")) %>% ungroup() 
  
  category = "AU"
  aver %<>% group_by(CellID) %>% mutate(Rolllul_n = unlist(tapply(Values, CellID, function(x) {(zoo::rollapply(as_vector(Values),width = window, FUN = parse_count, partial=T, align="center"))},simplify = TRUE)),Roll_AU = unlist(tapply(Values, CellID, function(x) {(zoo::rollapply(as_vector(Values),width = window, FUN = parse_freq, partial=T, align="center"))},simplify = TRUE)))
  
  category = "PF"
  aver %<>% group_by(CellID) %>% mutate(Roll_PF = unlist(tapply(Values, CellID, function(x) {(zoo::rollapply(as_vector(Values),width = window, FUN = parse_freq, partial=T, align="center"))},simplify = TRUE)))
  
  category = "SV"
  aver %<>% group_by(CellID) %>% mutate(Roll_SV = unlist(tapply(Values, CellID, function(x) {(zoo::rollapply(as_vector(Values),width = window, FUN = parse_freq, partial=T, align="center"))},simplify = TRUE)))
  
  category = "NF"
  aver %<>% group_by(CellID) %>% mutate(Roll_NF = unlist(tapply(Values, CellID, function(x) {(zoo::rollapply(as_vector(Values),width = window, FUN = parse_freq, partial=T, align="center"))},simplify = TRUE)))
  
  save(aver,file = here("interim",output))
}

not_ready_yet <- function(){
  variables = "Min_tmin"
  load(here("output_data/Historical_dataClimate_full.R"))
  reduct_five %>% filter(decimallatitude > lat_limits[1] & decimallongitude < lon_limits[1],
                         !(decimallatitude > lat_limits[2] | decimallongitude < lon_limits[2]),
                         !(decimallatitude > lat_limits[3] & decimallongitude > lon_limits[3]),
                         !CellID%in%all_of(non_id),AltRas<alt_limit) %>% dplyr::select(CellID,x,y,AltRas,starts_with(variables)) %>% distinct(CellID,.keep_all = T) %>% droplevels() %>% pivot_longer(cols = -c(1:4),names_to = "year",values_to = "Variable") %>% arrange(CellID,year) %>% mutate(Roll_est = unlist(tapply(Variable,CellID,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean, na.rm=T, partial=T, align="center")}),use.names = F))  %>% mutate(YearCat=as.numeric(sub(paste(variables,"_",sep=""),"",year))) %>% group_by(CellID) -> precs
  precs
  precs %>% nest_by() %>% mutate(Slopes = list(splines_mods(data$YearCat,data$Roll_est,jit=NULL,fit=NULL))) %>% mutate(SumSlopes = list(splines_sums(Slopes))) -> precs
  reduct_five %>% distinct(CellID,AltRas) %>% left_join(precs,.,by="CellID") -> precs
  precs %>% pull(AltRas) %>% cut(.,breaks=seq(0,alt_limit,belts),labels=seq(0,alt_limit,belts)[-1], include.lowest=T) -> alt_cats ### aggregate records into altitudinal belts
  precs <- precs %>% ungroup() %>% mutate(AltCats = alt_cats)
  precs %>% rowwise() %>% mutate(Predicted = list(splines_predict_bis(Slopes,der=deriv,pred=seq(start_year,2016,year_bins)))) -> preds
  #### bind_rows results
  lapply(preds$Predicted,"[",,2) -> fer 
  do.call(rbind,fer) -> fer 
  colnames(fer) <- paste("Year",seq(start_year,2016,year_bins),sep="_") 
  fer %>% as_tibble(.) ->fer
  preds %>% dplyr::select(CellID,AltCats) %>% bind_cols(.,fer) -> preds
  preds %>% pivot_longer(cols=-c(1:2),names_to = "YearCat",values_to = "Rate") %>% mutate(YearCat=as.numeric(sub("Year_","",YearCat)),AltCats=as.numeric(as.character(AltCats))) %>% group_by(AltCats,YearCat) %>% summarise(MeanRate=mean(Rate,na.rm=T)) -> preds
  if(prec) {direct=1; opcion_pal = "BrBG"}
  if(!prec) {direct=-1; opcion_pal="RdBu"}
  preds$MeanRate %>% range(na.rm=T) -> limits
  n_breaks=11
  brks <- seq(limits[1],limits[2], by= (limits[2]-limits[1])/n_breaks)[-c(1,n_breaks+1)]
  if(variables=="Seas_prec"){ brks <- c(-15,-10,-7,-4,0,4,7,15,25,35); direct=-1}
  if(variables=="Seas_tmean"){ brks <- c(-4.5,-3.5,-1.5,-0.5,0,0.5,1.5,2.5,3.5,4.5); direct=-1}
  if(variables=="Ann_tmean"){ brks <- c(-9,-7,-5,-1,0,5,10,15,20,25); direct=-1}
  if(variables=="Ann_prec"){ brks <- c(-135,-80,-30,-10,0,30,80,180,230,340); direct=1}
  preds %>% ggplot(aes(x=YearCat, y=AltCats, fill= MeanRate)) + 
    geom_tile() + labs(title = "Climate Change in cloud forests",
                       subtitle = paste("Rate of change (", "\u0394",") in ",code,sep=""),
                       y="Altitude",x="",caption = "Ramírez-Barahona et al.") + theme(panel.background = element_blank(),axis.ticks.y = element_blank()) + ylim(0,alt_limit+200) +
    scale_fill_fermenter(palette = opcion_pal,direction = direct,type="div",n.breaks=8) +
    # viridis::scale_fill_viridis(option=opcion_viridis,direction = direct,discrete=FALSE) + 
    NULL
  
  preds %>% rename_with(~c("AltCats","YearCat","Min_tmin")) %>% full_join(.,aa,by=c("AltCats","YearCat")) -> aa
  here("ClimateRAW_by_AltCat.R")
  save(aa,file=here("ClimateRAW_by_AltCat.R"))
  
  #### #####
  data="output_data/ClimateRAW_by_AltCat.R"
  inflect_point = 1976
  output="plots/ClimateSeries_TempRAW.pdf"
  variable = "tm"
  viridis_option = "F"
  load(here(data))
  aa %>% mutate(Shift=case_when(YearCat < inflect_point ~ "Pre",YearCat >= inflect_point ~ "Post")) %>% 
    group_by(Shift,AltCats) %>% summarize_at(2:9, median,na.rm=T) %>% 
    pivot_wider(names_from = c(Shift),values_from = c(Seas_tmean,Ann_prec,Ann_tmean,Seas_prec,Min_prec,Max_prec,Min_tmin,Max_tmax)) %>%
    pivot_longer(cols=-1,names_to = c("Variable","Shift","None"),names_pattern = "(.*)_(.*)_?(.*)",values_to = "Estimate") %>% select(-None) %>% pivot_wider(names_from = c(Variable,Shift),values_from = Estimate) %>% group_by(AltCats)
}
year_change <- function(x) {x %>% select(tme,Land_use) -> xx 
  index <- which(xx$Land_use != lag(xx$Land_use ,1))
  if(length(index)==0) return(0) else x$tme[index]
}
{
  chi %<>% mutate(year_change = year_change(data))
  
  reduct_five %>% filter(Accepted_Name %in% all_of(target_spp)) %>% select(CellID,Accepted_Name,Altitude,decimalLongitude,decimalLatitude) %>% distinct(CellID,Accepted_Name,.keep_all = TRUE) %>% left_join(.,chi %>% ungroup() %>% select(-data),by="CellID") %>% 
    #filter(Altitude >= 1500 & Altitude <=2500) %>% 
    mutate(Alt_cat = cut(Altitude,breaks=seq(0,4000,500),labels=seq(0,4000,500)[-1])) %>% 
    filter(!is.na(Alt_cat)) %>% 
    group_by(Alt_cat) %>% 
    summarise(Year = median(year_change,na.rm=TRUE))
  
  count(Manke_Cat) %>% mutate(Manke_Cat = ifelse(is.na(Manke_Cat),"STABLE",Manke_Cat))  %>% 
    mutate(Proportion = n / sum(n) * 100) %>% 
    filter(Manke_Cat%in% c("STABLE","WEAK")) %>% 
    mutate(Proportion = sum(Proportion)) %>% distinct(Alt_cat,.keep_all = TRUE) %>% 
    # filter(is.na(Manke_Cat)) %>% 
    ggplot(aes(y=Proportion,x=Alt_cat)) +
    geom_segment(aes(xend = Alt_cat,yend=0)) +
    geom_point(size=3) +
    NULL
  
  
  
  ggplot(aes(x=manken_lul,y=Alt_cat)) + 
    geom_density_ridges_gradient(rel_min_height = 0.001) +
    scale_fill_stepsn(colors=(MetBrewer::met.brewer("OKeeffe2",direction=1)),name="Mann-Kendall Z",n.breaks=10) + 
    theme(panel.background = element_blank(),panel.grid = element_blank(),
          axis.ticks.y  = element_blank(),
          panel.border = element_rect(fill=NA),
          legend.background = element_rect(fill = NA),
          legend.key.size = unit(0.7, "cm"),legend.key.width = unit(0.6,"cm"), legend.text = element_text(family="EB Garamond",size=8),legend.title = element_text(family="EB Garamond"),legend.spacing = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5)) +
    geom_vline(xintercept = c(1.96)) +
    NULL
}

### PLOT THINGS #####
plot_climatic_series_deprecated <- function(data="interim/ClimateSeries_v.2.RData",alt_limit = 4000){
  load(here(data))
  aqui %>% filter(!is.infinite(Roll_est)) %>%
    select(CellID,Var,YearCat,Roll_est) %>% filter(YearCat%in%1979:1989) %>% 
    mutate(Baseline=ifelse(YearCat%in%c(1979:1989),"Baseline",ifelse(YearCat%in%c(2000:2010),"Finish","Inter"))) %>% 
    group_by(CellID,Var,Baseline) %>% 
    summarise(Roll_est=mean(Roll_est,na.rm=TRUE)) %>% ungroup() %>% 
    pivot_wider(names_from = "Var",values_from = "Roll_est") -> pcas
  
  pcas %>% filter(Baseline=="Baseline") %>% drop_na() %>% select(where(is.double)) -> base.matrix
  cmdscale(dist(base.matrix))
  
  
  
  
  meancmi <- aqui %>% filter(!is.infinite(Roll_est)) %>% filter(Var=="mean.cmi") %>% mutate(AltCat=cut(Altitude,breaks=seq(0,alt_limit+100,250),labels=seq(0,alt_limit+100,250)[-1], include.lowest=T)) %>% group_by(YearCat,AltCat) %>% 
    summarise(mean=mean(Dif_Roll_est,na.rm=T)) %>% filter(!is.na(AltCat)) %>% mutate(All=1) %>% 
    ggplot(aes(x=as.numeric(YearCat),y=(AltCat),fill=mean*0.1)) + 
    labs(x="",y="Elevation (meters)",title="Mean Climate Moisture Index") +
    geom_tile() +
    scale_fill_gradientn(colours= colorspace::darken(MetBrewer::met.brewer("Demuth",direction=-1)[c(3:5,6:10)],amount=0),
                         n.breaks=10,
                         guide=guide_legend(reverse = FALSE,nrow = 1,title.position="top",label.position="bottom", title.hjust = 0),name="") +
    scale_color_gradientn(colours= colorspace::darken(MetBrewer::met.brewer("Demuth",direction=-1)[c(3:5,6:10)],amount=0),n.breaks=10,guide=guide_legend(reverse = FALSE,nrow = 1,title.position="top",label.position="bottom", title.hjust = 0),name="") +
    theme(legend.position = "bottom",
          panel.background = element_blank(),
          axis.line = element_line(),
          legend.background = element_blank(),
          #legend.key.height = unit(0.01,"cm"),legend.key.width = unit(0.01,"cm"),
          axis.title = element_text(family="EB Garamond"),
          axis.text = element_text(family="EB Garamond"),
          axis.ticks.y=element_blank(),
          #axis.text.y=element_blank(),
          #axis.line.y=element_blank(),
          legend.spacing.y = unit(c(0.2),"cm"),
          legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond"),
          legend.margin=margin(t=-25),
          title = element_text(family="EB Garamond")) +
    NULL
  
  meantas <- aqui %>% filter(!is.infinite(Roll_est)) %>% filter(Var=="mean.tas") %>% mutate(AltCat=cut(Altitude,breaks=seq(0,alt_limit+100,250),labels=seq(0,alt_limit+100,250)[-1], include.lowest=T)) %>% group_by(YearCat,AltCat) %>% 
    summarise(mean=mean(Dif_Roll_est,na.rm=T)) %>% filter(!is.na(AltCat)) %>% mutate(All=1) %>% 
    ggplot(aes(x=as.numeric(YearCat),y=(AltCat),fill=mean*0.1)) + 
    labs(x="",y="Elevation (meters)",title="Mean daily air temperature") +
    geom_tile() +
    scale_fill_gradientn(colours= colorspace::darken(MetBrewer::met.brewer("Demuth",direction=-1)[c(5,6:10)],amount=0),limits=c(-0.1,1.0),
                         #n.breaks=10,
                         breaks=c(-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
                         guide=guide_legend(reverse = FALSE,nrow = 1,title.position="top", label.position="bottom",title.hjust = 0),name="") +
    scale_color_gradientn(colours= colorspace::darken(MetBrewer::met.brewer("Demuth",direction=-1)[c(5,6:10)],amount=0),#n.breaks=10,
                          limits=c(-0.1,1.0),
                          breaks=c(-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
                          guide=guide_legend(reverse = FALSE,nrow = 1,title.position="top", label.position="bottom",title.hjust = 0),name="") +
    theme(legend.position = "bottom",
          panel.background = element_blank(),
          axis.line = element_line(),
          legend.background = element_blank(),
          #legend.key.height = unit(0.01,"cm"),legend.key.width = unit(0.01,"cm"),
          axis.title = element_text(family="EB Garamond"),
          axis.text = element_text(family="EB Garamond"),
          axis.ticks.y=element_blank(),
          #axis.text.y=element_blank(),
          #axis.line.y=element_blank(),
          legend.spacing.y = unit(c(0.2),"cm"),
          legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond"),
          legend.margin=margin(t=-25),
          title = element_text(family="EB Garamond")) +
    NULL
  
  meantasmax <- aqui %>% filter(!is.infinite(Roll_est)) %>% filter(Var=="mean.tasmax") %>% mutate(AltCat=cut(Altitude,breaks=seq(0,alt_limit+100,250),labels=seq(0,alt_limit+100,250)[-1], include.lowest=T)) %>% group_by(YearCat,AltCat) %>% 
    summarise(mean=mean(Dif_Roll_est,na.rm=T)) %>% filter(!is.na(AltCat)) %>% mutate(All=1) %>% 
    ggplot(aes(x=as.numeric(YearCat),y=(AltCat),fill=mean*0.1)) + 
    labs(x="",y="Elevation (meters)",title="Mean daily maximum air temperature") +
    geom_tile() +
    scale_fill_gradientn(colours= colorspace::darken(MetBrewer::met.brewer("Demuth",direction=-1)[c(5,6:10)],amount=0),limits=c(-0.1,1.0),
                         #n.breaks=10,
                         breaks=c(-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
                         guide=guide_legend(reverse = FALSE,nrow = 1,title.position="top",label.position="bottom", title.hjust = 0),name="") +
    scale_color_gradientn(colours= colorspace::darken(MetBrewer::met.brewer("Demuth",direction=-1)[c(5,6:10)],amount=0),#n.breaks=10,
                          limits=c(-0.1,1.0),
                          breaks=c(-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),guide=guide_legend(reverse = FALSE,nrow = 1,title.position="top",label.position="bottom", title.hjust = 0),name="") +
    theme(legend.position = "bottom",
          panel.background = element_blank(),
          axis.line = element_line(),
          legend.background = element_blank(),
          #legend.key.height = unit(0.01,"cm"),legend.key.width = unit(0.01,"cm"),
          axis.title = element_text(family="EB Garamond"),
          axis.text = element_text(family="EB Garamond"),
          axis.ticks.y=element_blank(),
          #axis.text.y=element_blank(),
          #axis.line.y=element_blank(),
          legend.spacing.y = unit(c(0.2),"cm"),
          legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond"),
          legend.margin=margin(t=-25),
          title = element_text(family="EB Garamond")) +
    NULL
  
  meanhurs <- aqui %>% filter(!is.infinite(Roll_est)) %>% filter(Var=="mean.hurs") %>% mutate(AltCat=cut(Altitude,breaks=seq(0,alt_limit+100,250),labels=seq(0,alt_limit+100,250)[-1], include.lowest=T)) %>% group_by(YearCat,AltCat) %>% 
    summarise(mean=mean(Dif_Roll_est,na.rm=T)) %>% filter(!is.na(AltCat)) %>%
    ggplot(aes(x=as.numeric(YearCat),y=(AltCat),fill=mean*0.01)) + 
    labs(x="",y="Elevation (meters)",title="Mean near-surface relative humidity") +
    geom_tile() +
    scale_fill_gradientn(colours= colorspace::darken(MetBrewer::met.brewer("Demuth",direction=-1)[c(3:5,6)],amount=0),
                         n.breaks=10,
                         guide=guide_legend(reverse = FALSE,nrow = 1,title.position="top",label.position="bottom", title.hjust = 0),name="") +
    scale_color_gradientn(colours= colorspace::darken(MetBrewer::met.brewer("Demuth",direction=-1)[c(2:5,6)],amount=0),n.breaks=10,guide=guide_legend(reverse = FALSE,nrow = 1,title.position="top",label.position="bottom", title.hjust = 0),name="") +
    theme(legend.position = "bottom",
          panel.background = element_blank(),
          axis.line = element_line(),
          legend.background = element_blank(),
          #legend.key.height = unit(0.01,"cm"),legend.key.width = unit(0.01,"cm"),
          axis.title = element_text(family="EB Garamond"),
          axis.text = element_text(family="EB Garamond"),
          axis.ticks.y=element_blank(),
          #axis.text.y=element_blank(),
          #axis.line.y=element_blank(),
          legend.spacing.y = unit(c(0.2),"cm"),
          legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond"),
          legend.margin=margin(t=-25),
          title = element_text(family="EB Garamond")) +
    NULL
  
  meanpet <- aqui %>% filter(!is.infinite(Roll_est)) %>% filter(Var=="mean.pet") %>% mutate(AltCat=cut(Altitude,breaks=seq(0,alt_limit+100,250),labels=seq(0,alt_limit+100,250)[-1], include.lowest=T)) %>% group_by(YearCat,AltCat) %>% 
    summarise(mean=mean(Dif_Roll_est,na.rm=T)) %>% filter(!is.na(AltCat)) %>%
    ggplot(aes(x=as.numeric(YearCat),y=(AltCat),fill=mean*0.01)) + 
    labs(x="",y="Elevation (meters)",title="Mean potential evapotranspiration") +
    geom_tile() +
    scale_fill_gradientn(colours= colorspace::darken(MetBrewer::met.brewer("Demuth",direction=-1)[c(5,6:9)],amount=0),
                         n.breaks=10,
                         guide=guide_legend(reverse = FALSE,nrow = 1,title.position="top",label.position="bottom", title.hjust = 0),name="") +
    scale_color_gradientn(colours= colorspace::darken(MetBrewer::met.brewer("Demuth",direction=-1)[c(5,6:9)],amount=0),n.breaks=10,guide=guide_legend(reverse = FALSE,nrow = 1,title.position="top",label.position="bottom", title.hjust = 0),name="") +
    theme(legend.position = "bottom",
          panel.background = element_blank(),
          axis.line = element_line(),
          legend.background = element_blank(),
          #legend.key.height = unit(0.01,"cm"),legend.key.width = unit(0.01,"cm"),
          axis.title = element_text(family="EB Garamond"),
          axis.text = element_text(family="EB Garamond"),
          axis.ticks.y=element_blank(),
          #axis.text.y=element_blank(),
          #axis.line.y=element_blank(),
          legend.spacing.y = unit(c(0.2),"cm"),
          legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond"),
          legend.margin=margin(t=-25),
          title = element_text(family="EB Garamond")) +
    NULL
  
  meanpr <- aqui %>% filter(!is.infinite(Roll_est)) %>% filter(Var=="mean.pr") %>% mutate(AltCat=cut(Altitude,breaks=seq(0,alt_limit+100,250),labels=seq(0,alt_limit+100,250)[-1], include.lowest=T)) %>% group_by(YearCat,AltCat) %>% 
    summarise(mean=mean(Dif_Roll_est,na.rm=T)) %>% filter(!is.na(AltCat)) %>%
    ggplot(aes(x=as.numeric(YearCat),y=(AltCat),fill=mean)) + 
    labs(x="",y="Elevation (meters)",title="Mean precipitation") +
    geom_tile() +
    scale_fill_gradientn(colours= colorspace::darken(MetBrewer::met.brewer("Demuth",direction=-1)[c(3:5,7:10)],amount=0),
                         n.breaks=10,
                         guide=guide_legend(reverse = FALSE,nrow = 1,title.position="top",label.position="bottom", title.hjust = 0),name="") +
    scale_color_gradientn(colours= colorspace::darken(MetBrewer::met.brewer("Demuth",direction=-1)[c(3:5,7:10)],amount=0),n.breaks=10,guide=guide_legend(reverse = FALSE,nrow = 1,title.position="top",label.position="bottom", title.hjust = 0),name="") +
    theme(legend.position = "bottom",
          panel.background = element_blank(),
          axis.line = element_line(),
          legend.background = element_blank(),
          #legend.key.height = unit(0.01,"cm"),legend.key.width = unit(0.01,"cm"),
          axis.title = element_text(family="EB Garamond"),
          axis.text = element_text(family="EB Garamond"),
          axis.ticks.y=element_blank(),
          #axis.text.y=element_blank(),
          #axis.line.y=element_blank(),
          legend.spacing.y = unit(c(0.2),"cm"),
          legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond"),
          legend.margin=margin(t=-25),
          title = element_text(family="EB Garamond")) +
    NULL
  
  minpr <- aqui %>% filter(!is.infinite(Roll_est)) %>% filter(Var=="min.pr") %>% mutate(AltCat=cut(Altitude,breaks=seq(0,alt_limit+100,250),labels=seq(0,alt_limit+100,250)[-1], include.lowest=T)) %>% group_by(YearCat,AltCat) %>% 
    summarise(mean=mean(Dif_Roll_est,na.rm=T)) %>% filter(!is.na(AltCat)) %>%
    ggplot(aes(x=as.numeric(YearCat),y=(AltCat),fill=mean)) + 
    labs(x="",y="Elevation (meters)",title="Minimum precipitation") +
    geom_tile() +
    scale_fill_gradientn(colours= colorspace::darken(MetBrewer::met.brewer("Demuth",direction=-1)[c(2:5,7:8)],amount=0),
                         n.breaks=10,
                         guide=guide_legend(reverse = FALSE,nrow = 1,title.position="top",label.position="bottom", title.hjust = 0),name="") +
    scale_color_gradientn(colours= colorspace::darken(MetBrewer::met.brewer("Demuth",direction=-1)[c(2:5,7:8)],amount=0),n.breaks=10,guide=guide_legend(reverse = FALSE,nrow = 1,title.position="top",label.position="bottom", title.hjust = 0),name="") +
    theme(legend.position = "bottom",
          panel.background = element_blank(),
          axis.line = element_line(),
          legend.background = element_blank(),
          #legend.key.height = unit(0.01,"cm"),legend.key.width = unit(0.01,"cm"),
          axis.title = element_text(family="EB Garamond"),
          axis.text = element_text(family="EB Garamond"),
          axis.ticks.y=element_blank(),
          #axis.text.y=element_blank(),
          #axis.line.y=element_blank(),
          legend.spacing.y = unit(c(0.2),"cm"),
          legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond"),
          legend.margin=margin(t=-25),
          title = element_text(family="EB Garamond")) +
    NULL
  
  
  layout <- "
AABBCC
DDEEFF
"
  
  pp<-wrap_plots(list(meantas,meantasmax,meanpr,minpr,meancmi,meanhurs)) +
    plot_layout(design = layout,byrow=T,nrow = 2,guides="keep",ncol=3,tag_level = 'new') + 
    plot_annotation(tag_levels = 'a',title="") & theme(plot.tag=element_text(family="EB Garamond",size=15,face="bold"),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))
  pp
  ggsave(filename = here("plots_tables/figure_3.pdf"),plot = pp,height = 7,width = 8.5)
  
}

plot_LUL_chords_deprecated <- function(data = "output/Historical_data_FULL_LUL.v.1.R",alt_limit = 4000,codes = c("PF","AU","NF","SV"),targetLU="PF") {
  load(here(data))
  reduct_five %>% select(CellID,Altitude,starts_with("LUC")) %>% distinct(CellID,.keep_all = T) %>% mutate(across(starts_with("LUC"), ~case_when(.x %in% c(1:6, 9, 12) ~"AU",.x %in% c(7)~"PF",.x %in% c(10,11)~"SV",.x %in% c(8)~"NF"))) %>% mutate(AltCat=cut(Altitude,breaks=seq(0,alt_limit+100,250),labels=seq(0,alt_limit+100,250)[-1], include.lowest=T),.before=Altitude) %>% group_by(LUC_1980,LUC_2015) %>% 
    summarise(Freq=n()) %>% filter(!is.na(LUC_1980)) %>% ungroup() %>%  filter(Freq  >= 0) %>% mutate(N=sum(Freq),Freq=Freq/sum(Freq)) %>% mutate_at(2,~paste(.x,"\nEnd",sep="")) %>% mutate_at(1,~paste(.x,"\nStart",sep="")) -> stats
  
  sub <- stats
  unique(c(sub[[1]],sub[[2]])) %>% length() -> n
  my_col <- rep(colorspace::darken(MetBrewer::met.brewer("Archambault",n=8)[c(2,8,8,8)],amount = 0.2),each=2)
  unique(c(sub[[1]],sub[[2]])) -> names
  lapply(codes,function(x) {grep(x,names)}) %>% unlist -> index
  names(my_col) <- names[index]
  my_col[!grep(targetLU,names(my_col))] <- colorspace::lighten(my_col[!grep(targetLU,names(my_col))],amount = 0.4)
  alphas <- rep(0.6,length(sub[[1]])); 
  alphas[grep(targetLU,sub[[1]])] <- 0.1
  circlize::chordDiagram(x = sub %>% select(-N), 
                         grid.col = my_col,transparency = alphas,directional = 1,
                         direction.type = c("arrows"), 
                         link.arr.type = "big.arrow",
                         link.arr.length = 0.1,link.sort = T,
                         link.largest.ontop = TRUE,grid.border="black",
                         link.border=adjustcolor("black",0.5),
                         link.lwd=0.2#,order=names(my_col)
  )
  
  
  
}

plot_lulseries_deprecated <- function(lul = "Historical_dataLULC.R",alt_breaks = c(0,300,1000,2000,4000),cat = c("Low","MidLow","MidHigh","High")){
  load(here(paste("output_data",lul,sep="/")))
  my_cols <- viridis(4,option="G",direction=1,begin = 0.2,end=0.8)
  lucs_tib %>% filter(decimallatitude > lat_limits[1] & decimallongitude < lon_limits[1],
                      !(decimallatitude > lat_limits[2] | decimallongitude < lon_limits[2]),
                      !(decimallatitude > lat_limits[3] & decimallongitude > lon_limits[3]),
                      !CellID%in%all_of(non_id),AltRas<alt_limit) %>% dplyr::select(CellID,x,y,AltRas,starts_with("LUC")) %>% distinct(CellID,.keep_all = T) %>% droplevels() %>% arrange(CellID) -> precs
  precs %>% pull(AltRas) %>% cut(.,breaks=seq(0,alt_limit+100,belts),labels=seq(0,alt_limit+100,belts)[-1], include.lowest=T) -> alt_cats ### aggregate records into altitudinal belts
  precs <- precs %>% ungroup() %>% mutate(AltCats = alt_cats)
  pp <- precs %>% mutate(Elevation=cut(AltRas,breaks=alt_breaks,labels=cat),.before=AltRas) %>% mutate(across(starts_with("LUC"), ~case_when(.x %in% c(1:6, 9, 12) ~"AU",.x %in% c(7)~"PF",.x %in% c(10,11)~"SV",.x %in% c(8)~"NF"))) %>% select(-AltRas,-x,-y) %>% pivot_longer(cols = starts_with("LUC"),names_to = "Year",values_to = "LandUse") %>% mutate(Year = as.numeric(sub("LUC_","",Year))) %>% group_by(Year,LandUse,Elevation) %>% summarize(Freq=n()) %>% filter(!is.na(LandUse)) %>%
    mutate(Category=paste(LandUse,Elevation,sep="_")) %>% 
    ggplot(aes(x=Year, y=Freq, fill=LandUse)) + 
    geom_area() + scale_fill_viridis_d(option="G",direction=-1,begin = 0.1,end=0.7) +
    theme(legend.position = "bottom",panel.background = element_blank(),panel.grid = element_blank(),axis.ticks = element_blank(),axis.text.x = element_text(size=10),axis.line = element_line()) + labs(y="Frequency",title="Land use in cloud forests through time") + geom_vline(xintercept = 1976) +
    facet_wrap(~ Category,strip.position = "top") +
    NULL
  ggsave(filename = here("plots/LUL_TimeSeries.pdf"),plot = pp)
}

plot_granger_deprecated <- function(data="output_data/GrangerFinal_Min.Rdata",output="plots/Granger_Min.pdf",TRY="output_data/growth_try.Rdata",  data_lul="output_data/Historical_dataLULC.R") {
  load(here(data))
  load(here(data_lul))
  load(here(TRY))
  load(here("output_data/alt_species.Rdata"))
  lucs_tib %>% distinct(Resolved_ACCEPTED) -> species_list
  growth_try %>% filter(Resolved_ACCEPTED %in% all_of(species_list$Resolved_ACCEPTED)) %>% group_by(Resolved_ACCEPTED) %>% nest() -> target_try
  target_try %>% rowwise() %>% mutate(LifeForm=(Mode(data$OrigValueStr))) -> target_try
  left_join(species_list,target_try,"Resolved_ACCEPTED") %>% select(-data) %>% mutate(LifeForm=replace_na(LifeForm,"Unknown")) %>% mutate(LifeForm=factor(LifeForm)) -> species_list
  final %>% unnest(gci_luc) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","gci_luc","prob_luc")) %>% inner_join(.,{final %>% unnest(gci_temps) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","gci_temps","prob_temps"))},by="Resolved_ACCEPTED") %>% 
    inner_join(.,{final %>% unnest(gci_precs) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","gci_precs","prob_precs"))},by="Resolved_ACCEPTED") %>% 
    inner_join(.,{final %>% unnest(gci_climate) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","gci_climate","prob_climate"))},by="Resolved_ACCEPTED") %>% 
    inner_join(.,{final %>% unnest(pgci_climate) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","pgci_climate","prob_p_climate"))},by="Resolved_ACCEPTED") %>% 
    inner_join(.,{final %>% unnest(gci_all) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","gci_all","prob_all"))},by="Resolved_ACCEPTED") %>% 
    inner_join(.,species_list,by="Resolved_ACCEPTED") %>% inner_join(.,alt_species,by="Resolved_ACCEPTED") %>% ungroup() %>% pivot_longer(cols = contains("gci"),names_to = "Series",values_to = "gci") %>% mutate(AltCat=fct_relevel(AltCat,"Low","Mid Low","Mid High","High")) %>% pivot_wider(names_from = "Series",values_from = "gci") %>% mutate(gci_luc = case_when(prob_luc > 0.05 ~ as.numeric(NA),prob_luc <= 0.05 ~ gci_luc),
                                                                                                                                                                                                                                                                                                                                                       gci_temps = case_when(prob_temps > 0.05 ~ as.numeric(NA),prob_temps <= 0.05 ~ gci_temps),
                                                                                                                                                                                                                                                                                                                                                       gci_precs = case_when(prob_precs > 0.05 ~ as.numeric(NA),prob_precs <= 0.05 ~ gci_precs),
                                                                                                                                                                                                                                                                                                                                                       gci_climate = case_when(prob_climate > 0.05 ~ as.numeric(NA),prob_climate <= 0.05 ~ gci_climate),
                                                                                                                                                                                                                                                                                                                                                       pgci_climate = case_when(prob_p_climate > 0.05 ~ as.numeric(NA),prob_p_climate <= 0.05 ~ pgci_climate),
                                                                                                                                                                                                                                                                                                                                                       gci_all = case_when(prob_all > 0.05 ~ as.numeric(NA),prob_all <= 0.05 ~ gci_all)
    ) %>% pivot_longer(contains("gci"),names_to = "Series",values_to = "gci") -> out_plot
  out_plot %>% group_by(Series) %>% filter(!is.na(gci)) %>% summarise(n=n()) %>% mutate(names=paste(Series," (",n,")",sep="")) %>% pull(names) -> names_series
  pp <- out_plot %>% #mutate(Series=fct_relevel(Series,"gci_all","gci_luc","gci_precs","gci_temps","gci_climate","pgci_climate")) %>%  
    ggplot(aes(x=LifeForm,y=gci,fill=Series)) + 
    geom_boxplot(outlier.size=0.3,position="dodge2",size=0.2) +
    labs(title="Causal influence of historical climate change and vegetation loss on species' elevational ranges",x="Lifeform",y="Granger Causality Index") + 
    theme(panel.background = element_blank(),panel.grid = element_blank(),axis.line = element_line(),legend.position = "right",axis.text.x = element_text(size = 10),legend.text = element_text()) + scale_fill_viridis_d(option="G",begin=0.4,name="Series",labels=names_series) + 
    ylim(c(0,4)) +
    NULL
  ggsave(filename = here(output),plot = pp)
}

plot_granger_alternative_deprecated <- function(data="output_data/GrangerFinal_Min.Rdata",output="plots/Granger_Mid_Alternative.pdf",TRY="output_data/growth_try.Rdata",  data_lul="output_data/Historical_dataLULC.R") {
  load(here(data))
  load(here(data_lul))
  load(here(TRY))
  load(here("output_data/alt_species.Rdata"))
  lucs_tib %>% distinct(Resolved_ACCEPTED) -> species_list
  growth_try %>% filter(Resolved_ACCEPTED %in% all_of(species_list$Resolved_ACCEPTED)) %>% group_by(Resolved_ACCEPTED) %>% nest() -> target_try
  target_try %>% rowwise() %>% mutate(LifeForm=(Mode(data$OrigValueStr))) -> target_try
  left_join(species_list,target_try,"Resolved_ACCEPTED") %>% select(-data) %>% mutate(LifeForm=replace_na(LifeForm,"Unknown")) %>% mutate(LifeForm=factor(LifeForm)) -> species_list
  final %>% unnest(gci_luc) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","gci_luc","prob_luc")) %>% inner_join(.,{final %>% unnest(gci_temps) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","gci_temps","prob_temps"))},by="Resolved_ACCEPTED") %>% 
    inner_join(.,{final %>% unnest(gci_precs) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","gci_precs","prob_precs"))},by="Resolved_ACCEPTED") %>% 
    inner_join(.,{final %>% unnest(gci_climate) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","gci_climate","prob_climate"))},by="Resolved_ACCEPTED") %>% 
    inner_join(.,{final %>% unnest(pgci_climate) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","pgci_climate","prob_p_climate"))},by="Resolved_ACCEPTED") %>% 
    inner_join(.,{final %>% unnest(gci_all) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","gci_all","prob_all"))},by="Resolved_ACCEPTED") %>% 
    inner_join(.,species_list,by="Resolved_ACCEPTED") %>% inner_join(.,alt_species,by="Resolved_ACCEPTED") %>% ungroup() %>% pivot_longer(cols = contains("gci"),names_to = "Series",values_to = "gci") %>% mutate(AltCat=fct_relevel(AltCat,"Low","Mid Low","Mid High","High")) %>% pivot_wider(names_from = "Series",values_from = "gci") %>% mutate(gci_luc = case_when(prob_luc > 0.05 ~ as.numeric(NA),prob_luc <= 0.05 ~ gci_luc),
                                                                                                                                                                                                                                                                                                                                                       gci_temps = case_when(prob_temps > 0.05 ~ as.numeric(NA),prob_temps <= 0.05 ~ gci_temps),
                                                                                                                                                                                                                                                                                                                                                       gci_precs = case_when(prob_precs > 0.05 ~ as.numeric(NA),prob_precs <= 0.05 ~ gci_precs),
                                                                                                                                                                                                                                                                                                                                                       gci_climate = case_when(prob_climate > 0.05 ~ as.numeric(NA),prob_climate <= 0.05 ~ gci_climate),
                                                                                                                                                                                                                                                                                                                                                       pgci_climate = case_when(prob_p_climate > 0.05 ~ as.numeric(NA),prob_p_climate <= 0.05 ~ pgci_climate),
                                                                                                                                                                                                                                                                                                                                                       gci_all = case_when(prob_all > 0.05 ~ as.numeric(NA),prob_all <= 0.05 ~ gci_all)
    ) %>% pivot_longer(contains("gci"),names_to = "Series",values_to = "gci") -> out_plot
  counts <- out_plot %>% filter(!is.na(gci)) %>% count(Series)
  meanas <- out_plot %>% filter(!is.na(gci)) %>% group_by(Series) %>% summarise(n=median(gci))
  names_series <- meanas %>% mutate(Names=str_to_sentence(sub("gci_","",Series)))
  names_series$Names[6] <- "Climate|LUC"
  names_series$Names[3] <- "LUC"
  pp <- out_plot %>% filter(!is.na(gci)) %>% 
    ggplot(aes(x = Series, y = gci,fill=Series)) + 
    ggdist::stat_halfeye(
      adjust = .9, 
      width = .6, 
      .width = 0, 
      alpha=0.7,
      slab_colour = "black",
      slab_size = .2,
      justification = -.4, 
      height = 2) + 
    geom_boxplot(width = .4, outlier.shape = NA,size=.2,fatten = 7,alpha=0.7) +
    geom_point(fill="grey50",size = 2,alpha = .4,shape = 21,stroke = 0.2,position = position_jitter(seed = 1, width = .2)) + 
    coord_flip() +
    scale_fill_viridis_d(option="G",begin=0.4,name="Series",labels=names_series)  + 
    scale_y_log10() + theme(panel.background = element_blank(),panel.grid = element_blank(),axis.line = element_line(),legend.position = "",axis.text.x = element_text(size = 10),legend.text = element_text()) + annotate(
      geom = "text", x = (1:6)+.2, y = 15, 
      label = paste0("n = ",counts$n), hjust = 0, vjust = -.8, size = 2,fontface = "bold",
    ) + annotate(
      geom = "text", x = (1:6)+.2, y = meanas$n, fontface = "bold",
      label = round(meanas$n,2), hjust = .5, vjust = -.8, size = 3
    ) + scale_x_discrete(labels= names_series$Names) + labs(x="",y="Granger Causality Index") +
    NULL
  ggsave(filename = here(output),plot = pp)
}

lmm_fit_deprecated <- function(data="output_data/Splines_v.4.Rdata",cual="Min",flex_point = 13,out_table="output_data/LMM_modelsMin.Rdata"){
  load(here(data))
  mean_alts %>% mutate(Type="OBS",.before=data) %>% rowwise(Resolved_ACCEPTED) %>% select(Resolved_ACCEPTED,Type,data) %>% unnest(cols = c(data)) %>% mutate(Quarter=case_when(YearCat < 1921 ~ "1",YearCat < 1926 & YearCat >= 1921 ~ "2",
                                                                                                                                                                               YearCat < 1931 & YearCat >= 1926 ~ "3",YearCat < 1936 & YearCat >= 1931 ~ "4",
                                                                                                                                                                               YearCat < 1941 & YearCat >= 1936 ~ "5",YearCat < 1946 & YearCat >= 1941 ~ "6",
                                                                                                                                                                               YearCat < 1951 & YearCat >= 1946 ~ "7",YearCat < 1956 & YearCat >= 1951 ~ "8",
                                                                                                                                                                               YearCat < 1961 & YearCat >= 1956 ~ "9",YearCat < 1966 & YearCat >= 1961 ~ "10",
                                                                                                                                                                               YearCat < 1971 & YearCat >= 1966 ~ "11",YearCat < 1976 & YearCat >= 1971 ~ "12",
                                                                                                                                                                               YearCat < 1981 & YearCat >= 1976 ~ "13",YearCat < 1986 & YearCat >= 1981 ~ "14",
                                                                                                                                                                               YearCat < 1991 & YearCat >= 1986 ~ "15",YearCat < 1996 & YearCat >= 1991 ~ "16",
                                                                                                                                                                               YearCat < 2001 & YearCat >= 1996 ~ "17",YearCat < 2006 & YearCat >= 2001 ~ "18",
                                                                                                                                                                               YearCat < 2011 & YearCat >= 2006 ~ "19",YearCat < 2016 & YearCat >= 2011 ~ "20")) %>% 
    filter(Type=="OBS") %>% group_by(Resolved_ACCEPTED,Quarter) %>% 
    summarize(across(c(RollAlt,RollAlt_Max,RollAlt_Min),~median(.x,na.rm=T))) -> usl
  usl %>% rename_with(~c("Resolved_ACCEPTED","Quarter","Mid","Max","Min")) %>% filter(!is.na(Quarter)) %>% pivot_wider(names_from = c(Quarter),values_from = c(Mid,Max,Min)) -> usl
  apply(usl[,-1],2,function(x) length(which(is.na(x))))
  usl %>% ungroup() %>% mutate_if(is.numeric,log) -> usl_log
  usl %>% ungroup() %>% pivot_longer(cols = 2:61,names_to = "Wave",values_to = "Alt") %>% separate(Wave,into = c("Var","Wave"),sep = "_") %>% mutate(Wave=as.numeric(Wave)) %>% mutate(wave0=Wave -1 ) %>% mutate(wave1 = wave0 - flex_point) %>% mutate(wave1=case_when(wave1 < 0 ~ 0,wave1>=0~wave1)) %>% mutate(wave01 = case_when(wave0 > flex_point ~ flex_point, wave0 <= flex_point ~ wave0))  %>% filter(Var==cual) -> usl_2
  lmer(data = usl_2, Alt ~ 1 + (1 | Resolved_ACCEPTED),control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4))) -> fit0
  lmer(data = usl_2, Alt ~ 1 + wave01 + (1 | Resolved_ACCEPTED),control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4))) -> fit1
  lmer(data = usl_2, Alt ~ 1 + wave01 + (1|Resolved_ACCEPTED) + (0 + wave01 | Resolved_ACCEPTED),control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4))) -> fit2
  lmer(data = usl_2, Alt ~  1 + wave01 + wave1 + (1|Resolved_ACCEPTED) + (0 + wave01| Resolved_ACCEPTED) +(0 + wave1| Resolved_ACCEPTED), control=lmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=1e4))) -> fit3
  models <- list(fit0,fit1,fit2,fit3)
  save(models,file = here(out_table))
}

plot_SppSeries <- function(data="output/Species_summaryRaw.RData",inflect_point = 1990,output="plots/Species_Series.TRY.pdf"){ 
  load(here(data))
  names(species_means) %>% grepl("trim",.) %>% !. -> cuales
  quien <- which(cuales)
  direct = c(-1,1,-1,-1)
  for (k in 1:length(quien)){
    cat("Plotting",k,"\r")
    pp <- species_means[quien][[k]] %>% group_by(Accepted_Name) %>% nest() %>% ungroup() %>% mutate(data = dif_clim(data,index = c(5,7),flex_pt = inflect_point)) %>% unnest(data) %>% #rename(Observed_baseline = Observed, Trend_baseline = Trend,YearCat2 = YearCat) %>% 
      inner_join(.,species_means[quien][[k]] %>% select(YearCat,Accepted_Name,LifeForm,SES_resamp_sig,SES_down_sig,SES_resamp,SES_down),by=c("Accepted_Name","YearCat")) %>% 
      #rename(LifeForm2  = LifeForm) %>% unnest(data,names_repair = "unique") %>% 
      group_by(YearCat,LifeForm) #%>%
    
    pp %>% group_by(YearCat,SES_1 = mean(abs(SES_resamp),na.rm=T),SES_2 = mean(abs(SES_down),na.rm=T)) %>% ungroup() %>% select(SES_2,SES_1) %>% pivot_longer(1:2) %>% summarise(range(value)) %>% pull(1) -> limites
    
    uno <- pp %>% group_by(YearCat) %>% #filter(SES_resamp_sig) %>%
      #mutate(SES_resamp=normalize(SES_resamp)) %>% 
      summarise(Observed_baseline=mean(Observed,na.rm=T),SES = mean(abs(SES_resamp),na.rm=T)) %>% 
      ggplot(aes(x=as.numeric(YearCat), y=Observed_baseline, fill= SES)) +
      geom_bar(stat="identity") + 
      labs(title = "Tendencias de elevación en plantas de bosque de niebla",
           subtitle = "Comparación con modelo nulo de remuestro", 
           y="Metros relativo a referencia",x="") + 
      theme(panel.background = element_blank(),axis.line = element_line(),axis.title = element_text(family="EB Garamond"),plot.title = element_text(family="EB Garamond",size=14,hjust = 0.5,face="bold"),
            plot.subtitle = element_text(family="EB Garamond",size=10,hjust = 0.5),axis.text = element_text(family="EB Garamond"),legend.text=element_text(family="EB Garamond"),
            legend.title = element_text(family="EB Garamond"),legend.position = "none",legend.key.size = unit(0.3, "cm"),legend.key.width = unit(1.5,"cm")) +
      scale_fill_stepsn(n.breaks=11,colors=MetBrewer::met.brewer("Hiroshige",direction=direct[k]),limits=limites) + 
      geom_vline(linetype="dashed",xintercept = inflect_point) + #facet_wrap(~LifeForm) +
      standardPrintOutput::watermark(show=T,lab="Ramirez-Barahona et al.") +
      NULL
    
    dos <- pp %>% group_by(YearCat) %>% #filter(SES_down_sig) %>%
      #mutate(SES_down=normalize(SES_down)) %>% 
      summarise(Observed_baseline=mean(Observed,na.rm=T),SES=mean(abs(SES_down),na.rm=T)) %>% 
      ggplot(aes(x=as.numeric(YearCat), y=Observed_baseline, fill= SES)) +
      geom_bar(stat="identity") + 
      labs(
        subtitle = "Comparación con modelo nulo de submuestreo", 
        y="Metros relativo a referencia",x="") + 
      theme(panel.background = element_blank(),axis.line = element_line(),axis.title = element_text(family="EB Garamond"),plot.title = element_text(family="EB Garamond",size=14,hjust = 0.5,face="bold"),
            plot.subtitle = element_text(family="EB Garamond",size=10,hjust = 0.5),axis.text = element_text(family="EB Garamond"),legend.text=element_text(family="EB Garamond"),
            legend.title = element_text(family="EB Garamond"),legend.position = "bottom",legend.key.size = unit(0.3, "cm"),legend.key.width = unit(1.5,"cm")) +
      scale_fill_stepsn(n.breaks=11,colors=MetBrewer::met.brewer("Hiroshige",direction=direct[k]),limits=limites) + 
      geom_vline(linetype="dashed",xintercept = inflect_point) + #facet_wrap(~LifeForm) +
      standardPrintOutput::watermark(show=T,lab="Ramirez-Barahona et al.") +
      NULL
    
    myplots <- "
AA
BB
"
    estos <- wrap_plots(list(uno,dos))  +
      plot_layout(design=myplots,guides="collect",tag_level = 'new') & 
      theme(legend.position = 'bottom')
    estos
    
    
    nombre <- names(species_means)[k]
    output2 <- output %>% sub(".pdf","",.) %>% paste0(.,"_",nombre,".pdf")
    ggsave(filename = here(output2),plot = estos)
  }
}

plot_raw_clim <- function (data="output/Historical_data_FULL_LUL.v.1.R",prec=F,variables="_cmi_"){
  load(here(data))
  
  reduct_five %>% select(-c(Accepted_Name:scientificName,year:BIOME)) %>% select(CellID:decimalLongitude,contains(all_of(variables))) %>% distinct(CellID,.keep_all = T) %>% droplevels() %>% pivot_longer(cols = -c(1:4),names_to = "year",values_to = "Variable") %>% 
    separate(year,into = c("summary","year"),sep=all_of(variables)) -> esta
  esta %>% filter(summary=="mean") %>% arrange(CellID,year) %>% mutate(Roll_est = unlist(tapply(Variable,CellID,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean, na.rm=T, partial=T, align="center")}),use.names = F)) %>% group_by(CellID) -> actual
  actual %>% pull(Altitude) %>% cut(.,breaks=seq(0,alt_limit+100,500),labels=seq(0,alt_limit+100,500)[-1], include.lowest=T) -> alt_cats ### aggregate records into altitudinal belts
  actual <- actual %>% ungroup() %>% mutate(AltCats = alt_cats)
  # Estimate means across years.
  actual$AltCats %>% unique() %>% as.character() %>% as.numeric() %>% .[!is.na(.)] %>% sort() -> catas
  catas_2 <- catas - 499
  #####
  #### PLOT ####
  plotas <- list()
  for (i in 1:length(catas)){
    plotas[[i]] <- actual %>% group_by(year,AltCats) %>% summarise(Roll_mean=mean(Roll_est,na.rm=T),seme=seme(Roll_est,na.rm=T)) %>% filter(AltCats==catas[i]) %>% ggplot(aes(x=year,y=Roll_mean,group=AltCats,color=Roll_mean)) + 
      geom_errorbar(aes(ymin = Roll_mean - seme,ymax = Roll_mean +  seme)) +
      geom_point() +
      viridis::scale_color_viridis(discrete=F,option="G",direction=-1) +
      
      theme(panel.background = element_blank(),legend.position = "none") + labs(y=code,x="",subtitle=paste(catas_2[i],"\u2013",catas[i]," masl",sep="")) + NULL
  }
  
  wrap_plots(plotas) +
    plot_layout(byrow=T,nrow = 3,ncol=4,guides = "collect",tag_level = 'new') + 
    plot_annotation(title = nombres,subtitle="1901\u20132016",
                    caption= "Ramírez-Barahona et al.") & theme(legend.position = '')
}

plot_Chords_deprecated <- function(lul = "../2.BMM/CF_Project/output_data/Historical_dataLULC.R",alt_breaks = c(0,300,1000,2000,4000),before = c("LUC_1936","LUC_1956","LUC_1976","LUC_1996"), after = c("LUC_2015","LUC_2015","LUC_2015","LUC_2015"),cat = c("Low","MidLow","MidHigh","High"),codes = c("PF","AU","NF","SV"),targetLU="PF",by_alt=T) {
  load(here(lul))
  lucs_tib %>% dplyr::select(CellID,x,y,AltRas,starts_with("LUC")) %>% distinct(CellID,.keep_all = T) %>% droplevels() %>% arrange(CellID) -> precs
  precs %>% pull(AltRas) %>% cut(.,breaks=seq(0,alt_limit+100,belts),labels=seq(0,alt_limit+100,belts)[-1], include.lowest=T) -> alt_cats ### aggregate records into altitudinal belts
  precs <- precs %>% ungroup() %>% mutate(AltCats = alt_cats)
  # Estimate means across years.
  for (j in 1:4){
    if(by_alt == T) { 
      pdf(here(paste("plots/LUL_transition_",paste(cat[i],"_",sub("LUC_","",before[j]),"-",sub("LUC_","",after[j]),sep=""),".pdf",sep="")))}
    if(by_alt == F) { 
      pdf(here(paste("plots/LUL_transition_",paste("Global","_",sub("LUC_","",before[j]),"-",sub("LUC_","",after[j]),sep=""),".pdf",sep="")))}
    if(by_alt == T){ 
      for (i in 1:4){
        precs %>% mutate(Elevation=cut(AltRas,breaks=alt_breaks,labels=cat),.before=AltRas) %>% mutate(across(starts_with("LUC"), ~case_when(.x %in% c(1:6, 9, 12) ~"AU",.x %in% c(7)~"PF",.x %in% c(10,11)~"SV",.x %in% c(8)~"NF"))) %>% group_by(across(all_of(c(before[j],after[j]))),Elevation) %>% summarise(Freq=n()) %>% ungroup() %>%  filter(across(1,~ !is.na(.))) %>% 
          filter(across(Freq,~ . >= 0)) %>% mutate(Elevation=as.factor(Elevation)) %>% 
          group_by(Elevation) %>%  mutate(N=sum(Freq),Freq=Freq/sum(Freq)) %>% mutate_at(2,~paste(.x,"\nEnd",sep="")) %>% mutate_at(1,~paste(.x,"\nStart",sep="")) -> stats
        stats %>% filter(Elevation == all_of(cat[i])) %>% ungroup() -> stats
        unique(c(stats[[1]],stats[[2]])) %>% length() -> n
        my_col <- viridis(n=n,option="G",direction=-1,end=0.7)
        unique(c(stats[[1]],stats[[2]])) -> names
        lapply(codes,function(x) {grep(x,names)}) %>% unlist -> index
        names(my_col) <- names[index]
        alphas <- rep(0.8,length(stats[[1]])); alphas[grep(targetLU,stats[[1]])] <- 0.2
        circlize::chordDiagram(x = stats %>% select(-N), 
                               grid.col = my_col,transparency = alphas,directional = 1,
                               direction.type = c("arrows"), link.arr.type = "big.arrow", link.arr.length = 0.1,
                               link.sort = T,link.largest.ontop = TRUE,grid.border="black",
                               link.border=adjustcolor("black",1),link.lwd=0.2,order=names(my_col))
        title(main=paste(cat[i]," (",sub("LUC_","",before[j]),"-",sub("LUC_","",after[j]),")",sep=""))
      }}
    if(by_alt == F){
      precs %>% mutate(Elevation=cut(AltRas,breaks=alt_breaks,labels=cat),.before=AltRas) %>% mutate(across(starts_with("LUC"), ~case_when(.x %in% c(1:6, 9, 12) ~"AU",.x %in% c(7)~"PF",.x %in% c(10,11)~"SV",.x %in% c(8)~"NF"))) %>% group_by(across(all_of(c(before[j],after[j])))) %>% summarise(Freq=n()) %>% ungroup() %>%  filter(across(1,~ !is.na(.))) %>% 
        filter(across(Freq,~ . >= 0)) %>% mutate(N=sum(Freq),Freq=Freq/sum(Freq)) %>% mutate_at(2,~paste(.x,"\nEnd",sep="")) %>% mutate_at(1,~paste(.x,"\nStart",sep="")) -> stats
      unique(c(stats[[1]],stats[[2]])) %>% length() -> n
      my_col <- viridis(n=n,option="G",direction=-1,end=0.7)
      unique(c(stats[[1]],stats[[2]])) -> names
      lapply(codes,function(x) {grep(x,names)}) %>% unlist -> index
      names(my_col) <- names[index]
      alphas <- rep(0.8,length(stats[[1]])); alphas[grep(targetLU,stats[[1]])] <- 0.2
      circlize::chordDiagram(x = stats %>% select(-N), 
                             grid.col = my_col,transparency = alphas,directional = 1,
                             direction.type = c("arrows"), link.arr.type = "big.arrow", link.arr.length = 0.1,
                             link.sort = T,link.largest.ontop = TRUE,grid.border="black",
                             link.border=adjustcolor("black",1),link.lwd=0.2,order=names(my_col))
      title(main=paste("Global"," (",sub("LUC_","",before[j]),"-",sub("LUC_","",after[j]),")",sep=""))
    }
    dev.off()
  }
}

plot_map_deprecated <- function(data="Historical_dataClimate_full.R",output_name="RichnessMap.pdf",palette = "Greens",title="Species richness") {
  load(here(paste("output_data/",data,sep="")))
  pp <- reduct_five %>% filter(decimallatitude > 12.5 & decimallongitude < -77,
                               !(decimallatitude > 26 | decimallongitude < -109),
                               !(decimallatitude > 17 & decimallongitude > -86),
                               !CellID%in%all_of(non_id),AltRas<alt_limit) %>% dplyr::select(CellID,x,y,AltRas) %>% group_by(CellID) %>% summarise(SR=n(),x=first(x),y=first(y),AltRas=mean(AltRas)) %>% ggplot() + geom_sf(data=ne_countries(scale=110,returnclass = "sf",type="countries"),colour="grey95",fill="grey95",size=0.5) + xlim(xlimits) + ylim(ylimits) + 
    geom_tile(aes(x=x,y=y,fill=log(SR)),color="black",size=0.05) + 
    scale_fill_fermenter(type = "seq",palette = palette,direction=1,n.breaks=9) + 
    theme(panel.background = element_blank(),panel.grid = element_blank(),axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.border = element_rect(fill = NA)) + labs(title=title)+
    NULL
  ggsave(filename=here(paste("plots/",output_name,sep="")),plot(pp))
}

### POSSIBLE DUPLICATES ####

create_weight_deprecated <- function(x, y, shape = 3, trend = 1, distortion = FALSE) {
  get <- x - start_year - 1
  # generate a random beta distribution
  if(distortion) {dist = 0.01 * get; dens <- rbeta(10000,shape, shape - dist)}
  if(!distortion) dens <- rbeta(10000,shape, shape)
  
  # transform beta distribution into elevation
  rr <- range_focal[,-1] + (trend * get)
  dens <- rr$Lower + (dens * (rr$Upper - rr$Lower))
  summary(dens)
  # get probability distribution and estimate probability for each elevation class
  dens_beta <- ecdf(dens)(seq(0,5000,100)) - lag(ecdf(dens)(seq(0,5000,100)),1)
  dens_beta[1] <- ecdf(dens)(seq(0,5000,100))[1]
  
  
  # Generate a random distribution with mean and sd from empirical baseline (pre 1985)
  dens <- rnorm(10000, mean = dist_focal$mean + (trend * get), sd = dist_focal$sd)
  dens <- dens[dens>=0]
  # get probability distribution and estimate probability for each elevation class
  dens_2 <- ecdf(dens)(seq(0,5000,100)) - lag(ecdf(dens)(seq(0,5000,100)),1)
  dens_2[1] <- ecdf(dens)(seq(0,5000,100))[1]
  points(seq(0,5000,100),dens_2,col="blue")
  
  
  # multiply beta distribution by normal distribution and standardize to sum = 1.
  dens_2 <- dens_2*dens_beta / sum(dens_2*dens_beta)
  points(seq(0,5000,100),dens_2,col="red")
  
  
  # Set sampling weights by defining elevation classes.
  y %>% mutate(Weight = cut(Altitude,breaks=seq(0,5000,50),labels = dens_2[-length(dens_2)]),.before=Altitude) %>% mutate(Weight = as.numeric(as.character(Weight)))
} 
LMMs_deprecated <- function(data_lul="output/Historical_data_FULL_GLASS.v.1.R",
                            output="LULSeriesGLASS_v.2.RData",
                            splines="output/Splines_v.3.Rdata",
                            start_year=1979,
                            end_year=2010,
                            pattern="Max"){
  
  load(here(data_lul))
  load(here(splines))
  aver <- mean_alts %>% mutate(Type="OBS",.before=data) %>% select(2:3) %>% unnest(data)
  
  if(is.null(pattern)) {
    pat = "RollAlt"
    pat_sd = paste0(pat,"_sd")
  }
  if(!is.null(pattern)) { 
    pat = paste0("RollAlt_",pattern)
    pat_sd = paste0(pat,"_sd")
  }
  aver %<>% ungroup() %>% 
    select(Accepted_Name,YearCat,ends_with(pat),Rolln,starts_with(pat_sd)) %>% 
    mutate(YearCat=as.numeric(YearCat)) %>% 
    mutate(RollAlt_error = unlist(across(starts_with(pat_sd), ~ .x / sqrt(Rolln)))) %>% 
    select(ends_with(pat),YearCat,Accepted_Name,RollAlt_error) %>% rename_with(~c("y","tme","ind","error"))
  
  
  
  
}
plot_ClimSeries_deprecated <- function(data="output_data/ClimateRAW_by_AltCat.R",inflect_point = 1976,output="plots/ClimateSeries_PrecRAW.pdf",variable = "pre",viridis_option = "G",ind=c(2:9)) {
  load(here(data))
  aa %>% mutate(Shift=case_when(YearCat < inflect_point ~ "Pre",YearCat >= inflect_point ~ "Post")) %>% group_by(Shift,AltCats) %>% summarise(across(c(2:9),mean)) -> means
  plotas <- list()
  colNames <- names(aa)[grep(variable,names(aa))]
  for (i in 1:length(colNames)){
    if(i %in% c(1,3)){
      plotas[[i]] <- aa %>% nest() %>% mutate(data = dif_clim(data,flex_pt = inflect_point,index=ind)) %>% unnest(cols = data) %>% 
        ggplot(aes_string(x=aa$YearCat, y=aa$AltCats, fill= colNames[i])) + 
        geom_tile() + labs(subtitle = colNames[i], y="Altitude",x="") + 
        theme(panel.background = element_blank(),axis.ticks.y = element_blank(),legend.position = "right",legend.key.height = unit(0.8,"cm"),legend.key.width = unit(0.2,"cm"),legend.text = element_text(size=6),legend.margin = margin(0,0,0,0)) + ylim(0,alt_limit+200) +
        scale_fill_viridis_b(n.breaks=20,option=viridis_option,direction= 1,name="")
    }
    if(i %in% c(2,4)){
      plotas[[i]] <- aa %>% nest() %>%  mutate(data = dif_clim(data,flex_pt = inflect_point,index=ind)) %>% unnest(cols = data) %>% 
        ggplot(aes_string(x=aa$YearCat, y=aa$AltCats, fill= colNames[i])) + 
        geom_tile() + labs(subtitle = colNames[i], y="Altitude",x="") + 
        theme(panel.background = element_blank(),axis.ticks.y = element_blank(),legend.position = "right",legend.key.height = unit(0.8,"cm"),legend.key.width = unit(0.2,"cm"),axis.text.y =  element_blank(),axis.title.y = element_blank(),legend.text = element_text(size=6),legend.margin = margin(0,0,0,0)) + ylim(0,alt_limit+200) +
        scale_fill_viridis_b(n.breaks=20,option=viridis_option,direction= 1,name="")
    }
  }
  myplots <- "
AABB
CCDD
"
  pp <- wrap_plots(plotas) +
    plot_layout(design=myplots,guides = "keep",tag_level = 'new') + 
    plot_annotation(title = "Historical climate trends across cloud forests",
                    subtitle="Deviations from baseline trend (per altitudinal belt)", 
                    caption= "Ramírez-Barahona et al.") & theme(legend.position = 'right')
  ggsave(filename = here(output),plot = pp)
}

plot_lulseries_deprecated <- function(lul = "Historical_dataLULC.R",alt_breaks = c(0,300,1000,2000,4000),cat = c("Low","MidLow","MidHigh","High")){
  load(here(paste("output_data",lul,sep="/")))
  my_cols <- viridis(4,option="G",direction=1,begin = 0.2,end=0.8)
  lucs_tib %>% filter(decimallatitude > lat_limits[1] & decimallongitude < lon_limits[1],
                      !(decimallatitude > lat_limits[2] | decimallongitude < lon_limits[2]),
                      !(decimallatitude > lat_limits[3] & decimallongitude > lon_limits[3]),
                      !CellID%in%all_of(non_id),AltRas<alt_limit) %>% dplyr::select(CellID,x,y,AltRas,starts_with("LUC")) %>% distinct(CellID,.keep_all = T) %>% droplevels() %>% arrange(CellID) -> precs
  precs %>% pull(AltRas) %>% cut(.,breaks=seq(0,alt_limit+100,belts),labels=seq(0,alt_limit+100,belts)[-1], include.lowest=T) -> alt_cats ### aggregate records into altitudinal belts
  precs <- precs %>% ungroup() %>% mutate(AltCats = alt_cats)
  pp <- precs %>% mutate(Elevation=cut(AltRas,breaks=alt_breaks,labels=cat),.before=AltRas) %>% mutate(across(starts_with("LUC"), ~case_when(.x %in% c(1:6, 9, 12) ~"AU",.x %in% c(7)~"PF",.x %in% c(10,11)~"SV",.x %in% c(8)~"NF"))) %>% select(-AltRas,-x,-y) %>% pivot_longer(cols = starts_with("LUC"),names_to = "Year",values_to = "LandUse") %>% mutate(Year = as.numeric(sub("LUC_","",Year))) %>% group_by(Year,LandUse,Elevation) %>% summarize(Freq=n()) %>% filter(!is.na(LandUse)) %>%
    mutate(Category=paste(LandUse,Elevation,sep="_")) %>% 
    ggplot(aes(x=Year, y=Freq, fill=LandUse)) + 
    geom_area() + scale_fill_viridis_d(option="G",direction=-1,begin = 0.1,end=0.7) +
    theme(legend.position = "bottom",panel.background = element_blank(),panel.grid = element_blank(),axis.ticks = element_blank(),axis.text.x = element_text(size=10),axis.line = element_line()) + labs(y="Frequency",title="Land use in cloud forests through time") + geom_vline(xintercept = 1976) +
    facet_wrap(~ Category,strip.position = "top") +
    NULL
  ggsave(filename = here("plots/LUL_TimeSeries.pdf"),plot = pp)
}

do_climate_deprecated <- function() {
  data="output/Historical_data_FULL_GLASS.v.1.R"
  load(here(data))
  reduct_five %<>% separate(Accepted_Name,into = c("G","S","extra"),extra = "merge",sep=" ") %>% 
    select(-extra) %>% unite("ind",G:S,sep="_") %>% right_join(.,esta,by="ind")
  
  reduct_five %>% select(CellID:ind,year,starts_with("max_tas_") ) %>%
    mutate(year=as.numeric(sub("Year_","",year))) %>% select(-Altitude,-year,-ind) %>% distinct(CellID,.keep_all = TRUE) %>% 
    pivot_longer(-c(1), names_to = c("Variable","year"),names_sep = "_tas_",values_to = "Values") %>% 
    group_by(CellID,year) %>% 
    #mutate(TSI = (abs(Values[2] - min(Values))) / (abs(Values[2] - max((Values))))) %>% 
    rename(Climate = Values) %>% 
    select(-Variable) %>% 
    pivot_wider(names_from = year,names_prefix = "Climate",values_from = Climate) %>% 
    right_join(reduct_five %>% select(CellID:ind),.,by="CellID") -> aqui 
  
  aqui %>% group_by(ind) %>% 
    select(-Altitude) %>% 
    pivot_longer(-c(ind,CellID),names_to = "year",values_to = "Climate",names_prefix="Climate") %>% 
    mutate(tme = as.numeric(year) - 1979) -> aqui
  
  left_join(aqui,nulas %>% select(ind,RSamp_SES,RSamp_SES_cat,DSamp_SES_cat,CAT,RSamp_SES,DSamp_SES,RSamp_ES,DSamp_ES,DSamp_Inter_Cat),by="ind") %>% 
    mutate(RSamp_SES_cat = cut(RSamp_SES,breaks=c(-10,-1.96,-1.28,0,1.28,1.96,10),labels=c("Negative","Negative","Weak","Weak","Positive","Positive"),include.lowest=TRUE)) %>% 
    mutate(DSamp_SES_cat = cut(DSamp_SES,breaks=c(-10,-1.96,-1.28,0,1.28,1.96,10),labels=c("Sensitive","Sensitive","Robust","Robust","Sensitive","Sensitive"),include.lowest=TRUE)) %>% 
    filter(DSamp_SES_cat %in% c("Robust")) -> fita
  lmer(Climate ~ 1 + tme + tme:RSamp_SES_cat + (1|ind), data=fita, control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5))) -> fit_aver
  lmer(Climate ~ 1 + tme + (1|ind), data=fita, control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5))) -> fit_sum
  broom.mixed::tidy(fit_aver)
  broom.mixed::tidy(fit_aver) %>% flextable::flextable() %>% flextable::save_as_docx(path = here("plots_tables/cmi_glmm.docx"))
  ggpredict(fit_aver, type = "random",terms=c("tme [0,5,10,15,20,30,40]","RSamp_SES_cat"),interval = "confidence") %>% tibble() -> predict_fit
  predict_fit %>% 
    ggplot(aes(x=x,y=predicted,color=group,fill=group)) +
    geom_line() +
    geom_ribbon(aes(ymin = conf.low,ymax=conf.high),alpha=0.1,color="NA") +
    scale_color_met_d(name="Hiroshige",direction = -1)+
    scale_fill_met_d(name="Hiroshige",direction = -1) +
    geom_line(data=ggpredict(fit_sum, type = "random",terms=c("tme [all]"),interval = "confidence") %>% tibble(),aes(x=x,y=predicted),color="black",inherit.aes = FALSE) +
    theme(legend.position = "bottom",panel.background = element_blank(),axis.line = element_line(),legend.background = element_blank(),
          axis.title = element_text(family="EB Garamond"),
          axis.text = element_text(family="EB Garamond"),
          #axis.text.x = element_blank(),
          legend.spacing.y = unit(c(1.2),"cm"),
          legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond")) +
    labs(x="Time",y="Predicted mean tas") +
    NULL
  
  ggsave(filename = "plots_tables/cmi_glmm.pdf",plot=last_plot())
  
  
}

do_Granger_deprecated <- function(data="output/GrangerInput_Mid.Rdata",output="output/GrangerFinal_Mid.Rdata") {
  load(here(data))
  to_plot %>% ungroup() %>% pivot_wider(names_from = "Variable",values_from = "Series") %>% group_by(Resolved_ACCEPTED) %>% nest() %>% rowwise() %>% 
    ## Land use
    mutate(gci_luc = list(granger_mult(x.1 = data$LUC,y = data$Observed,n_x = 1))) %>% 
    ## Temperature
    mutate(gci_temps = list(granger_mult(x.1 = data$AnnTmean,x.2 = data$SeasTmean,x.3=data$MaxTmax,x.4 = data$MinTmin, y = data$Observed,n_x = 4))) %>% 
    ## Precipitation
    mutate(gci_precs = list(granger_mult(x.1 = data$AnnPrec,x.2 = data$SeasPrec,x.3=data$MaxPrec,x.4 = data$MinPrec, y = data$Observed,n_x = 4))) %>% 
    ## Overall climate
    mutate(gci_climate = list(granger_mult(data$AnnTmean,x.2 = data$SeasTmean,x.3=data$MaxTmax,x.4 = data$MinTmin, x.5 = data$AnnPrec,x.6 = data$SeasPrec,x.7=data$MaxPrec,x.8 = data$MinPrec,y = data$Observed,n_x = 8))) %>% 
    ## Overall climate | land use
    mutate(pgci_climate = list(granger_mult(data$AnnTmean,x.2 = data$SeasTmean,x.3=data$MaxTmax,x.4 = data$MinTmin, x.5 = data$AnnPrec,x.6 = data$SeasPrec,x.7=data$MaxPrec,x.8 = data$MinPrec,z = data$LUC, y = data$Observed,n_x = 8))) %>% 
    ## Climate plus land use
    mutate(gci_all = list(granger_mult(data$AnnTmean,x.2 = data$SeasTmean,x.3=data$MaxTmax,x.4 = data$MinTmin, x.5 = data$AnnPrec,x.6 = data$SeasPrec,x.7=data$MaxPrec,x.8 = data$MinPrec,x.9 = data$LUC, y = data$Observed,n_x = 5))) -> final
  save(final,file=here(output))
}

plot_granger_deprecated <- function(data="output_data/GrangerFinal_Min.Rdata",output="plots/Granger_Min.pdf",TRY="output_data/growth_try.Rdata",  data_lul="output_data/Historical_dataLULC.R") {
  load(here(data))
  load(here(data_lul))
  load(here(TRY))
  load(here("output_data/alt_species.Rdata"))
  lucs_tib %>% distinct(Resolved_ACCEPTED) -> species_list
  growth_try %>% filter(Resolved_ACCEPTED %in% all_of(species_list$Resolved_ACCEPTED)) %>% group_by(Resolved_ACCEPTED) %>% nest() -> target_try
  target_try %>% rowwise() %>% mutate(LifeForm=(Mode(data$OrigValueStr))) -> target_try
  left_join(species_list,target_try,"Resolved_ACCEPTED") %>% select(-data) %>% mutate(LifeForm=replace_na(LifeForm,"Unknown")) %>% mutate(LifeForm=factor(LifeForm)) -> species_list
  final %>% unnest(gci_luc) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","gci_luc","prob_luc")) %>% inner_join(.,{final %>% unnest(gci_temps) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","gci_temps","prob_temps"))},by="Resolved_ACCEPTED") %>% 
    inner_join(.,{final %>% unnest(gci_precs) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","gci_precs","prob_precs"))},by="Resolved_ACCEPTED") %>% 
    inner_join(.,{final %>% unnest(gci_climate) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","gci_climate","prob_climate"))},by="Resolved_ACCEPTED") %>% 
    inner_join(.,{final %>% unnest(pgci_climate) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","pgci_climate","prob_p_climate"))},by="Resolved_ACCEPTED") %>% 
    inner_join(.,{final %>% unnest(gci_all) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","gci_all","prob_all"))},by="Resolved_ACCEPTED") %>% 
    inner_join(.,species_list,by="Resolved_ACCEPTED") %>% inner_join(.,alt_species,by="Resolved_ACCEPTED") %>% ungroup() %>% pivot_longer(cols = contains("gci"),names_to = "Series",values_to = "gci") %>% mutate(AltCat=fct_relevel(AltCat,"Low","Mid Low","Mid High","High")) %>% pivot_wider(names_from = "Series",values_from = "gci") %>% mutate(gci_luc = case_when(prob_luc > 0.05 ~ as.numeric(NA),prob_luc <= 0.05 ~ gci_luc),
                                                                                                                                                                                                                                                                                                                                                       gci_temps = case_when(prob_temps > 0.05 ~ as.numeric(NA),prob_temps <= 0.05 ~ gci_temps),
                                                                                                                                                                                                                                                                                                                                                       gci_precs = case_when(prob_precs > 0.05 ~ as.numeric(NA),prob_precs <= 0.05 ~ gci_precs),
                                                                                                                                                                                                                                                                                                                                                       gci_climate = case_when(prob_climate > 0.05 ~ as.numeric(NA),prob_climate <= 0.05 ~ gci_climate),
                                                                                                                                                                                                                                                                                                                                                       pgci_climate = case_when(prob_p_climate > 0.05 ~ as.numeric(NA),prob_p_climate <= 0.05 ~ pgci_climate),
                                                                                                                                                                                                                                                                                                                                                       gci_all = case_when(prob_all > 0.05 ~ as.numeric(NA),prob_all <= 0.05 ~ gci_all)
    ) %>% pivot_longer(contains("gci"),names_to = "Series",values_to = "gci") -> out_plot
  out_plot %>% group_by(Series) %>% filter(!is.na(gci)) %>% summarise(n=n()) %>% mutate(names=paste(Series," (",n,")",sep="")) %>% pull(names) -> names_series
  pp <- out_plot %>% #mutate(Series=fct_relevel(Series,"gci_all","gci_luc","gci_precs","gci_temps","gci_climate","pgci_climate")) %>%  
    ggplot(aes(x=LifeForm,y=gci,fill=Series)) + 
    geom_boxplot(outlier.size=0.3,position="dodge2",size=0.2) +
    labs(title="Causal influence of historical climate change and vegetation loss on species' elevational ranges",x="Lifeform",y="Granger Causality Index") + 
    theme(panel.background = element_blank(),panel.grid = element_blank(),axis.line = element_line(),legend.position = "right",axis.text.x = element_text(size = 10),legend.text = element_text()) + scale_fill_viridis_d(option="G",begin=0.4,name="Series",labels=names_series) + 
    ylim(c(0,4)) +
    NULL
  ggsave(filename = here(output),plot = pp)
}

plot_granger_alternative_deprecated <- function(data="output_data/GrangerFinal_Min.Rdata",output="plots/Granger_Mid_Alternative.pdf",TRY="output_data/growth_try.Rdata",  data_lul="output_data/Historical_dataLULC.R") {
  load(here(data))
  load(here(data_lul))
  load(here(TRY))
  load(here("output_data/alt_species.Rdata"))
  lucs_tib %>% distinct(Resolved_ACCEPTED) -> species_list
  growth_try %>% filter(Resolved_ACCEPTED %in% all_of(species_list$Resolved_ACCEPTED)) %>% group_by(Resolved_ACCEPTED) %>% nest() -> target_try
  target_try %>% rowwise() %>% mutate(LifeForm=(Mode(data$OrigValueStr))) -> target_try
  left_join(species_list,target_try,"Resolved_ACCEPTED") %>% select(-data) %>% mutate(LifeForm=replace_na(LifeForm,"Unknown")) %>% mutate(LifeForm=factor(LifeForm)) -> species_list
  final %>% unnest(gci_luc) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","gci_luc","prob_luc")) %>% inner_join(.,{final %>% unnest(gci_temps) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","gci_temps","prob_temps"))},by="Resolved_ACCEPTED") %>% 
    inner_join(.,{final %>% unnest(gci_precs) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","gci_precs","prob_precs"))},by="Resolved_ACCEPTED") %>% 
    inner_join(.,{final %>% unnest(gci_climate) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","gci_climate","prob_climate"))},by="Resolved_ACCEPTED") %>% 
    inner_join(.,{final %>% unnest(pgci_climate) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","pgci_climate","prob_p_climate"))},by="Resolved_ACCEPTED") %>% 
    inner_join(.,{final %>% unnest(gci_all) %>% select(orig,prob) %>% rename_with(~c("Resolved_ACCEPTED","gci_all","prob_all"))},by="Resolved_ACCEPTED") %>% 
    inner_join(.,species_list,by="Resolved_ACCEPTED") %>% inner_join(.,alt_species,by="Resolved_ACCEPTED") %>% ungroup() %>% pivot_longer(cols = contains("gci"),names_to = "Series",values_to = "gci") %>% mutate(AltCat=fct_relevel(AltCat,"Low","Mid Low","Mid High","High")) %>% pivot_wider(names_from = "Series",values_from = "gci") %>% mutate(gci_luc = case_when(prob_luc > 0.05 ~ as.numeric(NA),prob_luc <= 0.05 ~ gci_luc),
                                                                                                                                                                                                                                                                                                                                                       gci_temps = case_when(prob_temps > 0.05 ~ as.numeric(NA),prob_temps <= 0.05 ~ gci_temps),
                                                                                                                                                                                                                                                                                                                                                       gci_precs = case_when(prob_precs > 0.05 ~ as.numeric(NA),prob_precs <= 0.05 ~ gci_precs),
                                                                                                                                                                                                                                                                                                                                                       gci_climate = case_when(prob_climate > 0.05 ~ as.numeric(NA),prob_climate <= 0.05 ~ gci_climate),
                                                                                                                                                                                                                                                                                                                                                       pgci_climate = case_when(prob_p_climate > 0.05 ~ as.numeric(NA),prob_p_climate <= 0.05 ~ pgci_climate),
                                                                                                                                                                                                                                                                                                                                                       gci_all = case_when(prob_all > 0.05 ~ as.numeric(NA),prob_all <= 0.05 ~ gci_all)
    ) %>% pivot_longer(contains("gci"),names_to = "Series",values_to = "gci") -> out_plot
  counts <- out_plot %>% filter(!is.na(gci)) %>% count(Series)
  meanas <- out_plot %>% filter(!is.na(gci)) %>% group_by(Series) %>% summarise(n=median(gci))
  names_series <- meanas %>% mutate(Names=str_to_sentence(sub("gci_","",Series)))
  names_series$Names[6] <- "Climate|LUC"
  names_series$Names[3] <- "LUC"
  pp <- out_plot %>% filter(!is.na(gci)) %>% 
    ggplot(aes(x = Series, y = gci,fill=Series)) + 
    ggdist::stat_halfeye(
      adjust = .9, 
      width = .6, 
      .width = 0, 
      alpha=0.7,
      slab_colour = "black",
      slab_size = .2,
      justification = -.4, 
      height = 2) + 
    geom_boxplot(width = .4, outlier.shape = NA,size=.2,fatten = 7,alpha=0.7) +
    geom_point(fill="grey50",size = 2,alpha = .4,shape = 21,stroke = 0.2,position = position_jitter(seed = 1, width = .2)) + 
    coord_flip() +
    scale_fill_viridis_d(option="G",begin=0.4,name="Series",labels=names_series)  + 
    scale_y_log10() + theme(panel.background = element_blank(),panel.grid = element_blank(),axis.line = element_line(),legend.position = "",axis.text.x = element_text(size = 10),legend.text = element_text()) + annotate(
      geom = "text", x = (1:6)+.2, y = 15, 
      label = paste0("n = ",counts$n), hjust = 0, vjust = -.8, size = 2,fontface = "bold",
    ) + annotate(
      geom = "text", x = (1:6)+.2, y = meanas$n, fontface = "bold",
      label = round(meanas$n,2), hjust = .5, vjust = -.8, size = 3
    ) + scale_x_discrete(labels= names_series$Names) + labs(x="",y="Granger Causality Index") +
    NULL
  ggsave(filename = here(output),plot = pp)
}

plot_LMM_deprecated <- function(data="output_data/LMM_modelsMin.Rdata",cual="Min",outplot="plots/LMM_modelsMin.pdf",outtable="tables/LMM_modelsMin.docx"){
  load(here(data))
  print(anova(models[[1]],models[[2]],models[[3]],models[[4]]))
  as.numeric(readline("Select model")) -> quien
  target = models[[quien]]
  broom.mixed::tidy(target) %>% flextable::flextable() %>% flextable::save_as_docx(path = here(outtable))
  trends <- coef(target)$Resolved_ACCEPTED %>% as_tibble() %>% ggplot(aes(x=wave01)) + 
    geom_histogram( aes(x = wave01, y = ..density..), fill=rocket(10)[6],binwidth = 10) +
    geom_label(aes(x=max(wave01,wave1) - 20, y=0.005, label="Pre-shift trend"), color=rocket(10)[7]) +
    geom_histogram( aes(x = wave1, y = -..density..), fill= mako(10)[6],binwidth = 10) +
    geom_label( aes(x=max(wave01,wave1) - 20, y=-0.005, label="Post-shift trend"), color=mako(10)[7]) +
    labs(x="Trend (m)",y="density") + theme(panel.background = element_blank(),panel.grid = element_blank(), axis.line = element_line()) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    NULL
  
  pre <- coef(target)$Resolved_ACCEPTED %>% as_tibble() %>% ggplot(aes(x=wave01,y=`(Intercept)`)) + stat_density_2d(aes(fill = ..level..), geom = "polygon",bins=20) + labs(x="Pre-shift trend (m/time)",y="Initial elevation (m)") + scale_fill_viridis(option="F",direction=-1,begin=0.4) + theme(panel.background = element_blank(),panel.grid = element_blank(), axis.line = element_line()) + geom_vline(xintercept = 0,linetype="dashed")
  
  post <- coef(target)$Resolved_ACCEPTED %>% as_tibble() %>% ggplot(aes(x=wave1,y=`(Intercept)`)) + stat_density_2d(aes(fill = ..level..), geom = "polygon",bins=20) + labs(x="Post-shift trend (m/time)",y="Initial elevation (m)") + scale_fill_viridis(option="G",direction=-1,begin=0.4) + theme(panel.background = element_blank(),panel.grid = element_blank(), axis.line = element_line()) + geom_vline(xintercept = 0,linetype="dashed")
  
  layout <- "
AAAA
BBCC
"
  
  pp <- wrap_plots(list(trends,pre,post)) +
    plot_layout(design = layout,byrow=T,nrow = 2,ncol=2,guides = "collect",tag_level = 'new') + 
    plot_annotation(title = paste0("Species' elevational trends (",cual,"-point)"),subtitle="Turning point 1975-1976", caption= "Ramírez-Barahona et al.") & theme(legend.position = '')
  ggsave(filename = here(outplot),plot = pp)
  
}

alt_downsamp_nulls_deprecated <- function(data="output/Historical_data_FULL_LUL.v.1.R",n_sim=100,start_year=1979,records=100) {
  load(here(data))
  reduct_five %>% select(CellID:Accepted_Name,ID,class:family.x,decimalLatitude:BIOME) %>% 
    mutate(YearCat = as.numeric(sub("Year_","",year)),.after=1) %>% 
    filter(!is.na(Altitude),YearCat >= start_year) %>% 
    distinct(Accepted_Name,CellID,year,.keep_all = T) %>% group_by(Accepted_Name) -> esta
  esta %>% summarise(N=n()) -> tata
  target <- tata %>% filter(N >= records) %>% pull(Accepted_Name) %>% as.character(.)
  ## request parallelization if simulations are set to TRUE
  esta %>% filter(Accepted_Name %in% all_of(target)) %>% 
    mutate(AltCat=cut(Altitude,breaks=seq(0,5600,500)),.before=Altitude) %>% 
    mutate(FiveYear=cut(YearCat,breaks=seq(1979,2019,5),include.lowest = TRUE),.before=Altitude) %>% 
    count(AltCat,FiveYear,Accepted_Name) %>%
    mutate(n = floor(ifelse(n <= quantile(n,0.95), n, quantile(n,0.95))))-> N
  esta %>% filter(Accepted_Name %in% all_of(target)) %>% mutate(AltCat=cut(Altitude,breaks=seq(0,5600,500)),.before=Altitude) %>% 
    mutate(FiveYear=cut(YearCat,breaks=seq(1979,2019,5),include.lowest = TRUE),.before=Altitude) %>% 
    group_by(AltCat,FiveYear,Accepted_Name) %>% nest() %>% left_join(.,N,by=c("AltCat","Accepted_Name","FiveYear")) -> esta
  plan(multisession,workers=4)
  cat("     Performing",n_sim,"simulations","\n")
  ## simulate data by resampling
  simulations <- 1:n_sim %>% future_map_dfr(function(j){
    esta %>% mutate(data=map2(data,n,replace=F,sample_n)) %>% unnest(data) %>% mutate(Accepted_Name=paste(Accepted_Name,"_sim_",j,sep="")) %>% arrange(Accepted_Name,YearCat) -> rarified
    rarified %>% ungroup() %>% group_by(Accepted_Name,YearCat) %>% 
      summarise(Variable=mean(Altitude,na.rm=T),MaxVariable = quantile(Altitude,probs=0.95,na.rm=T),MinVariable=quantile(Altitude,probs=0.05,na.rm=T),n=n(),class=first(class)) %>% mutate_if(is.numeric, ~na_if(., -Inf)) %>% mutate_if(is.numeric, ~na_if(., Inf)) %>% arrange(Accepted_Name,YearCat) %>% 
      mutate(RollAlt=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=mean,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
      mutate(RollAlt_Max=unlist(tapply(MaxVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
      mutate(RollAlt_Min=unlist(tapply(MinVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
      mutate(RollAlt_error=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=seme,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
      mutate(RollAlt_Max_error=unlist(tapply(MaxVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = seme,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
      mutate(RollAlt_Min_error=unlist(tapply(MinVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = seme,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
      
      mutate(Rolln=unlist(tapply(n,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=sum,na.rm=T,partial=T,align="center")}),use.names = F)) %>% filter(Rolln >= min_by_year) -> r1
  },.progress = TRUE,.options = furrr_options(seed = 123)) ## set seed to avoid warning!
  cat("     Summarizing","simulations","\n")
  simulations %>% ungroup() %>% nest_by(Accepted_Name) -> nulls_m
  cat("Saving models (final)","\n")
  nombre_arch <- here("interim","Splines_Down_Alt_Partial_nulls_1.Rdata")
  save(nulls_m,file=nombre_arch) 
  rm(nulls_m)
  rm(simulations)
  #})
  plan(sequential)
  ## revert ot unparallelized session
}

### RANDOM CHUNKS OF CODE (NOT DELETING DUE TO OBSSESIVENESS) ####
## There is some testing grounds here....
{ 
  mero %>% separate(Accepted_Name,into=c("a","b","c"),sep = " ",extra="merge") %>% unite("ind",a:b,sep="_") -> aver
  aver %<>% 
    filter(c !="var. grandifolia (M.C.Johnst. & L.A.Johnst.) A.Pool") %>% filter(c!="var. ramirezii (Greenm.) Gonsoulin")
  dd %>% right_join(aver,.,by="ind") -> datos
  datos %<>% filter(type=="resamp") %>% group_by(ind)
  datos %>% select(estimate,mean,sd,SES,ES,type)
  fit_manken <- lmer(ES ~ 1 + (1|ind) + manken_tas + (0 + manken_tas|ind),data = datos, control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4)))
  broom.mixed::tidy(fit_manken)
  equatiomatic::extract_eq(fit_manken)
  
  
  una <- datos  %>% 
    ggplot(aes(y=manken_z,x=ind)) +
    #geom_segment(aes(y=quant10,yend=quant90,x=ind,xend=ind),alpha=0.3,size=0.2) +
    #geom_segment(aes(y=quant25,yend=quant75,x=ind,xend=ind),alpha=0.5) +
    geom_hline(yintercept = 0,linetype="dashed",color="grey80") +
    geom_point(aes(fill=ifelse(manken_z<=0,"1","2")),shape=21,stroke=0.1) +
    coord_flip()+
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line(),
          axis.line.y = element_blank(),
          legend.background = element_blank(),
          legend.key.height = unit(0.01,"cm"),
          legend.key.width = unit(0.01,"cm"),
          axis.title = element_text(family="EB Garamond"),
          axis.text = element_text(family="EB Garamond"),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.spacing.y = unit(c(1.2),"cm"),
          legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond")) +
    labs(y="Change in climate moisture index",x="")
  
  dos <- datos %>%  
    ggplot(aes(y=SES,x=as.numeric(ind))) +
    geom_boxplot(aes(fill=ifelse(manken_z<=0,"1","2"))) +
    geom_point(aes(fill=ifelse(manken_z<=0,"1","2")),shape=21,stroke=0.1) +
    geom_hline(yintercept = 0,linetype="dashed",color="grey70") +
    coord_flip() +
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line(),
          axis.line.y = element_blank(),
          legend.background = element_blank(),
          legend.key.height = unit(0.01,"cm"),
          legend.key.width = unit(0.01,"cm"),
          axis.title = element_text(family="EB Garamond"),
          axis.text = element_text(family="EB Garamond"),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.spacing.y = unit(c(1.2),"cm"),
          legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond")) +
    labs(y=bquote("SES of elevation shifts"),x="")
  
  lm(data=datos,SES~manken_z) %>% summary()
  datos %>%
    ggplot(aes(y=SES,x=(manken_z))) +
    geom_point() +
    scale_color_manual(values=MetBrewer::met.brewer("Hiroshige",n=16)) +
    geom_smooth(color="black",method="lm") +
    geom_hline(yintercept = 0,linetype="dashed",color="grey70") +
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line(),
          legend.background = element_blank(),
          legend.key.height = unit(0.001,"cm"),
          legend.key.width = unit(0.001,"cm"),
          axis.title = element_text(family="EB Garamond"),
          axis.text = element_text(family="EB Garamond"),
          #axis.ticks.y=element_blank(),
          #axis.text.x=element_blank(),
          legend.spacing.y = unit(c(1.2),"cm"),
          legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond")) +
    labs(y="Species elevation shifts (SES)",x="Multivariate MannKendall statistic (mean per species)")
  
  
  
  layout <- "
AABB
CCCC
"
  wrap_plots(list(una,dos,tres)) +
    plot_layout(design = layout,byrow=T,nrow = 2,guides="keep",ncol=2,tag_level = 'new') +
    plot_annotation(tag_levels = 'a',title="") & theme(legend.position = 'none',plot.tag=element_text(family="EB Garamond",size=15,face="bold"),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))
  
}
{
  splines="output/Splines_v.2.Rdata"
  load(here(splines))
  
  aver <- mean_alts %>% select(Accepted_Name,data) %>% unnest(data) %>% 
    select(Accepted_Name,YearCat,RollAlt) %>% 
    mutate(YearCat = paste0("year_",YearCat)) %>% 
    pivot_wider(names_from = "YearCat",values_from = "RollAlt") %>% 
    rowwise() 
  aver %>% mutate(baseline = mean(c_across(year_1980:year_1990),na.rm=T),.after="Accepted_Name") %>% 
    mutate(finish = mean(c_across(year_2000:year_2010),na.rm=T),.after="Accepted_Name") %>% 
    #mutate(across(starts_with("year"),~.x - baseline)) %>% 
    mutate(Overall = finish - baseline,.after="Accepted_Name") %>% 
    #pivot_longer(starts_with("year"),names_to = "year",values_to = "diff") %>% 
    #mutate(year=as.numeric(sub("year_","",year))) %>% 
    select(1:4) %>% 
    # pivot_longer(-1,names_to = "Where",values_to = "value") %>% 
    #mutate(Time=factor(ifelse(Where=="baseline",1,2))) %>% 
    # filter(Where!="Overall") %>% 
    mutate(Cat_Base = cut(baseline,seq(0,4000,200),labels=seq(0,4000,200)[-1], include.lowest=TRUE)) %>% 
    mutate(Cat_End = cut(finish,seq(0,4000,200),labels=seq(0,4000,200)[-1], include.lowest=TRUE)) %>% 
    group_by(Cat_Base) %>% 
    #filter(baseline >= 1000 & baseline < 2500) %>% 
    #filter(finish >= 1000 & finish < 2600) %>% 
    #filter(Overall > 50) %>% 
    filter(Overall < -50) %>% 
    ggplot() +
    geom_segment(aes(x=1,xend=2,y=baseline,yend=finish),alpha=0.5) +
    geom_point(aes(x=1,y=baseline)) +
    ylim(c(0,4000)) +
    geom_point(aes(x=2,y=finish)) +
    #geom_line() +
    #geom_boxplot() +
    theme(legend.position = "none") +
    NULL
  #####
  
  
  load(here("interim",lul))
  aver
  aver %<>% select(-YearCat) %>%  left_join(., reduct_five %>% select(CellID,Altitude) %>% distinct(CellID,.keep_all = T),by="CellID") %>% mutate(AltCat=cut(Altitude,breaks=seq(0,alt_limit+100,250),labels=seq(0,alt_limit+100,250)[-1], include.lowest=T))
  aver %<>% #filter(CellID==8977602) %>% 
    #rowwise() %>% 
    #mutate( Sum = sum(c_across(Roll_AU:Roll_NF))  ) %>% 
    # mutate(across(Roll_AU:Roll_NF,~.x/Sum)) %>% 
    mutate(across(Roll_AU:Roll_NF, ~.x/Rolllul_n*100))
  
  aver %>% filter(CellID==7014544) %>% 
    mutate(Decade = ifelse(year < 1990,"First",ifelse( year < 2000,"Second", ifelse(year < 2010 ,"Third", "Fourth" )))) %>% filter(year %in% c(1980,2010)) -> stats
  
  stats %>% select(year,Roll_AU:Roll_NF,  AltCat) %>% pivot_longer(Roll_AU:Roll_NF,names_to = "Var",values_to = "Val") %>% group_by(year,AltCat,Var) %>% summarize(Val=mean(Val)) %>% 
    #unite("x",c(year,Var)) %>% 
    pivot_wider(values_from = Val, names_from = year)
  
}
{
  
  ## 21 june I have the series for each grid-cell.
  library(lattice)
  library(sme)
  library(lme4)
  
  data = "interim/ClimateSeries_v.2.RData" #### series climáticas por gradilla
  load(data)
  aqui
  
  splines="output/Splines_v.2.Rdata" ### Mean altitude by year per species
  load(here(splines))
  mean_alts$data[[1]]
  all_models$resamp_models %>% ungroup() %>% filter(Derivative=="Raw")
  
  
  
  panelB.1 <- ggplot() +
    geom_segment(aes(y=0.5,yend=-1,x=-0.01,xend=-0.01),color="grey30") +
    geom_segment(aes(y=0,yend=0,x=-4,xend=5),color="grey30",size=0.1) +
    stat_dots(data  = dd,aes(x = SES,side=ifelse(type=="resamp","top","bottom"),
                             color=stat(x < 0),fill=stat(x<0)),scale=1) +
    scale_fill_manual(values= colorspace::darken(MetBrewer::met.brewer("Hiroshige",n=8,direction=1)[c(3,4)],amount=0.1),
                      guide=guide_legend(reverse = T,ncol = 1,title.position="top", title.hjust = 0, label.position = "top"),name="") +
    scale_color_manual(values= colorspace::darken(MetBrewer::met.brewer("Hiroshige",n=8,direction=1)[c(3,4)],amount=0.2),
                       guide=guide_legend(reverse = T,ncol = 1,title.position="top", title.hjust = 0, label.position = "top"),name="") +
    geom_text(inherit.aes = T,aes(label=paste0("Resampling"),x=c(3),y=0.3), 
              position=position_nudge(y=c(0)), size=3,family="EB Garamond") +
    geom_text(inherit.aes = T,aes(label=paste0("Downsampling"),x=c(3),y=-0.3), 
              position=position_nudge(y=c(0)), size=3,family="EB Garamond") +
    coord_cartesian(xlim = c(-4, 5),ylim=c(-1,0.3)) +
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line(),
          legend.background = element_blank(),
          legend.key.height = unit(0.01,"cm"),legend.key.width = unit(0.01,"cm"),
          axis.title = element_text(family="EB Garamond"),
          axis.text = element_text(family="EB Garamond"),
          axis.line.y=element_blank(),axis.ticks.y = element_blank(),
          axis.text.y=element_blank(),
          legend.spacing.y = unit(c(1.2),"cm"),
          legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond")) +
    labs(x=bquote("Standardized Effect Size"),y="") +
    scale_x_continuous(breaks=seq(-3,5,1))+
    
    NULL
  
  
  
  
  
  data="output/Historical_data_FULL_LUL.v.1.R"
  load(here(data))
  # reduct_five %>% distinct(CellID,Accepted_Name,.keep_all = T)
  ### ASSIGN TEMPERATURE TO EVERY RECORD BASED ON THE COLLECTION YEAR, THEN SUMMARIZE BY SPECIES
  variable="mean_tas_"
  reduct_five %>% 
    mutate(year=as.numeric(sub("Year_","",year))) %>% distinct(CellID,Accepted_Name,year,.keep_all = T) %>% 
    select(CellID,Accepted_Name,year,ends_with("_1980")) -> to_temp
  
  
  prcomp(to_temp[-c(1:4)])
  
  
  
  to_temp %>% pull(year) %>% as.character() -> target_year
  names(to_temp) %>% sub(variable,"",.) -> matching_year
  #baseline <-  to_temp %>% select(1:4,grep(paste(1979:1989,collapse="|"),names(.)))
  #baseline %<>% rowwise() %>% mutate(baseline = mean(c_across(starts_with(variable))))
  #baseline %>% select(baseline) %>% bind_cols(.,to_temp) %>% 
  #filter(year %in% all_of(matching_year[-c(1:3)])) -> to_temp
  #to_temp %<>% ungroup() %>% mutate(across(starts_with(variable), ~.x - baseline))
  target_year
  matching_year
  match_index <- match(target_year,matching_year)
  val <- c()
  for(i in 1:nrow(to_temp)){ 
    cat(i,"\r")
    to_temp[i,match_index[i]] %>% pull() -> val[i]
  }
  
  to_temp %>%  mutate(diff= mean_tas_2019 - mean_tas_1979,.after=year) %>% group_by(Accepted_Name) %>% summarise(var_year=mean(diff)) %>% rename("ind"=Accepted_Name) %>%  mutate(ind=as.character(ind)) %>% separate(ind,into = c("G","S","extra"),extra = "merge",sep=" ") %>% select(-extra) %>% unite("ind",G:S,sep="_") %>%  right_join(.,esta,by="ind") -> esta_otra
  
  to_temp %>% mutate(var_year = val,.after=year) %>% group_by(Accepted_Name) %>% summarise(var_year=mean(var_year)) %>% rename("ind"=Accepted_Name) %>%  mutate(ind=as.character(ind)) %>% separate(ind,into = c("G","S","extra"),extra = "merge",sep=" ") %>% select(-extra) %>% unite("ind",G:S,sep="_") %>%  right_join(.,esta,by="ind") -> esta_otra
  
  esta_otra %>% ggplot(aes(x=tme1,y=error.y,fill=var_year)) + 
    #stat_density_2d(aes(fill = ..level..), geom = "polygon",bins=20) + 
    geom_point(shape=21,size=1.5,stroke=.1)+
    scale_fill_stepsn(n.breaks=10,colors= rev(MetBrewer::met.brewer("Hiroshige")) ) + 
    theme(panel.background = element_blank(),panel.grid = element_blank(), axis.line = element_line()) + 
    geom_vline(xintercept = 0,linetype="dashed")
  
  
  
  
  
  
  broom.mixed::tidy(fit1) %>% select(estimate) %>% .[1,] %>% as_vector() -> intercept
  broom.mixed::tidy(fit1) %>% select(estimate) %>% .[2,] %>% as_vector() -> Btme1
  broom.mixed::tidy(fit1) %>% select(estimate) %>% .[3,] %>% as_vector() -> Berror
  broom.mixed::tidy(fit1) %>% select(estimate) %>% .[4,] %>% as_vector() -> Btme2
  broom.mixed::tidy(fit1) %>% select(estimate) %>% .[5,] %>% as_vector() -> Btme3
  broom.mixed::tidy(fit1) %>% select(estimate) %>% .[6,] %>% as_vector() -> Btme1_error
  broom.mixed::tidy(fit1) %>% select(estimate) %>% .[7,] %>% as_vector() -> Btme2_error
  broom.mixed::tidy(fit1) %>% select(estimate) %>% .[8,] %>% as_vector() -> Btme3_error
  
  tme1 = 0:10 #0:10 1989 - 1989
  tme2 = 0 #1:10 1990- 1999
  tme3 = 0 #1:16 2000-2015
  
  probas <- c(0.01,0.25,0.5,0.75,0.99)
  x = list()
  for(i in 1:length(probas)){ 
    error = otra %>% summarise(error=quantile(error,probs=probas[i])) %>% pull(error)
    x[[i]] <- intercept + Btme1*tme1 + Berror*error + Btme2*tme2 + Btme3*tme3 + Btme1_error*tme1*error + Btme2_error*tme2*error + Btme3_error*tme3*error
  }
  names(x) <- paste("prob",probas,sep="_")
  colores = MetBrewer::met.brewer("Hiroshige",n=length(probas))
  do.call(cbind,x) %>% as_tibble() %>% mutate(year=0:10) %>% 
    mutate(across(1:5, ~.x - lag(.x),.names="S{col}")) -> x 
  x %>% ggplot(aes(x=year,y=prob_0.99)) + 
    geom_line(color=colores[1]) +
    geom_line(aes(y=prob_0.75),color=colores[2]) +
    geom_line(aes(y=prob_0.5),color=colores[3]) +
    geom_line(aes(y=prob_0.25),color=colores[4]) +
    geom_line(aes(y=prob_0.01),color=colores[5]) +
    scale_x_continuous(labels=1979:1989,breaks=0:10) +
    #scale_x_continuous(labels=1990:1999,breaks=1:10) +
    #scale_x_continuous(labels=2000:2015,breaks=1:16) +
    annotate("text",x=nrow(x) +  1,y=apply(x[-1,1:length(probas)], FUN= function(y) tail(y,1),MARGIN=2),label=as_vector(round(x[2,(length(probas) +2 ):ncol(x)],2)),color=rev(colores),fontface="bold",family="EB Garamond") +
    theme(panel.background = element_blank(),axis.line = element_line(),
          axis.text.x = element_text(size=12,family="EB Garamond"),
          axis.text.y = element_text(family="EB Garamond"),
          axis.title = element_text(family="EB Garamond"),
          plot.title = element_text(size=18,face="bold",family="EB Garamond"),
          plot.subtitle = element_text(family="EB Garamond"),
          plot.caption = element_text(family="EB Garamond")
    ) +
    labs(x="",y="Altitude",title="Shifts in estimated altitude",
         caption="Note: values predicted under a linear mixed model with interaction terms",
         subtitle="Slopes estimated under different levels of uncertainty (standard error)"
    ) +
    NULL
}
{ 
  aver %>% select(y,ind,year) %>% group_by(ind) %>% arrange(year) %>% filter(!is.na(y)) %>% #summarise(start= nth(y,5L),end=nth(y,-5L))%>% 
    pivot_wider(names_from = year,values_from = y) %>% #rowwise() %>% mutate(start=`1985`,.after=1) %>% mutate(end=`2015`,.after=2) %>% select(1:3) %>% 
    ungroup() %>% #mutate(diff= end - start) %>% 
    summarise(across(-1,range,na.rm=T)) %>% mutate(variable=c("Min","Max"),.before=1) %>% pivot_longer(-1,names_to = "year",values_to = "sd") %>% ggplot(aes(x=year,y=sd)) + geom_point()
  
  aver %>% select(y,ind,year) %>% group_by(ind) %>% arrange(year) %>% filter(!is.na(y)) %>% 
    ggplot(aes(x=year,y=y,color=ind),alpha=0.3) + geom_line()
  
  
  aver %>% filter(year <= 1990) %>% group_by(Accepted_Name,CellID) %>% count(Values) %>% pivot_wider(names_from = Values,values_from = n) %>% mutate(across(everything(),replace_na,0)) %>% rowwise() %>% mutate(total=sum(c_across(everything()))) %>% rename_at(3:7,~paste0("Ref_",.x)) %>% group_by(CellID) %>% 
    #summarise(across(starts_with("Ref_"),sum)) %>% 
    arrange(desc(Ref_total))
  
  
  reduct_five %>% distinct(Accepted_Name) -> species_list
  all_models[[1]] %>% distinct(Accepted_Name) -> with_models
  with_models %>% mutate(Models="Model") %>% left_join(species_list,.,by="Accepted_Name") %>% mutate(Models=replace_na(Models,"No model")) -> species_list
  species_list %<>% separate(Accepted_Name,sep=" ",into=c("G","S","Author"),extra = "merge",remove = F) %>% select(-Author) %>% unite("binomial",G:S,sep="_")
  
  
  load(here(TRY))
  growth_try %>% rename("binomial"=Accepted_Name) %>% 
    filter(binomial %in% all_of(species_list$binomial)) %>% group_by(binomial) %>% nest() -> target_try
  target_try %>% rowwise() %>% mutate(LifeForm=(Mode(data$OrigValueStr))) -> target_try
  
  left_join(species_list,target_try,"binomial") %>% select(-data) %>% mutate(LifeForm=replace_na(LifeForm,"Unknown")) %>% mutate(LifeForm=factor(LifeForm)) -> species_list
  #levels(species_list$LifeForm) <- c("Aquatic","Climber","Epiphyte","Herb","Shrub","Shrub","Shrub","Tree","Unknown")
  
  
  
  
  
  
  select(Accepted_Name,CellID) %>% nest() -> tib_species
  tib_species %>% inner_join(.,species_list,by="Accepted_Name") -> tib_species
  
  ### This is to get the climate series in long format
  aqui %>% select(CellID,Var,YearCat,Roll_est) %>% pivot_wider(names_from = c(Var,YearCat),values_from = Roll_est) -> roll_ests
  names(roll_ests)[-1] = paste0("Rolling_",names(roll_ests)[-1])
  aqui %>% select(CellID,Var,YearCat,Values) %>% pivot_wider(names_from = c(Var,YearCat),values_from = Values) -> values
  names(values)[-1] = paste0("Raws_",names(values)[-1])
  aqui %>% select(CellID,Var,YearCat,Dif_Values) %>% pivot_wider(names_from = c(Var,YearCat),values_from = Dif_Values) -> Diffvalues
  names(Diffvalues)[-1] = paste0("RawsDif_",names(Diffvalues)[-1])
  aqui %>% select(CellID,Var,YearCat,Dif_Roll_est) %>% pivot_wider(names_from = c(Var,YearCat),values_from = Dif_Roll_est) -> Diffrolling
  names(Diffrolling)[-1] = paste0("RollingDif_",names(Diffrolling)[-1])
  roll_ests %>% 
    left_join(.,values,by="CellID")  %>% 
    left_join(.,Diffvalues,by="CellID")  %>% 
    left_join(.,Diffrolling,by="CellID") -> variables
  ####
  
  tib_species %<>% unnest(data) %>% group_by(CellID) %>% nest() %>% 
    left_join(.,variables,by="CellID")
  
  mutate(CellID=as.numeric(CellID)) %>% inner_join(tib_species %>% unnest(data),.,by="CellID") %>% summarise(across(.cols = -c(1:2), mean,na.rm=T)) -> sum_precs
  
  
  reduct_five %>% select(CellID:BIOME) %>% group_by(CellID) %>% nest() %>% left_join(.,variables,by="CellID") %>% unnest()
  
  
  #%>% left_join(.,variables,by="CellID")
  
  # %>% group_by(Resolved_ACCEPTED) %>% summarise(across(.cols = -c(1:3), freq_PF)) -> lul_una
  sum_precs
  lul_una %>% inner_join(sum_precs,.,by="Resolved_ACCEPTED") -> sum_precs
  
  all_models %>% ungroup() %>% filter(Type=="NULL",Derivative==quien) %>% group_by(Resolved_ACCEPTED) %>% summarise(across(starts_with("Year"), list(mean=mean,sd=sd))) -> sum_sims
  all_models %>% ungroup() %>% filter(Type=="OBS",Derivative==quien) -> observed
  sum_sims %>% pivot_longer(cols=-1,names_to = c("YearCat","Variable"),names_pattern = "Year_?(.*)_(.*)",values_to = "Estimate") %>% pivot_wider(names_from = Variable,values_from = Estimate) -> sims_vals
  observed %>% pivot_longer(cols=-c(Resolved_ACCEPTED,Type,Pvalue,Trend,Derivative),names_to = c("YearCat"),values_to = "Observed") %>% mutate(YearCat=as.numeric(sub("Year_","",YearCat))) %>% arrange(Resolved_ACCEPTED,YearCat) %>% 
    select(Observed) %>% bind_cols(sims_vals,.) %>% mutate(SES = (Observed - mean) / sd) -> to_plot
  sum_precs %>% pivot_longer(cols = -1,names_to = "YearCat",values_to = "Series" ) -> sum_precs
  
  to_plot %>% group_by(Resolved_ACCEPTED) %>% nest() %>% mutate(data=dif_clim(data,index = 4:5,flex_pt = 1976)) %>% unnest(data) -> to_plot
  sum_precs %>% separate(YearCat,into=c("Variable","YearCat"),sep="_") %>%inner_join(to_plot,.,by=c("Resolved_ACCEPTED","YearCat")) -> to_plot
  
  load(list_archs[4])
  all_models %>% ungroup() %>% filter(Type=="OBS",Derivative==quien) -> observed
  observed %>% pivot_longer(cols=-c(Resolved_ACCEPTED,Type,Pvalue,Trend,Derivative),names_to = c("YearCat"),values_to = "Observed") %>% mutate(YearCat=as.numeric(sub("Year_","",YearCat))) %>% arrange(Resolved_ACCEPTED,YearCat) %>% group_by(Resolved_ACCEPTED) %>% filter(YearCat < 1976) %>% summarise(Baseline=mean(Observed,na.rm=T)) %>% mutate(AltCat=case_when(Baseline < 600 ~ "Low", Baseline >=600 & Baseline < 1500 ~ "Mid Low", Baseline >= 1500 & Baseline < 2500 ~ "Mid High", Baseline >= 2500 ~ "High")) -> alt_species
  save(alt_species,file = here("output_data/alt_species.Rdata"))
  save(to_plot,file = here(output))
  #} 
  #"Observed"  "AnnPrec"   "AnnTmean"  
  #"MaxPrec"   "MaxTmax"   "MinPrec"   "MinTmin" 
  #"SeasPrec"  "SeasTmean" "LUC"
}
{
  data = "output/Splines_v.2.Rdata"
  load(data)
  predict_ci <- function(data,slopes,der,flex_point) {
    data %>% select(YearCat) %>% range -> years
    years = seq(years[1],years[2],1)
    una <- predict(slopes,x=years,se.fit=T,deriv=der)
    mean_ref <- una %>% slice(1:which(una$x == flex_point)) %>% summarise(mean_ref=mean(y)) %>% pull(mean_ref)
    una <- una %>% rowwise() %>% mutate(dif_mean_ref = y - mean_ref)
    return(una)
  }
  mean_alts %<>% rowwise(Accepted_Name) %>% 
    mutate(Slopes_seme = list(splines_mods(data$YearCat,data$RollAlt_error,jit = NULL,fit = NULL)))
  aver <- mean_alts %>% rowwise(Accepted_Name) %>%
    mutate(predict = list(predict_ci(data,Slopes,der=0,flex_point=1990)))
  aver %>% select(predict) %>% unnest(predict) %>% rename(YearCat=x) -> to_plot 
  
  aver %>% select(data) %>% unnest(data) %>% left_join(to_plot,.,by=c("Accepted_Name","YearCat")) -> to_plot
  
  uno <- to_plot %>% mutate(Year=((YearCat))) %>% 
    mutate(Alt=cut(y,breaks=seq(0,4600,200),labels=seq(0,4600,200)[-1], include.lowest=T)) %>% group_by(Alt,Year) %>% summarise(dif_mean_ref=mean(dif_mean_ref)) %>%  
    ggplot(aes(x=Year,y=Alt,fill=dif_mean_ref)) +
    geom_tile() +
    scale_fill_stepsn(colours=c(MetBrewer::met.brewer("Hiroshige",direction=1,n=20)),n.breaks=20)
  
  dos <- to_plot %>% mutate(Year=((YearCat))) %>% 
    mutate(Alt=cut(y,breaks=seq(0,4600,200),labels=seq(0,4600,200)[-1], include.lowest=T)) %>% 
    group_by(Alt,Year) %>% summarise(y=mean(y,na.rm=T),median_ERROR=quantile(mean_ERROR,probs=0.5,na.rm=T)) %>% 
    ggplot(aes(x=Year,y=Alt,fill=median_ERROR)) +
    geom_tile() +
    scale_fill_stepsn(colours=c(MetBrewer::met.brewer("Hiroshige",direction=-1)),n.breaks=10)
  
  to_plot %>% mutate(Year=((YearCat))) %>% 
    mutate(Alt=cut(y,breaks=seq(0,4600,200),labels=seq(0,4600,200)[-1], include.lowest=T)) %>% 
    group_by(Year,Alt) %>% summarise(y=mean(y,na.rm=T),median_ERROR=quantile(mean_ERROR,probs=0.5,na.rm=T),dif_mean_ref=mean(dif_mean_ref)) %>% ggplot(aes(x=abs(dif_mean_ref),y = median_ERROR,color=Year)) +
    geom_point() +
    scale_color_stepsn(colours=c(MetBrewer::met.brewer("Hiroshige",direction=-1)),n.breaks=10) +
    theme(panel.background = element_blank(),panel.grid=element_blank(),
          axis.line=element_line()) +
    labs(x="Altitudinal shift (relative to reference)",y="Standard error of the mean") +
    NULL
  
  to_plot %>% mutate(Year=((YearCat))) %>% 
    mutate(Alt=cut(y,breaks=seq(0,4600,200),labels=seq(0,4600,200)[-1], include.lowest=T)) %>% 
    group_by(Year,Alt) %>% summarise(y=mean(y,na.rm=T),median_ERROR=quantile(mean_ERROR,probs=0.5,na.rm=T),dif_mean_ref=mean(dif_mean_ref)) %>% ggplot(aes(y = median_ERROR)) +
    geom_histogram(fill=MetBrewer::met.brewer("Hiroshige",direction=1)[[10]]) +
    theme(panel.background = element_blank(),panel.grid=element_blank(),
          axis.line=element_line(),legend.position="") +
    labs(x="Altitudinal shift (relative to reference)",y="Standard error of the mean") +
    NULL
  
  
  library(patchwork)
  layout <- "
AA
BB
"
  wrap_plots(list(uno,dos)) +
    plot_layout(design = layout,byrow=T,nrow = 2,ncol=2,guides = "collect",tag_level = 'new') + 
    plot_annotation(title = paste0("Species' elevational trends (",cual,"-point)"),subtitle="Turning point 1975-1976", caption= "Ramírez-Barahona et al.")
  
  
  to_plot %>% filter(YearCat %in% c(1979:1985)) %>% 
    mutate(Year=factor((YearCat))) %>% 
    mutate(Alt=cut(y,breaks=seq(0,4500,500),labels=seq(0,4500,500)[-1], include.lowest=T)) %>% summarise(Alt_Cat = factor(Mode(Alt))) -> altas
  
  to_plot %>% mutate(Year=factor((YearCat))) %>% left_join(.,altas,by="Accepted_Name") %>% filter(!is.na(Alt_Cat)) %>% group_by(YearCat,Alt_Cat) %>% 
    summarise(All_Meanalt=mean(RollAlt,na.rm=TRUE),All_SDalt= sd(RollAlt,na.rm=TRUE)) %>% filter(YearCat<2016) %>% 
    ggplot(aes(x=YearCat,y=All_SDalt,color=fct_inseq(Alt_Cat))) + 
    geom_point() +
    geom_line() +
    #geom_smooth(method="loess") +
    scale_color_manual(values = MetBrewer::met.brewer("Tam"))
  
  species <- to_plot %>% distinct(Accepted_Name) %>% pull() %>% as.vector()
  centroids <- reduct_five %>% group_by(Accepted_Name) %>% summarise(across(decimalLatitude:decimalLongitude,mean))
  
  to_plot %>% left_join(.,centroids,by="Accepted_Name") %>% mutate(Year=factor((YearCat))) %>% left_join(.,altas,by="Accepted_Name") %>% filter(!is.na(Alt_Cat)) %>% mutate(LAT=cut(decimalLatitude,breaks=seq(5,30,1),labels=seq(5,30,1)[-1], include.lowest=T)) %>% mutate(LON=cut(decimalLongitude,breaks=seq(-110,-75,2),labels=seq(-110,-75,2)[-1], include.lowest=T)) %>% group_by(YearCat,LAT) %>% summarise(All_Meanalt=mean(RollAlt,na.rm=TRUE),All_SDalt= sd(RollAlt,na.rm=TRUE)) %>% 
    filter(YearCat<2016) %>% 
    ggplot(aes(x=YearCat,y=All_SDalt,color=LAT)) + 
    geom_point() +
    geom_line() + 
    #geom_smooth(method="loess",se = FALSE) +
    scale_color_manual(values = MetBrewer::met.brewer("Tam",n=20)) +
    NULL
  
  
  data="output/Historical_data_FULL_LUL.v.1.R"
  load(here(data))
  reduct_five %>% filter(Accepted_Name==species[2]) %>% group_by(year) %>% 
    filter(year %in% c("Year_1986","Year_1987","Year_1983","Year_1984","Year_1985")) %>% 
    ggplot(aes(x = Accepted_Name, y = Altitude,fill=Accepted_Name)) + 
    ggdist::stat_halfeye(adjust = .9, width = .6, .width = 0, alpha=0.7, 
                         slab_colour = "black", slab_size = .2, justification = -.4, height = 2) + 
    geom_boxplot(width = .4, outlier.shape = NA,size=.2,fatten = 7,alpha=0.7) +
    geom_point(fill="grey50", size = 2, alpha = .4, shape = 21, stroke = 0.2, 
               position = position_jitter(seed = 1, width = .2)) 
  
  reduct_five %>% select(CellID:Accepted_Name,ID,class:family.x,decimalLatitude:BIOME) %>% 
    mutate(YearCat = as.numeric(sub("Year_","",year)),.after=1) %>% 
    filter(!is.na(Altitude),YearCat >= start_year) %>% distinct(Accepted_Name,CellID,year,.keep_all = T) %>% group_by(Accepted_Name) -> esta
  esta %>% summarise(N=n()) -> tata
  target <- tata %>% filter(N >= records) %>% pull(Accepted_Name) %>% as.character(.)
  
  esta %>% #filter(Accepted_Name %in% all_of(target)) %>% 
    filter(Accepted_Name%in%c(species[2:3])) %>% droplevels() %>% group_by(Accepted_Name,YearCat) %>% 
    summarise(Variable=mean(Altitude,na.rm=T), MaxVariable = quantile(Altitude,probs=0.975,na.rm=T),MinVariable=quantile(Altitude,probs=0.025,na.rm=T),n=n(),class=first(class)) %>% ungroup() %>% 
    mutate_if(is.numeric, ~na_if(., -Inf)) %>% 
    mutate_if(is.numeric, ~na_if(., Inf)) %>% arrange(Accepted_Name,YearCat) -> esta_2
  
  esta_2 %>% rowwise() %>% mutate(Roll_dist=(unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width = window, FUN = is.unimodal, partial=T, align="center")}),use.names = F))) -> aver
}
{resamp <- readRDS("output/SES_shifting_v.3.RDS")
  resamp %<>% mutate(Samp_SES = cut(estimate,breaks=c(min(estimate,na.rm=TRUE),-1.96,-1.28,0,1.28,1.96,max(estimate,na.rm=TRUE)),labels=c("Strong_Neg","Moderate_Neg","Weak_Neg","Weak_Pos","Moderate_Pos","Strong_Pos"),include.lowest=TRUE))
  down <- readRDS("output/SES_shiftingSkew_v.3.RDS")
  down %>% distinct(ind,.keep_all = TRUE) %>% summarise(mean(obs_mid),seme(obs_mid))
  
  down %<>% mutate(Samp_SES = cut(estimate,breaks=c(min(estimate,na.rm=TRUE),-1.96,-1.28,0,1.28,1.96,max(estimate,na.rm=TRUE)),labels=c("Strong_Neg","Moderate_Neg","Weak_Neg","Weak_Pos","Moderate_Pos","Strong_Pos"),include.lowest=TRUE))
  resamp %>% filter(Samp_SES %in% c("Weak_Pos","Weak_Neg")) %>% 
    filter(SES=="SES_max") -> a
  down %>% filter(Samp_SES %in% c("Weak_Pos","Weak_Neg")) %>% 
    filter(SES=="SES_max") -> b
  nrow(a);nrow(b)}
{
  data="output/Historical_data_FULL_LUL.v.1.R"
  data_spp = "output/Species_summaryRaw.RData"
  load(here(data))  ### occurrences
  load(data_spp)    ### summary (vs. nulls)
  species_means[[1]] %>% count(Accepted_Name) %>% pull(Accepted_Name) -> targets
  reduct_five %<>% filter(Accepted_Name%in% targets)
  g <- raster::raster(nrows=180*24,ncols=360*24,xmn=-180,xmx=180,ymn=-90,ymx=90,vals=1) %>% 
    projectRaster(.,crs = CRS("+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m  +no_defs")) %>% as(., 'SpatialPixels')
  
  reduct_five %<>% 
    mutate(year=sub("Year_","",year)) %>% mutate(year = as.numeric(year)) %>% 
    mutate(Decade = cut(year,breaks=c(1979,1990,2000,2010,2020),include.lowest=TRUE,labels = c("80s","90s","00s","10s") ),.before=Altitude)
  
  reduct_five %>% distinct(CellID,Accepted_Name,.keep_all = TRUE) %>% 
    ggplot(aes(x=Decade, y=Altitude)) + geom_boxplot()
  
  reduct_five %>% distinct(Accepted_Name,CellID,Decade,.keep_all = TRUE) %>% group_by(Decade) %>% count(CellID) %>%
    pivot_wider(names_from = Decade,values_from = n,names_prefix = "SR_") %>% left_join(.,reduct_five %>% distinct(CellID,.keep_all = TRUE),by="CellID") -> SR
  
  g@coords %>% as_tibble() %>% mutate(CellID=1:nrow(.)) %>% right_join(.,SR,by="CellID") -> SR
  roads <- rnaturalearth::ne_countries(scale = 110,returnclass = "sf") %>% sf::st_transform(., crs = CRS("+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m  +no_defs"))
  xlimits = c(-11286442,-6932762)
  ylimits = c(988516,3849386)
  SR %>% mutate(n_std = (n - min(n)) / (max(n) - min(n))) -> SR
  reduct_five %>% count(CellID) %>% left_join(.,reduct_five %>% distinct(CellID,.keep_all = T),by="CellID") -> counts
  g@coords %>% as_tibble() %>% mutate(CellID=1:nrow(.)) %>% right_join(.,counts,by="CellID") -> counts
  counts %>% mutate(n_std = (n - min(n)) / (max(n) - min(n))) -> counts
  
  ggplot() + 
    geom_sf(data=roads,colour=NA,fill=adjustcolor(MetBrewer::met.brewer("Demuth",direction=-1)[1],alpha.f = 0.2),size=0.5) + xlim(xlimits) + ylim(ylimits) + 
    theme(panel.background = element_blank(),panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.border = element_rect(fill=NA),legend.position = c(0.9,0.35),legend.background = element_rect(fill = NA), legend.key.size = unit(0.7, "cm"),legend.key.width = unit(0.6,"cm"), legend.text = element_text(family="EB Garamond",size=8),legend.title = element_text(family="EB Garamond"),legend.spacing = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5)) +
    geom_tile(data = SR %>% filter(!is.na(SR_00s)),aes(x=x,y=y,fill=SR_00s,color=SR_00s)) +
    scale_fill_stepsn(colors=(MetBrewer::met.brewer("Tam",direction=1)),name="Species  \nrichness (log)",n.breaks=10,trans = "log",labels=function(x) sprintf("%.0f", (x))) +
    scale_color_stepsn(colors=(MetBrewer::met.brewer("Tam",direction=1)),name="Species  \nrichness (log)",n.breaks=10,trans = "log",labels=function(x) sprintf("%.0f", (x))) +
    labs(x="",y="",title="Occurrence records for cloud forests species", subtitle="SR_00s") +
    #annotate("text",x=-11132762,y=1288516,label=paste0("More than 50 occurrence records: ",formatC(nrow(plus50), big.mark=",")," species"," (",formatC(sum(plus50$n), big.mark=","),")"),family="EB Garamond",size=4,hjust=0) +
    #annotate("text",x=-11132762,y=1188516,label=paste0("More than 100 occurrence records: ",formatC(nrow(plus100), big.mark=",")," species"," (",formatC(sum(plus100$n), big.mark=","),")"),family="EB Garamond",size=4,hjust=0) +
    #annotate("text",x=-11132762,y=1088516,label=paste0("More than 500 occurrence records: ",formatC(nrow(plus500), big.mark=",")," species"," (",formatC(sum(plus500$n), big.mark=","),")"),family="EB Garamond",size=4,hjust=0) +
    #standardPrintOutput::watermark(show=T,lab="Ramirez-Barahona et al.") +
    annotation_custom(ggplotGrob(SR %>% filter(!is.na(SR_00s)) %>% 
                                   ggplot(aes(x=Altitude)) + 
                                   geom_histogram(stat="density",bins = 50,fill=MetBrewer::met.brewer("Hiroshige",direction=1)[10]) + theme(panel.background = element_blank(),panel.grid = element_blank(),axis.line = element_line(),legend.position = c(0.1,0.3),legend.background = element_rect(fill = NA), legend.key.size = unit(0.7, "cm"),legend.key.width = unit(0.6,"cm"), legend.text = element_text(family="EB Garamond",size=8),legend.title = element_text(family="EB Garamond"),legend.spacing = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=14,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),axis.text = element_text(family="EB Garamond"),axis.title = element_text(family="EB Garamond",size=11)
                                   ) + labs(y="No. cells",x="Elevation") ), xmin = -11432762, xmax = -9632762, 
                      ymin = 1388516, ymax = 2388516) +
    NULL
  
  SR %>% summarise(across(starts_with("SR_"),function(x) mean(x,na.rm=TRUE) ))
  SR %>% select(CellID,starts_with("SR_"),Altitude) %>% pivot_longer(-c(CellID,Altitude),names_to = "Decade",values_to = "SR") %>% mutate(Altitude=ifelse(is.na(SR),NA,Altitude)) %>% 
    summarise(mean(Altitude,na.rm=TRUE),.by=Decade)
  
  
  
  
  mean_alts
  mean_alts %<>% select(data) %>% unnest(data) %>% 
    mutate(Decade = cut(YearCat,breaks=c(1979,1990,2000,2010,2020),include.lowest=TRUE,labels = c("80s","90s","00s","10s") ),.before=YearCat) %>% mutate(tme = YearCat - 1979,.before=YearCat)
  
  mean_alts %<>% filter(Decade != "10s")
  lmer(RollAlt ~ 1 + tme + tme:RollAlt_sd + (1|Accepted_Name) + (0 + tme|Accepted_Name),data = mean_alts, control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4))) -> fit1
  lm(RollAlt ~ 1 + tme + tme:RollAlt_sd,data = mean_alts) %>% broom.mixed::tidy()
  
  ggeffects::hypothesis_test(fit1,terms="tme")
  ggeffects::ggpredict(fit1, type="random",terms="tme",interval="confidence")
  mean_alts %>%
    ggplot(aes(y=as.character(Decade), x=RollAlt,fill=after_stat(quantile))) +
    stat_density_ridges(quantile_lines = FALSE,
                        calc_ecdf = TRUE,
                        geom = "density_ridges_gradient",
                        quantiles = c(0.1,0.5, 0.9),
                        rel_min_height=0.05) +
    scale_fill_brewer(name = "") +
    coord_cartesian(xlim=c(0,4000))
  
  
  mean_alts %>% count(Accepted_Name) %>% pull(Accepted_Name) -> targets
  reduct_five %<>% filter(Accepted_Name%in% targets)
  reduct_five %<>%  mutate(year=sub("Year_","",year)) %>% mutate(year = as.numeric(year)) %>% 
    filter(year <= 2010)
  reduct_five %>% 
    ggplot(aes(y=as.character(year), x=Altitude,fill=stat(quantile))) +
    stat_density_ridges(quantile_lines = FALSE,
                        calc_ecdf = TRUE,
                        geom = "density_ridges_gradient",
                        quantiles = c(0.1,0.5, 0.9),
                        rel_min_height=0.05) +
    scale_fill_brewer(name = "") +
    coord_cartesian(xlim=c(0,4000))
}
{
  reduct_five %>% filter(Accepted_Name %in% all_of(target_spp)) %>% select(-Altitude,-Accepted_Name,-c(In_Villa:BIOME)) %>% 
    distinct(CellID,.keep_all = TRUE) %>% 
    pivot_longer(-CellID, names_to = "var",values_to = "val") -> mero
  mero$val[which(is.infinite(mero$val))] <- NA
  mero_mans <- mero %>% mutate(year=stringi::stri_sub(var,-4)) %>% mutate(var=substring(var,1,nchar(var)-5)) %>% drop_na() %>%
    group_by(CellID) %>% nest() %>% rowwise() %>% mutate(manken_tas = manken(data,tag="_pr")) 
  #mero_mans %<>%  mutate(manken_cmi = manken(data,tag="_cmi"))
  #mero_mans %<>% mutate(manken_pr = manken(data,tag="_pr"))
  #mero_mans %<>% mutate(manken_pet = manken(data,tag="_pet"))
  
  mero <- reduct_five %>% filter(Accepted_Name %in% all_of(target_spp)) %>% select(CellID,Accepted_Name,Altitude,decimalLongitude,decimalLatitude) %>% distinct(CellID,Accepted_Name,.keep_all = T) %>% left_join(.,mero_mans,by="CellID")
  
  mero
  aver <- reduct_five %>% filter(Accepted_Name%in% all_of(target_spp)) %>% 
    distinct(Accepted_Name,CellID,year,.keep_all = TRUE) %>%
    group_by(Accepted_Name) %>% select(CellID,year,Altitude,starts_with("d")) %>% 
    mutate(year=sub("Year_","",year)) %>% rename(YearCat=year) %>%
    ungroup() %>% select(-Accepted_Name) %>%
    distinct(CellID,.keep_all = TRUE) %>% 
    select(-decimalLatitude,-decimalLongitude) %>% 
    pivot_longer(cols = -c(1:3),names_to = "year",values_to = "Values") %>%  
    #mutate(year=sub("_","",year)) %>% 
    separate(year,sep="_",into=c("Var","year")) %>% arrange(CellID,year) %>% group_by(CellID) %>% select(-YearCat)
  
  reduct_five %>% select(CellID:Accepted_Name,year) %>% 
    #distinct(CellID,Accepted_Name,year) %>% 
    filter(Accepted_Name %in% all_of(target_spp)) %>% 
    mutate(year=sub("Year_","",year)) %>% 
    left_join(., aver,by=c("CellID","year")) %>% 
    rename(Altitude = Altitude.x) %>% 
    #filter(!is.na(Corrected_values)) %>% 
    group_by(Accepted_Name) %>%  
    mutate(year = as.numeric(sub("Year_","",year))) -> full_siono
  
  mero %>% select(CellID,starts_with("manken")) %>% left_join(full_siono,.,by="CellID") -> full_siono
  
  
  full_siono %>%
    rename(Mean_Temp = Values,Elevation=Altitude) %>% 
    select(Accepted_Name,year,Mean_Temp,Elevation) %>% 
    filter(!is.na(Mean_Temp)) %>%
    mutate(Mean_Temp = Mean_Temp*.1 - 273.15) %>% 
    rename(ind=Accepted_Name) %>% mutate(tme = year - min(year)) %>% 
    mutate(Elevation_Cat = cut(Elevation,breaks=seq(0,4000,1000),labels=seq(0,4000,1000)[-1], include.lowest=TRUE)) %>%
    separate(ind,into = c("G","S","extra"),extra = "merge",sep=" ") %>% select(-extra) %>% unite("ind",G:S,sep="_") %>% 
    left_join(., nulas,by="ind") %>%
    mutate(RSamp_SES_cat = cut(RSamp_SES,breaks=c(-10,-1.96,-1.28,0,1.28,1.96,10),labels=c("Negative","Negative","Weak","Weak","Positive","Positive"),include.lowest=TRUE)) %>% 
    mutate(DSamp_SES_cat = cut(DSamp_SES,breaks=c(-10,-1.96,-1.28,0,1.28,1.96,10),labels=c("Sensitive","Sensitive","Robust","Robust","Sensitive","Sensitive"),include.lowest=TRUE)) -> aqui # %>% #filter(DSamp_SES_cat %in% c("Robust")) 
  
  
  #otra %>% select(ind,y,year,error) %>% right_join(.,aqui,by=c("ind","year")) %>% 
  aqui %>% #mutate(Mean_Temp = scale(Mean_Temp)) %>% 
    filter(DSamp_SES_cat %in% c("Robust")) %>% 
    #glm(Land_use ~ 1 + tme,family =  binomial, data = .) -> fit_lul
    lmer(Mean_Temp ~ 1 + tme + tme:RSamp_SES_cat + (1|ind), data = .,
         control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4))) -> fit_lul
  broom.mixed::tidy(fit_lul)
  coef(fit_lul)$ind %>% dim()
  fit_lul_predict <- ggpredict(fit_lul, type = "random",terms=c("tme [all]","RSamp_SES_cat"),interval = "confidence") %>% tibble()
  #ggpredict(fit_lul,c("tme [all]","RSamp_SES_cat")) %>% hypothesis_test() %>% tibble() -> test_lul
  #test_lul
  fit_lul_predict %>% 
    ggplot(aes(x=x+1979,y=predicted,color=group,fill=group)) +
    geom_line() +
    geom_ribbon(aes(ymin = conf.low,ymax=conf.high),alpha=0.1,color="NA") +
    scale_color_met_d(name="Hiroshige",direction = -1)+
    scale_fill_met_d(name="Hiroshige",direction = -1) +
    theme(legend.position = "bottom",panel.background = element_blank(),axis.line = element_line(),legend.background = element_blank(),
          axis.title = element_text(family="EB Garamond"),
          axis.text = element_text(family="EB Garamond"),
          #axis.text.x = element_blank(),
          legend.spacing.y = unit(c(1.2),"cm"),
          legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond")) +
    labs(x="Time",y="Predicted probabilities (Forest land-cover class)") +
    NULL
  
  
  
  
  
  
  
  
  coef(fit_lul)$ind %>% as_tibble(rownames = "ind") %>% rename_with(~c("ind","Intercept_LUC","time_LUC","Interaction")) %>% left_join(.,nulas,by="ind") %>% 
    left_join(.,fit_lul_predict %>% filter(x==0) %>% rename(ind=group,year=x),by="ind") %>% 
    left_join(., esta %>% rename(time_elevation = estimate),by="ind") %>% 
    mutate(RSamp_SES_cat = cut(RSamp_SES,breaks=c(-6,-1.96,-1.28,0,1.28,1.96,6),labels=c("Strong_Neg","Moderate_Neg","Weak_Neg","Weak_Pos","Moderate_Pos","Strong_Pos"),include.lowest=TRUE)) %>% 
    mutate(DSamp_SES_cat = cut(DSamp_SES,breaks=c(-6,-1.96,-1.28,0,1.28,1.96,6),labels=c("Strong_Neg","Moderate_Neg","Weak_Neg","Weak_Pos","Moderate_Pos","Strong_Pos"),include.lowest=TRUE)) %>% 
    filter(!DSamp_SES_cat %in% c("Strong_Pos","Moderate_Pos")) %>% 
    ggplot(aes(x=RSamp_SES_cat, y = time_LUC)) +
    #geom_point() +
    geom_boxplot() +
    ggdist::geom_dotsinterval()+
    #geom_smooth(method="lm") +
    theme(legend.position = "none")+
    NULL
  
  
  fit_lul_predict %>% select(1:2,6) %>% pivot_wider(names_from = x,values_from = predicted,names_prefix = "Year_") %>% 
    pivot_longer(-1,names_to = "year",values_to = "prob") %>% rename(ind=group) %>% 
    mutate(year = factor(year,levels=paste0("Year_",0:33))) %>% 
    ggplot(aes(x=year,y=prob,color=ind,group=ind)) +
    geom_line() +
    theme(legend.position = "none")
  
  
  
  
  
  
  
  fit_lul_predict <- ggpredict(fit_lul, type = "random",terms=c("tme [all]"),interval = "confidence") %>% tibble()
  ggpredict(fit_lul,c("tme [all]")) %>% hypothesis_test() %>% tibble() -> test_lul
  fit_grid_predict <- ggpredict(fit_grid, type = "random",terms=c("tme [all]"),interval = "confidence") %>% tibble()
  ggpredict(fit_grid,c("tme [all]")) %>% hypothesis_test() %>% tibble() -> test_grid
  
  fit_lul_predict %>% 
    ggplot(aes(x=x+1979,y=predicted)) +
    geom_line() + 
    geom_ribbon(aes(ymin = conf.low,ymax=conf.high),alpha=0.1) +
    geom_line(data=fit_grid_predict,color="red") + 
    geom_ribbon(data=fit_grid_predict,aes(ymin = conf.low,ymax=conf.high),fill="red",alpha=0.1) +
    ylim(c(0,1)) +
    theme(legend.position = "none",panel.background = element_blank(),axis.line = element_line(),legend.background = element_blank(),legend.key.height = unit(0.01,"cm"),legend.key.width = unit(0.01,"cm"),
          axis.title = element_text(family="EB Garamond"),
          axis.text = element_text(family="EB Garamond"),
          #axis.text.x = element_blank(),
          legend.spacing.y = unit(c(1.2),"cm"),
          legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond")) +
    labs(x="Time",y="Predicted probabilities (Forest land-cover class)") +
    #equatiomatic::extract_eq(fit_lulPF) +
    #geom_text(x=1990, y= 0.6,label=paste(round(test_grid$Slope,3)," \U0028",round(test_grid$conf.low,3),", ",round(test_grid$conf.high,3),"\U0029",sep=""),color="red") +
    #geom_text(x=2000, y= 0.9,label=paste(round(test_lul$Slope,3)," \U0028",round(test_lul$conf.low,3),", ",round(test_lul$conf.high,3),"\U0029",sep=""),color="black")  +
    NULL
  
  
}
{
  reduct_five %>% select(CellID,year, contains("min_pr_")) %>% 
    mutate(year=sub("Year_","",year)) %>% rename(YearCat=year) %>%
    distinct(CellID,.keep_all = TRUE) %>% 
    pivot_longer(cols = -c(1:2),names_to = "year",values_to = "Values") %>%  
    separate(year,sep="_",into=c("stat","Var","year")) %>% 
    arrange(CellID,year) %>% group_by(CellID) %>% select(-YearCat) -> aver_clim
  
  siono2 <- aver_clim %>% 
    # mutate(Temperature = Values*.1 - 273.15) %>% 
    mutate(PP = Values) %>% 
    mutate(year = as.numeric(year)) %>% 
    mutate(tme = year - min(year))
  #mutate(CMI = as.vector(scale(CMI)))
  fit_cmi <- lmer(PP ~ 1 + tme + (1|CellID) +(0 + tme|CellID), data = siono2,control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4)))
  coef(fit_temp)$CellID %>% as_tibble(rownames = "CellID") %>% rename(Temp = `(Intercept)`, Change_Temp = tme) %>% 
    right_join(.,coef(fit_cmi)$CellID %>% as_tibble(rownames = "CellID") %>% rename(PP = `(Intercept)`, Change_PP = tme) ,by="CellID") %>% 
    ggplot(aes(x=Change_Temp,y=Change_PP)) +
    geom_hex(bins=90) +
    scale_fill_met_c(name="Hiroshige")
  
  
  testa_clim <- reduct_five %>% select(CellID:Accepted_Name,year) %>% 
    #distinct(CellID,Accepted_Name,year) %>% 
    filter(Accepted_Name %in% all_of(target_spp)) %>% 
    mutate(year=sub("Year_","",year)) %>% 
    left_join(., aver_clim,by=c("CellID")) %>% 
    filter(!is.na(Values)) %>% 
    group_by(Accepted_Name) %>% nest() %>% rowwise() %>% mutate(data = list({data %>% 
        mutate(year = as.numeric(sub("Year_","",year)))
    })) -> full_clim
  
  full_clim %>% unnest(data) %>% 
    mutate(Temperature = Values*.1 - 273.15) %>% 
    select(Accepted_Name:year,Temperature) %>% 
    filter(!is.na(Temperature)) %>%
    rename(ind=Accepted_Name) %>% 
    separate(ind,into = c("G","S","extra"),extra = "merge",sep=" ") %>% select(-extra) %>% unite("ind",G:S,sep="_") %>% ungroup() %>% left_join(.,otra,by=c("ind","year")) -> siono2
  siono2 %<>% filter(!is.na(tme)) %>% mutate(tme = tme - min(tme))
  
  coef(fit_temp)$ind %>% as_tibble(rownames = "ind") %>% rename(beta_temp = tme, beta0_temp = `(Intercept)`) %>% left_join(.,coef(fit_nested)$ind %>% as_tibble(rownames = "ind") %>%  rename(beta_time = tme, beta0_time = `(Intercept)`),by="ind") %>%
    left_join(.,otra %>% distinct(ind,.keep_all = TRUE),by="ind") %>% 
    ggplot(aes(y=beta_time,x=beta_temp)) + 
    geom_point() +
    geom_smooth(method = "lm") +
    NULL
  
  sumfun = function(.) {coef(.)$ind[,"tme"]}
  ## where the fun happens
  booty <- bootMer(fit_nested, sumfun, nsim = 1000, re.form= NULL)
  tibble(mean=apply(booty$t, 2, mean),sd=apply(booty$t, 2, sd),coef(fit_nested)$ind) %>% 
    mutate(bias=tme-mean) %>% ggplot(aes(x=bias)) + geom_histogram()
  
  booty_lul <- bootMer(fit_lulPF, sumfun, nsim = 1000, re.form= NULL)
  tibble(mean=apply(booty$t, 2, mean),sd=apply(booty$t, 2, sd),coef(fit_nested)$ind) %>% 
    mutate(bias=tme-mean) %>% ggplot(aes(x=bias)) + geom_histogram()
  
  
  
  
  ras <- raster("data/climate_data_CHELSA/Velocidad/dv50_tmean.tif")
  plot(ras) 
  
  #sims <- DHARMa::simulateResiduals(fittedModel = fit_lulPF, n = 100,re.form=NULL)
  #plot(y=sims$scaledResiduals,x=fitted(fit_lulPF))
  
  if(series=="LCCS") {
    reduct_five %>%  
      mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(10,11,12,20,30),"Human",.x))) %>% 
      mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(50,60:62,70:72,80:82,90,160,170),"Forest",.x))) %>%       mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(110,120:122,130,140,180,40,100),"Primary_Vegetation",.x))) %>% 
      mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(190),"Human",.x))) %>% 
      mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(150:153,200:202),"Bare land",.x))) %>% 
      mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(210),"Water bodies",.x))) %>% 
      mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(220),"IceSnow",.x))) %>% 
      group_by(Accepted_Name) %>% select(CellID,year,starts_with("LUC_")) %>% 
      mutate(year=sub("Year_","",year)) %>% rename(YearCat=year) %>% ungroup() %>% select(-Accepted_Name) %>% 
      distinct(CellID,.keep_all = TRUE) %>% pivot_longer(cols = -c(1:2),names_to = "year",values_to = "Values") %>%         separate(year,sep="_",into=c("Var","year")) %>% arrange(CellID,year) %>% 
      group_by(CellID) %>% select(-YearCat) -> aver
    
    aver %>% group_by(CellID,Values) %>% count()
    aver %>% mutate(New_values = ifelse(Values=="Vegetaion" & lag(Values,n=1) == "Human","Secondary",Values)) %>% 
      mutate(New_values = ifelse(Values=="Forest" & lag(Values,n=1) == "Human","Secondary",Values)) -> aver
    
    aver %>% nest() -> uno
    
    aver %>% nest() %>% mutate(New_data = map(data,~ lagas))
                               lagas(uno$data[[1]]) -> x
                               lagas <- function(x) {
                                 for(i in 2:nrow(x)){ 
                                   if(is.na(x$New_values[i])) next
                                   val <- ifelse(x$New_values[i-1] == "Secondary","Secondary",x$New_values[i])
                                   if(is.na(val)) next
                                   x$New_values[[i]] <- val
                                   return(x)
                                 }
                               }   
                               
                               
                               
                               
category = "Human"
aver %<>% group_by(CellID) %>% mutate(Rolllul_n = unlist(tapply(Values, CellID, function(x) {(zoo::rollapply(as_vector(Values),width = window, FUN = parse_count, partial=T, align="center"))},simplify = TRUE)),Roll_Human = unlist(tapply(Values, CellID, function(x) {(zoo::rollapply(as_vector(Values),width = window, FUN = parse_freq, partial=T, align="center"))},simplify = TRUE)))
                               
category = "Forest"
aver %<>% group_by(CellID) %>% mutate(Roll_Forest = unlist(tapply(Values, CellID, function(x) {(zoo::rollapply(as_vector(Values),width = window, FUN = parse_freq, partial=T, align="center"))},simplify = TRUE)))
                               
                               category = "Vegetation"
                               aver %<>% group_by(CellID) %>% mutate(Roll_Vegetation = unlist(tapply(Values, CellID, function(x) {(zoo::rollapply(as_vector(Values),width = window, FUN = parse_freq, partial=T, align="center"))},simplify = TRUE)))
                               
                               category = "Other"
                               aver %<>% group_by(CellID) %>% mutate(Roll_Other = unlist(tapply(Values, CellID, function(x) {(zoo::rollapply(as_vector(Values),width = window, FUN = parse_freq, partial=T, align="center"))},simplify = TRUE)))
                               #save(aver,file = here("interim",output))
                               
                               reduct_five %>% select(CellID:Accepted_Name,year) %>% 
                                 filter(Accepted_Name %in% all_of(target_spp)) %>% 
                                 mutate(year=sub("Year_","",year)) %>% 
                                 left_join(.,aver,by=c("CellID","year")) %>% filter(!is.na(Values)) %>% 
                                 group_by(Accepted_Name, year) %>% 
                                 group_by(year) %>% 
                                 summarise(across(starts_with("Roll"), ~ mean(.x))) %>%
                                 ggplot(aes(x=year,y=Roll_Vegetation/Rolllul_n)) + 
                                 geom_boxplot() +
                                 NULL
                               
                               load(here("output/Historical_data_FULL_LCCS.v.1.R"))
                               reduct_five %>%
                                 mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(10,11,12,20,30),"Human",.x))) %>% 
                                 mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(50,60:62,70:72,80:82,90,100,160,170),"Vegetation",.x))) %>% 
                                 mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(110,120:122,130,140,180,40),"Vegetation",.x))) %>% 
                                 mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(190),"Human",.x))) %>% 
                                 mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(150:153,200:202),"Bare land",.x))) %>% 
                                 mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(210),"Water bodies",.x))) %>% 
                                 mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(220),"IceSnow",.x))) %>% select(Accepted_Name,CellID,Altitude,year,starts_with("LUC")) %>% filter(Accepted_Name %in% all_of(target_spp)) %>%
                                 group_by(Accepted_Name) %>% nest() -> testa
                               

  library(fGarch)
  quantile(rsnorm(10000,xi=1,mean = dist_r1$mean,sd = 1),probs=c(0.05,0.5,0.95))
  quantile(rsnorm(10000,xi=-200,mean = dist_r1$mean,sd = 1),probs=c(0.05,0.5,0.95)) 
  quantile(rsnorm(10000,xi=0.8,mean = dist_r1$mean,sd = 1),probs=c(0.05,0.5,0.95))
  quantile(rsnorm(10000,xi=0.4,mean = dist_r1$mean,sd = 1),probs=c(0.05,0.5,0.95))
  quantile(rsnorm(10000,xi=0.1,mean = dist_r1$mean,sd = 1),probs=c(0.05,0.5,0.95))
  quantile(rsnorm(10000,xi=0.01,mean = dist_r1$mean,sd = 1),probs=c(0.05,0.5,0.95))
  
  hist(rsnorm(10000,xi=1,mean = dist_r1$mean,sd = 1),probs=c(0.05,0.5,0.95))
  hist(rsnorm(10000,xi=0.9,mean = dist_r1$mean,sd = 1),probs=c(0.05,0.5,0.95)) 
  hist(rsnorm(10000,xi=0.8,mean = dist_r1$mean,sd = 1),probs=c(0.05,0.5,0.95))
  hist(rsnorm(10000,xi=0.4,mean = dist_r1$mean,sd = 1),probs=c(0.05,0.5,0.95))
  hist(rsnorm(10000,xi=0.1,mean = dist_r1$mean,sd = 1),probs=c(0.05,0.5,0.95))
  hist(rsnorm(10000,xi=0.01,mean = dist_r1$mean,sd = 1),probs=c(0.05,0.5,0.95))
}}
{ 
  #spp_base_alts <- reduct_five %>% mutate(Year_num = as.numeric(sub("Year_","",year))) %>% filter(Year_num <= 1990) %>% 
  #summarise(across(starts_with("mean_tas_198"),~ mean(.x)),.by=Accepted_Name) %>% rowwise() %>% 
  #mutate(Base_temp = mean(c_across(mean_tas_1980:mean_tas_1989))) %>% select(Accepted_Name,Base_temp) %>% 
  #ungroup() %>% 
  #mutate(Temp_cat = cut(Base_temp,breaks=seq(2500,3500,50),labels=seq(2500,3450,50))) %>% rename(ind=Accepted_Name)
  #  summarise(Base_alt = mean(Altitude,na.rm=TRUE),.by=Accepted_Name) %>% 
  #mutate(Alt_cat = cut(Base_alt,breaks=seq(0,4600,200),labels=seq(0,4400,200))) %>% rename(ind=Accepted_Name)
  
  
  
  
  
  
  
  aver <- reduct_five %>% filter(Accepted_Name %in% target_spp) %>% 
    distinct(Accepted_Name,CellID,year,.keep_all = T) %>% 
    mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(10),"Agriculture",.x))) %>%
    mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(20),"Forest",.x))) %>%
    mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(30,40),"Grassland_Shrub",.x))) %>%
    mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(70,100,90),"Other",.x))) %>%
    mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(0),NA,.x))) %>% 
    group_by(Accepted_Name) %>% select(CellID,year,Altitude,starts_with("LUC_")) %>% 
    mutate(year=sub("Year_","",year)) %>% rename(YearCat=year) %>%
    ungroup() %>% select(-Accepted_Name) %>%
    distinct(CellID,.keep_all = TRUE) %>% 
    pivot_longer(cols = -c(1:3),names_to = "year",values_to = "Values") %>%  
    separate(year,sep="_",into=c("Var","year")) %>% arrange(CellID,year) %>% group_by(CellID) %>% select(-YearCat)
  aver
  
  aver %>% mutate(New_values = ifelse(Values == "Forest" & lag(Values,n=1) != "Forest","Altered_Forest",Values)) %>% mutate(New_values = ifelse(is.na(New_values),Values,New_values)) -> uno
  uno %<>% mutate(Corrected_values = New_values)
  uno %<>% nest()
  for(j in 1:nrow(uno)){
    cat(j,"\r")
    x <- uno$data[[j]]
    for(i in 2:nrow(x)){ 
      if(is.na(x$New_values[i])) next
      val <- ifelse(x$Corrected_values[i] == "Forest" & x$Corrected_values[i-1] != "Forest","Altered_Forest",x$Corrected_values[i])
      x$Corrected_values[i] <- val
    }
    x -> uno$data[[j]]
  }
  aver <- uno %>% unnest(data)
  rm(uno)
  aver %<>% mutate(Corrected_values = ifelse(is.na(Corrected_values),New_values,Corrected_values))
  
  aver %>% rename(Land_use = Corrected_values) %>% 
    mutate(Land_use = ifelse(Land_use %in% c("Agriculture","Grassland_Shrub"),"Altered_Forest",Land_use)) %>% 
    mutate(Land_use=factor(Land_use,levels=c("Forest","Altered_Forest"))) %>% 
    mutate(year=as.numeric(year)) %>% 
    mutate(tme = year - min(year)) %>% 
    mutate(Elevation_Cat = cut(Altitude,breaks=seq(0,5000,500),labels=seq(0,5000,500)[-1], include.lowest=TRUE)) -> test
  test %>% filter(year == 1982) %>% filter(Land_use=="Forest") %>% pull(CellID) %>% unique() -> target_cells
  chi <- test %>% filter(CellID %in% target_cells) %>% group_by(CellID) %>% nest() %>% rowwise()
  
  # manken_luc <- function(x) {
  #   x %>% select(tme,Land_use) %>% drop_na() -> xx
  #   if(nrow(xx)==0) return(NA) else xx$Land_use %>% ts(.) %>% 
  #     trend::mk.test(.) -> y
  #     y$statistic
  # }
  # 
  # manken_luc_error <- possibly(manken_luc,otherwise = NA)
  # 
  # chi %<>% mutate(manken_lul = manken_luc(data))
}