library(tidyverse)
library(ggplot2)
library(viridis)
library(future)
library(furrr)
library(patchwork)
library(TSA)
library("rnaturalearth")
library("rnaturalearthdata")
library(here)
library(lme4)
library(lavaan)

#### data with at least five columns named: 
#    decimallongitude
#    decimallatitude 
#    Resolved_ACCEPTED (species)
#    year
#    AltRas (altitude) 
### wgetCHELSA_extract: this takes a long time, depending on the connection speed
wgetCHELSA_extract <- function(hist_data="input_data/Historical_data.R",variable="tmax",start=1901,end=2017,file_get="input_data/envidatS3paths.txt"){ 
  load(here(hist_data)) # load data
  reduct %>% dplyr::filter(year >= start & year < end) %>% arrange(year) %>% mutate(YearCat=cut(.$year,end-start)) -> reduct_five ## setup the range of years
  levels(reduct_five$YearCat) %>% stringr::str_sub(.,-7) %>% sub("*.,","",.) %>% sub("\\..*","",.) %>% sub("]","",.) -> levels(reduct_five$YearCat)  
  reduct_five$YearCat <- as.numeric(as.character(reduct_five$YearCat))
  reduct_five$Resolved_ACCEPTED <- factor(reduct_five$Resolved_ACCEPTED)
  reduct_five <- reduct_five %>% mutate(year=paste("Year_",year,sep=""))
  
  # reduct_five %>% group_by(Resolved_ACCEPTED) %>% summarise(Records=n()) %>% filter(Records >= 50)
  
  files_to_get <- readLines(here(file_get), n = -1L)
  files_to_get %>% strsplit(.,variable) %>% lapply(.,"[",2) %>% unlist() %>% sub("_V.1.0.tif ","",.) -> to_get
  reduct_five %>% pull(year) %>% sub("Year_","",.) %>% unique() -> target_years
  to_get %>% strsplit(.,"_") %>% lapply(.,"[",2) %>% unlist() -> to_get
  files_to_get[which(to_get%in%target_years)] -> files_to_get
  files_to_get %>% strsplit(.,"/CHELSAcruts_") %>% lapply(.,"[",2) %>% unlist() %>% sub("_V.1.0.tif","",.) %>% str_trim() -> names
  xy <- reduct_five %>% dplyr::select(decimallongitude,decimallatitude)
  cual=paste("_",target_years,sep="")
  ## The loop downloads the monthly variables (12 per year) and extracts the data for the geographic occurrences. Then it saves a table of n occurrences X 12 variables.
  for (k in 1:length(cual)){
    one <- grep(cual[k],files_to_get)
    a <- list()
    for (i in 1:length(one)){
      cat("Downloading:",names[one][i],"\n")
      cat("------",files_to_get[one][i],"\n")
      system(paste("wget -nv", files_to_get[one][i],"-O temporal.tif",sep=" "))
      raster::raster("temporal.tif") -> r1
      cat("Extracting:",names[one][i],"\n")
      raster::extract(r1,xy) -> a[[i]]
    }
    cat("Mutating.....","\n")
    do.call(cbind,a) -> a_save
    save(a_save,file=paste(here("ClimateData"),sub("_","",variable),cual[k],".RData",sep="/"))
    system("rm temporal.tif")
  }
}
#### Aggregate monthly data into yearly data
month_to_year <- function(cual = "tmax", prec = F) {
focal_file = list.files(path=here("ClimateData"),full.names = T)
if(prec) focal_file <- focal_file[grep(cual,focal_file)]
if(!prec) focal_file <- focal_file[grep(cual,focal_file)]
for (i in 1:length(focal_file)){ 
  variables <- sub(".RData","",focal_file) %>% strsplit("/") %>% lapply(.,"[",2) %>% unlist() 
  cat(variables[i],"\r")
  load(focal_file[i])
  paste("Ann_",variables[i],sep="") -> ann
  paste("Seas_",variables[i],sep="") -> seas
  paste("Max_",variables[i],sep="") -> maxa
  paste("Min_",variables[i],sep="") -> mina
  if(!prec) { 
    if(cual=="tmin")  reduct_five <- reduct_five %>% mutate("{mina}" := a_save %>% apply(.,1,FUN = function(x) min(x)))
    if(cual=="tmax")  reduct_five <- reduct_five %>% mutate("{maxa}" := a_save %>% apply(.,1,FUN = function(x) max(x)))
  }
  if(prec){
    reduct_five <- reduct_five %>% mutate("{ann}" := a_save %>% rowSums()) %>% mutate("{seas}" := a_save %>% apply(.,1,FUN = function(x) sd(x)/mean(x) * 100),"{mina}" := a_save %>% apply(.,1,FUN = function(x) min(x)),"{maxa}" := a_save %>% apply(.,1,FUN = function(x) max(x)))
  }
}

## the code mutates the  yearly variables into the occurrence data base.
save(reduct_five,file = here("temporal_histdataClimate.R"))
load(here("temporal_histdataClimate.R"))
# the next has to be done to get mean temperature and seasonality, for which we need to merge data from tmin and tmax
focal_file = list.files(path=here("ClimateData"),full.names = T)
files_tmin <- focal_file[grep("tmin",focal_file)]
files_tmax <- focal_file[grep("tmax",focal_file)]
length(files_tmax) == length(files_tmin)
for (i in 1:length(files_tmax)){ 
  cat(i,"\r")
  variables <- sub(".RData","",files_tmax) %>% strsplit("_") %>% lapply(.,"[",2) %>% unlist() %>% paste("tmean",.,sep="_")
  variables_2 <- sub(".RData","",files_tmin) %>% strsplit("_") %>% lapply(.,"[",2) %>% unlist()
  (variables==variables_2) %>% unique()
  load(files_tmax[i])
  a_save -> a_tmax
  load(files_tmin[i])
  a_save -> a_tmin
  paste("Ann_",variables[i],sep="") -> ann
  paste("Seas_",variables[i],sep="") -> seas
  cbind(a_tmax,a_tmin) -> global
  reduct_five <- reduct_five %>% mutate("{ann}" := global %>% apply(.,1,FUN = function(x) mean(x)))
  (a_tmax + a_tmin)  / 2 -> mids
  reduct_five <- reduct_five %>% mutate("{seas}" := mids %>% apply(.,1,FUN = function(x) sd(x)/mean(x) * 100))
}

# again the code mutates the variables into the occurrence data base.
save(reduct_five,file=here("output_data/Historical_dataClimate_full.R"))
}
#### PROCESS TRY DATA #####
### Read try data
### Use the TRY database to explore rate of change as a function of growth form.
### I'm summarizing these as: trees, lianas, shrubs/herbs, epiphytes.
process_try <- function(data = here("input_data/TRY/11492.txt"), change_cats = here("input_data/TRY/CATS_TRY.csv")) { 
try_data <- data.table::fread(data)
#try_data %>% distinct(TraitID)
try_data %>% filter(!is.na(TraitID),TraitName=="Plant growth form") -> growth_try
growth_try %>% dplyr::select(SpeciesName,OrigValueStr) %>% mutate(Species=sub(" ","_",SpeciesName)) -> growth_try
growth_try %>% distinct() -> growth_try
#to_plot %>% distinct(OrigValueStr) %>% write.table(.,file="~/Desktop/CATS_TRY.csv",sep=",",row.names=F)
change <- read.table(change_cats,sep=",",header=T)
#unique(change$ChangeValue)
for (i in 1:nrow(change)){
  cat(i,"\n")
  growth_try$OrigValueStr[which(growth_try$OrigValueStr==change$OrigValueStr[i])] <- change$ChangeValue[i] 
}
names(growth_try)[3] <- "Resolved_ACCEPTED"
growth_try %>% as_tibble() %>% filter(OrigValueStr!="") -> growth_try
save(growth_try,file=here("output_data/growth_try.Rdata"))
}

#### FUNCTIONS AND SETTINGS ####
options(dplyr.summarise.inform = FALSE)
quantile_fun <- function(x,perc) ecdf(x)(perc)
splines_predict_bis <- function(mod,der,pred){
  una <- predict(mod,x=pred,deriv = der)
  return(una)
}
min_data_NA <- function(x,year_bins){
  x$YearCat -> z
  which.min(z) -> id
  return(seq(1901,2016,year_bins) >= z[id])
}
splines_mods <- function(x,y,jit,fit){
  npreg::ss(x=x,y=jitter(y,amount=jit),spar=fit,m=2,control.spar = list(lower=0.45,upper=0.7),all.knots=T,nknots = 20)
}
splines_sums <- function(mod){
  una <- summary(mod)
  una <- una$p.table
  return(una)
}
trim_dist <- function(x,lower=0.15,upper=0.85){ 
  xx <- x %>% lapply(.,function(una){ 
    una %>% select(AltRas) %>% summarize(AltRas=quantile(AltRas,na.rm=T,probs=c(lower,upper))) %>% pull() -> ok
    una %>% filter(AltRas >=  all_of(ok[1]) &  AltRas <= all_of(ok[2])) }
  )
}
sem <- function(x, na.rm = T) {
  out <-sd(x, na.rm = na.rm)/sqrt(length(x))
  return(out)
}
Mode_ant <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
Mode <- function(x) {
  ux <- unique(na.omit(x))
  if(length(ux) > 0) {names(which.max(table(x)))} else(NA)
}
dif_clim <- function(x,index,flex_pt) {lapply(x,function(x) {
  x %>% ungroup() %>% filter(YearCat < flex_pt) %>% summarise(across(all_of(index),~ mean(.x,na.rm=T))) -> xx
  simple <- function(z,y) {z-y}
  map2(x[,index],xx,simple) %>% as_tibble() %>% mutate(YearCat=x$YearCat)
})
}
# Experimental
new_ccfForlags <- function (DATA, lags=c(-10:-1), na.action = na.contiguous, na.max.fraction = 0.1) {
  vals <- rep(NA, length(lags))
  names(vals) <- paste("lag", lags)
  if ((sum(complete.cases(DATA)) <= 10) || (mean(complete.cases(DATA)) < na.max.fraction)) 
    return(vals)
  ans <- prewhiten(DATA[, 1], DATA[, 2], lag.max = max(abs(lags)),plot = F,order.max=4)
  vals[] <- drop((ans$ccf)[lags]$acf)
  vals
}
new_rollccf <- function (DATA, width = list(11), by = 1, lags=c(-10:-1), base.lag = estimateDelay(DATA, rises = rises, plot = FALSE), rises = FALSE, na.action = na.omit, na.max.fraction = 0.1){
  DATA <- as.zoo(DATA)
  stopifnot(NCOL(DATA) >= 2)
  if (rises) {
    if ("Q" %in% colnames(DATA)) {
      DATA[, "Q"] <- pmax(diff(DATA[, "Q"], na.pad = TRUE), 
                          0)
    }
    else {
      stop("Give an item 'Q' to use flow rises.")
    }
  }
  if ("Q" %in% colnames(DATA)) {
    whichQ <- which(colnames(DATA) == "Q")
    DATA <- DATA[, c(whichQ, (1:NCOL(DATA))[-whichQ])]
  }
  DATA <- DATA[, 1:2]
  width <- as.list(width)
  if (is.null(names(width))) 
    names(width) <- paste("width", unlist(width))
  obj <- list()
  obj$rolls <- lapply(width, function(wid) {
    tmp <- rollapply(DATA, width = wid, by = by, by.column = FALSE, 
                     FUN = new_ccfForlags, lags = lags, na.action = na.action, 
                     na.max.fraction = na.max.fraction)
    tmp
  })
  
  obj$data <- DATA
  obj$lags <- lags
  obj$width <- width
  obj$call <- match.call()
  class(obj) <- c("rollccf", class(obj))
  obj
}
granger_mult <- function(x.1,x.2=NULL,x.3=NULL,x.4=NULL,x.5=NULL,x.6=NULL,x.7=NULL,x.8=NULL,x.9=NULL,z=NULL,y,n_x=1,lag.order = 3){
  if(length(y[!is.na(y)]) <= 5) return(NULL)
  data.frame(x.1,y) -> m
  if(!is.null(z)) {add_column(.data = m,z,.after="y") -> m}
  if(!is.null(x.2)) {add_column(.data = m,x.2,.before="y") -> m}
  if(!is.null(x.3)) {add_column(.data = m,x.3,.before="y") -> m}
  if(!is.null(x.4)) {add_column(.data = m,x.4,.before="y") -> m}
  if(!is.null(x.5)) {add_column(.data = m,x.5,.before="y") -> m}
  if(!is.null(x.6)) {add_column(.data = m,x.6,.before="y") -> m}
  if(!is.null(x.7)) {add_column(.data = m,x.7,.before="y") -> m}
  if(!is.null(x.8)) {add_column(.data = m,x.8,.before="y") -> m}
  if(!is.null(x.9)) {add_column(.data = m,x.9,.before="y") -> m}
  m[!is.na(m$y),] -> m
  forecast::auto.arima(m$y,biasadj=T,approximation=F,ic="aicc",seasonal=F,stationary=F,max.order=5) -> fit
  m <- apply(m,2,function(x){x - fitted(forecast::Arima(x, model=fit))})
  return(FIAR::partGranger(m, nx=n_x, ny=1, order=lag.order,perm=F) %>% as_tibble(.rows = 1))
}
freq_PF <- function(x) {sum(x %in% c(7:8),na.rm=T)/length(x[!is.na(x)])}
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
substrLeft <- function(x, z, n) {
  substr(x, z,nchar(x)-n)
}

window = 10 ## size of moving window to estimate the rolling mean 
min_by_year = 10 ## minimum number of occurrences to consider a valid data point
min_data = 10 ## minimum number of data points (years) to do regression
limits = c(10,30) # latitudinal limit of the analyses
belts = 200 # the size (in meters) of altitudinal belts to visualize the trends
alt_limit = 3401
records = 50 ## minimum number of records per species
n_sim = 500 ## number of simulations
do_sims = T ## logical to indicate whether to perform simulations
# Ann_tmean Ann_prec Seas_prec Min_prec Max_prec Max_tmax Min_tmin Seas_tmean
start_year = 1901 # starting year to analyse
year_bins = 2
variables = "Min_tmin"
nombres = "Min_tmin"
code = "Min_tmin"
prec = F
seas = F
p_val = 0.01
cual = "Mean" # Mean, Max, Min
cual_n = 1 # mean: 1, max: 2, min: 3
quien = "Raw" # Raw, First
deriv = 0 # 0, 1
k = c(1:10) ### lags
ylimits = c(15,30)
xlimits = c(-115,-80)
lat_limits = c(12.5,26,17)
lon_limits = c(-77,-109,-86)
## getting a couple of Cells with weird estimates
non_id <- c(901353, 918489,918497,959612,968258,1026639)

#### LAND USE DATA ####
lulc_extract <- function(data="output_data/Historical_dataClimate_full.R",lulc_data = "LULC_850_2015") { 
load(here(data))
lista_luc <- list.files(path = here(paste("input_data/","LULC_850_2015",sep="")),pattern = ".tif$",full.names = T) %>% .[grep("1901_",.):grep("2015_",.)]
lucs <- raster::stack(lista_luc)
names(lucs) <- lista_luc %>% strsplit(.,"LULC_") %>% lapply(.,"[",3) %>% sub("_.*","",.) %>% paste("LUC",.,sep="_")
xy <- reduct_five %>% dplyr::select(decimallongitude,decimallatitude)
raster::extract(lucs,xy) -> lucs_ext
lucs_tib <- as_tibble(lucs_ext) %>% bind_cols(reduct_five,.)
save(lucs_tib,file = here("output_data/Historical_dataLULC.R"))
}

alt_sp_trim <- function(data="output_data/Historical_dataClimate_full.R",records=50,lat_limits = c(12.5,26,17),lon_limits = c(-77,-109,-86),output="output_data/Splines_v.4.trim.Rdata") {
  load(here(data))
  reduct_five %>% filter(decimallatitude > lat_limits[1] & decimallongitude < lon_limits[1],
                         !(decimallatitude > lat_limits[2] | decimallongitude < lon_limits[2]),
                         !(decimallatitude > lat_limits[3] & decimallongitude > lon_limits[3]),!is.na(AltRas),YearCat>=start_year,!CellID%in%all_of(non_id)) %>% group_by(Resolved_ACCEPTED) %>%  nest() %>% mutate(data = trim_dist(data,lower=0.01,upper = 0.99)) %>% unnest(cols=data) -> esta
esta %>% summarise(N=n()) -> tata
target <- tata %>% filter(N>=records) %>% pull(Resolved_ACCEPTED) %>% as.character()
esta %>% filter(Resolved_ACCEPTED %in% all_of(target)) %>% droplevels() %>% group_by(Resolved_ACCEPTED) %>%  nest() %>% mutate(data = trim_dist(data,lower=0.01,upper = 0.99)) %>% unnest(cols=data) %>% group_by(Resolved_ACCEPTED,YearCat) %>% summarise(Variable=median(AltRas,na.rm=T),MaxVariable = quantile(AltRas,probs=0.975,na.rm=T),MinVariable=quantile(AltRas,probs=0.025,na.rm=T),n=n(),class=first(class)) %>% mutate_if(is.numeric, ~na_if(., -Inf)) %>% mutate_if(is.numeric, ~na_if(., Inf)) %>% arrange(Resolved_ACCEPTED,YearCat) %>%
    mutate(RollAlt=unlist(tapply(Variable,Resolved_ACCEPTED,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = median ,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
    mutate(RollAlt_Max=unlist(tapply(MaxVariable,Resolved_ACCEPTED,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = quantile,probs=0.9,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
    mutate(RollAlt_Min=unlist(tapply(MinVariable,Resolved_ACCEPTED,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = quantile,probs=0.1,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
    mutate(Rolln=unlist(tapply(n,Resolved_ACCEPTED,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=sum,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
    filter(Rolln >= min_by_year) -> r1
r1 %>% nest_by() %>% mutate(Slopes = list(splines_mods(data$YearCat,data$RollAlt,jit = NULL,fit = NULL)), 
                              Slopes_Max = list(splines_mods(data$YearCat,data$RollAlt_Max,jit = NULL,fit = NULL)),
                              Slopes_Min = list(splines_mods(data$YearCat,data$RollAlt_Min,jit = NULL,fit = NULL))) %>% mutate(SumSlopes = list(splines_sums(Slopes)),SumSlopes_Max = list(splines_sums(Slopes_Max)),SumSlopes_Min = list(splines_sums(Slopes_Min))) -> mean_alts
  save(mean_alts,file=here(output))
}

alt_sims_trim <- function(data="output_data/Historical_dataClimate_full.R",do_sims = T,n_sim = 500,output="output_data/partial_NULLs") { 
  if(do_sims) {
load(here(data))
    ## request parallelization if simulations are set to TRUE
reduct_five %>% filter(decimallatitude > lat_limits[1] & decimallongitude < lon_limits[1],
                           !(decimallatitude > lat_limits[2] | decimallongitude < lon_limits[2]),
                           !(decimallatitude > lat_limits[3] & decimallongitude > lon_limits[3]),!is.na(AltRas),YearCat>=start_year,!CellID%in%all_of(non_id)) %>% group_by(Resolved_ACCEPTED) %>% 
  nest() %>% mutate(data = trim_dist(data,lower=0.01,upper = 0.99)) %>% unnest(cols=data) -> esta
esta %>%  summarise(N=n()) -> tata
target <- tata %>% filter(N>=records) %>% pull(Resolved_ACCEPTED) %>% as.character()
esta %>% filter(Resolved_ACCEPTED %in% all_of(target)) %>% droplevels() %>% group_by(Resolved_ACCEPTED) %>% summarise(Upper=quantile(AltRas,c(0.025),na.rm=T),Lower=quantile(AltRas,c(0.975),na.rm=T)) -> range_r1
# Estimate number of samples per year for all species
esta %>% filter(Resolved_ACCEPTED %in% all_of(target)) %>% droplevels() %>% group_by(YearCat,Resolved_ACCEPTED) %>% summarise(N=n()) %>% arrange(Resolved_ACCEPTED) -> n_years
    nulls_m <- list()
    count = 0
    plan(multisession,workers=8)
    system.time({ 
      for (i in 347:length(target)){ 
        # Identify records (global) within altitudinal range of focal species
        cat(i,toupper(target[i]),"\n")
        count <- count + 1
        range_r1 %>% filter(Resolved_ACCEPTED==target[i]) -> range_focal
        n_years %>% filter(Resolved_ACCEPTED==target[i]) -> years_focal
esta %>% filter(AltRas > range_focal[[2]],AltRas < range_focal[[3]]) %>% droplevels() %>% group_by(YearCat) %>% nest() %>% ungroup() %>% filter(YearCat %in% years_focal[[1]]) -> r_null
        cat("     Performing",n_sim,"simulations","\n")
        ## simulate data by resampling
simulations <- 1:n_sim %>% future_map_dfr(function(j){
          # Sample all records with the same per-year intensity as the observed for the focal species
    r_null %>% inner_join(.,years_focal,by="YearCat") %>% select(-Resolved_ACCEPTED) %>% mutate(samp = map2(data, N, replace=F,sample_n)) %>% dplyr::select(YearCat,N,samp) %>% unnest(samp) %>% droplevels() %>% group_by(YearCat) %>% summarise(Variable=median(AltRas,na.rm=T),MaxVariable = quantile(AltRas,probs=0.975,na.rm=T),MinVariable=quantile(AltRas,probs=0.025,na.rm=T),n=n(),class=first(class)) %>% mutate_if(is.numeric, ~na_if(., -Inf)) %>% mutate_if(is.numeric, ~na_if(., Inf)) %>% mutate(Resolved_ACCEPTED=paste(target[i],"_sim_",j,sep="")) %>% arrange(Resolved_ACCEPTED,YearCat) %>%
            mutate(RollAlt=unlist(tapply(Variable,Resolved_ACCEPTED,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=median,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
            mutate(RollAlt_Max=unlist(tapply(Variable,Resolved_ACCEPTED,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = quantile,probs=0.975 ,na.rm=T,partial=T,align="center")}),use.names = F)) %>% mutate(RollAlt_Min=unlist(tapply(Variable,Resolved_ACCEPTED,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = quantile,probs=0.025 ,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
            mutate(Rolln=unlist(tapply(n,Resolved_ACCEPTED,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=sum,na.rm=T,partial=T,align="center")}),use.names = F)) %>% filter(Rolln >= min_by_year) -> r1
        },.progress = TRUE,.options = furrr_options(seed = 123)) ## set seed to avoid warning!
        cat("     Summarizing","simulations","\n")
        simulations %>% nest_by(Resolved_ACCEPTED) %>% mutate(Slopes = list(splines_mods(data$YearCat,data$RollAlt,jit = NULL,fit = NULL)),Slopes_Max = list(splines_mods(data$YearCat,data$RollAlt_Max,jit = NULL,fit = NULL)),Slopes_Min = list(splines_mods(data$YearCat,data$RollAlt_Min,jit = NULL,fit = NULL))) %>% mutate(SumSlopes = list(splines_sums(Slopes)),SumSlopes_Max = list(splines_sums(Slopes_Max)),SumSlopes_Min = list(splines_sums(Slopes_Min))) -> nulls_m[[i]]
        
        names(nulls_m)[[i]] <- as.character(target[i])
        cat("------","\n")
        if(count == 50) { 
          cat("Saving models (partial)","\n")
          nulls_m %>% do.call(bind_rows,.) -> null_models
          nombre_arch <- paste(here(output),"Splines_Nulls_Partial_Trim_",i,".Rdata",sep="")
          save(null_models,file=nombre_arch) 
          count = 0
          rm(nulls_m,null_models)
          nulls_m <- list()
        }
        if(i==length(target)){
          cat("Saving models (final)","\n")
          nulls_m %>% do.call(bind_rows,.) -> null_models
          nombre_arch <- paste(here(output),"Splines_Nulls_Partial_Trim_",i,".Rdata",sep="")
          save(null_models,file=nombre_arch) 
          rm(nulls_m,null_models)
        }
        
      }
    })
    plan(sequential)
  }}

alt_splines <- function(data="output_data/Historical_dataClimate_full.R",records=50) {
load(here(data))
reduct_five %>% filter(decimallatitude > lat_limits[1] & decimallongitude < lon_limits[1],
                      !(decimallatitude > lat_limits[2] | decimallongitude < lon_limits[2]),
                      !(decimallatitude > lat_limits[3] & decimallongitude > lon_limits[3]),
                      !CellID%in%all_of(non_id),!is.na(AltRas),YearCat>=start_year) %>% group_by(Resolved_ACCEPTED) %>% summarise(N=n()) -> tata
target <- tata %>% filter(N>=records) %>% pull(Resolved_ACCEPTED) %>% as.character()
### The rolling mean is estimated by centering on the corresponding year and with partial estimates allowed (incomplete window)
reduct_five %>% filter(decimallatitude > lat_limits[1] & decimallongitude < lon_limits[1],
         !(decimallatitude > lat_limits[2] | decimallongitude < lon_limits[2]),
         !(decimallatitude > lat_limits[3] & decimallongitude > lon_limits[3]),
         Resolved_ACCEPTED %in% all_of(target),!is.na(AltRas),YearCat>=start_year,!CellID%in%all_of(non_id)) %>% droplevels() %>% group_by(Resolved_ACCEPTED,YearCat) %>% summarise(Variable=median(AltRas,na.rm=T),MaxVariable = quantile(AltRas,probs=0.9,na.rm=T),MinVariable=quantile(AltRas,probs=0.1,na.rm=T),n=n(),class=first(class)) %>% mutate_if(is.numeric, ~na_if(., -Inf)) %>% mutate_if(is.numeric, ~na_if(., Inf)) %>% arrange(Resolved_ACCEPTED,YearCat) %>%
  mutate(RollAlt=unlist(tapply(Variable,Resolved_ACCEPTED,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = median ,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
  mutate(RollAlt_Max=unlist(tapply(MaxVariable,Resolved_ACCEPTED,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = quantile,probs=0.9,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
  mutate(RollAlt_Min=unlist(tapply(MinVariable,Resolved_ACCEPTED,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = quantile,probs=0.1,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
  mutate(Rolln=unlist(tapply(n,Resolved_ACCEPTED,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=sum,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
  filter(Rolln >= min_by_year) -> r1
r1 %>% nest_by() %>% mutate(Slopes = list(splines_mods(data$YearCat,data$RollAlt,jit = NULL,fit = NULL)), 
                            Slopes_Max = list(splines_mods(data$YearCat,data$RollAlt_Max,jit = NULL,fit = NULL)),
                            Slopes_Min = list(splines_mods(data$YearCat,data$RollAlt_Min,jit = NULL,fit = NULL))) %>% mutate(SumSlopes = list(splines_sums(Slopes)),SumSlopes_Max = list(splines_sums(Slopes_Max)),SumSlopes_Min = list(splines_sums(Slopes_Min))) -> mean_alts
save(mean_alts,file=here("output_data/Splines_v.4.Rdata"))
}
alt_sims_nulls <- function(data="output_data/Historical_dataClimate_full.R",do_sims=F) { 
  if(do_sims) {
    load(here(data))
    reduct_five %>% filter(decimallatitude > lat_limits[1] & decimallongitude < lon_limits[1],
                           !(decimallatitude > lat_limits[2] | decimallongitude < lon_limits[2]),
                           !(decimallatitude > lat_limits[3] & decimallongitude > lon_limits[3]),
                           !CellID%in%all_of(non_id),!is.na(AltRas),YearCat>=start_year) %>% group_by(Resolved_ACCEPTED) %>% summarise(N=n()) -> tata
    target <- tata %>% filter(N>=records) %>% pull(Resolved_ACCEPTED) %>% as.character()
    ## request parallelization if simulations are set to TRUE
  reduct_five %>%  filter(decimallatitude > lat_limits[1] & decimallongitude < lon_limits[1],
                          !(decimallatitude > lat_limits[2] | decimallongitude < lon_limits[2]),
                          !(decimallatitude > lat_limits[3] & decimallongitude > lon_limits[3]),
                           !CellID%in%all_of(non_id),!is.na(AltRas)) %>% filter(Resolved_ACCEPTED %in% all_of(target)) %>% droplevels() %>% group_by(Resolved_ACCEPTED) %>% summarise(Upper=quantile(AltRas,c(0.025),na.rm=T),Lower=quantile(AltRas,c(0.975),na.rm=T)) -> range_r1
    # Estimate number of samples per year for all species
    reduct_five %>%  filter(decimallatitude > lat_limits[1] & decimallongitude < lon_limits[1],
                            !(decimallatitude > lat_limits[2] | decimallongitude < lon_limits[2]),
                            !(decimallatitude > lat_limits[3] & decimallongitude > lon_limits[3]),
                           !CellID%in%all_of(non_id),!is.na(AltRas)) %>% filter(Resolved_ACCEPTED %in% all_of(target)) %>% droplevels() %>% group_by(YearCat,Resolved_ACCEPTED) %>% summarise(N=n()) %>% arrange(Resolved_ACCEPTED) -> n_years
    nulls_m <- list()
    count=0
    plan(multisession,workers=8)
    system.time({ 
      for (i in 333:length(target)){ 
        # Identify records (global) within altitudinal range of focal species
        cat(i,toupper(target[i]),"\n")
        count <- count + 1
        range_r1 %>% filter(Resolved_ACCEPTED==target[i]) -> range_focal
        n_years %>% filter(Resolved_ACCEPTED==target[i]) -> years_focal
        reduct_five %>% filter(decimallatitude > 12.5 & decimallongitude < -77,
                               !(decimallatitude > 26 | decimallongitude < -109),
                               !(decimallatitude > 17 & decimallongitude > -86),
                               !CellID%in%all_of(non_id),AltRas<alt_limit) %>%
          filter(AltRas > range_focal[[2]],AltRas < range_focal[[3]]) %>% droplevels() %>% group_by(YearCat) %>% nest() %>% ungroup() %>% filter(YearCat %in% years_focal[[1]]) -> r_null
        cat("     Performing",n_sim,"simulations","\n")
        ## simulate data by resampling
        simulations <- 1:n_sim %>% future_map_dfr(function(j){
          # Sample all records with the same per-year intensity as the observed for the focal species
          r_null %>%  mutate(N=years_focal$N,samp = map2(data, N, replace=F,sample_n)) %>% dplyr::select(-2) %>% unnest(samp) %>% droplevels() %>% group_by(YearCat) %>% summarise(Variable=median(AltRas,na.rm=T),MaxVariable = quantile(AltRas,probs=0.9,na.rm=T),MinVariable=quantile(AltRas,probs=0.1,na.rm=T),n=n(),class=first(class)) %>% mutate_if(is.numeric, ~na_if(., -Inf)) %>% mutate_if(is.numeric, ~na_if(., Inf)) %>% mutate(Resolved_ACCEPTED=paste(target[i],"_sim_",j,sep="")) %>% arrange(Resolved_ACCEPTED,YearCat) %>%
            mutate(RollAlt=unlist(tapply(Variable,Resolved_ACCEPTED,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=median,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
            mutate(RollAlt_Max=unlist(tapply(Variable,Resolved_ACCEPTED,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = quantile,probs=0.9 ,na.rm=T,partial=T,align="center")}),use.names = F)) %>% mutate(RollAlt_Min=unlist(tapply(Variable,Resolved_ACCEPTED,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = quantile,probs=0.1 ,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
            mutate(Rolln=unlist(tapply(n,Resolved_ACCEPTED,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=sum,na.rm=T,partial=T,align="center")}),use.names = F)) %>% filter(Rolln >= min_by_year) -> r1
        },.progress = TRUE,.options = furrr_options(seed = 123)) ## set seed to avoid warning!
        
        cat("     Summarizing","simulations","\n")
        simulations %>% nest_by(Resolved_ACCEPTED) %>% mutate(Slopes = list(splines_mods(data$YearCat,data$RollAlt,jit = NULL,fit = NULL)),Slopes_Max = list(splines_mods(data$YearCat,data$RollAlt_Max,jit = NULL,fit = NULL)),Slopes_Min = list(splines_mods(data$YearCat,data$RollAlt_Min,jit = NULL,fit = NULL))) %>% mutate(SumSlopes = list(splines_sums(Slopes)),SumSlopes_Max = list(splines_sums(Slopes_Max)),SumSlopes_Min = list(splines_sums(Slopes_Min))) -> nulls_m[[i]]
        
        names(nulls_m)[[i]] <- as.character(target[i])
        cat("------","\n")
        if(count == 100) { 
          cat("Saving models (partial)","\n")
          nulls_m %>% do.call(bind_rows,.) -> null_models
          nombre_arch <- paste(here("Splines_Nulls_Partial_"),i,".Rdata",sep="")
          save(null_models,file=nombre_arch) 
          count=0
          rm(nulls_m,null_models)
          nulls_m <- list()
        }
        if(i==length(target)){
          cat("Saving models (final)","\n")
          nulls_m %>% do.call(bind_rows,.) -> null_models
          nombre_arch <- paste(here("Splines_Nulls_Partial_"),i,".Rdata",sep="")
          save(null_models,file=nombre_arch) 
          rm(nulls_m,null_models)
        }
        
      }
    })
    plan(sequential)
    ## revert ot unparallelized session
  }}

### READ NULL MODELS AND OBSERVED, PREDICT AND MERGE ####
##  CHANGE HERE: SumSlopes SumSlopes_Max SumSlopes_Min SumSlopes_Range
##  CHANGE HERE: Slopes Slopes_Max Slopes_Min  Slopes_Range
predicted_alts <- function(data = "output_data/Splines_v.4.trim.Rdata", path_partial = "output_data/partial_NULLs_trim/", start = 1901, end = 2016, year_bins = 2, var_name = "SumSlopes_Range",output="output_data/Predicted_SlopeRange_v4.trim.RData") { 
list_arch <- list.files(path=here(path_partial),pattern = "_Partial_",full.names = T)
load(data)
## Estimate models for the range: max - min 
mean_alts %>% mutate(Type="OBS",.before=data) %>% rowwise(Resolved_ACCEPTED) %>% mutate(Slopes_Range = list(splines_mods(data$YearCat,data$RollAlt_Max - data$RollAlt_Min,jit = NULL,fit = NULL)),.after=Slopes_Min) %>% mutate(SumSlopes_Range = list(splines_sums(Slopes_Range))) -> mean_alts
all_models <- list()
for (i in 1:length(list_arch)){
 cat("Start processing file",i,"\n")
  cat(" ........ Reading...... ","\n")
  null_models <- lapply(list_arch[i], function(x) mget(load(x)))
  null_models %>% lapply(function(x)x[[1]]) %>% .[[1]] %>% ungroup() %>% mutate(Resolved_ACCEPTED=sub("_sim_.*","",Resolved_ACCEPTED)) %>% mutate(Type="NULL",.before=data) %>% rowwise(Resolved_ACCEPTED) %>% mutate(Slopes_Range = list(splines_mods(data$YearCat,data$RollAlt_Max - data$RollAlt_Min,jit = NULL,fit = NULL)),.after=Slopes_Min) %>% mutate(SumSlopes_Range = list(splines_sums(Slopes_Range))) %>% ungroup() -> test_null
 rm(null_models)
  test_null %>% distinct(Resolved_ACCEPTED) %>% pull() -> whichas
  mean_alts %>% filter(Resolved_ACCEPTED %in% all_of(whichas)) -> test ## filter to fit nulls
  test %>% bind_rows(.,test_null) -> test
  # get significant trends and slopes
  lapply(test[[var_name]],"[",2,4) %>% unlist() -> sigs
  lapply(test[[var_name]],"[",2,1) %>% unlist() -> slopes
  test %>% ungroup() %>% mutate(Pvalue=sigs,Trend=slopes,.before=data) %>% rowwise() -> test
  lapply(test$data, min_data_NA,year_bins = year_bins) %>% do.call(rbind,.) %>% na_if(.,0) -> mininois_years ## get recorded years per species
  colnames(mininois_years) <- paste("Year",seq(start,end,year_bins),sep="_") 
  cat("...... calculating derivatives .......","\n")
  test %>% rowwise() %>% mutate(Predicted_Rates = list(splines_predict_bis(Slopes_Range,der=1,pred=seq(start,end,year_bins)))) -> preds ##  CHANGE HERE: Slopes Slopes_Max Slopes_Min  Slopes_Range
  preds %>% rowwise() %>% mutate(Predicted_Vals = list(splines_predict_bis(Slopes_Range,der=0,pred=seq(start,end,year_bins)))) -> preds ##  CHANGE HERE: Slopes Slopes_Max Slopes_Min  Slopes_Range
  #### bind_rows results
  lapply(preds$Predicted_Rates,"[",,2) -> fer 
  do.call(rbind,fer) -> fer
  colnames(fer) <- paste("Year",seq(start,end,year_bins),sep="_") 
  (fer*mininois_years) %>% as_tibble(.) -> fer
  preds %>% select(Resolved_ACCEPTED,Type) %>% bind_cols(.,fer) -> preds_rates
  preds_rates %>% bind_cols(.,test %>% ungroup() %>% select(c(3:4))) -> preds_rates
  preds_rates %>% mutate(Derivative="First",.after=Type) -> preds_rates
  
  lapply(preds$Predicted_Vals,"[",,2) -> fer 
  do.call(rbind,fer) -> fer
  colnames(fer) <- paste("Year",seq(start,end,year_bins),sep="_") 
  (fer*mininois_years) %>% as_tibble(.) -> fer
  preds %>% select(Resolved_ACCEPTED,Type) %>% bind_cols(.,fer) -> preds_vals
  preds_vals %>% bind_cols(.,test %>% ungroup() %>% select(c(3:4))) -> preds_vals
  preds_vals %>% mutate(Derivative="Raw",.after=Type)-> preds_vals
  preds <- rbind(preds_rates,preds_vals)
  cat("       DONE!     ",i,"\n")
  all_models[[i]] <- preds
}
all_models %>% do.call(rbind,.) -> all_models
save(all_models,file=here(output))
}

### SUMMARIZE ALTITUDE PER SPECIES AND MERGE WITH TRY #####
## First load and process try data
sum_species <- function(TRY="output_data/growth_try.Rdata",pattern="Predicted_Slope",quien = "Raw",output="output_data/Species_summaryRaw.RData"){
  load(here(TRY))
list_archs <- list.files(path=here("output_data"),pattern = pattern,full.names = T)
list_archs <- list_archs[grep("v4",list_archs)]
list_archs %>% sub(".*Predicted_","",.) %>% sub("_v4","",.) %>% sub(".RData","",.) -> nombres
species_means <- list()
for (k in 1:length(list_archs)){ 
  load(list_archs[k])
  cat(list_archs[k],"\n")
  all_models %>% ungroup() %>% filter(Type=="NULL",Derivative==quien) %>% group_by(Resolved_ACCEPTED) %>% summarise(across(starts_with("Year"), list(mean=mean,sd=sd))) -> sum_sims
  all_models %>% ungroup() %>% filter(Type=="OBS",Derivative==quien) -> observed
  sum_sims %>% pivot_longer(cols=-1,names_to = c("YearCat","Variable"),names_pattern = "Year_?(.*)_(.*)",values_to = "Estimate") %>% pivot_wider(names_from = Variable,values_from = Estimate) -> sims_vals
  observed %>% pivot_longer(cols=-c(Resolved_ACCEPTED,Type,Pvalue,Trend,Derivative),names_to = c("YearCat"),values_to = "Observed") %>% mutate(YearCat=as.numeric(sub("Year_","",YearCat))) %>% arrange(Resolved_ACCEPTED,YearCat) %>% 
    select(Observed) %>% bind_cols(sims_vals,.) %>% mutate(SES = (Observed - mean) / sd) -> to_plot
  all_models %>% ungroup() %>% filter(Type == "OBS",Derivative==quien) %>% pull(Resolved_ACCEPTED) -> target_spp
  to_plot$SES %>% range(na.rm=T) -> limits
  #### edit......
  to_plot %>% pull(Resolved_ACCEPTED) %>% unique() -> especies
  growth_try %>% filter(Resolved_ACCEPTED %in% all_of(especies)) %>% nest_by(Resolved_ACCEPTED) -> target_try
  target_try %>% rowwise() %>% mutate(LifeForm=(Mode(data$OrigValueStr))) -> target_try
  names(target_try)[1] <- "Resolved_ACCEPTED"
  left_join(to_plot,target_try,"Resolved_ACCEPTED") %>% select(-data) -> to_plot
  to_plot %>% mutate(LifeForm=as.factor(LifeForm)) -> to_plot
  
  levels(to_plot$LifeForm) <- c("Climbers","Epiphytes","Herbs","Shrubs","Shrubs","Trees")
  n_breaks=11
  brks <- seq(limits[1],limits[2], by= (limits[2]-limits[1])/n_breaks)[-c(1,n_breaks+1)]
  
  ### NEED TO ESTIMATE BASELINE TREND AND THEN DO THE DIFFERENCE AND PLOT..... PER SPECIES
  species_means[[k]] <- to_plot %>% filter(Resolved_ACCEPTED %in% all_of(target_spp)) %>% mutate(SES_sig = SES >= 1.96 | SES <= -1.96)
}
names(species_means) <- nombres
save(species_means,file = here(output))
}

##### CLIMATE TRENDS ##### 
# now to estimate the climate trends across cloud forest localities.
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

process_all_series <- function(data="output_data/Historical_dataClimate_full.R",TRY="output_data/growth_try.Rdata",data_lul="output_data/Historical_dataLULC.R",output="output_data/GrangerInput_MaxTrim.Rdata", quien="Raw", cual_n = 2){ 
list_archs <- list.files(path=here("output_data"),pattern="Predicted_Slope",full.names = T) %>% .[grep("v4",.)] 
list_archs %>% sub(".*Predicted_","",.) %>% sub(".RData","",.) %>% sub("_v4","",.) -> nombres
  load(here(data_lul))
  load(here(TRY))
  load(list_archs[cual_n])
  cat(nombres[cual_n],"\n")
lucs_tib %>% filter(decimallatitude > lat_limits[1] & decimallongitude < lon_limits[1],
                      !(decimallatitude > lat_limits[2] | decimallongitude < lon_limits[2]),
                      !(decimallatitude > lat_limits[3] & decimallongitude > lon_limits[3]),
                      !CellID%in%all_of(non_id),AltRas<alt_limit) %>% select(CellID,x,y,AltRas,contains("_prec"),contains("_tm")) %>% distinct(CellID,.keep_all = T) %>% droplevels() %>% pivot_longer(cols = -c(1:4),names_to = "year",values_to = "Values") %>% mutate(year=sub("_tm","Tm",year)) %>% mutate(year=sub("_pre","Pre",year)) %>% separate(year,sep="_",into=c("Var","year")) %>% unite(col="CellID",c(CellID,Var),sep="_",remove=T) %>% arrange(CellID,year) -> aqui
  
aqui %>% mutate(Roll_est = unlist(tapply(Values,CellID,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean, na.rm=T, partial=T, align="center")}),use.names = F)) %>% mutate(YearCat=year) %>% group_by(CellID) %>% nest() %>% mutate(data=dif_clim(data,index = 5:6,flex_pt = 1976)) %>% unnest(data) -> precs
  
  lucs_tib %>% distinct(Resolved_ACCEPTED) -> species_list
  all_models %>% distinct(Resolved_ACCEPTED) -> with_models
  with_models %>% mutate(Models="Model") %>% left_join(species_list,.,by="Resolved_ACCEPTED") %>% mutate(Models=replace_na(Models,"No model")) -> species_list
  growth_try %>% filter(Resolved_ACCEPTED %in% all_of(species_list$Resolved_ACCEPTED)) %>% group_by(Resolved_ACCEPTED) %>% nest() -> target_try
  target_try %>% rowwise() %>% mutate(LifeForm=(Mode(data$OrigValueStr))) -> target_try
  left_join(species_list,target_try,"Resolved_ACCEPTED") %>% select(-data) %>% mutate(LifeForm=replace_na(LifeForm,"Unknown")) %>% mutate(LifeForm=factor(LifeForm)) -> species_list
  levels(species_list$LifeForm) <- c("Climbers","Epiphytes","Herbs","Herbs","Shrubs","Shrubs","Trees","Unknown")
  lucs_tib %>% filter(Resolved_ACCEPTED %in% all_of(with_models %>% pull)) %>% group_by(Resolved_ACCEPTED) %>% select(Resolved_ACCEPTED,CellID) %>% nest() -> tib_species
  tib_species %>% inner_join(.,species_list %>% select(-Models),by="Resolved_ACCEPTED") -> tib_species
  
  precs %>% separate(CellID,into=c("CellID","Variable"),sep="_") -> precs
  precs %>% select(-Values) %>% pivot_wider(names_from = c(Variable,YearCat),values_from = Roll_est) %>% mutate(CellID=as.numeric(CellID)) %>% inner_join(tib_species %>% unnest(data),.,by="CellID") %>% summarise(across(.cols = -c(1:2), mean,na.rm=T)) -> sum_precs

lucs_tib %>% filter(decimallatitude > lat_limits[1] & decimallongitude < lon_limits[1],
                      !(decimallatitude > lat_limits[2] | decimallongitude < lon_limits[2]),
                      !(decimallatitude > lat_limits[3] & decimallongitude > lon_limits[3]),
                      !CellID%in%all_of(non_id),AltRas<alt_limit) %>% select(Resolved_ACCEPTED,x,y,AltRas,starts_with("LUC")) %>% filter(Resolved_ACCEPTED %in% all_of(with_models %>% pull)) %>% group_by(Resolved_ACCEPTED) %>% summarise(across(.cols = -c(1:3), freq_PF)) -> lul_una
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
} 
  #"Observed"  "AnnPrec"   "AnnTmean"  
  #"MaxPrec"   "MaxTmax"   "MinPrec"   "MinTmin" 
  #"SeasPrec"  "SeasTmean" "LUC"

do_Granger <- function(data="output_data/GrangerInput_MaxTrim.Rdata",output="output_data/GrangerFinal_MaxTrim.Rdata") {
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
  
plot_granger <- function(data="output_data/GrangerFinal_Min.Rdata",output="plots/Granger_Min.pdf",TRY="output_data/growth_try.Rdata",  data_lul="output_data/Historical_dataLULC.R") {
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

plot_granger_alternative <- function(data="output_data/GrangerFinal_Min.Rdata",output="plots/Granger_Mid_Alternative.pdf",TRY="output_data/growth_try.Rdata",  data_lul="output_data/Historical_dataLULC.R") {
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

lmm_fit <- function(data="output_data/Splines_v.4.Rdata",cual="Min",flex_point = 13,out_table="output_data/LMM_modelsMin.Rdata"){
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

plot_LMM <- function(data="output_data/LMM_modelsMin.Rdata",cual="Min",outplot="plots/LMM_modelsMin.pdf",outtable="tables/LMM_modelsMin.docx"){
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

#### PLOTS: STAND ALONE PLOTS, WHICH MAKES IT A BIT SLOW
plot_SppSeries <- function(data="output_data/Species_summaryRaw.RData",inflect_point = 1976,output="plots/Species_Series.TRY.pdf",trim=F, by_group=T){ 
load(here(data))
  plotas <- list()
  if (trim==T) {names(species_means) %>% grepl("trim",.) -> cuales} 
  if (trim==F) {names(species_means) %>% grepl("trim",.) %>% !. -> cuales}
  quien <- which(cuales)
  if(by_group) {
    for (k in 1:length(quien)){
      cat("Plotting",k,"\r")
      pp <- species_means[quien][[k]] %>% group_by(Resolved_ACCEPTED) %>% nest() %>% ungroup() %>% mutate(data = dif_clim(data,index = 4:5,flex_pt = inflect_point)) %>% unnest(data) %>% inner_join(.,species_means[quien][[k]] %>% select(Resolved_ACCEPTED,LifeForm),by="Resolved_ACCEPTED") %>% group_by(YearCat,LifeForm) %>% #mutate(SES_sig = SES >= 1.96 | SES <= -1.96) %>% filter(SES_sig) %>% 
        summarise(Observed=median(Observed,na.rm=T),SES=median(SES,na.rm=T)) %>% 
        ggplot(aes(x=as.numeric(YearCat), y=Observed, fill= SES)) +
        geom_bar(stat = "identity") + 
        labs(title = names(species_means)[quien][k],y="Meters",x="") + theme(panel.background = element_blank(),axis.line = element_line()) +
        scale_fill_viridis_b(n.breaks=10,option="G",na.value = "white",direction=-1,limits=c(-2,2)) + geom_vline(linetype="dashed",xintercept = inflect_point) + facet_wrap(~LifeForm)
      NULL
      output2 <- output %>% sub(".pdf","",.) %>% paste0(.,"_",k,".pdf")
ggsave(filename = here(output2),plot = pp)
    }}
  if(!by_group) { for (k in 1:length(quien)){
    cat("Plotting",k,"\r")
plotas[[k]] <- species_means[quien][[k]] %>% group_by(Resolved_ACCEPTED) %>% nest() %>% ungroup() %>% mutate(data = dif_clim(data,index = 4:5,flex_pt = inflect_point)) %>% unnest(data) %>% group_by(YearCat) %>% #mutate(SES_sig = SES >= 1.96 | SES <= -1.96) %>% filter(SES_sig) %>% 
  summarise(Observed=median(Observed,na.rm=T),SES=median(SES,na.rm=T)) %>% 
  ggplot(aes(x=as.numeric(YearCat), y=Observed, fill= SES)) +
  geom_bar(stat = "identity") + 
  labs(title = names(species_means)[quien][k],y="Meters",x="") + theme(panel.background = element_blank(),axis.line = element_line()) +
  scale_fill_viridis_b(n.breaks=10,option="G",na.value = "white",direction=-1,limits=c(-2,2)) + geom_vline(linetype="dashed",xintercept = inflect_point) + 
  NULL
  }
####
myplots <- "
AABB
CCDD
"
pp <- wrap_plots(plotas[-4])  +
  plot_layout(design=myplots,guides = "collect",tag_level = 'new') + 
  plot_annotation(title = "Historical climate trends across cloud forests",
                  subtitle="Deviation from baseline (pre-1976)", 
                  caption= "Ramírez-Barahona et al.") & theme(legend.position = 'right')

ggsave(filename = here(output),plot = pp)
}
}

plot_raw_clim <- function (data="output_data/Historical_dataClimate_full.R",prec=F,variables="Min_tmin"){load(here(data))
  reduct_five %>% filter(decimallatitude > lat_limits[1] & decimallongitude < lon_limits[1],
                         !(decimallatitude > lat_limits[2] | decimallongitude < lon_limits[2]),
                         !(decimallatitude > lat_limits[3] & decimallongitude > lon_limits[3]),
                         !CellID%in%all_of(non_id),AltRas<alt_limit) %>% dplyr::select(CellID,x,y,AltRas,starts_with(variables)) %>% distinct(CellID,.keep_all = T) %>% droplevels() %>% pivot_longer(cols = -c(1:4),names_to = "year",values_to = "Variable") %>% arrange(CellID,year) %>% mutate(Roll_est = unlist(tapply(Variable,CellID,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean, na.rm=T, partial=T, align="center")}),use.names = F))  %>% mutate(YearCat=as.numeric(sub(paste(variables,"_",sep=""),"",year))) %>% group_by(CellID) -> actual
  actual %>% pull(AltRas) %>% cut(.,breaks=seq(0,alt_limit+100,500),labels=seq(0,alt_limit+100,500)[-1], include.lowest=T) -> alt_cats ### aggregate records into altitudinal belts
  actual <- actual %>% ungroup() %>% mutate(AltCats = alt_cats)
  # Estimate means across years.
  actual$AltCats %>% unique() %>% as.character() %>% as.numeric() %>% .[!is.na(.)] %>% sort() -> catas
  catas_2 <- catas - 499
  #####
  #### PLOT ####
  if(prec) {direct=-1; opcion_viridis ="G"}
  if(!prec) {direct=1; opcion_viridis="F"}
  plotas <- list()
  for (i in 1:length(catas)){
    plotas[[i]] <- actual %>% group_by(YearCat,AltCats) %>% summarise(Roll_mean=mean(Roll_est),sem=sem(Roll_est)) %>% filter(AltCats==catas[i]) %>% ggplot(aes(x=YearCat,y=Roll_mean,group=AltCats,color=Roll_mean)) +  geom_point() +
      viridis::scale_color_viridis(discrete=F,option=opcion_viridis,direction=direct) +
      geom_errorbar(aes(ymin = Roll_mean - sem,ymax = Roll_mean +  sem)) + 
      theme(panel.background = element_blank(),legend.position = "none") + labs(y=code,x="",subtitle=paste(catas_2[i],"\u2013",catas[i]," masl",sep="")) + NULL
  }
  
  wrap_plots(plotas) +
    plot_layout(byrow=T,nrow = 3,ncol=4,guides = "collect",tag_level = 'new') + 
    plot_annotation(title = nombres,subtitle="1901\u20132016",
                    caption= "Ramírez-Barahona et al.") & theme(legend.position = '')
}

plot_Chords <- function(lul = "Historical_dataLULC.R",alt_breaks = c(0,300,1000,2000,4000),before = c("LUC_1936","LUC_1956","LUC_1976","LUC_1996"), after = c("LUC_2015","LUC_2015","LUC_2015","LUC_2015"),cat = c("Low","MidLow","MidHigh","High"),codes = c("PF","AU","NF","SV"),targetLU="PF",by_alt=T) {
  load(here(paste("output_data",lul,sep="/")))
lucs_tib %>% filter(decimallatitude > lat_limits[1] & decimallongitude < lon_limits[1],
      !(decimallatitude > lat_limits[2] | decimallongitude < lon_limits[2]),
      !(decimallatitude > lat_limits[3] & decimallongitude > lon_limits[3]),
      !CellID%in%all_of(non_id),AltRas<alt_limit) %>% dplyr::select(CellID,x,y,AltRas,starts_with("LUC")) %>% distinct(CellID,.keep_all = T) %>% droplevels() %>% arrange(CellID) -> precs
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

plot_lulseries <- function(lul = "Historical_dataLULC.R",alt_breaks = c(0,300,1000,2000,4000),cat = c("Low","MidLow","MidHigh","High")){
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

plot_map <- function(data="Historical_dataClimate_full.R",output_name="RichnessMap.pdf",palette = "Greens",title="Species richness") {
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

plot_ClimSeries <- function(data="output_data/ClimateRAW_by_AltCat.R",inflect_point = 1976,output="plots/ClimateSeries_PrecRAW.pdf",variable = "pre",viridis_option = "G",ind=c(2:9)) {
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

working_on <- function(){
data="GrangerFinal_Mid.Rdata"
TRY="output_data/growth_try.Rdata"
data_lul="output_data/Historical_dataLULC.R"
load(here(data))
load(here(data_lul))
load(here(TRY))
lucs_tib %>% distinct(Resolved_ACCEPTED) -> species_list
growth_try %>% filter(Resolved_ACCEPTED %in% all_of(species_list$Resolved_ACCEPTED)) %>% group_by(Resolved_ACCEPTED) %>% nest() -> target_try
target_try %>% rowwise() %>% mutate(LifeForm=(Mode(data$OrigValueStr))) -> target_try
left_join(species_list,target_try,"Resolved_ACCEPTED") %>% select(-data) %>% mutate(LifeForm=replace_na(LifeForm,"Unknown")) %>% mutate(LifeForm=factor(LifeForm)) -> species_list
data="GrangerFinal_Min.Rdata"
load(here(data))
min_final <- final
data="GrangerFinal_Max.Rdata"
load(here(data))
max_final <- final
min_final
min_final %>% mutate(Limit="minimum") -> min_final
min_final %>% unnest(data)
max_final %>% mutate(Limit="maximum") -> max_final
data="GrangerFinal_Mid.Rdata"
load(here(data))
mid_final <- final
mid_final
mid_final %>% mutate(Limit="mid") -> mid_final
max_final %>% bind_rows(.,min_final,mid_final) %>% select(gci_all,Limit) %>% unnest(gci_all) %>% select(orig,prob,Limit) %>% mutate(orig = case_when(prob > 0.05 ~ as.numeric(NA),prob <= 0.05 ~ orig)) %>% select(-prob) %>% 
  #pivot_wider(names_from = Limit,values_from = orig) %>% 
  inner_join(.,species_list,by="Resolved_ACCEPTED") %>% group_by(LifeForm) %>% #summarise(across(c(maximum,minimum),list(m=~mean(.x,na.rm=T),sd=~sd(.x,na.rm=T))))%>%
  ggplot(aes(x=LifeForm,y=orig,fill=Limit)) + 
  geom_boxplot() + scale_fill_viridis_d(option="G",begin=0.4,direction = -1) +
  ylim(c(0,4))



library(ggridges)
# Plot
max_final %>% bind_rows(.,min_final,mid_final) %>% select(gci_all,Limit) %>% unnest(gci_all) %>% select(orig,prob,Limit) %>% mutate(orig = case_when(prob > 0.05 ~ as.numeric(NA),prob <= 0.05 ~ orig)) %>% select(-prob) %>% 
  #pivot_wider(names_from = Limit,values_from = orig) %>% 
  inner_join(.,species_list,by="Resolved_ACCEPTED") %>% group_by(LifeForm) %>% filter(Limit=="mid") %>% 
ggplot(aes(x = orig, y = LifeForm, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.001) +
  scale_fill_viridis(option = "G",direction=-1) +
  scale_x_log10() +
 NULL
}


plot_range <- function(data_lul="output_data/Historical_dataLULC.R",data_spp="output_data/Species_summaryRaw.RData",data_lmm_min="output_data/LMM_modelsMin.Rdata",data_lmm_mid="output_data/LMM_modelsMin.Rdata",data_lmm_max="output_data/LMM_modelsMin.Rdata",which=4,output="plots/Range_Series.pdf"){ 
load(here(data_lul))
load(here(data_spp))
load(here(data_lmm_min))
models[[which]] -> uno
tibble(PRED=predict(uno),Resolved_ACCEPTED=uno@frame$Resolved_ACCEPTED,OBS=uno@frame$Alt,wave1=uno@frame$wave1,wave01=uno@frame$wave01) %>% mutate(Time=wave1+wave01,Range="lower") -> model_frame
data_lmm_mid="output_data/LMM_modelsMid.Rdata"
load(here(data_lmm_mid))
models[[4]] -> uno
tibble(PRED=predict(uno),Resolved_ACCEPTED=uno@frame$Resolved_ACCEPTED,OBS=uno@frame$Alt,wave1=uno@frame$wave1,wave01=uno@frame$wave01) %>% mutate(Time=wave1+wave01,Range="mid") %>% bind_rows(.,model_frame) -> model_frame
data_lmm_max="output_data/LMM_modelsMax.Rdata"
load(here(data_lmm_max))
models[[4]] -> uno

tibble(PRED=predict(uno),Resolved_ACCEPTED=uno@frame$Resolved_ACCEPTED,OBS=uno@frame$Alt,wave1=uno@frame$wave1,wave01=uno@frame$wave01) %>% mutate(Time=wave1+wave01,Range="upper") %>% bind_rows(.,model_frame) -> model_frame

model_frame %>% pivot_wider(names_from = Range, values_from = OBS) %>% group_by(Resolved_ACCEPTED,Time) %>% summarise(across(4:6,mean,na.rm=T)) %>% mutate(Range=abs(lower-upper)) %>% pivot_wider(names_from = Time,values_from = c(upper,mid,lower,Range)) %>% select(starts_with("Range")) %>% rowwise() %>% relocate(c(19,21,16,17),.before=2) %>% relocate(c(20,21),.before=8) %>% mutate(Baseline = mean(c_across(1:13),na.rm=T)) %>% relocate(Baseline,.before=2) %>% mutate(across(2:21, ~ .x - Baseline)) %>% pivot_longer(3:22,names_to = "Time",values_to = "Range") %>% mutate(Time=as.numeric(sub("Range_","",Time))) %>% group_by(Time) %>% summarise(N=median(Range,na.rm=T)) -> ns
pp <- model_frame %>% pivot_wider(names_from = Range, values_from = OBS) %>% group_by(Resolved_ACCEPTED,Time) %>% summarise(across(4:6,mean,na.rm=T)) %>% mutate(Range=abs(lower-upper)) %>% pivot_wider(names_from = Time,values_from = c(upper,mid,lower,Range)) %>% select(starts_with("Range")) %>% rowwise() %>% relocate(c(19,21,16,17),.before=2) %>% relocate(c(20,21),.before=8) %>% mutate(Baseline = mean(c_across(1:13),na.rm=T)) %>% relocate(Baseline,.before=2) %>% mutate(across(2:21, ~ .x - Baseline)) %>% pivot_longer(3:22,names_to = "Time",values_to = "Range") %>% mutate(Time=as.numeric(sub("Range_","",Time))) %>% filter(!is.na(Range))  %>%
  ggplot(aes(x=Time,y=Range)) +
  geom_hline(yintercept = 0,alpha=0.7) + 
  geom_jitter(alpha=0.4,size=0.8,width=0.2,shape=21,aes(fill=Time)) + 
  #ggdist::stat_halfeye(aes(fill=Time), adjust = .9, width = .6, .width = 0, alpha=0.7,slab_colour = "black", slab_size = .2,justification = -.4, height = 2) +
  stat_summary(aes(color=Time),fun = "median",size=1,geom = "crossbar",width=0.7,position = position_nudge(x=-.01),alpha=0.8) +
  stat_summary(aes(color=Time),fun = function(z) { quantile(z,0.10) }, geom = "crossbar",width=0.7,position = position_nudge(x=-.01),alpha=0.8) +
  stat_summary(aes(color=Time),fun = function(z) { quantile(z,0.90) }, geom = "crossbar",width=0.7,position = position_nudge(x=-.01),alpha=0.8) +
  scale_color_viridis_b(n.breaks=50,option="G",direction=1,end=.8) + 
  scale_fill_viridis_b(n.breaks=50,option="G",direction=1,end=.8)+ 
  theme(legend.position = "",panel.background = element_rect(fill="NA",color="black"),panel.grid = element_blank()) + 
  labs(x="Time",y="Deviation in range size") +
  geom_vline(xintercept = 13.5,linetype="dashed") + 
   annotate(geom = "text", x = ns$Time, y = ns$N, fontface = "bold",label = paste0(round(ns$N,1)), hjust = .5, size = 4,vjust=10) +
  NULL
pp
pp <- wrap_plots(list(pp)) +
  plot_layout(guides = "keep",tag_level = 'new') + 
  plot_annotation(title = "Historical altitudinal range for cloud forests",
                  subtitle="Deviations from baseline estimates", 
                  caption= "Ramírez-Barahona et al.") & theme(legend.position = '')
ggsave(filename = here(output),plot = pp)
}



trends <- coef(uno)$Resolved_ACCEPTED %>% as_tibble()
trends <- trends %>% mutate(Resolved_ACCEPTED=coef(uno)$Resolved_ACCEPTED %>% rownames())
lucs_tib %>% filter(decimallatitude > 12.5 & decimallongitude < -77,
       !(decimallatitude > 26 | decimallongitude < -109),
       !(decimallatitude > 17 & decimallongitude > -86),
       !CellID%in%all_of(non_id),AltRas<alt_limit) %>% distinct(CellID,Resolved_ACCEPTED,.keep_all = T) %>% select(CellID,Resolved_ACCEPTED,x,y,AltRas) %>% left_join(.,trends,by="Resolved_ACCEPTED") -> trends
trends %>% distinct(Resolved_ACCEPTED)
trends %>% pull(AltRas) %>% cut(.,breaks=seq(0,alt_limit,belts),labels=seq(0,alt_limit,belts)[-1], include.lowest=T) -> alt_cats ### aggregate records into altitudinal belts
trends <- trends %>% ungroup() %>% mutate(AltCats = alt_cats)
species_means$SlopeMin %>% distinct(Resolved_ACCEPTED) %>% pull() -> names
trends %>% filter(!is.na(wave01)) %>% #filter(Resolved_ACCEPTED %in% all_of(names)) %>% 
  distinct(Resolved_ACCEPTED)
trends %>% filter(!is.na(wave01)) %>% #filter(Resolved_ACCEPTED %in% all_of(names)) %>% 
  group_by(CellID) %>% summarise(across(5:7, mean),across(2:4,first)) -> to_plot
range(c(to_plot$wave01,to_plot$wave1),na.rm=T) -> limits_vals
p <- to_plot %>% ggplot() + geom_sf(data=ne_countries(scale=110,returnclass = "sf",type="countries"),colour="grey95",fill="grey95",size=0.5) + xlim(xlimits) + ylim(ylimits) + 
geom_tile(aes(x=x,y=y,fill=wave01),color="black",size=0.05) + 
  scale_fill_viridis_b(n.breaks=20,limits=limits_vals,name="Slope") + 
  #scale_fill_gradient2(limits=limits_vals,name="Slope") +
  theme(panel.background = element_rect(fill="NA",colour = "black"),panel.grid=element_blank(),axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())
pp <- to_plot %>% ggplot() + geom_sf(data=ne_countries(scale=110,returnclass = "sf",type="countries"),colour="grey95",fill="grey95",size=0.5) + xlim(xlimits) + ylim(ylimits) + 
geom_tile(aes(x=x,y=y,fill=wave1),color="black",size=0.05) + 
  scale_fill_viridis_b(n.breaks=20,limits=limits_vals,name="Slope") + 
  #scale_fill_gradient2(limits=limits_vals,name="Slope") +
  theme(panel.background = element_rect(fill="NA",colour = "black"),panel.grid=element_blank(),axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())
wrap_plots(list(p,pp)) +
  plot_layout(byrow=T,nrow = 2,ncol=1,guides = "collect",tag_level = 'new') + 
  plot_annotation(title = "Per-cell mean species' trends", subtitle="Breaking point 1975-1976", caption= "Ramírez-Barahona et al.") & theme(legend.position = 'bottom')




taget_spp = "Cyathea_divergens"
lucs_tib %>% filter(Resolved_ACCEPTED==taget_spp) -> test
test %>% distinct(CellID,.keep_all = T) %>%  select(CellID,x,y,YearCat,contains("tmean")) %>% pivot_longer(cols = -c(1:4),names_to = "Variable",values_to = "Values") %>% mutate(Year = as.numeric(substrRight(Variable,4))) %>% mutate(Match= YearCat==Year) %>% filter(Match) %>% mutate(Variable = (substrLeft(Variable,1,5))) %>% pivot_wider(names_from = "Variable",values_from = "Values") -> test
lucs_tib %>% filter(Resolved_ACCEPTED==taget_spp) %>% distinct(CellID,.keep_all = T) %>% select(CellID,contains("LUC")) %>% left_join(test,.,by="CellID") %>% filter(Ann_tmean!=max(Ann_tmean,na.rm=T),Ann_tmean!=min(Ann_tmean,na.rm=T),Seas_tmean!=max(Seas_tmean,na.rm=T),Seas_tmean!=min(Seas_tmean,na.rm=T)) -> test

library(gganimate)
aa <- test %>% mutate(across(contains("LUC"), ~ case_when(.x %in% c(7,8) ~ "Remaining",!.x %in% c(7,8) ~ "Lost"))) %>% pivot_longer(cols = -c(1:8),names_to = "Period",values_to = "Present") %>% mutate(Trans=as.integer(substrRight(Period,4))) %>% 
ggplot(aes(x=Ann_tmean,y=Seas_tmean,fill=Present,color=Present)) + 
geom_point(alpha=0.9,shape=21,size=3) + 
#stat_ellipse(size=1.2) +
#facet_wrap(~Period) + 
theme(panel.background = element_blank(),panel.grid=element_blank(),axis.line = element_line()) + scale_fill_manual(values= c("tomato","black"),name="Forest Cover") + scale_color_manual(values= c("tomato","black"),name="Forest Cover") + labs(title="Deforestation and species occurrences - {frame_time}",subtitle=taget_spp,caption="Ramírez-Barahona et al.",x="Annual mean temperature",y="Temperature Seasonality") + transition_time(Trans)

anim_save(aa,file="~/Desktop/Procrastination.gif",duration=30)

