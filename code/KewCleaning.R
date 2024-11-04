library(tidyverse)
library(magrittr)
library(ggdist)
library(patchwork)
library(here)
library(sp)
library(raster,exclude = "select")
library("rnaturalearth")
library("rnaturalearthdata")
library(ggplot2)
library(showtext)
library(furrr)
library(glue)
library(showtext)
library(ggridges)
library(ggeffects)
library(lme4)
library(MetBrewer)
font_add_google(name="EB Garamond")
showtext_auto()

###### PRE-PROCESSING: SETTING UP DATA FILES RETRIEVED FROM GBIF #####
###### USE BASH TO REDUCE FILE SIZE BY SELECTING COLUMNS
### bash commands are commented out!
# unzip 0188599-210914110416597.zip
# rename unzipped file to Traqueos_NeoTropics_raw.csv
# awk -F"\t" '{print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$13"\t"$22"\t"$23"\t"$33}' Traqueos_NeoTropics_raw.csv > Traqueos_NeoTropics_COR.csv
# wc -l TRAQUEOS_raw.csv
# wc -l Traqueos_NeoTropics_COR.csv

###### Keep records identified to species level or below ########

#### FUNCTIONS AND SETTINGS ####
#### SOME MIGHT BE DEPRECATED (27 MARCH 2023)
total <- function(x) { 
  x %>% summarise(n_tot = n()) %>% pull(n_tot) -> n_tot
  return(n_tot)
}
total_rbcL <- function(x) { 
  x %>% count(rbcL) %>% filter(!is.na(rbcL)) %>% pull(n) -> n_rbcl
  if(length(n_rbcl) == 0) n_rbcl <- 0
  return(n_rbcl)
}
font_add_google(name="EB Garamond")
showtext_auto()
options(dplyr.summarise.inform = FALSE)
quantile_fun <- function(x,perc) ecdf(x)(perc)
splines_predict_bis <- function(mod,der,pred){
  una <- predict(mod,x=pred,deriv = der)
  return(una)
}

min_data_NA <- function(x,year_bins){
  x$YearCat -> z
  which.min(z) -> id
  return(seq(1979,2019,year_bins) >= z[id])
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
seme <- function(x, na.rm = T) {
  out <-sd(x, na.rm = na.rm) / sqrt(length(which(!is.na(x))))
  return(out)
}
cv <- function(x, na.rm = T) {
  out <-sd(x, na.rm = na.rm)/mean(x, na.rm = na.rm) * 100
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
  x %>% ungroup() %>% filter(YearCat <= flex_pt) %>% summarise(across(all_of(index),~ mean(.x,na.rm=T))) -> xx
  simple <- function(z,y) {z-y}
  map2(x[,index],xx,simple) %>% as_tibble() %>% mutate(YearCat=x$YearCat)
})
}
dif_clim_2 <- function(x,index,flex_pt) {lapply(x,function(x){ 
  x %>% group_split(Var) -> x
  res <- lapply(x,FUN = function(i){ 
    i %>% ungroup() %>% filter(as.numeric(YearCat) <= flex_pt) %>% 
      summarise(across(all_of(index),~ mean(.x,na.rm=T))) -> xx
    simple <- function(z,y) {z-y}
    w <- map2(i[,index],xx,simple) %>% as_tibble() %>% 
      rename(Dif_Values=Values,Dif_Roll_est=Roll_est) %>% 
      bind_cols(i,.)
    return(w)
  })
  res %>% do.call(rbind,.)
})
}
parse_freq <- function(x,category) { length(x[which(x==category)])/ length(x) }
parse_count <- function(x) { length(x) }
parse_modal <- function(x) {x %>% tapply(.,.,length) %>% .[which.max(.)] %>% names(.)}
una_lm <- function(y,x) {lm(y ~ x,na.action = "na.omit") %>% coefficients(.) %>% .[2]}
man_ken_tau_deprecated <- function(x) {x %>% pivot_wider(names_from = var,values_from = val) %>% 
    select(contains("mean_tas_")) %>% drop_na()  -> xx
  if(ncol(xx)==0) return(NA) else xx  %>% as_vector() %>% Kendall::MannKendall() %>% .[["tau"]] %>% .[1]}
man_ken_s <- function(x) {Kendall::MannKendall(x)$S}
manken <- function(x,tag="_cmi") {x %>% pivot_wider(names_from = var,values_from = val) %>% 
    select(contains(tag)) %>% select(-contains("CV")) %>% 
    drop_na() -> xx
  if(ncol(xx)==0) return(NA) else xx %>% ts(.) %>% trend::mult.mk.test(.) %>% .[["statistic"]]
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
normalize <- function(x){ (x - min(x,na.rm=T)) / (max(x,na.rm=T)- min(x, na.rm=T)) }
downsampling <- function(){
  some4 %>% count(year,species) %>% filter(year >= 1950 & year < 2011) %>% 
    mutate(n = ifelse(n <= median(n), n, median(n)))-> N
  
  some4 %>% group_by(year,species) %>% filter(year >= 1950 & year < 2011) %>% nest() %>% left_join(.,N,by=c("year","species")) %>% mutate(data=map2(data,n,replace=FALSE,sample_n)) %>% unnest(data) -> rarified
}
rescale.coefs <- function(beta,mu,sigma) {
  beta2 <- beta ## inherit names etc.
  ## parameters other than intercept:
  beta2[-1] <- sigma[1]*beta[-1]/sigma[-1]
  ## intercept:
  beta2[1]  <- sigma[1]*beta[1] + mu[1]-sum(beta2[-1]*mu[-1])
  return(beta2)
}   
rescale.coefs_rand <- function(beta,sigma) {
  beta2 <- beta ## inherit names etc.
  ## parameters other than intercept:
  beta2 <- sigma[1]*beta / sigma[-1]
  return(beta2)
}


window = 10 ## size of moving window to estimate the rolling mean 
min_by_year = 10 ## minimum number of occurrences to consider a valid data point
min_data = 10
belts = 200 # the size (in meters) of altitudinal belts to visualize the trends
alt_limit = 3401
records = 100 ## minimum number of records per species
n_sim = 10 ## number of simulations
do_sims = T ## logical to indicate whether to perform simulations
# Ann_tmean Ann_prec Seas_prec Min_prec Max_prec Max_tmax Min_tmin Seas_tmean
start_year = 1979 # starting year to analyse
year_bins = 2
code = nombres = variables = "tas"
prec = F
seas = F
cual = "Mean" # Mean, Max, Min
cual_n = 1 # mean: 1, max: 2, min: 3
quien = "Raw" # Raw, First
deriv = 0 # 0, 1
k = c(1:10) ### lags
ylimits = c(15,30)
xlimits = c(-115,-80)
## getting a couple of Cells with weird estimates
non_id <- c(901353, 918489,918497,959612,968258,1026639)


#### FUNCTIONS FOR FILTERING AND PROCESSING ####
# These functions are from Ramírez-Barahona, Cuervo-Robayo, Magallón (2023)
gbif_pow_exactjoin <- function(gbif_data = "Traqueos_NeoTropics_COR.csv",
                               checklist = "wcvp_names.txt",
                               taxon_status_cats = c("Accepted","Orthographic","Synonym","Unplaced"),
                               output_file="Joined_base.Rdata") { 
aa <- data.table::fread(here("data",gbif_data)) %>% as_tibble()
aa %<>% mutate(ID=1:nrow(.)) %>% relocate(ID,.before=1) %>% separate(scientificName,into=c("G","S","extra"),extra = "merge",sep=" ") %>% unite("binomial",c(G,S),sep="_",remove=F) %>% unite("scientificName",c(G,S,extra),sep=" ",remove=T) %>% filter(species!="")
list_names <- data.table::fread(here("wcvp_2022",checklist),sep="|",quote="") %>% as_tibble() # %>% ## data provided by Kew (we are unable to provide it here)
list_names %<>% unite("scientificName",c(taxon_name,taxon_authors),sep=" ",remove = F) %>% filter(species!="") %>% mutate(Accepted_Name = scientificName[match(accepted_plant_name_id,plant_name_id,nomatch=NA)]) %>% as_tibble() %>%  filter(taxon_status %in% all_of(taxon_status_cats)) %>% unite("binomial",c(genus,species),sep="_",remove = F) %>% filter(accepted_plant_name_id!="")

## EXACT JOIN 
joined <- aa %>% #rename("taxon_name"=scientificName) %>% 
  left_join(.,list_names,by = "scientificName")
cat("Searching for duplicate IDs due to multiple matches!","\n")
joined %>% pull(ID) %>% duplicated() -> dups
joined$ID[which(dups)] -> id_dups
cat("Found",length(id_dups),"duplicated IDs","\n")
joined %>% filter(ID %in% all_of(id_dups)) -> to_correct
joined <- joined[which(!joined$ID %in% id_dups),]
pow_distributions="wcvp_distribution.txt"
pow_dist <- data.table::fread(here("wcvp_2022",pow_distributions),sep="|",quote="") %>% as_tibble() %>% #left_join(.,list_names,by="accepted_db_id") %>% 
  filter(introduced==0)
joined %<>% mutate(Dist_correction=T)
to_correct %<>% mutate(Dist_correction=T)
cat("First step: resolving by distribution","\n")
for(i in 1:nrow(to_correct)){
  cat("processing",i,"\r")
pow_dist %>% filter(plant_name_id==to_correct$accepted_plant_name_id[i]) %>% pull(continent_code_l1) -> tar_dist
tar_dist %in% c(1:6,9) -> res
if(sum(res)!=0) {to_correct$Dist_correction[i] <- FALSE}
}
to_correct %<>% filter(Dist_correction)
rbind(to_correct,joined) -> joined
joined %>% pull(ID) %>% duplicated() -> dups
joined$ID[which(dups)] -> id_dups
cat("Remaining",length(id_dups),"duplicated IDs","\n")
joined %>% filter(ID %in% all_of(id_dups)) -> to_correct
joined <- joined[which(!joined$ID %in% id_dups),]
to_correct %>% filter(homotypic_synonym=="") %>% group_split(ID) -> to_correct
cat("Second step: resolve manually","\n")
for(i in 1:length(to_correct)){
  cat("processing",i,"\n")
  if(nrow(to_correct[[i]])==1) next
to_correct[[i]] %>% select(scientificName,Accepted_Name,homotypic_synonym,accepted_plant_name_id,taxon_status,geographic_area) %>% print()
to_select <- as.numeric(readline("Which one to keep?"))
if(is.na(to_select)) next
to_correct[[i]] %<>% slice(-all_of(to_select))
}
data.table::rbindlist(to_correct) %>% as_tibble() -> to_correct
rbind(to_correct,joined) -> joined
joined %>% pull(ID) %>% duplicated() -> dups
joined$ID[which(dups)] -> id_dups
cat("Still",length(id_dups),"some unsolved records!","\n")
joined %>% filter(is.na(plant_name_id)) %>% distinct(scientificName) %>% nrow() -> spp
cat(spp,"unresolved names left!","\n")
save(joined,file=here("interim/",output_file))
}
gbif_pow_first_fuzzy <- function(data="Joined_base.Rdata",
                                 checklist = "wcvp_names.txt",
                                 taxon_status_cats = c("Accepted","Orthographic","Synonym","Unplaced"),
                                 output_file="Joined_fuzzy.Rdata"){ 
load(here("interim",data))
joined %>% filter(is.na(plant_name_id)) %>% distinct(scientificName) -> spp
  spp %>% filter(grepl("×",scientificName)) -> hybrids
  spp %>% filter(grepl("\\?",scientificName)) -> noidea
  spp %>% filter(!grepl("×",scientificName),!grepl("\\?",scientificName)) -> species

list_names <- data.table::fread(here("wcvp_2022",checklist),sep="|",quote="") %>% as_tibble() # %>% ## data provided by Kew (we are unable to provide it here)
  list_names %<>% unite("scientificName",c(taxon_name,taxon_authors),sep=" ",remove = F)
  list_names %<>% filter(species!="") %>% mutate(Accepted_Name = scientificName[match(accepted_plant_name_id,plant_name_id,nomatch=NA)]) %>% as_tibble() %>% filter(taxon_status %in% all_of(taxon_status_cats)) %>% unite("binomial",c(genus,species),sep="_",remove = F) %>% filter(accepted_plant_name_id!="")

joined_nana <- joined %>% filter(is.na(plant_name_id))
joined %<>% filter(!is.na(plant_name_id))

resolved_matches <- list()
non_matches <- list()
multiple_matches <- list()
max.dist=0.15
for (i in 1:nrow(species)){
  species %>% slice(i) %>% pull(scientificName) -> quien
  cat(i,"Attempting match on:", quien,"\n")
  joined_nana %>% filter(scientificName %in% all_of(quien)) -> target
target %>% slice(1) %>% select(binomial.x,scientificName,accepted_plant_name_id,Accepted_Name) %>% rename(binomial=binomial.x) -> sole_target

list_names %>% filter(binomial== sole_target$binomial) -> names_tomatch
fuzzyjoin::stringdist_left_join(species %>% slice(i),names_tomatch,by="scientificName",method="jw",distance_col="distance",max_dist=max.dist) -> result
if(nrow(result)>1){ 
  if(length(unique(result$accepted_plant_name_id))==1){
    cat("        Single match!","\n")
    result[1,] -> transfer
    target$accepted_plant_name_id <- transfer$accepted_plant_name_id
    target$Accepted_Name <- transfer$Accepted_Name
    resolved_matches[[i]] <- target
    } 
  if(length(unique(result$accepted_plant_name_id))>1) {
     if(min(result$distance)<0.01){
       cat("        Single match!","\n")
       result[which.min(result$distance),] -> transfer
       target$accepted_plant_name_id <- transfer$accepted_plant_name_id
       target$Accepted_Name <- transfer$Accepted_Name
       resolved_matches[[i]] <- target
     }
     if(min(result$distance)>0.01){
   cat("        Multiple matches!","\n")
   target$accepted_plant_name_id <- NA
   target$Accepted_Name <- NA
   multiple_matches[[i]] <- target
    }
  }
  next
}
if(is.na(result$plant_name_id)){ cat("       No match!","\n")
  target$accepted_plant_name_id <- NA
  target$Accepted_Name <- NA
  non_matches[[i]] <- target
  next}
if(!is.na(result$plant_name_id)){
  cat("        Single match!","\n")
result %>% select(scientificName.x,accepted_plant_name_id,Accepted_Name) -> transfer
target$accepted_plant_name_id <- transfer$accepted_plant_name_id
target$Accepted_Name <- transfer$Accepted_Name
resolved_matches[[i]] <- target
}
}

resolved_matches %>% data.table::rbindlist(.) %>% as_tibble() -> joined_fuzzy
joined_fuzzy %<>% bind_rows(joined,.)
joined_fuzzy %>% distinct(Accepted_Name)

save(joined_fuzzy,file=here("interim",output_file))
save(non_matches,file=here("interim","non_matches_fuzzy.Rdata"))
save(multiple_matches,file=here("interim","multiples_matches_fuzzy.Rdata"))

}

gbif_pow_second_fuzzy <- function(output_file="Joined_final.v1.Rdata"){ 
load(here("interim/non_matches_fuzzy.Rdata"))
  load("interim/Joined_fuzzy.Rdata")
  non_matches %<>% data.table::rbindlist() %>% as_tibble() 
  non_matches %>% distinct(scientificName) -> species
resolved_matches <- list()
non_matches_final <- list()
multiple_matches_final <- list()
max.dist=0.15
for (i in 1:nrow(species)){
  species %>% slice(i) %>% pull(scientificName) -> quien
  cat(i,"Attempting match on:", quien,"\n")
  non_matches %>% filter(scientificName %in% all_of(quien)) -> target
  target %>% slice(1) %>% select(binomial.x,scientificName,accepted_plant_name_id,Accepted_Name,genus.x) %>% rename(binomial=binomial.x,genus=genus.x) -> sole_target
  
  list_names %>% filter(genus == sole_target$genus) -> names_tomatch
  fuzzyjoin::stringdist_left_join(species %>% slice(i),names_tomatch,by="scientificName",method="jw",distance_col="distance",max_dist=max.dist) -> result
  if(nrow(result)>1){ 
    if(length(unique(result$accepted_plant_name_id))==1){
      cat("        Single match!","\n")
      result[1,] -> transfer
      target$accepted_plant_name_id <- transfer$accepted_plant_name_id
      target$Accepted_Name <- transfer$Accepted_Name
      resolved_matches[[i]] <- target
    } 
    if(length(unique(result$accepted_plant_name_id))>1) {
      if(min(result$distance)<0.01){
        cat("        Single match!","\n")
        result[which.min(result$distance),] -> transfer
        target$accepted_plant_name_id <- transfer$accepted_plant_name_id
        target$Accepted_Name <- transfer$Accepted_Name
        resolved_matches[[i]] <- target
      }
      if(min(result$distance)>0.01){
        cat("        Multiple matches!","\n")
        target$accepted_plant_name_id <- NA
        target$Accepted_Name <- NA
        multiple_matches_final[[i]] <- target
      }
    }
    next
  }
  if(is.na(result$plant_name_id)){ cat("       No match!","\n")
    target$accepted_plant_name_id <- NA
    target$Accepted_Name <- NA
    non_matches_final[[i]] <- target
    next}
  if(!is.na(result$plant_name_id)){
    cat("        Single match!","\n")
    result %>% select(scientificName.x,accepted_plant_name_id,Accepted_Name) -> transfer
    target$accepted_plant_name_id <- transfer$accepted_plant_name_id
    target$Accepted_Name <- transfer$Accepted_Name
    resolved_matches[[i]] <- target
  }
}
non_matches_final %>% data.table::rbindlist(.) %>% as_tibble() 
multiple_matches %>% data.table::rbindlist(.) %>% as_tibble() %>% distinct(scientificName)

resolved_matches %>% data.table::rbindlist(.) %>% as_tibble() -> joined_fuzzy_f
joined_fuzzy %<>% bind_rows(.,joined_fuzzy_f)
joined_fuzzy %>% distinct(Accepted_Name)
save(joined_fuzzy,file=here("output",output_file))
save(non_matches_final,file="output/no_matches_final.v1.Rdata")
save(multiple_matches_final,file="output/multiples_matches_final.v1.Rdata")

}
geographic_filter <- function(data="Joined_final.v1.Rdata",output_file = "Joined_finalPOW.v1.Rdata",
                              perform_tests=c("centroids","institutions", "equal", "gbif","capitals", "zeros","seas")) {
  load(here("output",data))
  ### APPLY FILTERS. NOT USING THE OUTLIER TEST, BECAUSE WE BASE THIS ON KEW'S DISTRIBUTIONS.
  cat("Filtering using CoordinateCleaner","\r")
  to_filter <- joined_fuzzy %>% select(decimalLongitude,decimalLatitude) %>% as.data.frame() 
  f <- CoordinateCleaner::clean_coordinates(to_filter, lon = "decimalLongitude",lat = "decimalLatitude",centroids_rad = 1000, centroids_detail = "both",tests=perform_tests,species=NULL,value="flagged",seas_scale = 110)
  joined_fuzzy %<>% mutate(Outlier_Test = f)
  save(joined_fuzzy,file=here("output",output_file))
}
powdist_filter <- function(data="Joined_finalPOW.v1.Rdata",output_file = "Joined_finalPOWdist.v1.Rdata",
                              pow_distributions="wcvp_distribution.txt",
                              wgsrpd = "level3/level3.shp") { 
  load(here("output",data))
  pow_dist <- data.table::fread(here("wcvp_2022",pow_distributions),sep="|",quote="") %>% as_tibble() %>% #left_join(.,list_names,by="accepted_db_id") %>% 
    filter(introduced==0)
  cat("....... reading polygons from WGSRPD","\r")
  poly <- rgdal::readOGR(here("wgsrpd-master",wgsrpd))
  #cat("....... buffering polygons from WGSRPD","\r")
  #spList = vector("list", length(poly))
  #for (i in 1:length(poly)) {
  #  cat(i,"\r")
  #  a <- rgeos::gBuffer(poly[i,], width = 0.5)
  #  a$LEVEL3_COD = poly[i,]$LEVEL3_COD
  #  a$LEVEL3_NAM = poly[i,]$LEVEL3_NAM
  #  spList[[i]] <- a
  ##}
  #poly <- do.call("rbind", spList)
  pointos <- SpatialPoints(as.data.frame(joined_fuzzy[,11:10]))
  proj4string(pointos)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  points_powo <- sp::over(pointos,poly)
  joined_fuzzy %<>% bind_cols(.,points_powo)
  cat("Filtering using KEW's distributions","\n")
  joined_fuzzy %>% distinct(Accepted_Name)
  joined_fuzzy %>% group_split(Accepted_Name) -> aver
  some <- list()
  for (i in 1:length(aver)){
      cat(i,"--","Getting POW data for", unique(aver[[i]]$Accepted_Name),"               ","\r")
    temp <- aver[[i]] %>% mutate(POW_distribution = FALSE)
    if(is.na(unique(aver[[i]]$Accepted_Name))) {
      some[[i]] <- temp %>% mutate(POW_distribution = TRUE)
      }
    matching <- !is.na(match(temp$LEVEL3_COD,pow_dist %>% filter(plant_name_id==unique(temp$accepted_plant_name_id)) %>% pull(area_code_l3),nomatch = NA))
    if(length(which(matching)) == 0) { 
      some[[i]] <- temp %>% mutate(POW_distribution = TRUE)
      next
    }
    if(length(which(matching)) != 0) some[[i]] <- temp %>% mutate(POW_distribution = matching)
  }
  some %<>% data.table::rbindlist(.) %>% as_tibble()
  cat("-------------- DONE --------------","\r")
  save(some,file=here("output",output_file))
}

#### These function are specific to the cloud forest analyses
## HARMONIZE AND QUERY VILLASEÑOR CHECKLIST
harmonize_bmm_checklist <- function(checklist = "wcvp_names.txt",
                                    taxon_status_cats = c("Accepted","Orthographic","Synonym","Unplaced"),
                                    bmm_checklist = "Lista_SPP.txt",
                                    output_file = "Lista_SPP_fuzzy.txt"){
list_names <- data.table::fread(here("wcvp_2022",checklist),sep="|",quote="") %>% as_tibble() # %>% ## data provided by Kew (we are unable to provide it here)
list_names %<>% unite("scientificName",c(taxon_name,taxon_authors),sep=" ",remove = F) %>% filter(species!="") %>% mutate(Accepted_Name = scientificName[match(accepted_plant_name_id,plant_name_id,nomatch=NA)]) %>% as_tibble() %>%  filter(taxon_status %in% all_of(taxon_status_cats)) %>% unite("binomial",c(genus,species),sep="_",remove = F) %>% filter(accepted_plant_name_id!="")

villa <- read.table(here("data",bmm_checklist),sep=";") %>% as_tibble() %>% select(1) %>% filter(V1!="") %>%
separate(V1,into=c("genus_villa","species_villa","author_villa"),remove=F,extra="merge") %>% filter(!is.na(author_villa)) %>% filter(!grepl("XX",genus_villa)) %>% rename("scientificName"=V1) %>% unite("binomial",c(genus_villa,species_villa),sep="_",remove = F)
villa
villa %<>% left_join(.,list_names,by="scientificName")
villa %>% filter(is.na(accepted_plant_name_id))
villa %>% filter(is.na(accepted_plant_name_id)) %>% distinct(scientificName) %>% filter(!grepl("×",scientificName)) -> species
villa_nana <- villa %>% filter(is.na(accepted_plant_name_id))
villa %<>% filter(!is.na(plant_name_id))

resolved_matches <- list()
non_matches <- list()
multiple_matches <- list()
max.dist=0.15
for (i in 1:nrow(species)){
  species %>% slice(i) %>% pull(scientificName) -> quien
  cat(i,"Attempting match on:", quien,"\n")
  villa_nana %>% filter(scientificName %in% all_of(quien)) -> target
  target %>% slice(1) %>% select(binomial.x,scientificName,accepted_plant_name_id,Accepted_Name) %>% rename(binomial=binomial.x) -> sole_target
  list_names %>% filter(binomial== sole_target$binomial) -> names_tomatch
  fuzzyjoin::stringdist_left_join(species %>% slice(i),names_tomatch,by="scientificName",method="jw",distance_col="distance",max_dist=max.dist) -> result
  if(nrow(result)>1){ 
    if(length(unique(result$accepted_plant_name_id))==1){
      cat("        Single match!","\n")
      result[1,] -> transfer
      target$accepted_plant_name_id <- transfer$accepted_plant_name_id
      target$Accepted_Name <- transfer$Accepted_Name
      resolved_matches[[i]] <- target
    } 
    if(length(unique(result$accepted_plant_name_id))>1) {
      if(min(result$distance)<0.01){
        cat("        Single match!","\n")
        result[which.min(result$distance),] -> transfer
        target$accepted_plant_name_id <- transfer$accepted_plant_name_id
        target$Accepted_Name <- transfer$Accepted_Name
        resolved_matches[[i]] <- target
      }
      if(min(result$distance)>0.01){
        cat("        Multiple matches!","\n")
        target$accepted_plant_name_id <- NA
        target$Accepted_Name <- NA
        multiple_matches[[i]] <- target
      }
    }
    next
  }
  if(is.na(result$plant_name_id)){ cat("       No match!","\n")
    target$accepted_plant_name_id <- NA
    target$Accepted_Name <- NA
    non_matches[[i]] <- target
    next}
  if(!is.na(result$plant_name_id)){
    cat("        Single match!","\n")
    result %>% select(scientificName.x,accepted_plant_name_id,Accepted_Name) -> transfer
    target$accepted_plant_name_id <- transfer$accepted_plant_name_id
    target$Accepted_Name <- transfer$Accepted_Name
    resolved_matches[[i]] <- target
  }
}

resolved_matches %>% data.table::rbindlist(.) %>% as_tibble() -> villa_fuzzy
villa_fuzzy %<>% bind_rows(villa,.)

non_matches %<>% data.table::rbindlist() %>% as_tibble() 
non_matches %>% distinct(scientificName) -> species
resolved_matches <- list()
non_matches_final <- list()
multiple_matches_final <- list()
max.dist=0.15
for (i in 1:nrow(species)){
  species %>% slice(i) %>% pull(scientificName) -> quien
  cat(i,"Attempting match on:", quien,"\n")
  non_matches %>% filter(scientificName %in% all_of(quien)) -> target
  target %>% slice(1) %>% select(binomial.x,scientificName,accepted_plant_name_id,Accepted_Name,genus_villa) %>% rename(binomial=binomial.x,genus=genus_villa) -> sole_target
  
  list_names %>% filter(genus == sole_target$genus) -> names_tomatch
  fuzzyjoin::stringdist_left_join(species %>% slice(i),names_tomatch,by="scientificName",method="jw",distance_col="distance",max_dist=max.dist) -> result
  if(nrow(result)>1){ 
    if(length(unique(result$accepted_plant_name_id))==1){
      cat("        Single match!","\n")
      result[1,] -> transfer
      target$accepted_plant_name_id <- transfer$accepted_plant_name_id
      target$Accepted_Name <- transfer$Accepted_Name
      resolved_matches[[i]] <- target
    } 
    if(length(unique(result$accepted_plant_name_id))>1) {
      if(min(result$distance)<0.01){
        cat("        Single match!","\n")
        result[which.min(result$distance),] -> transfer
        target$accepted_plant_name_id <- transfer$accepted_plant_name_id
        target$Accepted_Name <- transfer$Accepted_Name
        resolved_matches[[i]] <- target
      }
      if(min(result$distance)>0.01){
        cat("        Multiple matches!","\n")
        target$accepted_plant_name_id <- NA
        target$Accepted_Name <- NA
        multiple_matches_final[[i]] <- target
      }
    }
    next
  }
  if(is.na(result$plant_name_id)){ cat("       No match!","\n")
    target$accepted_plant_name_id <- NA
    target$Accepted_Name <- NA
    non_matches_final[[i]] <- target
    next}
  if(!is.na(result$plant_name_id)){
    cat("        Single match!","\n")
    result %>% select(scientificName.x,accepted_plant_name_id,Accepted_Name) -> transfer
    target$accepted_plant_name_id <- transfer$accepted_plant_name_id
    target$Accepted_Name <- transfer$Accepted_Name
    resolved_matches[[i]] <- target
  }
}
non_matches_final %>% data.table::rbindlist(.) %>% as_tibble() 
multiple_matches %>% data.table::rbindlist(.) %>% as_tibble() %>% distinct(scientificName)

resolved_matches %>% data.table::rbindlist(.) %>% as_tibble() -> villa_fuzzy_f
villa_fuzzy %<>% bind_rows(.,villa_fuzzy_f)
villa_fuzzy
write.table(villa_fuzzy,file=here("output",output_file),sep="|",quote=F,row.names = F,col.names = T)
}
extract_demMesoamerica <- function(data ="Joined_finalPOWdist.v1.Rdata",
                                   villa_fuzzy = "Lista_SPP_fuzzy.txt",
                                   dem = "EarthEnv_DEM90.tif"){ 
load(here("output",data))
villa_fuzzy <- data.table::fread(here("output",villa_fuzzy))
some %<>% select(ID:taxon_status,infraspecific_rank:POW_distribution)
villa_fuzzy %>% distinct(Accepted_Name) %>% mutate(In_Villa = TRUE) %>% 
  right_join(.,some, by="Accepted_Name") %>% as_tibble() -> some
dem90 <- raster(here("data",dem))
some %>% select(decimalLongitude,decimalLatitude) -> coords
coords %>% extract(dem90,.) %>% enframe(.) %>% rename_with(.,~all_of(c("name","Altitude"))) %>% select(Altitude) -> coords

some %>% bind_cols(coords,.) -> some

some %<>% filter(!is.na(Altitude)) %>% select(Altitude:family.x,scientificName:year,accepted_plant_name_id,Dist_correction:POW_distribution) 

save(some,file = "interim/Historical_dem.Rdata")

}

### FILTER BY VILLA AND THEN PROCESS
### Ecoregion spatial overlay
filter_bmm <- function(data ="Historical_dem.Rdata") {
load(here("interim",data))
g <- raster::raster(nrows=180*24,ncols=360*24,xmn=-180,xmx=180,ymn=-90,ymx=90,vals=1) %>% 
  projectRaster(.,crs = CRS("+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m  +no_defs")) %>% as(., 'SpatialPixels')
some_villa <- some %>% filter(In_Villa)
some_villa %>% select(decimalLongitude,decimalLatitude) %>% 
  sf::st_as_sf(x = .,coords = c("decimalLongitude","decimalLatitude"),crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>% 
  sf::st_transform(.,crs = proj4string(g)) %>%  
  sf::st_coordinates() %>% SpatialPoints(.,proj4string= CRS("+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m  +no_defs")) %>% sp::over(.,g) %>% enframe(.,name="name") %>% rename_with(.,~all_of(c("name","CellID"))) %>% 
  select(CellID) %>% bind_cols(.,some_villa) -> gridded

poly <- rgdal::readOGR("~/Documents/7.MAPOTECA/official/wwf_terr_ecos.shp")
pointos <- SpatialPoints(as.data.frame(gridded[,11:10]))
proj4string(pointos)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
points_ecos <- sp::over(pointos,poly)
gridded %<>% bind_cols(.,points_ecos)
gridded %<>% mutate(BIOME=as.factor(BIOME))
levels(gridded$BIOME) <- c("TropMoistForest","TropDryForest","TropConiferForest","TempConiferForest","TropGrasslands","MediterraneanForest","DesertsXeric","Mangroves",NA)
gridded %>% count(BIOME) %>% mutate(prop=n/sum(n)) ### Proportion of occurrences in each BIOME
##### Need to clean a bit more using biomes (although this is a very rough filter)
### The following biomes are 'suspect' (not necessarily incorrect): cleaning
bimas <- c("DesertsXeric","MediterraneanForest","Mangroves") ### These are the biomes to eliminate
spp <- unique(gridded$Accepted_Name[which(gridded$BIOME %in% bimas)])
non <- gridded[which(some$Accepted_Name %in% spp),]
nona <- tapply(non$BIOME,non$Accepted_Name, function(x)(round(100*(table(x) / length(x)),2)))
spp <- list()
for (i in 1:length(bimas)){
  spp[[i]] <- names(which(sapply(nona,"[[",bimas[i]) > 5)) ### Applying a threshold to every biome: species with more than this prop. of records are eliminated
}
spp <- unique(unlist(spp))
gridded2 <- gridded[-which(gridded$Accepted_Name %in% spp),]

gridded3 <- gridded2 %>% 
  filter(ECO_NAME%in%c("Belizian pine forests","Central American montane forests","Central American pine-oak forests","Central American Atlantic moist forests","Chiapas montane forests","Chimalapas montane forests","Oaxacan montane forests","Petén-Veracruz moist forests","Sierra de los Tuxtlas","Sierra Madre de Chiapas moist forests","Sierra Madre de Oaxaca pine-oak forests","Sierra Madre del Sur pine-oak forests","Veracruz montane forests","Veracruz moist forests","Sierra Madre Occidental pine-oak forests","Sierra Madre Oriental pine-oak forests","Trans-Mexican Volcanic Belt pine-oak forests",
"Isthmian-Pacific moist forests","Isthmian-Atlantic moist forests","Talamancan montane forests"))
gridded3 %>% filter(Outlier_Test,POW_distribution) -> reduct
save(reduct,file="interim/Historical_data.Rdata")
} 
wgetCHELSA_extract <- function(hist_data="interim/Historical_data.Rdata",
                               variable="tas",
                               start=1979,end=2019,
                               file_get="data/envidatS3paths_tas.txt"){ 
  load(here(hist_data)) # load data
  reduct %>% filter(year >= start & year <= end) %>% arrange(year) %>% mutate(YearCat=factor(year)) %>% 
  mutate(Accepted_Name=factor(Accepted_Name)) %>% mutate(year=paste("Year_",year,sep="")) %>% select(1:12,ECO_NAME:BIOME) -> reduct_five
  files_to_get <- readLines(here(file_get), n = -1L)
  files_to_get %>% strsplit(.,variable) %>% lapply(.,"[",3) %>% 
    unlist() %>% sub("_V.2.1.tif ","",.) %>% stringr::str_sub(.,-4) -> to_get
  reduct_five %>% pull(year) %>% sub("Year_","",.) %>% unique() -> target_years
  files_to_get[which(to_get %in% target_years)] -> files_to_get
  files_to_get %>% strsplit(.,"/CHELSA_") %>% lapply(.,"[",2) %>% unlist() %>% sub("_V.2.1.tif","",.) %>% str_trim() -> names
  xy <- reduct_five %>% select(decimalLongitude,decimalLatitude)
  cual = paste("_",target_years,sep="")
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
    save(a_save,file=paste(here("climate_data_CHELSA/"),sub("_","",variable),cual[k],".RData",sep=""))
    system("rm temporal.tif")
  }
}
month_to_year <- function(data="Historical_data.Rdata",start=1979,end=2019) {
  load(here("interim",data))
  reduct %>% filter(year >= start & year <= end) %>% arrange(year) %>% mutate(YearCat=factor(year)) %>% 
    mutate(Accepted_Name = factor(Accepted_Name)) %>% mutate(year=paste("Year_",year,sep="")) %>% 
    select(1:12,ECO_NAME:BIOME) -> reduct_five
  all_file = list.files(path=here("data/climate_data_CHELSA"),full.names = T)
  for (i in 1:length(all_file)){ 
    variables <- sub(".RData","",all_file) %>% strsplit("/") %>% lapply(.,"[",9) %>% unlist() 
    cat(variables[i],"\r")
    load(all_file[i])
    if(nrow(a_save)!=nrow(reduct_five)) {paste("Error in file ---",variables,"----");next}
    ok <-  a_save %>% apply(.,1,FUN = function(x) mean(x,na.rm=T))# %>% as_tibble() %>% rename(tas=value)
      name = paste("mean",variables[i],sep="_")
      reduct_five %<>% mutate("{name}" := ok)
    ok <-  a_save %>% apply(.,1,FUN = function(x) min(x, na.rm = T))# %>% as_tibble() %>% rename(tas=value)
      name = paste("min",variables[i],sep="_")
      reduct_five %<>% mutate("{name}" := ok)
    ok <-  a_save %>% apply(.,1,FUN = function(x) max(x, na.rm=T))# %>% as_tibble() %>% rename(tas=value)
      name = paste("max",variables[i],sep="_")
      reduct_five %<>% mutate("{name}" := ok)
    ok <-  a_save %>% apply(.,1,FUN = function(x) (sd(x,na.rm=T)/mean(x,na.rm=T)) * 100)# %>% as_tibble() %>% rename(tas=value)
     name = paste("CV",variables[i],sep="_")
     reduct_five %<>% mutate("{name}" := ok)
     
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
  save(reduct_five,file=here("output/Historical_data_FULL.v.1.R"))
}
process_try <- function(data = here("data/TRY/11492.txt"), change_cats = here("data/TRY/CATS_TRY.csv")) { 
  try_data <- data.table::fread(data)
  try_data %>% filter(!is.na(TraitID),TraitName=="Plant growth form") -> growth_try
  growth_try %>% select(SpeciesName,OrigValueStr) %>% mutate(Species=sub(" ","_",SpeciesName)) -> growth_try
  growth_try %>% distinct() -> growth_try
  #to_plot %>% distinct(OrigValueStr) %>% write.table(.,file="~/Desktop/CATS_TRY.csv",sep=",",row.names=F)
  change <- read.table(change_cats,sep=",",header=T)
  #unique(change$ChangeValue)
  for (i in 1:nrow(change)){
    cat(i,"\r")
    growth_try$OrigValueStr[which(growth_try$OrigValueStr==change$OrigValueStr[i])] <- change$ChangeValue[i] 
  }
  names(growth_try)[3] <- "Accepted_Name"
  growth_try %>% as_tibble() %>% filter(OrigValueStr!="") -> growth_try
  save(growth_try,file=here("output/growth_try.Rdata"))
 }
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
luc_extract <- function(data="output/Historical_data_FULL.v.1.R",lulc_data = "GLASS-GLC") { 
  load(here(data))
  lista_luc <- list.files(path = here("data",lulc_data),pattern = ".tif$",full.names = T)
  lucs <- raster::stack(lista_luc)
  names(lucs) <- lista_luc %>% strsplit(.,"7classes_") %>% lapply(.,"[",2) %>% sub(".tif$","",.) %>% paste("LUC",.,sep="_")
  xy <- reduct_five %>% select(decimalLongitude,decimalLatitude)
  raster::extract(lucs,xy) -> lucs_ext
  reduct_five <- as_tibble(lucs_ext) %>% bind_cols(reduct_five,.)
  save(reduct_five,file = here("output/Historical_data_FULL_GLASS.v.1.R"))
  
}
alt_splines <- function(data="output/Historical_data_FULL_LUL.v.1.R",records=100,output="output/Splines_v.2.Rdata") {
  cat("Loading data","\n")
  load(here(data))
reduct_five %>% select(CellID:Accepted_Name,ID,class:family.x,decimalLatitude:BIOME) %>% 
    mutate(YearCat = as.numeric(sub("Year_","",year)),.after=1) %>% 
    filter(!is.na(Altitude),YearCat >= start_year) %>% distinct(Accepted_Name,CellID,year,.keep_all = T) %>% group_by(Accepted_Name) -> esta
  esta %>% summarise(N=n()) -> tata
  target <- tata %>% filter(N >= records) %>% pull(Accepted_Name) %>% as.character(.)
  esta %>% filter(Accepted_Name %in% all_of(target)) %>% droplevels() %>% group_by(Accepted_Name,YearCat) %>% 
    summarise(Variable=mean(Altitude,na.rm=T), MaxVariable = quantile(Altitude,probs=0.975,na.rm=T),MinVariable=quantile(Altitude,probs=0.025,na.rm=T),n=n(),class=first(class)) %>% ungroup() %>% 
    mutate_if(is.numeric, ~na_if(., -Inf)) %>% mutate_if(is.numeric, ~na_if(., Inf)) %>% arrange(Accepted_Name,YearCat) %>%
    mutate(RollAlt=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean ,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
    mutate(RollAlt_Max=unlist(tapply(MaxVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
    mutate(RollAlt_Min=unlist(tapply(MinVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
    
    #mutate(RollAlt_error=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=seme,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
    #mutate(RollAlt_Max_error=unlist(tapply(MaxVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = seme,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
    #mutate(RollAlt_Min_error=unlist(tapply(MinVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = seme,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
    mutate(RollAlt_sd=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = sd,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
    mutate(RollAlt_Maxsd=unlist(tapply(MaxVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = sd,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
    mutate(RollAlt_Minsd=unlist(tapply(MinVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = sd,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
    
    mutate(Rolln=unlist(tapply(n,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=sum,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
    filter(Rolln >= min_by_year) -> r1
  cat("Summarising data","\n")
r1 %>% nest_by(Accepted_Name) %>% mutate(data_points = nrow(data)) %>% filter(data_points >= 10) %>% 
mutate(Slopes = list(splines_mods(data$YearCat,data$RollAlt,jit = NULL,fit = NULL)), 
       Slopes_Max = list(splines_mods(data$YearCat,data$RollAlt_Max,jit = NULL,fit = NULL)),
       Slopes_Min = list(splines_mods(data$YearCat,data$RollAlt_Min,jit = NULL,fit = NULL))) %>% mutate(SumSlopes = list(splines_sums(Slopes)),SumSlopes_Max = list(splines_sums(Slopes_Max)),SumSlopes_Min = list(splines_sums(Slopes_Min))) -> mean_alts
  
# mean_alts %>% unnest(data) %>% ggplot(aes(x=YearCat,y=RollAlt,group=Accepted_Name,color=Accepted_Name),
#                        alpha=0.2) + geom_path() + theme(legend.position = "")
  
save(mean_alts,file=here(output))

  }
alt_resamp_nulls <- function(data="output/Historical_data_FULL_LUL.v.1.R",n_sim=500,start_year=1979,records=100) {
  cat("Loading data","\n")
  load(here(data))
  reduct_five %>% select(CellID:Accepted_Name,ID,class:family.x,decimalLatitude:BIOME) %>% 
      mutate(YearCat = as.numeric(sub("Year_","",year)),.after=1) %>% 
      filter(!is.na(Altitude),YearCat >= start_year) %>% 
      distinct(Accepted_Name,CellID,year,.keep_all = T) %>% group_by(Accepted_Name) -> esta
    esta %>% summarise(N=n()) -> tata
    target <- tata %>% filter(N >= records) %>% pull(Accepted_Name) %>% as.character(.)
    ## request parallelization if simulations are set to TRUE
    esta %>% filter(Accepted_Name %in% all_of(target)) %>% droplevels() %>% group_by(Accepted_Name) %>% summarise(Lower=quantile(Altitude,c(0.025),na.rm=T),Upper=quantile(Altitude,c(0.975),na.rm=T)) -> range_r1
    # Estimate number of samples per year for all species
    esta %>% filter(Accepted_Name %in% all_of(target)) %>% group_by(YearCat,Accepted_Name) %>% summarise(N=n()) %>% arrange(Accepted_Name) -> n_years
    nulls_m <- list()
    count=0
    cat("Setting up multisessions","\n")
    plan(multisession,workers=8)
    #system.time({ 
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
          filter(Altitude >= range_focal[[2]],Altitude <= range_focal[[3]]) %>% droplevels() %>% group_by(YearCat) %>%
     nest() %>% ungroup() %>% filter(YearCat %in% years_focal[[1]]) -> r_null
        cat("     Performing",n_sim,"simulations","\n")
        ## simulate data by resampling
simulations <- 1:n_sim %>% future_map_dfr(function(j){
          # Sample all records with the same per-year intensity as the observed for the focal species
          r_null %>%  mutate(N=years_focal$N,samp = map2(data, N, replace=F,sample_n)) %>% select(-2) %>% unnest(samp) %>% droplevels() %>% group_by(YearCat) %>% summarise(Variable=mean(Altitude,na.rm=T),MaxVariable = quantile(Altitude,probs=0.975,na.rm=T),MinVariable=quantile(Altitude,probs=0.025,na.rm=T),n=n(),class=first(class)) %>% mutate_if(is.numeric, ~na_if(., -Inf)) %>% mutate_if(is.numeric, ~na_if(., Inf)) %>% mutate(Accepted_Name=paste(target[i],"_sim_",j,sep="")) %>% arrange(Accepted_Name,YearCat) %>%
      mutate(RollAlt=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=mean,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
      mutate(RollAlt_Max=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
      mutate(RollAlt_Min=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
            
      mutate(RollAlt_error=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=seme,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
      mutate(RollAlt_Max_error=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = seme,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
      mutate(RollAlt_Min_error=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = seme,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
            
      mutate(Rolln=unlist(tapply(n,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=sum,na.rm=T,partial=T,align="center")}),use.names = F)) %>% filter(Rolln >= min_by_year) -> r1
        },.progress = TRUE,.options = furrr_options(seed = 123)) ## set seed to avoid warning!
    cat("     Summarizing","simulations","\n")
    simulations %>% nest_by(Accepted_Name) %>% mutate(Slopes = list(splines_mods(data$YearCat,data$RollAlt,jit = NULL,fit = NULL)),Slopes_Max = list(splines_mods(data$YearCat,data$RollAlt_Max,jit = NULL,fit = NULL)),Slopes_Min = list(splines_mods(data$YearCat,data$RollAlt_Min,jit = NULL,fit = NULL))) %>% mutate(SumSlopes = list(splines_sums(Slopes)),SumSlopes_Max = list(splines_sums(Slopes_Max)),SumSlopes_Min = list(splines_sums(Slopes_Min))) -> nulls_m[[i]]
        
        names(nulls_m)[[i]] <- as.character(target[i])
        cat("------","\n")
        if(count == 100) { 
          cat("Saving models (partial)","\n")
          nulls_m %>% do.call(bind_rows,.) -> null_models
          nombre_arch <- paste(here("interim/splines_partial","Splines_Nulls_Partial_"),i,".Rdata",sep="")
          save(null_models,file=nombre_arch) 
          count=0
          rm(nulls_m,null_models)
          nulls_m <- list()
        }
        if(i==length(target)){
          cat("Saving models (final)","\n")
          nulls_m %>% do.call(bind_rows,.) -> null_models
          nombre_arch <- paste(here("interim/splines_partial","Splines_Nulls_Partial_"),i,".Rdata",sep="")
          save(null_models,file=nombre_arch) 
          rm(nulls_m,null_models)
        }
        
      }
    #})
    plan(sequential)
    ## revert ot unparallelized session
  }
alt_downsamp_nulls <- function(data="output/Historical_data_FULL_LUL.v.1.R",n_sim=100,start_year=1979,records=100) {
  load(here(data))
  reduct_five %>% select(CellID:Accepted_Name,ID,class:family.x,decimalLatitude:BIOME) %>% 
    mutate(YearCat = as.numeric(sub("Year_","",year)),.after=1) %>% 
    filter(!is.na(Altitude),YearCat >= start_year) %>% 
    distinct(Accepted_Name,CellID,year,.keep_all = T) %>% group_by(Accepted_Name) -> esta
  esta %>% summarise(N=n()) -> tata
  target <- tata %>% filter(N >= records) %>% pull(Accepted_Name) %>% as.character(.)
  ## request parallelization if simulations are set to TRUE
  esta %>% filter(Accepted_Name %in% all_of(target)) %>% count(year,Accepted_Name) %>%
    mutate(n = ifelse(n <= median(n), n, median(n)))-> N
  esta %>% filter(Accepted_Name %in% all_of(target)) %>% group_by(year,Accepted_Name) %>% 
    nest() %>% left_join(.,N,by=c("year","Accepted_Name")) -> esta
  plan(multisession,workers=8)
  #system.time({
      cat("     Performing",n_sim,"simulations","\n")
      ## simulate data by resampling
      simulations <- 1:n_sim %>% future_map_dfr(function(j){
      esta %>% mutate(data=map2(data,n,replace=F,sample_n)) %>% unnest(data) %>% mutate(Accepted_Name=paste(Accepted_Name,"_sim_",j,sep="")) %>% arrange(Accepted_Name,YearCat) -> rarified
rarified %>% ungroup() %>% group_by(Accepted_Name,YearCat) %>% summarise(Variable=mean(Altitude,na.rm=T),MaxVariable = quantile(Altitude,probs=0.975,na.rm=T),MinVariable=quantile(Altitude,probs=0.025,na.rm=T),n=n(),class=first(class)) %>% mutate_if(is.numeric, ~na_if(., -Inf)) %>% mutate_if(is.numeric, ~na_if(., Inf)) %>% arrange(Accepted_Name,YearCat) %>% 
          mutate(RollAlt=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=mean,na.rm=T,partial=T,align="center")}),use.names = F)) %>% mutate(RollAlt_Max=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean,na.rm=T,partial=T,align="center")}),use.names = F)) %>% mutate(RollAlt_Min=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
          mutate(RollAlt_error=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=seme,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
          mutate(RollAlt_Max_error=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = seme,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
          mutate(RollAlt_Min_error=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = seme,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
          
mutate(Rolln=unlist(tapply(n,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=sum,na.rm=T,partial=T,align="center")}),use.names = F)) %>% filter(Rolln >= min_by_year) -> r1
      },.progress = TRUE,.options = furrr_options(seed = 123)) ## set seed to avoid warning!
      cat("     Summarizing","simulations","\n")
      simulations %>% ungroup() %>% nest_by(Accepted_Name) %>% mutate(Slopes = list(splines_mods(data$YearCat,data$RollAlt,jit = NULL,fit = NULL)),Slopes_Max = list(splines_mods(data$YearCat,data$RollAlt_Max,jit = NULL,fit = NULL)),Slopes_Min = list(splines_mods(data$YearCat,data$RollAlt_Min,jit = NULL,fit = NULL))) %>% mutate(SumSlopes = list(splines_sums(Slopes)),SumSlopes_Max = list(splines_sums(Slopes_Max)),SumSlopes_Min = list(splines_sums(Slopes_Min))) -> nulls_m
        cat("Saving models (final)","\n")
        nombre_arch <- here("interim","Splines_Down_Partial_nulls_5.Rdata")
        save(nulls_m,file=nombre_arch) 
        rm(nulls_m)
        rm(simulations)
  #})
  plan(sequential)
  ## revert ot unparallelized session
}
alt_downsamp_nulls2 <- function(data="output/Historical_data_FULL_LUL.v.1.R",n_sim=100,start_year=1979,records=100) {
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
    rarified %>% ungroup() %>% group_by(Accepted_Name,YearCat) %>% summarise(Variable=mean(Altitude,na.rm=T),MaxVariable = quantile(Altitude,probs=0.975,na.rm=T),MinVariable=quantile(Altitude,probs=0.025,na.rm=T),n=n(),class=first(class)) %>% mutate_if(is.numeric, ~na_if(., -Inf)) %>% mutate_if(is.numeric, ~na_if(., Inf)) %>% arrange(Accepted_Name,YearCat) %>% 
      mutate(RollAlt=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=mean,na.rm=T,partial=T,align="center")}),use.names = F)) %>% mutate(RollAlt_Max=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean,na.rm=T,partial=T,align="center")}),use.names = F)) %>% mutate(RollAlt_Min=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
      mutate(RollAlt_error=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=seme,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
      mutate(RollAlt_Max_error=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = seme,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
      mutate(RollAlt_Min_error=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = seme,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
      
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

### READ NULL MODELS AND OBSERVED, PREDICT AND MERGE ####
##  CHANGE HERE: SumSlopes SumSlopes_Max SumSlopes_Min SumSlopes_Range
##  CHANGE HERE: Slopes Slopes_Max Slopes_Min  Slopes_Range
predicted_alts <- function(data = "output/Splines_v.2.Rdata", 
                           path_resample = "interim/splines_partial", 
                           path_downsample = "interim", start = 1979, end = 2019, year_bins = 1, 
                           var_name = "SumSlopes",output="output/Predicted_Slopes_v.2.RData") { 
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
  if(i!=length(list_arch_ds)){test <- test_null}
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

### SUMMARIZE ALTITUDE PER SPECIES AND MERGE WITH TRY #####
## First load and process try data
sum_species <- function(TRY="output/growth_try.Rdata",pattern="Predicted_Slope",
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


### Figures and tables
plot_bmm_things <- function(data="output/Historical_data_FULL_LUL.v.1.R",data_spp = "output/Species_summaryRaw.RData") {
  load(here(data))  ### occurrences
  load(data_spp)    ### summary (vs. nulls)
  species_means[[1]] %>% count(Accepted_Name) %>% pull(Accepted_Name) -> targets
  reduct_five %<>% filter(Accepted_Name%in% targets)
g <- raster::raster(nrows=180*24,ncols=360*24,xmn=-180,xmx=180,ymn=-90,ymx=90,vals=1) %>% 
    projectRaster(.,crs = CRS("+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m  +no_defs")) %>% as(., 'SpatialPixels')

reduct_five %>% 
  mutate(year=sub("Year_","",year)) %>% mutate(year = as.numeric(year)) %>% 
  mutate(Decade = cut(year,breaks=c(1979,1990,2000,2010,2020),include.lowest=TRUE,labels = c("80s","90s","00s","10s") ),.before=Altitude) %>% 
  distinct(Accepted_Name,CellID,.keep_all = T) %>% count(CellID) %>% left_join(.,reduct_five %>% distinct(CellID,.keep_all = T),by="CellID") -> SR
  g@coords %>% as_tibble() %>% mutate(CellID=1:nrow(.)) %>% right_join(.,SR,by="CellID") -> SR
  roads <- rnaturalearth::ne_countries(scale = 110,returnclass = "sf") %>% sf::st_transform(., crs = CRS("+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m  +no_defs"))
  xlimits = c(-11286442,-6932762)
  ylimits = c(988516,3849386)
  SR %>% mutate(n_std = (n - min(n)) / (max(n) - min(n))) -> SR
  reduct_five %>% count(CellID) %>% left_join(.,reduct_five %>% distinct(CellID,.keep_all = T),by="CellID") -> counts
  g@coords %>% as_tibble() %>% mutate(CellID=1:nrow(.)) %>% right_join(.,counts,by="CellID") -> counts
  counts %>% mutate(n_std = (n - min(n)) / (max(n) - min(n))) -> counts
  
  histo <- reduct_five %>% group_by(year) %>% mutate(year=sub("Year_","",year)) %>% mutate(year = as.numeric(year)) %>%
    ggplot(aes(x=year)) + 
    geom_histogram(stat="count",binwidth = 1,fill=MetBrewer::met.brewer("Hiroshige",direction=1)[10]) +  theme(panel.background = element_blank(),panel.grid = element_blank(),axis.line = element_line(),legend.position = c(0.1,0.3),legend.background = element_rect(fill = NA), legend.key.size = unit(0.7, "cm"),legend.key.width = unit(0.6,"cm"), legend.text = element_text(family="EB Garamond",size=8),legend.title = element_text(family="EB Garamond"),legend.spacing = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=14,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),axis.text = element_text(family="EB Garamond"),axis.title = element_text(family="EB Garamond",size=11)
  ) + labs(y="",x="") +
    NULL
  
reduct_five %>% group_by(Accepted_Name) %>% mutate(year=sub("Year_","",year)) %>% mutate(year = as.numeric(year)) %>%
    ggplot(aes(x=Accepted_Name)) + 
    geom_point(stat="count",fill=MetBrewer::met.brewer("Hiroshige",direction=1)[10]) +  
  theme(panel.background = element_blank(),panel.grid = element_blank(),axis.line = element_line(),legend.position = c(0.1,0.3),legend.background = element_rect(fill = NA), legend.key.size = unit(0.7, "cm"),legend.key.width = unit(0.6,"cm"), legend.text = element_text(family="EB Garamond",size=8),legend.title = element_text(family="EB Garamond"),legend.spacing = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=14,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),axis.text = element_text(family="EB Garamond"),axis.title = element_text(family="EB Garamond",size=11)
    ) + labs(y="",x="") +
    NULL
  
mappy <- ggplot() + 
    geom_sf(data=roads,colour=NA,fill=adjustcolor(MetBrewer::met.brewer("Demuth",direction=-1)[1],alpha.f = 0.2),size=0.5) + xlim(xlimits) + ylim(ylimits) + 
    theme(panel.background = element_blank(),panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.border = element_rect(fill=NA),legend.position = c(0.9,0.35),legend.background = element_rect(fill = NA), legend.key.size = unit(0.7, "cm"),legend.key.width = unit(0.6,"cm"), legend.text = element_text(family="EB Garamond",size=8),legend.title = element_text(family="EB Garamond"),legend.spacing = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5)) +
    geom_tile(data = SR ,aes(x=x,y=y,fill=(n),color=(n))) +
    scale_fill_stepsn(colors=(MetBrewer::met.brewer("Tam",direction=1)),name="Species  \nrichness (log)",n.breaks=10,trans = "log",labels=function(x) sprintf("%.0f", (x))) +
    scale_color_stepsn(colors=(MetBrewer::met.brewer("Tam",direction=1)),name="Species  \nrichness (log)",n.breaks=10,trans = "log",labels=function(x) sprintf("%.0f", (x))) +
    labs(x="",y="",title="Occurrence records for cloud forests species", subtitle="1979 \U2012 2019") +
    #annotate("text",x=-11132762,y=1288516,label=paste0("More than 50 occurrence records: ",formatC(nrow(plus50), big.mark=",")," species"," (",formatC(sum(plus50$n), big.mark=","),")"),family="EB Garamond",size=4,hjust=0) +
    #annotate("text",x=-11132762,y=1188516,label=paste0("More than 100 occurrence records: ",formatC(nrow(plus100), big.mark=",")," species"," (",formatC(sum(plus100$n), big.mark=","),")"),family="EB Garamond",size=4,hjust=0) +
    #annotate("text",x=-11132762,y=1088516,label=paste0("More than 500 occurrence records: ",formatC(nrow(plus500), big.mark=",")," species"," (",formatC(sum(plus500$n), big.mark=","),")"),family="EB Garamond",size=4,hjust=0) +
    annotation_custom(ggplotGrob(histo), xmin = -11432762, xmax = -9632762, 
                      ymin = 1388516, ymax = 2388516) +
    #standardPrintOutput::watermark(show=T,lab="Ramirez-Barahona et al.") +
    NULL

ggsave(filename = here("plots_tables/Supplementary_figure_S1.pdf"),plot = mappy,height = 7,width = 8.5)
  
}

data_descriptors <-  function(data="output/Historical_data_FULL_LUL.v.1.R",data_spp = "output/Species_summaryRaw.RData" ){
load(here(data))### occurrences
load(data_spp)### summary (vs. nulls)
species_means[[1]]
species_means[[1]] %>% count(Accepted_Name)

species_means[[1]] %>% count(Accepted_Name) %>% pull(Accepted_Name) -> targets
reduct_five %>% group_by(Accepted_Name) %>% filter(Accepted_Name%in%all_of(targets)) %>% 
  summarize(Class=first(class),Family=first(family.x),Occurrences=n()) %>% 
  flextable::flextable() -> ft
flextable::colformat_num(x = ft,decimal.mark=",") %>% flextable::save_as_docx(path = here("plots_tables/Supplementary_Table_S1.docx"))

reduct_five %>% group_by(ECO_NAME) %>% filter(Accepted_Name%in%all_of(targets)) %>% 
  summarize(Biome=first(BIOME),Occurrences=n()) %>% rename("Ecoregion"=ECO_NAME) %>% flextable::flextable() -> bt
flextable::colformat_num(x = bt,decimal.mark=",") %>% flextable::save_as_docx(path = here("plots_tables/Supplementary_Table_S2.docx"))

panelA <- reduct_five %>% select(1:15) %>% filter(Accepted_Name%in%all_of(targets)) %>% 
  mutate(year=sub("Year_","",year)) %>%  mutate(AltCat = cut(Altitude,breaks=seq(0,4500,500),include.lowest=T,labels=seq(0,4500,500)[-1]),.before=Altitude) %>% 
  mutate(AltYear = cut(as.numeric(year),breaks=seq(1979,2019,5),labels=c("1979-1984","1985-1989","1990-1994","1995-1999","2000-2004","2004-2009","2010-2014","2015-219"),include.lowest=T),.before=Altitude) %>% 
  filter(!is.na(AltCat)) %>% count(AltYear) %>% 
  ggplot(aes(y=n,x=AltYear,fill=AltYear)) +
  geom_col() + 
  theme(legend.position = "none",panel.background = element_blank(),axis.line = element_line(),
        axis.title = element_text(family="EB Garamond"),
        axis.text = element_text(family="EB Garamond"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond")) +
  labs(y="Number of occurrences",x="")  +
  scale_fill_grey(guide=guide_legend(reverse = F,title.position="top", title.hjust = 0.5),name="Five-year period")  +
  geom_text(aes(label=scales::comma(n)), vjust=-.5, color="black", size=3,family="EB Garamond") +
  guides(fill="none")+
  NULL

panelB <- reduct_five %>% select(1:15) %>% filter(Accepted_Name%in%all_of(targets)) %>% 
  mutate(year=sub("Year_","",year)) %>%  mutate(AltCat = cut(Altitude,breaks=seq(0,4500,500),include.lowest=T,labels=seq(0,4500,500)[-1]),.before=Altitude) %>% 
  mutate(AltYear = cut(as.numeric(year),breaks=seq(1979,2019,5),labels=c("1979-1984","1985-1989","1990-1994","1995-1999","2000-2004","2004-2009","2010-2014","2015-219"),include.lowest=T),.before=Altitude) %>% filter(!is.na(AltCat)) %>% count(AltCat,AltYear) %>% group_by(AltCat) %>% mutate(Prop = n/sum(n)) %>% 
ggplot(aes(y=Prop,x=AltCat,fill=AltYear,group=fct_reorder(.f=AltCat, .x = -as.numeric(AltYear),.fun = max))) +
  geom_col() + 
  theme(legend.position = "bottom",panel.background = element_blank(),axis.line = element_line(),
        axis.title = element_text(family="EB Garamond"),
        axis.text = element_text(family="EB Garamond"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond")) +
  labs(y="Proportion of occurrences",x="Altitudinal belt (upper limit in meters)")  +
  scale_fill_manual(values= MetBrewer::met.brewer("Tam",n=9,direction = 1),
                    guide=guide_legend(reverse = F,title.position="top", title.hjust = 0.5),name="Five-year period") +
  NULL
panelB <- reduct_five %>% select(1:15) %>% filter(Accepted_Name%in%all_of(targets)) %>% 
  mutate(year=sub("Year_","",year)) %>%  mutate(AltCat = cut(Altitude,breaks=seq(0,4500,500),include.lowest=T,labels=seq(0,4500,500)[-1]),.before=Altitude) %>% 
  mutate(AltYear = cut(as.numeric(year),breaks=seq(1979,2019,5),labels=c("1979-1984","1985-1989","1990-1994","1995-1999","2000-2004","2004-2009","2010-2014","2015-219"),include.lowest=T),.before=Altitude) %>% filter(!is.na(AltCat)) %>% count(AltCat,AltYear) %>% group_by(AltYear) %>% mutate(Prop = n/sum(n)) %>% 
  ggplot(aes(y=Prop,x=AltYear,fill=AltCat,group=fct_reorder(.f=AltCat, .x = -as.numeric(AltCat),.fun = max))) +
  geom_col() + 
  theme(legend.position = "bottom",panel.background = element_blank(),axis.line = element_line(),
        axis.title = element_text(family="EB Garamond"),
        axis.text = element_text(family="EB Garamond"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond")) +
  labs(y="Proportion of records",x="")  +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  scale_fill_manual(values= MetBrewer::met.brewer("Hiroshige",n=9,direction = -1),
                    guide=guide_legend(reverse = F,nrow = 1,title.position="top", title.hjust = 0.5,label.position = "bottom"),name="Altitudinal belt (upper limit in meters)") +
  NULL


panelC <- reduct_five %>% select(1:15) %>% filter(Accepted_Name%in%all_of(targets)) %>% 
  mutate(year=sub("Year_","",year)) %>%  mutate(AltCat = cut(Altitude,breaks=seq(0,4500,500),include.lowest=T,labels=seq(0,4500,500)[-1]),.before=Altitude) %>% 
  mutate(AltYear = cut(as.numeric(year),breaks=seq(1979,2019,5),labels=c("1979-1984","1985-1989","1990-1994","1995-1999","2000-2004","2004-2009","2010-2014","2015-219"),include.lowest=T),.before=Altitude) %>% filter(!is.na(AltCat)) %>% group_by(AltCat,AltYear) %>% 
  summarise(sr = n_distinct(Accepted_Name),n=n()) %>% group_by(AltYear) %>% 
  mutate(Prop = sr/sum(sr)) %>% 
  ggplot(aes(y=Prop,x=AltYear,fill=AltCat,group=fct_reorder(.f=AltCat, .x = -as.numeric(AltCat),.fun = max))) +
  geom_col(position="stack") + 
  theme(legend.position = "bottom",panel.background = element_blank(),axis.line = element_line(),
        axis.title = element_text(family="EB Garamond"),
        axis.text = element_text(family="EB Garamond"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond")) +
  labs(y="Proportion of Species",x="")  +
  
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  scale_fill_manual(values= MetBrewer::met.brewer("Hiroshige",n=9,direction = -1),
                    guide=guide_legend(reverse = F,nrow = 1,title.position="top", title.hjust = 0.5,
                                       label.position = "bottom"),name="Altitudinal belt (upper limit in meters)") +
  NULL

layout <- "
AAAA
BBCC
"
library(patchwork)
pp <- wrap_plots(list(panelA,panelB,panelC)) +
  plot_layout(design = layout,byrow=T,nrow = 2,guides="collect",ncol=2,tag_level = 'new') + 
  plot_annotation(tag_levels = 'a',title="") & theme(legend.position = 'bottom',plot.tag=element_text(family="EB Garamond",size=15,face="bold"),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),legend.key.height=unit(0.2,"cm"),legend.text=element_text(size=8),
                                                     legend.title=element_text(size=10))

ggsave(filename = here("plots_tables//Supplementary_figure_S1.pdf"),plot = pp,height = 7,width = 8.5)
pp

}

LMM_elevation <- function(TRY="output/growth_try.Rdata",splines="output/Splines_v.2.Rdata",
                     start_year=1979,end_year=2010, outtable="plots_tables/Supplementary_tableS3.docx",
                     ){
load(here(splines))
load(here(TRY))
growth_try %<>% rename("binomial"=Accepted_Name) %>% rename("growth"=OrigValueStr) %>% mutate(growth=as.factor(growth))
levels(growth_try$growth) <- c("Aquatic","Epiphyte/Climber","Epiphyte/Climber","Herb","Shrub","Shrub","Tree","Tree")
growth_try %>% group_by(binomial) %>% nest() %>% rowwise() %>% mutate(LifeForm=(Mode(data$growth))) -> target_try
aver <- mean_alts %>% mutate(Type="OBS",.before=data) %>% select(2:3) %>% unnest(data)
aver %<>% ungroup() %>% 
  select(Accepted_Name,YearCat,RollAlt,Rolln,RollAlt_sd) %>% 
  mutate(YearCat=as.numeric(YearCat)) %>% 
  mutate(RollAlt_error = RollAlt_sd / sqrt(Rolln)) %>% 
  select(RollAlt,YearCat,Accepted_Name,RollAlt_error) %>% rename_with(~c("y","tme","ind","error"))
#aver %>% right_join(.,spp_base_alts,by="ind")
aver %<>% #right_join(.,spp_base_alts,by="ind") %>% 
  filter(tme <= end_year) %>% mutate(year=tme) %>% mutate(tme = tme - start_year)
otra <- aver %>% separate(ind,into = c("G","S","extra"),extra = "merge",sep=" ") %>% select(-extra) %>% unite("ind",G:S,sep="_") #%>% left_join(.,target_try %>% rename("ind"=binomial) %>% select(-data),by="ind") #%>% mutate(LifeForm = replace_na(LifeForm,"Unknown"))
otra %<>% mutate(tme1 = ifelse(tme > 11,11,tme)) %>% mutate(tme2 = tme - 10) %>% mutate(tme2=ifelse(tme2 < 0,0,tme2)) %>% mutate(tme2 = ifelse(tme2 > 11,11,tme2)) %>% mutate(tme3 = tme - 20) %>% mutate(tme3=ifelse(tme3<0,0,tme3))
### rescale variables
otra <- transform(otra,tme_sc=scale(tme),error_sc=scale(error),tme1_sc=scale(tme1),tme2_sc=scale(tme2),tme3_sc=scale(tme3)) %>% as_tibble()
otra %<>% group_by(ind) %>% mutate(Quantiles=ecdf(y)(y)) #%>% mutate(y_cens = ifelse(Quantiles<0.975,y,quantile(y,0.975))) %>% mutate(y_cens = ifelse(Quantiles>0.025,y,quantile(y,0.0025)))
otra <- otra %>% filter(year<=2010) #%>% mutate(Alt_cat = as.factor(Alt_cat))
to_merge <- otra %>% group_by(ind) %>% summarise(error = median(error,na.rm=TRUE),#LifeForm=first(LifeForm)
                                                 ) #%>% mutate(LifeForm=replace_na(LifeForm,replace = "Unknown"))
fit_error <- lmer(error ~ 1 + tme + (1|ind) + (0 + tme|ind),data = otra, control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4)))
fit_pre <- lmer(y ~ 1 + tme + (1|ind) + (0 + tme|ind),data = otra, control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4)))
pre <- broom.mixed::tidy(fit_pre)
#plot(fit_pre)
### error as covariate
fit0 <- lmer(y ~ 1 + tme + error + (1|ind) + (0 + tme|ind),data = otra, control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4)))
cov <- broom.mixed::tidy(fit0)

#summary(fit0)
### error as covariate + interaction
fit1 <- lmer(y ~ 1 + tme + tme:error + (1|ind) + (0 + tme|ind),data = otra, control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4)))
covinter <- broom.mixed::tidy(fit1)

anova(fit_pre,fit0,fit1)
modelos <- list("time"=fit_pre,"cov"=fit0,"interaction"=fit1)
save(modelos,file="output/LMM_models.RData")
broom.mixed::tidy(fit0) %>% mutate(across(4:6,function(x) round(x,2) )) %>% flextable::flextable() %>% flextable::save_as_docx(path = here(outtable))
broom.mixed::tidy(fit1) %>% mutate(across(4:6,function(x) round(x,2) )) %>% flextable::flextable() %>% flextable::save_as_docx(path = here( "plots_tables/Supplementary_tableS2.docx"))
esta <- coef(fit0)$ind %>% as_tibble(rownames = "ind") %>% mutate(Inter_Cat=cut(`(Intercept)`,breaks=seq(0,4000,200),include.lowest=T,labels=seq(0,4000,200)[-1] )) %>% 
  left_join(.,to_merge,by="ind")
e <- esta %>% select(tme,ind) %>% pivot_longer(cols = 1,names_to = "trend",values_to = "estimate" ) %>% mutate(trend = sub("tme","1979-2019",trend)) %>% left_join(.,esta %>% select(ind,Inter_Cat),by="ind")
model_data <- list(e,to_merge)
save(model_data,file=here("output/LMM_modelData.Rdata"))

ggpredict(fit1,terms=c("tme","error [5,10,20,40,60,80,100]"),type="random",interval="confidence") %>%
  ggplot(aes(x=x+1979,y=predicted,group=group,color=group))  + geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill=group),color="NA", alpha = 0.1) +
  scale_color_met_d("Hiroshige") +  
  scale_fill_met_d("Hiroshige") +
  theme(legend.position = "right",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line()) +
  labs(x="",y="Predicted general elevation (meters above sea level)")+
  NULL
ggsave(filename = "plots_tables/Supplementary_figure_S1.pdf",plot = last_plot())
ggpredict(fit0,terms=c("tme[0,8,16,24,31]"),type="random",interval="confidence")
ggpredict(fit0,terms=c("tme[0,8,16,24,31]"),type="random",interval="confidence") %>% as.data.frame() -> az
az %>% ggplot(aes(x=as.character(x),y=predicted)) + geom_point() +
  geom_point(aes(y=conf.low),color="red") +
  geom_point(aes(y=conf.high),color="blue")

}

summary_LMMnulls_elevation <- function(nulls="output/SplinesNULLS_v.2.Rdata",
                              output="interim/resampSummary_LMM.Rdata",
                              output_fig = "plots_tables/resampSummary_LMM.pdf",
                              which=1,n_sim=500,end_year=2019,start_year=1979){
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

plot_LMM <- function(){
  data="output/LMM_modelData.Rdata"
  downsample="interim/DownSummary_LMM.Rdata"
  resample = "interim/resampSummary_LMM.Rdata"
  load(here(data))
  load(here(downsample))
    down <- e_nulls
    load(here(resample))
    resamp <- e_nulls
    dd <- down %>% group_by(ind) %>% summarise(mean=mean(estimate),sd=sd(estimate)) %>% mutate(type="down")
    dd <- resamp %>% group_by(ind) %>% summarise(mean=mean(estimate),sd=sd(estimate)) %>% mutate(type="resamp") %>% rbind(.,dd)
    dd %<>% filter(ind!="Frangula capreifolia var. grandifolia (M.C.Johnst. & L.A.Johnst.) A.Pool") %>% filter(ind!="Styrax argenteus var. ramirezii (Greenm.) Gonsoulin") 
    dd %<>% separate(ind,into = c("G","S","extra"),extra = "merge",sep=" ") %>% select(-extra) %>% unite("ind",G:S,sep="_")
    dd %<>% right_join(.,model_data[[1]],by="ind") %>% mutate(SES = (estimate - mean) / sd,ES=(estimate - mean))
    cat("Summary downsampling (only positives)","\n")
    dd %>% filter(type=="down") %>% filter(estimate > 0,SES > 0) %>% filter(SES > 1.28) %>% 
      summarise(mean(SES),mean(ES),sd(SES),sd(ES)) %>% print
    cat("Summary resampling (only positives)","\n")
    dd %>% filter(type=="resamp") %>% filter(estimate > 0,SES > 0) %>% filter(SES > 1.96) %>% 
      summarise(mean(SES),mean(ES),sd(SES),sd(ES)) %>% print
    
    dd %>% filter(type=="resamp") %>% select(-trend,-Inter_Cat) %>% setNames(paste0('RSamp_', names(.))) %>% 
      rename(ind=RSamp_ind) %>% 
      left_join(.,dd %>% filter(type=="down") %>% setNames(paste0('DSamp_', names(.))) %>% 
                  rename(ind=DSamp_ind),by="ind") %>% 
      mutate(RSamp_SES_cat = cut(RSamp_SES,breaks=c(-6,-1.96,-1.28,0,1.28,1.96,6),labels=c("Strong_Neg","Moderate_Neg","Weak_Neg","Weak_Pos","Moderate_Pos","Strong_Pos"),include.lowest=TRUE)) %>% 
      mutate(DSamp_SES_cat = cut(DSamp_SES,breaks=c(-6,-1.96,-1.28,0,1.28,1.96,6),labels=c("Strong_Neg","Moderate_Neg","Weak_Neg","Weak_Pos","Moderate_Pos","Strong_Pos"),include.lowest=TRUE)) %>% #unite("CAT",RSamp_SES_cat:DSamp_SES_cat) %>%
      mutate(CAT = ifelse(RSamp_estimate <= 0,"Negative","Positive")) -> nulas

  esta <- model_data[[1]] %>% left_join(.,model_data[[2]],by="ind")
data="output/Historical_data_FULL_GLASS.v.1.R"
load(here(data))
reduct_five %<>% separate(Accepted_Name,into = c("G","S","extra"),extra = "merge",sep=" ") %>% select(-extra) %>% unite("ind",G:S,sep="_") %>% right_join(.,esta,by="ind")

base_eco <- ecolandia %>% select(CellID:ind,year,starts_with("min_tas_"),starts_with("mean_tas_"),starts_with("max_tas_")) %>%  mutate(year=as.numeric(sub("Year_","",year))) %>% 
  filter(year<=1990) %>% 
  distinct(CellID,ind,.keep_all = TRUE) %>% 
  rowwise() %>% mutate(baseline_tas = mean(c_across(mean_tas_1979:mean_tas_1990),na.rm=TRUE),.before=year) %>% mutate(baseline_tasmin = mean(c_across(min_tas_1979:min_tas_1990),na.rm=TRUE),.before=year) %>% mutate(baseline_tasmax = mean(c_across(max_tas_1979:max_tas_1990),na.rm=TRUE),.before=year)

last_eco <- ecolandia %>% select(CellID:ind,year,starts_with("min_tas_"),starts_with("mean_tas_"),starts_with("max_tas_")) %>% mutate(year=as.numeric(sub("Year_","",year))) %>% 
 filter(year<=2010 & year >= 2000) %>% 
  distinct(CellID,ind,.keep_all = TRUE) %>% 
  rowwise() %>% mutate(last_tas = mean(c_across(mean_tas_2000:mean_tas_2010),na.rm=TRUE),.before=year) %>% mutate(last_tasmin = mean(c_across(min_tas_2000:min_tas_2010),na.rm=TRUE),.before=year) %>% mutate(last_tasmax = mean(c_across(max_tas_2000:max_tas_2010),na.rm=TRUE),.before=year)


base_eco %<>% ungroup() %>% select(1:6) %>% mutate(across(starts_with("baseline"), function(x) x*.1 -  273.15)) %>% 
  mutate(TSI = (abs(baseline_tas - baseline_tasmin)) / (abs(baseline_tas - baseline_tasmax))) %>% 
group_by(ind) %>% 
  summarise(across(baseline_tas:TSI,function(x) mean(x,na.rm=TRUE) ))

last_eco %<>% ungroup() %>% select(1:6) %>% mutate(across(starts_with("last"), function(x) x*.1 -  273.15)) %>% 
  mutate(TSI = (abs(last_tas - last_tasmin)) / (abs(last_tas - last_tasmax))) %>% 
  group_by(ind) %>% 
  summarise(across(last_tas:TSI,function(x) mean(x,na.rm=TRUE) ))

to_plot_tas <-  base_eco %>% rename(TSI_base = TSI) %>% left_join(.,last_eco %>% rename(TSI_last = TSI),by="ind") %>%
    left_join(., nulas %>% select(ind,RSamp_SES,RSamp_SES_cat,DSamp_SES_cat,CAT,RSamp_SES,DSamp_SES,RSamp_ES,DSamp_ES,DSamp_Inter_Cat),by="ind")

thres_alt = 1500
by_alt <- to_plot_tas %>%
 mutate(Alt_Cat = fct_rev(factor(ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) <= thres_alt ,"Lowlands", ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) > thres_alt ,"Highlands","Highlands") )))) %>% 
  ggplot(aes(x=RSamp_SES_cat,y = baseline_tas,group=Alt_Cat)) +
  ggdist::stat_interval(side="bottom",width=0,.width = seq(0,1,0.1),interval_size=3,
                        aes(color=after_stat(level),alpha=after_stat(level))) +
  ggdist::stat_dotsinterval(size=1,justification=-0.15,interval_alpha=0) +
  theme(legend.position = "none") +
  scale_color_manual(values= colorspace::darken(MetBrewer::met.brewer("OKeeffe2",n=11,direction=1),amount=0.3), guide=guide_legend(reverse = T,ncol = 2,title.position="top", title.hjust = 0, label.position = "top",direction = "vertical"),name="") +
  scale_alpha_manual(values=seq(0.2,1.0,.07)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title = element_text(size=14),
        axis.line = element_line(),
        text=element_text(family="EB Garamond")) +
  geom_text(data = to_plot_tas %>%  mutate(Alt_Cat = fct_rev(factor(ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) <= thres_alt ,"Lowlands",ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) > thres_alt ,"Highlands","NA") ))))  %>% group_by(RSamp_SES_cat,Alt_Cat) %>% summarise(med=median(baseline_tas)),family="EB Garamond",fontface="bold",size=5,
            aes(x=RSamp_SES_cat,y=med-7,group=Alt_Cat,label=paste0(round(med,2)))) +
  geom_text(data= to_plot_tas %>%  mutate(Alt_Cat = fct_rev(factor(ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) <= thres_alt ,"Lowlands",ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) > thres_alt ,"Highlands","NA") ))))  %>% 
              group_by(RSamp_SES_cat,Alt_Cat) %>% mutate(med=median(baseline_tas)) %>% 
              summarize(med=first(med),n=n()),family="EB Garamond",
            aes(x=RSamp_SES_cat,y=med-8,group=Alt_Cat,label=paste0("n=",round(n,2)))) +
  labs(x="",y="Baseline mean temperature (ºC)") +
  stat_summary(geom="crossbar",fun=function(x) quantile(x,probs=c(0.5)),width=0.2,size=0.4) +
  stat_summary(geom="crossbar",fun=function(x) quantile(x,probs=c(0.75,0.25)),width=0.18,size=0.2) +
  facet_wrap(~Alt_Cat) +
  NULL

general <- to_plot_tas %>% 
  ggplot(aes(x=RSamp_SES_cat,y = baseline_tas)) +
  ggdist::stat_interval(side="bottom",width=0,.width = seq(0,1,0.1),interval_size=3,
                        aes(color=after_stat(level),alpha=after_stat(level))) +
  ggdist::stat_dotsinterval(size=1,justification=-0.15,interval_alpha=0) +
  theme(legend.position = "none") +
  scale_color_manual(values= colorspace::darken(MetBrewer::met.brewer("OKeeffe2",n=11,direction=1),amount=0.3), guide=guide_legend(reverse = T,ncol = 2,title.position="top", title.hjust = 0, label.position = "top",direction = "vertical"),name="") +
  scale_alpha_manual(values=seq(0.2,1.0,.07)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title = element_text(size=14),
        axis.line = element_line(),
        text=element_text(family="EB Garamond")) +
  geom_text(data = to_plot_tas  %>% group_by(RSamp_SES_cat) %>% summarise(med=median(baseline_tas)),family="EB Garamond",fontface="bold",size=5,
            aes(x=RSamp_SES_cat,y=med-9,label=paste0(round(med,2)))) +
  geom_text(data= to_plot_tas  %>%
              group_by(RSamp_SES_cat) %>% mutate(med=median(baseline_tas)) %>% 
              summarize(med=first(med),n=n()),family="EB Garamond",
            aes(x=RSamp_SES_cat,y=med-10,label=paste0("n=",round(n,2)))) +
  labs(x="",y="Baseline mean temperature (ºC)") +
  stat_summary(geom="crossbar",fun=function(x) quantile(x,probs=c(0.5)),width=0.2,size=0.4) +
  stat_summary(geom="crossbar",fun=function(x) quantile(x,probs=c(0.75,0.25)),width=0.18,size=0.2) +
  NULL

layout <- "
AA
BB
"
pp<-wrap_plots(list(general,by_alt)) +
  plot_layout(design = layout,byrow=T,nrow = 2,guides="keep",ncol=2,tag_level = 'new') + 
  plot_annotation(tag_levels = 'a',title="") & theme(legend.position = 'none',plot.tag=element_text(family="EB Garamond",size=15,face="bold"),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))
pp

ggsave(filename=here("plots_tables/Shifts_by_baseline.pdf"),plot=pp)

thres_alt = 1500
to_plot_tas %>% mutate(Alt_Cat = fct_rev(factor(ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) <= thres_alt ,"Lowlands",ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) > thres_alt ,"Highlands","NA") )))) %>% group_by(Alt_Cat,RSamp_SES_cat) %>% 
  summarise(Median_temp=median(baseline_tas),n_species=n()) %>% mutate(Threshold = thres_alt,.before=1) -> a

to_plot_tas %>% mutate(Alt_Cat = fct_rev(factor(ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) <= 1200 ,"Lowlands",ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) > 1200 ,"Highlands","NA") )))) %>% group_by(Alt_Cat,RSamp_SES_cat) %>% 
  summarise(Median_temp=median(baseline_tas),n_species=n()) %>% mutate(Threshold = 1200,.before=1) -> b

to_plot_tas %>% mutate(Alt_Cat = fct_rev(factor(ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) <= 1000 ,"Lowlands",ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) > 1000 ,"Highlands","NA") )))) %>% group_by(Alt_Cat,RSamp_SES_cat) %>% 
  summarise(Median_temp=median(baseline_tas),n_species=n()) %>% mutate(Threshold = 1000,.before=1) -> c
bind_rows(a,b,c) %>% flextable::flextable() %>% flextable::save_as_docx(path = here( "plots_tables/Supplementary_tableS3.docx"))

to_plot_tas %>% select(ind,RSamp_SES_cat,DSamp_SES_cat,starts_with("TSI")) %>% 
  pivot_longer(-c(1:3),names_to = "Variable",values_to = "Values") %>% 
  mutate(RSamp_SES_cat = factor(RSamp_SES_cat,levels = c(
    "Strong_Neg", "Moderate_Neg","Weak_Neg"  ,"Weak_Pos", "Moderate_Pos","Strong_Pos" ))) %>% 
  ggplot(aes(y=Values,x=RSamp_SES_cat,group=Variable)) +
      ggdist::stat_interval(interval_size=10,.width = seq(0,1,0.1),aes(color=after_stat(level),interval_alpha=after_stat(level)),position = "dodge",justification = -1) +
      theme(legend.position = "none") +
      scale_color_manual(values= colorspace::darken(MetBrewer::met.brewer("OKeeffe2",n=11,direction=1),amount=0.1), guide=guide_legend(reverse = T,ncol = 2,title.position="top", title.hjust = 0, label.position = "top",direction = "vertical"),name="") +
  labs(y="TSI", x = "Elevation-shift SES") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(),
        axis.line = element_line()) +
  ggdist::stat_pointinterval(size=1,justification=-0.15,slab_colour="grey30",shape=19,point_interval="median_qi",position = "dodge") +
      NULL



### 
data="output/Historical_data_FULL_GLASS.v.1.R"
load(here(data))
reduct_five %<>% separate(Accepted_Name,into = c("G","S","extra"),extra = "merge",sep=" ") %>% select(-extra) %>% unite("ind",G:S,sep="_") %>% right_join(.,esta,by="ind")


reduct_five %>% select(CellID:ind,year,#starts_with("min_cmi_"),
                     starts_with("mean_tas_")
                    # ,starts_with("max_cmi_")
                     ) %>%
  mutate(year=as.numeric(sub("Year_","",year))) %>% select(-Altitude,-year,-ind) %>% distinct(CellID,.keep_all = TRUE) %>% 
  pivot_longer(-c(1), names_to = c("Variable","year"),names_sep = "_pr_",values_to = "Values") %>% 
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
  filter(DSamp_SES_cat %in% c("Robust")) %>% 
lmer(Climate ~ 1 + tme + tme:RSamp_SES_cat + (1|ind), data=., control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5))) -> fit_aver

broom.mixed::tidy(fit_aver)
#update(fit_aver)
ggpredict(fit_aver, type = "random",terms=c("tme [all]","RSamp_SES_cat"),interval = "confidence") %>% tibble() -> predict_fit
predict_fit %>% 
  ggplot(aes(x=x,y=predicted,color=group,fill=group)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low,ymax=conf.high),alpha=0.1,color="NA") +
  scale_color_met_d(name="Hiroshige",direction = -1)+
  scale_fill_met_d(name="Hiroshige",direction = -1) +
  #geom_line(data=ggpredict(fit_sum, type = "random",terms=c("tme [all]"),interval = "confidence") %>% tibble(),aes(x=x,y=predicted),color="black",inherit.aes = FALSE) +
  theme(legend.position = "bottom",panel.background = element_blank(),axis.line = element_line(),legend.background = element_blank(),
        axis.title = element_text(family="EB Garamond"),
        axis.text = element_text(family="EB Garamond"),
        #axis.text.x = element_blank(),
        legend.spacing.y = unit(c(1.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond")) +
  labs(x="Time",y="Predicted mean tas") +
  NULL
  

##### asdasd
  





















to_plot_pr <- base_eco %>% ungroup() %>% select(1:6) %>% mutate(across(starts_with("baseline"), function(x) x*.1 -  273.15)) %>% 
  mutate(TSI_base = (abs(baseline_pr - baseline_prmin)) / (abs(baseline_pr - baseline_prmax))) %>% left_join(., 
  last_eco %>% ungroup() %>% select(1:6) %>% mutate(across(starts_with("last"), function(x) x*.1 -  273.15)) %>% 
  mutate(TSI_last = (abs(last_pr - last_prmin)) / (abs(last_pr - last_prmax))),by=c("CellID","ind") ) %>% 
  mutate(TSI_delta = TSI_last - TSI_base) %>% select(-starts_with("Altitude")) %>% 
  group_by(ind) %>% 
  summarise(across(baseline_pr:TSI_delta,function(x) mean(x,na.rm=TRUE) )) %>% 
  left_join(., nulas %>% select(ind,RSamp_SES,RSamp_SES_cat,DSamp_SES_cat,CAT,RSamp_SES,DSamp_SES,RSamp_ES,DSamp_ES,DSamp_Inter_Cat),by="ind")

  
by_alt <- to_plot_pr %>%
  mutate(Alt_Cat = fct_rev(factor(ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) <= thres_alt ,"Lowlands", ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) > thres_alt ,"Highlands","Highlands") )))) %>% 
  ggplot(aes(x=RSamp_SES_cat,y = TSI_base,group=Alt_Cat)) +
  ggdist::stat_interval(side="bottom",width=0,.width = seq(0,1,0.1),interval_size=3,
                        aes(color=after_stat(level),alpha=after_stat(level))) +
  ggdist::stat_dotsinterval(size=1,justification=-0.15,interval_alpha=0) +
  theme(legend.position = "none") +
  scale_color_manual(values= colorspace::darken(MetBrewer::met.brewer("OKeeffe2",n=11,direction=1),amount=0.3), guide=guide_legend(reverse = T,ncol = 2,title.position="top", title.hjust = 0, label.position = "top",direction = "vertical"),name="") +
  scale_alpha_manual(values=seq(0.2,1.0,.07)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title = element_text(size=14),
        axis.line = element_line(),
        text=element_text(family="EB Garamond")) +
  geom_text(data = to_plot_pr %>%  mutate(Alt_Cat = fct_rev(factor(ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) <= thres_alt ,"Lowlands",ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) > thres_alt ,"Highlands","NA") ))))  %>% group_by(RSamp_SES_cat,Alt_Cat) %>% summarise(med=median(TSI_base)),family="EB Garamond",fontface="bold",size=5,
            aes(x=RSamp_SES_cat,y=med*.4,group=Alt_Cat,label=paste0(round(med,2)))) +
  geom_text(data= to_plot_pr %>%  mutate(Alt_Cat = fct_rev(factor(ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) <= thres_alt ,"Lowlands",ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) > thres_alt ,"Highlands","NA") ))))  %>% 
              group_by(RSamp_SES_cat,Alt_Cat) %>% mutate(med=median(TSI_base)) %>% 
              summarize(med=first(med),n=n()),family="EB Garamond",
            aes(x=RSamp_SES_cat,y=med*0.6,group=Alt_Cat,label=paste0("n=",round(n,2)))) +
  labs(x="",y="Baseline mean precipitation (mm)") +
  stat_summary(geom="crossbar",fun=function(x) quantile(x,probs=c(0.5)),width=0.2,size=0.4) +
  stat_summary(geom="crossbar",fun=function(x) quantile(x,probs=c(0.75,0.25)),width=0.18,size=0.2) +
  facet_wrap(~Alt_Cat) +
  NULL

general <- to_plot_pr %>% 
  ggplot(aes(x=RSamp_SES_cat,y = TSI_last)) +
  ggdist::stat_interval(side="bottom",width=0,.width = seq(0,1,0.1),interval_size=3,
                        aes(color=after_stat(level),alpha=after_stat(level))) +
  ggdist::stat_dotsinterval(size=1,justification=-0.15,interval_alpha=0) +
  theme(legend.position = "none") +
  scale_color_manual(values= colorspace::darken(MetBrewer::met.brewer("OKeeffe2",n=11,direction=1),amount=0.3), guide=guide_legend(reverse = T,ncol = 2,title.position="top", title.hjust = 0, label.position = "top",direction = "vertical"),name="") +
  scale_alpha_manual(values=seq(0.2,1.0,.07)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title = element_text(size=14),
        axis.line = element_line(),
        text=element_text(family="EB Garamond")) +
  geom_text(data = to_plot_pr  %>% group_by(RSamp_SES_cat) %>% summarise(med=median(TSI_last)),family="EB Garamond",fontface="bold",size=5,
            aes(x=RSamp_SES_cat,y=med*0.4,label=paste0(round(med,2)))) +
  geom_text(data= to_plot_pr  %>%
              group_by(RSamp_SES_cat) %>% mutate(med=median(TSI_last)) %>% 
              summarize(med=first(med),n=n()),family="EB Garamond",
            aes(x=RSamp_SES_cat,y=med*0.6,label=paste0("n=",round(n,2)))) +
  labs(x="",y="Baseline mean precipitation (mm)") +
  stat_summary(geom="crossbar",fun=function(x) quantile(x,probs=c(0.5)),width=0.2,size=0.4) +
  stat_summary(geom="crossbar",fun=function(x) quantile(x,probs=c(0.75,0.25)),width=0.18,size=0.2) +
  NULL

layout <- "
AA
BB
"
pp<-wrap_plots(list(general,by_alt)) +
  plot_layout(design = layout,byrow=T,nrow = 2,guides="keep",ncol=2,tag_level = 'new') + 
  plot_annotation(tag_levels = 'a',title="") & theme(legend.position = 'none',plot.tag=element_text(family="EB Garamond",size=15,face="bold"),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))
pp
ggsave(filename=here("plots_tables/Shifts_by_baseline_prec.pdf"),plot=pp)

thres_alt = 1500
to_plot_pr %>% mutate(Alt_Cat = fct_rev(factor(ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) <= thres_alt ,"Lowlands",ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) > thres_alt ,"Highlands","NA") )))) %>% group_by(Alt_Cat,RSamp_SES_cat) %>% 
  summarise(Median_prec=median(baseline_pr),n_species=n()) %>% mutate(Threshold = thres_alt,.before=1) -> a

to_plot_pr %>% mutate(Alt_Cat = fct_rev(factor(ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) <= 1200 ,"Lowlands",ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) > 1200 ,"Highlands","NA") )))) %>% group_by(Alt_Cat,RSamp_SES_cat) %>% 
  summarise(Median_prec=median(baseline_pr),n_species=n()) %>% mutate(Threshold = 1200,.before=1) -> b

to_plot_pr %>% mutate(Alt_Cat = fct_rev(factor(ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) <= 1000 ,"Lowlands",ifelse(as.numeric(as.vector(DSamp_Inter_Cat)) > 1000 ,"Highlands","NA") )))) %>% group_by(Alt_Cat,RSamp_SES_cat) %>% 
  summarise(Median_prec=median(baseline_pr),n_species=n()) %>% mutate(Threshold = 1000,.before=1) -> c
bind_rows(a,b,c) %>% flextable::flextable() %>% flextable::save_as_docx(path = here( "plots_tables/Supplementary_tableS3_prec.docx"))

to_plot_pr %>% select(ind,RSamp_SES_cat,DSamp_SES_cat,starts_with("TSI")) %>% 
  pivot_longer(-c(1:3),names_to = "Variable",values_to = "Values") %>% 
  mutate(RSamp_SES_cat = factor(RSamp_SES_cat,levels = c(
    "Strong_Neg", "Moderate_Neg","Weak_Neg"  ,"Weak_Pos", "Moderate_Pos","Strong_Pos" ))) %>% 
  mutate(Variable = factor(Variable,levels = c(
    "TSI_base", "TSI_last","TSI_delta"))) %>% 
  filter(Variable=="TSI_delta") %>% 
  ggplot(aes(y=Values,x=RSamp_SES_cat,group=Variable)) +
  ggdist::stat_interval(interval_size=10,.width = seq(0,1,0.1),aes(color=Variable,interval_alpha=after_stat(level)),position = "dodge",justification = -1) +
  theme(legend.position = "none") +
  #scale_color_manual(values= colorspace::darken(MetBrewer::met.brewer("VanGogh3",n=11,direction=1),amount=0.1), guide=guide_legend(reverse = T,ncol = 2,title.position="top", title.hjust = 0, label.position = "top",direction = "vertical"),name="") +
  labs(y="PSI", x = "Elevation-shift SES") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(),
        axis.line = element_line()) +
  ggdist::stat_pointinterval(size=1,justification=-0.15,slab_colour="grey30",shape=19,point_interval="median_qi",position = "dodge") +
  ggdist::stat_dist_dots() +
  geom_text(data =to_plot_pr %>% select(ind,RSamp_SES_cat,DSamp_SES_cat,starts_with("TSI")) %>% 
              mutate(TSI_delta = TSI_last - TSI_base) %>% 
              pivot_longer(-c(1:3),names_to = "Variable",values_to = "Values") %>% 
              mutate(RSamp_SES_cat = factor(RSamp_SES_cat,levels = c(
                "Strong_Neg", "Moderate_Neg","Weak_Neg"  ,"Weak_Pos", "Moderate_Pos","Strong_Pos" ))) %>% 
              mutate(Variable = factor(Variable,levels = c(
                "TSI_base", "TSI_last","TSI_delta"))) %>% 
              filter(Variable=="TSI_delta") %>% summarise(median = median(Values,na.rm=TRUE),.by = RSamp_SES_cat)
              ,family="EB Garamond",fontface="bold",size=5,inherit.aes = FALSE,
            aes(x=RSamp_SES_cat,y=median-.05,label=paste0(round(median,4))))+
  geom_hline(yintercept = 0) +
  NULL


esta_dos <- ecolandia %>% group_by(ind) %>% summarise(baseline_pr = mean(baseline_pr,na.rm=TRUE), baseline_tas = mean(baseline_tas,na.rm=TRUE)) %>% 
    left_join(esta,.,by="ind")

nulas %>% select(ind,RSamp_SES_cat,DSamp_SES_cat,CAT,RSamp_SES,DSamp_SES,RSamp_ES,DSamp_ES) %>% left_join(esta_dos,.,by="ind") %>% 
  #filter(RSamp_SES_cat %in% c("Strong_Pos","Moderate_Pos")) %>%
#filter(CAT=="Positive") %>% 
    ggplot(aes(y = RSamp_ES, x = baseline_tas,group=CAT)) +
    geom_point() +
    scale_fill_manual(values= colorspace::darken(MetBrewer::met.brewer("Hiroshige",n=9,direction=1)[c(9:5)],amount=0.1), guide=guide_legend(reverse = T,ncol = 1,title.position="top", title.hjust = 0, label.position = "top"),name="") +
    geom_smooth(method="lm") +
      NULL

model_data
  panelC <- ggplot() +
  #geom_segment(aes(y=6,yend=1,x=-0.01,xend=-0.01),color="grey30") +
  ggdist::stat_dots(data= esta,aes(x = estimate, y = LifeForm,fill=LifeForm,color=LifeForm),scale=0.9) +
  scale_fill_manual(values= colorspace::darken(MetBrewer::met.brewer("Hiroshige",n=9,direction=1)[c(9:5)],amount=0.1), guide=guide_legend(reverse = T,ncol = 1,title.position="top", title.hjust = 0, label.position = "top"),name="") +
  geom_text(inherit.aes = T,data = e %>% left_join(.,to_merge,by="ind") %>% filter(estimate >=0) %>% count(LifeForm),aes(label=paste0(n," species"),x=c(20,20,20,20,20),y=LifeForm), position=position_nudge(y=c(0.5,0.5,0.5,0.5,0.5)), size=3,family="EB Garamond") +
  geom_text(inherit.aes = T,data = e %>% left_join(.,to_merge,by="ind") %>% count(LifeForm),aes(label=paste0("\U0028",n,"\U0029"),x = 25,y=LifeForm), position=position_nudge(y=c(0.5,0.5,0.5,0.5,0.5)), size=3,family="EB Garamond") +
theme(legend.position = "none",panel.background = element_blank(),axis.line = element_line(),legend.background = element_blank(),legend.key.height = unit(0.01,"cm"),legend.key.width = unit(0.01,"cm"),
        axis.title.y = element_blank(),axis.title = element_text(family="EB Garamond"),
        axis.text = element_text(family="EB Garamond"),
        axis.text.y =  element_text(size=12,family="EB Garamond",margin=margin(t=-30,r=-15),vjust=-3),
        axis.ticks.y = element_blank(),
       axis.line.y = element_blank(),
        legend.spacing.y = unit(c(1.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond")) +
  labs(x=bquote("Elevation shift in meters"^-year)) +
  scale_y_discrete(expand = expansion(add=c(0.1, 1.5))) +
  scale_x_continuous(limits = c(-22,30)) +
  geom_segment(aes(y=c(1:5),yend=c(1:5),x=-20,xend=23),color="grey30",size=0.1) +
  geom_text(aes(label=c("downslope","upslope"),x = c(-14,10),y = 6.32),color="grey30", size=6,family="EB Garamond",fontface="bold",nudge_x=c(2,0)) +
  geom_segment(aes(y=6,yend=6,x=0,xend=-15),color="grey30",arrow = arrow(angle = 20,length=unit(0.2,"cm"))) +
  geom_segment(aes(y=6,yend=6,x=0,xend=15),color="grey30",arrow = arrow(angle = 20,length=unit(0.2,"cm"))) +
  scale_color_manual(values= colorspace::darken(MetBrewer::met.brewer("Hiroshige",n=9,direction=1)[c(9:5)],amount=0.2), guide=guide_legend(reverse = T,ncol = 1,title.position="top", title.hjust = 0, label.position = "top"),name="")  +
  NULL

panelA <- ggplot() +
  geom_segment(aes(y=1,yend=0.0,x=-0.01,xend=-0.01),color="grey30") +
  ggdist::stat_dots(data= esta,aes(x = estimate,fill=stat(x <  0),color=stat(x < 0)),size=0.5) + scale_fill_manual(values= colorspace::darken(MetBrewer::met.brewer("Hiroshige",n=9,direction=1)[c(1,2)],amount=0.1),guide=guide_legend(reverse = T,ncol = 1,title.position="top", title.hjust = 0, label.position = "top"),name="") +
  scale_color_manual(values= colorspace::darken(MetBrewer::met.brewer("Hiroshige",n=9,direction=1)[c(1,2)],amount=0.2),
                    guide=guide_legend(reverse = T,ncol = 1,title.position="top", title.hjust = 0, label.position = "top"),name="") +
  geom_text(inherit.aes = T,data = e %>% left_join(.,to_merge,by="ind") %>% filter(estimate > 0) %>% summarise(n=n()),
            aes(label=paste0(n," species"),x=c(20),y=0.5), 
            position=position_nudge(y=c(0)), size=3,family="EB Garamond") +
  geom_text(inherit.aes = T,data = e %>% left_join(.,to_merge,by="ind")%>% filter(estimate <= 0) %>% summarise(n=n()),
            aes(label=paste0(n," species"),x = c(-15),y=0.5), 
            position=position_nudge(y=c(0)), size=3,family="EB Garamond") +
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(),
        legend.background = element_blank(),
        legend.key.height = unit(0.01,"cm"),
        legend.key.width = unit(0.01,"cm"),
        axis.title = element_text(family="EB Garamond"),
        axis.text = element_text(family="EB Garamond"),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.line.y=element_blank(),
        legend.spacing.y = unit(c(1.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond")) +
labs(x=bquote("Elevation shift in meters"^-year),y="") +
  scale_x_continuous(limits = c(-25,30)) +
  geom_text(aes(label=c("downslope","upslope"),x = c(-14,10),y = 1.1),color="grey30", 
            size=6,family="EB Garamond",fontface="bold",nudge_x=c(2,0)) +
  geom_segment(aes(y=1,yend=1,x=0,xend=-15),color="grey30",arrow = arrow(angle = 20,length=unit(0.2,"cm"))) +
  geom_segment(aes(y=1,yend=1,x=0,xend=15),color="grey30",arrow = arrow(angle = 20,length=unit(0.2,"cm"))) +
  NULL

panel_B.2 <- ggplot() +
  geom_segment(aes(y=16,yend=1,x=-0.01,xend=-0.01),color="grey30") +
  ggdist::stat_dots(data=esta,aes(x = estimate, y = Inter_Cat,fill=Inter_Cat,color=Inter_Cat),scale=0.9) +
  scale_fill_manual(values= colorspace::darken(MetBrewer::met.brewer("Hiroshige",n=15,direction=1),amount=0.1),
                    guide=guide_legend(reverse = T,ncol = 1,title.position="top", title.hjust = 0, label.position = "top"),name="") +
  geom_text(inherit.aes = T,data = e %>% left_join(.,to_merge,by="ind") %>% filter(estimate >=0) %>% count(Inter_Cat),
            aes(label=paste0(n," species"),x=rep(25,1),y=Inter_Cat), 
            position=position_nudge(y=rep(0.5,1)), size=3,family="EB Garamond") +
  geom_text(inherit.aes = T,data = e %>% left_join(.,to_merge,by="ind") %>% count(Inter_Cat),
            aes(label=paste0("\U0028",n,"\U0029"),x = rep(29,1),y=Inter_Cat), 
            position=position_nudge(y=rep(0.5,1)), size=3,family="EB Garamond") +
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(),
        legend.background = element_blank(),
        legend.key.height = unit(0.01,"cm"),legend.key.width = unit(0.01,"cm"),
        axis.title.y = element_blank(),axis.title = element_text(family="EB Garamond"),
        axis.text = element_text(family="EB Garamond"),
        axis.text.y =  element_text(size=12,family="EB Garamond",margin=margin(t=-30,r=-15),vjust=-1),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.spacing.y = unit(c(1.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond")) +
  labs(x=bquote("Elevation shift in meters"^-year)) +
  scale_y_discrete(expand = expansion(add=c(0.1, 1.5))) +
  scale_x_continuous(limits = c(-20,30)) +
  geom_segment(aes(y=c(1:15),yend=c(1:15),x=-20,xend=23),color="grey30",size=0.1) +
  geom_text(aes(label=c("downslope","upslope"),x = c(-14,10),y = 16.32),color="grey30", 
            size=6,family="EB Garamond",fontface="bold",nudge_x=c(2,0))+
  geom_segment(aes(y=16,yend=16,x=0,xend=-15),color="grey30",arrow = arrow(angle = 20,length=unit(0.2,"cm"))) +
  geom_segment(aes(y=16,yend=16,x=0,xend=15),color="grey30",arrow = arrow(angle = 20,length=unit(0.2,"cm"))) +
  scale_color_manual(values= colorspace::darken(MetBrewer::met.brewer("Hiroshige",n=15,direction=1),amount=0.2),
                     guide=guide_legend(reverse = T,ncol = 1,title.position="top", title.hjust = 0, label.position = "top"),name="")  +
  NULL

layout <- "
AACC
AACC
"

pp<-wrap_plots(list(panelA,panel_B.2)) +
  plot_layout(design = layout,byrow=T,nrow = 2,guides="keep",ncol=2,tag_level = 'new') + 
  plot_annotation(tag_levels = 'a',title="") & theme(legend.position = 'none',plot.tag=element_text(family="EB Garamond",size=15,face="bold"),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))
pp
ggsave(filename = here("plots_tables/figure_1.pdf"),plot = pp,height = 7,width = 8.5)
ggsave(filename = here("plots_tables/Supplementary_figure_S2.pdf"),plot = panelC,height = 7,width = 8.5)

  }


plot_null <- function(data="output/LMM_modelData.Rdata",
                      downsample="interim/DownSummary_LMM.Rdata",
                      resample= "interim/resampSummary_LMM.Rdata"){
load(here(data))
load(here(downsample))
down <- e_nulls
load(here(resample))
resamp <- e_nulls
dd <- down %>% group_by(ind) %>% summarise(mean=mean(estimate),sd=sd(estimate)) %>% mutate(type="down")
dd <- resamp %>% group_by(ind) %>% summarise(mean=mean(estimate),sd=sd(estimate)) %>% mutate(type="resamp") %>% rbind(.,dd)
dd %<>% filter(ind!="Frangula capreifolia var. grandifolia (M.C.Johnst. & L.A.Johnst.) A.Pool") %>% filter(ind!="Styrax argenteus var. ramirezii (Greenm.) Gonsoulin") 
dd %<>% separate(ind,into = c("G","S","extra"),extra = "merge",sep=" ") %>% select(-extra) %>% unite("ind",G:S,sep="_")
dd %<>% right_join(.,model_data[[1]],by="ind") %>% mutate(SES = (estimate - mean) / sd,ES=(estimate - mean))
cat("Summary downsampling (only positives)","\n")
dd %>% filter(type=="down") %>% filter(estimate > 0,SES > 0) %>% filter(SES > 1.28) %>% 
summarise(mean(SES),mean(ES),sd(SES),sd(ES)) %>% print
cat("Summary resampling (only positives)","\n")
dd %>% filter(type=="resamp") %>% filter(estimate > 0,SES > 0) %>% filter(SES > 1.96) %>% 
  summarise(mean(SES),mean(ES),sd(SES),sd(ES)) %>% print

dd %>% filter(type=="resamp") %>% select(-trend,-Inter_Cat) %>% setNames(paste0('RSamp_', names(.))) %>% 
  rename(ind=RSamp_ind) %>% 
  left_join(.,dd %>% filter(type=="down") %>% setNames(paste0('DSamp_', names(.))) %>% 
                                                         rename(ind=DSamp_ind),by="ind") %>% 
  mutate(RSamp_SES_cat = cut(RSamp_SES,breaks=c(-6,-1.96,-1.28,0,1.28,1.96,6),labels=c("Strong_Neg","Moderate_Neg","Weak_Neg","Weak_Pos","Moderate_Pos","Strong_Pos"),include.lowest=TRUE)) %>% 
  mutate(DSamp_SES_cat = cut(DSamp_SES,breaks=c(-6,-1.96,-1.28,0,1.28,1.96,6),labels=c("Strong_Neg","Moderate_Neg","Weak_Neg","Weak_Pos","Moderate_Pos","Strong_Pos"),include.lowest=TRUE)) %>% #unite("CAT",RSamp_SES_cat:DSamp_SES_cat) %>%
  mutate(CAT = ifelse(RSamp_estimate <= 0,"Negative","Positive")) -> nulas

  nulas %>% #filter(RSamp_estimate > 0) %>% 
  ggplot(aes(x=RSamp_SES,y=DSamp_SES,fill=CAT,color=CAT)) + 
  geom_vline(xintercept = c(0,-1.28,-1.96,1.28,1.96),size=0.2,linetype="dashed") +
  geom_hline(yintercept = c(0,-1.28,-1.96,1.28,1.96),size=0.2,linetype="dashed") +
  geom_point(alpha=0.6,shape=22,stroke=0.1,size=2) +
  ggside::geom_ysidehistogram(fill="black",stat="identity") +
  ggside::geom_xsidehistogram(fill="black",stat="identity") +
  scale_color_manual(values=MetBrewer::met.brewer("Hiroshige",direction=-1,n=6)) +
  scale_fill_manual(values=MetBrewer::met.brewer("Hiroshige",direction=-1,n=6)) +
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(),
        legend.background = element_blank(),
        legend.key.height = unit(0.01,"cm"),
        legend.key.width = unit(0.01,"cm"),
        axis.title = element_text(family="EB Garamond"),
        axis.text = element_text(family="EB Garamond"),
        #axis.ticks.y=element_blank(),
        #axis.text.y=element_blank(),
        legend.spacing.y = unit(c(1.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond")) +
  labs(x="SES (Resampled null models)",y="SES (Dowsampled null models)") +
  
  NULL



dd %>% filter(type=="resamp") %>% select(-trend,-Inter_Cat) %>% setNames(paste0('RSamp_', names(.))) %>% 
  rename(ind=RSamp_ind) %>% 
  left_join(.,dd %>% filter(type=="down") %>% setNames(paste0('DSamp_', names(.))) %>% 
              rename(ind=DSamp_ind),by="ind") %>% 
  mutate(across(contains("SES"),~cut(.x,breaks=c(-5,-1.96,-1.28,0,1.28,1.96,5),labels=c("Strong_Neg","Moderate_Neg","Weak_Neg","Weak_Pos","Moderate_Pos","Strong_Pos"),include.lowest=TRUE))) %>% 
  filter(DSamp_SES %in% c("Moderate_Neg","Moderate_Pos"),DSamp_estimate > 0)
  group_by(RSamp_SES,DSamp_SES) %>% count() %>% view



dddown <- dd %>% filter(type=="down") %>% filter(estimate > 0)

  panel_B <- ggplot() + 
  geom_segment(aes(y=0,yend=1,x=-0.01,xend=-0.01),color="grey30") +
  stat_dots(data=dddown,aes(x = ES,fill=stat(x <  0),color=stat(x < 0)),size=0.5) + 
  scale_fill_manual(values= colorspace::darken(MetBrewer::met.brewer("Hiroshige",n=9,direction=1)[c(8,1)],amount=0.1),guide=guide_legend(reverse = T,ncol = 1,title.position="top", title.hjust = 0, label.position = "top"),name="") +
  scale_color_manual(values= colorspace::darken(MetBrewer::met.brewer("Hiroshige",n=9,direction=1)[c(8,1)],amount=0.1),guide=guide_legend(reverse = T,ncol = 1,title.position="top", title.hjust = 0, label.position = "top"),name="") +
  geom_text(inherit.aes = T,data = dddown %>% filter(estimate > 0,ES > 0) %>% summarize(n=n()),
            aes(label=paste0(n," species"),x=c(4),y=0.5), 
            position=position_nudge(y=c(0)), size=3,family="EB Garamond") +
  geom_text(inherit.aes = T,data = dddown %>% filter(estimate > 0,ES <= 0) %>% summarise(n=n()),
            aes(label=paste0(n," species"),x = c(-2),y=0.5), 
            position=position_nudge(y=c(0)), size=3,family="EB Garamond") +
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(),
        legend.background = element_blank(),
        legend.key.height = unit(0.01,"cm"),legend.key.width = unit(0.01,"cm"),
        axis.title = element_text(family="EB Garamond"),
        axis.text = element_text(family="EB Garamond"),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.line.y=element_blank(),
        legend.spacing.y = unit(c(1.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond")) +
  labs(x=bquote("Elevation shift in meters"^-year),y="") +
  scale_x_continuous(limits = c(-2.5,4.5)) +
  geom_text(aes(label=c("downslope","upslope"), x = c(-1,1),y = 1.05),color="grey30", size=6,family="EB Garamond",fontface="bold")+
  geom_segment(aes(y=1,yend=1,x=0,xend=-2),color="grey30",arrow = arrow(angle = 20,length=unit(0.2,"cm"))) +
  geom_segment(aes(y=1,yend=1,x=0,xend=2.3),color="grey30",arrow = arrow(angle = 20,length=unit(0.2,"cm"))) +
  NULL

  ddresamp <- dd %<>% filter(type=="resamp") %>% filter(estimate > 0)
  panel_A <- ggplot() + 
    geom_segment(aes(y=0,yend=1,x=-0.01,xend=-0.01),color="grey30") +
    stat_dots(data=dd,aes(x = ES,fill=stat(x <  0),color=stat(x < 0)),size=0.5) + 
    scale_fill_manual(values= colorspace::darken(MetBrewer::met.brewer("Hiroshige",n=9,direction=1)[c(8,1)],amount=0.1),guide=guide_legend(reverse = T,ncol = 1,title.position="top", title.hjust = 0, label.position = "top"),name="") +
    scale_color_manual(values= colorspace::darken(MetBrewer::met.brewer("Hiroshige",n=9,direction=1)[c(8,1)],amount=0.1),guide=guide_legend(reverse = T,ncol = 1,title.position="top", title.hjust = 0, label.position = "top"),name="") +
    geom_text(inherit.aes = T,data = ddresamp %>% filter(estimate > 0,ES > 0) %>% summarize(n=n()),
              aes(label=paste0(n," species"),x=c(20),y=0.5), 
              position=position_nudge(y=c(0)), size=3,family="EB Garamond") +
    geom_text(inherit.aes = T,data = ddresamp %>% filter(estimate > 0,ES <= 0) %>% summarise(n=n()),
              aes(label=paste0(n," species"),x = c(-4),y=0.5), 
              position=position_nudge(y=c(0)), size=3,family="EB Garamond") +
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line(),
          legend.background = element_blank(),
          legend.key.height = unit(0.01,"cm"),legend.key.width = unit(0.01,"cm"),
          axis.title = element_text(family="EB Garamond"),
          axis.text = element_text(family="EB Garamond"),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.line.y=element_blank(),
          legend.spacing.y = unit(c(1.2),"cm"),
          legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond")) +
    labs(x=bquote("Elevation shift in meters"^-year),y="") +
    scale_x_continuous(limits = c(-6,30)) +
    geom_text(aes(label=c("downslope","upslope"), x = c(-3.5,4),y = 1.05),color="grey30", size=6,family="EB Garamond",fontface="bold")+
    geom_segment(aes(y=1,yend=1,x=0,xend=-6),color="grey30",arrow = arrow(angle = 20,length=unit(0.2,"cm"))) +
    geom_segment(aes(y=1,yend=1,x=0,xend=10),color="grey30",arrow = arrow(angle = 20,length=unit(0.2,"cm"))) +
    NULL
  
panel_A


layout <- "
AABB
"

pp<-wrap_plots(list(panel_A,panel_B)) +
  plot_layout(design = layout,byrow=T,nrow = 2,guides="keep",ncol=2,tag_level = 'new') + 
  plot_annotation(tag_levels = 'a',title="") & theme(legend.position = 'none',plot.tag=element_text(family="EB Garamond",size=15,face="bold"),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))
ggsave(filename = here("plots_tables/figure_2.pdf"),plot = pp,height = 7,width = 8.5)

} 
plot_climatic_series <- function(data="interim/ClimateSeries_v.2.RData",alt_limit = 4000){
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

# PLOT LUL SERIES
plot_LUL_chords <- function(data = "output/Historical_data_FULL_LUL.v.1.R",alt_limit = 4000,codes = c("PF","AU","NF","SV"),targetLU="PF") {
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

manken_maps <- function(){ 
data="output/Historical_data_FULL_GLASS.v.1.R"
load(here(data))


load(here("output/LMM_modelData.Rdata"))
load(here("interim/DownSummary_LMM.Rdata"))
down <- e_nulls
load(here("interim/resampSummary_LMM.Rdata"))
resamp <- e_nulls
dd <- down %>% group_by(ind) %>% summarise(mean=mean(estimate),sd=sd(estimate)) %>% mutate(type="down")
dd <- resamp %>% group_by(ind) %>% summarise(mean=mean(estimate),sd=sd(estimate)) %>% mutate(type="resamp") %>% rbind(.,dd)
dd %<>% filter(ind!="Frangula capreifolia var. grandifolia (M.C.Johnst. & L.A.Johnst.) A.Pool") %>% filter(ind!="Styrax argenteus var. ramirezii (Greenm.) Gonsoulin") 
dd %<>% separate(ind,into = c("G","S","extra"),extra = "merge",sep=" ") %>% select(-extra) %>% unite("ind",G:S,sep="_")
dd %<>% right_join(.,model_data[[1]],by="ind") %>% mutate(SES = (estimate - mean) / sd,ES=(estimate - mean))
down %>% distinct(ind) %>% pull -> target_spp
reduct_five %>% filter(Accepted_Name %in% all_of(target_spp)) %>% select(-Altitude,-Accepted_Name,-c(In_Villa:BIOME)) %>% 
  distinct(CellID,.keep_all = TRUE) %>% 
  pivot_longer(-CellID, names_to = "var",values_to = "val") -> mero
mero$val[which(is.infinite(mero$val))] <- NA
mero_mans <- mero %>% mutate(year=stringi::stri_sub(var,-4)) %>% mutate(var=substring(var,1,nchar(var)-5)) %>% drop_na() %>%
  group_by(CellID) %>% nest() %>% rowwise() %>% mutate(manken_cmi = manken(data,tag="_cmi"))
mero_mans %<>% mutate(manken_tas = manken(data,tag="_tas"))
mero_mans %<>% mutate(manken_pr = manken(data,tag="_pr"))
#mero_mans %<>% mutate(manken_pet = manken(data,tag="_pet"))

mero <- reduct_five %>% filter(Accepted_Name %in% all_of(target_spp)) %>% select(CellID,Accepted_Name,Altitude,decimalLongitude,decimalLatitude) %>% distinct(CellID,Accepted_Name,.keep_all = T) %>% left_join(.,mero_mans,by="CellID")

mero %>% distinct(CellID,.keep_all = TRUE) %>%
  mutate(Alt_cat = cut(Altitude,breaks=seq(0,4600,200),labels=seq(0,4400,200))) %>% 
  filter(!is.na(Alt_cat)) %>% 
ggplot(aes(x=(manken_pr),y=(manken_tas),fill=Alt_cat,color=Alt_cat)) +
  geom_point(size=1,stroke=0.05,shape=21) +
  scale_fill_manual(values=MetBrewer::met.brewer("Hiroshige",n=23),name="Altitude") +
  scale_color_manual(values=MetBrewer::met.brewer("Hiroshige",n=23),name="Altitude") +
  stat_ellipse(size=1.5,level=0.9)+
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),axis.line = element_line(),
        legend.position="right",legend.background = element_rect(fill = NA), legend.key.size = unit(0.7, "cm"),legend.key.width = unit(0.6,"cm"), legend.text = element_text(family="EB Garamond",size=8),legend.title = element_text(family="EB Garamond"),legend.spacing = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=14,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5),axis.text = element_text(family="EB Garamond"),axis.title = element_text(family="EB Garamond",size=11)
  )

roads <- rnaturalearth::ne_countries(scale = 110,returnclass = "sf")
xlimits = c(-110,-75)
ylimits = c(5,30)

prec_mk <- ggplot() + 
  geom_sf(data=roads,colour=NA,fill=adjustcolor(MetBrewer::met.brewer("Demuth",direction=-1)[1],alpha.f = 0.2),size=0.5) + xlim(xlimits) + ylim(ylimits) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.border = element_rect(fill=NA),
        legend.position = c(0.2,0.25),
        legend.background = element_rect(fill = NA), legend.key.size = unit(0.7, "cm"),legend.key.width = unit(0.6,"cm"), legend.text = element_text(family="EB Garamond",size=8),legend.title = element_text(family="EB Garamond"),legend.spacing = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5)) +
  geom_point(data = mero %>% filter(abs(manken_pr) >= 1.96),aes(x=decimalLongitude,y=decimalLatitude,fill=manken_pr,color=manken_pr),size=0.5,shape=15) +
  scale_fill_stepsn(colors=(MetBrewer::met.brewer("Isfahan1",direction=1)),name="Mann-Kendall Z",n.breaks=10,labels=function(x) sprintf("%.2f", (x))) +
  scale_color_stepsn(colors=(MetBrewer::met.brewer("Isfahan1",direction=1)),name="Mann-Kendall Z",n.breaks=10,labels=function(x) sprintf("%.2f", (x))) +
  labs(x="",y="") +
  NULL

tas_mk  <- ggplot() + 
  geom_sf(data=roads,colour=NA,fill=adjustcolor(MetBrewer::met.brewer("Demuth",direction=-1)[1],alpha.f = 0.2),size=0.5) + xlim(xlimits) + ylim(ylimits) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.border = element_rect(fill=NA),
        legend.position = c(0.2,0.25),
        legend.background = element_rect(fill = NA), legend.key.size = unit(0.7, "cm"),legend.key.width = unit(0.6,"cm"), legend.text = element_text(family="EB Garamond",size=8),legend.title = element_text(family="EB Garamond"),legend.spacing = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5)) +
  geom_point(data = mero %>% filter(abs(manken_tas) >= 1.96) ,aes(x=decimalLongitude,y=decimalLatitude,fill=manken_tas,color=manken_tas),size=0.5,shape=15) +
  scale_fill_stepsn(colors=(MetBrewer::met.brewer("OKeeffe2",direction=1)),name="Mann-Kendall Z",n.breaks=10,labels=function(x) sprintf("%.2f", (x))) +
  scale_color_stepsn(colors=(MetBrewer::met.brewer("OKeeffe2",direction=1)),name="Mann-Kendall Z",n.breaks=10,labels=function(x) sprintf("%.2f", (x))) +
  labs(x="",y="") +
  NULL

cmi_mk <- ggplot() + 
  geom_sf(data=roads,colour=NA,fill=adjustcolor(MetBrewer::met.brewer("Demuth",direction=-1)[1],alpha.f = 0.2),size=0.5) + xlim(xlimits) + ylim(ylimits) +
  theme(panel.background = element_blank(),panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.border = element_rect(fill=NA),
        legend.position = c(0.2,0.25),
        legend.background = element_rect(fill = NA), legend.key.size = unit(0.7, "cm"),legend.key.width = unit(0.6,"cm"), legend.text = element_text(family="EB Garamond",size=8),legend.title = element_text(family="EB Garamond"),legend.spacing = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5)) +
  geom_point(data = mero %>% filter(abs(manken_cmi) >= 1.96),aes(x=decimalLongitude,y=decimalLatitude,fill=manken_cmi,color=manken_cmi),size=0.5,shape=15) +
  scale_fill_stepsn(colors=(MetBrewer::met.brewer("VanGogh3",direction=1)),name="Mann-Kendall Z",n.breaks=10,labels=function(x) sprintf("%.2f", (x))) +
  scale_color_stepsn(colors=(MetBrewer::met.brewer("VanGogh3",direction=1)),name="Mann-Kendall Z",n.breaks=10,labels=function(x) sprintf("%.2f", (x))) +
  labs(x="",y="") +
  NULL

tas_mk
ggsave(filename = here("plots_tables/figure_3a.pdf"),plot = last_plot(),height = 7,width = 8.5)
prec_mk
ggsave(filename = here("plots_tables/figure_3b.pdf"),plot = last_plot(),height = 7,width = 8.5)
cmi_mk
ggsave(filename = here("plots_tables/figure_3c.pdf"),plot = last_plot(),height = 7,width = 8.5)

ggsave(filename = here("plots_tables//Supplementary_figure_S1.pdf"),plot = pp,height = 7,width = 8.5)
pp





mero %>% mutate(Alt_cat = cut(Altitude,breaks=seq(0,4500,500),labels=seq(0,4000,500))) %>% 
  filter(!is.na(Alt_cat)) %>% 
  filter(abs(manken_tas) >= 1.96) %>% 
  ggplot(aes(x=manken_tas,y=Alt_cat,fill=..x..)) + 
  #ggridges::geom_density_ridges() +
  geom_density_ridges_gradient(size = NA,rel_min_height = 0.001) +
  scale_fill_stepsn(colors=(MetBrewer::met.brewer("OKeeffe2",direction=1)),name="Mann-Kendall Z",n.breaks=10) + 
  theme(panel.background = element_blank(),panel.grid = element_blank(),
       axis.ticks.y  = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.background = element_rect(fill = NA),
        legend.key.size = unit(0.7, "cm"),legend.key.width = unit(0.6,"cm"), legend.text = element_text(family="EB Garamond",size=8),legend.title = element_text(family="EB Garamond"),legend.spacing = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5)) +
  geom_vline(xintercept = c(1.96)) +
  NULL


mero %>% mutate(Alt_cat = cut(Altitude,breaks=seq(0,4500,500),labels=seq(0,4000,500))) %>% 
  filter(!is.na(Alt_cat)) %>% 
  filter(abs(manken_pr) >= 1.96) %>% 
  ggplot(aes(x=manken_pr,y=Alt_cat,fill=..x..)) + 
  #ggridges::geom_density_ridges() +
  geom_density_ridges_gradient(size = NA,rel_min_height = 0.001) +
  scale_fill_stepsn(colors=(MetBrewer::met.brewer("Isfahan1",direction=1)),name="Mann-Kendall Z",n.breaks=10) + 
  theme(panel.background = element_blank(),panel.grid = element_blank(),
       axis.ticks.y  = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.background = element_rect(fill = NA),
        legend.key.size = unit(0.7, "cm"),legend.key.width = unit(0.6,"cm"), legend.text = element_text(family="EB Garamond",size=8),legend.title = element_text(family="EB Garamond"),legend.spacing = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5)) +
  geom_vline(xintercept = c(-1.96,1.96)) +
  NULL


mero %>% mutate(Alt_cat = cut(Altitude,breaks=seq(0,4500,500),labels=seq(0,4000,500))) %>% 
  filter(!is.na(Alt_cat)) %>% 
  filter(abs(manken_cmi) >= 1.96) %>% 
  ggplot(aes(x=manken_cmi,y=Alt_cat,fill=..x..)) + 
  #ggridges::geom_density_ridges() +
  geom_density_ridges_gradient(size = NA,rel_min_height = 0.001) +
  scale_fill_stepsn(colors=(MetBrewer::met.brewer("VanGogh3",direction=1)),name="Mann-Kendall Z",n.breaks=10) + 
  theme(panel.background = element_blank(),panel.grid = element_blank(),
        axis.ticks.y  = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.background = element_rect(fill = NA),
        legend.key.size = unit(0.7, "cm"),legend.key.width = unit(0.6,"cm"), legend.text = element_text(family="EB Garamond",size=8),legend.title = element_text(family="EB Garamond"),legend.spacing = unit(0.1,"cm"),plot.title = element_text(family="EB Garamond",face="bold",size=18,hjust = 0.5),plot.subtitle = element_text(family="EB Garamond",size=12,hjust = 0.5)) +
  geom_vline(xintercept = c(-1.96,1.96)) +
  NULL



}


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
        legend.key.height = unit(0.01,"cm"),
        legend.key.width = unit(0.01,"cm"),
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


do_Granger <- function(data="output/GrangerInput_Mid.Rdata",output="output/GrangerFinal_Mid.Rdata") {
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
  
  
  null_alts <- function(data = "output/Splines_v.2.Rdata", 
                        path_resample = "interim/splines_partial", 
                        path_downsample = "interim",
                        start = 1979, end = 2019, year_bins = 1,
                        var_name = "SumSlopes",output="output/SplinesNULLS_v.2.Rdata") { 
    load(data)
    ## Estimate models for the range: max - min 
    mean_alts <- mean_alts %>% mutate(Type="OBS",.before=data)
    
    list_arch <- list.files(path=here(path_resample),pattern = "_Partial_",full.names = T)
    list_arch_ds <- list.files(path=here(path_downsample),pattern = "_Partial_",full.names = T)
    list_arch_ds <- list_arch_ds[grep("_Alt_",list_arch_ds)]
    
    all_models <- list()
    for (i in 1:length(list_arch_ds)){
      cat("Start processing file",i,"\n")
      cat(" ........ Reading...... ","\n")
      null_models <- mget(load(list_arch_ds[i]))
      test_null <- null_models %>% .[[1]] %>% 
        ungroup() %>% select(Accepted_Name,data) %>% mutate(Accepted_Name=sub("_sim_.*","",Accepted_Name)) %>% mutate(Type="NULL_downsamp",.before=data)
      rm(null_models,nulls_m)
      all_models[[i]] <- test_null %>% select(Accepted_Name,data)
      
    }
    all_models %>% do.call(rbind,.) -> down_models
    down_models %>% distinct(Accepted_Name)
    
    all_models <- list()
    for (i in 1:length(list_arch)){
      cat("Start processing file",i,"\n")
      cat(" ........ Reading...... ","\n")
      null_models <- mget(load(list_arch[i]))
      test_null <- null_models %>% .[[1]] %>% ungroup() %>% select(Accepted_Name,data) %>% mutate(Accepted_Name=sub("_sim_.*","",Accepted_Name)) %>% mutate(Type="NULL_resamp",.before=data)
      rm(null_models,nulls_m)
      all_models[[i]] <- test_null %>% select(Accepted_Name,data)
    }
    all_models %>% do.call(rbind,.) -> resamp_models
    all_models <- list(down_models,resamp_models)
    names(all_models) = c("down_models","resamp_models")
    save(all_models,file=here(output))
    
  }
  
  
  
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
#####

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


do_Granger <- function(data="output/GrangerInput_Mid.Rdata",output="output/GrangerFinal_Mid.Rdata") {
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

#### PLOTS: STAND ALONE PLOTS, WHICH MAKES IT A BIT SLOW
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

plot_Chords <- function(lul = "../2.BMM/CF_Project/output_data/Historical_dataLULC.R",alt_breaks = c(0,300,1000,2000,4000),before = c("LUC_1936","LUC_1956","LUC_1976","LUC_1996"), after = c("LUC_2015","LUC_2015","LUC_2015","LUC_2015"),cat = c("Low","MidLow","MidHigh","High"),codes = c("PF","AU","NF","SV"),targetLU="PF",by_alt=T) {
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



