library(tidyverse)
library(magrittr)
library(here)
library(furrr)
library(glue)
library(lme4)
library(sp)
library(raster,exclude = "select")
library("rnaturalearth")
library("rnaturalearthdata")
library(ggplot2)
library(ggdist)
library(patchwork)
library(showtext)
library(ggridges)
library(ggeffects)
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
ifelse(length(ux) > 0, names(which.max(table(x))), NA)
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
}}
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
start_year = 1979 # starting year to analyse

## Need to check if this one is the latest version
wgetCHELSA_extract <- function(hist_data="interim/Historical_data.Rdata",
                               variable="tas",
                               start=1979,end=2019,
                               file_get="data/envidatS3paths_tas.txt"){ 
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

process_try <- function(data = here("data/TRY/11492.txt"), trait = "Dispersal syndrome",
                        change_cats = here("data/TRY/TRY_CATS_Disp.csv")) { 
try <- data.table::fread(data)
try %>% as_tibble() %>% filter(!DatasetID %in% 
                                 # Exclude restricted datasets
                                 c(429,368,384,398,312,410,75,71,412,214)) %>%
  filter(!is.na(TraitID),TraitName %in% c("Dispersal syndrome","Plant growth form")) %>% 
  count(DatasetID) %>% saveRDS(.,file="output/TRY_datasets.RDS")
  growth_try <- try %>% as_tibble() %>% filter(!DatasetID %in% 
# Exclude restricted datasets
  c(429,368,384,398,312,410,75,71,412,214)) %>%
    filter(!is.na(TraitID),TraitName==trait) %>% select(SpeciesName,OrigValueStr) %>% distinct()

  change <- read.table(change_cats,sep=",",header=T)
  
# This is very crude....
if(trait == "Dispersal syndrome") {
  growth_try %<>% mutate(OrigValueStr = gsub(",","",OrigValueStr)) %>% 
  mutate(OrigValueStr = gsub("\\(","",OrigValueStr)) %>% 
  mutate(OrigValueStr = gsub("\\)","",OrigValueStr)) %>% 
  mutate(OrigValueStr = gsub('"',"",OrigValueStr)) %>% 
  mutate(OrigValueStr = gsub("'","",OrigValueStr)) %>% 
  mutate(OrigValueStr = gsub("  "," ",OrigValueStr))
}

  for (i in 1:nrow(change)){
    cat(i,"\r")
    growth_try$OrigValueStr[which(growth_try$OrigValueStr == change$OrigValueStr[i])] <- change$ChangeValue[i] 
  }
  
  names(growth_try)[1] <- "Accepted_Name"
  growth_try %>% as_tibble() %>% filter(OrigValueStr!="") -> growth_try
  growth_try %>% mutate(Accepted_Name=sub(" ","_",Accepted_Name)) -> growth_try
  save(growth_try,file=here("output/dispersal_try.Rdata"))
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

############## NULL MODELS GENERATION ################
## RUN THE NULL MODELS AND ESTIMATE ROLLING ESTIMATES ## THIS TAKES A WHILE!
# Generate null models through sampling and estimate rolling mean estimates for each null model
alt_resamp_nulls <- function(data="output/Historical_data_FULL_LUL.v.1.R",n_sim=500,start_year=1979,records=100,workers=8) {
  cat("Loading data","\n")
  load(here(data))
  reduct_five %>% select(CellID:Accepted_Name,ID,class:family.x,decimalLatitude:BIOME) %>% 
    mutate(YearCat = as.numeric(sub("Year_","",year)),.after=1) %>% 
    filter(!is.na(Altitude),YearCat >= start_year) %>% 
    distinct(Accepted_Name,CellID,year,.keep_all = T) %>% group_by(Accepted_Name) -> esta
  esta %>% summarise(N=n()) -> tata
  target <- tata %>% filter(N >= records) %>% pull(Accepted_Name) %>% as.character(.)
  ## Estimate elevational ranges for targeted sampling
  esta %>% filter(Accepted_Name %in% all_of(target)) %>% droplevels() %>% group_by(Accepted_Name) %>% 
    summarise(Lower=quantile(Altitude,c(0.025),na.rm=T),Upper=quantile(Altitude,c(0.975),na.rm=T)) -> range_r1
  # Estimate number of samples per year for all species
  esta %>% filter(Accepted_Name %in% all_of(target)) %>% group_by(YearCat,Accepted_Name) %>% summarise(N=n()) %>% arrange(Accepted_Name) -> n_years
  nulls_m <- list()
  count=0
  cat("Setting up multisessions","\n")
  plan(multisession,workers=workers)
  #system.time({ 
  for (i in 1:length(target)){ 
    cat(i,toupper(target[i]),"\n")
    count <- count + 1
    # Identify records (global) within altitudinal range of focal species, and for the specific years it has been sampled.
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
      r_null %>%  mutate(N=years_focal$N,samp = map2(data, N, replace=F,sample_n)) %>% select(-2) %>% 
        unnest(samp) %>% droplevels() %>% group_by(YearCat) %>% 
        summarise(Variable=mean(Altitude,na.rm=T),
                  MaxVariable = quantile(Altitude,probs=0.975,na.rm=T),
                  MinVariable=quantile(Altitude,probs=0.025,na.rm=T),
                  n=n(),class=first(class)) %>% 
        mutate_if(is.double, ~na_if(., -Inf)) %>% 
        mutate_if(is.double, ~na_if(., Inf)) %>% 
        mutate(Accepted_Name=paste(target[i],"_sim_",j,sep="")) %>% 
        arrange(Accepted_Name,YearCat) %>%
        mutate(RollAlt=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=mean,na.rm=TRUE,partial=TRUE,align="center")}),use.names = FALSE)) %>% 
        mutate(RollAlt_Max=unlist(tapply(MaxVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
        mutate(RollAlt_Min=unlist(tapply(MinVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
        
        mutate(RollAlt_sd=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=sd,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
        mutate(RollAlt_Max_sd=unlist(tapply(MaxVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = sd,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
        mutate(RollAlt_Min_sd=unlist(tapply(MinVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = sd,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
        mutate(Rolln=unlist(tapply(n,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=sum,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
        filter(Rolln >= min_by_year) -> r1
    },.progress = TRUE,.options = furrr_options(seed = 123)) ## set seed to avoid warning!
    cat("     Summarizing","simulations","\n")
    simulations %>% ungroup() %>% nest_by(Accepted_Name) -> nulls_m[[i]] 
    names(nulls_m)[[i]] <- as.character(target[i])
    cat("------","\n")
    if(count == 100) { 
      cat("Saving models (partial)","\n")
      nulls_m %>% do.call(bind_rows,.) -> null_models
      nombre_arch <- paste(here("interim/splines_partial_v3","Splines_Resample_Partial_"),i,".Rdata",sep="")
      save(null_models,file=nombre_arch) 
      count=0
      rm(nulls_m,null_models)
      nulls_m <- list()
    }
    if(i==length(target)){
      cat("Saving models (final)","\n")
      nulls_m %>% do.call(bind_rows,.) -> null_models
      nombre_arch <- paste(here("interim/splines_partial_v3","Splines_Resample_Partial_"),i,".Rdata",sep="")
      save(null_models,file=nombre_arch) 
      rm(nulls_m,null_models)
    }
    
  }
  #})
  plan(sequential)
  ## revert ot unparallelized session
}

### HAD TO DO THIS BY BITS DUE TO MEMORY ISSUES (CAN'T TRACE THE PROBLEM), BUT I NEED TO CHANGE THE SEED FOR EACH BIT, OTHERWISE ALL BITS ARE EXACTLY THE SAME.
### THIS A RE THE SEEDS USED: 123, 234, 345, 456, 567 
alt_downsamp_nulls <- function(data="output/Historical_data_FULL_LUL.v.1.R",n_sim=100,start_year=1979,records=100, output = "Splines_Down_Partial_nulls5.Rdata", workers=6, seed=567) {
  cat("Loading data","\n")
  load(here(data))
  reduct_five %>% select(CellID:Accepted_Name,ID,class:family.x,decimalLatitude:BIOME) %>% 
    mutate(YearCat = as.numeric(sub("Year_","",year)),.after=1) %>% 
    filter(!is.na(Altitude),YearCat >= start_year) %>% 
    distinct(Accepted_Name,CellID,year,.keep_all = T) %>% group_by(Accepted_Name) -> esta
  rm(reduct_five)
  
  esta %>% summarise(N=n()) -> tata
  target <- tata %>% filter(N >= records) %>% pull(Accepted_Name) %>% as.character(.)
  ## request parallelization if simulations are set to TRUE
  esta %>% filter(Accepted_Name %in% all_of(target)) %>% count(year,Accepted_Name) %>%
    group_by(year) %>% mutate(n = ifelse(n <= median(n), n, median(n))) -> N
  esta %>% filter(Accepted_Name %in% all_of(target)) %>% group_by(year,Accepted_Name) %>% 
    nest() %>% left_join(.,N,by=c("year","Accepted_Name")) -> esta
  cat("Setting up multisessions","\n")
  plan(multisession,workers=workers)
  #system.time({
  cat("     Performing",n_sim,"simulations","\n")
  ## simulate data by resampling
  simulations <- 1:n_sim %>% future_map_dfr(function(j){
    esta %>% mutate(data=map2(data,n,replace=FALSE,sample_n)) %>% unnest(data) %>% 
      mutate(Accepted_Name=paste(Accepted_Name,"_sim_",j,sep="")) %>% 
      arrange(Accepted_Name,YearCat) -> rarified
    rarified %>% ungroup() %>% group_by(Accepted_Name,YearCat) %>% 
      summarise(Variable=mean(Altitude,na.rm=TRUE),
                MaxVariable = quantile(Altitude,probs=0.975,na.rm=T),
                MinVariable=quantile(Altitude,probs=0.025,na.rm=T),
                n=n(),class=first(class)) %>% 
      mutate_if(is.double, ~na_if(., -Inf)) %>% 
      mutate_if(is.double, ~na_if(., Inf)) %>% 
      arrange(Accepted_Name,YearCat) %>% 
      mutate(RollAlt=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=mean,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
      mutate(RollAlt_Max=unlist(tapply(MaxVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
      mutate(RollAlt_Min=unlist(tapply(MinVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
      mutate(RollAlt_sd=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=sd,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
      mutate(RollAlt_Max_sd=unlist(tapply(MaxVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = sd,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
      mutate(RollAlt_Min_sd=unlist(tapply(MinVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = sd,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
      mutate(Rolln=unlist(tapply(n,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=sum,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
      filter(Rolln >= min_by_year) -> r1
  },.progress = TRUE,.options = furrr_options(seed = seed)) ## set seed to avoid warning!
  cat("     Summarizing","simulations","\n")
  simulations %>% ungroup() %>% nest_by(Accepted_Name) -> nulls_m
  cat("Saving models (final)","\n")
  nombre_arch <- here("interim/splines_partial_v3",output)
  save(nulls_m,file=nombre_arch) 
  rm(nulls_m)
  rm(simulations)
  #})
  plan(sequential)
  ## revert ot unparallelized session
}

## If trend = 0, then there is no shift in the distribution, and the only difference with the resampling would be that occurrences are sampled with weights according to baseline elevations.
## If trend != 0, then the normal distribution is shifted every year.
## If distorion = TRUE, the normal distribution is multiplied by a right-skewing beta distribution.
create_weight <- function(x, y, shape = 3, trend = 2.72 * 0.75, skewed = TRUE) {
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
    kappa = 2 + (0.011*get)
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
                               model_prefix = "Splines_Nulls_ShiftSkewed",
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
  esta %>% filter(Accepted_Name %in% all_of(target)) %>% droplevels() %>% group_by(Accepted_Name) %>%
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
        summarise(Variable=mean(Altitude,na.rm=T),MaxVariable = quantile(Altitude,probs=0.95,na.rm=T),MinVariable=quantile(Altitude,probs=0.05,na.rm=T),n=n(),class=first(class)) %>% 
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

############## ESTIMATION OF ROLLING ESTIMATES ##################
## Estimate rolling mean estimates for observed data and combine null models
# Note: null models already contain rolling estimates
obs_elevation <- function(data="output/Historical_data_FULL_LUL.v.1.R", records=100, output="output/Splines_v.3.Rdata") {
  cat("Loading data","\n")
  load(here(data))
reduct_five %>% select(CellID:Accepted_Name,ID,class:family.x,decimalLatitude:BIOME) %>% 
    mutate(YearCat = as.numeric(sub("Year_","",year)),.after=1) %>% 
    filter(!is.na(Altitude),YearCat >= start_year) %>% distinct(Accepted_Name,CellID,year,.keep_all = T) %>% group_by(Accepted_Name) -> esta
  esta %>% summarise(N=n()) -> tata
  target <- tata %>% filter(N >= records) %>% pull(Accepted_Name) %>% as.character(.)
  cat("Estimating rolling means....","\n")
  esta %>% filter(Accepted_Name %in% all_of(target)) %>% droplevels() %>% group_by(Accepted_Name,YearCat) %>% 
    summarise(Variable=mean(Altitude,na.rm=T), 
              MaxVariable = quantile(Altitude,probs=0.975,na.rm=T),
              MinVariable=quantile(Altitude,probs=0.025,na.rm=T),
              n=n(),class=first(class)) %>% ungroup() %>% 
    mutate_if(is.double, ~na_if(., -Inf)) %>% mutate_if(is.double, ~na_if(., Inf)) %>% arrange(Accepted_Name,YearCat) %>%
    mutate(RollAlt=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean ,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
    mutate(RollAlt_Max=unlist(tapply(MaxVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
    mutate(RollAlt_Min=unlist(tapply(MinVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = mean,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
    mutate(RollAlt_sd=unlist(tapply(Variable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=sd,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
    mutate(RollAlt_Max_sd=unlist(tapply(MaxVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = sd,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
    mutate(RollAlt_Min_sd=unlist(tapply(MinVariable,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window, FUN = sd,na.rm=T,partial=T,align="center")}),use.names = F)) %>%
    mutate(Rolln=unlist(tapply(n,Accepted_Name,function(x) {zoo::rollapply(as_vector(x),width=window,FUN=sum,na.rm=T,partial=T,align="center")}),use.names = F)) %>% 
    filter(Rolln >= min_by_year) -> r1
  summary(r1$Rolln)

cat("Summarising data...... ","\n")
mean_alts <- r1 %>% nest_by(Accepted_Name) %>% mutate(data_points = nrow(data)) %>% filter(data_points >= 10)
save(mean_alts,file=here(output))
cat("Done!!","\n")
}

############## SUMARISE AND ESTIMATE TRENDS #################
# Fit LMMs to observed data
LMM_models <- function(splines="output/Splines_v.3.Rdata",
                       fit_try = TRUE,
                       pattern="max",start_year=1979,
                       end_year=2010,outtable="plots_tables/Supplementary_tableS3.docx"){
  
  load(here(splines))
  pattern = NULL
  aver <- mean_alts %>% mutate(Type="OBS",.before=data) %>% select(2:3) %>% unnest(data)
  if(is.null(pattern)) {
    pat = "RollAlt"
    pat_sd = paste0(pat,"_sd")
    names <- c("ind","intercept_mid","trend_mid")
  }
  if(!is.null(pattern)) { 
    pat = paste0("RollAlt_",pattern)
    pat_sd = paste0(pat,"_sd")
    names <- c("ind",paste0(c("intercept_","trend_"),pattern))
  }
  aver %<>% ungroup() %>% 
    select(Accepted_Name,YearCat,ends_with(pat),Rolln,starts_with(pat_sd)) %>% 
    mutate(YearCat=as.numeric(YearCat)) %>% 
    mutate(RollAlt_error = unlist(across(starts_with(pat_sd), ~ .x / sqrt(Rolln)))) %>% 
    select(ends_with(pat),YearCat,Accepted_Name,RollAlt_error) %>% rename_with(~c("y","tme","ind","error"))
  aver %<>% filter(tme <= end_year) %>% mutate(year=tme) %>% mutate(tme = tme - start_year)
  otra <- aver %>% separate(ind,into = c("G","S","extra"),extra = "merge",sep=" ") %>% select(-extra) %>% unite("ind",G:S,sep="_")
  otra <- otra %>% filter(year<=2010)
to_merge <- otra %>% group_by(ind) %>% summarise(error = median(error,na.rm=TRUE))

if(fit_try){ 
target_try <- readRDS("output/TRY_processed.RDS")
otra %<>% left_join(., target_try %>% rename("ind"=binomial) %>% select(ind,Disp_Syndrome,Disp_Syndrome_Inf,LifeForm),by="ind")

# Do a bit of wrangling of the try data to get inferred dispersal syndrome by genus
target_try %<>% select(binomial,Disp_Syndrome_Inf) %>% 
  rename("Disp_Syndrome_Genus" = Disp_Syndrome_Inf) %>% 
  separate(binomial,c("Genus","spp"),sep="_",remove = FALSE) %>% 
  distinct(Genus,.keep_all = TRUE) %>% select(-binomial,-spp)
otra %<>% separate(ind,c("Genus","spp"),sep="_",remove=FALSE) %>% 
  left_join(., target_try,by="Genus" ) %>% select(-Disp_Syndrome_Inf)

otra %<>% mutate(Disp_Syndrome = replace_na(Disp_Syndrome,"Unknown")) %>% 
  mutate(LifeForm = replace_na(LifeForm,"Unknown")) %>% 
  mutate(Disp_Syndrome = ifelse(Disp_Syndrome == "Water","Unknown",Disp_Syndrome))
otra %<>% mutate(LifeForm= factor(LifeForm,levels=c("Unknown","Herb","Shrub","Epiphyte/Climber","Tree"))) %>%
  mutate(Disp_Syndrome= factor(Disp_Syndrome,levels=c("Unknown","Unassisted","Animal","Wind")))

fit1 <- lmer(y ~ 1 + tme + tme:LifeForm + tme:error + error + (1|ind) + (0 + tme|ind),data = otra, control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4)))

broom.mixed::tidy(fit1) %>% mutate(across(4:6, \(x) round(x,3))) %>% 
  flextable::flextable() %>% 
  flextable::colformat_num(x = .,decimal.mark=".") %>% flextable::save_as_docx(path = here("plots_tables/Supplementary_Table_S2a_max.docx"))

ggpredict(fit1,terms=c("tme","LifeForm"),type="random",interval="confidence") %>% 
  hypothesis_test(test= "pairwise") %>% as_tibble() %>% 
  mutate(across(3:6, \(x) round(x,3))) %>% 
  flextable::flextable() %>% 
  flextable::colformat_num(x = .,decimal.mark=".") %>% flextable::save_as_docx(path = here("plots_tables/Supplementary_Table_S3a_max.docx"))

a <- ggpredict(fit1,terms=c("tme","LifeForm"),type="random",interval="confidence") %>%
  ggplot(aes(x=x+1979,y=predicted,group=group,color=group)) + 
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill=group),color="NA", alpha = 0.1) +
  scale_color_met_d("Hiroshige") +  
  scale_fill_met_d("Hiroshige") +
  theme(legend.position = "right",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line()) +
  labs(x="", y = "Mid-point elevation (masl)")+
  guides(color = guide_legend(title="Growth Form",reverse = TRUE),
         fill = guide_legend(title="Growth Form",reverse = TRUE)) +
  NULL

fit1 <- lmer(y ~ 1 + tme + tme:Disp_Syndrome + tme:error + error + (1|ind) + (0 + tme|ind),data = otra, control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4)))
broom.mixed::tidy(fit1) %>% mutate(across(4:6, \(x) round(x,3))) %>% 
  flextable::flextable() %>% 
  flextable::colformat_num(x = .,decimal.mark=".") %>% flextable::save_as_docx(path = here("plots_tables/Supplementary_Table_S2b_max.docx"))

ggpredict(fit1,terms=c("tme","Disp_Syndrome"),type="random",interval="confidence") %>% 
  hypothesis_test(test= "pairwise") %>% as_tibble() %>% 
  mutate(across(3:6, \(x) round(x,3))) %>% 
  flextable::flextable() %>% 
  flextable::colformat_num(x = .,decimal.mark=".") %>% flextable::save_as_docx(path = here("plots_tables/Supplementary_Table_S3b_max.docx"))

b <- ggpredict(fit1,terms=c("tme","Disp_Syndrome"),type="random",interval="confidence") %>%
  ggplot(aes(x=x+1979,y=predicted,group=group,color=group)) + 
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill=group),color="NA", alpha = 0.1) +
  scale_color_met_d("Hiroshige") +  
  scale_fill_met_d("Hiroshige") +
  theme(legend.position = "right",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line()) +
  labs(x="", y = "Mid-point elevation (masl)")+
  guides(color = guide_legend(title="Dispersal",reverse = TRUE),
         fill = guide_legend(title="Dispersal",reverse = TRUE)) +
  NULL

layout <- "AAAC\nBBBC"
a + b + guide_area() +
  plot_layout(design = layout, byrow=T, nrow = 2, guides = "collect", ncol=2,
              tag_level = 'new') +
  plot_annotation(tag_levels = 'a') &
  theme(
    plot.title = element_text(family="EB Garamond",size=16,face="bold"),
    plot.tag = element_text(family="EB Garamond",size=15,face="bold"),
    plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
    legend.key.height=unit(0.7,"cm"),
    legend.text=element_text(size=8,family="EB Garamond"), 
    legend.title=element_text(size=10,family="EB Garamond"),
    axis.title = element_text(family="EB Garamond",size=10))

ggsave(filename = "plots_tables/Supplementary_figure_Max_Traits.pdf",plot = last_plot())


}

# Assess estimation error through time
  fit_error <- lmer(error ~ 1 + tme + (1|ind) + (0 + tme|ind),data = otra, control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4)))
  err <- broom.mixed::tidy(fit_error)
  err
  
  # Assess elevation through time
  fit_pre <- lmer(y ~ 1 + tme + (1|ind) + (0 + tme|ind),data = otra, control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4)))
  pre <- broom.mixed::tidy(fit_pre)
  pre
  
  ### Add error as covariate
  fit0 <- lmer(y ~ 1 + tme + error + (1|ind) + (0 + tme|ind),data = otra, control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4)))
  cov <- broom.mixed::tidy(fit0)
  cov
  
  ### Add error as covariate + interaction
  fit1 <- lmer(y ~ 1 + tme + tme:error + error + (1|ind) + (0 + tme|ind),data = otra, control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4)))
  covinter <- broom.mixed::tidy(fit1)
  covinter
  
  ggpredict(fit1,terms=c("tme"),type="random",interval="confidence") %>% hypothesis_test()
  
  anova(fit_pre,fit0,fit1)
  
  modelos <- list("error" = fit_error, "time"=fit_pre,"cov"=fit0,"interaction"=fit1)
  save(modelos,file="output/LMM_models_Mid_v3.RData")
  
  
  if(is.null(pattern)){  
    esta <- coef(fit0)$ind %>% as_tibble(rownames = "ind") %>% #mutate(Inter_Cat=cut(`(Intercept)`,breaks=seq(0,4000,200),include.lowest=T,labels=seq(0,4000,200)[-1] )) %>% 
      select(ind,`(Intercept)`,tme) %>% rename_with(~names) %>% left_join(.,to_merge,by="ind")}
  
  if(!is.null(pattern)){  
    esta <- coef(fit0)$ind %>% as_tibble(rownames = "ind") %>% #mutate(Inter_Cat=cut(`(Intercept)`,breaks=seq(0,4000,200),include.lowest=T,labels=seq(0,4000,200)[-1] )) %>% 
      select(ind,`(Intercept)`,tme) %>% rename_with(~names) %>% left_join(.,esta,by="ind")}
  
  saveRDS(esta,file = "output/LMM_trends_v.3.RDS")
}

## Combine simulations for each null model
null_models <- function(data = "output/Splines_v.3.Rdata", 
                      path = "interim/splines_partial_v3",
                      start = 1979, end = 2019, year_bins = 1, 
                      output_prefix="SplinesNULLS") { 

# Resampling  
  list_arch <- list.files(path=here(path),pattern = "Resample_",full.names = TRUE)
  all_models <- list()
for (i in 1:length(list_arch)){
    cat("Start processing file",i,"\n")
    cat(" ........ Reading...... ","\n")
    null_models <- mget(load(list_arch[i]))
    test_null <- null_models %>% .[[1]] %>% ungroup() %>% select(Accepted_Name,data) %>% mutate(Accepted_Name=sub("_sim_.*","",Accepted_Name)) %>% mutate(Type="NULL_resamp",.before=data)
    rm(null_models)
    all_models[[i]] <- test_null %>% select(Accepted_Name,data)
  }
  all_models %>% do.call(rbind,.) -> resamp_models
saveRDS(resamp_models,file=here("output",paste0(output_prefix,"_resamp.RDS")))
  rm(resamp_models)
  
# Downsample
  list_arch_ds <- list.files(path=here(path),pattern = "Down_Partial_",full.names = T)
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
  saveRDS(down_models,file=here("output",paste0(output_prefix,"_downsample.RDS")))
  rm(down_models)
  
# Shifting sample
#  list_arch_sh <- list.files(path=here(path),pattern = "Nulls_Shift_",full.names = T)
#  all_models <- list()
#  for (i in 1:length(list_arch_sh)){
 #   cat("Start processing file",i,"\n")
#    cat(" ........ Reading...... ","\n")
#    null_models <- mget(load(list_arch_sh[i]))
#    test_null <- null_models %>% .[[1]] %>% ungroup() %>% select(Accepted_Name,data) %>% mutate(Accepted_Name=sub("_sim_.*","",Accepted_Name)) %>% mutate(Type="NULL_Shifting",.before=data)
#    rm(null_models)
#    all_models[[i]] <- test_null %>% select(Accepted_Name,data)
#  }
#  all_models %>% do.call(rbind,.) -> shift_models
  # saveRDS(shift_models,file=here("output",paste0(output_prefix,"_shifting.RDS")))
  # rm(shift_models)
  # 
  # # Skewed sample 
  # list_arch_di <- list.files(path=here(path),pattern = "Nulls_DistortedNoShift_",full.names = T)
  # all_models <- list()
  # for (i in 1:length(list_arch_di)){
  #   cat("Start processing file",i,"\n")
  #   cat(" ........ Reading...... ","\n")
  #   null_models <- mget(load(list_arch_di[i]))
  #   test_null <- null_models %>% .[[1]] %>% ungroup() %>% select(Accepted_Name,data) %>% mutate(Accepted_Name=sub("_sim_.*","",Accepted_Name)) %>% mutate(Type="NULL_Skewed",.before=data)
  #   rm(null_models)
  #   all_models[[i]] <- test_null %>% select(Accepted_Name,data)
  # }
  # all_models %>% do.call(rbind,.) -> skewed_models
  # saveRDS(skewed_models,file=here("output",paste0(output_prefix,"_skewed.RDS")))
  # rm(skewed_models)
  
  
  # # Skewed-shifting sample
  # list_arch_dish <- list.files(path=here(path),pattern = "Nulls_ShiftSkewed",full.names = T)
  #  all_models <- list()
  # for (i in 1:length(list_arch_dish)){
  #   cat("Start processing file",i,"\n")
  #   cat(" ........ Reading...... ","\n")
  #   null_models <- mget(load(list_arch_dish[i]))
  #   test_null <- null_models %>% .[[1]] %>% ungroup() %>% select(Accepted_Name,data) %>% mutate(Accepted_Name=sub("_sim_.*","",Accepted_Name)) %>% mutate(Type="NULL_SkewedShifting",.before=data)
  #   rm(null_models)
  #   all_models[[i]] <- test_null %>% select(Accepted_Name,data)
  # }
  # all_models %>% do.call(rbind,.) -> skewshift_models
  # saveRDS(skewshift_models,file=here("output",paste0(output_prefix,"_shiftingSkew.RDS")))
  # rm(skewshift_models)
  # 
}

## Fit LMMs to every null model
summary_LMMnulls_elevation <- function(){
  nulls="output/SplinesNULLS_downsample.RDS"
  output_prefix="SummaryLMM_downsamp"
  var="mid"
  n_sim=500
  end_year=2010
  start_year=1979
  models <- readRDS(here(nulls)) ### null models
  models <- models %>% group_by(Accepted_Name) %>% 
    mutate(Type=paste("SIMS",1:n_sim,sep="_"),.before=data) %>% 
    group_by(Accepted_Name,Type) %>% unnest(data)
  pat_sd = paste0("RollAlt","_sd")
  
  models <- models %>% ungroup() %>%
    select(starts_with("RollAlt"),YearCat,Accepted_Name,Type,starts_with(pat_sd),Rolln) %>% 
    mutate(YearCat=as.numeric(YearCat)) %>% 
    rename_with(~c("y","y.max","y.min","error","error.max","error.min","tme","ind","sim","Rolln")) %>%
    mutate(across(starts_with("error"), ~ .x / sqrt(Rolln))) %>% 
    filter(tme <= end_year) %>% mutate(year=tme) %>% mutate(tme = tme - start_year)

  models %<>% group_by(sim) %>% group_split()
  e_nulls <- list()
  sum_nulls <- list()
  trends <- c()
  for(i in 1:length(models)){
    poruna <- models[[i]]
    cat("Simulation -- ",i,"\n")
    cat("Fitting model","\r")
    if(var == "min"){ fit1 <- lmer(y.min ~ 1 + tme + error.min + (1|ind) + (0 + tme|ind),data = poruna, control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4)))}
    if(var == "max"){ fit1 <- lmer(y.max ~ 1 + tme + error.max + (1|ind) + (0 + tme|ind),data = poruna, control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4)))}
    if(var == "mid"){ fit1 <- lmer(y ~ 1 + tme + error + (1|ind) + (0 + tme|ind),data = poruna, control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4)))}
    #sum_nulls[[i]] <- fit1
    trends[[i]] <- broom.mixed::tidy(fit1)[2,4]
    ef <-  coef(fit1)$ind %>% as_tibble(rownames = "ind") %>% select(tme,ind) %>% pivot_longer(cols = 1,names_to = "trend",values_to = "estimate" ) %>% mutate(trend = sub("tme","1979-2019",trend)) %>% mutate(sim=paste0("SIMS_",i))
    e_nulls[[i]] <- ef
  }
  do.call(rbind,e_nulls) %>% as_tibble() -> e_nulls
  data.table::rbindlist(trends) %>% summary()
  mean(e_nulls$estimate)
  e_nulls <- list(e_nulls,#sum_nulls,
                  data.table::rbindlist(trends))
  save(e_nulls,file= here("output",paste0(output_prefix,"_",var,"_v3.RData")) ) 

  
}

## Plot observed trends

merge_TRY <- function(TRY="output/growth_try.Rdata",  TRY2="output/dispersal_try.Rdata") {
  ## Read the TRY data and merge
###### There was a mysterious bug here with the old data file, re-working with new file (13 march 2024) ########
load(here(TRY))
  growth_try %<>% rename("binomial"=Accepted_Name) %>% rename("growth"=OrigValueStr) %>% mutate(growth=as.factor(growth))
  levels(growth_try$growth) <- c("Aquatic","Epiphyte/Climber","Epiphyte/Climber","Herb","Shrub","Shrub","Tree","Tree")
  growth_try %>% group_by(binomial) %>% nest() %>% rowwise() %>% mutate(LifeForm=Mode(data$growth)) -> target_try
  
load(here(TRY2))
growth_try %<>% rename("binomial"=Accepted_Name) %>% rename("dispersal"=OrigValueStr) %>% mutate(dispersal=as.factor(dispersal))
growth_try %<>% separate(binomial,into=c("Genus","spp"),sep="_",remove=FALSE) %>% group_by(Genus) %>% nest() %>% rowwise() %>% 
  mutate(dispersal_genus = Mode(data$dispersal)) %>% unnest(data) %>% select(-spp)

growth_try %>% group_by(binomial) %>% nest() %>% rowwise() %>% mutate(Disp_Syndrome=Mode(data$dispersal),Disp_Syndrome_Inf = Mode(data$dispersal_genus)) %>% full_join(.,target_try,by="binomial") -> target_try
target_try %<>% select(-data.x,-data.y) 

saveRDS(target_try,here("output","TRY_processed.RDS"))
}

plot_observed <- function(){
model_title = "Observed_trends"
mids = "output/LMM_models_Mid_v3.RData"
mins = "output/LMM_models_Min_v3.RData"
maxs = "output/LMM_models_Max_v3.RData"
target_try <- readRDS("output/TRY_processed.RDS")
esta <- readRDS(here("output/LMM_trends_v.3.RDS"))
e <- esta %>%
  mutate(Inter_Cat=cut(intercept_mid,breaks=seq(0,4600,100),include.lowest=T,labels=seq(0,4600,100)[-1] )) %>% 
  pivot_longer(cols = -c(ind,Inter_Cat),names_to = "trend",values_to = "estimate" ) %>% 
  left_join(., target_try %>% rename("ind"=binomial) %>% select(ind,Disp_Syndrome,Disp_Syndrome_Inf,LifeForm),by="ind")

# Do a bit of wrangling of the try data to get inferred dispersal syndrome by genus
target_try %<>% select(binomial,Disp_Syndrome_Inf) %>% 
  rename("Disp_Syndrome_Genus" = Disp_Syndrome_Inf) %>% 
  separate(binomial,c("Genus","spp"),sep="_",remove = FALSE) %>% 
  distinct(Genus,.keep_all = TRUE) %>% select(-binomial,-spp)
e %<>% separate(ind,c("Genus","spp"),sep="_",remove=FALSE) %>% 
  left_join(., target_try,by="Genus" ) %>% select(-Disp_Syndrome_Inf)


  as <- e %>% filter(trend=="error") %>% reframe(summary(estimate)) %>% pull(1)
  as <- round(as,1)
  as <- paste(as,collapse = ",")
  as <- paste0("error ","[",as,"]")

  ## Read linear mixed models, predict and then plot
load(here(mids))
  mid <- ggpredict(modelos$interaction,terms=c("tme",as),type="random",interval="confidence") %>%
    ggplot(aes(x=x+1979,y=predicted,group=group,color=group)) + 
    geom_line() +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill=group),color="NA", alpha = 0.1) +
    scale_color_met_d("Hiroshige") +  
    scale_fill_met_d("Hiroshige") +
    theme(legend.position = "right",
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line()) +
    labs(x="", y = "Mid-point elevation (masl)")+
    guides(color = guide_legend(title="Error",reverse = TRUE),
           fill = guide_legend(title="Error",reverse = TRUE)) +
    coord_cartesian(ylim=c(1000,1500)) +
    NULL
  
  load(here(mins))
  min <- ggpredict(modelos$interaction,terms=c("tme",as),type="random",interval="confidence") %>%
    ggplot(aes(x=x+1979,y=predicted,group=group,color=group))  + geom_line() +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill=group),color="NA", alpha = 0.1) +
    scale_color_met_d("Hiroshige") +  
    scale_fill_met_d("Hiroshige") +
    theme(legend.position = "right",
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line()) +
    labs(x="",y="Lower-point elevation (masl)")+
    guides(color = guide_legend(title="Error",reverse = TRUE),
           fill = guide_legend(title="Error",reverse = TRUE)) +
    coord_cartesian(ylim=c(400,1000)) +
    NULL
  
  load(here(maxs))
  max <- ggpredict(modelos$interaction,terms=c("tme",as),type="random",interval="confidence") %>%
    ggplot(aes(x=x+1979,y=predicted,group=group,color=group))  + geom_line() +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill=group),color="NA", alpha = 0.1) +
    scale_color_met_d("Hiroshige") +  
    scale_fill_met_d("Hiroshige") +
    theme(legend.position = "right",
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line()) +
    labs(x="",y="Upper-point elevation (masl)")+
    guides(color = guide_legend(title="Error",reverse = TRUE),
           fill = guide_legend(title="Error",reverse = TRUE)) +
    coord_cartesian(ylim=c(1500,2000)) +
    NULL
  
layout <- "AAAD\nBBCC"

  mid + min + max + guide_area() +
    plot_layout(design = layout, byrow=T, nrow = 2, guides = "collect", ncol=2,
                tag_level = 'new') +
    plot_annotation(tag_levels = 'a',title=model_title) &
    theme(
      plot.title = element_text(family="EB Garamond",size=16,face="bold"),
plot.tag = element_text(family="EB Garamond",size=15,face="bold"),
plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
legend.key.height=unit(0.7,"cm"),
legend.text=element_text(size=8,family="EB Garamond"), 
legend.title=element_text(size=10,family="EB Garamond"),
axis.title = element_text(family="EB Garamond",size=10))
  
ggsave(filename = "plots_tables/Supplementary_figure_SX.pdf",plot = last_plot())


}

## Plot null trends
plot_SES <- function(){
target_try <- readRDS("output/TRY_processed.RDS")
target_try_2 <- target_try %>% select(binomial,Disp_Syndrome_Inf) %>% 
    rename("Disp_Syndrome_Genus" = Disp_Syndrome_Inf) %>% 
    separate(binomial,c("Genus","spp"),sep="_",remove = FALSE) %>% 
    distinct(Genus,.keep_all = TRUE) %>% select(-binomial,-spp)
  
model_labs = c("Resampling","Downsampling")
  model_titles = paste0(model_labs," null models")
  listMods <- list.files(here("output"),pattern = "SummaryLMM",full.names = TRUE)
  cols <- met.brewer("Kandinsky",n=3)
  model <- c("resamp_","downsamp_")
  plots <- list()
  data <- list()
  obs <- readRDS(here("output/LMM_trends_v.3.RDS"))
  obs %<>% select(ind,starts_with("trend")) %>% rename_with(.cols = starts_with("trend"),
                                                           ~sub("trend","obs",.x))
  lines <- tibble(y = rep(1:3,each=3), yend = rep(1.8:3.8,each=3), 
                  x = rep(c(-1.28,0,1.28),times = 3), 
                  xend = rep(c(-1.28,0,1.28),times = 3))
  
for(i in 1:length(model)){
  cat("Reading and plotting null model",model[i],"\n")
  s <- listMods[grep(model[i],listMods)]
  mids = s[grep("_mid_",s)]
  mins = s[grep("_min_",s)]
  maxs = s[grep("_max_",s)]
load(here(mids))
esta <- e_nulls[[1]] %>% group_by(ind) %>% 
  summarise(trend_mid=mean(estimate),sd_mid=sd(estimate))
load(here(mins))
esta <- e_nulls[[1]] %>% group_by(ind) %>% 
    summarise(trend_min=mean(estimate),sd_min=sd(estimate)) %>% 
    left_join(.,esta,by="ind")
load(here(maxs))
esta <- e_nulls[[1]] %>% group_by(ind) %>% 
  summarise(trend_max=mean(estimate),sd_max=sd(estimate)) %>% 
  left_join(.,esta,by="ind")

# Estimate SES and ES
e <- esta %>% separate(ind,into = c("G","S"),sep = " ") %>% 
  unite(col = "ind",c(G,S), sep = "_") %>% 
  left_join(.,obs,by="ind") %>% 
  mutate(SES_mid = (obs_mid - trend_mid) / sd_mid, ES_mid=(obs_mid - trend_mid)) %>% 
  mutate(SES_min = (obs_min - trend_min) / sd_min, ES_min=(obs_min - trend_min)) %>% 
  mutate(SES_max = (obs_max - trend_max) / sd_max, ES_max=(obs_max - trend_max)) #%>% 

e %<>% left_join(., target_try %>% rename("ind"=binomial) %>% select(ind,Disp_Syndrome,Disp_Syndrome_Inf,LifeForm),by="ind") %>% 
  separate(ind,c("Genus","spp"),sep="_",remove=FALSE) %>% 
   left_join(., target_try_2,by="Genus" ) %>% select(-Disp_Syndrome_Inf,-Genus,-spp)

e %<>% pivot_longer(starts_with("SES"),names_to = "SES",values_to = "estimate")
name <- here("output",paste0("SES_",model,"v.3.RDS"))[i]
saveRDS(e,file = name)
e %<>% mutate(SES = factor(SES,levels=c("SES_min","SES_mid","SES_max")))
data[[i]] <- e

numbers <- e %>% mutate(Samp_SES = cut(estimate,breaks=c(min(estimate,na.rm=TRUE),-1.96,-1.28,0,1.28,1.96,max(estimate,na.rm=TRUE)),labels=c("Strong_Neg","Moderate_Neg","Weak_Neg","Weak_Pos","Moderate_Pos","Strong_Pos"),include.lowest=TRUE)) %>% 
  group_by(SES,Samp_SES) %>% count() %>% group_by(SES) %>% mutate(prop = n / sum(n)) %>%  filter(!grepl("Weak",Samp_SES)) %>% summarise(n = sum(n),prop = sum(prop))

plots[[i]] <-  ggplot()  + 
    ggdist::geom_dots(data = e , aes(x = abs(estimate), y= SES,
                                     fill=SES,
                                     color = SES,
                                     #side = ifelse(estimate < 0, "bottom","top")
                                     ),
                                     inherit.aes = FALSE,shape=21) +
    scale_color_met_d("Kandinsky") +  
    scale_fill_met_d("Kandinsky") +
    geom_hline(yintercept = 1:3,linewidth=0.1) +
    geom_vline(xintercept = c(1.28),color="grey30", linetype = c("dashed"),linewidth=0.4) +
    geom_segment(aes(x=1.28,xend=10,y=3.8,yend=3.8),arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "open")) +
    geom_text(data = numbers %>% arrange(desc(SES)),aes(x=rep(10,3),y=seq(1.5,3.5,1),label=paste0(n," species")),family="EB Garamond") +
    #geom_vline(xintercept = 0) +
    theme(legend.position = "none",
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line.x = element_line(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(family="EB Garamond",size=14,hjust = 0.1,face="bold.italic"),
          axis.text.x = element_text(family="EB Garamond",size=12),
          axis.title.x = element_text(family="EB Garamond",size=12,face="bold"),
          axis.text.y = element_text(family="EB Garamond",size=12,face="bold")) +
  scale_y_discrete(expand = c(0, 0.1),labels=c("Lower-point", "Mid-point", "Upper-point")) +
  scale_x_continuous(trans="sqrt",breaks = c(0,1.28,10,20)) +
    labs(y="",x="Standardised Effect Sizes (absolute value)") +
    guides(color = guide_legend(title="Elevation point",reverse = TRUE),
           fill = guide_legend(title="Elevation point",reverse = TRUE)) +
    NULL
  
}
    
wrap_plots(plots) +
    plot_layout(byrow=T, nrow = 2, guides = "collect", ncol=1,tag_level = 'new',axes = "collect") +
    plot_annotation(tag_levels = 'a') & theme(
      plot.title = element_text(family="EB Garamond",size=10,face="bold"),
      plot.tag = element_text(family="EB Garamond",size=16,face="bold"),
      plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
      legend.key.height=unit(0.7,"cm"),
      legend.text=element_text(size=8,family="EB Garamond"), 
      legend.title=element_text(size=10,family="EB Garamond"),
      axis.title = element_text(family="EB Garamond",size=10))

ggsave(filename = "plots_tables/Supplementary_figure_SES.pdf",plot = last_plot())
  

re <- data[[1]] %>% select(ind,SES,estimate,Disp_Syndrome_Genus,LifeForm,obs_max ,obs_min ,obs_mid) %>% 
  rename("RS_SES" = estimate) %>%
  filter(!grepl("Frangula_capreifolia|Styrax_argenteus",ind)) %>%
  left_join(.,data[[2]] %>% select(ind,SES,estimate) %>% rename("DS_SES" = estimate) %>%  
              filter(!grepl("Frangula_capreifolia|Styrax_argenteus",ind)),by=c("ind","SES"))

re <- data[[1]] %>% select(ind,SES,estimate,Disp_Syndrome_Genus,LifeForm,obs_max ,obs_min ,obs_mid) %>% 
  rename("RS_SES" = estimate) %>% bind_cols(.,data[[2]] %>% select(estimate) %>% rename("DS_SES" = estimate))


mid_p <- re %>% filter(SES=="SES_mid") %>% mutate(Trend = ifelse(obs_mid < 0,"Downslope","Upslope")) %>% 
  ggplot(aes(x=RS_SES,y=DS_SES,col=Trend,fill=Trend)) + 
  geom_vline(xintercept = c(0,-1.28,1.28),linewidth=0.2,linetype="dashed") +
  geom_hline(yintercept = c(0,-1.28,1.28),linewidth=0.2,linetype="dashed") +
  geom_point(alpha=0.6,shape=22,stroke=0.1,size=2) +
  ggside::geom_ysidehistogram(stat="identity") +
  ggside::geom_xsidehistogram(stat="identity") +
  scale_color_manual(values=MetBrewer::met.brewer("Hiroshige",direction=-1,n=6)) +
  scale_fill_manual(values=MetBrewer::met.brewer("Hiroshige",direction=-1,n=6)) +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        legend.background = element_blank(),
        legend.key.height = unit(0.01,"cm"),
        legend.key.width = unit(0.01,"cm"),
        axis.title = element_text(family="EB Garamond"),
        axis.text = element_text(family="EB Garamond"),
        legend.spacing.y = unit(c(1.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond")) +
  labs(x="SES (Resampled null models)",y="SES (Dowsampled null models)") +
  NULL

min_p <- re %>% filter(SES=="SES_min") %>% mutate(Trend = ifelse(obs_min < 0,"Downslope","Upslope")) %>% 
ggplot(aes(x=RS_SES,y=DS_SES,col=Trend,fill=Trend)) + 
  geom_vline(xintercept = c(0,-1.28,1.28),linewidth=0.2,linetype="dashed") +
  geom_hline(yintercept = c(0,-1.28,1.28),linewidth=0.2,linetype="dashed") +
  geom_point(alpha=0.6,shape=22,stroke=0.1,size=2) +
  ggside::geom_ysidehistogram(stat="identity") +
  ggside::geom_xsidehistogram(stat="identity") +
  scale_color_manual(values=MetBrewer::met.brewer("Hiroshige",direction=-1,n=6)) +
  scale_fill_manual(values=MetBrewer::met.brewer("Hiroshige",direction=-1,n=6)) +
  theme(panel.background = element_blank(),
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

max_p <- re %>% filter(SES=="SES_max") %>% mutate(Trend = ifelse(obs_max < 0,"Downslope","Upslope")) %>% 
ggplot(aes(x=RS_SES,y=DS_SES,col=Trend,fill=Trend)) + 
  geom_vline(xintercept = c(0,-1.28,1.28),linewidth=0.2,linetype="dashed") +
  geom_hline(yintercept = c(0,-1.28,1.28),linewidth=0.2,linetype="dashed") +
  geom_point(alpha=0.6,shape=22,stroke=0.1,size=2) +
  scale_color_manual(values=MetBrewer::met.brewer("Hiroshige",direction=-1,n=6)) +
  scale_fill_manual(values=MetBrewer::met.brewer("Hiroshige",direction=-1,n=6)) +
  ggside::geom_ysidehistogram(stat="identity") +
  ggside::geom_xsidehistogram(stat="identity") +
  theme(panel.background = element_blank(),
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

layout <- "AAAD\nBBCC"
mid_p + min_p + max_p + guide_area() +
  plot_layout(design = layout, byrow=T, nrow = 2, guides = "collect", ncol=2,
              tag_level = 'new',axes = "collect") +
  plot_annotation(tag_levels = 'a') &
  theme(
    plot.title = element_text(family="EB Garamond",size=10,face="bold"),
    plot.tag = element_text(family="EB Garamond",size=16,face="bold"),
    plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
    legend.key.height=unit(0.7,"cm"),
    legend.text=element_text(size=8,family="EB Garamond"), 
    legend.title=element_text(size=10,family="EB Garamond"),
    axis.title = element_text(family="EB Garamond",size=10))
ggsave(filename = "plots_tables/Supplementary_figure_SES_2.pdf",plot = last_plot())

re %>% filter(SES == "SES_mid") %>% count(Random = abs(RS_SES) >= 1.28, Robust = abs(DS_SES) <= 1.28 ) %>% 
  mutate(prop = n / sum(n) * 100) %>% group_by(Random) %>% mutate(sum_rand = sum(n),prop_rand = sum(prop)) %>% 
  filter(Random)
re %>% filter(SES == "SES_min") %>% count(Random = abs(RS_SES) >= 1.28, Robust = abs(DS_SES) <= 1.28 ) %>% 
  mutate(prop = n / sum(n) * 100) %>% group_by(Random) %>% mutate(sum_rand = sum(n),prop_rand = sum(prop)) %>% 
  filter(Random)
re %>% filter(SES == "SES_max") %>% count(Random = abs(RS_SES) >= 1.28, Robust = abs(DS_SES) <= 1.28 ) %>% 
  mutate(prop = n / sum(n) * 100) %>% group_by(Random) %>% mutate(sum_rand = sum(n),prop_rand = sum(prop)) %>% 
  filter(Random)


re %>% mutate(LifeForm = sub("Epiphyte/Climber","Epiphyte\nClimber",LifeForm)) %>% 
  filter(!is.na(LifeForm)) %>% 
  mutate(Trend = ifelse(obs_mid <= 0,"Downslope","Upslope")) %>% 
  ggplot(aes(x = obs_mid,y=LifeForm,alpha = ifelse(RS_SES < 1.28,"|SES| < 1.28","|SES| > 1.28"),side = ifelse(RS_SES < 1.28,"bottom","top"))) +
ggdist::geom_dots(shape=21,linewidth=0.4,col="black",fill="black") +
    geom_hline(yintercept = 1:4,linewidth=0.2) +
    geom_vline(xintercept = c(1.28),color="grey30", linetype = c("dashed"),linewidth=0.4) +
  scale_alpha_manual(values=c(0.4,1),name="Resampling")+
  facet_wrap(~SES) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        #strip.text.x = element_blank(),
    plot.title = element_text(family="EB Garamond",size=14,hjust = 0.1,face="bold.italic"),
    axis.text.x = element_text(family="EB Garamond",size=12),
    axis.title.x = element_text(family="EB Garamond",size=12,face="bold"),
    axis.text.y = element_text(family="EB Garamond",size=12,face="bold",angle=90,hjust = 0.5),
    legend.position="bottom",
    legend.title.position = "left",
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
        legend.background = element_blank(),
        legend.key.height = unit(0.01,"cm"),
        legend.key.width = unit(0.01,"cm"),
        axis.title = element_text(family="EB Garamond"),
        axis.text = element_text(family="EB Garamond"),
        legend.spacing.y = unit(c(1.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond")) +
  labs(y="",x="Observed elevational shift") +
  scale_y_discrete(expand = c(0, 0.1)) +
  scale_x_continuous(expand = c(0, 0.1),breaks=c(-20,-10,0,10,20)) +
  NULL

ggsave(filename = "plots_tables/Supplementary_figure_LifeForm.pdf",plot = last_plot())
  
re %>% mutate(LifeForm = sub("Epiphyte/Climber","Epiphyte\nClimber",LifeForm)) %>% 
  filter(!is.na(LifeForm)) %>% 
  mutate(Trend = ifelse(obs_mid <= 0,"Downslope","Upslope")) %>% 
  ggplot(aes(x = obs_mid,y=LifeForm, alpha = ifelse(DS_SES < 1.28,"|SES| < 1.28","|SES| > 1.28"),side = ifelse(DS_SES < 1.28,"bottom","top"))) +
  ggdist::geom_dots(shape=21,linewidth=0.4,col="black",fill="black") +
  geom_hline(yintercept = 1:4,linewidth=0.2) +
  geom_vline(xintercept = c(1.28),color="grey30", linetype = c("dashed"),linewidth=0.4) +
  scale_alpha_manual(values=c(0.4,1),name="Downsampling")+
  facet_wrap(~SES) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        #strip.text.x = element_blank(),
        plot.title = element_text(family="EB Garamond",size=14,hjust = 0.1,face="bold.italic"),
        axis.text.x = element_text(family="EB Garamond",size=10),
        axis.title.x = element_text(family="EB Garamond",size=12,face="bold"),
        axis.text.y = element_text(family="EB Garamond",size=12,face="bold",angle=90,hjust = 0.5),
        legend.position="bottom",
        legend.title.position = "left",
        strip.text = element_text(family="EB Garamond",size=12,face="bold"),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.background = element_blank(),
        legend.key.height = unit(0.01,"cm"),
        legend.key.width = unit(0.01,"cm"),
        axis.title = element_text(family="EB Garamond"),
        axis.text = element_text(family="EB Garamond"),
        legend.spacing.y = unit(c(1.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond")) +
  labs(y="",x="Observed elevational shift") +
  scale_y_discrete(expand = c(0, 0.1)) +
  scale_x_continuous(expand = c(0, 0.1),breaks=c(-20,-10,0,10,20)) +
  NULL

ggsave(filename = "plots_tables/Supplementary_figure_LifeForm2.pdf",plot = last_plot())

re %>%
  filter(!is.na(Disp_Syndrome_Genus)) %>% 
  mutate(Trend = ifelse(obs_mid <= 0,"Downslope","Upslope")) %>% 
  ggplot(aes(x = obs_mid,y=Disp_Syndrome_Genus,alpha = ifelse(RS_SES < 1.28,"|SES| < 1.28","|SES| > 1.28"),side = ifelse(RS_SES < 1.28,"bottom","top"))) +
  ggdist::geom_dots(shape=21,linewidth=0.4,col="black",fill="black") +
  geom_hline(yintercept = 1:4,linewidth=0.2) +
  geom_vline(xintercept = c(1.28),color="grey30", linetype = c("dashed"),linewidth=0.4) +
  scale_alpha_manual(values=c(0.4,1),name="Resampling")+
  facet_wrap(~SES) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        #strip.text.x = element_blank(),
        plot.title = element_text(family="EB Garamond",size=14,hjust = 0.1,face="bold.italic"),
        axis.text.x = element_text(family="EB Garamond",size=12),
        axis.title.x = element_text(family="EB Garamond",size=12,face="bold"),
        axis.text.y = element_text(family="EB Garamond",size=12,face="bold",angle=90,hjust = 0.5),
        legend.position="bottom",
        legend.title.position = "left",
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.background = element_blank(),
        legend.key.height = unit(0.01,"cm"),
        legend.key.width = unit(0.01,"cm"),
        axis.title = element_text(family="EB Garamond"),
        axis.text = element_text(family="EB Garamond"),
        legend.spacing.y = unit(c(1.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond")) +
  labs(y="",x="Observed elevational shift") +
  scale_y_discrete(expand = c(0, 0.1)) +
  scale_x_continuous(expand = c(0, 0.1),breaks=c(-20,-10,0,10,20)) +
  NULL

ggsave(filename = "plots_tables/Supplementary_figure_Dispersal.pdf",plot = last_plot())

re %>% 
  filter(!is.na(Disp_Syndrome_Genus)) %>% 
  mutate(Trend = ifelse(obs_mid <= 0,"Downslope","Upslope")) %>% 
  ggplot(aes(x = obs_mid,y=Disp_Syndrome_Genus, alpha = ifelse(DS_SES < 1.28,"|SES| < 1.28","|SES| > 1.28"),side = ifelse(DS_SES < 1.28,"bottom","top"))) +
  ggdist::geom_dots(shape=21,linewidth=0.4,col="black",fill="black") +
  geom_hline(yintercept = 1:4,linewidth=0.2) +
  geom_vline(xintercept = c(1.28),color="grey30", linetype = c("dashed"),linewidth=0.4) +
  scale_alpha_manual(values=c(0.4,1),name="Downsampling")+
  facet_wrap(~SES) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        #strip.text.x = element_blank(),
        plot.title = element_text(family="EB Garamond",size=14,hjust = 0.1,face="bold.italic"),
        axis.text.x = element_text(family="EB Garamond",size=10),
        axis.title.x = element_text(family="EB Garamond",size=12,face="bold"),
        axis.text.y = element_text(family="EB Garamond",size=12,face="bold",angle=90,hjust = 0.5),
        legend.position="bottom",
        legend.title.position = "left",
        strip.text = element_text(family="EB Garamond",size=12,face="bold"),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.background = element_blank(),
        legend.key.height = unit(0.01,"cm"),
        legend.key.width = unit(0.01,"cm"),
        axis.title = element_text(family="EB Garamond"),
        axis.text = element_text(family="EB Garamond"),
        legend.spacing.y = unit(c(1.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond")) +
  labs(y="",x="Observed elevational shift") +
  scale_y_discrete(expand = c(0, 0.1)) +
  scale_x_continuous(expand = c(0, 0.1),breaks=c(-20,-10,0,10,20)) +
  NULL

ggsave(filename = "plots_tables/Supplementary_figure_Dispersal2.pdf",plot = last_plot())





i=1
cat("Reading and plotting null model",model[i],"\n")
s <- listMods[grep(model[i],listMods)]
mids = s[grep("_mid_",s)]
mins = s[grep("_min_",s)]
maxs = s[grep("_max_",s)]
load(here(mids))
esta <- e_nulls[[2]] %>% as_tibble() %>% rename("Mid_tme"=estimate)
load(here(mins))
esta <- e_nulls[[2]] %>% as_tibble() %>% rename("Min_tme"=estimate) %>% 
  bind_cols(.,esta)
load(here(maxs))
esta <- esta <- e_nulls[[2]] %>% as_tibble() %>% rename("Max_tme"=estimate) %>% 
  bind_cols(.,esta)
summary(esta$Mid_tme)
sd(esta$Mid_tme) / sqrt(500)
nulls <- esta %>% summarise(across(1:3, mean)) %>% pivot_longer(1:3,names_to = "Term",values_to = "Mean") %>% left_join(.,esta %>% summarise(across(1:3, sd)) %>% 
                              pivot_longer(1:3,names_to = "Term",values_to = "sd"),by="Term")

esta %>% ggplot(aes(x=Mid_tme)) +
  geom_histogram() +
  geom_vline(xintercept = 2.72)


load("output/LMM_models_Mid_v3.RData")
mid <- broom.mixed::tidy(modelos$cov)
load("output/LMM_models_Min_v3.RData")
min <- broom.mixed::tidy(modelos$cov)
load("output/LMM_models_Max_v3.RData")
max <- broom.mixed::tidy(modelos$cov)

mid %<>% filter(term=="tme") %>% select(estimate) %>% mutate(Term = "Mid_tme") 
min %<>% filter(term=="tme") %>% select(estimate) %>% mutate(Term = "Min_tme") 
max %<>% filter(term=="tme") %>% select(estimate) %>% mutate(Term = "Max_tme") 

DSES <- bind_rows(mid,min,max) %>% left_join(.,nulls,by="Term") %>% 
  mutate(DSES_tme = (estimate - Mean) / sd )

RSES <- bind_rows(mid,min,max) %>% left_join(.,nulls,by="Term") %>% 
  mutate(DSES_tme = (estimate - Mean) / sd )

DSES
RSES

}


## Data wrangling for land-cover and climate to ready Linear Mixed Models
process_LandCover <- function(edge= "SES_max",output="Ready_LMM_cover_max.RDS") { 
  data_lul="output/Historical_data_FULL_GLASS.v.1.R"
  nulas <- readRDS("output/SES_resamp_v.3.RDS")
  target_spp <- nulas %>% distinct(ind) %>% pull(.)
  nulas <- nulas %>% 
    mutate(SES_cat = cut(estimate,breaks=c(-Inf,-1.96,-1.28,0,1.28,1.96,Inf),labels=c("Strong","Moderate","Weak","Weak","Moderate","Strong"),include.lowest=TRUE)) 
  nulas %<>% mutate(Shift = ifelse(obs_mid < 0, "Negative","Positive")) %>% 
    mutate(SES_cat = paste0(SES_cat,"_",Shift))

load(here(data_lul))
aver <- reduct_five %>% 
  separate(Accepted_Name,into=c("G","S"),sep=" ",extra = "drop",remove=FALSE) %>% 
    unite("ind", c(G,S),sep = "_",remove=TRUE) 

aver %<>% 
    filter(ind %in% target_spp) %>% ## filter to target species
    distinct(CellID,.keep_all = TRUE) %>% # we want only one row per grid-cell
    mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(10),"Agriculture",.x))) %>%
    mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(20),"Forest",.x))) %>%
    mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(30,40),"Grassland_Shrub",.x))) %>%
    mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(70,100,90),"Other",.x))) %>%
    mutate(across(starts_with("LUC"), ~ ifelse(.x %in% c(0),NA,.x))) %>% 
    select(CellID,year,Altitude,starts_with("LUC_")) %>% 
    mutate(year=sub("Year_","",year)) %>% rename(YearCat=year) %>%
    pivot_longer(cols = -c(1:3),names_to = "year",values_to = "Values") %>%  
    separate(year,sep="_",into=c("Var","year")) %>% arrange(CellID,year) %>% 
    group_by(CellID) %>% select(-YearCat)

aver %>% mutate(New_values = ifelse(Values == "Forest" & lag(Values,n=1) != "Forest","Altered_Forest",Values)) %>% 
    mutate(New_values = ifelse(is.na(New_values),Values,New_values)) -> uno
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
  filter(Land_use %in% c("Altered_Forest","Agriculture","Forest")) %>%
  mutate(Land_use = ifelse(Land_use == "Agriculture","Altered_Forest",Land_use)) %>% 
  mutate(Land_use=factor(Land_use,levels=c("Altered_Forest","Forest"))) %>% 
  mutate(Alt_cat = cut(Altitude,breaks=seq(0,4000,1000), labels = seq(0,4000,1000)[-1],include.lowest = TRUE )) %>% filter(!is.na(Alt_cat)) %>% 
  mutate(year = as.numeric(year)) %>% 
  #mutate(year_bin = cut(year,breaks=seq(1979,2010,1), labels = seq(1979,2010,1)[-1],include.lowest = TRUE )) %>% 
  #distinct(CellID,year,.keep_all = TRUE) %>% 
  group_by(year,Alt_cat) %>%
  count(Land_use) %>% 
  ggplot(aes(x=year,y=n,fill=Land_use)) +
  geom_bar(position="fill", stat="identity") +
  facet_wrap(~Alt_cat,nrow = 1)
# We have a table with time series for every grid-cell. Now we need to join this with species occurrences

aqui <- reduct_five %>% select(CellID:Accepted_Name,year) %>% 
    #distinct(CellID,Accepted_Name,year) %>% 
    separate(Accepted_Name,into=c("G","S"),sep=" ",extra = "drop",remove=FALSE) %>% 
    unite("ind", c(G,S),sep = "_",remove=TRUE) %>% 
    filter(ind %in% target_spp) %>% 
    mutate(year=sub("Year_","",year)) %>% 
  # Join with time series
    left_join(., aver,by=c("CellID","year")) %>% 
    rename(Altitude = Altitude.x) %>% 
    mutate(year = as.numeric(sub("Year_","",year))) %>%
    rename(Land_use = Corrected_values,Elevation=Altitude) %>% 
    select(CellID,ind,year,Land_use,Elevation) %>% # Change here to keep CellID
    filter(!is.na(Land_use)) %>%
  ## Every occurrence has data on Land_use
    filter(Land_use %in% c("Altered_Forest","Agriculture","Forest")) %>%
    mutate(Land_use = ifelse(Land_use == "Agriculture","Altered_Forest",Land_use)) %>% 
    mutate(Land_use=factor(Land_use,levels=c("Altered_Forest","Forest"))) %>%
    mutate(tme = year - min(year)) %>% 
    mutate(Elevation_Cat = cut(Elevation,breaks=seq(0,4000,1000),labels=seq(0,4000,1000)[-1], include.lowest=TRUE)) %>% 
  ## This is the join with species trend data
  left_join(., nulas %>% filter(SES == edge) %>% distinct(ind,.keep_all = TRUE),by="ind")
## nulas == estimates by species
## reduct_five == all occurrences ~550k
## aver == 14743 cells, all years ~500k
## aqui == time series for all species with trend data ~315k
 saveRDS(aqui,file = here("output",output))
}
Process_climate_LMM <- function(  var = "mean_pr_", edge = "SES_min", output = "Ready_LMM_pr.RDS") {
  data="output/Historical_data_FULL_GLASS.v.1.R"
  load(here(data))
  colnames(reduct_five)
  nulas <- readRDS("output/SES_resamp_v.3.RDS")
  target_spp <- nulas %>% distinct(ind) %>% pull(.)
  nulas <- nulas %>% 
    mutate(SES_cat = cut(estimate,breaks=c(-Inf,-1.96,-1.28,0,1.28,1.96,Inf),labels=c("Strong","Moderate","Weak","Weak","Moderate","Strong"),include.lowest=TRUE)) 
  nulas %<>% mutate(Shift = ifelse(obs_mid < 0, "Negative","Positive")) %>% 
    mutate(SES_cat = paste0(SES_cat,"_",Shift))
  
  red <- reduct_five %>% 
    separate(Accepted_Name,into = c("G","S","extra"),extra = "merge",sep=" ") %>%
    select(-extra) %>% unite("ind",G:S,sep="_",remove=TRUE)
  
  aqui <- red %>% filter(ind %in% target_spp) %>%
    distinct(CellID,.keep_all = TRUE) %>% 
    select(CellID,year,starts_with(var)) %>%
    mutate(year=as.numeric(sub("Year_","",year))) %>%
    rename(YearCat=year) %>%
    pivot_longer(-c(1:2), names_to = c("Variable","year"),names_sep = sub("mean","",var),values_to = "Values") %>% 
    arrange(CellID,year) %>% 
    group_by(CellID) %>% select(-YearCat) %>% 
    rename(Climate = Values) 
  
  ## nulas == estimates by species
  ## red == all occurrences ~550k
  ## aqui == time series for all species with trend data ~ 392k
  aqui <- red %>% select(CellID:ind,year) %>%
    filter(ind %in% target_spp) %>% 
    mutate(year= sub("Year_","",year)) %>% 
    # Join with time series
    left_join(. , aqui,by=c("CellID","year")) %>% 
    mutate(year = as.numeric(year)) %>% 
    mutate(tme = year - min(year)) %>% 
    select(-Variable) %>% 
    left_join(., nulas %>% filter(SES == edge) %>% distinct(ind,.keep_all = TRUE), by="ind")
  
  saveRDS(aqui,file = here("output",output))
  
}

## Fit the LMM and plot. The function returns the plot and saves all other things. 
LandCover_LMM <- function(){
  aqui <- readRDS(here("output","Ready_LMM_cover_max.RDS"))
  aqui %<>% rename("SES_val" = estimate)
## General trend in Forest alteration
aqui %>% glmer(Land_use ~ 1 + tme + (1|ind) + (0 + tme|ind),family =  binomial, data = ., nAGQ = 0,control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4))) -> fit_lul_sum
  broom.mixed::tidy(fit_lul_sum)
  
## General trend in Forest alteration, with interaction by group of species (we can do this directly by SES or with SES_cat)
aqui %>% glmer(Land_use ~ 1 + tme + tme:SES_val + SES_val + (1|ind) + (0 + tme|ind),family =  binomial, data = ., nAGQ = 0,control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4))) -> fit_lul
  broom.mixed::tidy(fit_lul)

broom.mixed::tidy(fit_lul) %>% flextable::flextable() %>% flextable::save_as_docx(path = here("plots_tables/LMM_cover.docx"))

ggpredict(fit_lul,type = "random",terms = c("tme [all]","SES_val [-1.96,-1.28,0,1.28,1.96]"),interval = "confidence") %>% hypothesis_test() %>% tibble()

fit_lul_predict <- ggpredict(fit_lul, type = "random",terms=c("tme [all]","SES_val [-1.96,-1.28,-1.05,0,1.05,1.28,1.96]"),interval = "confidence") %>% tibble()
cover <- fit_lul_predict %>% 
    ggplot(aes(x=x+1982,y=predicted,color=group,fill=group)) +
    geom_line() +
    geom_ribbon(aes(ymin = conf.low,ymax=conf.high),alpha=0.1,color="NA") +
    scale_color_manual(name="SES",values= MetBrewer::met.brewer("Hiroshige",n=7,direction=-1)) +
    scale_fill_manual(name="SES",values= MetBrewer::met.brewer("Hiroshige",n=7,direction=-1)) +
    guides(color = guide_legend(nrow=1)) +
  guides(fill = guide_legend(nrow=1)) +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.byrow = TRUE,
          panel.background = element_blank(),
          axis.line = element_line(),
          legend.background = element_blank(),
          axis.title = element_text(family="EB Garamond"),
          axis.text = element_text(family="EB Garamond"),
          #axis.text.x = element_blank(),
          legend.spacing.y = unit(c(1.2),"cm"),
          legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond")) +
    labs(x="",y="Forest cover (probability)") +
    NULL
return(cover)
  
#ggsave(filename = "plots_tables/Figure_land_use_glmm.pdf",plot=plot_glmm)

  
}
Climate_LMM <- function(variable = "cmi", infile = "output/Ready_LMM_cmi.RDS",ylabel="XX") {
  aqui <- readRDS(here(infile))
aqui %<>% rename("SES_val" = estimate)
if(grepl("temp",variable)) {aqui %<>% mutate(across(starts_with("Climate"), function(x) x*.1 - 273.15))}
if(grepl("prec",variable)) {aqui %<>% mutate(across(starts_with("Climate"), function(x) x * .1))}
if(grepl("cmi",variable)) {aqui %<>% mutate(across(starts_with("Climate"), function(x) x * 1))}

# General climate trend
lm(Climate ~ 1 + tme, data = aqui %>% distinct(CellID,year,.keep_all = TRUE) ) -> fit1
summary(fit1)

# Linear mixed model with random effects by species
lmer(Climate ~ 1 + tme + tme:SES_val + SES_val + (1|ind), 
     data = aqui, control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5))) -> fit_2
broom.mixed::tidy(fit_2)
name <- sub("output/Ready_","",infile) %>% sub("RDS","",.)
name <- paste0("plots_tables/",name,"docx")
broom.mixed::tidy(fit_2) %>% flextable::flextable() %>% 
  flextable::save_as_docx(path = here(name))

ggpredict(fit_2, type = "random",terms=c("tme [all]", "SES_val [-1.96,-1.28,-1.05,0,1.05,1.28,1.96]"),interval = "confidence") %>% tibble() -> predict_fit
plot <- predict_fit %>% 
    ggplot(aes(x=x + 1979,y=predicted,color=group,fill=group)) +
    geom_line() +
    geom_ribbon(aes(ymin = conf.low,ymax=conf.high),alpha=0.1,color="NA") +
    scale_color_manual(values = MetBrewer::met.brewer("Hiroshige",n=7,direction = -1,),name="SES") +
  scale_fill_manual(values = MetBrewer::met.brewer("Hiroshige",n=7,direction = -1 ),name="SES") +
  guides(color = guide_legend(nrow=1)) +
  guides(fill = guide_legend(nrow=1)) +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.byrow = TRUE,
          panel.background = element_blank(),
          axis.line = element_line(),
          legend.background = element_blank(),
          axis.title = element_text(family="EB Garamond"),
          axis.text = element_text(family="EB Garamond"),
          #axis.text.x = element_blank(),
          legend.spacing.y = unit(c(1.2),"cm"),
          legend.title = element_text(family="EB Garamond"),
          legend.text = element_text(family="EB Garamond")) +
    labs(x="",y=ylabel) +
    NULL
return(plot)




}

## We need to save all plots. Uses objects generated above.
plot_all_factors <- function(A = tas, B = prec, C = cmi, D = cover) { 
  layout = "AABB\n#EE#\nCCDD"
  A + B + C + D + guide_area() +
  plot_layout(design = layout, nrow = 3, guides = "collect", ncol=2, tag_level = 'new',axes = "collect", heights = c(10,3,10)) +
    #guides(SES = guide_legend(nrow=7,ncol=1,direction = "horizontal")) +
    plot_annotation(tag_levels = 'a') & 
      theme(
      legend.direction = "horizontal",
            legend.text.position = "top",
            legend.byrow = TRUE,
      plot.title = element_text(family="EB Garamond",size=10,face="bold"),
      plot.tag = element_text(family="EB Garamond",size=16,face="bold"),
      plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
      legend.key.height=unit(0.5,"cm"),
      legend.key.width = unit(0.3,"cm"),
      legend.text=element_text(size=8,family="EB Garamond"), 
      legend.title=element_text(size=10,family="EB Garamond"),
      axis.title = element_text(family="EB Garamond",size=10))
    
    
  ggsave(filename = "plots_tables/Supplementary_figure_LAST.pdf",plot = last_plot())
  
}  
##### EXPERIMENTAL #####

ClimbyCover <- function(cover_ready = "Ready_LMM_cover_max.RDS", climate_ready = "Ready_LMM_cmi.RDS",paleta = "Isfahan1"){ #Ingres Greek
alla <- readRDS(here("output",cover_ready))
aqui <- readRDS(here( "output",climate_ready))
aculla <- aqui %>% left_join(.,alla %>% distinct(CellID,year,.keep_all = TRUE) %>% select(CellID,year,Land_use),by=c("CellID","year")) 
var_id <- sub("Ready_LMM_","",climate_ready)
aculla %<>% select(Climate,tme, CellID,year,ind,Land_use) %>% filter(!is.na(Land_use))
aculla %>% distinct(CellID,year,.keep_all = TRUE) %>% 
  lmer(Climate ~ 1 + tme:Land_use  + Land_use + (1|ind), 
       data = ., control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5))) -> fit_2


ggpredict(fit_2, type = "random",terms=c("tme [1979,1980,1990,2000,2010]", "Land_use"),interval = "confidence") %>% tibble() -> predict_fit

a <- predict_fit %>% 
  ggplot(aes(x=x + 1979,y=predicted,color=group,fill=group)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low,ymax=conf.high),alpha=0.1,color="NA") +
  scale_color_manual(values = MetBrewer::met.brewer(paleta,n=7,direction = 1,)[c(2,6)],name="SES") +
  scale_fill_manual(values = MetBrewer::met.brewer(paleta,n=7,direction = 1 )[c(2,6)],name="SES") +
  guides(color = guide_legend(nrow=1)) +
  guides(fill = guide_legend(nrow=1)) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.byrow = TRUE,
        panel.background = element_blank(),
        axis.line = element_line(),
        legend.background = element_blank(),
        axis.title = element_text(family="EB Garamond"),
        axis.text = element_text(family="EB Garamond"),
        #axis.text.x = element_blank(),
        legend.spacing.y = unit(c(1.2),"cm"),
        legend.title = element_text(family="EB Garamond"),
        legend.text = element_text(family="EB Garamond")) +
  labs(x="",y=var_id) +
  NULL

return(a)
}

plot_ClimbyCover <- function(A = temp, B = prec, C = cmi){
  layout = "AA\nBB\nCC"
A + B + C +
  plot_layout(design = layout, nrow = 3, guides = "collect", 
              ncol=1, tag_level = 'new',axes = "collect") +
  plot_annotation(tag_levels = 'a') & 
  theme(legend.position = "right",
    legend.direction = "vertical",
    legend.text.position = "top",
    legend.byrow = TRUE,
    plot.title = element_text(family="EB Garamond",size=10,face="bold"),
    plot.tag = element_text(family="EB Garamond",size=16,face="bold"),
    plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
    legend.key.height=unit(0.5,"cm"),
    legend.key.width = unit(0.3,"cm"),
    legend.text=element_text(size=8,family="EB Garamond"), 
    legend.title=element_text(size=10,family="EB Garamond"),
    axis.title = element_text(family="EB Garamond",size=10))

ggsave(filename = "plots_tables/Supplementary_figure_4.pdf",plot = last_plot())
}

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
  load(here(data)) ### occurrences
  load(data_spp) ### summary (vs. nulls)
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

## Not sure if this one is seriated with other functions, need to check and modify to make it stand alone
do_venn <- function(){
  esta <- readRDS("output/LMM_trends_v.3.RDS")
  data="output/Historical_data_FULL_GLASS.v.1.R"
  load(here(data))
  reduct_five %<>% separate(Accepted_Name,into = c("G","S","extra"),extra = "drop",sep=" ") %>% 
    unite("ind",G:S,sep="_") %>% right_join(.,esta,by="ind")
  
  reduct_five %>% select(CellID:ind,year,starts_with("mean_tas_")) %>%
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
  
  nulas <- readRDS("output/SES_resamp_v.3.RDS")
  nulas <- nulas %>% 
    mutate(SES_cat = cut(estimate,breaks=c(-Inf,-1.96,-1.28,0,1.28,1.96,Inf),labels=c("Strong","Moderate","Weak","Weak","Moderate","Strong"),include.lowest=TRUE)) 
  nulas %<>% mutate(Shift = ifelse(obs_mid < 0, "Negative","Positive")) %>% 
    mutate(SES_cat = paste0(SES_cat,"_",Shift))
  nulas <- readRDS("output/SES_resamp_v.3.RDS")
  nulas <- nulas %>% 
    mutate(SES_cat = cut(estimate,breaks=c(-Inf,-1.96,-1.28,0,1.28,1.96,Inf),labels=c("Strong","Moderate","Weak","Weak","Moderate","Strong"),include.lowest=TRUE)) 
  nulas %<>% mutate(Shift = ifelse(obs_mid < 0, "Negative","Positive")) %>% 
    mutate(SES_cat = paste0(SES_cat,"_",Shift))
  nulas %<>% distinct(ind,.keep_all = TRUE)
  
  left_join(aqui,nulas %>% select(ind,SES,estimate,SES_cat,Shift),by="ind") %>% 
    lmer(Climate ~ 1 + tme + tme:SES_cat + (1|ind), data=., control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5))) -> fit_aver
  
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
  
  nulas
  
  nulas %>% select(SES,ind,estimate) %>% distinct(SES, ind,.keep_all = TRUE) %>% 
    pivot_wider(id_cols = "ind",values_from = "estimate",names_from = "SES") %>% 
    mutate(across(2:4,\(x) cut(x,c(-Inf,-1.96,-1.28,0,1.28,1.96,Inf),labels=c("Negative","Negative","No","No","Positive","Positive"),include.lower =TRUE))) %>% group_by(SES_mid,SES_min,SES_max) %>% 
    mutate(Pre = 1) %>% 
    pivot_wider(id_cols = ind,names_from = starts_with("SES_"),values_from = "Pre",values_fill = 0)
  
  # Upslope shifts
  x <- list(
    MidPoint = a %>% filter(SES_mid == "Positive") %>% pull(ind),
    ColdEdge = a %>% filter(SES_max != "No") %>% pull(ind),
    WarmEdge = a %>% filter(SES_min != "No") %>% pull(ind))
  ggvenn(x)
  
  # Upslope with retreating warm edged
  x <- list(
    MidPoint = a %>% filter(SES_mid == "Positive") %>% pull(ind),
    ColdEdge = a %>% filter(SES_max != "No") %>% pull(ind),
    WarmEdge = a %>% filter(SES_min == "Positive") %>% pull(ind))
  ggvenn(x)
  
  # Download shifts
  x <- list(
    MidPoint = a %>% filter(SES_mid == "Negative") %>% pull(ind),
    ColdEdge = a %>% filter(SES_max != "No") %>% pull(ind),
    WarmEdge = a %>% filter(SES_min != "No") %>% pull(ind))
  ggvenn(x)
  
  # Download with retracting cold edged
  x <- list(
    MidPoint = a %>% filter(SES_mid == "Positive") %>% pull(ind),
    ColdEdge = a %>% filter(SES_max == "Negative") %>% pull(ind),
    WarmEdge = a %>% filter(SES_min != "No") %>% pull(ind))
  ggvenn(x)   
  
  # Contractions
  x <- list(
    ColdEdge = a %>% filter(SES_max == "Negative") %>% pull(ind),
    WarmEdge = a %>% filter(SES_min == "Positive") %>% pull(ind))
  ggvenn(x)
  
  # Expansions
  x <- list(
    #MidPoint = a %>% filter(SES_mid == "Negative") %>% pull(ind),
    ColdEdge = a %>% filter(SES_max == "Positive") %>% pull(ind),
    WarmEdge = a %>% filter(SES_min == "Negative") %>% pull(ind))
  ggvenn(x)
}

plot_baselines <- function(){
  data="output/LMM_trends_v.3.RDS"
  downsample="output/SummaryLMM_downsamp_mid_v3.RData"
  resample = "output/SummaryLMM_resamp_mid_v3.RData"
  model_data <- readRDS(here(data))
  load(here(downsample))
  down <- e_nulls[[1]]
  load(here(resample))
  resamp <- e_nulls[[1]]
  dd <- down %>% group_by(ind) %>% summarise(mean=mean(estimate),sd=sd(estimate)) %>% mutate(type="down")
  dd <- resamp %>% group_by(ind) %>% summarise(mean=mean(estimate),sd=sd(estimate)) %>% mutate(type="resamp") %>% rbind(.,dd)
  dd %<>% filter(ind!="Frangula capreifolia var. grandifolia (M.C.Johnst. & L.A.Johnst.) A.Pool") %>% filter(ind!="Styrax argenteus var. ramirezii (Greenm.) Gonsoulin") 
  dd %<>% separate(ind,into = c("G","S","extra"),extra = "merge",sep=" ") %>% select(-extra) %>% unite("ind",G:S,sep="_")
  dd %<>% right_join(.,model_data,by="ind")
  
  data="output/Historical_data_FULL_GLASS.v.1.R"
  load(here(data))
  reduct_five %<>% separate(Accepted_Name,into = c("G","S","extra"),extra = "merge",sep=" ") %>% select(-extra) %>% unite("ind",G:S,sep="_")
  base_eco <- reduct_five %>% select(CellID:ind,year,starts_with("min_tas_"),starts_with("mean_tas_"),starts_with("max_tas_")) %>%  mutate(year=as.numeric(sub("Year_","",year))) %>% 
    filter(year<=1990) %>% 
    distinct(CellID,ind,.keep_all = TRUE) %>% 
    rowwise() %>% mutate(baseline_tas = mean(c_across(mean_tas_1979:mean_tas_1990),na.rm=TRUE),.before=year) %>% mutate(baseline_tasmin = mean(c_across(min_tas_1979:min_tas_1990),na.rm=TRUE),.before=year) %>% mutate(baseline_tasmax = mean(c_across(max_tas_1979:max_tas_1990),na.rm=TRUE),.before=year)
  
  last_eco <- reduct_five %>% select(CellID:ind,year,starts_with("min_tas_"),starts_with("mean_tas_"),starts_with("max_tas_")) %>% mutate(year=as.numeric(sub("Year_","",year))) %>% 
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
  
  
  nulas <- dd %>% mutate(SES = (trend_mid - mean) / sd,.before="mean")
  nulas %<>% select(ind, SES, type,trend_mid) %>% pivot_wider(id_cols = c(ind,trend_mid),values_from = "SES",names_from = "type",names_prefix = "SES_") %>% mutate(Shift = ifelse(trend_mid < 0, "Downslope","Upslope"))
  to_plot_tas <-  base_eco %>% rename(TSI_base = TSI) %>% left_join(.,last_eco %>% rename(TSI_last = TSI),by="ind") %>% 
    left_join(., nulas %>% select(ind,SES_resamp,SES_down,Shift),by="ind")
  to_plot_tas %<>% 
    mutate(SES_resamp_cat = cut(SES_resamp,breaks=c(-Inf,-1.96,-1.28,0,1.28,1.96,Inf),labels=c("Strong","Moderate","Weak","Weak","Moderate","Strong"),include.lowest=TRUE)) %>% 
    unite(col = "SES_resamp_cat", c("SES_resamp_cat","Shift"), sep="\n",remove = FALSE) %>% 
    mutate(SES_resamp_cat = factor(SES_resamp_cat,levels= c("Strong\nDownslope","Moderate\nDownslope","Weak\nDownslope","Weak\nUpslope","Moderate\nUpslope","Strong\nUpslope")))
  to_plot_tas %<>% filter(!is.na(SES_resamp))
  
  cont <- to_plot_tas %>% 
    ggplot(aes(x= SES_resamp,y = baseline_tas,color=Shift, fill=Shift)) +
    geom_smooth(method="lm", color="grey50",fill="grey50") +
    geom_point(alpha=0.7,size=2) +
    theme(legend.position = "top") +
    scale_color_manual(values = MetBrewer::met.brewer("Ingres",5)[c(1,3)],name="") +
    scale_fill_manual(values = MetBrewer::met.brewer("Ingres",5)[c(1,3)],name="") +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text.y = element_text(size=10),
          axis.title = element_text(size=14),
          axis.line = element_line(),
          text=element_text(family="EB Garamond")) +
    labs(x="Standardised Effect Size (SES)",y="Baseline mean temperature (ºC)") +
    geom_vline(xintercept = c(-1.28,0,1.28),type="dashed",alpha=0.5) +
    NULL
  
  splines="output/Splines_v.3.Rdata"
  load(here(splines))
  pattern = NULL
  end_year = 2010
  aver <- mean_alts %>% mutate(Type="OBS",.before=data) %>% select(2:3) %>% unnest(data)
  if(is.null(pattern)) {
    pat = "RollAlt"
    pat_sd = paste0(pat,"_sd")
    names <- c("ind","intercept_mid","trend_mid")
  }
  if(!is.null(pattern)) { 
    pat = paste0("RollAlt_",pattern)
    pat_sd = paste0(pat,"_sd")
    names <- c("ind",paste0(c("intercept_","trend_"),pattern))
  }
  aver %<>% ungroup() %>% 
    select(Accepted_Name,YearCat,ends_with(pat),Rolln,starts_with(pat_sd)) %>% 
    mutate(YearCat=as.numeric(YearCat)) %>% 
    mutate(RollAlt_error = unlist(across(starts_with(pat_sd), ~ .x / sqrt(Rolln)))) %>% 
    select(ends_with(pat),YearCat,Accepted_Name,RollAlt_error) %>% rename_with(~c("y","tme","ind","error"))
  aver %<>% filter(tme <= end_year) %>% mutate(year=tme) %>% mutate(tme = tme - start_year)
  otra <- aver %>% separate(ind,into = c("G","S","extra"),extra = "merge",sep=" ") %>% select(-extra) %>% unite("ind",G:S,sep="_")
  otra <- otra %>% filter(year<=2010)
  to_merge <- otra %>% group_by(ind) %>% summarise(error = median(error,na.rm=TRUE))
  fit_try = TRUE
  if(fit_try){ 
    target_try <- readRDS("output/TRY_processed.RDS")
    otra %<>% left_join(., target_try %>% rename("ind"=binomial) %>% select(ind,Disp_Syndrome,Disp_Syndrome_Inf,LifeForm),by="ind")
    
    # Do a bit of wrangling of the try data to get inferred dispersal syndrome by genus
    target_try %<>% select(binomial,Disp_Syndrome_Inf) %>% 
      rename("Disp_Syndrome_Genus" = Disp_Syndrome_Inf) %>% 
      separate(binomial,c("Genus","spp"),sep="_",remove = FALSE) %>% 
      distinct(Genus,.keep_all = TRUE) %>% select(-binomial,-spp)
    otra %<>% separate(ind,c("Genus","spp"),sep="_",remove=FALSE) %>% 
      left_join(., target_try,by="Genus" ) %>% select(-Disp_Syndrome_Inf)
    
    otra %<>% mutate(Disp_Syndrome = replace_na(Disp_Syndrome,"Unknown")) %>% 
      mutate(LifeForm = replace_na(LifeForm,"Unknown")) %>% 
      mutate(Disp_Syndrome = ifelse(Disp_Syndrome == "Water","Unknown",Disp_Syndrome))
    otra %<>% mutate(LifeForm= factor(LifeForm,levels=c("Unknown","Herb","Shrub","Epiphyte/Climber","Tree"))) %>%
      mutate(Disp_Syndrome= factor(Disp_Syndrome,levels=c("Unknown","Unassisted","Animal","Wind")))
  }
  
  fit0 <-lm(SES_resamp ~ 1 + baseline_tas + baseline_tas:LifeForm + LifeForm,data = to_plot_tas %>% left_join(.,otra %>% distinct(ind,.keep_all = TRUE),by="ind"))
  broom::tidy (fit0) %>% flextable::flextable() %>% flextable::save_as_docx(path = here("plots_tables/baseline_lm_gf.docx"))
  
  growth <- to_plot_tas %>% left_join(.,otra %>% distinct(ind,.keep_all = TRUE),by="ind") %>% 
    ggplot(aes(x= SES_resamp,y = baseline_tas, color= LifeForm,fill=LifeForm)) +
    geom_point(alpha=0.4,size=1) +
    theme(legend.position = "top") +
    scale_fill_manual(values = MetBrewer::met.brewer("Archambault",5),name="") +
    scale_color_manual(values = MetBrewer::met.brewer("Archambault",5),name="") +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text.y = element_text(size=10),
          axis.title = element_text(size=14),
          axis.line = element_line(),
          legend.key.size = unit(0.01,"cm"),
          legend.spacing.x = unit(0.01,"cm"),
          text=element_text(family="EB Garamond")) +
    geom_smooth(method="lm", alpha=0.5) +
    labs(x="Standardised Effect Size (SES)",y="Baseline mean temperature (ºC)") +
    geom_vline(xintercept = c(-1.28,0,1.28),type="dashed",alpha=0.5) +
    NULL
  fit0 <-lm(SES_resamp ~ 1 + baseline_tas + baseline_tas:Disp_Syndrome + Disp_Syndrome,data = to_plot_tas %>% left_join(.,otra %>% distinct(ind,.keep_all = TRUE),by="ind"))
  broom::tidy (fit0) %>% flextable::flextable() %>% flextable::save_as_docx(path = here("plots_tables/baseline_lm_ds.docx"))
  
  disp <- to_plot_tas %>% left_join(.,otra %>% distinct(ind,.keep_all = TRUE),by="ind") %>% 
    ggplot(aes(x= SES_resamp,y = baseline_tas, color= Disp_Syndrome,fill=Disp_Syndrome)) +
    geom_point(alpha=0.4,size=1) +
    theme(legend.position = "top") +
    scale_fill_manual(values = MetBrewer::met.brewer("Egypt",4),name="") +
    scale_color_manual(values = MetBrewer::met.brewer("Egypt",4),name="") +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text.y = element_text(size=10),
          axis.title = element_text(size=14),
          axis.line = element_line(),
          legend.key.size = unit(0.01,"cm"),
          legend.spacing.x = unit(0.01,"cm"),
          text=element_text(family="EB Garamond")) +
    geom_smooth(method="lm", alpha=0.5) +
    labs(x="Standardised Effect Size (SES)",y="Baseline mean temperature (ºC)") +
    geom_vline(xintercept = c(-1.28,0,1.28),type="dashed",alpha=0.5) +
    NULL
  
  fit0 <-lm(SES_resamp ~ 1 + baseline_tas,data = to_plot_tas %>% left_join(.,otra %>% distinct(ind,.keep_all = TRUE),by="ind"))
  broom::tidy (fit0) %>% flextable::flextable() %>% flextable::save_as_docx(path = here("plots_tables/baseline_lm.docx"))
  
  
  general  <- to_plot_tas %>%  
    ggplot(aes(x= SES_resamp_cat,y = baseline_tas)) +
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
    geom_text(data = to_plot_tas  %>% group_by(SES_resamp_cat) %>% 
                summarise(med=median(baseline_tas)),family="EB Garamond",fontface="bold",size=5,
              aes(x=SES_resamp_cat,y=med-9,label=paste0(round(med,2)))) +
    geom_text(data= to_plot_tas  %>%
                group_by(SES_resamp_cat) %>% mutate(med=median(baseline_tas)) %>% 
                summarize(med=first(med),n=n()),family="EB Garamond",
              aes(x=SES_resamp_cat,y=med-10,label=paste0("n=",round(n,2)))) +
    labs(x="",y="Baseline mean temperature (ºC)") +
    stat_summary(geom="crossbar",fun=function(x) quantile(x,probs=c(0.5)),width=0.2,size=0.4) +
    stat_summary(geom="crossbar",fun=function(x) quantile(x,probs=c(0.75,0.25)),width=0.18,size=0.2) +
    NULL
  
  
  layout <- "
AAAAAA
BBCCDD
"
  pp <- general + cont + growth + disp +
    plot_layout(design = layout,byrow=T,nrow = 2,guides="keep",axes="collect_x",tag_level = 'new') + 
    plot_annotation(tag_levels = 'a',title="") & 
    theme(
      plot.tag=element_text(family="EB Garamond",size=15,face="bold"),
      plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"), 
      plot.title = element_text(family="EB Garamond",size=16,face="bold"),
      legend.key.height=unit(0.7,"cm"),
      legend.key.size = unit(0.01,"cm"),
      legend.spacing.x = unit(0.01,"cm"),
      legend.text=element_text(size=8,family="EB Garamond"), 
      legend.title=element_text(size=10,family="EB Garamond"),
      axis.title = element_text(family="EB Garamond",size=10))
  pp
  
  ggsave(filename=here("plots_tables/Shifts_by_baseline.pdf"),plot=pp, device = "pdf",width = 25,units = "cm")
  
  
  
  vnulas <- dd %>% mutate(ES = (trend_mid - mean) / sd,.before="mean")
  nulas %<>% select(ind, ES, type,trend_mid) %>% pivot_wider(id_cols = c(ind,trend_mid),values_from = "ES",names_from = "type",names_prefix = "ES_") %>% mutate(Shift = ifelse(trend_mid < 0, "Downslope","Upslope"))
  to_plot_tas %<>% select(ind,SES_resamp) %>%   left_join(., nulas %>% select(ind,ES_resamp,trend_mid,Shift),by="ind")
  to_plot_tas %<>% 
    mutate(SES_resamp_cat = cut(SES_resamp,breaks=c(-Inf,-1.96,-1.04,0,1.04,1.96,Inf),labels=c("Strong","Strong","Weak","Weak","Strong","Strong"),include.lowest=TRUE)) %>% 
    #mutate(SES_down_SES_cat = cut(SES_down,breaks=c(-Inf,-1.96,-1.28,0,1.28,1.96,Inf),labels=c("Sensitive","Sensitive","Robust","Robust","Sensitive","Sensitive"),include.lowest=TRUE)) %>% 
    unite(col = "SES_resamp_cat", c("SES_resamp_cat","Shift"), sep="\n",remove = FALSE) #%>% 
  # mutate(SES_resamp_cat = factor(SES_resamp_cat,levels= c("Strong\nDownslope","Moderate\nDownslope","Weak\nDownslope","Weak\nUpslope","Moderate\nUpslope","Strong\nUpslope")))
  
  to_plot_tas %>% group_by(SES_resamp_cat) %>% 
    summarise(magnitude = mean(ES_resamp),n=n(),sd=sd(ES_resamp))
  
  to_plot_tas %>% filter(!grepl("Weak",SES_resamp_cat)) %>% 
    aov(abs(SES) ~ SES_resamp_cat,data=.) %>% TukeyHSD()
  
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
    mutate(DSamp_SES_cat = cut(DSamp_SES,breaks=c(-6,-1.96,-1.28,0,1.28,1.96,6),labels=c("Strong_Neg","Moderate_Neg","Weak_Neg","Weak_Pos","Moderate_Pos","Strong_Pos"),include.lowest=TRUE)) %>%
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

manken_maps <- function(){ 
  data="output/Historical_data_FULL_GLASS.v.1.R"
  load(here(data))
  load(here("output/LMM_modelData.Rdata"))
  load(here("interim/downSummary_LMM_Mid_v3.Rdata"))
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















