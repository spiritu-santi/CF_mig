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

#### FUNCTIONS AND SETTINGS
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
dif_clim <- function(x) {lapply(x,function(x) {
  x %>% filter(YearCat < 1976) %>% summarise(across(c(-1),mean,na.rm=T)) -> xx
  simple <- function(z,y) {z-y}
  map2(x[,-1],xx,simple) %>% as_tibble() %>% mutate(YearCat=x$YearCat)
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

lulc_extract <- function(data="Historical_dataClimate_full.R",lulc_data = "LULC_850_2015") { 
load(here(paste("output_data",data,sep="/")))
lista_luc <- list.files(path = here(paste("input_data/","LULC_850_2015",sep="")),pattern = ".tif$",full.names = T) %>% .[grep("1901_",.):grep("2015_",.)]
lucs <- raster::stack(lista_luc)
names(lucs) <- lista_luc %>% strsplit(.,"LULC_") %>% lapply(.,"[",3) %>% sub("_.*","",.) %>% paste("LUC",.,sep="_")
xy <- reduct_five %>% dplyr::select(decimallongitude,decimallatitude)
raster::extract(lucs,xy) -> lucs_ext
lucs_tib <- as_tibble(lucs_ext) %>% bind_cols(reduct_five,.)
save(lucs_tib,file = here("output_data/Historical_dataLULC.R"))
}


##### CLIMATE TRENDS ##### 
# now to estimate the climate trends across cloud forest localities.
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

#### PLOT ####
# now plot the actual values 
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
if(prec) {direct=-1; opcion_viridis ="D"}
if(!prec) {direct=1; opcion_viridis="A"}
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
#####





#### PLOTS: STAND ALONE PLOTS, WHICH MAKES IT A BIT SLOW
plot_Chords <- function(lul = "Historical_dataLULC.R",alt_breaks = c(0,300,1000,2000,4000),before = c("LUC_1936","LUC_1956","LUC_1976","LUC_1996"), after = c("LUC_1956","LUC_1976","LUC_1996","LUC_2015"),cat = c("Low","MidLow","MidHigh","High"),codes = c("PF","AU","NF","SV"),targetLU="PF") {
  load(here(paste("output_data",lul,sep="/")))
lucs_tib %>% filter(decimallatitude > lat_limits[1] & decimallongitude < lon_limits[1],
      !(decimallatitude > lat_limits[2] | decimallongitude < lon_limits[2]),
      !(decimallatitude > lat_limits[3] & decimallongitude > lon_limits[3]),
      !CellID%in%all_of(non_id),AltRas<alt_limit) %>% dplyr::select(CellID,x,y,AltRas,starts_with("LUC")) %>% distinct(CellID,.keep_all = T) %>% droplevels() %>% arrange(CellID) -> precs
precs %>% pull(AltRas) %>% cut(.,breaks=seq(0,alt_limit+100,belts),labels=seq(0,alt_limit+100,belts)[-1], include.lowest=T) -> alt_cats ### aggregate records into altitudinal belts
precs <- precs %>% ungroup() %>% mutate(AltCats = alt_cats)
# Estimate means across years.
for (j in 1:4){ 
pdf(here(paste("plots/LUL_transition_",paste(cat[i],"_",sub("LUC_","",before[j]),"-",sub("LUC_","",after[j]),sep=""),".pdf",sep="")))
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
  theme(legend.position = "bottom",panel.background = element_blank(),panel.grid = element_blank(),axis.ticks = element_blank(),axis.text.x = element_text(size=10),axis.line = element_line()) + labs(y="Frequency",title="Land use in cloud forests through time") +
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
plot_ClimSeries <- function(data="output_data/ClimateRAW_by_AltCat.R",inflect_point = 1976,output="plots/ClimateSeries_TempRAW.pdf",variable = "tm",viridis_option = "F") {
  load(here(data))
  aa %>% mutate(Shift=case_when(YearCat < inflect_point ~ "Pre",YearCat >= inflect_point ~ "Post")) %>% group_by(Shift,AltCats) %>% summarise(across(c(2:9),mean)) -> means
  plotas <- list()
  colNames <- names(aa)[grep("tm",names(aa))]
  for (i in 1:length(colNames)){
    if(i %in% c(1,3)){
      plotas[[i]] <- aa %>% nest() %>% mutate(data = (dif_clim(data))) %>% unnest(cols = data) %>% 
        ggplot(aes_string(x=aa$YearCat, y=aa$AltCats, fill= colNames[i])) + 
        geom_tile() + labs(subtitle = colNames[i], y="Altitude",x="") + 
        theme(panel.background = element_blank(),axis.ticks.y = element_blank(),legend.position = "right",legend.key.height = unit(0.8,"cm"),legend.key.width = unit(0.2,"cm"),legend.text = element_text(size=6),legend.margin = margin(0,0,0,0)) + ylim(0,alt_limit+200) +
        scale_fill_viridis_b(n.breaks=20,option=viridis_option,direction= 1,name="")
    }
    if(i %in% c(2,4)){
      plotas[[i]] <- aa %>% nest() %>% mutate(data = (dif_clim(data))) %>% unnest(cols = data) %>% 
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
  ggsave(filename = output,plot = pp)
}



