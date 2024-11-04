library(forecast)
bld.mbb.bootstrap.stationary <- function (x, num, n_b = 20, span=0.5, block_size = NULL, mean_block=10, stationary = TRUE) {
  if (length(x) <= 1L) 
    return(rep(list(x), num))
  freq <- frequency(x)
  if (length(x) <= 2 * freq) 
    freq <- 1L
  if (is.null(block_size)) {
    block_size <- ifelse(freq > 1, 2 * freq, min(8, floor(length(x)/2)))
  }
  xs <- list()
  xs[[1]] <- x
  if (num > 1) {
    if (min(x) > 1e-06) {
      lambda <- BoxCox.lambda(x, lower = 0, upper = 2)
    }
    else {
      lambda <- 1
    }
    x.bc <- BoxCox(x, lambda)
    lambda <- attr(x.bc, "lambda")
    if (freq > 1) {
      x.stl <- stl(ts(x.bc, frequency = freq), "per")$time.series
      seasonal <- x.stl[, 1]
      trend <- x.stl[, 2]
      remainder <- x.stl[, 3]
    }
    else {
      trend <- 1:length(x)
      suppressWarnings(x.loess <- loess(ts(x.bc, frequency = 1) ~ trend, span = span, degree = 1))
      seasonal <- rep(0, length(x))
      trend <- x.loess$fitted
      remainder <- x.loess$residuals
    }
  }
  if(stationary == TRUE){
    for (i in 2:num) {
      xs[[i]] <- ts(InvBoxCox(trend + seasonal + tsbootstrap(remainder, nb = n_b, b =mean_block, type="stationary"), lambda))
      tsp(xs[[i]]) <- tsp(x)
    }
  }
  if(stationary == FALSE){
    for (i in 2:num) {
      xs[[i]] <- ts(InvBoxCox(trend + seasonal + MBB(remainder, block_size), lambda))
      tsp(xs[[i]]) <- tsp(x)
    }
  }   
  
  xs
}
loess.apply <- function(x) {
  loess(ts(x) ~ c(1:length(x)), span=span,degree=1,family = "gaussian")$fitted
}
do_lms <- function(data) { 
  data %>% summarise(across(c(y,roll_y,SBB_mean), \(x) lm( x ~ tme, data=.) %>% 
                              broom::tidy() %>% filter(term=="tme") %>% pull(2), .names = "{.col}")) %>% 
    pivot_longer(cols = everything(),names_to = "LM_estimate", values_to = "Value")
}
stationary_bootstrap <- function(){
  splines="output/Splines_v.3.Rdata"
  load(here(splines))
  pattern = NULL
  start_year=1979
  end_year=2010
  plot = FALSE
  span = 0.17
  num_bs = 1000
  mean_alts %<>% mutate(Boostraps = list(NULL))
  species <- mean_alts %>% mutate(Accepted_Name = as.character(Accepted_Name)) %>% pull(Accepted_Name)
  for(i in 1:length(species)) {
    cat("Bootstrapping for",i,"--",species[i],"\n")
    aver <- mean_alts %>% filter(Accepted_Name == species[i]) %>% unnest(data)
    if(is.null(pattern)) {
      pat = "Variable"
      pat2 = "RollAlt"
    }
    if(!is.null(pattern)) { 
      pat = paste0(pattern,"Variable")
      pat2 = paste0("RollAlt_",pattern)
    }
    aver %<>% ungroup() %>% 
      select(YearCat,starts_with(pat),pat2,n) %>% 
      mutate(YearCat=as.numeric(YearCat)) %>%
      select(starts_with(pat),YearCat,pat2,n) %>% rename_with(~c("y","tme","roll_y","n")) %>% 
      #filter(tme <= end_year) %>% 
      mutate(year=tme) %>% mutate(tme = tme - start_year)
    
    # There a couple of things to consider before bootstraping a time series:
    # 1) Stationary boostrapping (and other techniques) work on stationary series, but here we have non-stationary series with a trend.
    # 2) We can decompose the time series into trend and remainder components, and then apply boostrapping to the remainders (which we assume are stationary).
    # 3) The time series for species are actually 'smoothed versions' of the raw series, using the moving window mean estimates. Thus, applying the bootstrap to these series seems not appropriate.
    # 4) We have 'gappy' series for most species, where some years have no data, and others have more than one occurrence (this is the rationale behind the moving window smoothing we applied). So, given the data structure we can either re-sample/bootstrap: the original occurrence data (something similar to the downsampled models); the 'smoothed' mean series or something in between using means per year. 
    # The procedure is as follows: 
    # - Estimate annuals means of elevation from occurrences, since we need a single data point per time unit. 
    # - Fill in the gaps in the 'annual mean' time series using a 'last observation carried forward' approach. Here, missing values are given the preceding value in the series. For instance, the following series 'ts(1980 = 10, 1981 = NA, 1982 = 13)' will look like this 'ts(1980 = 10, 1981 = 10, 1982 = 13)'.
    # - Box-Cox transfor the time series and decompose by fitting a 'loess' model: trend + remainder. 
    # - Extract the residuals/remainders for each time unit, which are the values to be boostrapped.
    # - Resample the remainders using the 'Stationary Bootstrap' approach with a mean block length of 10.
    # - Reverse the Box-Cox transformation using the trend + the resampled remainders to obtain boostrapped time series.
    # - Fit 'loess' models to the boostrapped time series to obtain confidence intervals using the mean and standard deviation.
    
    aver %<>% droplevels() %>% mutate(YearDate = lubridate::ymd(year, truncated = 2L)) 
    x <- padr::pad(aver) %>% fill(y,roll_y) %>% mutate(tme = year(YearDate) - 1979)
    boots_station <- bld.mbb.bootstrap.stationary(ts(x$y), num = 2, n_b = num_bs, span = span, mean_block = 10, stationary = TRUE) %>% 
      .[[2]] %>% as_tibble(.name_repair = "universal")
    
    confs <- boots_station %>% 
      mutate(across(starts_with("Series"), \(x) loess.apply(x))) %>% 
      rowwise() %>% 
      mutate(SBB_sd = sd(c_across(starts_with("Series"))),.before=1) %>% 
      mutate(SBB_mean = mean(c_across(starts_with("Series"))),.before=1) %>% 
      #mutate(Upper = quantile(c_across(starts_with("Series")),probs=0.975)) %>% 
      #mutate(Lower = quantile(c_across(starts_with("Series")),probs=0.025)) %>% 
      ungroup() %>% 
      mutate(SBB_sem = SBB_sd / sqrt(num_bs),.before=1) %>% 
      mutate(SBB_upper = SBB_mean + (SBB_sd * 1.96 ),.before=1) %>% 
      mutate(SBB_lower = SBB_mean - (SBB_sd * 1.96 ),.before=1) %>%
      mutate(tme = x$tme,.before=1) %>% 
      mutate(n = x$n,.before=1)
    
    complete <- right_join(aver,confs,by="tme")
    mean_alts$Boostraps[[i]]<- complete
    
  }
  
  station_boot = mean_alts
  station_boot %<>% mutate(LMs = list(do_lms(Boostraps)))
  save(station_boot,file="interim/StationaryBoot_Sampling.RData")
}

bias <- function(x) {
  x %>% pivot_longer(cols = -c(y, roll_y),names_to = "Series",values_to = "Value") %>% 
    mutate(Bias = y - Value) %>% summarise(Bias = sum(Bias) * (1/500)) %>% pull(Bias)
}
bias_roll <- function(x) {
  x %>% pivot_longer(cols = -c(y, roll_y),names_to = "Series",values_to = "Value") %>% 
    mutate(Bias = roll_y - Value) %>% summarise(Bias = sum(Bias) * (1/500)) %>% pull(Bias)
}
sum_diff <- function(data) { 
  data %>% #filter(!grepl("Series",LM_estimate)) %>% 
  pivot_wider(names_from = LM_estimate, values_from = Value) %>% 
  mutate(Diff = roll_y - SBB_mean) %>% pull(Diff)
}


load("interim/StationaryBoot_Sampling.RData")
station_boot
roll_boots <- function(data) {
  data %>% 
    mutate(across(starts_with("Series"), \(x) ifelse(is.na(n.x),NA,x))) %>% 
    mutate(across(starts_with("Series"), \(x) zoo::rollapply(as_vector(x),width=10, FUN = mean ,na.rm=T,partial=T,align="center"), .names= "{'Roll'}_{.col}" ),.after=2) %>% 
    mutate(across(starts_with("Series"), \(x) zoo::rollapply(as_vector(x),width=10, FUN = sd ,na.rm=T,partial=T,align="center"), .names= "{'Roll_sd'}_{.col}" ),.after=2) %>% 
    select(-starts_with("Series"))
  }
station_boot %<>% rowwise() %>% mutate(Roll_Boots = list(roll_boots(Boostraps)))
save(station_boot,file="interim/StationaryBoot_Sampling.RData")

# LOOK AT BIAS AND CONFIDENCE IN BOOTSTRAPS AS A FUNCTION OF OCCURRENCE NUMBERS
# LOOK AT WHAT HAPPENS WITH SPECIES WITH 50 < OCCURRENCES > 100 ?!
load("interim/StationaryBoot_Sampling.RData")
otra <- station_boot %>% unnest(Roll_Boots) %>% 
separate(Accepted_Name,into = c("G","S","extra"),extra = "merge",sep=" ") %>% select(-extra) %>% unite("ind",G:S,sep="_") %>% 
  select(-data,-data_points,-Boostraps)
series_n <- paste0("Series.",1:1000)
est_boot <- list()
coeffs <- list()
for (i in 1:length(series_n)){
  cat("Processing",i,"\n")
una <- otra %>% select(ind,tme,n.x,ends_with(series_n[i])) %>% 
  mutate(across(starts_with("Roll_sd"), \(x) x / sqrt(n.x), .names = "{'Roll_error'}")) %>% 
  rename_with(~c("ind", "tme", "n", "Roll_sd", "Roll_Series", "Roll_error"))
cat("       Fitting model","\n")
fit1 <- lmer(Roll_Series ~ 1 + tme + Roll_error + (1|ind) + (0 + tme|ind),data = una, control=lmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e4)))
est_boot[[i]] <- broom.mixed::tidy(fit1) %>% mutate(across(4:6, \(x) round(x,3))) %>% mutate(Series = series_n[i])
coeffs[[i]] <- coef(fit1)$ind %>% as_tibble(rownames = "ind") %>% mutate(Series = series_n[i])
}
boot_list <- list(est_boot,coeffs)
save(boot_list,file = "output/boostrap_processed.RData")

load("output/boostrap_processed.RData")
est_boot <- boot_list[[1]]
coeffs <- boot_list[[2]]

mat <- lapply(est_boot, FUN = function(x){as.matrix(x[,-c(1:3,7)])})
mean_boot <- apply(abind::abind(mat, along=3), MARGIN=c(1,2), FUN = mean)
matad <- lapply(est_boot, FUN = "[", 2,) %>% data.table::rbindlist() %>% as_tibble()

load("output/LMM_models_Mid_v3.RData")
obs <- broom.mixed::tidy(modelos$cov) %>% mutate(across(4:6, \(x) round(x,3))) 
mean_boot %<>% as_tibble() %>%  mutate(term = c('(Intercept)', "tme","error","sd__(Intercept)","sd__tme","sd__Observation"))
matad %>% 
ggplot() + 
  ggdist::stat_dotsinterval(aes(y = term, x = estimate),point_interval = "mean_qi") +
  geom_point(data = obs %>% filter(term=="tme"),aes(y = 1.8, x = estimate), size=4,
               color=MetBrewer::met.brewer("Hiroshige",direction=-1,n=6)[6]) +
  geom_segment(data = obs %>% filter(term=="tme"),aes(y = 1.8,yend=term, x = estimate,xend=estimate),
               linewidth=1, color=MetBrewer::met.brewer("Hiroshige",direction=-1,n=6)[6]) +
    coord_cartesian(xlim=c(1.5,3),expand=TRUE) +
  scale_y_discrete(expand=expansion(0,.1)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(family="EB Garamond",size=14),
        axis.text = element_text(family="EB Garamond")) +
  labs(x="Estimated general elevational shift (meters per year)") +
  scale_y_discrete(expand=c(0.02,0)) +
  NULL
  

coe_mat <- lapply(coeffs, FUN = function(x){as.matrix(x[,c(2:3)])})
coe_boot <- apply(abind::abind(coe_mat, along=3), MARGIN=c(1,2), FUN = mean)
coe_boot_sd <- apply(abind::abind(coe_mat, along=3), MARGIN=c(1,2), FUN = sd)
coe_boot %<>% as_tibble() %>% mutate(ind = coeffs[[1]]$ind,.before=1) %>% 
  mutate(tme_sd = coe_boot_sd[,2], '(Intercept)_sd' = coe_boot_sd[,1])

coe_boot %>% mutate(UCI = tme + 2.58*(tme_sd/sqrt(1000))) %>% mutate(LCI = tme - 2.58*(tme_sd/sqrt(1000))) %>% 
  mutate(SIG = ifelse(UCI < 0 & LCI < 0, TRUE, ifelse(UCI > 0 & LCI > 0, TRUE, FALSE))) %>% filter(SIG)

ns <- station_boot %>% unnest(data) %>% group_by(Accepted_Name) %>% summarise(ns = sum(n, na.rm=TRUE))
range <- station_boot %>% unnest(data) %>% group_by(Accepted_Name) %>% 
  summarise(up = max(Variable,na.rm=TRUE), down = min(Variable,na.rm=TRUE)) %>% 
  mutate(range = up - down)
ns %>% separate(Accepted_Name,into = c("G","S","extra"),extra = "merge",sep=" ") %>% select(-extra) %>% unite("ind",G:S,sep="_") %>% left_join(.,coe_boot,by="ind") %>% 
  ggplot() +
  #geom_point(aes(x= range, y= tme_sd)) +
  #geom_smooth(aes(x= range, y= tme_sd),method="lm") +
  geom_point(aes(x= ns, y= tme_sd)) +
  geom_smooth(aes(x= ns, y= tme_sd),method="lm") +
  scale_x_continuous(transform="log10")


esta <- readRDS(here("output/LMM_trends_v.3.RDS"))
coe_boot %>% right_join(.,esta %>% select(ind,trend_mid),by="ind") %>% 
  mutate(Upper = tme + 2.58*(tme_sd/sqrt(1000))) %>% 
  mutate(Lower = tme - 2.58*(tme_sd/sqrt(1000))) %>% 
  mutate(Within = trend_mid < 0 & Upper < 0 & Lower < 0) %>% 
  mutate(Within = ifelse(!Within, trend_mid > 0 & Upper  > 0 & Lower > 0,Within)) %>% 
  count(Within)

ggplot() +
  #geom_point(aes(x=trend_mid,y=fct_reorder(ind,trend_mid)),size=0.5) +
  geom_segment(aes(x = Upper, xend = Lower,y= fct_reorder(ind,trend_mid),yend=fct_reorder(ind,trend_mid),
                   color=Within),alpha=0.9,linewidth=0.5) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.line.x = element_line()) +
  geom_vline(xintercept = 0) +
  NULL




una <- station_boot %>% select(LMs) %>% unnest(LMs) %>% 
  #filter(LM_estimate %in% c("y","roll_y","SBB_mean","SBB_upper","SBB_lower")) %>% 
  pivot_wider(values_from = Value, names_from = LM_estimate) %>% 
  mutate(Bias = y - SBB_mean)
  #select(Accepted_Name:roll_y,starts_with("Series")) %>% 
  #rowwise() %>% 
  #mutate(SBB_mean = mean(c_across(starts_with("Series")))) %>% 
  #mutate(SBB_sd = sd(c_across(starts_with("Series")))) %>% 
  #nest(.by = Accepted_Name) %>% 
  #rowwise() %>% 
  #mutate(Bias_y = bias(data),.after=1) %>% 
  #mutate(Bias_roll = bias_roll(data),.after=1) %>% 
  #unnest(data) %>% select(Accepted_Name:roll_y,SBB_mean,SBB_sd)

to_plot <- una %>% 
  #mutate(SBB_upper = SBB_mean + (SBB_sd * 2.32 )) %>% 
#  mutate(SBB_lower = SBB_mean - (SBB_sd * 2.32 )) %>%
#  mutate(sd = SBB_upper - SBB_mean) %>% 
#  mutate(Diff = roll_y - SBB_mean) %>% 
#  mutate(Within = abs(y - SBB_mean) < abs((SBB_sd * 2.32 ))) %>%
  mutate(SameSign = y < 0 & SBB_mean < 0 | y > 0 & SBB_mean > 0 ) %>% 
 # mutate(Significant = SBB_lower < 0 & SBB_upper < 0 | SBB_lower > 0 & SBB_upper > 0 ) %>% 
  #ungroup() %>% count(Within)
  # 921/1021 species with the same tredn direction as mean bootstraps #898/1021 for roll means
  # 931/1021 species within confidence intervals of bootstraps # 891/1021 for roll means
  arrange(y)

#### 

esta <- readRDS(here("output/LMM_trends_v.3.RDS"))
una %<>% separate(Accepted_Name,into=c("G","S"),sep=" ",extra = "drop",remove=FALSE) %>% 
  unite("ind", c(G,S),sep = "_",remove=TRUE) 
esta %>% left_join(.,una,by="ind") %>% 
  mutate(Bias_Cat = cut(Bias, breaks = c(-Inf,-15,-10,-5,0,5,10,15,Inf)),
         include.lowest=TRUE) %>% 
  ggplot() +
  geom_point(aes(y=Bias, x=trend_mid, color=error, fill=error)) +
 geom_segment(aes(x = Bias, xend = Bias,  y = Bias, yend = trend_mid, color=error))
  #facet_wrap(~Bias_Cat) +
   geom_hline(yintercept = 0, linewidth=0.5, linetype="dashed")
  geom_abline(slope = 1,intercept = 0, linewidth=0.5, linetype="dashed") +
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(family="EB Garamond",size=14,hjust = 0.1,face="bold.italic"),
        axis.text.x = element_text(family="EB Garamond",size=12,face="bold"),
        axis.title = element_text(family="EB Garamond",size=12,face="bold"),
        axis.text.y = element_text(family="EB Garamond",size=12,face="bold")) +
  NULL
  


ns <- to_plot %>% ungroup() %>% count(SameSign) %>% pull(n)
facert_label = c("TRUE" = paste0("Same Sign\n",ns[2], " species"), "FALSE" = paste0("Different sign\n",ns[1], " species"))
to_plot %>%  ggplot()  + 
  geom_point(aes(y = y, x = Bias, color= SameSign, fill = SameSign),alpha=0.9) +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap( ~ SameSign, labeller = as_labeller(facert_label)) +
  #ggdist::geom_dots(aes(x = SBB_mean),inherit.aes = FALSE,shape=21) +
  scale_color_met_d("Kandinsky") +  
  scale_fill_met_d("Kandinsky") 
  geom_hline(yintercept = prom, 
             color = rep(met.brewer(name = "Kandinsky",n=2),2), 
             linetype = c("dashed"),linewidth=0.6) +

  geom_hline(yintercept = 0,color="grey30",linewidth=0.6) +
  geom_segment(aes(x=y,xend= y, y=abs(Diff),yend=0,color=Within),linewidth=0.3,alpha=0.5) +
#geom_text(data = numbers %>% arrange(desc(SES)),aes(x=rep(10,3),y=seq(1.5,3.5,1),label=paste0(n," species")),family="EB Garamond") +
  #geom_vline(xintercept = 0) +
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(family="EB Garamond",size=14,hjust = 0.1,face="bold.italic"),
        axis.text.x = element_text(family="EB Garamond",size=12,face="bold"),
        axis.title = element_text(family="EB Garamond",size=12,face="bold"),
        axis.text.y = element_text(family="EB Garamond",size=12,face="bold")) +
  scale_y_continuous(expand = c(0, 0)) +
 # scale_x_continuous(trans="sqrt",breaks = c(0,1.28,10,20)) +
  labs(x="Observed linear trend",y="Difference of observed to mean boostrap trend") +
  guides(color = guide_legend(title="Elevation point",reverse = TRUE),
         fill = guide_legend(title="Elevation point",reverse = TRUE)) +
  NULL



 to_plot <- una %>% mutate(SBB_upper = SBB_mean + (SBB_sd * 2.32 )) %>% 
  mutate(SBB_lower = SBB_mean - (SBB_sd * 2.32 )) %>% 
  mutate(Diff = y - SBB_mean ) %>% 
  mutate(sd = SBB_upper - SBB_mean) %>% 
  mutate(Within = abs(Diff) < abs((SBB_sd * 2.32 ))) %>% 
  mutate(SameSign = y < 0 & SBB_mean < 0 | y > 0 & SBB_mean > 0 ) %>% 
  mutate(Significant = SBB_lower < 0 & SBB_upper < 0 | SBB_lower > 0 & SBB_upper > 0 ) %>% 
   # 591/1021 species with boostrap intervals not spanning zero
  arrange(y)
ns <- to_plot %>% ungroup() %>% count(Significant) %>% pull(n)
facert_label = c("TRUE" = paste0("Significant\n",ns[1], " species"), "FALSE" = paste0("Non significant\n",ns[2], " species"))
to_plot %>% ggplot()  + 
  geom_point(aes(y = SBB_sd, x = y, color= Significant, fill = Significant),alpha=0.9) +
  facet_wrap(~Significant,labeller = as_labeller(facert_label))+
  #ggdist::geom_dots(aes(x = SBB_mean),inherit.aes = FALSE,shape=21) +
  scale_color_met_d("Kandinsky") +  
  scale_fill_met_d("Kandinsky") +
  #geom_hline(yintercept = 0,color="grey30",linewidth=0.6) +
  geom_segment(aes(x = y,xend = y, y=SBB_sd ,yend=0,color=Significant),linewidth=0.3,alpha=0.5) +
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(family="EB Garamond",size=14,hjust = 0.1,face="bold.italic"),
        axis.text.x = element_text(family="EB Garamond",size=12,face="bold"),
        axis.title = element_text(family="EB Garamond",size=12,face="bold"),
        axis.text.y = element_text(family="EB Garamond",size=12,face="bold")) +
  scale_y_continuous(expand = c(0, 0)) +
  # scale_x_continuous(trans="sqrt",breaks = c(0,1.28,10,20)) +
  labs(x="Observed linear trends",y="Standard deviation bootstrap replicates") +
  guides(color = guide_legend(title="Elevation point",reverse = TRUE),
         fill = guide_legend(title="Elevation point",reverse = TRUE)) +
  NULL


to_plot %>% ggplot()  + 
  geom_point(aes(y = SBB_sd, x = Diff),alpha=0.9) +
  geom_segment(aes(x = Diff,xend = Diff, y=SBB_sd ,yend=0),linewidth=0.3,alpha=0.5) +
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(family="EB Garamond",size=14,hjust = 0.1,face="bold.italic"),
        axis.text.x = element_text(family="EB Garamond",size=12,face="bold"),
        axis.title = element_text(family="EB Garamond",size=12,face="bold"),
        axis.text.y = element_text(family="EB Garamond",size=12,face="bold")) +
  scale_y_continuous(expand = c(0, 0)) +
  # scale_x_continuous(trans="sqrt",breaks = c(0,1.28,10,20)) +
  labs(x="Difference of observed to mean boostrap trend",y="Standard deviation bootstrap replicates") +
  guides(color = guide_legend(title="Elevation point",reverse = TRUE),
         fill = guide_legend(title="Elevation point",reverse = TRUE)) +
  NULL

to_plot %>% arrange(Accepted_Name) %>% 
  #pivot_longer(cols = c(Accepted_Name:SBB_mean),names_to = "Estimate", values_to = "Values") %>% 
ggplot() +
 #geom_segment(aes(y = 1:1021, yend = Accepted_Name, x = y, xend = SBB_mean)) +
  geom_point(aes(y = 1:1021, x=SBB_mean),color="red") +
  geom_point(aes(y = 1:1021, x=y),color="blue") 
  NULL
