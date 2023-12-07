library(dygraphs)
library(GGally)
library(viridis)
data="interim/ClimateSeries_v.2.RData"
load(here(data))
aqui %>% group_by(YearCat,Var) %>% summarise(mean=median(Values,na.rm=TRUE),sd=mad(Values,na.rm=TRUE)) %>% filter(Var == "mean.tas") %>% ggplot(aes(x=YearCat),y=mean)
  
load(here("output/Historical_data_FULL_LUL.v.1.R"))
files_list <- list.files("data/climate_data_CHELSA",pattern = "tas_",full.names = TRUE)
a <- list()
for(k in 1:length(files_list)){
  cat("Processing: ",sub(".RData","",files_list[k]) %>% str_sub(.,start=-4L),"\n")
load(files_list[k])
if(dim(a_save)[2]!=12) {
  cat("----- Skipping year","\n"); next}
colnames(a_save) <- paste0("M_",1:ncol(a_save))
a_save %>% as_tibble() %>% bind_cols(reduct_five %>% select(CellID,Altitude),.) %>% 
  distinct(CellID,.keep_all = TRUE) %>% 
mutate(Altitude = ifelse(Altitude <= 1500, "Lowland","Highland")) %>% 
group_by(Altitude) %>% 
  summarise(across(-c(1),mean,na.rm=TRUE)) %>% mutate(Year = sub(".RData","",files_list[k]) %>% str_sub(.,start=-4L)) -> a[[k]]

}
tas_monthly <- do.call(rbind,a)
tas_monthly %<>% mutate(Decade = case_when(Year < 1990 ~ "1980s",
                    Year >= 1990 & Year < 2000 ~ "1990s",
                    Year >= 2000 & Year < 2010 ~ "2000s",
                     Year >= 2010 ~ "2010s"))
tas_monthly %<>% mutate(across(starts_with("M_"),~ .x * 0.1 - 273.15))

tas_monthly %>% filter(Altitude!="Highland") %>% group_by(Year) %>% summarise(across(starts_with("M_"), mean)) %>% 
  mutate(across(starts_with("M_"), ~ .x - mean(head(.x,10)))) %>% 
  #mutate(across(starts_with("M_"), ~ .x - mean(.x))) %>% 
  rowwise() %>% 
  mutate(Close = mean(c_across(starts_with("M_")))) %>% 
  mutate(High = quantile(c_across(starts_with("M_")),probs = 1)) %>% 
  mutate(Low = quantile(c_across(starts_with("M_")),probs = 0)) %>% 
  ungroup() %>% 
  mutate(Open = lag(Close)) %>% mutate(Open = coalesce(Open, Close)) %>% 
  select(1,Open,High,Low,Close) -> data


#data %<>% mutate(Year = ymd(Year,truncated = 2L))

to_plot <- data %>% slice(-1) %>% 
  ggplot(aes(x=Year,xend=Year,color=ifelse(Close > Open,"Warmer","Cooler"))) +
  geom_hline(yintercept = 0,alpha=0.6) +
  geom_segment(aes(y=High,yend=Low),color="black") +
  geom_segment(aes(y=Open,yend=Close),lwd=5) +
  scale_color_manual(values=MetBrewer::met.brewer("Hiroshige",n=10)[c(8,1)],name=NULL) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        legend.direction = "horizontal", 
        axis.text.x = element_blank(),
        axis.text.y = element_text(family="EB Garamond",size=16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(family="EB Garamond",size=16),
        axis.line.y = element_line(),
        axis.line.x = element_blank(),
        axis.ticks = element_blank(),
        text=element_text(family="EB Garamond"))  +
  labs(x="",y="Temperature deviation from baseline mean (ºC)") +
  coord_cartesian(ylim=c(-1.5,2.1)) +
  geom_text(aes(x=1,y = -1.65),color="black",label="1980",family="EB Garamond",size=6,hjust=0,angle=90) +
  geom_text(aes(x=39,y = -0.8),color="black",label="2019",family="EB Garamond",size=6,hjust=1,angle=90) +
NULL


inset <- tibble(x=c(1,1.2),open=1:2,close=2:1,low=0.5,high=2.5) %>% 
ggplot(aes(x=x,xend=x,color=ifelse(close > open,"Warmer","Cooler"))) + 
  scale_color_manual(values=MetBrewer::met.brewer("Hiroshige",n=10)[c(1,8)],name=NULL) +
  geom_segment(aes(y=high,yend=low),color="black",lwd=1) +
  geom_segment(aes(y=open,yend=close),lwd=10) +
  #xlim(c(0.9,1.3)) +
  coord_cartesian(xlim=c(0.98,1.3)) +
  geom_text(aes(x=1.01,y = 2.09),color="black",label="Previous mean",family="EB Garamond",size=3,hjust=0) +
  geom_text(aes(x=1.01,y = 0.91),color="black",label="Mean",family="EB Garamond",size=3,hjust=0) +
  geom_text(aes(x=1.01,y = 2.5),color="black",label="Maximum",family="EB Garamond",size=3,hjust=0) +
  geom_text(aes(x=1.01,y = 0.5),color="black",label="Minimum",family="EB Garamond",size=3,hjust=0) +
  
  geom_text(aes(x=1.21,y = 2.09),color="black",label="Mean",family="EB Garamond",size=3,hjust=0,fontface="plain") +
  geom_text(aes(x=1.21,y = 0.91),color="black",label="Previous mean",family="EB Garamond",size=3,hjust=0) +
  geom_text(aes(x=1.21,y = 2.5),color="black",label="Maximum",family="EB Garamond",size=3,hjust=0) +
  geom_text(aes(x=1.21,y = 0.5),color="black",label="Minimum",family="EB Garamond",size=3,hjust=0) +
  
  geom_segment(aes(x=1.06, xend=1.06, y = 1.95, yend = 1.05),color=MetBrewer::met.brewer("Hiroshige",n=10)[c(8)],arrow = arrow(type="open",length=unit(10,"points")),lwd=1) +
  geom_segment(aes(x=1.26, xend=1.26, y = 1.05, yend = 1.95),color=MetBrewer::met.brewer("Hiroshige",n=10)[c(1)],arrow = arrow(type="open",length=unit(10,"points")),lwd=1) +
  geom_text(aes(x=1,y = 3.2),color="black",label="Cooler",family="EB Garamond",size=6,hjust=0,fontface="bold") +
  geom_text(aes(x=1.2,y = 3.2),color="black",label="Warmer",family="EB Garamond",size=6,hjust=0,fontface="bold") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_blank(),
         axis.title = element_blank(),
        axis.line =  element_blank(),
        axis.ticks = element_blank(),
        #plot.background = element_rect(fill = NA,color="black"),
        text = element_text(family="EB Garamond")) +
NULL

to_plot + 
  annotation_custom(ggplotGrob(inset),xmin=1,xmax=10,ymin=1.2,ymax=2) +
NULL




files_list <- list.files("data/climate_data_CHELSA",pattern = "pr_",full.names = TRUE)
a <- list()
for(k in 1:length(files_list)){
  cat("Processing: ",sub(".RData","",files_list[k]) %>% str_sub(.,start=-4L),"\n")
  load(files_list[k])
  if(dim(a_save)[2]!=12) {
    cat("----- Skipping year","\n"); next}
  colnames(a_save) <- paste0("M_",1:ncol(a_save))
  a_save %>% as_tibble() %>% bind_cols(reduct_five %>% select(CellID,Altitude),.) %>% 
    distinct(CellID,.keep_all = TRUE) %>% 
    mutate(Altitude = ifelse(Altitude <= 1500, "Lowland",ifelse(Altitude > 1500 & Altitude <= 2500,"Highland","Alpine"))) %>% 
    group_by(Altitude) %>% summarise(across(-c(1),mean,na.rm=TRUE)) %>% mutate(Year = sub(".RData","",files_list[k]) %>% str_sub(.,start=-4L)) -> a[[k]]
  
}
pr_monthly <- do.call(rbind,a)
pr_monthly %<>% mutate(Decade = case_when(Year < 1990 ~ "1980s",
                                           Year >= 1990 & Year < 2000 ~ "1990s",
                                           Year >= 2000 & Year < 2010 ~ "2000s",
                                           Year >= 2010 ~ "2010s"))
pr_monthly %<>% mutate(across(starts_with("M_"),~ .x * 0.1))

pr_monthly %>% filter(Altitude=="Lowland") %>% 
  group_by(Year) %>% summarise(across(starts_with("M_"), mean)) %>% mutate(across(starts_with("M_"), ~ .x - mean(head(.x,10)))) %>% 
  rowwise() %>% 
  mutate(Close = mean(c_across(starts_with("M_")))) %>% 
  mutate(High = quantile(c_across(starts_with("M_")),probs = 1)) %>% 
  mutate(Low = quantile(c_across(starts_with("M_")),probs = 0)) %>% 
  ungroup() %>% 
  mutate(Open = lag(Close)) %>% mutate(Open = coalesce(Open, Close)) %>% 
  select(1,Open,High,Low,Close) -> data


#data %<>% mutate(Year = ymd(Year,truncated = 2L))

to_plot <- data %>% slice(-1) %>% 
  ggplot(aes(x=Year,xend=Year,color=ifelse(Close > Open,"Wetter","Dryer"))) +
  geom_hline(yintercept = 0,alpha=0.6) +
  geom_segment(aes(y=High,yend=Low),color="black") +
  geom_segment(aes(y=Open,yend=Close),lwd=5) +
  scale_color_manual(values=MetBrewer::met.brewer("Isfahan1",n=10)[c(3,7)],name=NULL) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        legend.direction = "horizontal", 
        axis.text.x = element_blank(),
        axis.text.y = element_text(family="EB Garamond",size=16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(family="EB Garamond",size=16),
        axis.line.y = element_line(),
        axis.line.x = element_blank(),
        axis.ticks = element_blank(),
        text=element_text(family="EB Garamond"))  +
  labs(x="",y="Precipitation deviation from mean (mm)") +
  coord_cartesian(ylim=c(-1100,1600)) +
  geom_text(aes(x=1,y = -1000),color="black",label="1980",family="EB Garamond",size=6,hjust=0,angle=90) +
  geom_text(aes(x=39,y = -1000),color="black",label="2019",family="EB Garamond",size=6,hjust=1,angle=90) +
  NULL


inset <- tibble(x=c(1,1.2),open=1:2,close=2:1,low=0.5,high=2.5) %>% 
  ggplot(aes(x=x,xend=x,color=ifelse(close > open,"Warmer","Cooler"))) + 
  scale_color_manual(values=MetBrewer::met.brewer("Isfahan1",n=10)[c(7,3)],name=NULL) +
  geom_segment(aes(y=high,yend=low),color="black",lwd=1) +
  geom_segment(aes(y=open,yend=close),lwd=10) +
  #xlim(c(0.9,1.3)) +
  coord_cartesian(xlim=c(0.98,1.3),ylim=c(0.5,3.4)) +
  geom_text(aes(x=1.01,y = 2.09),color="black",label="Previous mean",family="EB Garamond",size=3,hjust=0) +
  geom_text(aes(x=1.01,y = 0.91),color="black",label="Mean",family="EB Garamond",size=3,hjust=0) +
  geom_text(aes(x=1.01,y = 2.5),color="black",label="Maximum",family="EB Garamond",size=3,hjust=0) +
  geom_text(aes(x=1.01,y = 0.5),color="black",label="Minimum",family="EB Garamond",size=3,hjust=0) +
  
  geom_text(aes(x=1.21,y = 2.09),color="black",label="Mean",family="EB Garamond",size=3,hjust=0,fontface="plain") +
  geom_text(aes(x=1.21,y = 0.91),color="black",label="Previous mean",family="EB Garamond",size=3,hjust=0) +
  geom_text(aes(x=1.21,y = 2.5),color="black",label="Maximum",family="EB Garamond",size=3,hjust=0) +
  geom_text(aes(x=1.21,y = 0.5),color="black",label="Minimum",family="EB Garamond",size=3,hjust=0) +
  
  geom_segment(aes(x=1.06, xend=1.06, y = 1.95, yend = 1.05),color=MetBrewer::met.brewer("Isfahan1",n=10)[c(3)],arrow = arrow(type="open",length=unit(10,"points")),lwd=1) +
  geom_segment(aes(x=1.26, xend=1.26, y = 1.05, yend = 1.95),color=MetBrewer::met.brewer("Isfahan1",n=10)[c(7)],arrow = arrow(type="open",length=unit(10,"points")),lwd=1) +
  geom_text(aes(x=1,y = 3.2),color="black",label="Dryer",family="EB Garamond",size=6,hjust=0,fontface="bold") +
  geom_text(aes(x=1.2,y = 3.2),color="black",label="Wetter",family="EB Garamond",size=6,hjust=0,fontface="bold") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line =  element_blank(),
        axis.ticks = element_blank(),
        #plot.background = element_rect(fill = NA,color="black"),
        text = element_text(family="EB Garamond")) +
  NULL


to_plot + 
  annotation_custom(ggplotGrob(inset),xmin=1,xmax=10,ymin=1000,ymax=1500) +
  NULL


files_list <- list.files("data/climate_data_CHELSA",pattern = "tas_",full.names = TRUE)
a <- list()
for(k in 1:length(files_list)){
  cat("Processing: ",sub(".RData","",files_list[k]) %>% str_sub(.,start=-4L),"\n")
  load(files_list[k])
  if(dim(a_save)[2]!=12) {
    cat("----- Skipping year","\n"); next}
  colnames(a_save) <- paste0("M_",1:ncol(a_save))
  a_save %>% as_tibble() %>% bind_cols(reduct_five %>% select(CellID,Altitude,Accepted_Name),.) %>% 
    distinct(CellID,Accepted_Name,.keep_all = TRUE) %>% 
    group_by(Accepted_Name) %>% summarise(across(-c(1),mean,na.rm=TRUE)) %>% mutate(Year = sub(".RData","",files_list[k]) %>% str_sub(.,start=-4L)) -> a[[k]]
  
}

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
  dd %>% filter(type=="resamp") %>% select(-trend,-Inter_Cat) %>% setNames(paste0('RSamp_', names(.))) %>% 
    rename(ind=RSamp_ind) %>% 
    left_join(.,dd %>% filter(type=="down") %>% setNames(paste0('DSamp_', names(.))) %>% 
                rename(ind=DSamp_ind),by="ind") %>% 
    mutate(RSamp_SES_cat = cut(RSamp_SES,breaks=c(-6,-1.96,-1.28,0,1.28,1.96,6),labels=c("Strong_Neg","Moderate_Neg","Weak_Neg","Weak_Pos","Moderate_Pos","Strong_Pos"),include.lowest=TRUE)) %>% 
    mutate(DSamp_SES_cat = cut(DSamp_SES,breaks=c(-6,-1.96,-1.28,0,1.28,1.96,6),labels=c("Strong_Neg","Moderate_Neg","Weak_Neg","Weak_Pos","Moderate_Pos","Strong_Pos"),include.lowest=TRUE)) %>% mutate(CAT = ifelse(RSamp_estimate <= 0,"Negative","Positive")) -> nulas
  
  esta <- model_data[[1]] %>% left_join(.,model_data[[2]],by="ind")
esta %<>% left_join(.,nulas,by="ind")

monthly <- do.call(rbind,a)
monthly %<>% mutate(across(starts_with("M_"),~ .x * 0.1))

monthly %<>% separate(Accepted_Name,into = c("G","S","extra"),extra = "merge",sep=" ") %>% select(-extra) %>% unite("ind",G:S,sep="_") %>% 
  left_join(esta,.,by="ind")
monthly %>% select(ind, Year, RSamp_SES_cat, starts_with("M_")) %>% 
  group_by(ind) %>% arrange(Year) %>%  mutate(across(starts_with("M_"), ~ .x - mean(head(.x,10)))) %>%
  rowwise() %>% 
  mutate(Close = mean(c_across(starts_with("M_")))) %>% 
  mutate(High = quantile(c_across(starts_with("M_")),probs = 1)) %>% 
  mutate(Low = quantile(c_across(starts_with("M_")),probs = 0)) %>% 
  ungroup() %>% 
  mutate(Open = lag(Close)) %>% mutate(Open = coalesce(Open, Close)) %>% 
  select(1,Year,RSamp_SES_cat,Open,High,Low,Close) -> data

data %>% group_by(Year,RSamp_SES_cat) %>% summarise(across(Open:Close,mean)) %>% 
  ggplot(aes(x=Year,xend=Year,group=RSamp_SES_cat,color=RSamp_SES_cat)) +
  geom_hline(yintercept = 0,alpha=0.6) +
  geom_line(aes(y=Close)) +
  #geom_segment(aes(y=High,yend=Low),color="black") +
  #geom_segment(aes(y=Open,yend=Close),lwd=3) +
  scale_color_manual(values=MetBrewer::met.brewer("Benedictus",n=6),name=NULL)  +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal", 
        axis.text.x = element_blank(),
        axis.text.y = element_text(family="EB Garamond",size=16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(family="EB Garamond",size=16),
        axis.line.y = element_line(),
        axis.line.x = element_blank(),
        axis.ticks = element_blank(), q
        text=element_text(family="EB Garamond"))  +
  #labs(x="",y="Climate mositure index deviation from baseline") +
  #facet_wrap(~RSamp_SES_cat) + 
  #coord_cartesian(ylim=c(-1.5,2.1)) +
  #geom_text(aes(x=1,y = -1.65),color="black",label="1980",family="EB Garamond",size=6,hjust=0,angle=90) +
  #geom_text(aes(x=39,y = -0.8),color="black",label="2019",family="EB Garamond",size=6,hjust=1,angle=90) +
  NULL
