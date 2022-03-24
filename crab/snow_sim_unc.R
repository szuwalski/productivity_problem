#======================================================
# Inputs
#======================================================
library(reshape)
library(dplyr)
library(ggplot2)
library(ggridges)
library(gridExtra)
library(png)
library(grid)
source("functions/pop_dy.R")
source("functions/calc_SBPR.R")

.THEME    = theme_bw(base_size = 12, base_family = "") +
  theme(strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(color="white",fill="white")) 
.COL =  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

sizes<-seq(27.5,132.5,5)
size_trans<-as.matrix(read.csv("crab/size_trans.csv",header=F))
specs<-read.csv("crab/snow_specs.csv")
time_specs<-read.csv("crab/time_specs.csv")

#================================================================
# simulate for historical status
#================================================================
#==recruitment
log_avg_rec<-13.2956704191
in_rec<-exp(log_avg_rec+time_specs$rec_dev)/2

#==natural and fishing mortality
nat_m<-matrix(c(rep(0.31,nrow(time_specs)),rep(0.30,nrow(time_specs))),ncol=2)
in_f<-time_specs$fmort

orig<-pop_dy(size_trans=size_trans,
             specs=specs,
             rec=in_rec,
             fmort=in_f,
             sizes=sizes,
             nat_m=nat_m)

no_f<-pop_dy(size_trans=size_trans,
             specs=specs,
             rec=in_rec,
             fmort=rep(0,nrow(time_specs)),
             sizes=sizes,
             nat_m=nat_m)

#==reference point beginnings
SBPR<-calc_SBPR(size_trans,specs,sizes,in_m=c(0.31,0.30))
SBPRF35<-SBPR[[2]]
F35<-SBPR[[1]]

#==calculate status for proxies over time as recruitment updates
B_35<-rep(NA,nrow(time_specs))
B_35[1:10]<-SBPRF35 * mean(in_rec[1:10])
for(x in 11:nrow(time_specs))
  B_35[x]<-SBPRF35 * mean(in_rec[1:x])

#==calculate mature biomass
MMB<-matrix(ncol=nrow(time_specs),nrow=2)
for(x in 1:ncol(MMB))
{
  MMB[1,x]<-sum(orig$n_at_l_mate[x,,2]*specs$weight,na.rm=T)
  MMB[2,x]<-sum(no_f$n_at_l_mate[x,,2]*specs$weight,na.rm=T)
}  

#==plot yields
#==calculate mature biomass
Yields<-matrix(ncol=nrow(time_specs),nrow=2)
for(x in 1:ncol(MMB))
{
  Yields[1,x]<-sum(orig$c_at_l[x,,1]*specs$weight,na.rm=T) + sum(orig$c_at_l[x,,2]*specs$weight,na.rm=T)
  Yields[2,x]<-sum(no_f$c_at_l[x,,1]*specs$weight,na.rm=T) + sum(no_f$c_at_l[x,,2]*specs$weight,na.rm=T)
}  


#===============================================================
# Project under different productivity and management
#===============================================================
#===========================================
# Recruitment changes scenario
#===========================================
alpha<-0.1 # Harvest control rule parameter
set.seed(10)
proj_yr<-2099
nsim<-100
#==historical recruitment
log_avg_rec<-13.2956704191
in_rec_1<-exp(log_avg_rec+time_specs$rec_dev)/2

#========================================================
#==relationship to environmental indices
#==based off of Szuwalski et al. 2020
env_data<-read.csv("climate/Bering sea environmental variables.csv")
snow_lag<-5
snow_yrs<-seq(1982,1981+length(in_rec_1))[snow_lag:length(in_rec_1)]-snow_lag
use_snow_rec<-in_rec_1[snow_lag:length(in_rec_1)]
use_snow_SB<-MMB[1,1:(ncol(MMB)-snow_lag+1)]
use_env_snow<-env_data[which(!is.na(match(env_data$Year,snow_yrs))),]
tmp<-data.frame(logRS = log(use_snow_rec/use_snow_SB),
                AO = use_env_snow$Arctic_oscillation,
                Ice = use_env_snow$Ice_cover,
                SB = use_snow_SB)
rownames(tmp)<-snow_yrs

#==removed because initial and final years estimates are not good
tmp<-tmp[-nrow(tmp),]
tmp<-tmp[-c(1,2),]
mod<-lm(logRS~.,data=tmp)

#==projected environmental indices
#==ice cover
proj_ice_raw<-unlist(read.csv("climate/proj_ice_raw.csv"))
#==Arctic oscillation
AO_proj<-read.csv("climate/AO_proj.csv")
colnames(AO_proj)<-"AO"

#==predict future recruitment with constant mean SB
in_dat<-data.frame(SB=rep(mean(tmp$SB),length(proj_ice_raw)),
                   Ice=proj_ice_raw,
                   AO=AO_proj)
rownames(in_dat)<-seq(2013,proj_yr)

#==natural mortality
nat_m<-matrix(c(rep(0.31,length(in_rec_1)),rep(0.30,length(in_rec_1))),ncol=2)
in_f<-time_specs$fmort
proj_rec_SQ<-list(list())
proj_rec_chng<-list(list())

big_proj_sq<-list(list())
big_proj_ch<-list(list())

#==set up initial year
proj_rec_SQ[[1]]<-pop_dy(size_trans=size_trans,
                         specs=specs,
                         rec=in_rec_1,
                         fmort=in_f,
                         sizes=sizes,
                         nat_m=nat_m)

proj_rec_chng[[1]]<-pop_dy(size_trans=size_trans,
                         specs=specs,
                         rec=in_rec_1,
                         fmort=in_f,
                         sizes=sizes,
                         nat_m=nat_m)

proj_yrs<-nrow(in_dat)
b35_sq<-matrix(0,nrow=nsim,ncol=proj_yrs)
b35_chg<-matrix(0,nrow=nsim,ncol=proj_yrs)
in_f_sq<-matrix(0,nrow=nsim,ncol=proj_yrs)
for(x in 1:nrow(in_f_sq))
 in_f_sq[x,1:length(time_specs$fmort)]<-time_specs$fmort

in_f_chg<-matrix(0,nrow=nsim,ncol=proj_yrs)
for(x in 1:nrow(in_f_chg))
  in_f_chg[x,1:length(time_specs$fmort)]<-time_specs$fmort

in_rec_in_sq<-matrix(0,nrow=nsim,ncol=proj_yrs)
for(x in 1:nrow(in_rec_in_sq))
  in_rec_in_sq[x,1:length(time_specs$fmort)]<-in_rec_1

in_rec_in_chg<-matrix(0,nrow=nsim,ncol=proj_yrs)
for(x in 1:nrow(in_rec_in_chg))
  in_rec_in_chg[x,1:length(time_specs$fmort)]<-in_rec_1

for(m in 1:nsim)
{
 for(z in 2:(proj_yrs-length(time_specs$fmort)))
 {
  #==update natural mortality
  nat_m<-matrix(c(rep(0.31,length(in_rec_1)+z-1),rep(0.30,length(in_rec_1)+z-1)),ncol=2)
  #=====================================
  # status quo, doesn't change reference
  #=====================================
  #==harvest control rule
  #==calculate MMB
  temp<-proj_rec_SQ[[z-1]]$n_at_l_mate[,,2]
  MMB<-rep(NA,nrow(temp))
  for(x in 1:length(MMB))
    MMB[x]<-sum(temp[x,]*specs$weight,na.rm=T)
    
  #==calculate B35
  b35_sq[m,z]<-SBPRF35*mean(in_rec_in_sq[m,1:(length(in_rec_1)+z-2)])
  
  #==find F
  new_f<-F35
  if(MMB[length(MMB)-1]/b35_sq[m,z] < 0.25)
    new_f <- 0
  if(MMB[length(MMB)-1]/b35_sq[m,z] >= 0.25 & MMB[length(MMB)-1]/b35_sq[m,z]<1)
    new_f <- (F35*(MMB[length(MMB)-1]/b35_sq[m,z] - alpha))/(1-alpha)

  in_f_sq[m,z+length(time_specs$fmort)-1]<-new_f
  
  
  #==CALCULATE RECRUITMENT FOR YEAR BASED ON SSB 5 YEARS IN PAST + ENV INDICES
  pred_dat<-in_dat[(z-1),]
  pred_dat[1]<-MMB[length(MMB)-5]
  tmp_pred<-  predict(mod,newdata=pred_dat,se.fit=TRUE)
  pred_rec_in<-MMB[length(MMB)-5]*exp(rnorm(1,tmp_pred$fit,tmp_pred$se.fit))

  in_rec_in_sq[m,z+length(time_specs$fmort)-1] <- pred_rec_in
  
  proj_rec_SQ[[z]]<-pop_dy(size_trans=size_trans,
                        specs=specs,
                        rec=in_rec_in_sq[m,1:(z+length(time_specs$fmort)-1)],
                        fmort=in_f_sq[m,1:(z+length(time_specs$fmort)-1)],
                        sizes=sizes,
                        nat_m=nat_m)
  
  
  #=====================================
  # change reference points
  #=====================================
  #==harvest control rule
  #==calculate MMB
  temp<-proj_rec_chng[[z-1]]$n_at_l_mate[,,2]
  MMB<-rep(NA,nrow(temp))
  for(x in 1:length(MMB))
    MMB[x]<-sum(temp[x,]*specs$weight,na.rm=T)
  
  #==calculate B35
  if(z<10)
    b35_chg[m,z]<-SBPRF35*mean(in_rec_in_chg[m,1:(length(in_rec_1)+z-2)])
  if(z>=10)
    b35_chg[m,z]<-SBPRF35*mean(in_rec_in_chg[m,(length(in_rec_1)+5):(length(in_rec_1)+z-2)])
  
  #==find F
  new_f<-F35
  if(MMB[length(MMB)-1]/b35_chg[m,z] < 0.25)
    new_f <- 0
  if(MMB[length(MMB)-1]/b35_chg[m,z] >= 0.25 & MMB[length(MMB)-1]/b35_chg[m,z]<1)
    new_f <- (F35*(MMB[length(MMB)-1]/b35_chg[m,z] - alpha))/(1-alpha)
  
  in_f_chg[m,z+length(time_specs$fmort)-1]<-new_f
  
  #==predict recruitment
  pred_dat<-in_dat[(z-1),]
  pred_dat[1]<-MMB[length(MMB)-5]
  tmp_pred<-  predict(mod,newdata=pred_dat,se.fit=TRUE)
  pred_rec_in<-MMB[length(MMB)-5]*exp(rnorm(1,tmp_pred$fit,tmp_pred$se.fit))

  in_rec_in_chg[m,z+length(time_specs$fmort)-1] <- pred_rec_in
  
  
  proj_rec_chng[[z]]<-pop_dy(size_trans=size_trans,
                           specs=specs,
                           rec=in_rec_in_chg[m,1:(z+length(time_specs$fmort)-1)],
                           fmort=in_f_chg[m,1:(z+length(time_specs$fmort)-1)],
                           sizes=sizes,
                           nat_m=nat_m)
 } #years
  big_proj_sq[[m]]<-  proj_rec_SQ[[z]]
  big_proj_ch[[m]]<-  proj_rec_chng[[z]]
} # sim loop


#==plot yields
yields_sq<-matrix(ncol=proj_yrs,nrow=nsim)
yields_ch<-matrix(ncol=proj_yrs,nrow=nsim)
for(x in 1:nsim)
 for(y in 1:(proj_yrs-1))
  {
   yields_sq[x,y]<-sum(big_proj_sq[[x]]$c_at_l[y,,1]*specs$weight,na.rm=T) + sum(big_proj_sq[[x]]$c_at_l[y,,2]*specs$weight,na.rm=T)
   yields_ch[x,y]<-sum(big_proj_ch[[x]]$c_at_l[y,,1]*specs$weight,na.rm=T) + sum(big_proj_ch[[x]]$c_at_l[y,,2]*specs$weight,na.rm=T)
 }  

#==plot MMB
  mmb_sq<-matrix(ncol=proj_yrs,nrow=nsim)
  mmb_ch<-matrix(ncol=proj_yrs,nrow=nsim)
  for(x in 1:nsim)
    for(y in 1:(proj_yrs-1))
    {
      mmb_sq[x,y]<- sum(big_proj_sq[[x]]$n_at_l[y,,2]*specs$weight,na.rm=T)
      mmb_ch[x,y]<- sum(big_proj_ch[[x]]$n_at_l[y,,2]*specs$weight,na.rm=T)
    }  
  
   #==plot eMMB
  emmb_sq<-matrix(ncol=proj_yrs,nrow=nsim)
  emmb_ch<-matrix(ncol=proj_yrs,nrow=nsim)
  for(x in 1:nsim)
    for(y in 1:(proj_yrs-1))
    {
      emmb_sq[x,y]<- sum(big_proj_sq[[x]]$n_at_l[y,,2]*specs$weight*specs$fish_sel,na.rm=T)
      emmb_ch[x,y]<- sum(big_proj_ch[[x]]$n_at_l[y,,2]*specs$weight*specs$fish_sel,na.rm=T)
    }  

  #==plot F
  in_u_ch<-1-exp(-in_f_chg)
  in_u_sq<-1-exp(-in_f_sq)


#===============================================================
# Plot figure that shows catches under productivity shift
# and change in reference points vs. no change
# do one for both a change in recruitment and a change in M
# remember to change SBPR35
#===============================================================

make_dat<-function(input,input_q,input_hcr,adj=1,begin_HCR=NA)
 {
  tmp<-apply(input,2,sort)
  med_ch<-as.numeric(tmp[50,])/adj
  up_ch<-as.numeric(tmp[1,])/adj
  dn_ch<-as.numeric(tmp[100,])/adj
  in_HCR<-rep(input_hcr,length(med_ch))
  if(!is.na(begin_HCR))
    in_HCR[1:begin_HCR]<-"Historical"
  out<-data.frame(median=med_ch,upper=up_ch,lower=dn_ch,Year=seq(from=1982,length.out=length(med_ch)),
             Quantity=input_q,HCR=in_HCR)
  return(out)
 }

use_colors<-c("blue","dark grey","green")

  .THEME    = theme_bw(base_size = 12, base_family = "") +
    theme(strip.text.x = element_text(margin= margin(1,0,1,0)),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          strip.background = element_rect(color="white",fill="white")) 
  .COL =  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
    
  tmp<-make_dat(input=in_u_ch,input_q="Exp_rate",input_hcr="Climate adaptive",begin_HCR=38)
  tmp2<-make_dat(input=in_u_sq,input_q="Exp_rate",input_hcr="Status quo",begin_HCR=38)
  in_exp<-as.data.frame(rbind(tmp2,tmp))
  
  exp_proj<-ggplot() +
    geom_ribbon(data=in_exp,aes(x=Year,ymin = lower, ymax = upper, fill = HCR),alpha=0.1) +
    geom_line(data=in_exp,aes(x=Year,y=median,col=HCR)) +
    .THEME +
  scale_y_continuous(position = "right")+
    theme(legend.position = "none",
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y=element_blank()) + 
  annotate("text",x=1993,y=1,label="Harvest rate") +
    geom_vline(xintercept=2019.5,lty=2)+ 
    xlim(1982,2065)+
    scale_color_manual(values=use_colors) +
    scale_fill_manual(values=use_colors)
  
  tmp<-make_dat(input=in_rec_in_sq,input_q="Recruits",input_hcr="Status quo",adj=10000,begin_HCR=38)
  tmp2<-make_dat(input=in_rec_in_chg,input_q="Recruits",input_hcr="Climate adaptive",adj=10000,begin_HCR=38)
  in_exp<-as.data.frame(rbind(tmp,tmp2))

  in_png<-readPNG('crab/silhouette_crab.png')
  
  rec_proj<-ggplot() +
    geom_ribbon(data=in_exp,aes(x=Year,ymin = lower, ymax = upper, fill = HCR),alpha=0.1) +
    geom_line(data=in_exp,aes(x=Year,y=median,col=HCR)) +
    .THEME +
    scale_y_continuous(position = "right")+
    theme(legend.position = 'none',
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y=element_blank()) + 
  annotate("text",x=1994,y=300,label="Recruits (10000s)") +
    annotation_raster(in_png,
                      ymin=150,
                      ymax=350,
                      xmin=2035,
                      xmax=2065)+
    geom_vline(xintercept=2019.5,lty=2)+ 
    xlim(1982,2065)+
    scale_color_manual(values=use_colors) +
    scale_fill_manual(values=use_colors)

  tmp<-make_dat(input=yields_sq[,1:(ncol(yields_sq)-2)],input_q="Yields",input_hcr="Status quo",adj=1,begin_HCR=38)
  tmp2<-make_dat(input=yields_ch[,1:(ncol(yields_sq)-2)],input_q="Yields",input_hcr="Climate adaptive",adj=1,begin_HCR=38)
  in_exp<-as.data.frame(rbind(tmp,tmp2))
  
  yield_proj<-ggplot() +
    geom_ribbon(data=in_exp,aes(x=Year,ymin = lower, ymax = upper, fill = HCR),alpha=0.1) +
    geom_line(data=in_exp,aes(x=Year,y=median,col=HCR)) +
    .THEME +
    scale_y_continuous(position = "right")+
    theme(legend.position = c(.75,1.05),
          axis.title.y=element_blank(),axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x = element_blank()) + 
    annotate("text",x=1990,y=140,label="Yield (1000t)")+
    geom_vline(xintercept=2019.5,lty=2) + 
    xlim(1982,2065)+
    scale_color_manual(values=use_colors) +
    scale_fill_manual(values=use_colors)
  
  
  tmp<-make_dat(input=emmb_sq[,1:(ncol(emmb_sq)-1)],input_q="MMB",input_hcr="Status quo",begin_HCR=38)
  tmp2<-make_dat(input=emmb_ch[,1:(ncol(emmb_sq)-1)],input_q="MMB",input_hcr="Climate adaptive",begin_HCR=38)
  in_exp<-as.data.frame(rbind(tmp,tmp2))
  
  emmb_proj<-ggplot() +
    geom_ribbon(data=in_exp,aes(x=Year,ymin = lower, ymax = upper, fill = HCR),alpha=0.1) +
    geom_line(data=in_exp,aes(x=Year,y=median,col=HCR)) +
    .THEME +
    scale_y_continuous(position = "right")+
    theme(legend.position = 'none',
          axis.title.y=element_blank()) +
    annotate("text",x=1998,y=180,label="Mature male biomass")+
    scale_color_manual(values=use_colors) +
    geom_vline(xintercept=2019.5,lty=2) +
    scale_fill_manual(values=use_colors)
  


  #======================================================================
  # make figure for cold pool
  #=========================================================================
  proj_dat_85<-read.csv('climate/BT_sta_gfdl_85.csv')
  proj_dat_45<-read.csv('climate/BT_sta_gfdl_45.csv')
  station_meta<-read.csv('climate/station_meta.csv')
  BT <- as_tibble(read.csv("climate/BT.csv"))
  
  cold_pool_85<-colSums(proj_dat_85[14:ncol(proj_dat_85)]<2)
  cold_pool_45<-colSums(proj_dat_45[14:ncol(proj_dat_45)]<2)
  past_cold<-colSums(BT[14:ncol(BT)]<2)
  # plot(past_cold~seq(1970,2018),type='l',ylim=c(0,450),xlim=c(1970,2099))
  # lines(cold_pool_85~seq(2006,2100),col=2)
  # lines(cold_pool_45~seq(2006,2100),col=3)
  
  in3<-data.frame(Scenario="Historical",Cold_pool = past_cold,Year=seq(1970,2018))
  hind_u<-mean(in3[37:49,]$Cold_pool)
  hind_var<-var(in3[37:49,]$Cold_pool)
  
  in1<-data.frame(Scenario="RCP 8.5",Cold_pool = cold_pool_85,Year=seq(2006,2100))
  fut_u<-mean(in1[1:12,]$Cold_pool)
  fut_var<- var(in1[1:12,]$Cold_pool)
  
  in2<-data.frame(Scenario="RCP 4.5",Cold_pool = cold_pool_45,Year=seq(2006,2100))
  mean(in2[1:12,]$Cold_pool)
  fut_u<-mean(in1[1:12,]$Cold_pool)
  fut_var<- var(in1[1:12,]$Cold_pool)
  
  in_dat<-rbind(in1[13:nrow(in1),],in2[13:nrow(in2),],in3)
  
  use_colors<-c("dark grey","cornflowerblue","salmon")
  
  proj_cold_pool<-ggplot(data=in_dat,aes(x=Year,y=Cold_pool,color=Scenario,fill=Scenario)) +
    geom_line() +
    geom_smooth() +
    .THEME +
    theme(legend.background = element_rect(fill = "transparent"), 
          axis.title.x=element_blank()) +
    labs(y="Cold pool") + 
    scale_x_continuous(position = 'top')+
    scale_color_manual(values=use_colors) +
    scale_fill_manual(values=use_colors)

  plot1 <- readPNG('plots/coldpool_hind.png')
  
  png("plots/fig_snow_big.png",height=5.5,width=8.5,res=600,units='in')
  grid.arrange(rasterGrob(plot1),proj_cold_pool,rec_proj,exp_proj,yield_proj,emmb_proj,layout_matrix = cbind(c(2,1,1,1),c(3,4,5,6)) )
  dev.off()
  

  