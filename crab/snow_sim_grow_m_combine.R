#======================================================================================
# Produce a figure that has both the output of changes in growth and natural mortality
#======================================================================================
sizes<-seq(27.5,132.5,5)
size_trans<-as.matrix(read.csv("crab/size_trans.csv",header=F))
specs<-read.csv("crab/snow_specs_m.csv")
time_specs<-read.csv("crab/time_specs.csv")
source("functions/pop_dy.R")
source("functions/calc_SBPR.R")

library(reshape)
library(dplyr)
library(ggplot2)
library(ggridges)

.THEME    = theme_bw(base_size = 12, base_family = "") +
  theme(strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(color="white",fill="white")) 
.COL =  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))


make_size_trans<-function(alpha,beta,growth_sd,out_plot=FALSE,no_x=FALSE)
{
  avg_post_molt<-alpha + beta*sizes
  size_trans_mat<-matrix(ncol=length(sizes),nrow=length(sizes))
  
  for(x in 1:nrow(size_trans_mat))
  {
    tmp<-dnorm(sizes,avg_post_molt[x],growth_sd)  
    tmp[seq(1,nrow(size_trans_mat))<x]<-0
    size_trans_mat[x,]<-round(tmp/sum(tmp) ,3)
    
  }

    rownames(size_trans_mat)<-sizes
    colnames(size_trans_mat)<-sizes
    in_g<-data.frame(melt(t(size_trans_mat)))
    colnames(in_g)<-c("Postmolt","Premolt","Density")
    
    p <- ggplot(in_g)
    p <- p + geom_density_ridges(aes(x=Postmolt, y=Premolt, height = Density, group=Premolt,
                                     fill=stat(y),alpha=.9999),stat = "identity",scale=3) +
      scale_fill_viridis_c()+
      .THEME +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 90)) 
    
    if(no_x==TRUE)
      p<-p+theme(axis.text.x=element_blank(),
                  axis.title.x=element_blank(),
                  axis.ticks.x = element_blank())
    
    if(out_plot==TRUE)
      print(p)
    
  list(size_trans_mat,p)
}

size_trans_1<-make_size_trans(alpha=4,beta=1.25,growth_sd=5,out_plot=F,no_x=FALSE)
size_trans_2<-make_size_trans(alpha=4,beta=1.1,growth_sd=5,out_plot=F)

png("plots/change_size.png",height=10,width=5,res=400,units='in')
grid.arrange(size_trans_1[[2]],size_trans_2[[2]])
dev.off()
#====================================
# plot comparison of HCRs for growth
#===================================
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

#==calculate mature biomass
MMB<-matrix(ncol=nrow(time_specs),nrow=2)
for(x in 1:ncol(MMB))
{
  MMB[1,x]<-sum(orig$n_at_l_mate[x,,2]*specs$weight,na.rm=T)
  MMB[2,x]<-sum(no_f$n_at_l_mate[x,,2]*specs$weight,na.rm=T)
}  

#==reference point beginnings
SBPR<-calc_SBPR(size_trans_1[[1]],specs,sizes,in_m=c(0.31,0.30))
SBPRF35<-SBPR[[2]]
F35<-SBPR[[1]]

#==reference point beginnings
SBPR_2<-calc_SBPR(size_trans_2[[1]],specs,sizes,in_m=c(0.31,0.30))
SBPRF35_2<-SBPR_2[[2]]
F35_2<-SBPR_2[[1]]

tmp_rec<-exp(log_avg_rec)
tmp_b35_1_grow<-SBPRF35*tmp_rec
tmp_b35_2_grow<-SBPRF35_2*tmp_rec
F35_1_grow<-F35
F35_2_grow<-F35_2


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

snow_lag<-5
snow_yrs<-seq(1982,1981+length(in_rec_1))[snow_lag:length(in_rec_1)]-snow_lag
use_snow_rec<-in_rec_1[snow_lag:length(in_rec_1)]
use_snow_SB<-MMB[1,1:(ncol(MMB)-snow_lag+1)]
tmp<-data.frame(logRS = log(use_snow_rec/use_snow_SB),
                SB = use_snow_SB)
rownames(tmp)<-snow_yrs

#==removed because initial and final years estimates are not good
tmp<-tmp[-nrow(tmp),]
tmp<-tmp[-c(1,2),]
mod<-lm(logRS~.,data=tmp)

in_dat<-data.frame(SB=seq(0,max(use_snow_SB),length.out=100))

proj_snow<-in_dat*exp(predict(mod,newdata=in_dat,interval="prediction"))
in_f<-time_specs$fmort
proj_rec_SQ<-list(list())
proj_rec_chng<-list(list())

big_proj_sq<-list(list())
big_proj_ch<-list(list())

#==set up initial year
proj_rec_SQ[[1]]<-pop_dy(size_trans=size_trans_1[[1]],
                         specs=specs,
                         rec=in_rec_1,
                         fmort=in_f,
                         sizes=sizes,
                         nat_m=nat_m)

proj_rec_chng[[1]]<-pop_dy(size_trans=size_trans_1[[1]],
                           specs=specs,
                           rec=in_rec_1,
                           fmort=in_f,
                           sizes=sizes,
                           nat_m=nat_m)

proj_yrs<-100
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
    nat_m<-matrix(c(rep(0.31,length(in_rec_1)+z),rep(0.30,length(in_rec_1)+z)),ncol=2)
    # if(z>6)
    #  {
    #   nat_m[(length(time_specs$fmort)+6):(length(time_specs$fmort)+z),1]<-change_m[1]
    #   nat_m[(length(time_specs$fmort)+6):(length(time_specs$fmort)+z),2]<-change_m[2]
    #   }
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
    pred_dat<-data.frame(SB=MMB[length(MMB)-5])
    tmp_pred<-  predict(mod,newdata=pred_dat,se.fit=TRUE,interval="prediction")
    #==use prediction inteval to generate recruitment
    imp_sd<-((tmp_pred$fit[1]-tmp_pred$fit[2] )/1.96)
    pred_rec_in<-MMB[length(MMB)-5]*exp(rnorm(1,tmp_pred$fit[1],imp_sd))
    #pred_rec_in<-MMB[length(MMB)-5]*exp(rnorm(1,tmp_pred$fit,tmp_pred$se.fit))
    
    in_rec_in_sq[m,z+length(time_specs$fmort)-1] <- pred_rec_in
    
    proj_rec_SQ[[z]]<-pop_dy(size_trans=size_trans_1[[1]],
                             specs=specs,
                             rec=in_rec_in_sq[m,1:(z+length(time_specs$fmort)-1)],
                             fmort=in_f_sq[m,1:(z+length(time_specs$fmort)-1)],
                             sizes=sizes,
                             nat_m=nat_m,
                             size_trans_2=size_trans_2[[1]],
                             switch_growth=length(time_specs$fmort)+6)
    
    
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
      b35_chg[m,z]<-SBPRF35_2*mean(in_rec_in_chg[m,(length(in_rec_1)+5):(length(in_rec_1)+z-2)])
    
    #==find F
    new_f<-F35
    #==after F change
    if(z>=10)
      new_f<-F35_2
    
    if(MMB[length(MMB)-1]/b35_chg[m,z] < 0.25)
      new_f <- 0
    if(MMB[length(MMB)-1]/b35_chg[m,z] >= 0.25 & MMB[length(MMB)-1]/b35_chg[m,z]<1)
      new_f <- (F35*(MMB[length(MMB)-1]/b35_chg[m,z] - alpha))/(1-alpha)
    
    in_f_chg[m,z+length(time_specs$fmort)-1]<-new_f
    
    #==predict recruitment
    pred_dat<-data.frame(SB=MMB[length(MMB)-5])
    tmp_pred<-  predict(mod,newdata=pred_dat,se.fit=TRUE,interval="prediction")
    #==use prediction inteval to generate recruitment
    imp_sd<-((tmp_pred$fit[1]-tmp_pred$fit[2] )/1.96)
    pred_rec_in<-MMB[length(MMB)-5]*exp(rnorm(1,tmp_pred$fit[1],imp_sd)) 
    in_rec_in_chg[m,z+length(time_specs$fmort)-1] <- pred_rec_in
    
    
    proj_rec_chng[[z]]<-pop_dy(size_trans=size_trans_1[[1]],
                               specs=specs,
                               rec=in_rec_in_chg[m,1:(z+length(time_specs$fmort)-1)],
                               fmort=in_f_chg[m,1:(z+length(time_specs$fmort)-1)],
                               sizes=sizes,
                               nat_m=nat_m,
                               size_trans_2=size_trans_2[[1]],
                               switch_growth=length(time_specs$fmort)+6)
  } #years
  big_proj_sq[[m]]<-  proj_rec_SQ[[z]]
  big_proj_ch[[m]]<-  proj_rec_chng[[z]]
} # sim loop


in_yr<-seq(1982,length.out=nrow(nat_m))

#==plot catches and biomass over time
plot_ind<-1
#==yields
yields_sq<-matrix(ncol=proj_yrs,nrow=nsim)
yields_ch<-matrix(ncol=proj_yrs,nrow=nsim)
for(x in 1:nsim)
  for(y in 1:(proj_yrs-1))
  {
    yields_sq[x,y]<-sum(big_proj_sq[[x]]$c_at_l[y,,1]*specs$weight,na.rm=T) + sum(big_proj_sq[[x]]$c_at_l[y,,2]*specs$weight,na.rm=T)
    yields_ch[x,y]<-sum(big_proj_ch[[x]]$c_at_l[y,,1]*specs$weight,na.rm=T) + sum(big_proj_ch[[x]]$c_at_l[y,,2]*specs$weight,na.rm=T)
  }  


#== MMB
mmb_sq<-matrix(ncol=proj_yrs,nrow=nsim)
mmb_ch<-matrix(ncol=proj_yrs,nrow=nsim)
for(x in 1:nsim)
  for(y in 1:(proj_yrs-1))
  {
    mmb_sq[x,y]<- sum(big_proj_sq[[x]]$n_at_l[y,,2]*specs$weight,na.rm=T)
    mmb_ch[x,y]<- sum(big_proj_ch[[x]]$n_at_l[y,,2]*specs$weight,na.rm=T)
  }  

#== eMMB
emmb_sq<-matrix(ncol=proj_yrs,nrow=nsim)
emmb_ch<-matrix(ncol=proj_yrs,nrow=nsim)
for(x in 1:nsim)
  for(y in 1:(proj_yrs-1))
  {
    emmb_sq[x,y]<- sum(big_proj_sq[[x]]$n_at_l[y,,2]*specs$weight*specs$fish_sel,na.rm=T)
    emmb_ch[x,y]<- sum(big_proj_ch[[x]]$n_at_l[y,,2]*specs$weight*specs$fish_sel,na.rm=T)
  }  


in_u_ch<-1-exp(-in_f_chg)
in_u_sq<-1-exp(-in_f_sq)

#===============================================================
# Plot figure 
#===============================================================
library(ggplot2)
make_dat<-function(input,input_q,input_hcr,adj=1,begin_HCR=NA)
{
  tmp<-apply(input,2,sort)
  med_ch<-as.numeric(tmp[50,])/adj
  up_ch<-as.numeric(tmp[10,])/adj
  dn_ch<-as.numeric(tmp[90,])/adj
  in_HCR<-rep(input_hcr,length(med_ch))
  if(!is.na(begin_HCR))
    in_HCR[1:begin_HCR]<-"Historical"
  out<-data.frame(median=med_ch,upper=up_ch,lower=dn_ch,Year=seq(from=1982,length.out=length(med_ch)),
                  Quantity=input_q,HCR=in_HCR)
  return(out)
}

use_colors<-c("dark grey","blue","green")

.THEME    = theme_bw(base_size = 12, base_family = "") +
  theme(strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(color="white",fill="white")) 
.COL =  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

tmp<-make_dat(input=in_u_ch,input_q="Exp_rate",input_hcr="Reference_change",begin_HCR=38)
tmp2<-make_dat(input=in_u_sq,input_q="Exp_rate",input_hcr="Status_quo",begin_HCR=38)
in_exp<-as.data.frame(rbind(tmp2,tmp))

exp_proj_grow<-ggplot() +
  geom_ribbon(data=in_exp,aes(x=Year,ymin = lower, ymax = upper, fill = HCR),alpha=0.1) +
  geom_line(data=in_exp,aes(x=Year,y=median,col=HCR)) +
  .THEME +
  scale_y_continuous(position = "left")+
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x = element_blank()) + 
  #annotate("text",x=1993,y=1,label="Exploitation rate") +
  geom_vline(xintercept=2019.5,lty=2)+ 
  xlim(1982,2065)+
  scale_color_manual(values=use_colors) +
  scale_fill_manual(values=use_colors) + 
  labs(y="Exploitation rate")

tmp<-make_dat(input=in_rec_in_sq,input_q="Recruits",input_hcr="Status_quo",adj=10000,begin_HCR=38)
tmp2<-make_dat(input=in_rec_in_chg,input_q="Recruits",input_hcr="Reference_change",adj=10000,begin_HCR=38)
in_exp<-as.data.frame(rbind(tmp,tmp2))

library(png)
in_png<-readPNG('crab/silhouette_crab.png')

rec_proj_grow<-ggplot() +
  geom_ribbon(data=in_exp,aes(x=Year,ymin = lower, ymax = upper, fill = HCR),alpha=0.1) +
  geom_line(data=in_exp,aes(x=Year,y=median,col=HCR)) +
  .THEME +
  scale_y_continuous(position = "left")+
  theme(legend.position = 'none',
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x = element_blank()) + 
 # annotate("text",x=1994,y=300,label="Recruits (10000s)") +
  geom_vline(xintercept=2019.5,lty=2)+ 
  xlim(1982,2065)+
  scale_color_manual(values=use_colors) +
  scale_fill_manual(values=use_colors) +
  labs(title = "Growth") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(y="Recruits (10000s)")
#print(rec_proj)

tmp<-make_dat(input=yields_sq[,1:(ncol(yields_sq)-2)],input_q="Yields",input_hcr="Status_quo",adj=1,begin_HCR=38)
tmp2<-make_dat(input=yields_ch[,1:(ncol(yields_sq)-2)],input_q="Yields",input_hcr="Reference_change",adj=1,begin_HCR=38)
in_exp<-as.data.frame(rbind(tmp,tmp2))

yield_proj_grow<-ggplot() +
  geom_ribbon(data=in_exp,aes(x=Year,ymin = lower, ymax = upper, fill = HCR),alpha=0.1) +
  geom_line(data=in_exp,aes(x=Year,y=median,col=HCR)) +
  .THEME +
  scale_y_continuous(position = "left")+
  theme(legend.position = c(.75,.75)) + 
  #annotate("text",x=1990,y=340,label="Yield (1000t)")+
  geom_vline(xintercept=2019.5,lty=2) + 
  xlim(1982,2065)+
  scale_color_manual(values=use_colors) +
  scale_fill_manual(values=use_colors) + 
  labs(y="Yield (1000t)")


tmp<-make_dat(input=emmb_sq[,1:(ncol(emmb_sq)-1)],input_q="MMB",input_hcr="Status_quo",begin_HCR=38)
tmp2<-make_dat(input=emmb_ch[,1:(ncol(emmb_sq)-1)],input_q="MMB",input_hcr="Reference_change",begin_HCR=38)
in_exp<-as.data.frame(rbind(tmp,tmp2))

emmb_proj_grow<-ggplot() +
  geom_ribbon(data=in_exp,aes(x=Year,ymin = lower, ymax = upper, fill = HCR),alpha=0.1) +
  geom_line(data=in_exp,aes(x=Year,y=median,col=HCR)) +
  .THEME +
  theme(legend.position = 'none',
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(position = "left")+
  theme(legend.position = 'none') + 
  geom_vline(xintercept=2019.5,lty=2) + 
  xlim(1982,2065)+
  scale_color_manual(values=use_colors) +
  scale_fill_manual(values=use_colors)+ 
  labs(y="Mature biomass (1000t)")

library(gridExtra)
# grid.arrange(rec_proj_grow,exp_proj_grow,yield_proj_grow)


#=========================================
# Close up
#=======================================
tmp<-make_dat(input=in_rec_in_sq,input_q="Recruits",input_hcr="Status_quo",adj=10000,begin_HCR=38)
tmp2<-make_dat(input=in_rec_in_chg,input_q="Recruits",input_hcr="Reference_change",adj=10000,begin_HCR=38)
in_exp<-as.data.frame(rbind(tmp,tmp2))
in_exp<-filter(in_exp,HCR!="Historical")
library(png)

rec_proj_grow_sm<-ggplot() +
  geom_ribbon(data=in_exp,aes(x=Year,ymin = lower, ymax = upper, fill = HCR),alpha=0.1) +
  geom_line(data=in_exp,aes(x=Year,y=median,col=HCR)) +
  .THEME +
  scale_y_continuous(position = "right")+
  theme(legend.position = 'none',
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y=element_blank())+ 
  xlim(2025,2065)+
  ylim(0,75) +
  scale_color_manual(values=use_colors[2:3]) +
  scale_fill_manual(values=use_colors[2:3])+ 
  annotate("text",x=2032,y=65,label="Recruits (10000s)")
#print(rec_proj)

tmp<-make_dat(input=yields_sq[,1:(ncol(yields_sq)-2)],input_q="Yields",input_hcr="Status_quo",adj=1,begin_HCR=38)
tmp2<-make_dat(input=yields_ch[,1:(ncol(yields_sq)-2)],input_q="Yields",input_hcr="Reference_change",adj=1,begin_HCR=38)
in_exp<-as.data.frame(rbind(tmp,tmp2))
in_exp<-filter(in_exp,HCR!="Historical")
yield_proj_grow_sm<-ggplot() +
  geom_ribbon(data=in_exp,aes(x=Year,ymin = lower, ymax = upper, fill = HCR),alpha=0.1) +
  geom_line(data=in_exp,aes(x=Year,y=median,col=HCR)) +
  .THEME +
  scale_y_continuous(position = "right")+
  theme(legend.position = c(.75,.8),
        axis.title.y=element_blank()) + 
  annotate("text",x=2030,y=90,label="Yield (1000t)")+
  geom_vline(xintercept=2019.5,lty=2) + 
  xlim(2025,2065) +
  ylim(0,90) +
  scale_color_manual(values=use_colors[2:3]) +
  scale_fill_manual(values=use_colors[2:3])




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

#==reference point beginnings
SBPR_2<-calc_SBPR(size_trans,specs,sizes,in_m=c(0.41,0.4))
SBPRF35_2<-SBPR_2[[2]]
F35_2<-SBPR_2[[1]]


#==========================================
# recruitment change for HCR
#==========================================
#==recruitment
log_avg_rec<-13.2956704191 # from assessment
log_avg_rec_2<-11.89279 # from projections
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
F35_rec<-SBPR[[1]]
B35_highrec<-SBPRF35*exp(log_avg_rec)
B35_lowrec<-SBPRF35*exp(log_avg_rec_2)

#====================================
# plot comparison of HCRs
#===================================
png("plots/comp_hcr.png",units='in',res=500,height=3.5,width=8)
xmax<-1.5*max(tmp_b35_1_grow,tmp_b35_2_grow)
par(mfrow=c(1,3),mar=c(.1,.1,.1,.1),oma=c(4,4,1,4))
plot(-10,xlim=c(0,xmax),ylim=c(0,0.7),las=1,ylab='',xlab='')
lines(in_y_1_grow~in_x_1_grow,lty=4,lwd=2)
lines(in_y_2_grow~in_x_2_grow,lty=2,lwd=2,col=2)
mtext(side=2,expression(paste("Fishing mortality  "," yr"^"-1")),line=2.2)
legend("bottomright",bty='n',lty=c(2,4),legend=c("Decrease growth","Status quo"),lwd=2,col=c(2,1))


tmp_rec<-exp(log_avg_rec)
tmp_b35<-SBPRF35*tmp_rec
tmp_b35_2<-SBPRF35_2*tmp_rec
plot(-10,xlim=c(0,max(2*tmp_b35,tmp_b35_2)),ylim=c(0,0.7),las=1,ylab='',xlab='',yaxt='n')
in_x<-c(.1,seq(.1*tmp_b35,tmp_b35,length.out=100),seq(tmp_b35,tmp_b35*2,length.out=100))
in_y<-c(0,seq(0,F35,length.out=100),rep(F35,100))
lines(in_y~in_x,lty=4,lwd=2)

in_x_2<-c(.1,seq(.1*tmp_b35_2,tmp_b35_2,length.out=100),seq(tmp_b35_2,tmp_b35*2,length.out=100))
in_y_2<-c(0,seq(0,F35_2,length.out=100),rep(F35_2,100))
lines(in_y_2~in_x_2,lty=2,lwd=2,col=2)
legend("bottomright",bty='n',lty=c(2,4),legend=c("Increase natural mortality","Status quo"),lwd=2,col=c(2,1))


F35_rec<-SBPR[[1]]
B35_highrec<-SBPRF35*exp(log_avg_rec)
B35_lowrec<-SBPRF35*exp(log_avg_rec_2)
plot(-10,xlim=c(0,400),ylim=c(0,0.7),las=1,ylab='',xlab='',yaxt='n')
in_x<-c(.1,seq(.1*B35_highrec,B35_highrec,length.out=100),seq(B35_highrec,B35_highrec*2,length.out=100))
in_y<-c(0,seq(0,F35_rec,length.out=100),rep(F35_rec,100))
lines(in_y~in_x,lty=4,lwd=2)

in_x_2<-c(.1,seq(.1*B35_lowrec,B35_lowrec,length.out=100),seq(B35_lowrec,B35_highrec*2,length.out=100))
in_y_2<-c(0,seq(0,F35_rec,length.out=100),rep(F35_rec,100))
lines(in_y_2~in_x_2,lty=2,lwd=2,col=2)
legend("bottomright",bty='n',lty=c(2,4),legend=c("Decrease recruitment","Status quo"),lwd=2,col=c(2,1))
mtext(side=1,"Biomass",line=2.2,outer=TRUE)
dev.off()



#==calculate mature biomass
MMB<-matrix(ncol=nrow(time_specs),nrow=2)
for(x in 1:ncol(MMB))
{
  MMB[1,x]<-sum(orig$n_at_l_mate[x,,2]*specs$weight,na.rm=T)
  MMB[2,x]<-sum(no_f$n_at_l_mate[x,,2]*specs$weight,na.rm=T)
}  

#===========================================
# Recruitment changes scenario
#===========================================
alpha<-0.1 # Harvest control rule parameter
set.seed(1)
proj_yr<-2099
nsim<-100

#==historical recruitment
log_avg_rec<-13.2956704191
in_rec_1<-exp(log_avg_rec+time_specs$rec_dev)/2

snow_lag<-5
snow_yrs<-seq(1982,1981+length(in_rec_1))[snow_lag:length(in_rec_1)]-snow_lag
use_snow_rec<-in_rec_1[snow_lag:length(in_rec_1)]
use_snow_SB<-MMB[1,1:(ncol(MMB)-snow_lag+1)]
tmp<-data.frame(logRS = log(use_snow_rec/use_snow_SB),
                SB = use_snow_SB)
rownames(tmp)<-snow_yrs

#==removed because initial and final years estimates are not good
tmp<-tmp[-nrow(tmp),]
tmp<-tmp[-c(1,2),]
mod<-lm(logRS~.,data=tmp)

in_dat<-data.frame(SB=seq(0,max(use_snow_SB),length.out=100))

proj_snow<-in_dat*exp(predict(mod,newdata=in_dat,interval="prediction"))

#==initial natural mortality
nat_m<-matrix(c(rep(0.31,length(in_rec_1)),rep(0.30,length(in_rec_1))),ncol=2)
change_m<-c(0.41,0.40)
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

proj_yrs<-100
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
    nat_m<-matrix(c(rep(0.31,length(in_rec_1)+z),rep(0.30,length(in_rec_1)+z)),ncol=2)
    if(z>6)
    {
      nat_m[(length(time_specs$fmort)+6):(length(time_specs$fmort)+z),1]<-change_m[1]
      nat_m[(length(time_specs$fmort)+6):(length(time_specs$fmort)+z),2]<-change_m[2]
    }#=====================================
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
    # pred_dat<-data.frame(SB=MMB[length(MMB)-5])
    # tmp_pred<-  predict(mod,newdata=pred_dat,se.fit=TRUE)
    # pred_rec_in<-MMB[length(MMB)-5]*exp(rnorm(1,tmp_pred$fit,tmp_pred$se.fit))
    
    #==predict recruitment
    pred_dat<-data.frame(SB=MMB[length(MMB)-5])
    tmp_pred<-  predict(mod,newdata=pred_dat,se.fit=TRUE,interval="prediction")
    #==use prediction inteval to generate recruitment
    imp_sd<-((tmp_pred$fit[1]-tmp_pred$fit[2] )/1.96)
    pred_rec_in<-MMB[length(MMB)-5]*exp(rnorm(1,tmp_pred$fit[1],imp_sd)) 
    
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
      b35_chg[m,z]<-SBPRF35_2*mean(in_rec_in_chg[m,(length(in_rec_1)+5):(length(in_rec_1)+z-2)])
    
    #==find F
    new_f<-F35
    #==after F change
    if(z>=10)
      new_f<-F35_2
    
    if(MMB[length(MMB)-1]/b35_chg[m,z] < 0.25)
      new_f <- 0
    if(MMB[length(MMB)-1]/b35_chg[m,z] >= 0.25 & MMB[length(MMB)-1]/b35_chg[m,z]<1)
      new_f <- (F35*(MMB[length(MMB)-1]/b35_chg[m,z] - alpha))/(1-alpha)
    
    in_f_chg[m,z+length(time_specs$fmort)-1]<-new_f
    
    #==predict recruitment
    pred_dat<-data.frame(SB=MMB[length(MMB)-5])
    tmp_pred<-  predict(mod,newdata=pred_dat,se.fit=TRUE,interval="prediction")
    #==use prediction inteval to generate recruitment
    imp_sd<-((tmp_pred$fit[1]-tmp_pred$fit[2] )/1.96)
    pred_rec_in<-MMB[length(MMB)-5]*exp(rnorm(1,tmp_pred$fit[1],imp_sd)) 
    
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


in_yr<-seq(1982,length.out=nrow(nat_m))



#==plot catches and biomass over time
plot_ind<-1
#==plot yields
yields_sq<-matrix(ncol=proj_yrs,nrow=nsim)
yields_ch<-matrix(ncol=proj_yrs,nrow=nsim)
for(x in 1:nsim)
  for(y in 1:(proj_yrs-1))
  {
    yields_sq[x,y]<-sum(big_proj_sq[[x]]$c_at_l[y,,1]*specs$weight,na.rm=T) + sum(big_proj_sq[[x]]$c_at_l[y,,2]*specs$weight,na.rm=T)
    yields_ch[x,y]<-sum(big_proj_ch[[x]]$c_at_l[y,,1]*specs$weight,na.rm=T) + sum(big_proj_ch[[x]]$c_at_l[y,,2]*specs$weight,na.rm=T)
  }  


#== MMB
mmb_sq<-matrix(ncol=proj_yrs,nrow=nsim)
mmb_ch<-matrix(ncol=proj_yrs,nrow=nsim)
for(x in 1:nsim)
  for(y in 1:(proj_yrs-1))
  {
    mmb_sq[x,y]<- sum(big_proj_sq[[x]]$n_at_l[y,,2]*specs$weight,na.rm=T)
    mmb_ch[x,y]<- sum(big_proj_ch[[x]]$n_at_l[y,,2]*specs$weight,na.rm=T)
  }  


#== eMMB
emmb_sq<-matrix(ncol=proj_yrs,nrow=nsim)
emmb_ch<-matrix(ncol=proj_yrs,nrow=nsim)
for(x in 1:nsim)
  for(y in 1:(proj_yrs-1))
  {
    emmb_sq[x,y]<- sum(big_proj_sq[[x]]$n_at_l[y,,2]*specs$weight*specs$fish_sel,na.rm=T)
    emmb_ch[x,y]<- sum(big_proj_ch[[x]]$n_at_l[y,,2]*specs$weight*specs$fish_sel,na.rm=T)
  }  

#== F
in_u_ch<-1-exp(-in_f_chg)
in_u_sq<-1-exp(-in_f_sq)


#===============================================================
# Plot figure 
#===============================================================
library(ggplot2)
make_dat<-function(input,input_q,input_hcr,adj=1,begin_HCR=NA)
{
  tmp<-apply(input,2,sort)
  med_ch<-as.numeric(tmp[50,])/adj
  up_ch<-as.numeric(tmp[10,])/adj
  dn_ch<-as.numeric(tmp[90,])/adj
  in_HCR<-rep(input_hcr,length(med_ch))
  if(!is.na(begin_HCR))
    in_HCR[1:begin_HCR]<-"Historical"
  out<-data.frame(median=med_ch,upper=up_ch,lower=dn_ch,Year=seq(from=1982,length.out=length(med_ch)),
                  Quantity=input_q,HCR=in_HCR)
  return(out)
}

use_colors<-c("dark grey","blue","green")

.THEME    = theme_bw(base_size = 12, base_family = "") +
  theme(strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(color="white",fill="white")) 
.COL =  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

tmp<-make_dat(input=in_u_ch,input_q="Exp_rate",input_hcr="Reference_change",begin_HCR=38)
tmp2<-make_dat(input=in_u_sq,input_q="Exp_rate",input_hcr="Status_quo",begin_HCR=38)
in_exp<-as.data.frame(rbind(tmp2,tmp))

exp_proj_m<-ggplot() +
  geom_ribbon(data=in_exp,aes(x=Year,ymin = lower, ymax = upper, fill = HCR),alpha=0.1) +
  geom_line(data=in_exp,aes(x=Year,y=median,col=HCR)) +
  .THEME +
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y=element_blank()) +
  geom_vline(xintercept=2019.5,lty=2)+ 
  xlim(1982,2065)+
  scale_color_manual(values=use_colors) +
  scale_fill_manual(values=use_colors)

tmp<-make_dat(input=in_rec_in_sq,input_q="Recruits",input_hcr="Status_quo",adj=10000,begin_HCR=38)
tmp2<-make_dat(input=in_rec_in_chg,input_q="Recruits",input_hcr="Reference_change",adj=10000,begin_HCR=38)
in_exp<-as.data.frame(rbind(tmp,tmp2))

library(png)
in_png<-readPNG('crab/silhouette_crab.png')

rec_proj_m<-ggplot() +
  geom_ribbon(data=in_exp,aes(x=Year,ymin = lower, ymax = upper, fill = HCR),alpha=0.1) +
  geom_line(data=in_exp,aes(x=Year,y=median,col=HCR)) +
  .THEME +
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y=element_blank()) +
  annotation_raster(in_png,
                    ymin=150,
                    ymax=300,
                    xmin=2035,
                    xmax=2065)+
  geom_vline(xintercept=2019.5,lty=2)+ 
  xlim(1982,2065)+
  scale_color_manual(values=use_colors) +
  scale_fill_manual(values=use_colors) +
  labs(title = "Natural mortality")+
  theme(plot.title = element_text(hjust = 0.5))
#print(rec_proj)

tmp<-make_dat(input=yields_sq[,1:(ncol(yields_sq)-2)],input_q="Yields",input_hcr="Status_quo",adj=1,begin_HCR=38)
tmp2<-make_dat(input=yields_ch[,1:(ncol(yields_sq)-2)],input_q="Yields",input_hcr="Reference_change",adj=1,begin_HCR=38)
in_exp<-as.data.frame(rbind(tmp,tmp2))

yield_proj_m<-ggplot() +
  geom_ribbon(data=in_exp,aes(x=Year,ymin = lower, ymax = upper, fill = HCR),alpha=0.1) +
  geom_line(data=in_exp,aes(x=Year,y=median,col=HCR)) +
  .THEME +
  scale_y_continuous(position = "right")+
  theme(axis.title.y=element_blank()) + 
  theme(legend.position = "none",
        axis.title.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y=element_blank()) +
  geom_vline(xintercept=2019.5,lty=2) + 
  xlim(1982,2065)+
  scale_color_manual(values=use_colors) +
  scale_fill_manual(values=use_colors)

tmp<-make_dat(input=emmb_sq[,1:(ncol(emmb_sq)-1)],input_q="MMB",input_hcr="Status_quo",begin_HCR=38)
tmp2<-make_dat(input=emmb_ch[,1:(ncol(emmb_sq)-1)],input_q="MMB",input_hcr="Reference_change",begin_HCR=38)
in_exp<-as.data.frame(rbind(tmp,tmp2))

emmb_proj_m<-ggplot() +
  geom_ribbon(data=in_exp,aes(x=Year,ymin = lower, ymax = upper, fill = HCR),alpha=0.1) +
  geom_line(data=in_exp,aes(x=Year,y=median,col=HCR)) +
  .THEME +
  theme(legend.position = 'none',
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y=element_blank()) +
  geom_vline(xintercept=2019.5,lty=2) + 
  xlim(1982,2065)+
  theme(legend.position = 'none') + 
  scale_color_manual(values=use_colors) +
  scale_fill_manual(values=use_colors)

library(gridExtra)
png("plots/fig_snow_combo.png",height=8,width=8,res=600,units='in')
grid.arrange(rec_proj_grow,rec_proj_m,exp_proj_grow,exp_proj_m,emmb_proj_grow,emmb_proj_m,yield_proj_grow,yield_proj_m,ncol=2)
dev.off()

