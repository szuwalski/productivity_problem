#==read in parameters
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(plotrix)
library(reshape)
library(gridExtra)
source('functions/color_legend2.R')

ram_fits_MM<-read.csv("RAM/RAM_sp_fits.csv")
ram_sp<-read.csv("RAM/RAM_surplus_prod_inputs.csv")
.THEME    = theme_bw(base_size = 12, base_family = "") +
  theme(strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(color="white",fill="white")) 
.COL =  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))


tmp <- ram_fits_MM %>%
       filter(Fit.Summary == "[Sufficient Data][Converged][Passed Filters][Usable]" &
                !is.na(phi.SPmod))

max_BMSY<-max(tmp$Bmsy.SPmod,tmp$Bmsy.assess,na.rm=T)
colrange<-seq(0,log(max_BMSY),length.out=1000)
cols 		<-colorRampPalette(brewer.pal(9,'RdYlBu'))(length(colrange))
in_leg  <-c(0,20,700,18000,480000,12600000)

#==plot the production fits
plot_out<-1
if(plot_out==1)
{
unq_stk<-unique(ram_sp$stocklong)

png("plots/RAM_fits.png",height=10,width=8,res=600,units='in')
par(mfrow=c(23,23),mar=c(0.1,.1,.1,.1),oma=c(4,4,1,9))
for(x in 1:length(unq_stk))
{
temp<-filter(ram_sp,stocklong==unq_stk[x]) 

#==get SP fit pars
temp2<-filter(tmp,stocklong==unq_stk[x])
temp2<-temp2[nrow(temp2),]

if(length(temp2$phi.SPmod)>0)
{
in_b<-seq(0,max(temp$B,na.rm=T),length.out=100)

if(!is.na(temp2$ERmsy.SPmod) & !is.na(temp2$Bmsy.SPmod))
 in_p<-((temp2$phi.SPmod/(temp2$phi.SPmod-1))*in_b*temp2$ERmsy.SPmod) - ((temp2$ERmsy.SPmod*in_b^temp2$phi.SPmod)/((temp2$phi.SPmod-1)*(temp2$Bmsy.SPmod^(temp2$phi.SPmod-1))))

if(is.na(temp2$ERmsy.SPmod) & !is.na(temp2$Bmsy.SPmod))
  in_p<-((temp2$phi.SPmod/(temp2$phi.SPmod-1))*in_b*temp2$Umsy.assess) - ((temp2$Umsy.assess*in_b^temp2$phi.SPmod)/((temp2$phi.SPmod-1)*(temp2$Bmsy.SPmod^(temp2$phi.SPmod-1))))

if(!is.na(temp2$ERmsy.SPmod) & is.na(temp2$Bmsy.SPmod))
  in_p<-((temp2$phi.SPmod/(temp2$phi.SPmod-1))*in_b*temp2$ERmsy.SPmod) - ((temp2$ERmsy.SPmod*in_b^temp2$phi.SPmod)/((temp2$phi.SPmod-1)*(temp2$Bmsy.assess^(temp2$phi.SPmod-1))))

if(is.na(temp2$ERmsy.SPmod) & is.na(temp2$Bmsy.SPmod))
  in_p<-((temp2$phi.SPmod/(temp2$phi.SPmod-1))*in_b*temp2$Umsy.assess) - ((temp2$Umsy.assess*in_b^temp2$phi.SPmod)/((temp2$phi.SPmod-1)*(temp2$Bmsy.assess^(temp2$phi.SPmod-1))))


in_BMSY<-temp2$Bmsy.assess
if(is.na(in_BMSY))
  in_BMSY<-temp2$Bmsy.SPmod

incols<-cols[which(abs(colrange-log(in_BMSY))==min(abs(colrange-log(in_BMSY))))]  


in_pred<-data.frame(pred_P=in_p,B=in_b)

plot(temp$P~temp$B,xlim=c(0,max(temp$B,na.rm=T)),ylim=c(min(c(0,temp$P),na.rm=T),max(c(temp$P,in_pred$pred_P),na.rm=T)),
     bty='n',yaxt='n',xaxt='n',pch=16,cex=0.8,col=incols)
 lines(in_pred$pred_P~in_pred$B)
 abline(h=0,lty=2)
 
}
}
mtext("Surplus production",side=2,outer=T,line=1)
mtext("Biomass",side=1,outer=T,line=1)

par(xpd=NA)
color.legend2(2*max(temp$B),max(temp$P)*3,3.1*max(temp$B),max(temp$P)*30,rect.col=cols,align="rb",
              legend=in_leg,cex=1,gradient='y')
text(x=2.5*max(temp$B),y=max(temp$P)*30.5,"BMSY (t)")



dev.off()
}


#==projection specifications + function
proj_RAM<-function(years,ERmsy,Binit,Finit,TBmsy,p,change_ref,alpha_HCR)
{
  biomass<-rep(0,years)
  f_t<-rep(0,years)
  yield<-rep(0,years)
  biomass[1]<-Binit		
  f_t[1] <- Finit
  
  for(i in 1:years)
  {
      f_t[i]<-ERmsy
        BMSY<-TBmsy[1]
      if(change_ref==1)
        BMSY<-TBmsy[i]

      if(biomass[i]/(BMSY) < 0.2)
        f_t[i] <- 0
      if(biomass[i]/(BMSY) >= 0.20 & biomass[i]/(BMSY)<1)
        f_t[i] <- (f_t[i]*((biomass[i]/BMSY) - alpha_HCR))/(1-alpha_HCR)

      yield[i] <-    (biomass[i]) * f_t[i]
    biomass[i+1] <-  biomass[i] + ((p/(p-1))*biomass[i]*ERmsy) - ((ERmsy*biomass[i]^p)/((p-1)*(TBmsy[i]^(p-1)))) - yield[i]

    if(biomass[i+1]<=0) biomass[i+1]<-0.001
  }
  
  list(f_t=f_t,yield=yield,biomass=biomass)
}


#============================
# Project
#=============================
pr_yrs<-50
nsim<-100
RAM_sq<-list(list())
RAM_ch<-list(list())

stocks<-unique(tmp$stocklong)
k_err<-matrix(nrow=nsim,ncol=length(stocks))
for(x in 1:nsim)
{
k_err_up<-rnorm(length(stocks),1.5,.1)
k_err_dn<-rnorm(length(stocks),0.5,.1)
k_err_in<-k_err_up
ind<-sample(x=seq(1,length(stocks)),size=length(stocks)/2,replace=FALSE)
k_err_in[ind]<-k_err_dn[ind]
k_err[x,]<-k_err_in
}

biomass<-array(dim=c(length(stocks),pr_yrs,nsim,2))
yield  <-array(dim=c(length(stocks),pr_yrs,nsim,2))
f_t    <-array(dim=c(length(stocks),pr_yrs,nsim,2))

tmp <- ram_fits_MM %>%
  filter(Fit.Summary == "[Sufficient Data][Converged][Passed Filters][Usable]" &
           !is.na(phi.SPmod))

for(x in 1:length(stocks))
{
  #==get SP fit pars
  temp<-filter(ram_sp, stocklong == stocks[x]) 
  temp<-temp[complete.cases(temp),]
  temp2<-filter(tmp,stocklong==unq_stk[x])
  temp2<-temp2[nrow(temp2),]
  
  if(length(temp2$phi.SPmod)>0)
  {
  in_phi<-temp2$phi.SPmod
  in_umsy<-temp2$ERmsy.SPmod
  in_bmsy<-temp2$Bmsy.SPmod
  if(is.na(in_umsy))
    in_umsy<-temp2$Umsy.assess
  if(is.na(in_bmsy))
    in_bmsy<-temp2$Bmsy.assess
  in_init_b<-in_bmsy
  for(z in 1:nsim)  
  {
    RAM_sq<-proj_RAM(years=pr_yrs,
                 ERmsy=in_umsy,
                 Binit=in_init_b,
                 Finit=in_umsy,
                 TBmsy=c(rep(in_bmsy,pr_yrs/2),rep(in_bmsy*k_err[z,x],pr_yrs/2)),
                 p=in_phi,
                 change_ref=0,
                 alpha_HCR=0.2)

   RAM_ch<-proj_RAM(years=pr_yrs,
                    ERmsy=in_umsy,
                   Binit=in_init_b,
                   Finit=in_umsy,
                   TBmsy=c(rep(in_bmsy,pr_yrs/2),rep(in_bmsy*k_err[z,x],pr_yrs/2)),
                   p=in_phi,
                   change_ref=1,
                   alpha_HCR=0.2)
   
   yield[x,,z,1]<-RAM_sq$yield
   yield[x,,z,2]<-RAM_ch$yield
   biomass[x,,z,1]<-RAM_sq$biomass[-1]
   biomass[x,,z,2]<-RAM_ch$biomass[-1]
   f_t[x,,z,1]<-RAM_sq$f_t
   f_t[x,,z,2]<-RAM_ch$f_t
  }
  }
}


#==sum over stocks
tot_y_sq<-apply(yield[,,,1],2:3,sum,na.rm=T)
tot_y_ch<-apply(yield[,,,2],2:3,sum,na.rm=T)
tot_b_sq<-apply(biomass[,,,1],2:3,sum,na.rm=T)
tot_b_ch<-apply(biomass[,,,2],2:3,sum,na.rm=T)


use_colors<-c("blue","green")
adj<-1000000
tmp<-apply(tot_y_sq,1,sort)
med_ch<-as.numeric(tmp[50,])/adj
years<-seq(from=2015,length.out=length(med_ch))
.THEME    = theme_bw(base_size = 12, base_family = "") +
  theme(strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(color="white",fill="white")) 
.COL =  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

#====================================
# plot change up vs. change down
#====================================
long_yield<-melt(yield)
names(long_yield)<-c("Stock","Year","Sim","HCR","Yield")

long_kerr<-melt(k_err)
names(long_kerr)<-c("Sim","Stock","Mult")
mrg_data<-merge(long_yield,long_kerr)

#==pull based on HCR and Mult
#==take top, bottom, median
tmp2_y_tot<-mrg_data %>%
  group_by(Year,Sim,HCR) %>%
  summarize(tot_y=sum(Yield,na.rm=T))
tmp3_y_tot<-tmp2_y_tot %>%
  group_by(Year,HCR) %>%
  summarize(upper=max(tot_y)/adj,lower=min(tot_y)/adj,med=median(tot_y)/adj)

tmp1_y_dn<-filter(mrg_data,Mult<1)
tmp2_y_dn<-tmp1_y_dn %>%
  group_by(Year,Sim,HCR) %>%
  summarize(tot_y=sum(Yield,na.rm=T))
tmp3_y_dn<-tmp2_y_dn %>%
  group_by(Year,HCR) %>%
  summarize(upper=max(tot_y)/adj,lower=min(tot_y)/adj,med=median(tot_y)/adj)

tmp1_y_up<-filter(mrg_data,Mult>=1)
tmp2_y_up<-tmp1_y_up %>%
  group_by(Year,Sim,HCR) %>%
  summarize(tot_y=sum(Yield,na.rm=T))
tmp3_y_up<-tmp2_y_up %>%
  group_by(Year,HCR) %>%
  summarize(upper=max(tot_y)/adj,lower=min(tot_y)/adj,med=median(tot_y)/adj)

#==biomass
long_bio<-melt(biomass)
names(long_bio)<-c("Stock","Year","Sim","HCR","Biomass")
mrg_data_b<-merge(long_bio,long_kerr)

#==pull based on HCR and Mult
#==take top, bottom, median

tmp2_b_tot<-mrg_data_b %>%
  group_by(Year,Sim,HCR) %>%
  summarize(tot_b=sum(Biomass,na.rm=T))
tmp3_b_tot<-tmp2_b_tot %>%
  group_by(Year,HCR) %>%
  summarize(upper=max(tot_b)/adj,lower=min(tot_b)/adj,med=median(tot_b)/adj)

tmp1_b_dn<-filter(mrg_data_b,Mult<1)
tmp2_b_dn<-tmp1_b_dn %>%
  group_by(Year,Sim,HCR) %>%
  summarize(tot_b=sum(Biomass,na.rm=T))
tmp3_b_dn<-tmp2_b_dn %>%
  group_by(Year,HCR) %>%
  summarize(upper=max(tot_b)/adj,lower=min(tot_b)/adj,med=median(tot_b)/adj)

tmp1_b_up<-filter(mrg_data_b,Mult>=1)
tmp2_b_up<-tmp1_b_up %>%
  group_by(Year,Sim,HCR) %>%
  summarize(tot_b=sum(Biomass,na.rm=T))
tmp3_b_up<-tmp2_b_up %>%
  group_by(Year,HCR) %>%
  summarize(upper=max(tot_b)/adj,lower=min(tot_b)/adj,med=median(tot_b)/adj)

in_HCR<-c("Status quo","Climate adaptive")

tmp3_y_tot$HCR<-as.factor(in_HCR[tmp3_y_tot$HCR])
tmp3_b_tot$HCR<-as.factor(in_HCR[tmp3_b_tot$HCR])
tmp3_y_dn$HCR<-as.factor(in_HCR[tmp3_y_dn$HCR])
tmp3_y_up$HCR<-as.factor(in_HCR[tmp3_y_up$HCR])
tmp3_b_dn$HCR<-as.factor(in_HCR[tmp3_b_dn$HCR])
tmp3_b_up$HCR<-as.factor(in_HCR[tmp3_b_up$HCR])

tmp3_y_tot$Year<-years[tmp3_y_tot$Year]
tmp3_b_tot$Year<-years[tmp3_b_tot$Year]
tmp3_y_dn$Year<-years[tmp3_y_dn$Year]
tmp3_y_up$Year<-years[tmp3_y_up$Year]
tmp3_b_dn$Year<-years[tmp3_b_dn$Year]
tmp3_b_up$Year<-years[tmp3_b_up$Year]


both_tot<-ggplot() +
  geom_ribbon(data=tmp3_y_tot,aes(x=Year,ymin = lower, ymax = upper, fill = HCR),alpha=0.1) +
  geom_line(data=tmp3_y_tot,aes(x=Year,y=med,col=HCR))+
  geom_ribbon(data=tmp3_b_tot,aes(x=Year,ymin = lower, ymax = upper, fill = HCR),alpha=0.1) +
  geom_line(data=tmp3_b_tot,aes(x=Year,y=med,col=HCR)) +
  .THEME +
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x = element_blank()) +
  annotate("text",x=2030,y=145,label="Biomass") +
  annotate("text",x=2029,y=45,label="Yield") +
  annotate("text",x=2029,y=162,label="TOTAL") +
  xlim(2025,2065)+
  ylim(0,165)+
  ylab(" ")+
  scale_color_manual(values=use_colors) +
  scale_fill_manual(values=use_colors) 

both_down<-  ggplot() +
  geom_ribbon(data=tmp3_y_dn,aes(x=Year,ymin = lower, ymax = upper, fill = HCR),alpha=0.1) +
  geom_line(data=tmp3_y_dn,aes(x=Year,y=med,col=HCR))+
  geom_ribbon(data=tmp3_b_dn,aes(x=Year,ymin = lower, ymax = upper, fill = HCR),alpha=0.1) +
  geom_line(data=tmp3_b_dn,aes(x=Year,y=med,col=HCR)) +
  .THEME +
  theme(legend.position = c(.7,.75),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        legend.background=element_blank()) + 
  annotate("text",x=2030,y=70,label="Biomass") +
  annotate("text",x=2029,y=25,label="Yield") +
  annotate("text",x=2029,y=150,label="DECLINE") +
  xlim(2025,2065)+
  ylim(0,150)+
  ylab("Million metric tons")+
  scale_color_manual(values=use_colors) +
  scale_fill_manual(values=use_colors) 

both_up<-  ggplot() +
  geom_ribbon(data=tmp3_y_up,aes(x=Year,ymin = lower, ymax = upper, fill = HCR),alpha=0.1) +
  geom_line(data=tmp3_y_up,aes(x=Year,y=med,col=HCR))+
  geom_ribbon(data=tmp3_b_up,aes(x=Year,ymin = lower, ymax = upper, fill = HCR),alpha=0.1) +
  geom_line(data=tmp3_b_up,aes(x=Year,y=med,col=HCR)) +
  .THEME +
  theme(legend.position = "none") + 
  ylab(" ")+
  annotate("text",x=2030,y=68,label="Biomass") +
  annotate("text",x=2029,y=25,label="Yield") +
  annotate("text",x=2029,y=150,label="INCREASE") +
  xlim(2025,2065)+
  ylim(0,150)+
  scale_color_manual(values=use_colors) +
  scale_fill_manual(values=use_colors) 


png("plots/tot_fig_RAM.png",height=8,width=4,res=600,units='in')
grid.arrange(both_tot,both_down,both_up)
dev.off()

#==differences between scenarios
sq_y<-filter(tmp3_y_tot,HCR=="Status quo"&Year==2060)$med
ch_y<-filter(tmp3_y_tot,HCR=="Climate adaptive"&Year==2060)$med

(sq_y-ch_y)/ch_y

sq_b<-filter(tmp3_b_tot,HCR=="Status quo"&Year==2060)$med
ch_b<-filter(tmp3_b_tot,HCR=="Climate adaptive"&Year==2060)$med

(sq_b-ch_b)/ch_b

sq_y<-filter(tmp3_y_dn,HCR=="Status quo"&Year==2060)$med
ch_y<-filter(tmp3_y_dn,HCR=="Climate adaptive"&Year==2060)$med

(sq_y-ch_y)/ch_y

sq_b<-filter(tmp3_b_dn,HCR=="Status quo"&Year==2060)$med
ch_b<-filter(tmp3_b_dn,HCR=="Climate adaptive"&Year==2060)$med

(sq_b-ch_b)/ch_b

sq_y<-filter(tmp3_y_up,HCR=="Status quo"&Year==2060)$med
ch_y<-filter(tmp3_y_up,HCR=="Climate adaptive"&Year==2060)$med

(sq_y-ch_y)/ch_y

