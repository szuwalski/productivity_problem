#==read in parameters
library(dplyr)
library(ggplot2)
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
unq_stk<-unique(ram_sp$stocklong)

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


#===========================================================================
# Project under different proportions of stocks decreasing vs. increasing
#===========================================================================
pr_yrs<-50
nsim<-100
RAM_sq<-list(list())
RAM_ch<-list(list())
stocks<-unique(tmp$stocklong)

up_down<-seq(0.1,0.9,0.05)
biomass<-array(dim=c(length(stocks),pr_yrs,nsim,2,length(up_down)))
yield  <-array(dim=c(length(stocks),pr_yrs,nsim,2,length(up_down)))
f_t    <-array(dim=c(length(stocks),pr_yrs,nsim,2,length(up_down)))


for(p in 1:length(up_down))
{
  
k_err<-matrix(nrow=nsim,ncol=length(stocks))
for(x in 1:nsim)
{
k_err_up<-rnorm(length(stocks),1.5,.1)
k_err_dn<-rnorm(length(stocks),0.5,.1)
k_err_in<-k_err_up
ind<-sample(x=seq(1,length(stocks)),size=length(stocks)*up_down[p],replace=FALSE)
k_err_in[ind]<-k_err_dn[ind]
k_err[x,]<-k_err_in
}

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
   
   yield[x,,z,1,p]<-RAM_sq$yield
   yield[x,,z,2,p]<-RAM_ch$yield
   biomass[x,,z,1,p]<-RAM_sq$biomass[-1]
   biomass[x,,z,2,p]<-RAM_ch$biomass[-1]
   f_t[x,,z,1,p]<-RAM_sq$f_t
   f_t[x,,z,2,p]<-RAM_ch$f_t
  }
  }
}
}

tot_y_sq<-list(list())
tot_y_ch<-list(list())
tot_b_sq<-list(list())
tot_b_ch<-list(list())

#==sum over stocksf
for(x in 1:length(up_down))
{
tot_y_sq[[x]]<-apply(yield[,,,1,x],2:3,sum,na.rm=T)
tot_y_ch[[x]]<-apply(yield[,,,2,x],2:3,sum,na.rm=T)
tot_b_sq[[x]]<-apply(biomass[,,,1,x],2:3,sum,na.rm=T)
tot_b_ch[[x]]<-apply(biomass[,,,2,x],2:3,sum,na.rm=T)
}

#==get ready for plot
plot_dat_b_sq<-data.frame(upper=rep(NA,length(up_down)),
                     lower=rep(NA,length(up_down)),
                     mid=rep(NA,length(up_down)),
                     Quantity=rep("Biomass",length(up_down)),
                     HCR=rep("Status quo",length(up_down)),
                     Proportion=rep(NA,length(up_down)))

plot_dat_b_ch<-plot_dat_b_sq
plot_dat_b_ch$HCR<-"Climate adaptive"
plot_dat_y_sq<-plot_dat_b_sq
plot_dat_y_sq$Quantity<-"Yield"
plot_dat_y_ch<-plot_dat_y_sq
plot_dat_y_ch$HCR<-"Climate adaptive"

for(x in 1:length(up_down))
{
  tmp_b_sq<-sort(tot_b_sq[[x]][pr_yrs,] )  / 1000000
  tmp_b_ch<-sort(tot_b_ch[[x]][pr_yrs,] )  / 1000000
  tmp_y_sq<-sort(tot_y_sq[[x]][pr_yrs,] )  / 1000000
  tmp_y_ch<-sort(tot_y_ch[[x]][pr_yrs,] )  / 1000000

  plot_dat_b_sq$upper[x]<-tmp_b_sq[length(tmp_b_sq)] 
  plot_dat_b_sq$lower[x]<-tmp_b_sq[1]  
  plot_dat_b_sq$mid[x]<-median(tmp_b_sq) 
   
  plot_dat_b_ch$upper[x]<-tmp_b_ch[length(tmp_b_ch)] 
  plot_dat_b_ch$lower[x]<-tmp_b_ch[1]  
  plot_dat_b_ch$mid[x]<-median(tmp_b_ch)
  
  plot_dat_y_sq$upper[x]<-tmp_y_sq[length(tmp_y_sq)] 
  plot_dat_y_sq$lower[x]<-tmp_y_sq[1]  
  plot_dat_y_sq$mid[x]<-median(tmp_y_sq)
  
  plot_dat_y_ch$upper[x]<-tmp_y_ch[length(tmp_y_ch)] 
  plot_dat_y_ch$lower[x]<-tmp_y_ch[1]  
  plot_dat_y_ch$mid[x]<-median(tmp_y_ch)
  
  plot_dat_y_ch$Proportion[x]<-up_down[x]
  plot_dat_b_ch$Proportion[x]<-up_down[x]
  plot_dat_y_sq$Proportion[x]<-up_down[x]
  plot_dat_b_sq$Proportion[x]<-up_down[x]
}

in_bio<-rbind(plot_dat_b_ch,plot_dat_b_sq)
in_yield<-rbind(plot_dat_y_ch,plot_dat_y_sq)

#==do one where sum within a region and plot those
library(reshape)
library(ggplot2)
library(gridExtra)
#==plot
use_colors<-c("blue","green")

.THEME    = theme_bw(base_size = 12, base_family = "") +
  theme(strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(color="white",fill="white")) 
.COL =  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

both_proj<-ggplot() +
  geom_ribbon(data=in_bio,aes(x=Proportion,ymin = lower, ymax = upper, fill = HCR),alpha=0.1) +
  geom_line(data=in_bio,aes(x=Proportion,y=mid,col=HCR)) +
  geom_ribbon(data=in_yield,aes(x=Proportion,ymin = lower, ymax = upper, fill = HCR),alpha=0.1) +
  geom_line(data=in_yield,aes(x=Proportion,y=mid,col=HCR)) +
  .THEME +
  theme(legend.position = c(.2,.5)) + 
  scale_color_manual(values=use_colors) +
  scale_fill_manual(values=use_colors) + ylim(0,205) + xlim(0.1,0.9) + 
  annotate("text",x=0.75,y=160,label="Biomass") +
  annotate("text",x=0.75,y=50,label="Yield") +
  labs(y="Biomass (1,000,000 t)") +
  xlab("Proportion of populations decreasing in productivity")+
  geom_vline(xintercept=0.5,lty=2)

png("plots/prop_unc_proj.png",height=4,width=7.5,res=600,units='in')
print(both_proj)
dev.off()



