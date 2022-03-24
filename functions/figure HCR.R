#====================================
# Effort dynamics toy model
#====================================
alpha<-0.9 			  # sets equilibrium relative to BMSY
chi<-0.3			  # makes the exploitation rate smaller
years<-100       # years in simulation
r<-0.25         # growth rate of population
Binit<-1000     # initial bioamss
Finit<-0.01     # initial fishing mortality
alpha_HCR<-0.2
K<-c(rep(1000,30),rep(500,years-30))
#===========================================
# POpulation dynamics
#========================================
proj_pop<-function(alpha,chi,years,r,Binit,Finit,K,change_ref,use_HCR=NA,HCR_yr=NA)
{
biomass<-rep(0,years)
f_t<-rep(0,years)
yield<-rep(0,years)
exp_rate<-rep(0,years)

biomass[1]<-Binit		
f_t[1] <- Finit

for(i in 2:years)
{
  if(change_ref==1)
   f_t[i] <- f_t[i-1] * (biomass[i-1]/(alpha*K[i]/2))^chi
  if(change_ref==0)
   f_t[i] <- f_t[i-1] * (biomass[i-1]/(alpha*K[1]/2))^chi
  if(use_HCR==1 & i>HCR_yr)
  {
    f_t[i]<-r/2
    BMSY<-K[1]/2
    if(change_ref==1)
      BMSY<-K[i-1]/2
    if(biomass[i-1]/(BMSY) < 0.2)
      f_t[i] <- 0
    if(biomass[i-1]/(BMSY) >= 0.20 & biomass[i-1]/(BMSY)<1)
      f_t[i] <- (f_t[i]*((biomass[i-1]/BMSY) - alpha_HCR))/(1-alpha_HCR)
  } 
  #(F35*(MMB[length(MMB)]/b35_chg[m,z] - alpha))/(1-alpha)
  yield[i] <- ( biomass[i-1]+r*biomass[i-1]*(1-biomass[i-1]/K[i]) ) * (1-exp(-f_t[i]))
  biomass[i] <- ( biomass[i-1]+r*biomass[i-1]*(1-biomass[i-1]/K[i]) ) * exp(-f_t[i])
  exp_rate[i] <- yield[i]/biomass[i]
}

list(f_t=f_t,yield=yield,biomass=biomass,exp_rate=exp_rate)
}

#============================
# Project
#=============================
chg_sq<-proj_pop(alpha,chi,years,r,Binit,Finit,K=c(rep(1000,years/2),rep(500,years/2)),change_ref=0,use_HCR=1,HCR_yr=years/2)
chg_chg<-proj_pop(alpha,chi,years,r,Binit,Finit,K=c(rep(1000,years/2),rep(500,years/2)),change_ref=1,use_HCR=1,HCR_yr=years/2)

incex<-1.5
in_tag<-106
sm_cex<-0.7
#==use color to represent 'perceived status'
png("plots/fig_hcr.png",height=4,width=8,res=600,units='in')
incol<-c("dark grey","green","blue")
layout(cbind(c(1,2,3),c(4,4,4),c(4,4,4)))
par(mar=c(.1,.1,.1,.1),oma=c(4,4,1,4))

plot(chg_sq$biomass[(years/2):years]~seq((years/2),years),
     xlim=c(1,years+10),type="l",xaxt='n',ylim=c(0,K[1]),las=1,lty=4,col=incol[2],lwd=2)
lines(chg_sq$biomass[1:(years/2)]~seq(1,(years/2)),col=incol[1],lwd=2)
lines(chg_chg$biomass[(years/2):years]~seq((years/2),years),lty=2,col=incol[3],lwd=2)
legend("top",bty='n',"Biomass",cex=incex)
text(x=in_tag,y=K[1]/2,"BMSY",col=incol[2],cex=sm_cex)
text(x=in_tag,y=K[length(K)]/2,"BMSY",col=incol[3],cex=sm_cex)

plot(chg_sq$yield[(years/2):years]~seq((years/2),years),
     xlim=c(1,years+10),type="l",xaxt='n',ylim=c(0,110),las=1,lty=4,col=incol[2],lwd=2)
lines(chg_sq$yield[1:(years/2)]~seq(1,(years/2)),col=incol[1],lwd=2)
lines(chg_chg$yield[(years/2):years]~seq((years/2),years),lty=2,col=incol[3],lwd=2)
legend("top",bty='n',"Harvest",cex=incex)
text(x=in_tag,y=(r*K[1])/4,"MSY",col=incol[2],cex=sm_cex)
text(x=in_tag,y=(r*K[length(K)])/4,"MSY",col=incol[3],cex=sm_cex)

plot(chg_sq$exp_rate[(years/2):years]~seq((years/2),years),
     xlim=c(1,years+10),type="l",ylim=c(0,0.3),las=1,lty=4,col=incol[2],lwd=2)
lines(chg_sq$exp_rate[1:(years/2)]~seq(1,(years/2)),col=incol[1],lwd=2)
lines(chg_chg$exp_rate[(years/2):years]~seq((years/2),years),lty=2,col=incol[3],lwd=2)
legend("top",bty='n',"Harvest rate",cex=incex)
mtext(side=1,line=2.7,"Year")
text(x=in_tag,y=r/2 + .01,"UMSY",col=incol[2],cex=sm_cex)
text(x=in_tag,y=r/2 - .01,"UMSY",col=incol[3],cex=sm_cex)

plot(-10,ylab='',xlab='',ylim=c(0,r),xlim=c(0,K[1]),las=1,yaxt='n')
axis(side=4,las=1)
in_x<-c(.1,seq(.1*K[1],K[1]/2,length.out=100),seq(K[1]/2,K[1],length.out=100))
in_y<-c(0,seq(0,r/2,length.out=100),rep(r/2,100))
lines(in_y~in_x,lty=2,lwd=2)

in_x<-c(.1,seq(.1*K[length(K)],K[length(K)]/2,length.out=100),seq(K[length(K)]/2,K[length(K)],length.out=100))
in_y<-c(0,seq(0,r/2,length.out=100),rep(r/2,100))
lines(in_y~in_x,lty=4,lwd=2)

points(x=K[1]/2,y=r/2,cex=3,pch=16,col=incol[2])
points(x=K[length(K)]/2,y=r/2,cex=3,pch=16,col=incol[3])
text("A",x=K[1]/2,y=0.02+r/2,cex=3,col=incol[2])
text("B",x=K[length(K)]/2,y=0.02+r/2,cex=3,col=incol[3])
mtext(side=1,"Biomass",line=2.5)
mtext(side=4,"Harvest rate",line=2.75)

text(expression("K"^"A"),x=K[1]*.985,y=0.01,cex=3,col=incol[2])
text(expression("K"^"B"),x=K[length(K)],y=0.01,cex=3,col=incol[3])

#lines(chg_sq$exp_rate~chg_sq$biomass,lwd=3,col=3,lty=2)
lines(chg_sq$exp_rate[1:(50)]~chg_sq$biomass[1:(50)],lwd=3,col=incol[1],lty=1)
#lines(chg_chg$exp_rate~chg_chg$biomass,lwd=3,col=3,lty=4)
arrows(y0=chg_sq$exp_rate[4],
       y1=chg_sq$exp_rate[6],
       x0=chg_sq$biomass[4],
       x1=chg_sq$biomass[6],lwd=3,col=incol[1])
arrows(y0=chg_sq$exp_rate[14],
       y1=chg_sq$exp_rate[16],
       x0=chg_sq$biomass[14],
       x1=chg_sq$biomass[16],lwd=3,col=incol[1])
arrows(y0=chg_sq$exp_rate[24],
       y1=chg_sq$exp_rate[26],
       x0=chg_sq$biomass[24],
       x1=chg_sq$biomass[26],lwd=3,col=incol[1])

text('5',y=chg_sq$exp_rate[4],
     x=chg_sq$biomass[4],cex=incex)
text('15',y=chg_sq$exp_rate[15],
     x=chg_sq$biomass[15],cex=incex)
text('25',y=chg_sq$exp_rate[24],
     x=chg_sq$biomass[24],cex=incex)

text('50',y=chg_sq$exp_rate[50],
     x=chg_sq$biomass[50],cex=incex)

legend("topleft",lty=c(1,4,2),bty='n',c("Historical","Climate adaptive","Status quo"),col=incol[c(1,3,2)],lwd=2)

dev.off()



