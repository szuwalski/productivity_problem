#=================================================================
#==simluate a size structured population similar to snow crab
#=================================================================

pop_dy<-function(size_trans,
                 specs,
                 rec,
                 fmort,
                 sizes,
                 nat_m,
                 size_trans_2=NA,
                 switch_growth=NA)
{
  #==storage
  n_at_l      <-array(dim=c(length(rec),length(sizes),2))
  surv_n_at_l <-array(dim=c(length(rec),length(sizes),2))
  c_at_l      <-array(dim=c(length(rec),length(sizes),2))
  maturing_crab<-matrix(nrow=length(rec),ncol=length(sizes))
  n_at_l_mate <-array(dim=c(length(rec),length(sizes),2))
  
  #==set up initial size structure
  n_at_l[1,,1]<-specs$init_immature
  n_at_l[1,,2]<-specs$init_mature
  
  for(x in 1:(length(rec)-1))
  {
    imm_to_mat<-NULL
    for(y in 1:2)
    {
      temp<-n_at_l[x,,y]
      #==survey observations
      surv_n_at_l[x,,y]<-temp*specs$surv_sel
      #==natural mortality (part 1)
      temp<-temp*exp(-nat_m[x,y]*0.625)
      
      #==fishery (0.625 years after survey)
      c_at_l[x,,y] <- temp - temp*exp(-fmort[x]*specs$fish_sel)
      temp<- temp*exp(-fmort[x]*specs$fish_sel)
      n_at_l_mate[x,,y] <- temp
      
      #==growth & rec (immediately after fishery)
      if(y==1) # immature
      {
        #==growth only for immature 
        use_size_trans<-size_trans
        if(!is.na(switch_growth) & x>=switch_growth)
          use_size_trans<-size_trans_2
        temp<-temp%*%use_size_trans
        imm_to_mat <- temp*specs$prob_maturing
        #imm_to_mat <- temp*specs$prob_maturing_alt   # test knifeedge prob of maturing
        
        #==recruitment only to immature
        temp <- temp + specs$prop_rec*rec[x]
        
        #==update numbers at length 
        n_at_l[x+1,,y] <- temp*(1-specs$prob_maturing)
      }
      
      #==update numbers at length (mature)
      if(y == 2)
      {
        maturing_crab[x,] <- imm_to_mat  
        n_at_l[x+1,,y] <- temp + imm_to_mat
      }
      #==natural mortality (part 2)
      temp<-temp*exp(-nat_m[x,y]*0.375)
      
    }
  }
  
  list(n_at_l=n_at_l,c_at_l=c_at_l,mat_crab=maturing_crab,n_at_l_mate=n_at_l_mate)
}
