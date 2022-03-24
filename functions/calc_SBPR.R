calc_SBPR<-function(size_trans,specs,sizes,in_m)
{
  proj_yr<-100
  const_rec<-rep(100000,proj_yr)
  nat_m<-matrix(c(rep(in_m[1],length(const_rec)),rep(in_m[2],length(const_rec))),ncol=2)
  
  vir_bio<-pop_dy(size_trans=size_trans,
                  specs=specs,
                  rec=const_rec,
                  fmort=rep(0,length(const_rec)),
                  sizes=sizes,
                  nat_m=nat_m)
  
  vir_MMB<-sum(vir_bio$n_at_l_mate[length(const_rec)-1,,2]*specs$weight,na.rm=T)
  
  min_f<-0.001
  max_f<-10
  for(x in 1:20)
  {
    in_f<-(max_f+min_f)/2  
    test<-pop_dy(size_trans=size_trans,
                 specs=specs,
                 rec=const_rec,
                 fmort=rep(in_f,length(const_rec)),
                 sizes=sizes,
                 nat_m=nat_m)
    temp_bio<-sum(test$n_at_l_mate[length(const_rec)-1,,2]*specs$weight,na.rm=T)
    if(temp_bio/vir_MMB > 0.35)
      min_f<-in_f
    if(temp_bio/vir_MMB < 0.35)
      max_f<-in_f
    
    print(temp_bio/vir_MMB)
  }
  F35<-in_f
  SBPRF35<-temp_bio/const_rec[1]
  list(F35,SBPRF35)
}
