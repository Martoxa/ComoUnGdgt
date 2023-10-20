bayspR_tex<-function(raw,lon,lat,prior_std,runname,varargin=1000,data="area",complete=TRUE){
  
  if(data=="area"){dats<-TEX86(raw,complete=complete,data=data)
  } else if (data=="fa"){dats<-TEX86(raw,complete=complete,data=data)
  } else if (data=="tex"){dats<-raw}
  
  if(length(varargin)>2){stop("varargin accepts no more than two arguments")
  } else if(length(varargin)==2){
    Nsamps<-varargin[1]
    ens_sel<-varargin[2]
  } else if(length(varargin)==1) {
    Nsamps<-varargin
    ens_sel<-0
  }
  
  if(runname=="SST"){
    Params<-readRDS(file="SST_param_std.rds")
    Obs<-readRDS(file = "obsSST.rds")
  } else if(runname=="subT"){
    Params<-readRDS(file="SubT_param_std.rds")
    Obs<-readRDS(file = "obssubT.rds")
  }
  
  grid_half_space<-10
  min_num<-1
  max_dist<-500
  
  ind_s<-round(seq(1,length(Params$tau2.samples),length.out=Nsamps))
  alpha_samples_comp<-Params$alpha.samples.comp[,ind_s]
  beta_samples_comp<-Params$beta.samples.comp[,ind_s]
  tau2_samples<-Params$tau2.samples[ind_s]
  
  Nd<-length(dats)
  
  Output_Preds<-matrix(NaN,nrow=Nd,ncol=3)
  Output_SiteLoc<-as.matrix(cbind(lon,lat))
  Output_GridLoc<-matrix(NaN,nrow = 2,ncol=1)
  Output_PriorMean<-matrix(NaN,nrow = 1,ncol=1)
  Output_PriorStd<-prior_std
  
  if(ens_sel==1){Output_PredsEns<-matrix(NaN,nrow = Nd,ncol=Nsamps)} else {Output_PredsEns<-NA}
  
  dists_prior<-RthChordDistances(Output_SiteLoc,Obs$locs.st.obs)
  sorted_dists_prior<-if(length(lon)==1){sort(dists_prior)} else if(length(lon)>1){apply(dists_prior,2,sort)}
  inds_dists_prior<-if(length(lon)==1){match(sorted_dists_prior,dists_prior)} else if(length(lon)>1){which(sorted_dists_prior==dists_prior,arr.ind = TRUE)[,1]}
  num_below_dist<-which(sorted_dists_prior == tail(c(sorted_dists_prior)[c(sorted_dists_prior)<max_dist],n=1),arr.ind = FALSE)
  
  if(num_below_dist>min_num){prior_mean<-mean(Obs$st.obs.ave.vec[inds_dists_prior[seq(from=1,to=num_below_dist,by=1)]])
  } else {prior_mean<-mean(Obs$st.obs.ave.vec[inds_dists_prior[seq(from=1,to=min_num,by=1)]])}
  Output_PriorMean<-prior_mean
  
  inder_g<- which(abs(Params$Locs.Comp[,1]-lon)<=grid_half_space & abs(Params$Locs.Comp[,2]-lat)<= grid_half_space)
  alpha_samples_comp<-as.data.frame(alpha_samples_comp)[inder_g,]
  beta_samples_comp<-as.data.frame(beta_samples_comp)[inder_g,]
  Output_GridLoc<-Params$Locs.Comp[inder_g,]
  
  pmu<-kronecker(matrix(1,1,Nsamps),matrix(1,Nd,1)*prior_mean)
  pinv_cov<-kronecker(matrix(1,Nd,Nsamps),prior_std)^(-2)
  sigmaS<-sqrt(tau2_samples)
  
  post_mean_num<-pinv_cov*pmu + kronecker(matrix(1,Nd,1),t(sigmaS))^(-2)*kronecker(matrix(1,Nd,1),as.matrix(beta_samples_comp))*(dats- kronecker(matrix(1,Nd,1),as.matrix(alpha_samples_comp)))
  post_mean_den<-pinv_cov+kronecker(matrix(1,Nd,1),as.matrix(beta_samples_comp))^2 * kronecker(matrix(1,Nd,1),t(sigmaS))^(-2)
  post_mean<-post_mean_num / post_mean_den
  post_sig<-sqrt(post_mean_den^(-1))
  Preds<-post_mean+matrix(rnorm(Nd*Nsamps,0,1),nrow = Nd,ncol=Nsamps)*post_sig
  
  Output_Preds<- if(length(dats)==1){quantile(sort(Preds),probs=c(0.05,0.5,0.95))} else if(length(dats)>1){t(apply(apply(Preds,1,sort),2,quantile,probs=c(0.05,0.5,0.95)))}
  
  if(ens_sel==1){
    Output_PredsEns<-Preds
  }
  
  Out<-list("Predictions"=Output_Preds,"SiteLoc"=Output_SiteLoc,"GridLoc"=Output_GridLoc,"PriorMean"=Output_PriorMean,"PriorStd"=Output_PriorStd,"Ensemble"=Output_PredsEns)

}