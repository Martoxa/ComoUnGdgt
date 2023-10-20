bayspR_tex_analog<-function(raw,prior_mean,prior_std,search_tol,runname,varargin=1000,data="area",complete=TRUE){
  
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
    Params<-readRDS(file="SST_param_ana.rds")
    DataIn<-readRDS(file = "SST_DataIn.rds")
  } else if(runname=="subT"){
    Params<-readRDS(file="SubT_param_ana.rds")
    DataIn<-readRDS(file = "subT_DataIn.rds")
  }

  ind_s<-round(seq(1,length(Params$tau2.samples),length.out=Nsamps))
  alpha_samples<-Params$alpha.samples[,ind_s]
  beta_samples<-Params$beta.samples[,ind_s]
  tau2_samples<-Params$tau2.samples[ind_s]
  
  Nd<-length(dats)
  
  Output_Preds<-matrix(NaN,nrow=Nd,ncol=3)
  Output_AnLocs<-NA
  Output_PriorMean<-prior_mean
  Output_PriorStd<-prior_std
  
  if(ens_sel==1){Output_PredsEns<-NA} else {Output_PredsEns<-NA} #maybe move to later
  
  N_bg<-dim(matrix(unlist(DataIn[[1]][2]),ncol=2))[1]
  
  spatialMean<-unlist(lapply(1:N_bg,function(a){mean(unlist(DataIn[[1]][8])[unlist(DataIn[[1]][5])==a])}))
  inder_g<- spatialMean >= mean(dats)-search_tol & spatialMean <= mean(dats)+search_tol
  if(sum(inder_g)==0){stop("Your search tolerance is too narrow")}
 
  alpha_samples<-as.data.frame(alpha_samples)[inder_g,]
  beta_samples<-as.data.frame(beta_samples)[inder_g,]
  tau2_samples<-kronecker(matrix(1,dim(alpha_samples)[1],1),tau2_samples)
  alpha_samples<-matrix(as.matrix(alpha_samples),nrow=1,ncol=dim(alpha_samples)[1]*dim(alpha_samples)[2])
  beta_samples<-matrix(as.matrix(beta_samples),nrow=1,ncol=dim(beta_samples)[1]*dim(beta_samples)[2])
  
  Output_AnLocs<-matrix(unlist(DataIn[[1]][2]),ncol=2)[inder_g,]
  
  
  pmu<-kronecker(matrix(1,1,dim(alpha_samples)[2]),matrix(1,Nd,1)*prior_mean)
  pinv_cov<-kronecker(matrix(1,Nd,dim(alpha_samples)[2]),prior_std)^(-2)
  sigmaS<-sqrt(tau2_samples)
  
  post_mean_num<-pinv_cov*pmu + kronecker(matrix(1,Nd,1),t(sigmaS))^(-2)*kronecker(matrix(1,Nd,1),as.matrix(beta_samples))*(dats- kronecker(matrix(1,Nd,1),as.matrix(alpha_samples)))
  post_mean_den<-pinv_cov+kronecker(matrix(1,Nd,1),as.matrix(beta_samples))^2 * kronecker(matrix(1,Nd,1),t(sigmaS))^(-2)
  post_mean<-post_mean_num / post_mean_den
  post_sig<-sqrt(post_mean_den^(-1))
  Preds<-post_mean+matrix(rnorm(Nd*dim(alpha_samples)[2],0,1),nrow = Nd,ncol=dim(alpha_samples)[2])*post_sig
  
  Output_Preds<- if(length(dats)==1){quantile(sort(Preds),probs=c(0.05,0.5,0.95))} else if(length(dats)>1){t(apply(apply(Preds,1,sort),2,quantile,probs=c(0.05,0.5,0.95)))}
  
  if(ens_sel==1){
    Output_PredsEns<-Preds
  }
  
  Out<-list("Predictions"=Output_Preds,"AnLocs"=Output_AnLocs,"PriorMean"=Output_PriorMean,"PriorStd"=Output_PriorStd,"Ensemble"=Output_PredsEns)
  
}