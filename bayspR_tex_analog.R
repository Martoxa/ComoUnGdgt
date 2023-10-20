bayspR_tex_analog<-function(raw,prior_mean,prior_std,search_tol,runname,varargin=1000,data="area",complete=TRUE){
#
# INPUTS:
# raw      - a matrix of HPLC-MS peak area values  (6 x N) or (N x 6) 
#if data="area" (default), or a scalar or vector of TEX86 values (1 x N) or (N x 1). 
# prior_mean - prior mean on temperature. assumed same for all if dats is a
# vector
# prior_std  - prior std on temperature. assumed same for all. 
# search_tol - tolerance for finding analog locations. comparison is
# between the mean of dats and the mean tex value within each large
# gridcell. 
# runname   - enter 'SST' for sea-surface temperature or 'subT' for subsurface T
#
# varargin  - used to set the number of posterior draws of parameters to
# mix across, and whether to save the 5th/50th/95th percentiles of the
# predictions, or the whole ensemble as well.
# - if left empty, then 1000 draws are used and the ensemble of predictions
# is not saved. 
# - if only one argument, it is the number of draws (cannot exceed 15000),
# and the ensemble is not saved. Note that the first N_samps are used.
# - if two arguments, the first gives the number of draws, while the second
# is an indicator:
# 0: (default) save only the 5th/50th/95th percentiles of the predictions. 
# 1: save the whole ensemble as well. 
#
# data: A string corresponding the type of data to run in the model.
# "area" = Peak areas for each isoGDGT (default option)
# "fa" = Fractional abundance of each isoGDGT
# "tex" = Calculated TEX86 values  
# complete: Determines how is TEX86 calculated
# TRUE = Calculated exactly as per Schouten, et al., 2002 (default option)
# FALSE = Omits missing compounds. A warning is generated and the ommitted compounds are shown.
#
#Output structure:
#Output (list)
#Predictions      - Nd by 3 array, where Nd=length(dats), giving the 5/50/95
#percentiles of the predictions. 
#AnLocs=- Centroids of the grid cells selected as analogs, as [lon,
#lat] pairs 
#PriorMean - The prior mean as calculated above. 
#PriorStd  - The prior std as input. 
#PredsEns  - Nd by No. analog locatioins by Nsamps array of predictions. Only included if
#ens_sel==1. Note that the second dimension corresponds to the locations
#in .AnLocs.
    
  #Deal with the input parameters and decide how to handle the input data
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
  
  #Load data files needed for the analysis
  if(runname=="SST"){
    Params<-readRDS(file="SST_param_ana.rds")
    DataIn<-readRDS(file = "SST_DataIn.rds")
  } else if(runname=="subT"){
    Params<-readRDS(file="SubT_param_ana.rds")
    DataIn<-readRDS(file = "subT_DataIn.rds")
  }
  
  #Thin the samples to the right number 
  #(so as to use the full span of the ensemble even if few samples are used)
  ind_s<-round(seq(1,length(Params$tau2.samples),length.out=Nsamps))
  alpha_samples<-Params$alpha.samples[,ind_s]
  beta_samples<-Params$beta.samples[,ind_s]
  tau2_samples<-Params$tau2.samples[ind_s]
  
  #Get the number of observations
  Nd<-length(dats)
  
  #Build the output structures. Probably not needed here, 
  #mostly as to keep similar to Matlab code
  Output_Preds<-matrix(NaN,nrow=Nd,ncol=3)
  Output_AnLocs<-NA
  Output_PriorMean<-prior_mean
  Output_PriorStd<-prior_std
  
  if(ens_sel==1){Output_PredsEns<-NA} else {Output_PredsEns<-NA} #maybe move to later
  
  #Find the analogs
  #Cylce through the alpha/beta grid cells, find those that feature mean
  #modern TEX obs within the tolerance
  
  #Number of big grids
  N_bg<-dim(matrix(unlist(DataIn[[1]][2]),ncol=2))[1]
  
  #Calculate mean SST across spatial grid cells
  spatialMean<-unlist(lapply(1:N_bg,function(a){mean(unlist(DataIn[[1]][8])[unlist(DataIn[[1]][5])==a])}))
  
  #Identify mean values within the tolerance
  inder_g<- spatialMean >= mean(dats)-search_tol & spatialMean <= mean(dats)+search_tol
  if(sum(inder_g)==0){stop("Your search tolerance is too narrow")}
 
  alpha_samples<-as.data.frame(alpha_samples)[inder_g,]
  beta_samples<-as.data.frame(beta_samples)[inder_g,]
  #Tile tau2 to match
  tau2_samples<-kronecker(matrix(1,dim(alpha_samples)[1],1),tau2_samples)
  #Reshape
  alpha_samples<-matrix(as.matrix(alpha_samples),nrow=1,ncol=dim(alpha_samples)[1]*dim(alpha_samples)[2])
  beta_samples<-matrix(as.matrix(beta_samples),nrow=1,ncol=dim(beta_samples)[1]*dim(beta_samples)[2])
  
  Output_AnLocs<-matrix(unlist(DataIn[[1]][2]),ncol=2)[inder_g,]
  
  #Solve
  
  #Prior mean and inverse covariance matrix
  pmu<-kronecker(matrix(1,1,dim(alpha_samples)[2]),matrix(1,Nd,1)*prior_mean)
  pinv_cov<-kronecker(matrix(1,Nd,dim(alpha_samples)[2]),prior_std)^(-2)
  sigmaS<-sqrt(tau2_samples)
  
  #Posterior calculations
  post_mean_num<-pinv_cov*pmu + kronecker(matrix(1,Nd,1),t(sigmaS))^(-2)*kronecker(matrix(1,Nd,1),as.matrix(beta_samples))*(dats- kronecker(matrix(1,Nd,1),as.matrix(alpha_samples)))
  post_mean_den<-pinv_cov+kronecker(matrix(1,Nd,1),as.matrix(beta_samples))^2 * kronecker(matrix(1,Nd,1),t(sigmaS))^(-2)
  post_mean<-post_mean_num / post_mean_den
  post_sig<-sqrt(post_mean_den^(-1))
  Preds<-post_mean+matrix(rnorm(Nd*dim(alpha_samples)[2],0,1),nrow = Nd,ncol=dim(alpha_samples)[2])*post_sig
  
  #Get the percentiles
  Output_Preds<- if(length(dats)==1){quantile(sort(Preds),probs=c(0.05,0.5,0.95))} else if(length(dats)>1){t(apply(apply(Preds,1,sort),2,quantile,probs=c(0.05,0.5,0.95)))}
  
  #Get the ensemble if asked
  if(ens_sel==1){
    Output_PredsEns<-Preds
  }
  
  Out<-list("Predictions"=Output_Preds,"AnLocs"=Output_AnLocs,"PriorMean"=Output_PriorMean,"PriorStd"=Output_PriorStd,"Ensemble"=Output_PredsEns)
  
}