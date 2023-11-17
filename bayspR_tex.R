bayspR_tex<-function(raw,lon,lat,prior_std,runname,varargin=1000,data="area",complete=TRUE,na.ignore=FALSE){
#
# INPUTS:
# raw      - a matrix of HPLC-MS peak area values  (6 x N) or (N x 6) 
#if data="area" (default), or a scalar or vector of TEX86 values (1 x N) or (N x 1). 
# lon       - longitude, from -180 to 180
# lat       - latitude, from -90 to 90
# prior_std - prior standard deviation.
# runname   - specify which model to use. Enter either SST or subT
#
# varargin  - used to set the number of posterior draws of parameters to
# mix across, and whether to save the 5th/50th/95th percentiles of the
# predictions, or the whole ensemble as well.
# - if left empty, then 1000 draws are used and the ensemble of predictions
# is not saved. 
# - if only one argument, it is the number of draws (cannot exceed 15000),
# and the ensemble is not saved. Note that samples are thinned to use the
# full span of the ensmeble. 
# - if two arguments, the first gives the number of draws, while the second
# is an indicator:
# 0: (default) save only the 5th/50th/95th percentiles of the predictions. 
# 1: save the whole ensemble as well. 
#
# in all of the above, the total number of ensemble members is capped by
# the number of draws in the model output.
#
# data: A string corresponding the type of data to run in the model.
# "area" = Peak areas for each isoGDGT (default option)
# "fa" = Fractional abundance of each isoGDGT
# "tex" = Calculated TEX86 values  
# complete: Determines how is TEX86 calculated
# TRUE = Calculated exactly as per Schouten, et al., 2002 (default option)
# FALSE = Omits missing compounds. A warning is generated and the ommitted compounds are shown.
#
#prior mean is set as the mean over all instrumental SST observations
#within max_dist (distance is chordal), or the closest min_num points:
#whichever has the larger number of observations.  
#To use the closest K points, set min_num=K and max_dist=0. 
#to use all points within L km, set min_num=1 (ensure at least one) and
#max_dist=L. Results from the paper are with min_num=1; max_dist=500; 
#
#
#Output structure:
#Output (list)
#Predictions      - Nd by 3 array, where Nd=length(dats), giving the 5/50/95
#percentiles of the predictions. 
#SiteLoc   - Location of dats, as the inputted [lon, lat];
#GridLoc 	- Location of the grid centroid closest to dats, used to pull the alpha and beta
#samples. 
#PriorMean - The prior mean as calculated above. 
#PriorStd  - The prior std as input. 
#PredsEns  - Nd by Nsamps array of predictions. Only include if
#ens_sel==1. 

  #Deal with the input parameters and decide how to handle the input data
  
  if(data=="area"){dats<-TEX86(raw,complete=complete,data=data, na.ignore = na.ignore)
  } else if (data=="fa"){dats<-TEX86(raw,complete=complete,data=data, na.ignore = na.ignore)
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
    Params<-readRDS(file="SST_param_std.rds")
    Obs<-readRDS(file = "obsSST.rds")
  } else if(runname=="subT"){
    Params<-readRDS(file="SubT_param_std.rds")
    Obs<-readRDS(file = "obssubT.rds")
  }

  #Grid spacing is hard-coded here  
  grid_half_space<-10
  #Minimum number of grid cells used to calculate modern prior SST
  min_num<-1
  #Maximum distance to search over for modern prior SST
  max_dist<-500
  
  #Thin the samples to the right number 
  #(so as to use the full span of the ensemble even if few samples are used)
  ind_s<-round(seq(1,length(Params$tau2.samples),length.out=Nsamps))
  alpha_samples_comp<-Params$alpha.samples.comp[,ind_s]
  beta_samples_comp<-Params$beta.samples.comp[,ind_s]
  tau2_samples<-Params$tau2.samples[ind_s]
  
  #Get the number of observations
  Nd<-length(dats)
  
  #Build the output structures. Probably not needed here, 
  #mostly as to keep similar to Matlab code 
  Output_Preds<-matrix(NaN,nrow=Nd,ncol=3)
  Output_SiteLoc<-as.matrix(cbind(lon,lat))
  Output_GridLoc<-matrix(NaN,nrow = 2,ncol=1)
  Output_PriorMean<-matrix(NaN,nrow = 1,ncol=1)
  Output_PriorStd<-prior_std
  
  if(ens_sel==1){Output_PredsEns<-matrix(NaN,nrow = Nd,ncol=Nsamps)} else {Output_PredsEns<-NA}
  
  #Get the prior mean: constant across the calibration cases,
  #as we just average min_num sst obs, on the original 1 degree scale,
  #that are closest to the timeseries, or all that are within a max_dist search radius.
  
  #Get the distance to each of the sst gridcells, on the original resolution 
  dists_prior<-RthChordDistances(Output_SiteLoc,Obs$locs.st.obs)
  #order by distance. Slight deviation from Matlab code, made this in two steps
  sorted_dists_prior<-if(length(lon)==1){sort(dists_prior)} else if(length(lon)>1){apply(dists_prior,2,sort)}
  inds_dists_prior<-if(length(lon)==1){match(sorted_dists_prior,dists_prior)} else if(length(lon)>1){which(sorted_dists_prior==dists_prior,arr.ind = TRUE)[,1]}
  #Get the numbers that are bellow the distance cutoff
  num_below_dist<-which(sorted_dists_prior == tail(c(sorted_dists_prior)[c(sorted_dists_prior)<max_dist],n=1),arr.ind = FALSE)
  
  #If this is larger than the min number, use it to select them:
  if(num_below_dist>min_num){prior_mean<-mean(Obs$st.obs.ave.vec[inds_dists_prior[seq(from=1,to=num_below_dist,by=1)]])
  } else {prior_mean<-mean(Obs$st.obs.ave.vec[inds_dists_prior[seq(from=1,to=min_num,by=1)]])} #else use the min_num
  Output_PriorMean<-prior_mean
  
  #Figure out the alpha and beta series to draw from.
  
  #Just find the index of the Locs.Comp entry that contains the inputted location: 
  inder_g<- which(abs(Params$Locs.Comp[,1]-lon)<=grid_half_space & abs(Params$Locs.Comp[,2]-lat)<= grid_half_space)
  #Extract the alpha, beta series
  alpha_samples_comp<-as.data.frame(alpha_samples_comp)[inder_g,]
  beta_samples_comp<-as.data.frame(beta_samples_comp)[inder_g,]
  Output_GridLoc<-Params$Locs.Comp[inder_g,]
  
  #Solve
  
  #Prior mean and inverse covariance matrix
  pmu<-kronecker(matrix(1,1,Nsamps),matrix(1,Nd,1)*prior_mean)
  pinv_cov<-kronecker(matrix(1,Nd,Nsamps),prior_std)^(-2)
  sigmaS<-sqrt(tau2_samples)
  
  #posterior calculations
  post_mean_num<-pinv_cov*pmu + kronecker(matrix(1,Nd,1),t(sigmaS))^(-2)*kronecker(matrix(1,Nd,1),as.matrix(beta_samples_comp))*(dats- kronecker(matrix(1,Nd,1),as.matrix(alpha_samples_comp)))
  post_mean_den<-pinv_cov+kronecker(matrix(1,Nd,1),as.matrix(beta_samples_comp))^2 * kronecker(matrix(1,Nd,1),t(sigmaS))^(-2)
  post_mean<-post_mean_num / post_mean_den
  post_sig<-sqrt(post_mean_den^(-1))
  Preds<-post_mean+matrix(rnorm(Nd*Nsamps,0,1),nrow = Nd,ncol=Nsamps)*post_sig
  
  #Take the percentiles
  Output_Preds<- if(length(dats)==1){quantile(sort(Preds),probs=c(0.05,0.5,0.95))} else if(length(dats)>1){t(apply(apply(Preds,1,sort),2,quantile,probs=c(0.05,0.5,0.95)))}
  
  #Get the ensemble if asked
  if(ens_sel==1){
    Output_PredsEns<-Preds
  }
  
  Out<-list("Predictions"=Output_Preds,"SiteLoc"=Output_SiteLoc,"GridLoc"=Output_GridLoc,"PriorMean"=Output_PriorMean,"PriorStd"=Output_PriorStd,"Ensemble"=Output_PredsEns)

}