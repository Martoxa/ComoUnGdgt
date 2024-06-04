#' BaySPR TEX86
#'
#' Predict T from TEX86 using the standard model. The assumption in this case is that oceanographic conditions are sufficiently similar to today so that it is reasonable to use the spatial distribution of regression parameters as fitted from the core tops.
#' @param gdgt A matrix of HPLC-MS peak area values (15xN) or (Nx15), or a
#' scalar vector of MBT'5Me values (1xN) or (Nx1).
#' @param data A string corresponding the type of data to run in the model.
#' Options are "raw" or "tex", which use peak areas or fractional abundances,
#' or TEX86 values, respectively.
#' @param lon Longitude from -180 to 180
#' @param lat Latitude from -90 to 90
#' @param prior_std Prior standard deviation
#' @param runname Specify which model to use. Enter either "SST" or "subT"
#' @param varargin Used to set the number of posterior draws of parameters to
#' mix across, and whether to save the 5th/50th/95th percentiles of the
#' predictions, or the whole ensemble as well. If left empty, then 1000 draws
#' are used and the ensemble of predictions is not saved. If only one argument,
#' it is the number of draws (cannot exceed 15000), and the ensemble is not
#' saved. Nothe that samples are thinned to use the full span of the ensemble.
#' If two arguments, the first gives the number of draws, while the second is
#' an indicator: 0 (default) or 1, will save the 5th/50th/95th percentiles
#' only, or the full ensemble, respectively. In all the above, the total
#' number of ensemble members is capped by the number of draws in the model
#' output.
#' @param complete Determines how is TEX86 calculated. TRUE calculates exactly
#' as per Schouten, et al., 2002 (default). FALSE ommits missing compounds. A
#' warning is generated and the ommitted compounds are shown.
#' prior mean is set as the mean over all instrumental SST observations
#' within max_dist (distance is chordal), or the closest min_num points:
#' whichever has the larger number of observations.
#' To use the closest K points, set min_num=K and max_dist=0.
#' to use all points within L km, set min_num=1 (ensure at least one) and
#' max_dist=L. Results from the paper are with min_num=1; max_dist=500;
#'
#' @return Predictions - Nd by 3 array, where Nd=length(dats), giving the 5/50/95
#' percentiles of the predictions.
#' @return SiteLoc   - Location of dats, as the inputted [lon, lat];
#' @return GridLoc 	- Location of the grid centroid closest to dats, used to pull the alpha and beta
#' samples.
#' @return PriorMean - The prior mean as calculated above.
#' @return PriorStd  - The prior std as input.
#' @return PredsEns  - Nd by Nsamps array of predictions. Only include if
#' ens_sel==1.
#' @references Tierney, J. E., & Tingley, M. P. (2014). A Bayesian, spatially-varying calibration model for the TEX86 proxy. Geochimica et Cosmochimica Acta, 127, 83-106.
#' @references Schouten, S., Hopmans, E. C., Schefuß, E., & Damste, J. S. S. (2002). Distributional variations in marine crenarchaeotal membrane lipids: a new tool for reconstructing ancient sea water temperatures?. Earth and Planetary Science Letters, 204(1-2), 265-274.


bayspR_tex<-function(gdgt,lon,lat,prior_std,runname,varargin=1000,data="raw",complete=TRUE,na.ignore=FALSE){

  #Deal with the input parameters and decide how to handle the input data

  if(data=="raw"){dats<-TEX86(gdgt,complete=complete,na.ignore = na.ignore)
  } else if (data=="tex"){dats<-gdgt}

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
    Params<-SST_param_std
    Obs<-obsSST
  } else if(runname=="subT"){
    Params<-subT_param_std
    Obs<-obssubT
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
