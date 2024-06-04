#' BaySPR TEX86 (Analogue)
#'
#' The assumption in this case is that oceanographic conditions are sufficiently dissimilar to the modern configuration (e.g. due to tectonic changes) so as to preclude the use of the regression parameters that correspond to the current location of the study site (see Tierney and Tingley, 2014, GCA for details).
#' @param gdgt A matrix of HPLC-MS peak area values (15xN) or (Nx15), or a
#' scalar vector of MBT'5Me values (1xN) or (Nx1).
#' @param prior_mean Prior mean on temperature. Assumed same for all.
#' @param prior_std Prior standard deviation for temperature. Assumed same for all.
#' @param search_tol Tolerance for finding analog locations. Comparison is between the mean of dats and the mean TEX value within each large gridcell.
#' @param runname Enter 'SST' for sea-surface temperature or 'subT' for subsurface T.
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
#' @param data A string corresponding the type of data to run in the model.
#' Options are "raw" or "tex", which use peak areas or fractional abundances,
#' or TEX86 values, respectively.
#' @param complete Determines how is TEX86 calculated. TRUE calculates exactly
#' as per Schouten, et al., 2002 (default). FALSE ommits missing compounds. A
#' warning is generated and the ommitted compounds are shown.
#' prior mean is set as the mean over all instrumental SST observations
#' within max_dist (distance is chordal), or the closest min_num points:
#' whichever has the larger number of observations.
#' To use the closest K points, set min_num=K and max_dist=0.
#' to use all points within L km, set min_num=1 (ensure at least one) and
#' max_dist=L. Results from the paper are with min_num=1; max_dist=500;
#' @return Predictions      - Nd by 3 array, where Nd=length(dats), giving the 5/50/95
#' percentiles of the predictions.
#' @return AnLocs=- Centroids of the grid cells selected as analogs, as [lon,
#' lat] pairs
#' @return PriorMean - The prior mean as calculated above.
#' @return PriorStd  - The prior std as input.
#' @return PredsEns  - Nd by No. analog locatioins by Nsamps array of predictions. Only included if
#' ens_sel==1. Note that the second dimension corresponds to the locations in .AnLocs.
#' @references Tierney, J. E., & Tingley, M. P. (2014). A Bayesian, spatially-varying calibration model for the TEX86 proxy. Geochimica et Cosmochimica Acta, 127, 83-106.



bayspR_tex_analog<-function(gdgt,prior_mean,prior_std,search_tol,runname,varargin=1000,data="raw",complete=TRUE){
  #Deal with the input parameters and decide how to handle the input data
  if(data=="raw"){dats<-TEX86(gdgt,complete=complete)
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
