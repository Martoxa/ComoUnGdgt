#' BAYMBT prediction model
#'
#' Predicts Mean annual air temperature or mean temperatures above zero from MBT5Me measured in soils and peats.
#' @param gdgt A matrix of HPLC-MS peak area values (15xN) or (Nx15), or a
#' scalar vector of MBT'5Me values (1xN) or (Nx1).
#' @param prior_mean A scalar prior mean value of T in degrees C.
#' @param prior_std A scalar prior standard deviation value of T in degrees C.
#' @param Tmodel A string corresponding the temperature model you want to
#' calculate. Options are: "T" or "T0", they calculate the mean annual air
#' temperature (BayMBT) or the mean annual temperatures above zero (BayMBT0),
#' respectively.
#' @param Type A string corresponding to the data type. Options are "soil" or
#' "lake". The function uses the respective calibration.
#' @param data A string corresponding the type of data to run in the model.
#' Options are "raw" or "mbt5", which use peak areas or fractional abundances,
#' or MBT'5Me values, respectively.
#' @param complete Logic variable. Indicates whether all 15 brGDGTs are present in the data frame (TRUE, default value) or if some are missing. If TRUE the function will not proceed if any variable is missing, it will produce an error. If FALSE, the function will remove the missing variables from the equation.
#' @param na.ignore Logic variable. If FALSE (default) the function will turn any NA values into 0 for the calculations.
#' Example: bayRmbt_predict(mbt5me,10,10,Tmodel="T0",Type="lake",data="area",complete=TRUE)
#'
#' @references Dearing Crampton-Flood, E., Tierney, J. E., Peterse, F., Kirkels, F. M.,
#' & Sinninghe Damsté, J. S. (2020). BayMBT: A Bayesian calibration model
#' for branched glycerol dialkyl glycerol tetraethers in soils and peats.
#' Geochimica et Cosmochimica Acta, 268, 142-159.
#'
#' @references Martinez-Sosa, P., Tierney, J. E., Stefanescu, I. C., Dearing
#' Crampton-Flood, E., Shuman, B. N., Routson, C. (2021) A global Bayesian
#' temperature calibration for lacustrine brGDGTs.
#' Geochimica et Cosmochimica Acta, 305, 87-105.
#'
#' @return output$prior_mean: User choice of prior mean
#' @return output$prior_std: User choice of prior standard deviation
#' @return output$T: 2.5%, 50%, and 97.5% confidence intervals of posterior SST
#' @return output$ens: full ensemble of posterior SST (N x 1000);

bayrmbt_predict<-function(gdgt,prior_mean,prior_std,Tmodel,Type,data="raw",complete=TRUE,na.ignore=FALSE){
  #adjust data type
  if(data=="raw"){mbt5me<-MBT5(gdgt,complete=complete, na.ignore = na.ignore)
  } else if (data=="mbt5"){mbt5me<-gdgt}
  #ensure vector
  if(is.vector(mbt5me) == FALSE){mbt5me<-as.vector(mbt5me)}
  #load appropriate model
  if(Type == "soil" & Tmodel == "T"){
    params<-baymbt_params_soil
  } else if(Type == "soil" & Tmodel == "T0"){
    params<-baymbt0_params_soil
  } else if(Type == "lake" & Tmodel == "T"){
    params<-baymbt_params_lake
  } else if(Type == "lake" & Tmodel == "T0"){
    params<-baymbt0_params_lake
  } else print("Type or TModel not recognized")

  b_draws_final<-params[[1]]
  tau2_draws_final<-params[[2]]

  #get dimensions of time series and draws
  nd<-length(mbt5me)
  n_draws<-length(tau2_draws_final)

  #parse parameters
  alpha<-b_draws_final[,2]
  betaT<-b_draws_final[,1]
  sigmab<-sqrt(tau2_draws_final)

  #prior mean and inverse covariance matrix
  pmu<- replicate(n_draws,rep(1,nd)*prior_mean)
  pinv_cov<-(replicate(n_draws,rep(prior_std,nd)))^(-2)

  #posterior calculation
  post_mean_num<- pinv_cov*pmu+(t(replicate(nd,sigmab)))^(-2)*t(replicate(nd,betaT))*(mbt5me-t(replicate(nd,alpha)))
  post_mean_den<- pinv_cov+(t(replicate(nd,betaT)))^2*(t(replicate(nd,sigmab)))^(-2)
  post_mean<- post_mean_num/post_mean_den
  post_sig<-sqrt(post_mean_den^(-1))
  output.ens<-post_mean+matrix(rnorm(nd*n_draws),nrow = nd,ncol = n_draws)*post_sig
  output.prior_mean<-prior_mean
  output.prior_std<-prior_std

  #if using BayMBT0, truncate at T<0
  if(Tmodel == "T0"){
    output.ens[output.ens<0]<-0
  }

  output.T<-t(apply(output.ens,1,FUN = quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975)))
  output<-list(output.ens,output.prior_mean,output.prior_std,output.T)
  names(output)<-c("ens","prior_mean","prior_std","T")
  return(output)

}
