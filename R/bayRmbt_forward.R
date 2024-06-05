#' BAYMBT forward model
#'
#' Applies the BAYMBT forward model for MBT5Me measured in soils and peats.
#' Predicts MBT5Me based on mean annual air temperature or mean temperatures
#' above zero.
#' @param t A scalar or vector of temperature values (1 x N) or (N x 1)
#' @param Tmodel A string corresponding the model you want to use. Options are
#' "T" = assumes mean annual air temperature (BayMBT); "T0" = assumes mean
#' annual temperature of months above zero (BayMBT0).
#' @param Type A string corresponding to the data type. Options are: "soil" or
#' "lake". Function will use the respective calibration.
#'
#' Example: bayRmbt_forward(t,Tmodel="T0",Type="lake")
#'
#' @return A 1000-member ensemble of MBT'5Me values (N x 1000)
#' @references Dearing Crampton-Flood, E., Tierney, J. E., Peterse, F., Kirkels, F. M.,
#' & Sinninghe Damst?, J. S. (2020). BayMBT: A Bayesian calibration model
#' for branched glycerol dialkyl glycerol tetraethers in soils and peats.
#' Geochimica et Cosmochimica Acta, 268, 142-159.
#'
#' @references Martinez-Sosa, P., Tierney, J. E., Stefanescu, I. C., Dearing
#' Crampton-Flood, E., Shuman, B. N., Routson, C. (2021) A global Bayesian
#' temperature calibration for lacustrine brGDGTs.
#' Geochimica et Cosmochimica Acta, 305, 87-105.
#'
bayrmbt_forward<-function(t,Tmodel,Type){

  #ensure vector
  if(is.vector(t) == FALSE){t<-as.vector(t)}
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

  #define dimensions
  Nobs<-length(t)

  #parse parameters
  alpha<-b_draws_final[,2]
  betaT<-b_draws_final[,1]
  sigmab<-sqrt(tau2_draws_final)

  #vectorized calculation of ensemble
  mbt5me<-matrix(rnorm(Nobs*1000,mean=alpha+t%*%t(as.matrix(betaT)),sd=t(replicate(Nobs,sigmab))),ncol = 1000,byrow = TRUE)
  #any MBT values outside 0 to 1 are forced to be in that range
  mbt5me[mbt5me<0]<-0
  mbt5me[mbt5me>1]<-1

  return(mbt5me)
}
