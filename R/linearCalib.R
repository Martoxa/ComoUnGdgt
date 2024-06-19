#' Application of linear calibrations for GDGTs
#'
#' Applies published forward models based on linear functions for reconstructing environmental parameters using GDGTs
#' @param gdgt Data frame with the peak areas or fractional abundances of GDGTs. Each row should correspond to a sample and each column to each GDGT with the appropriated name format.
#' @param env Selection of the environmental parameter to reconstruct. Options currently are "Temperature" and "pH"
#' @param calibration Selection of the forward model to be used. Currently options are "Russell" (Russell, 2018), "Martinez" (Martinez-Sosa, 2022), "DeJonge" (De Jonge, 2014), and "Naafs" (Naafs, 2017) for temperature; "Russell" (Russell, 2018), and "DeJonge" (De Jonge, 2014) for pH.
#' @param data Data is presented as either peak areas ("raw") or a vector with the corresponding index values for the calibration ("index").
#' @param complete Logic variable. Indicates whether all brGDGTs in the original formulation are present in the data frame (TRUE, default value) or if some are missing. If TRUE the function will not proceed if any variable is missing, it will produce an error. If FALSE, the function will remove the missing variables from the equation.
#' @param na.ignore Logic variable. If FALSE (default) the function will turn any NA values into 0 for the calculations.
#' @return The function returns the MBT'5Me value for each sample provided.
#' @references Russell, J. M., Hopmans, E. C., Loomis, S. E., Liang, J., & Damsté, J. S. S. (2018). Distributions of 5-and 6-methyl branched glycerol dialkyl glycerol tetraethers (brGDGTs) in East African lake sediment: Effects of temperature, pH, and new lacustrine paleotemperature calibrations. Organic Geochemistry, 117, 56-69.
#' @references Martínez-Sosa, P., Tierney, J. E., Stefanescu, I. C., Crampton-Flood, E. D., Shuman, B. N., & Routson, C. (2021). A global Bayesian temperature calibration for lacustrine brGDGTs. Geochimica et Cosmochimica Acta, 305, 87-105.
#' @references De Jonge, C., Hopmans, E. C., Zell, C. I., Kim, J. H., Schouten, S., & Damsté, J. S. S. (2014). Occurrence and abundance of 6-methyl branched glycerol dialkyl glycerol tetraethers in soils: Implications for palaeoclimate reconstruction. Geochimica et Cosmochimica Acta, 141, 97-112.
#' @references Naafs, B. D. A., Inglis, G. N., Zheng, Y., Amesbury, M. J., Biester, H., Bindler, R., ... & Pancost, R. D. (2017). Introducing global peat-specific temperature and pH calibrations based on brGDGT bacterial lipids. Geochimica et Cosmochimica Acta, 208, 285-301.

linearCalib<-function(gdgt,env,calibration,data="raw",complete=TRUE,na.ignore=FALSE){
  ## Calibrations
  #Russell et al., 2018
  AfricaRussellMAAT<-alist(MBT5=,m=32.42,b=-1.21,(MBT5*m)+b)
  AfricaRussellpH<-alist(CBTp=,m=2.65,b=8.95,(CBTp*m)+b)
  #Martinez-Sosa et al., 2022
  GlobalMartinezMAF<-alist(MBT5=,m=1/0.03,b=0.075/0.03,(MBT5*m)+b)
  #DeJonge et al., 2014
  GlobalDeJongeMAAT<-alist(MBT5=,m=31.45,b=-8.57,(MBT5*m)+b)
  GlobalDeJongepH<-alist(CBTp=,m=1.59,b=7.15,(CBTp*m)+b)
  #Naafs et al., 2017
  GlobalNaafsMAAT<-alist(MBT5=,m=52.18,b=-23.05,(MBT5*m)+b)

  Temperature<-list(Russell=AfricaRussellMAAT,Martinez=GlobalMartinezMAF,DeJonge=GlobalDeJongeMAAT,Naafs=GlobalNaafsMAAT)
  pH<-list(Russell=AfricaRussellpH,DeJonge=GlobalDeJongepH)
  Calibrations<-list(Temperature=Temperature,pH=pH)

  if(data=="raw"){
    if(env=="Temperature"){prediction<-as.function(Calibrations[[env]][[calibration]])(MBT5(gdgt,complete = complete,na.ignore=na.ignore))
    } else if(env=="pH"){prediction<-as.function(Calibrations[[env]][[calibration]])(CBTp(gdgt,complete = complete,na.ignore=na.ignore))
    } else stop("Select one of the available environmental parameters: Temperature or pH")
  } else if(data=="index"){
    if(env=="Temperature"){prediction<-as.function(Calibrations[[env]][[calibration]])(gdgt)
    } else if(env=="pH"){prediction<-as.function(Calibrations[[env]][[calibration]])(gdgt)
    } else stop("Select one of the available environmental parameters: Temperature or pH")
  }
  prediction
}
