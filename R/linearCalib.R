#' Application of calibrations for GDGTs
#'
#' Applies published forward models based on linear functions (Eq. 1 to 15) for reconstructing temperature using GDGTs
#' @param gdgt Data frame with the peak areas or fractional abundances of GDGTs. Each row should correspond to a sample and each column to each GDGT with the appropriated name format.
#' @param depth Water depth in meters for the Stefanescu temperature calibration (Eq. 10). 1 by default. Depth must be a vector with equal length as the number of samples provided.
#' @param calibration Selection of the forward model to be used. Current options are presented in Suplementary Table 1 of the (Terrestrial GDGT review).
#' @param data Data is presented as peak areas ("raw") or fractional abundances ("frac"). Default is peak areas.
#' @param complete Logic variable. Indicates whether all brGDGTs in the original formulation are present in the data frame (TRUE, default value) or if some are missing. If TRUE the function will not proceed if any variable is missing, it will produce an error. If FALSE, the function will remove the missing variables from the equation.
#' @param na.ignore Logic variable. If FALSE (default) the function will turn any NA values into 0 for the calculations.
#' @return The function returns the MBT'5Me value for each sample provided.
#' @export

linearCalib<-function(gdgt,depth=NA,calibration,data="raw",complete=TRUE,na.ignore=FALSE){
  Calimx<-matrix(c(-8.57,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,31.45,0,0,
                   5.05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,14.86,0,
                   7.17,17.1,25.9,34.4,-28.6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   5.58,17.91,0,0,-18.77,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   -23.05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,52.18,0,0,
                   12.22,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,18.79,0,
                   -1.21,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32.42,0,0,
                   -1.353,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,29.412,0,0,
                   -1.82,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,56.06,0,0,
                   -3.846,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,38.462,0,2.5,
                   -2.19,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,31.19,0,0,
                   4.81,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15.64,0,0,
                   23.81,0,68.8,0,-24.7,0,-41.91,51.59,0,0,-31.02,0,0,0,0,0,0,0,0,
                   7.11,0,67.66,0,-13.54,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   5.19,16.22,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    ),ncol = 19,nrow=15,byrow = TRUE) #This is a matrix with the parameters for all the equations allowed here.

  if(data=="frac"){gdgt<-gdgt} else if(data=="raw") {gdgt<-FracA(gdgt,group="br",coerce=TRUE)} else {print("Data must be peak area or fractional abundance.")}

  if(calibration!=10){
    gdata<-matrix(c(rep(1,dim(gdgt)[1]),gdgt$Ia,gdgt$Ib,gdgt$Ic,gdgt$IIa5,gdgt$IIa6,gdgt$IIb5,gdgt$IIb6,gdgt$IIc5,gdgt$IIc6,gdgt$IIIa5,gdgt$IIIa6,gdgt$IIIb5,gdgt$IIIb6,gdgt$IIIc5,gdgt$IIIc6,MBT5(gdgt),log10((gdgt$Ia+gdgt$Ib+gdgt$Ic+gdgt$IIa6+gdgt$IIIa6)/(gdgt$Ic+gdgt$IIa5+gdgt$IIc5+gdgt$IIIa5+gdgt$IIIa6))),ncol=18)
                  out<-rowSums(t(apply(gdata,1,function(x) Calimx[calibration,1:18]*x)),na.rm = TRUE)
                  out
  } else if (length(depth)==dim(gdgt)[1]){
    gdata<-matrix(c(rep(1,dim(gdgt)[1]),gdgt$Ia,gdgt$Ib,gdgt$Ic,gdgt$IIa5,gdgt$IIa6,gdgt$IIb5,gdgt$IIb6,gdgt$IIc5,gdgt$IIc6,gdgt$IIIa5,gdgt$IIIa6,gdgt$IIIb5,gdgt$IIIb6,gdgt$IIIc5,gdgt$IIIc6,MBT5(gdgt),log10((gdgt$Ia+gdgt$Ib+gdgt$Ic+gdgt$IIa6+gdgt$IIIa6)/(gdgt$Ic+gdgt$IIa5+gdgt$IIc5+gdgt$IIIa5+gdgt$IIIa6)),log(depth)),ncol=19)
                  out<-rowSums(t(apply(gdata,1,function(x) Calimx[calibration,]*x)),na.rm = TRUE)
                  out
  } else {print("Depth must be a vector with water depth for each sample.")}
}
