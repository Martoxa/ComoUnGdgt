#' Degrees of Cyclization' index
#'
#' Calculates the DC' index (De Jonge, 2021) for brGDGTs presented as either peak areas or fractional abundances.
#' @param gdgt Data frame with the peak areas or fractional abundances of brGDGTs. Each row should correspond to a sample and each column to each brGDGT with the appropriated name format.
#' @param complete Logic variable. Indicates whether all brGDGTs in the original formulation are present in the data frame (TRUE, default value) or if some are missing. If TRUE the function will not proceed if any variable is missing, it will produce an error. If FALSE, the function will remove the missing variables from the equation.
#' @param na.ignore Logic variable. If FALSE (default) the function will turn any NA values into 0 for the calculations.
#' @return The function returns the DC' value for each sample provided. Any value that blows up (<-10000 or >10000) is turned into NA.
#' @references De Jonge, C., Kuramae, E. E., Radujković, D., Weedon, J. T., Janssens, I. A., & Peterse, F. (2021). The influence of soil chemistry on branched tetraether lipids in mid-and high latitude soils: Implications for brGDGT-based paleothermometry. Geochimica et Cosmochimica Acta, 310, 95-112.

DegCycp<-function(gdgt,complete=TRUE,na.ignore=FALSE){
  if(na.ignore==TRUE){gdgt[is.na(gdgt)]<-0}
  formula<-"( gdgt$Ib + gdgt$IIb5 + gdgt$IIb6 )/( gdgt$Ia + gdgt$IIa5 + gdgt$IIa6 + gdgt$Ib + gdgt$IIb5 + gdgt$IIb6 )"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(gdgt))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(gdgt,formula)}
  out<-eval(parse(text = formula))
  out[-10000 > out]<-NA
  out[out>10000]<-NA
  out
}
