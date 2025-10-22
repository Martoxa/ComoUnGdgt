#' Isomer 7Me Ratio
#'
#' Calculates the IR6Me ratio (Wang et al., 2021) for brGDGTs presented as either peak areas or fractional abundances.
#' @param gdgt Data frame with the peak areas or fractional abundances of brGDGTs, extended version. Each row should correspond to a sample and each column to each brGDGT with the appropriated name format.
#' @param complete Logic variable. Indicates whether all brGDGTs in the original formulation are present in the data frame (TRUE, default value) or if some are missing. If TRUE the function will not proceed if any variable is missing, it will produce an error. If FALSE, the function will remove the missing variables from the equation.
#' @param na.ignore Logic variable. If FALSE (default) the function will turn any NA values into 0 for the calculations.
#' @return The function returns the IR7Me value for each sample provided.
#' @references Wang, H., Liu, W., He, Y., Zhou, A., Zhao, H., Liu, H., ... & Liu, Z. (2021). Salinity-controlled isomerization of lacustrine brGDGTs impacts the associated MBT5ME'terrestrial temperature index. Geochimica et Cosmochimica Acta, 305, 33-48.
#' @export

IsoRatsevn<-function(gdgt,complete=TRUE,na.ignore=FALSE){
  if(na.ignore==TRUE){gdgt[is.na(gdgt)]<-0}
  formula<-"( gdgt$IIa7 + gdgt$IIIa7 )/( gdgt$IIa5 + gdgt$IIIa5 + gdgt$IIa6 + gdgt$IIIa6 + gdgt$IIa7 + gdgt$IIIa7 )"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(gdgt))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(gdgt,formula)}
  eval(parse(text = formula))
}
