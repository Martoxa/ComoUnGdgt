#' Community index
#'
#' Calculates the Community Index (De Jonge et al., 2019) for brGDGTs presented as either peak areas or fractional abundances.
#' @param gdgt Data frame with the peak areas or fractional abundances of brGDGTs. Each row should correspond to a sample and each column to each brGDGT with the appropriated name format.
#' @param complete Logic variable. Indicates whether all brGDGTs in the original formulation are present in the data frame (TRUE, default value) or if some are missing. If TRUE the function will not proceed if any variable is missing, it will produce an error. If FALSE, the function will remove the missing variables from the equation.
#' @param na.ignore Logic variable. If FALSE (default) the function will turn any NA values into 0 for the calculations.
#' @return The function returns the CI value for each sample provided.
#' @references De Jonge, C., Radujković, D., Sigurdsson, B. D., Weedon, J. T., Janssens, I., & Peterse, F. (2019). Lipid biomarker temperature proxy responds to abrupt shift in the bacterial community composition in geothermally heated soils. Organic Geochemistry, 137, 103897.
#' @export

ComIdx<-function(gdgt,complete=TRUE,na.ignore=FALSE){
  if(na.ignore==TRUE){gdgt[is.na(gdgt)]<-0}
  formula<-"( gdgt$Ia )/( gdgt$Ia + gdgt$IIa5 + gdgt$IIIa5 )"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(gdgt))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(gdgt,formula)}
  eval(parse(text = formula))
}
