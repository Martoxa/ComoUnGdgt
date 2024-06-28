#' Isomer 6Me Ratio
#'
#' Calculates the IR6Me ratio (De Jonge et al., 2015) for brGDGTs presented as either peak areas or fractional abundances.
#' @param gdgt Data frame with the peak areas or fractional abundances of brGDGTs. Each row should correspond to a sample and each column to each brGDGT with the appropriated name format.
#' @param complete Logic variable. Indicates whether all brGDGTs in the original formulation are present in the data frame (TRUE, default value) or if some are missing. If TRUE the function will not proceed if any variable is missing, it will produce an error. If FALSE, the function will remove the missing variables from the equation.
#' @param na.ignore Logic variable. If FALSE (default) the function will turn any NA values into 0 for the calculations.
#' @return The function returns the IR6Me value for each sample provided.
#' @references De Jonge, C., Stadnitskaia, A., Fedotov, A., & Damsté, J. S. S. (2015). Impact of riverine suspended particulate matter on the branched glycerol dialkyl glycerol tetraether composition of lakes: The outflow of the Selenga River in Lake Baikal (Russia). Organic Geochemistry, 83, 241-252.

IsoRat<-function(gdgt,complete=TRUE,na.ignore=FALSE){
  if(na.ignore==TRUE){gdgt[is.na(gdgt)]<-0}
  formula<-"( gdgt$IIa6 + gdgt$IIb6 + gdgt$IIc6 + gdgt$IIIa6 + gdgt$IIIb6 + gdgt$IIIc6 )/( gdgt$IIa5 + gdgt$IIb5 + gdgt$IIc5 + gdgt$IIIa5 + gdgt$IIIb5 + gdgt$IIIc5 + gdgt$IIa6 + gdgt$IIb6 + gdgt$IIc6 + gdgt$IIIa6 + gdgt$IIIa6 + gdgt$IIIb6 + gdgt$IIIc6 )"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(gdgt))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(gdgt,formula)}
  eval(parse(text = formula))
}
