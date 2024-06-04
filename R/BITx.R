#' BIT index
#'
#' Calculates BIT index (Hopmans, 2004) for GDGTs presented as either peak areas or fractional abundances.
#' @param gdgt Data frame with the peak areas or fractional abundances of GDGTs. Each row should correspond to a sample and each column to each GDGT with the appropriated name format.
#' @param complete Logic variable. Indicates whether all GDGTs in the original formulation are present in the data frame (TRUE, default value) or if some are missing. If TRUE the function will not proceed if any variable is missing, it will produce an error. If FALSE, the function will remove the missing variables from the equation.
#' @param na.ignore Logic variable. If FALSE (default) the function will turn any NA values into 0 for the calculations.
#' @return The function returns the BIT value for each sample provided.
#' @references Hopmans, E. C., Weijers, J. W., Schefuß, E., Herfort, L., Damsté, J. S. S., & Schouten, S. (2004). A novel proxy for terrestrial organic matter in sediments based on branched and isoprenoid tetraether lipids. Earth and Planetary Science Letters, 224(1-2), 107-116.

BITx<-function(gdgt,complete=TRUE,out="all",na.ignore=FALSE){
  if(na.ignore==TRUE){gdgt[is.na(gdgt)]<-0}
  formula<-"( gdgt$Ia + gdgt$IIa5 + gdgt$IIa6 + gdgt$IIIa5 + gdgt$IIIa6 )/( gdgt$Ia + gdgt$IIa5 + gdgt$IIa6 + gdgt$IIIa5 + gdgt$IIIa6 + gdgt$Cren )"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(gdgt))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(gdgt,formula)}
  eval(parse(text = formula))
}
