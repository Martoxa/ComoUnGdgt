#' Percentage of GDGT-2
#'
#' Calculates the percentage of GDGT-2 to other isoGDGTs (Sinninghe Damsté, 2012) for isoGDGTs presented as either peak areas or fractional abundances.
#' @param gdgt Data frame with the peak areas or fractional abundances of isoGDGTs. Each row should correspond to a sample and each column to each isoGDGT with the appropriated name format.
#' @param complete Logic variable. Indicates whether all isoGDGTs in the original formulation are present in the data frame (TRUE, default value) or if some are missing. If TRUE the function will not proceed if any variable is missing, it will produce an error. If FALSE, the function will remove the missing variables from the equation.
#' @param na.ignore Logic variable. If FALSE (default) the function will turn any NA values into 0 for the calculations.
#' @return The function returns the %GDGT2 value for each sample provided.
#' @references Sinninghe Damsté, J. S. S., Ossebaar, J., Schouten, S., & Verschuren, D. (2012). Distribution of tetraether lipids in the 25-ka sedimentary record of Lake Challa: extracting reliable TEX86 and MBT/CBT palaeotemperatures from an equatorial African lake. Quaternary Science Reviews, 50, 43-54.
#' @export

pGDGT2<-function(gdgt,complete=TRUE,na.ignore=FALSE){
  if(na.ignore==TRUE){gdgt[is.na(gdgt)]<-0}
  formula<-"( gdgt$GDGT1 + gdgt$GDGT2 + gdgt$GDGT3 )/( gdgt$GDGT1 + gdgt$GDGT2 + gdgt$GDGT3 + gdgt$Cren + gdgt$Creni )"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(gdgt))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(gdgt,formula)}
  eval(parse(text = formula))
}
