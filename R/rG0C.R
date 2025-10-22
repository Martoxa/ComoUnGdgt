#' Ratio of GDGT0 to Cren
#'
#' Calculates the ratio of GDGT0 to Crenarchaeol (Blaga, 2009) fractional abundances of isoGDGTs. Each row should correspond to a sample and each column to each isoGDGT with the appropriated name format.
#' @param gdgt Data frame with the peak areas or fractional abundances of isoGDGTs. Each row should correspond to a sample and each column to each brGDGT with the appropriated name format.
#' @param complete Logic variable. Indicates whether all isoGDGTs in the original formulation are present in the data frame (TRUE, default value) or if some are missing. If TRUE the function will not proceed if any variable is missing, it will produce an error. If FALSE, the function will remove the missing variables from the equation.
#' @param na.ignore Logic variable. If FALSE (default) the function will turn any NA values into 0 for the calculations.
#' @return The function returns the GDGT0/Cren value for each sample provided.
#' @references Blaga, C. I., Reichart, G. J., Heiri, O., & Sinninghe Damsté, J. S. (2009). Tetraether membrane lipid distributions in water-column particulate matter and sediments: a study of 47 European lakes along a north–south transect. Journal of Paleolimnology, 41, 523-540.
#' @export

rG0C<-function(gdgt,complete=TRUE,na.ignore=FALSE){
  if(na.ignore==TRUE){gdgt[is.na(gdgt)]<-0}
  formula<-"( gdgt$GDGT0 /( gdgt$GDGT0 + gdgt$Cren ) )*100"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(gdgt))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(gdgt,formula)}
  eval(parse(text = formula))
}
