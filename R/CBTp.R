#' CBT' index
#'
#' Calculates the CBT' index (De Jonge et al., 2014) for brGDGTs presented as either peak areas or fractional abundances.
#' @param gdgt Data frame with the peak areas or fractional abundances of brGDGTs. Each row should correspond to a sample and each column to each brGDGT with the appropriated name format.
#' @param complete Logic variable. Indicates whether all brGDGTs in the original formulation are present in the data frame (TRUE, default value) or if some are missing. If TRUE the function will not proceed if any variable is missing, it will produce an error. If FALSE, the function will remove the missing variables from the equation.
#' @param na.ignore Logic variable. If FALSE (default) the function will turn any NA values into 0 for the calculations.
#' @return The function returns the CBT' value for each sample provided. Any value that blows up (<-10000 or >10000) is turned into NA.
#' @references De Jonge, C., Hopmans, E. C., Zell, C. I., Kim, J. H., Schouten, S., & Damsté, J. S. S. (2014). Occurrence and abundance of 6-methyl branched glycerol dialkyl glycerol tetraethers in soils: Implications for palaeoclimate reconstruction. Geochimica et Cosmochimica Acta, 141, 97-112.

CBTp<-function(gdgt,complete=TRUE,na.ignore=FALSE){
  if(na.ignore==TRUE){gdgt[is.na(gdgt)]<-0}
  formula<-"log10(( gdgt$Ic + gdgt$IIa6 + gdgt$IIb6 + gdgt$IIc6 + gdgt$IIIa6 + gdgt$IIIb6 + gdgt$IIIc6 )/( gdgt$Ia + gdgt$IIa5 + gdgt$IIIa5 ))"
  if(complete==TRUE){if(NA %in% match(correctGs(formula),colnames(gdgt))){stop("Missing variables")}}
  if(complete==FALSE){formula<-partialEq(gdgt,formula)}
  out<-eval(parse(text = formula))
  out[-10000 > out]<-NA
  out[out>10000]<-NA
  out
}
