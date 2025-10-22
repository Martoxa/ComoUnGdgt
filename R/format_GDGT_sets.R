#' Format GDGT data to calculate FAs by sets
#'
#' Formats a dataset formatted as a table where columns are GDGT concentrations
#' or peak areas and rows samples, with a first column containing unique ids,
#' and transforms the GDGT columns to long format with structural information.
#' @param raw A matrix of HPLC-MS peak area values (Nx15+1). Column 1 must have
#' unique identifiers for each sample. GDGT names should be formatted as per
#' the publication.
#' @references Raberg, J. H., Harning, D. J., Crump, S. E., de Wet, G., Blumm, A., Kopf, S., ... & Sepúlveda, J. (2021). Revised fractional abundances and warm-season temperatures substantially improve brGDGT calibrations in lake sediments. Biogeosciences, 18(12), 3579-3603.
#' @return output$prior_mean: User choice of prior mean
#' @return raw_data: A data frame with N*15 rows, and 5 additional columns
#' to the original dataset. Columns are: Compound (brGDGT name),
#' Amount (concentration or peak area), Methyl (I to III, depending on the
#' number of methyl groups), Cycle (a to c depending on the number of cyclo
#' pentanes), and Pos (5 or 6 depending on the position of the methyl groups)
#' @author Jonathan H. Raberg
#' @note Function structure developed with the aid of Claude AI
#' @note Here we're assigning Ia, Ib, and Ic to be "5-methyl" compounds here
#' in accordance with Raberg et al. (2021)
#' @export

format_GDGT_sets <- function(raw) {

  ##Required packages
  #if (!require(dplyr)){cat("informative error message")}
  #if (!require(tidyr)){cat("informative error message")}
  #if (!require(stringr)){cat("informative error message")}

  ##Check data type
  raw_data<-raw |>
    # parse compound names
    tidyr::pivot_longer(tidyr::matches("I+[a-c][5-6]?"), names_to = "compound", values_to = "amount") |>
    dplyr::mutate(methyl = dplyr::case_when(stringr::str_detect(compound, "III") ~ "III",
                              stringr::str_detect(compound, "II") ~ "II",
                              stringr::str_detect(compound, "I") ~ "I"),
           cycle = dplyr::case_when(stringr::str_detect(compound, "a") ~ "a",
                             stringr::str_detect(compound, "b") ~ "b",
                             stringr::str_detect(compound, "c") ~ "c"),
           pos = dplyr::if_else(stringr::str_detect(compound, "6"), "6", "5")
    )|>
    dplyr::rename_with(.cols=1,~"sample_id") #First column in dataset must be sample ID or unique identifier
  raw_data
}



