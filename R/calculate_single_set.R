#' Calculate FAs by sets
#'
#' Calculates the FAs of brGDGTs acording to the sets as defined by
#' Raberg et al., (2022).Data must be formatted with structural information
#' of the GDGTs using the format_GDGT_sets function.
#' @param data A data frame of HPLC-MS peak area values (Nx15+1). Column 1 must
#' have unique identifiers for each sample. GDGT names should be formatted as per
#' the publication.
#' @param set_name Indicates the set to be calculated. Options are:
#' "Meth" (methylation set), "Cyc" (cyclization set), "Isom" (Isomer set),
#' "MI" (Meth-Isom set), "CI" (Cyc-Isom set), "MC" (Meth-Cyc set),
#' and "Full" (Full set).
#' @param group_vars Indicates the structural variables that will be used to
#' select the brGDGTs for each set. Options are "methyl", "cycle", and "pos",
#' following the formatting structure of the format_GDGT_sets function.
#' @references Raberg, J. H., Harning, D. J., Crump, S. E., de Wet, G., Blumm, A., Kopf, S., ... & Sepúlveda, J. (2021). Revised fractional abundances and warm-season temperatures substantially improve brGDGT calibrations in lake sediments. Biogeosciences, 18(12), 3579-3603.
#' @return result: A data frame with an additional column where the indicated
#' FA for a set is calculated for each brGDGT for each sample.
#' @author Jonathan H. Raberg
#' @note Function structure developed with the aid of Claude AI
#' @export

# Function to calculate FAs for a single set
calculate_single_set <- function(data, set_name, group_vars) {

  #Required packages
  #if (!require(dplyr)){cat("informative error message")}

  # Create FA column name
  fa_col <- paste0(set_name, "_set_FA")

  # Handle the grouping
  if (length(group_vars) == 0) {
    # For Full set (no grouping)
    result <- data %>%
      dplyr::group_by(sample_id, .add = TRUE) %>%
      dplyr::mutate(!!fa_col := amount / sum(amount)) %>%
      dplyr::ungroup()
  } else {
    # For sets with grouping variables
    result <- data %>%
      dplyr::group_by(dplyr::across(tidyselect::all_of(c(group_vars, "sample_id"))), .add = TRUE) %>%
      dplyr::mutate(!!fa_col := amount / sum(amount)) %>%
      dplyr::ungroup()
  }

  return(result)
}
