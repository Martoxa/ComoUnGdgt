#' Structural Set Fractional Abundance (FA) Calculation Function
#'
#' Calculate fractional abundances for the structural sets defined in
#' Raberg et al. (2021): (Meth, Cyc, Isom, MI, CI, MC, Full)
#' @param raw A data frame of HPLC-MS peak area values (Nx15+1). Column 1 must
#' have unique identifiers for each sample. GDGT names should be formatted as per
#' the publication.
#' @param set_type Indicates the set to be calculated. Options are:
#' "Meth" (methylation set), "Cyc" (cyclization set), "Isom" (Isomer set),
#' "MI" (Meth-Isom set), "CI" (Cyc-Isom set), "MC" (Meth-Cyc set),
#' and "Full" (Full set).
#' @param data Indicates the type of data being processed. The only accepted
#' data types are peak areas or concentrations ("raw").
#' @references Raberg, J. H., Harning, D. J., Crump, S. E., de Wet, G., Blumm, A., Kopf, S., ... & Sepúlveda, J. (2021). Revised fractional abundances and warm-season temperatures substantially improve brGDGT calibrations in lake sediments. Biogeosciences, 18(12), 3579-3603.
#' @return result: A data frame with additional columns where the indicated
#' FAs for the selected sets are calculated for each brGDGT for each sample.
#' @author Jonathan H. Raberg
#' @note Function structure developed with the aid of Claude AI
#' @export


calculate_FAs_by_set <- function(raw, set_type = "all", data="raw") {

  ##Required packages
  #if (!require(dplyr)){cat("informative error message")}

  ##Check data is concentration or peak area
  if(data!="raw"){stop("Data should be in concentration or peak area format.")}

  df<-format_GDGT_sets(raw)

  # Define the grouping variables for each set type
  grouping_vars <- list(
    "Meth" = c("pos", "cycle"),
    "Cyc" = c("methyl", "pos"),
    "Isom" = c("methyl", "cycle"),
    "MI" = c("cycle"),
    "CI" = c("methyl"),
    "MC" = c("pos"),
    "Full" = character(0)  # empty character vector for no grouping
  )

  # Validate set_type input
  valid_sets <- c(names(grouping_vars), "all")
  if (!set_type %in% valid_sets) {
    stop("set_type must be one of: ", paste(valid_sets, collapse = ", "))
  }

  # Calculate for specific set or all sets
  if (set_type == "all") {
    # Calculate all sets and join them
    result <- df

    for (set_name in names(grouping_vars)) {
      group_vars <- grouping_vars[[set_name]]
      temp_result <- calculate_single_set(df, set_name, group_vars)

      # Extract just the new FA column
      fa_col <- paste0(set_name, "_set_FA")
      temp_fa <- temp_result %>%
        dplyr::select(tidyselect::all_of(c(names(df), fa_col)))

      # Join with main result
      if (set_name == names(grouping_vars)[1]) {
        result <- temp_fa
      } else {
        # Join by all original columns
        join_cols <- names(df)
        result <- result %>%
          dplyr::left_join(temp_fa, by = join_cols)
      }
    }

  } else {
    # Calculate for single set
    group_vars <- grouping_vars[[set_type]]
    result <- calculate_single_set(df, set_type, group_vars)
  }

  return(result)
}
