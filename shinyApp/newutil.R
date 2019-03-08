library(dplyr)
library(stringr)
library(purrr)


uniq_vals_description <- function(ff1) {
  ## remove NA's
  ff <- ff1[!is.na(ff1)]
  uu <- unique(ff)
  if(class(uu)=="numeric") {
    return(str_c("#unique: ", length(uu)))
  }
  ss <- str_c(uu, collapse = ",")
  flog.info(str_c("uniq_vals_description: ", ss))
  if(!length(ss)) {
    flog.error("nothing in this string")
    flog.error(uu)
    return(NULL)
  }
  if (length(ss) && nchar(ss) > 20) {
    ss <- str_c( str_sub(ss, 1, 16), "...")
  }
  ss
}

update_grouping_factor <- function(nm, f, labels, true_label) {
}

handle_factor <- function(nm, f) {

  u_vals <- unique(f)
  num_unique <- length(u_vals)
  ready <- (num_unique == 2)
  type <- class(u_vals)

  labels <- "FALSE,TRUE"
  true_label <- "TRUE"

  if(type == "numeric") {
    description <- sprintf("Numeric: min %4.3f max %4.3f median %4.3f mean %4.3f",
            min(f, na.rm = TRUE),
            max(f, na.rm = TRUE),
            median(f, na.rm = TRUE),
            mean(f, na.rm = TRUE))
  } else {
    description <- sprintf("character, num values %d", num_unique)
    if(num_unique == 2) {
      labels <- str_c(u_vals, collapse = ",")
      true_label <- u_vals[2]
    }
  }

  uu <- uniq_vals_description(f)
  new_df <- data_frame(name = nm,
                       num_unique = num_unique,
                       unique_values = uu,
                       method_applied = "none",
                       ready = (num_unique == 2),
                       type = type,
                       description = description,
                       labels = labels,
                       true_label = true_label)
}

#
identify_factors <- function(df) {
  rows <- map(colnames(df), ~ handle_factor(., pull(df,.)))
  dplyr::bind_rows(rows)
}
