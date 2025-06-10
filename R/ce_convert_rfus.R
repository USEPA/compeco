#' Converts RFUs to concentration in μg/L
#' 
#' This functions needs input data from a fluorometer of RFU's and then
#' uses a standard curve to make the conversion.  The standard curves are are
#' calculated and stored internally in this package and may be updated.  
#' Fluorometer and spec files needed for the standard curves are detailed in 
#' ADDSOPHERE and will work on both blanked and unblanked spec measurements.  
#' Note that the ADDSOPHERE SOP sepcifies blanked mesaurements.  The unblanked 
#' measurements are include for older standard curve methods but should not be 
#' used going forward. The module, fluorometer, and date arguments will specify 
#' which curve data file to use.  
#' 
#' @param rfus A data frame with sample metadata and RFU.  See ADDEXAMPLEFILE
#'                for an example of the expected file format.
#' @param module One of two options: "ext_chla", "invivo_chla", or "phyco", 
#'                default is "ext_chla".
#' @param fluorometer One of two options: "g04" or "m07", default is 
#'                     "g04". These represent the two Turner Triology's at
#'                     ACESD.  The "g04" one was purchased by the compeco lab, 
#'                     the "m07", has been in use for sometime and currently 
#'                     resides in M07. The default is "g04".
#' @param year Year of the standard curve.  If more than one curve caluclated in
#'             a year add "a", "b", ...  Acceptable values are:  "2021"
#' @param output An optional output path and csv file name for converted values.
#' @param std_check Logical to determine if a solid standard check should be 
#'                  done.  Defaults to TRUE.
#' @param conversion_slope This argument is used to overide the standard curves 
#'                         that are generated internally.  If you are using a 
#'                         different fluorometer than specified above and 
#'                         already know the slope, specify that here as a 
#'                         numeric.
#' @param blank_correction Optional arguement to correct (or not) RFU values  
#'                         with blank RFUs.  Default is TRUE
#' @return returns a tibble with proper metadata, input RFUs and converted
#'          concentrations                     
#' @note While it is possible to mix and match modules and fluorometers, don't 
#'        do this.  Keep the "g04" modules with the "g04" fluorometer 
#'        and vice-versa. 
#' @importFrom readr read_csv write_csv
#' @importFrom dplyr group_by
#' @export
#' @examples
#' examp_chla_data <- readr::read_csv(system.file("extdata/chla_2021_7_14.csv", package = "compeco"), na = c("", "NA", "na"))
#' examp_phyco_data <- readr::read_csv(system.file("extdata/phyco_2021_06_28.csv", package = "compeco"), na = c("", "NA", "na"))
#' surge_chla <- readr::read_csv(system.file("extdata/surge_file.csv", package = "compeco"), na = c("", "NA", "na"))
#' 
#' ce_convert_rfus(rfus = examp_chla_data, module = "ext_chla", year = "2021", fluorometer = "m07")
#' ce_convert_rfus(rfus = examp_phyco_data, module = "phyco", fluorometer = "g04")
#' ce_convert_rfus(rfus = surge_chla, year = "2022", module = "ext_chla", 
#'                 fluorometer = "g04", std_check = FALSE)
ce_convert_rfus <- function(rfus, 
                            module = c("ext_chla", "invivo_chla", "phyco"),
                            year = c("2021", "2022"),
                            fluorometer = c("g04", "m07"),
                            output = NULL,
                            std_check = TRUE,
                            conversion_slope = NULL,
                            blank_correction = TRUE){
  
  # Drop rows with NA in value - nothing to convert so dont need to keep
  # Drop rows with no fluorometer listed - can't convert without that info.
  
  rfus <- filter(rfus, !is.na(value))
  if("fluorometer" %in% names(rfus)){
    rfus <- filter(rfus, !is.na(fluorometer))
    fluorometer <- tolower(unique(rfus$fluorometer))  
  } else {
    fluorometer <- match.arg(fluorometer)
  }
  
  day_na <- any(is.na(rfus$day))
  month_na <- any(is.na(rfus$month))
  year_na <- any(is.na(rfus$year))
  if(any(c(day_na, month_na, year_na))){
    stop("At least one of the dates is missing.  Correct data file and note change.")
  }
  module <- match.arg(module)
  year <- match.arg(year)
  
  if(is.null(conversion_slope)){
    std_curve <- ce_create_std_curve(module, year, fluorometer)
  } else {
    std_curve <- list(std_curve = conversion_slope, solid_std = NA)
  }
  
  
  if("dup" %in% names(rfus)){
    rfus <- dplyr::rename(rfus, field_dups = dup)
  }
  if("rep" %in% names(rfus)){ 
    rfus <- dplyr::rename(rfus, lab_reps = rep)
  }
  if("dups" %in% names(rfus)){
    rfus <- dplyr::rename(rfus, field_dups = dups)
  }
  if("reps" %in% names(rfus)){ 
    rfus <- dplyr::rename(rfus, lab_reps = reps)
  }
  
  # Add missing columns
  names_to_check <- c("waterbody", "site", "depth", "field_dups", "lab_reps", "units", "notes")
  miss_names <- setdiff(names_to_check, names(rfus))
  rfus[miss_names] <- NA
  
  sample_solid_std <- mean(rfus$value[grepl("solid std", rfus$site)], na.rm = TRUE)
  
  # Check solid standard drift
  if(all(is.na(std_curve))){
    perc_diff <- 0
  } else {
    perc_diff <- abs(1-(sample_solid_std/std_curve$solid_std))
  }
  
  perc_diff <- round(perc_diff, 3) * 100
  
  if(is.na(perc_diff) & !is.null(conversion_slope)){
    message("Standard check not possible with conversion_slope argument")
  } else if(is.na(perc_diff)){
    message(paste0("Calculation of percent diff did not succeed and returned ", perc_diff))
  } else if(perc_diff >= 10 & std_check){
    stop(paste0("Sample solid standard is ", perc_diff, "% off from standard curve solid standard."))
  } else if(perc_diff >= 5 & std_check){
    warning(paste0("Sample solid standard is ", perc_diff, "% off from standard curve solid standard."))
  } else {
    message(paste0("Solid standard is off by ", perc_diff, "%."))
  }
  
  # Convert RFU to Concentration
  conc <- ce_convert_to_conc(module, rfus, std_curve, blank_correction)
  
  # Add missing columns
  names_to_check <- c("waterbody", "site", "depth", "field_dups", "lab_reps", "notes", "time")
  miss_names <- setdiff(names_to_check, names(conc))
  conc[miss_names] <- NA
  
  
  
  # Clean up output concentrations
  conc1 <- dplyr::select(conc, date, waterbody, site, depth, field_dups, lab_reps, variable, 
                         units, value = value_rfu_dilute_correct, notes)
  conc1 <- dplyr::mutate(conc1, variable = module, units = "rfu")
  conc2 <- dplyr::select(conc, date, waterbody, site, depth, field_dups, lab_reps, variable, 
                         units, value = value_field_conc_dilute_correct, notes)
  conc2 <- dplyr::mutate(conc2, variable = module, units = "µg/L")
  conc <- dplyr::bind_rows(conc1, conc2)
  conc <- dplyr::mutate(conc, 
                        waterbody = as.character(waterbody),
                        site = as.character(site),
                        depth = as.character(depth),
                        field_dups = as.character(field_dups),
                        lab_reps = as.character(lab_reps),
                        variable = as.character(variable),
                        units = as.character(units),
                        value = as.double(value),
                        notes = as.character(notes))
  
  if(module == "invivo_chla"){
    slope <- NA
  } else if(class(std_curve$std_curve) == "lm"){
    slope <- coef(std_curve$std_curve)
  } else {
    slope <- std_curve$std_curve
  }
  if(!is.null(output)) write_csv(conc, output)
  message(paste0("Std Curve - year: ", year, ", fluorometer: ", fluorometer, 
                 ", slope: ", round(slope, 4)))
  
  # Add perc_diff to Notes
  conc <- conc |>
    mutate(notes = case_when(is.na(notes)~
                              paste0("solid std drift: ", perc_diff),
                            !is.na(notes)~
                              paste0(notes, "; solid std drift: ", perc_diff, "%"),
                            TRUE ~ notes))
  conc
}

#' Create a standard curve
#'
#' This function creates a standard curve from input fluorometer and 
#' spectrophotometer .csv files.  Curve is used to convert RFUs to µg/L.
#' 
#' @param ... Arguments for module, fluorometer and date of input standard curve
#'            csv files.  Passed from \code{\link{ce_convert_rfus}}
#' @keywords internal
ce_create_std_curve <- function(...){
  args <- list(...)
  module <- args[[1]]
  year <- args[[2]]
  fluorom <- args[[3]]
  if(module != "invivo_chla"){
    ffile <- paste0(system.file("extdata", package = "compeco"), "/", 
                    module, "_",  year, "_", fluorom, "_fluorometer.csv")
    sfile <-  paste0(system.file("extdata", package = "compeco"), "/", 
                     module, "_",  year, "_", fluorom, "_spec.xlsx")
    fluoro <- suppressMessages(readr::read_csv(ffile,
                                      na = c("","NA","na")))
  }
  
  if(module == "ext_chla"){
    
    fluoro <- dplyr::rename_all(fluoro, tolower)
    fluoro <- dplyr::group_by(fluoro, standard)
    fluoro <- dplyr::summarize(fluoro, avg_value = mean(value))
    fluoro <- dplyr::ungroup(fluoro)
    blank <- fluoro[fluoro$standard == "blank",]$avg_value
    fluoro <- dplyr::mutate(fluoro, blanked_rfus = 
                              dplyr::case_when(standard != "solid" ~ 
                                                 avg_value - blank,
                                               TRUE ~ avg_value),
                            standard = tolower(standard))
    fluoro <- dplyr::filter(fluoro, .data$standard != "blank")
    sheets <- readxl::excel_sheets(sfile)
    specs <- purrr::map(sheets, function(x) {
      spec <- readxl::read_excel(sfile, sheet = x, skip = 4, 
                                 na = c("","NA","na"))
      spec <- dplyr::mutate(spec, standard = x)})
    specs <- do.call(rbind, specs)
    specs <- dplyr::filter(specs, nm == 750 | nm == 664)
    specs <- dplyr::rename_all(specs, tolower)
    # Blank correction
    blank750 <- specs[tolower(specs$standard) == "blank.sample" & specs$nm == 750,]$a
    blank664 <- specs[tolower(specs$standard) == "blank.sample" & specs$nm == 664,]$a
    
    blanked <- dplyr::near(0, blank750, tol = 0.001) | 
                             dplyr::near(0, blank664, tol = 0.001)
    if(!blanked){
      stop("Find Jeff and Stephen.  The code to deal with un-blanked spec data has not yet been written!")
    } else {
      specs <- tidyr::pivot_wider(specs,id_cols = standard,names_from = nm, 
                                  names_prefix = "nm", values_from = a)
      specs <- dplyr::mutate(specs, corrected_abs = nm664 - nm750)
    }
    
    specs <- dplyr::mutate(specs, conc = (11.4062*(.data$corrected_abs/10))*1000,
                           standard = gsub(".sample", "", 
                                           tolower(.data$standard)))
    specs <- dplyr::filter(specs, .data$standard != "blank")
    fluoro <- dplyr::left_join(fluoro, specs, by = "standard")
    
    std_curve <- lm(conc ~ 0 + blanked_rfus, 
                    data = fluoro[fluoro$standard != "solid",])
  } else if(module == "phyco"){
    fluoro <- dplyr::group_by(fluoro, standard)
    fluoro <- dplyr::summarize(fluoro, avg_value = mean(value), 
                               conc = mean(concentration))
    fluoro <- dplyr::ungroup(fluoro)
    blank <- fluoro[fluoro$standard == "blank",]$avg_value
    fluoro <- dplyr::mutate(fluoro, blanked_rfus = 
                              dplyr::case_when(standard != "solid" ~ 
                                                 avg_value - blank,
                                               TRUE ~ avg_value),
                            standard = tolower(standard))
    fluoro <- dplyr::filter(fluoro, .data$standard != "blank")
    std_curve <- lm(conc ~ 0 + blanked_rfus, 
                    data = fluoro[fluoro$standard != "solid",])
    
  } else if(module == "invivo_chla"){
    return(NA)
  }
  
  list(std_curve = std_curve, solid_std = fluoro$blanked_rfus[fluoro$standard == "solid"])
}


#' Convert RFUS to concentration
#'
#' This function creates takes an compeco input rfu object and standard curve 
#' object and converts the rfus to concentration
#' 
#' @param module which module/variable
#' @param rfus rfu object from \code{ce_convert_rfus}
#' @param std_curve std_curve object from \code{ce_convert_rfus}
#' @keywords internal
ce_convert_to_conc <- function(module = c("ext_chla", "invivo_chla", "phyco"), 
                               rfus, std_curve, blank_correction = TRUE){
  
  module <- match.arg(module)
  rfus <- dplyr::mutate(rfus, 
                        variable = dplyr::case_when(lab_reps == "blank" ~ "blank", 
                                                    TRUE ~ variable))
  
  if(module != "invivo_chla"){
    
    rfus <- dplyr::mutate(rfus, date = lubridate::ymd(paste(year, month, day)))
    blanks <- dplyr::filter(rfus, variable == "blank")
    blanks <- dplyr::group_by(blanks, waterbody, site, date)
    blanks <- dplyr::summarize(blanks , blank_cor = mean(value))
    blanks <- dplyr::ungroup(blanks)
    blanks <- dplyr::select(blanks, waterbody, site, date, blank_cor)
    blanks <- unique(blanks)
    conc <- dplyr::filter(rfus, variable != "blank")
    conc <- dplyr::left_join(conc, blanks)
    if(blank_correction){
      if(any(is.na(conc$blank_cor))){
        conc <- mutate(conc, blank_cor = case_when(is.na(blank_cor) ~
                                                     0,
                                                   TRUE ~ blank_cor))
        message("Blank correction is set to true, but some blanks are NA.  Changing blanks to 0.")
      }
      conc <- dplyr::mutate(conc, value_cor = value - blank_cor)
    } else {
      conc <- dplyr::mutate(conc, value_cor = value)
    }
  } else if(module == "invivo_chla"){
    
    rfus <- dplyr::mutate(rfus, date = lubridate::ymd(paste(year, month, day)))
    blanks <- dplyr::filter(rfus, variable == "blank")
    blanks <- dplyr::group_by(blanks, waterbody, site, field_dups, date)
    blanks <- dplyr::summarize(blanks , blank_cor = mean(value))
    blanks <- dplyr::ungroup(blanks)
    # Do I need all these columns or just waterbody and date.  
    blanks <- dplyr::select(blanks, waterbody, site, date, field_dups, blank_cor)
    conc <- dplyr::left_join(rfus, blanks)
    conc <- dplyr::filter(conc, variable != "blank")
    if(blank_correction){
      if(any(is.na(conc$blank_cor))){
        conc <- mutate(conc, blank_cor = case_when(is.na(blank_cor) ~
                                                     0,
                                                   TRUE ~ blank_cor))
        message("Blank correction is set to true, but some blanks are NA.  Changing blanks to 0.")
      }
      conc <- dplyr::mutate(conc, value_cor = value - blank_cor)
    } else {
      conc <- dplyr::mutate(conc, value_cor = value)
    }
  }
  
  if(all(is.na(std_curve))){
    std_curve_slope <- 1
  } else if(is.na(std_curve$solid_std)){
    std_curve_slope <- std_curve$std_curve
  } else {
    std_curve_slope <- coef(std_curve$std_curve)
  }
  
  conc <- dplyr::mutate(conc, value_cuvette_conc = value_cor * std_curve_slope)
  if(module == "ext_chla"){
    conc <- dplyr::mutate(conc, value_field_conc = value_cuvette_conc * 
                            (0.01/(filter_vol/1000)))
    conc <- dplyr::mutate(conc, value_rfu = value_cor * 
                            (0.01/(filter_vol/1000)))
  } else if(module == "invivo_chla"){
    message("Invivo chla not currently calculating concentration, still just RFUs")
    conc <- dplyr::mutate(conc, value_field_conc = value_cuvette_conc)
    conc <- dplyr::mutate(conc, value_rfu = value_cuvette_conc)
  } else if(module == "phyco"){
    conc <- dplyr::mutate(conc, value_field_conc = value_cuvette_conc * 
                            (0.02/(filter_vol/1000)))
    conc <- dplyr::mutate(conc, value_rfu = value_cor * 
                            (0.02/(filter_vol/1000)))
  }
  idx_site <- is.na(conc$site)
  conc$site[idx_site] <- "missing damnit"
  conc <- dplyr::filter(conc, site != "solid std")

  # Add dilution column if missing the correct for dilutions
  dilution_col <- c("dilution")
  miss_names <- setdiff(dilution_col, names(conc))
  conc[miss_names] <- 1
  
  
  # Fix any dilutions with "x"
  # Seeing problems with fedsdata/turner/chl_2021_07_28.csv
  
  conc <- dplyr::mutate(conc,
                        dilution = as.character(dilution),
                        dilution = 
                          dplyr::case_when(grepl("x", dilution, ignore.case = TRUE) ~ 
                                      gsub("x", "", dilution, 
                                           ignore.case = TRUE),
                                    TRUE ~ dilution),
                        dilution = as.numeric(dilution))
 
  conc <- dplyr::mutate(conc,
                        dilution = dplyr::case_when(is.na(dilution) ~ 1,
                                             TRUE ~ dilution),
                        value_cor_dilute_correct = value_cor * dilution,
                        value_field_conc_dilute_correct = value_field_conc * 
                          dilution,
                        value_rfu_dilute_correct = value_rfu * dilution)
  conc <- dplyr::group_by(conc, waterbody, site, depth, field_dups, lab_reps, variable, day, 
                   month, year)
  
  conc <- dplyr::filter(conc, dilution == max(.data$dilution))
  conc <- dplyr::ungroup(conc)
  conc
}

#' Export Standard Curve
#' 
#' This function exports the linear model for the different standard curves.
#' 
#' @param module One of two options: "ext_chla", "invivo_chla", or "phyco", 
#'                default is "ext_chla".
#' @param fluorometer One of two options: "g04" or "m07", default is 
#'                     "g04". These represent the two Turner Triology's at
#'                     ACESD.  The "g04" one was purchased by the compeco lab, 
#'                     the "m07", has been in use for sometime and currently 
#'                     resides in M07. The default is "g04".
#' @param year Year of the standard curve.  If more than one curve caluclated in
#'             a year add "a", "b", ...  Acceptable values are:  "2021" for
#'             all but ext_chla and g04.  That has "2021" and "2022".
#'             
#' @export
#' @examples 
#' ce_export_std_curve("ext_chla", "g04", "2022")
ce_export_std_curve <- function(module, fluorometer, year){
  std_curve <- ce_create_std_curve(module, year, fluorometer)
  list(std_curve = std_curve$std_curve, 
       data = tibble::tibble(std_curve$std_curve$model), 
       r2 = summary(std_curve$std_curve)$r.squared,
       solid_std = std_curve$solid_std)
}
