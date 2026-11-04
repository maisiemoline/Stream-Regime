# ==============================================================================
# File: R/fetch_streamflow.R
# Description: Streamflow Data Ingestion & Harmonization Module with Robust Caching
# Author: Municipal Hydrology & Watershed Science Pair
# ==============================================================================

library(tidyverse)
library(dataRetrieval)
library(httr)
library(jsonlite)
library(lubridate)

#' Fetch Daily Streamflow Data from USGS NWIS with Caching & Retry
#'
#' Retrieves daily mean discharge (parameter code 00060, stat 00003) for one or more
#' USGS streamgage stations with automated per-site query fallback and disk caching.
#'
#' @param site_numbers Character vector of USGS station numbers (e.g. "06696980", "06700000")
#' @param start_date Character string "YYYY-MM-DD"
#' @param end_date Character string "YYYY-MM-DD"
#' @param cache_dir Directory to store/load cached RDS data (default: "data_cache")
#' @return A tidy tibble with standardized columns
#' @export
fetch_usgs_streamflow <- function(site_numbers, 
                                  start_date = "1990-01-01", 
                                  end_date = as.character(Sys.Date()),
                                  cache_dir = "data_cache") {
  
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  
  site_results <- list()
  
  for (site_no in site_numbers) {
    cache_file <- file.path(cache_dir, paste0("usgs_", site_no, "_", start_date, "_", end_date, ".rds"))
    
    if (file.exists(cache_file)) {
      message(sprintf("Loading cached USGS data for site: %s", site_no))
      site_results[[site_no]] <- readRDS(cache_file)
      next
    }
    
    message(sprintf("Fetching USGS streamflow for site: %s", site_no))
    
    # Metadata lookup
    site_meta <- tryCatch({
      readNWISsite(site_no) %>%
        select(site_no, station_nm, dec_lat_va, dec_long_va, drain_area_va)
    }, error = function(e) {
      tibble(site_no = site_no, station_nm = paste("USGS Site", site_no),
             dec_lat_va = NA_real_, dec_long_va = NA_real_, drain_area_va = NA_real_)
    })
    
    # Query with retry
    raw_df <- NULL
    for (attempt in 1:3) {
      raw_df <- tryCatch({
        readNWISdv(
          siteNumbers = site_no,
          parameterCd = "00060",
          startDate = start_date,
          endDate = end_date,
          statCd = "00003"
        )
      }, error = function(e) {
        if (attempt < 3) Sys.sleep(2 * attempt)
        NULL
      })
      if (!is.null(raw_df)) break
    }
    
    if (is.null(raw_df) || nrow(raw_df) == 0) {
      warning("No daily records returned or query failed for USGS site: ", site_no)
      next
    }
    
    clean_site <- raw_df %>%
      renameNWISColumns() %>%
      left_join(site_meta, by = "site_no") %>%
      transmute(
        agency = "USGS",
        site_id = as.character(site_no),
        station_name = if_else(is.na(station_nm), paste("USGS", site_no), station_nm),
        latitude = dec_lat_va,
        longitude = dec_long_va,
        drainage_area_sqmi = drain_area_va,
        date = as.Date(Date),
        doy = yday(date),
        year = year(date),
        month = month(date),
        water_year = if_else(month >= 10, year + 1L, year),
        flow_cfs = as.numeric(Flow),
        qualifier = if ("Flow_cd" %in% names(.)) Flow_cd else NA_character_
      ) %>%
      filter(!is.na(flow_cfs) & flow_cfs >= 0) %>%
      arrange(date)
    
    saveRDS(clean_site, cache_file)
    site_results[[site_no]] <- clean_site
  }
  
  if (length(site_results) == 0) return(tibble())
  return(bind_rows(site_results))
}

#' Fetch Daily Streamflow Data from Colorado DWR OpenCDSS REST API v2 with Caching
#'
#' Retrieves daily surface water streamflow records directly from the State of Colorado's
#' Division of Water Resources (CDSS) REST API for specified station abbreviations.
#'
#' @param abbrevs Character vector of DWR station abbreviations (e.g. "PLAHARCO", "PLAANTCO")
#' @param start_date Character string "YYYY-MM-DD"
#' @param end_date Character string "YYYY-MM-DD"
#' @param cache_dir Directory to store/load cached RDS data (default: "data_cache")
#' @return A tidy tibble with standardized columns matching USGS output
#' @export
fetch_dwr_streamflow <- function(abbrevs, 
                                 start_date = "1990-01-01", 
                                 end_date = as.character(Sys.Date()),
                                 cache_dir = "data_cache") {
  
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  
  results <- list()
  s_fmt <- format(as.Date(start_date), "%m-%d-%Y")
  e_fmt <- format(as.Date(end_date), "%m-%d-%Y")
  
  for (abbrev in abbrevs) {
    cache_file <- file.path(cache_dir, paste0("dwr_", abbrev, "_", start_date, "_", end_date, ".rds"))
    
    if (file.exists(cache_file)) {
      message(sprintf("Loading cached DWR data for site: %s", abbrev))
      results[[abbrev]] <- readRDS(cache_file)
      next
    }
    
    message(sprintf("Fetching Colorado DWR streamflow for site: %s", abbrev))
    
    # 1. Station metadata
    meta_url <- paste0(
      "https://dwr.state.co.us/Rest/GET/api/v2/surfacewater/surfacewaterstations/?format=json&abbrev=", 
      abbrev
    )
    meta_res <- GET(meta_url, user_agent("Hydrology-Regime-Analysis-R"))
    
    st_name <- abbrev
    lat <- NA_real_
    lon <- NA_real_
    
    if (status_code(meta_res) == 200) {
      meta_json <- fromJSON(content(meta_res, as = "text", encoding = "UTF-8"))
      if (!is.null(meta_json$ResultList) && nrow(meta_json$ResultList) > 0) {
        m <- meta_json$ResultList[1, ]
        st_name <- if (!is.null(m$stationName) && !is.na(m$stationName)) m$stationName else abbrev
        lat <- if (!is.null(m$latitude)) as.numeric(m$latitude) else NA_real_
        lon <- if (!is.null(m$longitude)) as.numeric(m$longitude) else NA_real_
      }
    }
    
    # 2. Daily time series data
    ts_url <- paste0(
      "https://dwr.state.co.us/Rest/GET/api/v2/surfacewater/surfacewatertsday/?format=json",
      "&abbrev=", abbrev,
      "&min-measDate=", s_fmt,
      "&max-measDate=", e_fmt
    )
    
    ts_res <- GET(ts_url, user_agent("Hydrology-Regime-Analysis-R"))
    
    if (status_code(ts_res) == 200) {
      ts_json <- fromJSON(content(ts_res, as = "text", encoding = "UTF-8"))
      if (!is.null(ts_json$ResultList) && nrow(ts_json$ResultList) > 0) {
        df_site <- ts_json$ResultList %>%
          transmute(
            agency = "CO DWR",
            site_id = abbrev,
            station_name = st_name,
            latitude = lat,
            longitude = lon,
            drainage_area_sqmi = NA_real_,
            date = as.Date(substr(measDate, 1, 10)),
            doy = yday(date),
            year = year(date),
            month = month(date),
            water_year = if_else(month >= 10, year + 1L, year),
            flow_cfs = as.numeric(value),
            qualifier = if ("flagA" %in% names(.)) as.character(flagA) else NA_character_
          ) %>%
          filter(!is.na(flow_cfs) & flow_cfs >= 0) %>%
          arrange(date)
        
        saveRDS(df_site, cache_file)
        results[[abbrev]] <- df_site
      } else {
        warning("No daily records returned from DWR for: ", abbrev)
      }
    } else {
      warning(sprintf("DWR REST API returned HTTP %d for station: %s", status_code(ts_res), abbrev))
    }
  }
  
  if (length(results) == 0) return(tibble())
  return(bind_rows(results))
}

#' Fetch Streamflow from Any Mix of USGS and DWR Stations
#'
#' Unified interface to pull and standardize daily streamflow data from both USGS
#' and Colorado DWR.
#'
#' @param usgs_sites Character vector of USGS site numbers (or NULL)
#' @param dwr_sites Character vector of DWR station abbreviations (or NULL)
#' @param start_date Character string "YYYY-MM-DD"
#' @param end_date Character string "YYYY-MM-DD"
#' @param cache_dir Directory to store/load cached RDS data (default: "data_cache")
#' @return A unified, harmonized tidy tibble
#' @export
fetch_all_streamflow <- function(usgs_sites = NULL, 
                                 dwr_sites = NULL, 
                                 start_date = "2000-01-01", 
                                 end_date = "2023-12-31",
                                 cache_dir = "data_cache") {
  
  dfs <- list()
  
  if (!is.null(usgs_sites) && length(usgs_sites) > 0) {
    dfs$usgs <- fetch_usgs_streamflow(usgs_sites, start_date, end_date, cache_dir = cache_dir)
  }
  
  if (!is.null(dwr_sites) && length(dwr_sites) > 0) {
    dfs$dwr <- fetch_dwr_streamflow(dwr_sites, start_date, end_date, cache_dir = cache_dir)
  }
  
  combined <- bind_rows(dfs) %>%
    arrange(agency, site_id, date)
  
  return(combined)
}
