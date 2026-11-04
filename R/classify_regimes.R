# ==============================================================================
# File: R/classify_regimes.R
# Description: Streamflow Regime Classification Engine (Excel Demarcation & Dynamic)
# Author: Municipal Hydrology & Watershed Science Pair
# ==============================================================================

library(tidyverse)
library(lubridate)

#' Calculate Day-of-Year (DOY) Hydrograph Envelopes & Percentiles
#'
#' Aggregates multi-year daily streamflow by Day-of-Year (1 to 365/366) and calculates
#' statistical quantiles across historical years. Optional rolling window smoothing
#' eliminates day-to-day weather noise while preserving seasonal inflection points.
#'
#' @param df Tidy streamflow data frame with columns: site_id, date, flow_cfs
#' @param smooth_window Integer width of rolling smoothing window in days (default: 7). Set to 1 for raw.
#' @param p_high Numeric percentile for high-flow envelope curve (default: 0.20, representing 20th percentile)
#' @return Tibble containing DOY hydrograph envelopes for each station
#' @export
calculate_doy_envelopes <- function(df, smooth_window = 7, p_high = 0.20) {
  
  df_clean <- df %>%
    filter(!is.na(flow_cfs) & flow_cfs >= 0) %>%
    mutate(doy = yday(date))
  
  # Compute DOY statistics across all years
  doy_stats <- df_clean %>%
    group_by(agency, site_id, station_name, doy) %>%
    summarise(
      n_years = n_distinct(year(date)),
      n_obs = n(),
      q_min = min(flow_cfs, na.rm = TRUE),
      q_p10 = quantile(flow_cfs, 0.10, na.rm = TRUE),
      q_p20 = quantile(flow_cfs, p_high, na.rm = TRUE),
      q_p50 = median(flow_cfs, na.rm = TRUE),
      q_mean = mean(flow_cfs, na.rm = TRUE),
      q_p80 = quantile(flow_cfs, 0.80, na.rm = TRUE),
      q_p90 = quantile(flow_cfs, 0.90, na.rm = TRUE),
      q_max = max(flow_cfs, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(site_id, doy)
  
  # Apply circular / centered rolling smoothing if requested
  if (smooth_window > 1) {
    k <- smooth_window
    half_k <- floor(k / 2)
    
    doy_stats <- doy_stats %>%
      group_by(site_id) %>%
      group_modify(~ {
        d_sub <- .x %>% arrange(doy)
        n_d <- nrow(d_sub)
        
        smooth_col <- function(vec) {
          # Circular wrap for rolling filter to avoid edge distortion at Jan 1 / Dec 31
          extended <- c(tail(vec, half_k), vec, head(vec, half_k))
          res <- stats::filter(extended, rep(1 / k, k), sides = 2)
          clean_res <- as.numeric(res[(half_k + 1):(half_k + n_d)])
          return(clean_res)
        }
        
        d_sub %>%
          mutate(
            q_min_smooth = smooth_col(q_min),
            q_p10_smooth = smooth_col(q_p10),
            q_p20_smooth = smooth_col(q_p20),
            q_p50_smooth = smooth_col(q_p50),
            q_mean_smooth = smooth_col(q_mean),
            q_p80_smooth = smooth_col(q_p80),
            q_p90_smooth = smooth_col(q_p90),
            q_max_smooth = smooth_col(q_max)
          )
      }) %>%
      ungroup()
  } else {
    doy_stats <- doy_stats %>%
      mutate(
        q_min_smooth = q_min,
        q_p10_smooth = q_p10,
        q_p20_smooth = q_p20,
        q_p50_smooth = q_p50,
        q_mean_smooth = q_mean,
        q_p80_smooth = q_p80,
        q_p90_smooth = q_p90,
        q_max_smooth = q_max
      )
  }
  
  return(doy_stats)
}

#' Classify Streamflow Regimes Using Excel Demarcation Envelope Methodology
#'
#' Implements the Excel quantile-intersection demarcation workflow:
#' 1. Demarcation threshold (T_demarcation) is computed as the p_baseflow (default 85th) percentile of DOY minimums.
#' 2. Runoff season boundaries (d_run_start, d_run_end) are where the maximum envelope curve intersects T_demarcation (or upward seasonal inflection for perennial streams).
#' 3. High flow season boundaries (d_high_start, d_high_end) are where the p_high (default 20th) percentile curve intersects T_demarcation.
#' 4. Days are partitioned into four regimes:
#'    - 'low flow season'
#'    - 'shoulder: low to high'
#'    - 'high flow season'
#'    - 'shoulder: high to low'
#'
#' @param doy_envelopes Output from calculate_doy_envelopes()
#' @param p_baseflow Percentile for demarcation threshold (default: 0.85)
#' @param p_high Percentile for high flow curve (default: 0.20)
#' @param use_smoothed Logical whether to use smoothed envelopes for intersection calculation (default: TRUE)
#' @return A list containing: threshold_summary, seasonal_boundaries, doy_classified
#' @export
classify_climatological_regime <- function(doy_envelopes, 
                                           p_baseflow = 0.85, 
                                           p_high = 0.20,
                                           use_smoothed = TRUE) {
  
  sites <- unique(doy_envelopes$site_id)
  
  boundaries_list <- list()
  doy_class_list <- list()
  thresholds_list <- list()
  
  regime_levels <- c(
    "low flow season",
    "shoulder: low to high",
    "high flow season",
    "shoulder: high to low"
  )
  
  for (s in sites) {
    site_doy <- doy_envelopes %>% 
      filter(site_id == s) %>% 
      arrange(doy)
    
    st_name <- unique(site_doy$station_name)[1]
    agency_nm <- unique(site_doy$agency)[1]
    
    # Select working vectors
    q_min_vec <- if (use_smoothed) site_doy$q_min_smooth else site_doy$q_min
    q_p20_vec <- if (use_smoothed) site_doy$q_p20_smooth else site_doy$q_p20
    q_p50_vec <- if (use_smoothed) site_doy$q_p50_smooth else site_doy$q_p50
    q_max_vec <- if (use_smoothed) site_doy$q_max_smooth else site_doy$q_max
    
    # 1. Calculate demarcation threshold (85th %ile of daily minimums)
    t_demarcation <- as.numeric(quantile(q_min_vec, p_baseflow, na.rm = TRUE))
    
    # Peak timing from median hydrograph
    doy_peak <- site_doy$doy[which.max(q_p50_vec)]
    q_peak_med <- max(q_p50_vec, na.rm = TRUE)
    
    # Winter baseline (Jan-Feb and late Dec)
    winter_sub <- site_doy %>% filter(doy %in% c(1:60, 320:365))
    winter_base <- if (nrow(winter_sub) > 0) median(winter_sub$q_p50_smooth, na.rm = TRUE) else min(q_p50_vec, na.rm = TRUE)
    
    # 2. Runoff season onset and termination
    if (min(q_max_vec, na.rm = TRUE) >= t_demarcation) {
      # For perennial regulated rivers where maximums exceed baseflow all year,
      # runoff onset/end is when median rises above winter base + 15% amplitude
      t_runoff <- winter_base + 0.15 * (q_peak_med - winter_base)
      rising_days <- site_doy %>% filter(doy <= doy_peak, q_p50_smooth >= t_runoff)
      falling_days <- site_doy %>% filter(doy >= doy_peak, q_p50_smooth >= t_runoff)
      
      d_run_start <- if (nrow(rising_days) > 0) min(rising_days$doy) else 90
      d_run_end   <- if (nrow(falling_days) > 0) max(falling_days$doy) else 290
    } else {
      # Natural snowmelt: exact intersection of Q_max with T_demarcation
      rising_days <- site_doy %>% filter(doy <= doy_peak, q_max_vec >= t_demarcation)
      falling_days <- site_doy %>% filter(doy >= doy_peak, q_max_vec >= t_demarcation)
      
      d_run_start <- if (nrow(rising_days) > 0) min(rising_days$doy) else 90
      d_run_end   <- if (nrow(falling_days) > 0) max(falling_days$doy) else 290
    }
    
    # 3. High flow season window
    # Intersection of Q_p20 with T_demarcation on rising and falling limbs
    high_rising <- site_doy %>% filter(doy <= doy_peak, doy >= d_run_start, q_p20_vec >= t_demarcation)
    high_falling <- site_doy %>% filter(doy >= doy_peak, doy <= d_run_end, q_p20_vec >= t_demarcation)
    
    effective_p_used <- p_high
    if (nrow(high_rising) == 0 || nrow(high_falling) == 0) {
      # Fallback to median envelope (45% peak amplitude) if Q_p20 has high variability
      t_high <- winter_base + 0.45 * (q_peak_med - winter_base)
      high_rising <- site_doy %>% filter(doy <= doy_peak, doy >= d_run_start, q_p50_vec >= t_high)
      high_falling <- site_doy %>% filter(doy >= doy_peak, doy <= d_run_end, q_p50_vec >= t_high)
      effective_p_used <- 0.50
    }
    
    d_high_start <- if (nrow(high_rising) > 0) min(high_rising$doy) else d_run_start + 18
    d_high_end   <- if (nrow(high_falling) > 0) max(high_falling$doy) else d_run_end - 18
    
    # Enforce strict chronological ordering and spacing
    if (d_high_start <= d_run_start) d_high_start <- d_run_start + 12
    if (d_high_end >= d_run_end) d_high_end <- d_run_end - 12
    if (d_high_end <= d_high_start) d_high_end <- d_high_start + 20
    
    # Record threshold metadata
    thresholds_list[[s]] <- tibble(
      agency = agency_nm,
      site_id = s,
      station_name = st_name,
      t_demarcation_cfs = round(t_demarcation, 2),
      p_baseflow = p_baseflow,
      p_high_target = p_high,
      p_high_effective = effective_p_used,
      min_envelope_min_cfs = round(min(q_min_vec, na.rm = TRUE), 2),
      max_envelope_max_cfs = round(max(q_max_vec, na.rm = TRUE), 2),
      median_peak_cfs = round(q_peak_med, 2)
    )
    
    ref_year <- 2023
    doy_to_date_str <- function(d) {
      format(as.Date(d - 1, origin = paste0(ref_year, "-01-01")), "%b %d")
    }
    
    boundaries_list[[s]] <- tibble(
      agency = agency_nm,
      site_id = s,
      station_name = st_name,
      runoff_start_doy = d_run_start,
      runoff_start_date = doy_to_date_str(d_run_start),
      high_flow_start_doy = d_high_start,
      high_flow_start_date = doy_to_date_str(d_high_start),
      high_flow_end_doy = d_high_end,
      high_flow_end_date = doy_to_date_str(d_high_end),
      runoff_end_doy = d_run_end,
      runoff_end_date = doy_to_date_str(d_run_end),
      t_demarcation_cfs = round(t_demarcation, 2)
    )
    
    # Classify each DOY (1 to 366)
    classified_doy <- site_doy %>%
      mutate(
        t_demarcation = t_demarcation,
        flow_regime = case_when(
          doy < d_run_start | doy > d_run_end ~ "low flow season",
          doy >= d_run_start & doy < d_high_start ~ "shoulder: low to high",
          doy >= d_high_start & doy <= d_high_end ~ "high flow season",
          doy > d_high_end & doy <= d_run_end ~ "shoulder: high to low",
          TRUE ~ "low flow season"
        ),
        flow_regime = factor(flow_regime, levels = regime_levels)
      )
    
    doy_class_list[[s]] <- classified_doy
  }
  
  return(list(
    thresholds = bind_rows(thresholds_list),
    boundaries = bind_rows(boundaries_list),
    doy_classified = bind_rows(doy_class_list)
  ))
}

#' Classify Daily Real-Time Hydrologic State
#'
#' Provides an operational daily classification based on actual observed flow magnitude
#' relative to baseline/high flow thresholds and 7-day rate-of-change (rising vs recession).
#'
#' @param df Daily streamflow data frame with site_id, date, flow_cfs
#' @param thresholds Output table from classify_climatological_regime()$thresholds
#' @return Data frame with daily dynamic flow regime classifications
#' @export
classify_daily_dynamic_regime <- function(df, thresholds) {
  
  regime_levels <- c(
    "low flow season",
    "shoulder: low to high",
    "high flow season",
    "shoulder: high to low"
  )
  
  df_dyn <- df %>%
    left_join(thresholds %>% select(site_id, t_demarcation_cfs, median_peak_cfs), by = "site_id") %>%
    group_by(site_id) %>%
    arrange(date) %>%
    mutate(
      dq_7d = (lead(flow_cfs, 3) - lag(flow_cfs, 3)) / 6,
      dq_7d = if_else(is.na(dq_7d), (flow_cfs - lag(flow_cfs, 1)), dq_7d),
      
      high_threshold = pmax(t_demarcation_cfs * 1.5, median_peak_cfs * 0.55),
      
      dynamic_regime = case_when(
        flow_cfs < t_demarcation_cfs ~ "low flow season",
        flow_cfs >= high_threshold ~ "high flow season",
        flow_cfs >= t_demarcation_cfs & flow_cfs < high_threshold & dq_7d >= 0 ~ "shoulder: low to high",
        flow_cfs >= t_demarcation_cfs & flow_cfs < high_threshold & dq_7d < 0 ~ "shoulder: high to low",
        TRUE ~ "low flow season"
      ),
      dynamic_regime = factor(dynamic_regime, levels = regime_levels)
    ) %>%
    ungroup()
  
  return(df_dyn)
}

#' Master Workflow: Diagnose and Classify Flow Regimes for Any Station
#'
#' End-to-end wrapper that ingests raw streamflow data, applies the Excel demarcation
#' envelope method, maps daily regime classifications, and produces tidy analysis tables.
#'
#' @param df Tidy daily streamflow tibble (from fetch_all_streamflow or custom)
#' @param smooth_window Rolling window smoothing width (default: 7 days)
#' @param p_baseflow Percentile for demarcation line (default: 0.85)
#' @param p_high Target percentile for high-flow curve (default: 0.20)
#' @return A comprehensive list containing: daily_classified, doy_envelopes, boundaries, thresholds
#' @export
diagnose_station_regimes <- function(df, 
                                     smooth_window = 7, 
                                     p_baseflow = 0.85, 
                                     p_high = 0.20) {
  
  message("Computing Day-of-Year hydrograph envelopes...")
  envelopes <- calculate_doy_envelopes(df, smooth_window = smooth_window, p_high = p_high)
  
  message("Determining demarcation thresholds and seasonal regime boundaries...")
  regime_res <- classify_climatological_regime(
    envelopes, 
    p_baseflow = p_baseflow, 
    p_high = p_high,
    use_smoothed = TRUE
  )
  
  # Map regime classification back to full daily time series
  message("Mapping regime classifications to daily date records...")
  doy_map <- regime_res$doy_classified %>%
    select(site_id, doy, flow_regime, t_demarcation)
  
  daily_classified <- df %>%
    mutate(doy = yday(date)) %>%
    left_join(doy_map, by = c("site_id", "doy")) %>%
    arrange(agency, site_id, date)
  
  # Also attach dynamic real-time classification for comparison
  daily_classified <- classify_daily_dynamic_regime(daily_classified, regime_res$thresholds)
  
  return(list(
    daily_classified = daily_classified,
    doy_envelopes = regime_res$doy_classified,
    boundaries = regime_res$boundaries,
    thresholds = regime_res$thresholds
  ))
}
