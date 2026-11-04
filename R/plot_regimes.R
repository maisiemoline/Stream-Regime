# ==============================================================================
# File: R/plot_regimes.R
# Description: Publication-Grade Visualizations for Streamflow Regime Diagnostics
# Author: Municipal Hydrology & Watershed Science Pair
# ==============================================================================

library(tidyverse)
library(lubridate)
library(scales)

# Curated, accessible color palette for the 4 flow regimes
REGIME_COLORS <- c(
  "low flow season"       = "#3182CE",  # Solid cool blue (baseflow/winter)
  "shoulder: low to high" = "#38A169",  # Fresh green (spring rising freshet)
  "high flow season"      = "#DD6B20",  # Warm amber/orange (peak runoff)
  "shoulder: high to low" = "#805AD5"   # Royal purple (autumn recession)
)

REGIME_FILLS <- c(
  "low flow season"       = "#EBF8FF",  # Soft blue fill
  "shoulder: low to high" = "#F0FFF4",  # Soft green fill
  "high flow season"      = "#FFFAF0",  # Soft orange fill
  "shoulder: high to low" = "#FAF5FF"   # Soft purple fill
)

#' Plot Day-of-Year Hydrograph Envelopes with Demarcation Line and Regime Zones
#'
#' Generates the core diagnostic envelope plot showing historical minimum, 20th percentile,
#' median, and maximum hydrograph envelopes along with the 85th percentile minimum demarcation line
#' and shaded seasonal regime zones.
#'
#' @param doy_classified Data frame from diagnose_station_regimes()$doy_envelopes
#' @param boundaries Data frame from diagnose_station_regimes()$boundaries
#' @param thresholds Data frame from diagnose_station_regimes()$thresholds
#' @param site_filter Optional character string for specific site_id
#' @param log_scale Logical whether to display discharge on log10 scale (default: TRUE)
#' @return ggplot object
#' @export
plot_regime_envelope <- function(doy_classified, 
                                 boundaries, 
                                 thresholds, 
                                 site_filter = NULL,
                                 log_scale = TRUE) {
  
  if (!is.null(site_filter)) {
    doy_sub <- doy_classified %>% filter(site_id == site_filter)
    bound_sub <- boundaries %>% filter(site_id == site_filter)
    thresh_sub <- thresholds %>% filter(site_id == site_filter)
  } else {
    doy_sub <- doy_classified
    bound_sub <- boundaries
    thresh_sub <- thresholds
  }
  
  # Month breaks for X-axis (DOY of first day of each month)
  month_breaks <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
  month_labels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  
  # Create background regime bands
  regime_bands <- doy_sub %>%
    group_by(agency, site_id, station_name, flow_regime) %>%
    summarise(
      start_doy = min(doy),
      end_doy = max(doy),
      .groups = "drop"
    )
  
  p <- ggplot(doy_sub, aes(x = doy)) +
    # Background shaded regime vertical rectangles
    geom_rect(
      data = doy_sub,
      aes(xmin = doy - 0.5, xmax = doy + 0.5, ymin = -Inf, ymax = Inf, fill = flow_regime),
      alpha = 0.22,
      inherit.aes = FALSE
    ) +
    # 10th to 90th percentile envelope ribbon
    geom_ribbon(
      aes(ymin = pmax(0.01, q_p10), ymax = q_p90),
      fill = "#CBD5E0",
      alpha = 0.45
    ) +
    # Min to Max historical range ribbon
    geom_ribbon(
      aes(ymin = pmax(0.01, q_min), ymax = q_max),
      fill = "#E2E8F0",
      alpha = 0.30
    ) +
    # Smoothed Envelope curves
    geom_line(aes(y = q_max_smooth, color = "Historical Max (Q_max)"), linewidth = 1.0, linetype = "solid") +
    geom_line(aes(y = q_p50_smooth, color = "Historical Median (Q_50)"), linewidth = 1.1, linetype = "solid") +
    geom_line(aes(y = q_p20_smooth, color = "High Flow Baseline (Q_20)"), linewidth = 0.9, linetype = "dotdash") +
    geom_line(aes(y = q_min_smooth, color = "Historical Min (Q_min)"), linewidth = 0.8, linetype = "dotted") +
    # Demarcation Line (85th percentile of daily minimums)
    geom_hline(
      data = thresh_sub,
      aes(yintercept = t_demarcation_cfs),
      color = "#E53E3E",
      linewidth = 1.0,
      linetype = "dashed"
    ) +
    # Vertical boundary demarcation lines
    geom_vline(
      data = bound_sub,
      aes(xintercept = runoff_start_doy),
      color = "#2D3748",
      linetype = "dotted",
      linewidth = 0.7
    ) +
    geom_vline(
      data = bound_sub,
      aes(xintercept = high_flow_start_doy),
      color = "#2D3748",
      linetype = "dotted",
      linewidth = 0.7
    ) +
    geom_vline(
      data = bound_sub,
      aes(xintercept = high_flow_end_doy),
      color = "#2D3748",
      linetype = "dotted",
      linewidth = 0.7
    ) +
    geom_vline(
      data = bound_sub,
      aes(xintercept = runoff_end_doy),
      color = "#2D3748",
      linetype = "dotted",
      linewidth = 0.7
    ) +
    facet_wrap(~ paste0(agency, ": ", site_id, " - ", station_name), scales = "free_y", ncol = 1) +
    scale_x_continuous(
      breaks = month_breaks,
      labels = month_labels,
      limits = c(1, 366),
      expand = c(0, 0)
    ) +
    scale_fill_manual(
      name = "Flow Regime Zone",
      values = REGIME_COLORS
    ) +
    scale_color_manual(
      name = "Hydrograph Quantile Curves",
      values = c(
        "Historical Max (Q_max)" = "#742A2A",
        "Historical Median (Q_50)" = "#1A202C",
        "High Flow Baseline (Q_20)" = "#2C5282",
        "Historical Min (Q_min)" = "#4A5568"
      )
    ) +
    labs(
      title = "Streamflow Regime Classification: Day-of-Year Envelope & Demarcation Threshold",
      subtitle = "Red Dashed Line: Demarcation Threshold (85th %ile of DOY Mins) | Colored Shading: Classified Seasonal Regimes",
      x = "Day of Year (Calendar Month)",
      y = if (log_scale) "Discharge (cfs) [Log10 Scale]" else "Discharge (cfs)",
      caption = "Analysis methodology: Excel Quantile Demarcation Workflow Formalized in R | Watershed & Municipal Hydrology"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13, color = "#1A202C"),
      plot.subtitle = element_text(size = 10, color = "#4A5568", margin = margin(b = 10)),
      strip.background = element_rect(fill = "#EDF2F7", color = NA),
      strip.text = element_text(face = "bold", size = 10, color = "#2D3748"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E2E8F0", linewidth = 0.4),
      panel.spacing = unit(1.2, "lines")
    )
  
  if (log_scale) {
    p <- p + scale_y_log10(labels = label_comma(), expand = expansion(mult = c(0.05, 0.15)))
  } else {
    p <- p + scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0.05, 0.15)))
  }
  
  return(p)
}

#' Plot Multi-Year Daily Time Series Colored by Flow Regime
#'
#' Plots the continuous daily hydrograph across historical years with points/lines colored
#' by their classified seasonal regime.
#'
#' @param daily_classified Data frame from diagnose_station_regimes()$daily_classified
#' @param site_filter Optional site_id
#' @param year_range Optional numeric vector c(start_year, end_year)
#' @return ggplot object
#' @export
plot_regime_timeseries <- function(daily_classified, 
                                   site_filter = NULL, 
                                   year_range = NULL) {
  
  df_plot <- daily_classified
  if (!is.null(site_filter)) {
    df_plot <- df_plot %>% filter(site_id == site_filter)
  }
  if (!is.null(year_range)) {
    df_plot <- df_plot %>% filter(year(date) >= year_range[1] & year(date) <= year_range[2])
  }
  
  p <- ggplot(df_plot, aes(x = date, y = flow_cfs)) +
    geom_line(color = "#A0AEC0", linewidth = 0.4, alpha = 0.7) +
    geom_point(aes(color = flow_regime), size = 0.8, alpha = 0.85) +
    geom_hline(
      aes(yintercept = t_demarcation),
      color = "#E53E3E",
      linetype = "dashed",
      linewidth = 0.7,
      alpha = 0.8
    ) +
    facet_wrap(~ paste0(agency, ": ", site_id, " - ", station_name), scales = "free_y", ncol = 1) +
    scale_color_manual(
      name = "Classified Flow Regime",
      values = REGIME_COLORS
    ) +
    scale_y_log10(labels = label_comma()) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(
      title = "Multi-Year Daily Streamflow Time Series by Classified Hydrologic Regime",
      subtitle = "Red Dashed Line: Baseflow Demarcation Threshold (85th %ile of DOY Mins)",
      x = "Date",
      y = "Daily Mean Discharge (cfs) [Log10 Scale]",
      caption = "Municipal Water Supply & Watershed Quality Regime Classification"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13, color = "#1A202C"),
      plot.subtitle = element_text(size = 10, color = "#4A5568", margin = margin(b = 8)),
      strip.background = element_rect(fill = "#EDF2F7", color = NA),
      strip.text = element_text(face = "bold", size = 10, color = "#2D3748"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 9),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1.0, "lines")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  return(p)
}

#' Plot Annual Regime Calendar Heatmap
#'
#' Displays Day of Year on X-axis and Year on Y-axis with tiles colored by daily regime classification.
#'
#' @param daily_classified Data frame from diagnose_station_regimes()$daily_classified
#' @param site_filter Specific site_id to plot
#' @return ggplot object
#' @export
plot_regime_annual_heatmap <- function(daily_classified, site_filter) {
  
  df_sub <- daily_classified %>%
    filter(site_id == site_filter) %>%
    mutate(
      year_fac = factor(year(date), levels = sort(unique(year(date)), decreasing = TRUE)),
      doy = yday(date)
    )
  
  st_name <- unique(df_sub$station_name)[1]
  agency_nm <- unique(df_sub$agency)[1]
  
  month_breaks <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
  month_labels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  
  p <- ggplot(df_sub, aes(x = doy, y = year_fac, fill = flow_regime)) +
    geom_tile(width = 1, height = 0.9) +
    scale_fill_manual(
      name = "Hydrologic Regime",
      values = REGIME_COLORS
    ) +
    scale_x_continuous(
      breaks = month_breaks,
      labels = month_labels,
      limits = c(1, 366),
      expand = c(0, 0)
    ) +
    labs(
      title = paste0("Annual Flow Regime Progression Matrix: ", agency_nm, " ", site_filter),
      subtitle = paste0(st_name, " | Day-of-Year Seasonal Regime Distribution Across Historical Years"),
      x = "Day of Year (Month)",
      y = "Year",
      caption = "Municipal Water Resources Hydrologic Regime Analysis"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12, color = "#1A202C"),
      plot.subtitle = element_text(size = 9.5, color = "#4A5568", margin = margin(b = 8)),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 9),
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 8.5)
    )
  
  return(p)
}
