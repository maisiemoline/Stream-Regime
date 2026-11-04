# ==============================================================================
# Script: test_classification.R
# Description: Automated test script verifying USGS & DWR regime classification
# ==============================================================================

source("R/fetch_streamflow.R")
source("R/classify_regimes.R")
source("R/plot_regimes.R")

cat("=== 1. Ingesting Streamflow Data ===\n")
test_usgs <- c("06696980", "06700000")
test_dwr  <- c("PLAHARCO", "PLAANTCO")

flow_data <- fetch_all_streamflow(
  usgs_sites = test_usgs,
  dwr_sites = test_dwr,
  start_date = "2000-01-01",
  end_date = "2023-12-31"
)

cat("Total records fetched:", nrow(flow_data), "\n")
print(table(flow_data$site_id, flow_data$agency))

cat("\n=== 2. Running Regime Classification Pipeline ===\n")
results <- diagnose_station_regimes(flow_data, smooth_window = 7, p_baseflow = 0.85, p_high = 0.20)

cat("\n--- Thresholds Summary Table ---\n")
print(results$thresholds)

cat("\n--- Seasonal Regime Boundaries Table ---\n")
print(results$boundaries)

cat("\n--- Daily Classification Breakdown by Site ---\n")
print(
  results$daily_classified %>%
    group_by(agency, site_id, flow_regime) %>%
    summarise(
      days = n(),
      pct_days = round(100 * n() / nrow(.), 2),
      mean_flow = round(mean(flow_cfs, na.rm = TRUE), 1),
      median_flow = round(median(flow_cfs, na.rm = TRUE), 1),
      .groups = "drop"
    )
)

cat("\n=== 3. Testing Plot Generation ===\n")
p_env <- plot_regime_envelope(results$doy_envelopes, results$boundaries, results$thresholds)
ggsave("output_envelope_test.png", plot = p_env, width = 11, height = 12, dpi = 150)
cat("Saved envelope plot to output_envelope_test.png\n")

p_ts <- plot_regime_timeseries(results$daily_classified, year_range = c(2018, 2023))
ggsave("output_timeseries_test.png", plot = p_ts, width = 11, height = 10, dpi = 150)
cat("Saved timeseries plot to output_timeseries_test.png\n")

cat("\n=== All automated tests completed successfully! ===\n")
