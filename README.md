# Streamflow Regime Classification Engine

An R-based hydrologic analysis toolkit developed for municipal water utilities and watershed scientists to diagnose, classify, and visualize streamflow regimes across **USGS** and **Colorado Division of Water Resources (DWR)** streamgage networks.

🌐 **Live Interactive Web Report**: [https://maisiemoline.github.io/Stream-Regime/](https://maisiemoline.github.io/Stream-Regime/)  
📊 **Upper South Platte Assessment**: [https://maisiemoline.github.io/Stream-Regime/stream_regime_classification.html](https://maisiemoline.github.io/Stream-Regime/stream_regime_classification.html)

---

## 🌊 Core Flow Regimes

The classification engine partitions streamflow into four hydrologic seasons:

1. **`low flow season` (Winter Baseflow & Storage Drawdown):**
   Groundwater-dominated winter period; crucial for water treatment plant dilution monitoring, 7Q10 low-flow permit compliance, and winter reservoir storage management.
2. **`shoulder: low to high` (Spring Freshet / Rising Limb):**
   Onset of low-elevation snowmelt and spring runoff; triggers increased coagulant dosing (alum / ferric chloride) for water treatment plants and water rights storage diversion initiation.
3. **`high flow season` (Peak Runoff / Snowmelt Window):**
   High volume snowmelt runoff; peak municipal reservoir storage capture period, requiring management of elevated turbidity, suspended sediment, and TOC.
4. **`shoulder: high to low` (Late Summer / Autumn Recession Limb):**
   Snowpack depletion and hydrograph recession; transition into autumn baseflow, late-summer temperature monitoring, and algal bloom management.

---

## 📁 Repository Structure

```
├── index.html                 # Main live webpage (GitHub Pages) with full multi-basin analysis
├── app.R                      # Interactive Shiny dashboard application
├── stream_regime_classifier.R # Multi-station classification & visualization engine
├── streamflow_daily_cache.rds # Pre-cached daily streamflow dataset for instant app launch
├── colorado_streamflow_regime_seasonal_timing.csv # Multi-station seasonal timing results
├── stream_flow_regime_analysis.Rmd # Multi-basin seasonal timing analysis report
├── stream_flow_regime_analysis.html # Rendered multi-basin analysis report
├── StreamRegime.Rproj         # RStudio Project configuration
├── R/
│   ├── fetch_streamflow.R     # Ingestion for USGS NWIS & CO DWR OpenCDSS REST API v2
│   ├── classify_regimes.R     # Quantile envelope & demarcation classification engine
│   └── plot_regimes.R         # Publication-grade ggplot2 visualizations
├── stream_regime_classification.Rmd   # Comprehensive R Markdown report
├── stream_regime_classification.html  # Rendered standalone HTML report
├── render_report.R            # Script to re-render the HTML report
├── test_classification.R      # Automated test script for pilot sites
├── daily_flow_regimes_classified.csv  # Full daily classified records (2000–2023)
├── station_regime_seasonal_boundaries.csv # Seasonal transition calendar dates
└── station_demarcation_thresholds.csv     # Demarcation line thresholds
```

---

## 🚀 Quickstart Guide

### 1. View Directly in Any Browser
Anyone can view the full interactive report immediately without installing R or any software:  
👉 **[https://maisiemoline.github.io/Stream-Regime/](https://maisiemoline.github.io/Stream-Regime/)**

### 2. Launch the Shiny App from GitHub
If someone has R/RStudio installed, they can launch the interactive Shiny app directly from your GitHub repo with a single command:
```r
shiny::runGitHub("Stream-Regime", "maisiemoline")
```
Or open `StreamRegime.Rproj` in RStudio and click **Run App** (or run `shiny::runApp()`).

### 3. Ingest Data and Classify Any Streamgage


```r
source("R/fetch_streamflow.R")
source("R/classify_regimes.R")
source("R/plot_regimes.R")

# Fetch streamflow for any USGS and DWR stations
flow_data <- fetch_all_streamflow(
  usgs_sites = c("06696980", "06700000"),
  dwr_sites  = c("PLAHARCO", "PLAANTCO"),
  start_date = "2000-01-01",
  end_date   = "2023-12-31"
)

# Run classification engine
results <- diagnose_station_regimes(
  df = flow_data,
  smooth_window = 7,     # 7-day circular moving average
  p_baseflow = 0.85,     # 85th percentile of daily minimums for demarcation
  p_high = 0.20          # 20th percentile for high-flow curve
)

# View seasonal transition dates
print(results$boundaries)

# Generate diagnostic envelope plot
p <- plot_regime_envelope(results$doy_envelopes, results$boundaries, results$thresholds)
print(p)
```

### 4. Re-render the HTML Report
```bash
Rscript render_report.R
```

---

## 📊 Summary of Pilot Stations (Upper South Platte Basin, CO)

| Station ID | Source | Station Name | Demarcation Line | Runoff Onset | High Flow Window | Recession End |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **06696980** | USGS | Tarryall Creek at Upper Station nr Como | **8.76 cfs** | Apr 18 | May 06 – Aug 27 | Oct 17 |
| **06700000** | USGS | South Platte River above Cheesman Lake | **100.88 cfs** | May 10 | May 22 – Sep 06 | Oct 18 |
| **PLAHARCO** | DWR | South Platte River above Elevenmile Res. | **58.75 cfs** | Apr 13 | Apr 28 – Oct 06 | Oct 21 |
| **PLAANTCO** | DWR | South Platte River below Antero Reservoir | **3.27 cfs** | May 17 | Jun 01 – Jul 25 | Aug 09 |
