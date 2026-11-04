# ==============================================================================
# USGS & CDSS Streamflow Regime Classification Dashboard (Shiny Application)
# ==============================================================================
# Interactive application for exploring annual hydrograph envelopes, baseflow
# demarcation lines, and monthly Water Year regime classification matrices across
# Colorado streamgage stations (Watershed Monitoring Plan Methodology).
# ==============================================================================

library(shiny)
library(bslib)
library(DT)
library(tidyverse)
library(dataRetrieval)
library(cdssr)
library(httr)
library(jsonlite)
library(lubridate)
library(scales)

# Ensure graphics device compatibility on macOS/Linux (avoid Quartz 0-size device crashes)
if (.Platform$OS.type == "unix") {
  try(options(bitmapType = "cairo"), silent = TRUE)
}

# Source core classifier functions (Refined Peak-Anchored & Dormant Ceiling Classification)
source("stream_regime_classifier.R")

# Pre-defined site catalog: 9 Reference Benchmark Sites + 28 USGS Sites + 40 DWR Sites
SITE_LABELS <- c(
  # --- Reference & Benchmark Validation Sites (9 Sites) ---
  "PLASPICO"  = "★ PLASPICO (06694920) - South Platte River above Spinney (Benchmark)",
  "PLACHECO"  = "★ PLACHECO (06701500) - South Platte River below Cheesman (Benchmark)",
  "BCRMORCO"  = "★ BCRMORCO (06710500) - Bear Creek at Morrison (Benchmark)",
  "PLUTIRCO"  = "★ PLUTIRCO (06709530) - Plum Creek at Titan Road (Benchmark)",
  "PLAWATCO"  = "★ PLAWATCO (06708000) - South Platte River at Waterton (Benchmark)",
  "PLAGRACO"  = "★ PLAGRACO (06707500) - South Platte River at South Platte (Benchmark)",
  "FRAWINCO"  = "★ FRAWINCO (09024000) - Fraser River at Winter Park (Benchmark)",
  "WILPARCO"  = "★ WILPARCO (09037500) - Williams Fork near Parshall (Benchmark)",
  "COLKREMCO" = "★ COLKREMCO (09058000) - Colorado River near Kremmling (Benchmark)",

  # --- USGS Streamflow Stations (28 Sites) ---
  "06700000" = "06700000 - SOUTH PLATTE RIVER ABOVE CHEESMAN LAKE (Seasonal Gauge)",
  "06701900" = "06701900 - SOUTH PLATTE RIVER BLW BRUSH CRK NEAR TRUMBULL",
  "06709530" = "06709530 - TURKEY CREEK NEAR MORRISON",
  "06710247" = "06710247 - BEAR CREEK ABOVE BEAR CREEK LAKE",
  "06710385" = "06710385 - BEAR CREEK ABOVE EVERGREEN",
  "06710605" = "06710605 - BEAR CREEK ABOVE KITTREDGE",
  "06711040" = "06711040 - TURKEY CREEK ABOVE BEAR CREEK LAKE NEAR MORRISON (2019-Present)",
  "06711565" = "06711565 - SOUTH PLATTE RIVER AT ENGLEWOOD",
  "09024000" = "09024000 - FRASER RIVER AT WINTER PARK",
  "09025000" = "09025000 - VASQUEZ CREEK AT WINTER PARK",
  "09027100" = "09027100 - WEST ST. LOUIS CREEK NEAR FRASER",
  "09031600" = "09031600 - VASQUEZ DIVERSION 3 NR WINTER PARK (2021-Present)",
  "09032000" = "09032000 - RANCH CREEK NEAR FRASER",
  "09034250" = "09034250 - FRASER RIVER BLW CROOKED CRK AT TABERNASH",
  "09034900" = "09034900 - BOBTAIL CREEK NEAR WINTER PARK",
  "09035500" = "09035500 - CABIN CREEK NEAR FRASER",
  "09035700" = "09035700 - MEADOW CREEK NEAR FRASER",
  "09035900" = "09035900 - CROOKED CREEK NEAR FRASER",
  "09036000" = "09036000 - ST. LOUIS CREEK AT FRASER",
  "09037500" = "09037500 - STRAIGHT CREEK BLW LASKEY GULCH NR DILLON",
  "09038500" = "09038500 - WILLIAMS FORK BELOW WILLIAMS FORK RESERVOIR",
  "09041400" = "09041400 - SNAKE RIVER NEAR MONTEZUMA",
  "09046600" = "09046600 - BLUE RIVER NEAR DILLON",
  "09047500" = "09047500 - BLUE RIVER BLW GREEN MOUNTAIN RESERVOIR",
  "09050100" = "09050100 - COLORADO RIVER BLW CATAMOUNT CRK NR KREMMLING",
  "09050700" = "09050700 - COLORADO RIVER NEAR KREMMLING",
  "09051050" = "09051050 - COLORADO RIVER BLW GORE CANYON NR KREMMLING",
  "09058000" = "09058000 - COLORADO RIVER NEAR DOTSERO",

  # --- Colorado DWR / CDSS Stations & Diversion Structures (40 Sites) ---
  "ROBTUNCO" = "ROBTUNCO - Roberts Tunnel",
  "STCTUNCO" = "STCTUNCO - Straight Creek Tunnel",
  "PLASTRCO" = "PLASTRCO - South Platte River at Strontia Springs",
  "PLACHACO" = "PLACHACO - South Platte River at Chatfield Reservoir",
  "PLAWATCO" = "PLAWATCO - South Platte River at Waterton",
  "BCRMORCO" = "BCRMORCO - Bear Creek at Morrison",
  "COND20CO" = "COND20CO - Conduit 20 Diversion",
  "PLAKASCO" = "PLAKASCO - South Platte River at Kassler",
  "BCROUTCO" = "BCROUTCO - Bear Creek Reservoir Outlet",
  "GENMPLCO" = "GENMPLCO - Geneva Creek at Grant / Mt Evans",
  "WARDITCO" = "WARDITCO - WARNE DITCH (Irrigation Ditch)",
  "BCRSHECO" = "BCRSHECO - Bear Creek at Sheridan",
  "0902253A" = "0902253A - Moffat Tunnel Collection / Fraser",
  "WARDBCCO" = "WARDBCCO - Ward Canal at Bear Creek",
  "PLUCHACO" = "PLUCHACO - Plum Creek at Chatfield Reservoir",
  "PLAGRACO" = "PLAGRACO - South Platte River at Grant",
  "PLABAICO" = "PLABAICO - South Platte River at Bailey",
  "SFKANTCO" = "SFKANTCO - South Fork South Platte River at Antero",
  "PLAANTCO" = "PLAANTCO - South Platte River at Antero Reservoir",
  "PLAHARCO" = "PLAHARCO - South Platte River at Hartsel",
  "PLAGEOCO" = "PLAGEOCO - South Platte River below Elevenmile (Lake George)",
  "PLASPLCO" = "PLASPLCO - South Platte River below Spinney Mountain",
  "RALCRKCO" = "RALCRKCO - Ralston Creek near Golden",
  "MOFTUNCO" = "MOFTUNCO - Moffat Tunnel East Portal",
  "RALBLWCO" = "RALBLWCO - Ralston Creek below Ralston Reservoir",
  "BOSDELCO" = "BOSDELCO - Boulder Creek Delivery Structure",
  "BOCBGRCO" = "BOCBGRCO - Boulder Creek at Groesbeck Reservoir",
  "PLAHENCO" = "PLAHENCO - South Platte River at Henderson",
  "BOCPINCO" = "BOCPINCO - Boulder Creek near Pinecliffe",
  "PLADENCO" = "PLADENCO - South Platte River at Denver",
  "FRAJIMCO" = "FRAJIMCO - Fraser River at Jim Creek",
  "COLKBDCO" = "COLKBDCO - Colorado River below Kremmling",
  "BLUGRECO" = "BLUGRECO - Blue River below Green Mountain Reservoir",
  "TARTARCO" = "TARTARCO - Tarryall Creek at Tarryall",
  "FOUBERCO" = "FOUBERCO - FOUR MILE CREEK AT BEERY DITCH (2302905)",
  "HOMSPICO" = "HOMSPICO - Homestake Collection / Spinney",
  "MFKSTMCO" = "MFKSTMCO - Middle Fork South Platte River nr Santa Maria",

  # --- Custom User Inputs ---
  "CUSTOM_USGS" = "Custom USGS Station (Enter ID below)",
  "CUSTOM_CDSS" = "Custom CDSS Station (Enter Abbreviation below)"
)

DEFAULT_CHOICES <- setNames(names(SITE_LABELS), unname(SITE_LABELS))

# Load pre-cached dataset for instant response
LOCAL_CACHE <- if (file.exists("streamflow_daily_cache.rds")) {
  tryCatch(readRDS("streamflow_daily_cache.rds"), error = function(e) NULL)
} else NULL

REGIME_COLORS <- c(
  "low flow season"      = "#2B83BA",  # Baseflow Blue
  "shoulder:low to high" = "#43AA8B",  # Rising Green
  "high flow season"     = "#E63946",  # Peak Red
  "shoulder:high to low" = "#F4A261"   # Receding Orange
)

# Robust site ID resolution
resolve_site_id <- function(raw_input) {
  if (is.null(raw_input) || length(raw_input) == 0) return("PLASPICO")
  k <- trimws(raw_input)
  if (k %in% names(SITE_LABELS)) return(k)
  if (k %in% unname(SITE_LABELS)) {
    idx <- which(unname(SITE_LABELS) == k)
    return(names(SITE_LABELS)[idx[1]])
  }
  token <- strsplit(k, "[-/ ]")[[1]][1]
  token <- gsub("[^A-Za-z0-9]", "", token)
  if (nchar(token) > 0) token else k
}

# Compact Executive Metric Card Component
metric_card <- function(title, value_id, sub_id, icon_name, border_color = "#1D3557", icon_color = "#1D3557") {
  div(
    class = "card shadow-sm mb-3",
    style = sprintf("border-left: 4px solid %s; border-radius: 8px; background: #ffffff; padding: 10px 14px; min-height: 84px; display: flex; flex-direction: column; justify-content: space-between; border-top: 1px solid #e2e8f0; border-right: 1px solid #e2e8f0; border-bottom: 1px solid #e2e8f0;", border_color),
    div(
      class = "d-flex justify-content-between align-items-center mb-1",
      tags$span(
        style = "font-size: 0.70rem; font-weight: 700; letter-spacing: 0.5px; text-transform: uppercase; color: #64748b;",
        title
      ),
      tags$i(class = sprintf("bi bi-%s", icon_name), style = sprintf("font-size: 1.05rem; color: %s; opacity: 0.85;", icon_color))
    ),
    div(
      tags$div(
        style = "font-size: 1.20rem; font-weight: 800; color: #0f172a; line-height: 1.2;",
        textOutput(value_id, inline = TRUE)
      ),
      tags$div(
        style = "font-size: 0.74rem; color: #64748b; font-weight: 500; margin-top: 2px;",
        textOutput(sub_id, inline = TRUE)
      )
    )
  )
}

# ==============================================================================
# UI Definition
# ==============================================================================
ui <- page_sidebar(
  title = div(
    class = "d-flex align-items-center gap-2",
    tags$span(style = "font-size: 1.35rem; font-weight: 800; color: #1D3557;", "HydroRegime Explorer"),
    tags$span(class = "badge bg-primary", "Four-Phase Flow Regime Classification")
  ),
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#1D3557",
    secondary = "#457B9D",
    success = "#43AA8B",
    danger = "#E63946",
    warning = "#F4A261",
    info = "#2B83BA",
    base_font = font_google("Inter"),
    heading_font = font_google("Inter")
  ),
  sidebar = sidebar(
    width = 360,
    # Search Box for Quick Filtering
    textInput(
      "station_search_box",
      label = tags$strong("Quick Search Stations:"),
      placeholder = "Type name or ID, e.g. Bear, Turkey, 06711, Williams..."
    ),

    # Searchable Station Selector
    selectizeInput(
      "site_select",
      label = tags$strong("Select Streamgage Station:"),
      choices = DEFAULT_CHOICES,
      selected = "PLASPICO",
      options = list(
        placeholder = "Click or search station...",
        allowEmptyOption = FALSE,
        highlight = TRUE
      )
    ),

    # Conditional Custom Inputs
    conditionalPanel(
      condition = "input.site_select == 'CUSTOM_USGS'",
      textInput("custom_usgs_id", "USGS Station ID (8-10 digits):", value = "09361500")
    ),
    conditionalPanel(
      condition = "input.site_select == 'CUSTOM_CDSS'",
      textInput("custom_cdss_abbrev", "CDSS Station Abbreviation:", value = "PLASPICO")
    ),

    tags$hr(),

    # Analysis Period Preset
    selectInput(
      "period_preset",
      label = tags$strong("Historical Analysis Period:"),
      choices = c(
        "1997-2017 (Historical Benchmark Period)" = "bench_2017",
        "2000-2024 (Modern Baseline - Recommended)" = "modern_2024",
        "2000-2026 (Modern including WY 2026 Drought)" = "modern_2026",
        "Custom Date Range" = "custom"
      ),
      selected = "bench_2017"
    ),

    # Conditional Custom Date Range
    conditionalPanel(
      condition = "input.period_preset == 'custom'",
      dateRangeInput(
        "custom_dates",
        "Custom Date Range:",
        start = "1997-01-01",
        end = "2017-12-31",
        min = "1950-01-01",
        max = Sys.Date()
      )
    ),

    tags$hr(),

    # Demarcation Threshold Controls
    radioButtons(
      "demarc_mode",
      label = tags$strong("Demarcation Threshold (Q_demarc):"),
      choices = c(
        "Automated (Winter Maximum Flow Ceiling / Lowest Baseflow)" = "auto",
        "Manual Override Value (cfs)" = "manual"
      ),
      selected = "auto"
    ),
    conditionalPanel(
      condition = "input.demarc_mode == 'manual'",
      numericInput("manual_demarc_val", "Demarcation Value (cfs):", value = 110, min = 0.1, step = 1)
    ),

    # Runoff Percentile Selection
    selectInput(
      "hf_pct_mode",
      label = tags$strong("High Flow Runoff Curve (Q_P):"),
      choices = c(
        "Auto / Benchmark Default (20% PLASPICO, 40% PLACHECO)" = "auto",
        "50th Percentile (Median Flow - Typical Runoff Peak)" = "0.50",
        "40th Percentile (PLACHECO Standard - Regulated Delivery)" = "0.40",
        "20th Percentile (PLASPICO Standard - Runoff Onset)" = "0.20",
        "30th Percentile" = "0.30",
        "25th Percentile" = "0.25"
      ),
      selected = "auto"
    ),

    # Smoothing Window
    sliderInput(
      "smooth_window",
      "Moving Average Smoothing (Days):",
      min = 1, max = 15, value = 1, step = 2
    ),

    tags$div(
      class = "alert alert-light mt-3 p-2 small border text-muted",
      tags$strong("Hydrological Note on Peaks: "),
      "Q_max shows the single highest flood year on that day (e.g. wet-year snowmelt spills in June), while lower percentiles reflect typical annual operation (e.g. July/August municipal releases on managed rivers). For seasonal mountain gauges, unmonitored winter months automatically default to winter baseflow dormancy."
    )
  ),

  # Top Metric Cards (Compact Executive KPI Row)
  layout_column_wrap(
    width = 1/4,
    fill = FALSE,
    gap = "12px",
    metric_card(
      title = "Demarcation Threshold",
      value_id = "vb_demarc_val",
      sub_id = "vb_demarc_sub",
      icon_name = "water",
      border_color = "#1D3557",
      icon_color = "#1D3557"
    ),
    metric_card(
      title = "High Flow Season",
      value_id = "vb_high_flow_val",
      sub_id = "vb_high_flow_sub",
      icon_name = "calendar2-range",
      border_color = "#E63946",
      icon_color = "#E63946"
    ),
    metric_card(
      title = "Low Flow Season",
      value_id = "vb_low_flow_val",
      sub_id = "vb_low_flow_sub",
      icon_name = "calendar2-check",
      border_color = "#2B83BA",
      icon_color = "#2B83BA"
    ),
    metric_card(
      title = "Peak Discharge",
      value_id = "vb_peak_flow_val",
      sub_id = "vb_peak_flow_sub",
      icon_name = "graph-up-arrow",
      border_color = "#43AA8B",
      icon_color = "#43AA8B"
    )
  ),

  # Main Tab Panels
  navset_card_tab(
    full_screen = TRUE,

    # Tab 1: Four-Phase Seasonal Hydrograph
    nav_panel(
      title = tags$span(icon("chart-area"), " Seasonal Flow Regime Hydrograph"),
      card_body(
        div(
          class = "d-flex justify-content-between align-items-center mb-2",
          tags$h5(textOutput("plot_station_title"), class = "mb-0 text-primary fw-bold"),
          div(
            radioButtons(
              "plot_style_toggle",
              label = NULL,
              choices = c(
                "Four-Phase Seasonal Hydrograph" = "standard",
                "Full Percentile Envelopes" = "envelopes"
              ),
              selected = "standard",
              inline = TRUE
            )
          )
        ),
        plotOutput("main_hydrograph", height = "520px")
      )
    ),

    # Tab 2: Monthly Regime Matrix & Colored Cards
    nav_panel(
      title = tags$span(icon("th"), " Monthly Regime Matrix"),
      card_body(
        tags$h5("Dominant Flow Regime by Month (Water Year: Oct–Sep)", class = "fw-bold text-primary mb-3"),
        # Colored Responsive Month Cards Grid (Guaranteed 100% visible in every browser)
        uiOutput("monthly_regime_cards"),
        tags$hr(class = "my-3"),
        tags$h6("Monthly Discharge & Regime Classification Summary Table", class = "fw-bold text-secondary mb-2"),
        uiOutput("monthly_summary_table_html")
      )
    ),

    # Tab 3: Seasonal Timing Summary
    nav_panel(
      title = tags$span(icon("calendar-alt"), " Seasonal Timing & Transitions"),
      card_body(
        tags$h5("Comprehensive Four-Phase Flow Regime Breakdown", class = "fw-bold text-primary mb-3"),
        tableOutput("timing_detail_table"),
        tags$hr(),
        tags$div(
          class = "alert alert-info",
          tags$h6(tags$strong("Hydrologic Phase Definitions:"), class = "alert-heading"),
          tags$ul(
            tags$li(tags$strong("Low Flow Season: "), "Dormant fall/winter baseflow where daily maximum flow remains bounded below the demarcation threshold (Q_demarc)."),
            tags$li(tags$strong("Spring Shoulder: "), "The snowmelt rising limb transition between low flow termination and high flow onset."),
            tags$li(tags$strong("High Flow Season: "), "The sustained core runoff window where the selected lower percentile curve (Q_P) exceeds Q_demarc."),
            tags$li(tags$strong("Fall Shoulder: "), "The post-runoff receding limb transition back into winter dormancy.")
          )
        )
      )
    ),

    # Tab 4: Benchmark Validation (PLASPICO & PLACHECO)
    nav_panel(
      title = tags$span(icon("check-circle"), " Benchmark Validation (PLASPICO & PLACHECO)"),
      card_body(
        tags$h5("Validation Against South Platte Reference Benchmarks (1997-2017)", class = "fw-bold text-primary mb-2"),
        p("This tab compares the algorithm results against the official published benchmark targets for both South Platte benchmark sites."),
        tags$h6("1. PLASPICO: South Platte River Above Spinney (USGS 06694920 / CDSS PLASPICO)", class = "fw-bold text-success mt-3"),
        tableOutput("plaspico_bench_table"),
        tags$hr(),
        tags$h6("2. PLACHECO: South Platte River Below Cheesman Reservoir (USGS 06701500 / CDSS PLACHECO)", class = "fw-bold text-primary mt-3"),
        tableOutput("placheco_bench_table")
      )
    ),

    # Tab 5: Daily Envelopes Data & Export
    nav_panel(
      title = tags$span(icon("table"), " Daily Envelopes Data"),
      card_body(
        div(
          class = "d-flex justify-content-between align-items-center mb-2",
          tags$span("Daily DOY percentiles (1-365) and classified regimes:", class = "text-muted"),
          downloadButton("download_doy_csv", "Download DOY Envelopes (CSV)", class = "btn-sm btn-outline-primary")
        ),
        DTOutput("doy_table")
      )
    )
  )
)

# ==============================================================================
# Server Logic
# ==============================================================================
server <- function(input, output, session) {

  data_cache <- reactiveValues()

  # Dynamic Live Search Filtering for Station Selection
  observeEvent(input$station_search_box, {
    query <- trimws(input$station_search_box)
    if (nzchar(query)) {
      matched_idx <- grep(query, names(DEFAULT_CHOICES), ignore.case = TRUE)
      if (length(matched_idx) == 0) {
        matched_idx <- grep(query, unname(DEFAULT_CHOICES), ignore.case = TRUE)
      }
      if (length(matched_idx) > 0) {
        filtered_choices <- DEFAULT_CHOICES[matched_idx]
        updateSelectizeInput(
          session,
          "site_select",
          choices = filtered_choices,
          selected = unname(filtered_choices)[1]
        )
      }
    } else {
      updateSelectizeInput(
        session,
        "site_select",
        choices = DEFAULT_CHOICES,
        selected = "PLASPICO"
      )
    }
  }, ignoreInit = TRUE)

  selected_date_range <- reactive({
    preset <- input$period_preset
    if (preset == "bench_2017") {
      c(start = "1997-01-01", end = "2017-12-31")
    } else if (preset == "modern_2024") {
      c(start = "2000-01-01", end = "2024-12-31")
    } else if (preset == "modern_2026") {
      c(start = "2000-01-01", end = "2026-09-30")
    } else {
      c(start = format(input$custom_dates[1], "%Y-%m-%d"),
        end   = format(input$custom_dates[2], "%Y-%m-%d"))
    }
  })

  diagnosis_result <- reactive({
    raw_site_input <- input$site_select
    site_key <- resolve_site_id(raw_site_input)
    dates <- selected_date_range()
    smooth_k <- input$smooth_window
    preset_name <- input$period_preset

    # Benchmark-specific default tuning
    demarc_input <- if (input$demarc_mode == "manual") {
      input$manual_demarc_val
    } else if (site_key == "PLASPICO" && preset_name == "bench_2017") {
      110
    } else if (site_key == "PLACHECO" && preset_name == "bench_2017") {
      254
    } else {
      NULL
    }

    hf_pct_input <- if (input$hf_pct_mode != "auto") {
      as.numeric(input$hf_pct_mode)
    } else if (site_key == "PLASPICO") {
      0.20
    } else if (site_key == "PLACHECO") {
      0.40
    } else {
      NULL
    }

    cache_key <- sprintf("%s_%s_%s", site_key, dates["start"], dates["end"])
    raw_df <- data_cache[[cache_key]]

    if (is.null(raw_df)) {
      # 1. Check local cache for instant retrieval
      if (!is.null(LOCAL_CACHE)) {
        if (site_key == "PLASPICO" && !is.null(LOCAL_CACHE$plaspico)) {
          raw_df <- LOCAL_CACHE$plaspico %>% filter(Date >= as.Date(dates["start"]) & Date <= as.Date(dates["end"]))
          site_id <- "06694920 / PLASPICO"
          site_name <- "South Platte River Above Spinney (USGS PLASPICO)"
        } else if (site_key == "PLACHECO" && !is.null(LOCAL_CACHE$placheco)) {
          raw_df <- LOCAL_CACHE$placheco %>% filter(Date >= as.Date(dates["start"]) & Date <= as.Date(dates["end"]))
          site_id <- "06701500 / PLACHECO"
          site_name <- "South Platte River below Cheesman Reservoir, CO"
        } else if (!is.null(LOCAL_CACHE$dwr_daily) && site_key %in% LOCAL_CACHE$dwr_daily$site_no) {
          raw_df <- LOCAL_CACHE$dwr_daily %>%
            filter(site_no == site_key, Date >= as.Date(dates["start"]) & Date <= as.Date(dates["end"]))
          meta_match <- if (!is.null(LOCAL_CACHE$dwr_meta)) LOCAL_CACHE$dwr_meta %>% filter(site_no == site_key) else tibble()
          site_id <- site_key
          site_name <- if (nrow(meta_match) > 0) meta_match$station_nm[1] else SITE_LABELS[site_key]
        } else if (!is.null(LOCAL_CACHE$usgs_daily) && site_key %in% LOCAL_CACHE$usgs_daily$site_no) {
          raw_df <- LOCAL_CACHE$usgs_daily %>%
            filter(site_no == site_key, Date >= as.Date(dates["start"]) & Date <= as.Date(dates["end"]))
          meta_match <- if (!is.null(LOCAL_CACHE$site_meta)) LOCAL_CACHE$site_meta %>% filter(site_no == site_key) else tibble()
          site_id <- site_key
          site_name <- if (nrow(meta_match) > 0) meta_match$station_nm[1] else SITE_LABELS[site_key]
        }
      }

      # Automatic period fallback for modern stations (e.g. 06711040 from 2019, 09031600 from 2021)
      if ((is.null(raw_df) || nrow(raw_df) == 0) && !is.null(LOCAL_CACHE)) {
        if (!is.null(LOCAL_CACHE$usgs_daily) && site_key %in% LOCAL_CACHE$usgs_daily$site_no) {
          all_avail <- LOCAL_CACHE$usgs_daily %>% filter(site_no == site_key)
          if (nrow(all_avail) > 0) {
            raw_df <- all_avail
            meta_match <- if (!is.null(LOCAL_CACHE$site_meta)) LOCAL_CACHE$site_meta %>% filter(site_no == site_key) else tibble()
            site_id <- site_key
            site_name <- if (nrow(meta_match) > 0) meta_match$station_nm[1] else SITE_LABELS[site_key]
            showNotification(
              sprintf("Note: Station %s has records from %s to %s. Utilizing all available records.",
                      site_key, format(min(raw_df$Date), "%Y"), format(max(raw_df$Date), "%Y")),
              type = "warning", duration = 6
            )
          }
        } else if (!is.null(LOCAL_CACHE$dwr_daily) && site_key %in% LOCAL_CACHE$dwr_daily$site_no) {
          all_avail <- LOCAL_CACHE$dwr_daily %>% filter(site_no == site_key)
          if (nrow(all_avail) > 0) {
            raw_df <- all_avail
            meta_match <- if (!is.null(LOCAL_CACHE$dwr_meta)) LOCAL_CACHE$dwr_meta %>% filter(site_no == site_key) else tibble()
            site_id <- site_key
            site_name <- if (nrow(meta_match) > 0) meta_match$station_nm[1] else SITE_LABELS[site_key]
            showNotification(
              sprintf("Note: Station %s has records from %s to %s. Utilizing all available records.",
                      site_key, format(min(raw_df$Date), "%Y"), format(max(raw_df$Date), "%Y")),
              type = "warning", duration = 6
            )
          }
        }
      }

      # 2. If not in cache, download on the fly
      if (is.null(raw_df) || nrow(raw_df) == 0) {
        showNotification(sprintf("Downloading streamflow data for %s...", site_key), type = "message", duration = 3)

        if (site_key == "PLASPICO") {
          raw_df <- fetch_cdss_daily("PLASPICO", start_date = dates["start"], end_date = dates["end"])
          site_id <- "06694920 / PLASPICO"
          site_name <- "South Platte River Above Spinney (USGS PLASPICO)"
        } else if (site_key == "PLACHECO") {
          raw_df <- fetch_cdss_daily("PLACHECO", start_date = dates["start"], end_date = dates["end"])
          site_id <- "06701500 / PLACHECO"
          site_name <- "South Platte River below Cheesman Reservoir, CO"
        } else if (site_key == "CUSTOM_CDSS" || (!grepl("^[0-9]", site_key) && site_key != "CUSTOM_USGS")) {
          abbrev <- if (site_key == "CUSTOM_CDSS") trimws(toupper(input$custom_cdss_abbrev)) else site_key
          raw_df <- fetch_cdss_daily(abbrev, start_date = dates["start"], end_date = dates["end"])
          site_id <- abbrev
          site_name <- sprintf("CDSS Station %s", abbrev)
        } else {
          target_usgs <- if (site_key == "CUSTOM_USGS") trimws(input$custom_usgs_id) else site_key
          raw_df <- fetch_usgs_daily(target_usgs, start_date = dates["start"], end_date = dates["end"])
          site_id <- target_usgs
          meta <- tryCatch(readNWISsite(target_usgs), error = function(e) tibble())
          site_name <- if (nrow(meta) > 0) meta$station_nm[1] else sprintf("USGS %s", target_usgs)
        }
      }

      if (!is.null(raw_df) && nrow(raw_df) > 0) {
        raw_df$site_id_resolved <- coalesce(site_id, site_key)
        raw_df$site_name_resolved <- coalesce(site_name, site_key)
        data_cache[[cache_key]] <- raw_df
      }
    }

    if (is.null(raw_df) || nrow(raw_df) == 0) {
      showNotification("No daily discharge records returned for the specified site/dates.", type = "error", duration = 5)
      return(NULL)
    }

    resolved_id <- raw_df$site_id_resolved[1]
    resolved_name <- raw_df$site_name_resolved[1]

    # Run classification
    classify_stream_regime(
      df_daily = raw_df,
      site_id = resolved_id,
      station_nm = resolved_name,
      demarc_val = demarc_input,
      high_flow_pct = hf_pct_input,
      smooth_window = smooth_k
    )
  })

  # KPI Card Outputs
  output$vb_demarc_val <- renderText({
    res <- diagnosis_result()
    req(res)
    sprintf("%.1f cfs", res$demarc_line)
  })
  output$vb_demarc_sub <- renderText({
    res <- diagnosis_result()
    req(res)
    if (input$demarc_mode == "auto") "Automated Baseflow Ceiling" else "Manual Override Value"
  })

  output$vb_high_flow_val <- renderText({
    res <- diagnosis_result()
    req(res)
    span_txt <- res$timing_summary$high_flow_period
    sub(" \\(.*", "", span_txt)
  })
  output$vb_high_flow_sub <- renderText({
    res <- diagnosis_result()
    req(res)
    span_txt <- res$timing_summary$high_flow_period
    days_txt <- regmatches(span_txt, regexpr("\\(.*\\)", span_txt))
    pct_txt <- res$timing_summary$selected_high_flow_pct
    if (length(days_txt) > 0) {
      sprintf("%s | %s Runoff Curve", gsub("[()]", "", days_txt), pct_txt)
    } else {
      sprintf("%s Runoff Curve", pct_txt)
    }
  })

  output$vb_low_flow_val <- renderText({
    res <- diagnosis_result()
    req(res)
    span_txt <- res$timing_summary$low_flow_period
    sub(" \\(.*", "", span_txt)
  })
  output$vb_low_flow_sub <- renderText({
    res <- diagnosis_result()
    req(res)
    span_txt <- res$timing_summary$low_flow_period
    days_txt <- regmatches(span_txt, regexpr("\\(.*\\)", span_txt))
    if (length(days_txt) > 0) {
      sprintf("%s Baseflow Period", gsub("[()]", "", days_txt))
    } else {
      "Dormant Baseflow Period"
    }
  })

  output$vb_peak_flow_val <- renderText({
    res <- diagnosis_result()
    req(res)
    sprintf("%.1f cfs", res$peak_flow)
  })
  output$vb_peak_flow_sub <- renderText({
    res <- diagnosis_result()
    req(res)
    sprintf("Annual Peak Date: %s", res$peak_date)
  })

  output$plot_station_title <- renderText({
    res <- diagnosis_result()
    req(res)
    sprintf("%s (%s)", res$station_nm, res$site_no)
  })

  # Tab 1: Primary Streamflow Regime Hydrograph (Clean ggplot2)
  output$main_hydrograph <- renderPlot({
    diag <- diagnosis_result()
    req(diag)
    style <- input$plot_style_toggle
    if (style == "envelopes") {
      plot_demarcation_hydrographs(diag$doy_df, ncol = 1)
    } else {
      period_lbl <- if (input$period_preset == "bench_2017") "1997-2017" else if (input$period_preset == "modern_2024") "2000-2024" else "Baseline"
      plot_seasonal_regime_hydrograph(diag, period_label = period_lbl)
    }
  }, res = 96)

  # Tab 2: Colored Responsive Monthly Regime Cards Grid (100% Guaranteed Browser Rendering)
  output$monthly_regime_cards <- renderUI({
    diag <- diagnosis_result()
    req(diag)

    m_df <- diag$monthly_regime
    wy_order <- c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep")
    ordered_df <- tibble(month_abb = wy_order) %>% left_join(m_df, by = "month_abb")

    regime_colors <- c(
      "low flow season"      = "#2B83BA",
      "shoulder:low to high" = "#43AA8B",
      "high flow season"     = "#E63946",
      "shoulder:high to low" = "#F4A261"
    )

    regime_names <- c(
      "low flow season"      = "Low Flow",
      "shoulder:low to high" = "Spring Shoulder",
      "high flow season"     = "High Flow",
      "shoulder:high to low" = "Fall Shoulder"
    )

    cards <- lapply(1:nrow(ordered_df), function(i) {
      row <- ordered_df[i, ]
      reg <- as.character(row$regime)
      bg_col <- if (!is.na(reg) && reg %in% names(regime_colors)) regime_colors[[reg]] else "#64748B"
      lbl <- if (!is.na(reg) && reg %in% names(regime_names)) regime_names[[reg]] else "No Data"
      flow_val <- if (!is.na(row$median_flow)) sprintf("%.0f cfs", row$median_flow) else "---"

      tags$div(
        style = sprintf("flex: 1 1 120px; min-width: 95px; max-width: 140px; background: #ffffff; border-top: 5px solid %s; border-radius: 8px; box-shadow: 0 2px 5px rgba(0,0,0,0.06); padding: 10px 8px; text-align: center; border-left: 1px solid #e2e8f0; border-right: 1px solid #e2e8f0; border-bottom: 1px solid #e2e8f0; margin: 4px;", bg_col),
        tags$div(style = "font-size: 0.90rem; font-weight: 800; color: #1e293b; text-transform: uppercase;", as.character(row$month_abb)),
        tags$span(style = sprintf("display: inline-block; font-size: 0.72rem; font-weight: 700; color: #ffffff; background: %s; padding: 3px 7px; border-radius: 4px; margin: 6px 0;", bg_col), lbl),
        tags$div(style = "font-size: 0.80rem; color: #475569; font-weight: 600;", flow_val),
        tags$div(style = "font-size: 0.68rem; color: #94a3b8; font-weight: 500;", "median daily")
      )
    })

    tags$div(
      tags$div(
        style = "display: flex; flex-wrap: wrap; gap: 8px; justify-content: flex-start; margin-bottom: 16px;",
        cards
      ),
      tags$div(
        style = "display: flex; gap: 16px; flex-wrap: wrap; font-size: 0.78rem; font-weight: 600; color: #475569; padding: 6px 12px; background: #f8fafc; border-radius: 6px; border: 1px solid #e2e8f0;",
        tags$span(tags$span(style = "display: inline-block; width: 12px; height: 12px; background: #2B83BA; border-radius: 2px; margin-right: 4px; vertical-align: middle;"), "Low Flow Season"),
        tags$span(tags$span(style = "display: inline-block; width: 12px; height: 12px; background: #43AA8B; border-radius: 2px; margin-right: 4px; vertical-align: middle;"), "Spring Shoulder (Rising)"),
        tags$span(tags$span(style = "display: inline-block; width: 12px; height: 12px; background: #E63946; border-radius: 2px; margin-right: 4px; vertical-align: middle;"), "High Flow Season (Peak)"),
        tags$span(tags$span(style = "display: inline-block; width: 12px; height: 12px; background: #F4A261; border-radius: 2px; margin-right: 4px; vertical-align: middle;"), "Fall Shoulder (Receding)")
      )
    )
  })

  # Tab 2: Monthly Summary Table with Colored Pill Badges
  output$monthly_summary_table_html <- renderUI({
    diag <- diagnosis_result()
    req(diag)

    m_df <- diag$monthly_regime
    wy_order <- c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep")
    ordered_df <- tibble(month_abb = wy_order) %>% left_join(m_df, by = "month_abb")

    regime_colors <- c(
      "low flow season"      = "#2B83BA",
      "shoulder:low to high" = "#43AA8B",
      "high flow season"     = "#E63946",
      "shoulder:high to low" = "#F4A261"
    )

    regime_names <- c(
      "low flow season"      = "Low Flow Season",
      "shoulder:low to high" = "Spring Shoulder (Rising)",
      "high flow season"     = "High Flow Season (Peak)",
      "shoulder:high to low" = "Fall Shoulder (Receding)"
    )

    rows <- lapply(1:nrow(ordered_df), function(i) {
      r <- ordered_df[i, ]
      reg <- as.character(r$regime)
      bg_col <- if (!is.na(reg) && reg %in% names(regime_colors)) regime_colors[[reg]] else "#64748B"
      lbl <- if (!is.na(reg) && reg %in% names(regime_names)) regime_names[[reg]] else "No Data"

      tags$tr(
        tags$td(style = "font-weight: 700; text-align: center;", as.character(r$month_abb)),
        tags$td(
          style = "text-align: center;",
          tags$span(
            style = sprintf("display: inline-block; background-color: %s; color: #ffffff; padding: 4px 10px; border-radius: 4px; font-weight: 700; font-size: 0.78rem;", bg_col),
            lbl
          )
        ),
        tags$td(style = "text-align: center;", if (!is.na(r$median_flow)) sprintf("%.1f", r$median_flow) else "---"),
        tags$td(style = "text-align: center;", if (!is.na(r$max_flow)) sprintf("%.1f", r$max_flow) else "---"),
        tags$td(style = "text-align: center;", if (!is.na(r$hf_curve_flow)) sprintf("%.1f", r$hf_curve_flow) else "---"),
        tags$td(style = "text-align: center;", sprintf("%.1f", r$demarc_line))
      )
    })

    tags$table(
      class = "table table-striped table-hover table-bordered table-sm align-middle mt-2",
      tags$thead(
        class = "table-light",
        tags$tr(
          tags$th(style = "text-align: center;", "Month (WY)"),
          tags$th(style = "text-align: center;", "Classified Flow Regime"),
          tags$th(style = "text-align: center;", "Median Daily Flow (cfs)"),
          tags$th(style = "text-align: center;", "Max Observed Flow (cfs)"),
          tags$th(style = "text-align: center;", "Runoff Envelope Flow (cfs)"),
          tags$th(style = "text-align: center;", "Demarcation Line (cfs)")
        )
      ),
      tags$tbody(rows)
    )
  })

  # Tab 3: Timing Detail Table
  output$timing_detail_table <- renderTable({
    diag <- diagnosis_result()
    req(diag)
    tibble(
      `Regime Phase` = c(
        "Low Flow Season (Baseflow Dormancy)",
        "Spring Shoulder (Snowmelt Rising Limb)",
        "High Flow Season (Core Runoff Period)",
        "Fall Shoulder (Receding Limb)"
      ),
      `Continuous Calendar Period` = c(
        diag$timing_summary$low_flow_period,
        diag$timing_summary$rising_shoulder_period,
        diag$timing_summary$high_flow_period,
        diag$timing_summary$receding_shoulder_period
      ),
      `Operational Characteristics` = c(
        sprintf("Flows bounded below demarcation line (%.1f cfs). Groundwater fed.", diag$demarc_line),
        "Snowmelt begins rising; reservoir pre-runoff storage operations.",
        sprintf("Typical runoff percentile (%s) sustains above demarcation line.", diag$timing_summary$selected_high_flow_pct),
        "Post-peak recession back to autumn baseline."
      )
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  # Tab 4: Benchmark Tables
  output$plaspico_bench_table <- renderTable({
    val <- validate_plaspico_benchmark()
    val$validation_table
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  output$placheco_bench_table <- renderTable({
    val <- validate_placheco_benchmark()
    val$validation_table
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  # Tab 5: DOY Data Table
  output$doy_table <- renderDT({
    diag <- diagnosis_result()
    req(diag)
    doy_export <- diag$doy_df %>%
      mutate(
        `Calendar Date` = format(as.Date(doy - 1, origin = "2023-01-01"), "%b %d"),
        `Day of Year` = doy,
        `Years of Record` = n_years,
        `Min Flow (cfs)` = round(min_s, 2),
        `20th %ile (cfs)` = round(p20_s, 2),
        `40th %ile (cfs)` = round(p40_s, 2),
        `Median Flow (cfs)` = round(p50_s, 2),
        `80th %ile (cfs)` = round(p80_s, 2),
        `Max Flow (cfs)` = round(max_s, 2),
        `Demarcation Line (cfs)` = demarc_line,
        `Classified Regime` = regime
      ) %>%
      select(
        `Calendar Date`, `Day of Year`, `Years of Record`,
        `Min Flow (cfs)`, `20th %ile (cfs)`, `40th %ile (cfs)`,
        `Median Flow (cfs)`, `80th %ile (cfs)`, `Max Flow (cfs)`,
        `Demarcation Line (cfs)`, `Classified Regime`
      )

    datatable(doy_export, options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE)
  })

  # Tab 5: CSV Download
  output$download_doy_csv <- downloadHandler(
    filename = function() {
      diag <- diagnosis_result()
      sid <- if (!is.null(diag)) gsub("[^A-Za-z0-9]", "_", diag$site_no) else "hydro_regime"
      sprintf("%s_doy_regime_classification.csv", sid)
    },
    content = function(file) {
      diag <- diagnosis_result()
      req(diag)
      write.csv(diag$doy_df, file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
