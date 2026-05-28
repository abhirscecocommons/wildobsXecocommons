# =============================================================================
# WildObs_occurrence_prep.R
# Self-contained | CamtrapDP v1.0 | Alectura lathami | Mt Isa QLD
#
# ── ECOLOGICAL CONTEXT ────────────────────────────────────────────────────────
# Mt Isa (-20.336385, 137.527452) = WESTERN RANGE MARGIN for Alectura lathami.
# Core range: eastern coastal QLD/NSW. These records are either a genuine
# range extension or require AI-classification validation. Declare in methods.
#
# ── FOUR CASES ────────────────────────────────────────────────────────────────
# Case 1: WildObs site + ALA records (100 km, auto-expands to 500 km)
#         → range-wide SDM with local anchor
# Case 2: WildObs API query for same species across all available projects
#         → multi-site data from broader WildObs network
# Case 3: Single-site dispersal jitter → >= 10 points for MaxEnt/EcoCommons
#         → home range buffer (Göth & Vogel 2002)
# Case 4: Detection/non-detection matrix for occupancy modelling
#         → directly feeds EcoCommons EC_Occupancy_practical notebook
#         → unmarkedFrameOccu() ready: y matrix + siteCovs + obsCovs CSVs
#
# ── FOLDER STRUCTURE ──────────────────────────────────────────────────────────
# wildobsxecocommons/
# ├── wildobs_local_data/          ← your 3 CamtrapDP input CSVs (read-only)
# │   ├── observations.csv
# │   ├── deployments.csv
# │   └── media.csv
# └── outputs/                     ← all script outputs (CSVs + PNGs)
#     ├── occurrences_case1.csv
#     ├── occurrences_case2.csv
#     ├── occurrences_case3.csv
#     ├── occupancy_y_matrix.csv
#     ├── occupancy_site_covs.csv
#     ├── occupancy_obs_covs_effort.csv
#     ├── plot01_detection_counts.png
#     ├── plot02_confidence_histogram.png
#     ├── plot03_case1_ala_map.png
#     ├── plot04_case3_jitter_map.png
#     ├── plot05_occupancy_heatmap.png
#     └── plot06_all_cases_summary.png
#
# ── REFS ──────────────────────────────────────────────────────────────────────
# Phillips et al. (2006) Ecol. Model. 190:231     MaxEnt >= 10 points
# Göth & Vogel (2002) Emu 102:37                  home range / daily movement
# MacKenzie et al. (2002) Ecology 83(8):2248       occupancy model theory
# Guillera-Arroita et al. (2014) MEE 5:914        SDM detection bias
# EcoCommons EC_Occupancy_practical notebook (Feb 2026)
# =============================================================================


# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
# Run install block once, then comment out.

# options(repos = c(CRAN = "https://cran.rstudio.com/"))
# packages_cran <- c("tidyverse","sf","galah","janitor","lubridate",
#                    "ggplot2","scales","unmarked","frictionless")
# for (pkg in packages_cran) {
#   if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
# }
# # WildObsR — install directly from GitHub (no token needed)
# install.packages(
#   "https://github.com/WildObs/WildObsR/archive/refs/heads/main.tar.gz",
#   repos = NULL, type = "source"
# )

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(galah)
  library(janitor)
  library(lubridate)
  library(ggplot2)
  library(scales)
  library(WildObsR)
  library(frictionless)
  library(unmarked)
})


# ── 1. CONFIGURATION ──────────────────────────────────────────────────────────

ALA_EMAIL      <- "support@ecocommons.org.au"
TARGET_SPECIES <- "Alectura lathami"
EXCLUDE_TAXA   <- c("Mammalia", "Aves", "Homo sapiens")

# Case 1 — ALA radius (metres). Auto-expands if no records found.
ALA_RADIUS_M   <- 100000    # start at 100 km
ALA_RADIUS_MAX <- 500000    # expand to 500 km if 100 km returns nothing

# Case 3 — buffer jitter
# 500 m = conservative daily movement radius (Göth & Vogel 2002, Emu 102:37)
BUFFER_M       <- 500
MIN_POINTS     <- 10        # MaxEnt / EcoCommons minimum
SET_SEED       <- 42

# Case 4 — occupancy matrix
OCC_WINDOW_DAYS <- 7        # survey occasion = 7-day window
OCC_MIN_SURVEYS <- 3        # min survey occasions per site (MacKenzie et al. 2002)
OCC_SPECIES     <- TARGET_SPECIES

CONF_THRESHOLD  <- 0.7      # classification confidence flag

# ── Folder paths ──────────────────────────────────────────────────────────────
# Set DATA_DIR to where your 3 CamtrapDP CSVs live.
# Set OUT_DIR  to where all outputs (CSVs + PNGs) will be saved.
# Both folders must exist before running — create them manually first.
# Paths are relative to the WILDOBSXECOCOMMONS folder — run this script
# from inside that folder (i.e. set working directory to WILDOBSXECOCOMMONS).
# Folder names match what is visible in your repo:
#   WILDOBSXECOCOMMONS/wildobs_raw_data/   ← your 3 CamtrapDP CSVs
#   WILDOBSXECOCOMMONS/outputs/            ← all script outputs (CSVs + PNGs)
DATA_DIR <- "wildobs_raw_data"
OUT_DIR  <- "outputs"

# Validate folders exist at startup — fail early with a clear message
if (!dir.exists(DATA_DIR))
  stop("DATA_DIR not found: ", DATA_DIR,
       "\n  Expected: WILDOBSXECOCOMMONS/wildobs_raw_data/ — check working directory is set to WILDOBSXECOCOMMONS.")
if (!dir.exists(OUT_DIR))
  stop("OUT_DIR not found: ", OUT_DIR,
       "\n  Expected: WILDOBSXECOCOMMONS/outputs/ — check working directory is set to WILDOBSXECOCOMMONS.")

# Helper: build output file path
out <- function(filename) file.path(OUT_DIR, filename)

# galah 2.2.0 requires download_reason_id for atlas_occurrences() — without
# it the ALA API returns HTTP 400. Reason 10 = "education/research" (valid
# for all ecological research use). See show_all(reasons) for full list.
galah_config(
  email              = ALA_EMAIL,
  download_reason_id = 10
)

# Load WildObs API key (ships with WildObsR package)
data(wildobsr_api_key)

cat("══════════════════════════════════════════════════════════════\n")
cat("WildObs Occurrence Prep — Alectura lathami\n")
cat("Site: Mt Isa QLD (-20.336385, 137.527452)\n")
cat("⚠  Western range margin — validate AI classifications\n")
cat("══════════════════════════════════════════════════════════════\n\n")


# ── 2. LOAD LOCAL CAMTRAP DP DATA ─────────────────────────────────────────────
cat("── Loading CamtrapDP CSVs ────────────────────────────────────\n")

obs_raw <- read_csv(file.path(DATA_DIR, "observations.csv"), show_col_types = FALSE) |> clean_names()
dep_raw <- read_csv(file.path(DATA_DIR, "deployments.csv"),  show_col_types = FALSE) |> clean_names()
med_raw <- read_csv(file.path(DATA_DIR, "media.csv"),        show_col_types = FALSE) |> clean_names()

cat("observations:", nrow(obs_raw), "| deployments:", nrow(dep_raw),
    "| media:", nrow(med_raw), "\n\n")


# ── 3. VERIFICATION PLOT 1: Species detection summary ─────────────────────────
cat("── Verification Plot 1: Detection counts ─────────────────────\n")

det_summ <- obs_raw |>
  filter(observation_type == "animal") |>
  count(scientific_name, sort = TRUE) |>
  mutate(scientific_name = fct_reorder(scientific_name, n),
         target = scientific_name == TARGET_SPECIES)

p1 <- ggplot(det_summ, aes(scientific_name, n, fill = target)) +
  geom_col() +
  geom_text(aes(label = n), hjust = -0.2, size = 3.2) +
  scale_fill_manual(values = c("FALSE" = "#7FBADC", "TRUE" = "#E84B23"),
                    labels = c("Other species", TARGET_SPECIES)) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Verification 1: Detections by species",
       subtitle = "Red = target | CHECK: counts plausible for Mt Isa?",
       x = NULL, y = "Detections", fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(colour = "grey40", size = 9))
print(p1)
ggsave(out("plot01_detection_counts.png"), p1, width=9, height=6, dpi=150)
cat("CHECK: 38 Alectura lathami at Mt Isa — western range margin.\n\n")


# ── 4. VERIFICATION PLOT 2: Classification confidence ─────────────────────────
cat("── Verification Plot 2: Classification confidence ────────────\n")

conf_data <- obs_raw |>
  filter(observation_type == "animal",
         scientific_name == TARGET_SPECIES,
         !is.na(classification_probability))

p2 <- ggplot(conf_data, aes(classification_probability)) +
  geom_histogram(binwidth = 0.05, fill = "#E84B23", colour = "white", alpha = 0.85) +
  geom_vline(xintercept = CONF_THRESHOLD, linetype = "dashed",
             colour = "navy", linewidth = 0.8) +
  annotate("text", x = CONF_THRESHOLD + 0.02, y = Inf,
           label = paste0("Flag threshold: ", CONF_THRESHOLD),
           vjust = 1.8, hjust = 0, colour = "navy", size = 3.2) +
  labs(title = paste0("Verification 2: Confidence — ", TARGET_SPECIES),
       subtitle = paste0("Records below ", CONF_THRESHOLD,
                         " need manual image review"),
       x = "Classification probability", y = "Count") +
  theme_minimal(base_size = 11)
print(p2)
ggsave(out("plot02_confidence_histogram.png"), p2, width=8, height=5, dpi=150)

n_low <- sum(conf_data$classification_probability < CONF_THRESHOLD, na.rm = TRUE)
cat(sprintf("CHECK: %d / %d records below %.0f%% confidence.\n\n",
            n_low, nrow(conf_data), CONF_THRESHOLD * 100))


# ── 5. BUILD BASE OCCURRENCE TABLE ────────────────────────────────────────────
occ_base <- obs_raw |>
  filter(observation_type == "animal",
         !scientific_name %in% EXCLUDE_TAXA) |>
  left_join(dep_raw |> select(deployment_id, latitude, longitude,
                               location_name, deployment_start, deployment_end),
            by = "deployment_id") |>
  select(species = scientific_name,
         decimal_latitude = latitude, decimal_longitude = longitude,
         location_name, event_start, event_end,
         classification_probability, deployment_id) |>
  mutate(event_start = ymd_hms(event_start, tz = "Australia/Brisbane"),
         event_end   = ymd_hms(event_end,   tz = "Australia/Brisbane"),
         source = "WildObs_CamtrapDP") |>
  drop_na(decimal_latitude, decimal_longitude)

al_cam <- occ_base |> filter(species == TARGET_SPECIES)
site_lat <- unique(al_cam$decimal_latitude)
site_lon <- unique(al_cam$decimal_longitude)

cat(sprintf("Base table: %d records | %d species | %d unique coords\n",
            nrow(occ_base), n_distinct(occ_base$species),
            n_distinct(paste(occ_base$decimal_latitude, occ_base$decimal_longitude))))
cat(sprintf("%s: %d records at (%.4f, %.4f)\n\n",
            TARGET_SPECIES, nrow(al_cam), site_lat, site_lon))


# =============================================================================
# CASE 1: WildObs site + ALA within 100–500 km
# =============================================================================
cat("══ CASE 1: WildObs + ALA (100–500 km) ═══════════════════════\n\n")

cat("Querying ALA for", TARGET_SPECIES, "...\n")
galah_ver <- tryCatch(as.character(packageVersion("galah")), error = function(e) "unknown")
cat("galah version:", galah_ver, "\n")

# galah 2.2.0 confirmed syntax (from galah_apply_profile.R source):
#   - galah_filter() takes field conditions ONLY — profile = "ALA" inside
#     galah_filter() always throws "named input" error in 2.x.
#   - galah_apply_profile("ALA") is a SEPARATE pipe step after galah_filter().
#   - download_reason_id = 10 must be in galah_config() for atlas_occurrences().
ala_raw <- tryCatch({
  galah_call() |>
    galah_identify(TARGET_SPECIES) |>
    galah_filter(
      coordinateUncertaintyInMeters <= 1000,
      year >= 2000
    ) |>
    galah_apply_profile("ALA") |>
    galah_select(group = "basic") |>
    atlas_occurrences()
}, error = function(e) {
  cat("\u26a0 ALA failed:", conditionMessage(e), "\n")
  cat("  Check: galah_config() has your personal ALA-registered email\n")
  cat("  and download_reason_id = 10. Test with atlas_counts() first.\n")
  NULL
})

if (!is.null(ala_raw)) {
  ala_raw <- ala_raw |> clean_names()
  cat("ALA columns:", paste(names(ala_raw), collapse = ", "), "\n")
  cat("ALA records:", nrow(ala_raw), "\n")
  # Post-filter in R as safety net
  if ("coordinate_uncertainty_in_meters" %in% names(ala_raw))
    ala_raw <- filter(ala_raw, is.na(coordinate_uncertainty_in_meters) |
                        coordinate_uncertainty_in_meters <= 1000)
  if ("year" %in% names(ala_raw))
    ala_raw <- filter(ala_raw, is.na(year) | year >= 2000)
}
ala_df <- NULL
ALA_RADIUS_USED <- NA

if (!is.null(ala_raw) && nrow(ala_raw) > 0) {
  site_sf_m <- st_sfc(st_point(c(site_lon, site_lat)), crs = 4326) |>
    st_transform(7854)
  # Column names now snake_case after clean_names():
  # decimalLatitude -> decimal_latitude, decimalLongitude -> decimal_longitude
  ala_sf_m  <- ala_raw |>
    drop_na(decimal_latitude, decimal_longitude) |>
    st_as_sf(coords = c("decimal_longitude", "decimal_latitude"), crs = 4326) |>
    st_transform(7854)

  # Try 100 km first; auto-expand to 500 km if empty
  for (r in c(ALA_RADIUS_M, ALA_RADIUS_MAX)) {
    buf    <- st_buffer(site_sf_m, r)
    inside <- ala_sf_m[st_within(ala_sf_m, buf, sparse = FALSE), ] |>
      st_transform(4326)
    if (nrow(inside) > 0) {
      ALA_RADIUS_USED <- r
      cat(sprintf("ALA records within %d km: %d\n", r / 1000, nrow(inside)))
      break
    }
    cat(sprintf("No ALA records within %d km (expected for range margin) — expanding.\n",
                r / 1000))
  }

  if (!is.na(ALA_RADIUS_USED) && nrow(inside) > 0) {
    # Bug fix: extract coordinates BEFORE st_drop_geometry() — after drop,
    # geometry column no longer exists and st_coordinates() would error.
    coords_inside <- st_coordinates(inside)
    # Column names after clean_names(): scientific_name, event_date
    ala_df <- inside |>
      st_drop_geometry() |>
      mutate(decimal_latitude  = coords_inside[, 2],
             decimal_longitude = coords_inside[, 1],
             species     = scientific_name,
             event_start = as.Date(event_date),
             source      = paste0("ALA_", ALA_RADIUS_USED / 1000, "km")) |>
      select(species, decimal_latitude, decimal_longitude, event_start, source)
  }
}

case1 <- bind_rows(
  al_cam |> select(species, decimal_latitude, decimal_longitude,
                   event_start, source),
  ala_df
) |> distinct(decimal_latitude, decimal_longitude, .keep_all = TRUE)

cat(sprintf("Case 1: %d records | %d unique coords\n", nrow(case1),
            n_distinct(paste(case1$decimal_latitude, case1$decimal_longitude))))

# Verification plot — Case 1 map
p3 <- ggplot(case1, aes(decimal_longitude, decimal_latitude, colour = source)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_point(data = filter(case1, source == "WildObs_CamtrapDP"),
             size = 6, colour = "#E84B23", shape = 18) +
  labs(title = "Verification 3: Case 1 — WildObs + ALA",
       subtitle = paste0("Diamond = WildObs site (Mt Isa) | Dots = ALA within ",
                         ifelse(is.na(ALA_RADIUS_USED), "n/a",
                                paste0(ALA_RADIUS_USED / 1000, "km"))),
       x = "Longitude", y = "Latitude", colour = "Source") +
  theme_minimal(base_size = 11) + coord_equal()
print(p3)
ggsave(out("plot03_case1_ala_map.png"), p3, width=9, height=7, dpi=150)
cat("CHECK: WildObs point isolated from eastern ALA cluster? (Expected)\n\n")

write_csv(
  case1 |> select(species, decimalLatitude = decimal_latitude,
                  decimalLongitude = decimal_longitude),
  out("occurrences_case1.csv"))
message("✓ Exported: occurrences_case1.csv (", nrow(case1), " records)")


# =============================================================================
# CASE 2: WildObs API — multi-site query for same species
# =============================================================================
cat("\n══ CASE 2: WildObs API multi-site query ═════════════════════\n\n")
cat("Querying WildObs MongoDB for projects with", TARGET_SPECIES, "...\n")
cat("Spatial extent: broad Australia bounding box\n\n")

# Query WildObs for any project with Alectura lathami, Australia-wide
wildobs_ids <- tryCatch({
  wildobs_mongo_query(
    api_key  = wildobsr_api_key,
    spatial  = list(xmin = 110.0, xmax = 155.0, ymin = -45.0, ymax = -10.0),
    temporal = list(minDate = as.Date("2018-01-01"),
                    maxDate = as.Date("2025-12-31")),
    tabularSharingPreference = c("open", "partial")
  )
}, error = function(e) {
  cat("⚠ WildObs API query failed:", conditionMessage(e), "\n")
  cat("  Check API key and internet connection.\n")
  NULL
})

# Bug fix 6: wildobs_mongo_query() has no taxonomic filter — it returns ALL
# projects matching spatial/temporal/sharing criteria. Downloading everything
# at once could be huge. Fix: metadata_only=TRUE pass first to identify which
# projects contain the target species, then full download of those only.

case2_placeholder <- function(reason) {
  cat(reason, "\n  Writing local WildObs data as placeholder.\n")
  write_csv(
    al_cam |> select(species, decimalLatitude = decimal_latitude,
                      decimalLongitude = decimal_longitude),
    out("occurrences_case2.csv"))
  message("✓ Exported: occurrences_case2.csv (placeholder)")
}

if (!is.null(wildobs_ids) && length(wildobs_ids) > 0) {
  cat("Projects found:", length(wildobs_ids), "\n")
  cat("Pass 1: metadata-only download to identify species-relevant projects...\n")

  # Pass 1: lightweight metadata scan to find which projects have target species
  meta_list <- tryCatch({
    wildobs_dp_download(
      api_key       = wildobsr_api_key,
      project_ids   = wildobs_ids,
      media         = FALSE,
      metadata_only = TRUE    # fast — no tabular data
    )
  }, error = function(e) {
    cat("⚠ Metadata download failed:", conditionMessage(e), "\n"); NULL
  })

  # Filter to projects whose taxonomic metadata mentions target species
  if (!is.null(meta_list) && length(meta_list) > 0) {
    relevant_ids <- keep(names(meta_list), function(nm) {
      dp <- meta_list[[nm]]
      # CamtrapDP metadata: taxonomic list is in dp$taxonomic
      taxa <- tryCatch(dp$taxonomic, error = function(e) NULL)
      if (is.null(taxa)) return(FALSE)
      any(grepl(TARGET_SPECIES, unlist(taxa), ignore.case = TRUE))
    })
    cat(sprintf("Projects containing %s: %d / %d\n",
                TARGET_SPECIES, length(relevant_ids), length(wildobs_ids)))
  } else {
    relevant_ids <- character(0)
  }

  if (length(relevant_ids) == 0) {
    case2_placeholder(paste0("⚠ No WildObs projects found with ", TARGET_SPECIES, "."))
  } else {
    cat("Pass 2: full download of", length(relevant_ids), "relevant project(s)...\n")

    dp_list <- tryCatch({
      wildobs_dp_download(
        api_key     = wildobsr_api_key,
        project_ids = relevant_ids,
        media       = FALSE,
        metadata_only = FALSE
      )
    }, error = function(e) {
      cat("⚠ Full download failed:", conditionMessage(e), "\n"); NULL
    })

    if (is.null(dp_list) || length(dp_list) == 0) {
      case2_placeholder("⚠ Download returned empty list.")
    } else {
      cat("Data packages downloaded:", length(dp_list), "\n")

      # Extract obs + deployments from all packages, filter to target species
      all_obs_list <- map(dp_list, function(dp) {
        tryCatch({
          obs_dp <- frictionless::read_resource(dp, "observations") |> clean_names()
          dep_dp <- frictionless::read_resource(dp, "deployments")  |> clean_names()
          obs_dp |>
            filter(observation_type == "animal",
                   scientific_name == TARGET_SPECIES) |>
            left_join(dep_dp |> select(deployment_id, latitude, longitude,
                                        location_name),
                      by = "deployment_id") |>
            select(species = scientific_name,
                   decimal_latitude = latitude, decimal_longitude = longitude,
                   location_name, event_start) |>
            drop_na(decimal_latitude, decimal_longitude)
        }, error = function(e) NULL)
      }) |> compact()

      if (length(all_obs_list) == 0) {
        case2_placeholder("⚠ No animal records extracted from downloaded packages.")
      } else {
        case2_api <- bind_rows(all_obs_list) |>
          mutate(source = "WildObs_API") |>
          distinct(decimal_latitude, decimal_longitude, event_start, .keep_all = TRUE)

        n_sites_api <- n_distinct(paste(case2_api$decimal_latitude,
                                         case2_api$decimal_longitude))
        cat(sprintf("%s from WildObs API: %d records | %d sites\n",
                    TARGET_SPECIES, nrow(case2_api), n_sites_api))

        # Verification plot — Case 2 API map
        p_c2 <- ggplot(case2_api, aes(decimal_longitude, decimal_latitude)) +
          geom_point(colour = "#9B59B6", size = 3, alpha = 0.75) +
          labs(title = "Verification: Case 2 — WildObs API multi-site",
               subtitle = paste0(n_sites_api, " unique sites | ",
                                 nrow(case2_api), " records"),
               x = "Longitude", y = "Latitude") +
          theme_minimal(base_size = 11) + coord_equal()
        print(p_c2)
ggsave(out("plot_case2_multisite.png"), p_c2, width=9, height=7, dpi=150)

        case2_out <- case2_api |>
          select(species, decimalLatitude = decimal_latitude,
                 decimalLongitude = decimal_longitude)
        write_csv(case2_out, out("occurrences_case2.csv"))
        message("✓ Exported: occurrences_case2.csv (", nrow(case2_out),
                " records, ", n_sites_api, " sites)")
      }
    }
  }
} else {
  case2_placeholder("⚠ No WildObs projects returned for Australia-wide query.")
}


# =============================================================================
# CASE 3: Dispersal jitter — >= 10 unique points
# =============================================================================
cat("\n══ CASE 3: Dispersal jitter (>= 10 points) ══════════════════\n\n")
cat(sprintf("Buffer: %dm | Justification: daily movement, Göth & Vogel (2002)\n",
            BUFFER_M))

site_sf_c3 <- st_sfc(st_point(c(site_lon, site_lat)), crs = 4326) |>
  st_sf(species = TARGET_SPECIES)

buffer_c3 <- site_sf_c3 |>
  st_transform(7854) |>
  st_buffer(BUFFER_M) |>
  st_transform(4326)

n_gen <- max(MIN_POINTS, nrow(al_cam))
set.seed(SET_SEED)

# Bug fix: st_coordinates(geometry) inside mutate() fails — dplyr evaluates
# column names as symbols, not geometry objects, so 'geometry' is not found.
# Fix: extract coordinates into a matrix BEFORE mutate(), reference by object.
jitter_sf     <- st_sample(buffer_c3, size = n_gen, type = "random") |> st_sf()
jitter_coords <- st_coordinates(jitter_sf)   # extract before mutate/drop

jitter_pts <- jitter_sf |>
  mutate(species           = TARGET_SPECIES,
         decimal_latitude  = jitter_coords[, 2],
         decimal_longitude = jitter_coords[, 1],
         source            = paste0("jitter_", BUFFER_M, "m")) |>
  st_drop_geometry()

case3 <- bind_rows(
  al_cam |> select(species, decimal_latitude, decimal_longitude) |>
    mutate(source = "WildObs_CamtrapDP"),
  jitter_pts |> select(species, decimal_latitude, decimal_longitude, source)
)

# ALA fallback if still < MIN_POINTS
n_uniq_c3 <- n_distinct(paste(case3$decimal_latitude, case3$decimal_longitude))
if (n_uniq_c3 < MIN_POINTS && !is.null(ala_df) && nrow(ala_df) > 0) {
  case3 <- bind_rows(
    case3,
    ala_df |> mutate(source = "ALA_supplement") |>
      select(species, decimal_latitude, decimal_longitude, source)
  ) |> distinct(decimal_latitude, decimal_longitude, .keep_all = TRUE)
  cat("ALA supplement added. Total:", nrow(case3), "\n")
}

cat(sprintf("Case 3: %d records | %d unique coords\n", nrow(case3),
            n_distinct(paste(case3$decimal_latitude, case3$decimal_longitude))))

# Verification plot — buffer map
buf_coords <- buffer_c3 |> st_coordinates() |>
  as.data.frame() |> rename(lon = X, lat = Y)

p4 <- ggplot() +
  geom_polygon(data = buf_coords, aes(lon, lat),
               fill = "#FFD700", alpha = 0.2, colour = "#FFD700", linewidth = 0.8) +
  geom_point(data = case3,
             aes(decimal_longitude, decimal_latitude, colour = source,
                 size = source == "WildObs_CamtrapDP"), alpha = 0.85) +
  scale_colour_manual(values = c("WildObs_CamtrapDP" = "#E84B23",
                                  "jitter_500m"        = "#2ECC71",
                                  "ALA_supplement"     = "#3A86FF")) +
  scale_size_manual(values = c("TRUE" = 5, "FALSE" = 2.5), guide = "none") +
  labs(title = "Verification 4: Case 3 — dispersal jitter",
       subtitle = paste0(n_gen, " points within ", BUFFER_M,
                         "m buffer | Diamond = original site"),
       x = "Longitude", y = "Latitude", colour = "Source") +
  theme_minimal(base_size = 11) + coord_equal()
print(p4)
ggsave(out("plot04_case3_jitter_map.png"), p4, width=8, height=7, dpi=150)
cat("CHECK: All jittered points inside yellow buffer?\n")
cat("CHECK: 500m ecologically justified (daily movement, Göth & Vogel 2002)\n\n")

write_csv(
  case3 |> select(species, decimalLatitude = decimal_latitude,
                  decimalLongitude = decimal_longitude),
  out("occurrences_case3.csv"))
message("✓ Exported: occurrences_case3.csv (", nrow(case3), " records)")


# =============================================================================
# CASE 4: Detection/non-detection matrix for EcoCommons Occupancy notebook
# =============================================================================
# This section prepares data directly for:
#   EC_Occupancy_practical.html (EcoCommons/WildObs notebook, Feb 2026)
# Output feeds unmarkedFrameOccu():
#   y        = occupancy_y_matrix.csv   (sites × occasions, 1/0/NA)
#   siteCovs = occupancy_site_covs.csv  (1 row per site)
#   obsCovs  = occupancy_obs_covs_effort.csv (sites × occasions, effort in days)
#
# SINGLE-SEASON ASSUMPTIONS (MacKenzie et al. 2002 Ecology 83:2248):
#   - Closure: occupancy state constant within season (no colonisation/extinction)
#   - Sites spatially independent (distance > home range)
#   - Detections independent across occasions conditional on occupancy
#   - No false positives
#
# CURRENT DATA LIMITATION:
#   Your data has 1 unique lat/lon across 15 deployments — all are the same
#   spatial location. Occupancy model needs spatially independent sites.
#   This section builds the matrix structure correctly from your data and
#   flags the limitation explicitly. It will produce a valid matrix when
#   multi-site data is loaded (via Case 2 WildObs API or new deployments).
# =============================================================================
cat("\n══ CASE 4: Occupancy matrix for EcoCommons notebook ═════════\n\n")
cat("Target: EC_Occupancy_practical notebook (Feb 2026)\n")
cat("Method: detection/non-detection per site per 7-day window\n")
cat("Output: y matrix + siteCovs + obsCovs CSVs for unmarkedFrameOccu()\n\n")

# Use all animal observations for Case 4 (not just target species)
# suppressMessages() wraps ymd_hms() calls to silence ISO8601 timezone
# conversion messages — these are informational only, not errors.
occ_all <- obs_raw |>
  filter(observation_type == "animal") |>
  left_join(dep_raw |> select(deployment_id, latitude, longitude,
                               location_name, deployment_start, deployment_end),
            by = "deployment_id") |>
  mutate(
    event_start      = suppressMessages(ymd_hms(event_start,      tz = "Australia/Brisbane")),
    deployment_start = suppressMessages(ymd_hms(deployment_start, tz = "Australia/Brisbane")),
    deployment_end   = suppressMessages(ymd_hms(deployment_end,   tz = "Australia/Brisbane")),
    detected         = as.integer(scientific_name == OCC_SPECIES)
  ) |>
  drop_na(latitude, longitude, event_start)

# ── Step 1: Define sites ─────────────────────────────────────────────────────
sites <- dep_raw |>
  select(site_id = deployment_id, latitude, longitude,
         location_name, deployment_start, deployment_end) |>
  mutate(
    deployment_start = suppressMessages(ymd_hms(deployment_start, tz = "Australia/Brisbane")),
    deployment_end   = suppressMessages(ymd_hms(deployment_end,   tz = "Australia/Brisbane")),
    survey_days      = as.numeric(difftime(deployment_end, deployment_start, units = "days"))
  )

# Compute n_occ safely using base R ifelse (not dplyr if_else) — avoids strict
# type checking that causes if_else to propagate NA when survey_days is NA.
# Case 1: NA survey_days (malformed date row 11) → 0 occasions
# Case 2: survey_days = 0 (start == end, rows 1/12/13) → 0 occasions
# Case 3: survey_days < OCC_WINDOW_DAYS (rows 3/7, 4 days) → 0 occasions
sites$n_occ <- ifelse(
  is.na(sites$survey_days) | sites$survey_days < OCC_WINDOW_DAYS,
  0L,
  as.integer(floor(sites$survey_days / OCC_WINDOW_DAYS))
)

n_sites <- nrow(sites)
max_occ <- max(sites$n_occ, 1L)

cat(sprintf("Sites (deployments)              : %d\n", n_sites))
cat(sprintf("Unique lat/lon pairs             : %d\n",
            n_distinct(paste(sites$latitude, sites$longitude))))
cat(sprintf("Survey occasion window           : %d days\n", OCC_WINDOW_DAYS))
cat(sprintf("Sites with 0 occasions (excluded): %d\n", sum(sites$n_occ == 0L)))
cat(sprintf("Sites with >= 1 occasion         : %d\n", sum(sites$n_occ > 0L)))
cat(sprintf("Max occasions per site           : %d\n", max_occ))
cat("\nSite summary:\n")
print(sites[, c("site_id", "survey_days", "n_occ")], n = Inf)

# ── Step 2: Build y matrix and effort matrix ──────────────────────────────────
occ_rows    <- vector("list", n_sites)
effort_rows <- vector("list", n_sites)

for (i in seq_len(n_sites)) {
  # Use [[ for safe scalar extraction from tibble row — avoids 1-length vector
  # issues when comparing with ==
  n_occ   <- sites$n_occ[[i]]
  site_id <- sites$site_id[[i]]

  if (is.na(n_occ) || n_occ == 0L) {
    occ_rows[[i]]    <- rep(NA_integer_, max_occ)
    effort_rows[[i]] <- rep(NA_real_,    max_occ)
    next
  }

  start_d <- as.Date(sites$deployment_start[[i]])

  det_vec <- map_int(seq_len(n_occ), function(k) {
    win_start <- start_d + (k - 1L) * OCC_WINDOW_DAYS
    win_end   <- win_start + OCC_WINDOW_DAYS
    det <- occ_all |>
      filter(deployment_id == site_id,
             as.Date(event_start) >= win_start,
             as.Date(event_start) <  win_end) |>
      pull(detected)
    if (length(det) == 0L) return(0L)
    return(as.integer(any(det == 1L)))
  })

  occ_rows[[i]]    <- c(det_vec,    rep(NA_integer_, max_occ - n_occ))
  effort_rows[[i]] <- c(rep(as.numeric(OCC_WINDOW_DAYS), n_occ),
                        rep(NA_real_, max_occ - n_occ))
}

# Bug fix 4: do.call(rbind, args = _) — the _ placeholder requires R >=4.2 and
# a named pipe; plain do.call(rbind, list) is safer across R versions.
y_matrix      <- do.call(rbind, occ_rows)    |> as.data.frame()
effort_matrix <- do.call(rbind, effort_rows) |> as.data.frame()

col_names <- paste0("occ_", seq_len(max_occ))
colnames(y_matrix)      <- col_names
colnames(effort_matrix) <- col_names

# Bug fix 5: rownames after do.call(rbind) on a list of vectors are NULL.
# Set explicitly from sites$site_id.
rownames(y_matrix)      <- sites$site_id
rownames(effort_matrix) <- sites$site_id

cat(sprintf("y matrix: %d sites × %d occasions\n", nrow(y_matrix), ncol(y_matrix)))
cat(sprintf("Detections (1s) : %d\n", sum(y_matrix == 1L, na.rm = TRUE)))
cat(sprintf("Non-NA cells    : %d\n", sum(!is.na(y_matrix))))
cat(sprintf("NA cells (no survey): %d\n", sum(is.na(y_matrix))))

naive_occ <- mean(rowSums(y_matrix, na.rm = TRUE) > 0)
cat(sprintf("Naive occupancy (>= 1 detection): %.1f%%\n\n", naive_occ * 100))
# ── Step 3: Site covariates (siteCovs) ───────────────────────────────────────
# One row per site. These go into ~habitat in occu(~p ~psi).
# Available from your data: lat, lon, survey_days.
# Ideal additions: elevation, habitat type, IBRA region (via WildObsR).

site_covs <- sites |>
  select(site_id, latitude, longitude, location_name, survey_days) |>
  mutate(
    # Placeholder habitat — replace with real NVIS/habitat data
    # Use WildObsR::ibra_classification() for IBRA bioregion
    # Use WildObsR::locationName_buffer_CAPAD() for protected area name
    habitat_placeholder = "unknown",
    note = "Add real habitat covariates via WildObsR or raster extract"
  )

# ── Step 4: Observation-level covariates (obsCovs) ───────────────────────────
# One value per site × occasion. Survey effort in days per window.
# All windows = OCC_WINDOW_DAYS except possibly the last (truncate at end).

# Effort matrix — same fixes as y_matrix above.
# Uses pre-computed sites$n_occ so no NULL rows possible.
effort_rows <- map(seq_len(n_sites), function(i) {
  n_occ <- sites$n_occ[i]
  if (n_occ == 0L) return(rep(NA_real_, max_occ))
  c(rep(as.numeric(OCC_WINDOW_DAYS), n_occ),
    rep(NA_real_, max_occ - n_occ))
})

effort_matrix <- do.call(rbind, effort_rows) |> as.data.frame()
colnames(effort_matrix) <- paste0("occ_", seq_len(max_occ))
rownames(effort_matrix) <- sites$site_id

# ── Verification Plot 5: Detection matrix heatmap ────────────────────────────
cat("── Verification Plot 5: Detection history heatmap ────────────\n")

y_long <- y_matrix |>
  rownames_to_column("site") |>
  pivot_longer(-site, names_to = "occasion", values_to = "detection") |>
  mutate(occasion = factor(occasion, levels = paste0("occ_", seq_len(max_occ))))

p5 <- ggplot(y_long, aes(occasion, site, fill = factor(detection))) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = c("0" = "#BDC3C7", "1" = "#E84B23", "NA" = "#F0F0F0"),
    na.value = "#F0F0F0",
    labels  = c("0" = "Not detected", "1" = "Detected", "NA" = "Not surveyed"),
    na.translate = FALSE
  ) +
  labs(title = paste0("Verification 5: Detection history — ", OCC_SPECIES),
       subtitle = paste0(nrow(y_matrix), " sites × ", ncol(y_matrix),
                         " occasions (", OCC_WINDOW_DAYS, "-day windows)"),
       x = "Survey occasion", y = "Site (deploymentID)",
       fill = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        plot.subtitle = element_text(colour = "grey40", size = 9))
print(p5)
ggsave(out("plot05_occupancy_heatmap.png"), p5, width=11, height=7, dpi=150)

cat("CHECK: Are there enough detections (1s) across occasions?\n")
cat("CHECK: Is there spatial variation across sites?\n")
cat(sprintf("⚠ WARNING: All %d sites share 1 lat/lon — violates spatial\n", n_sites))
cat("  independence assumption of occupancy models.\n")
cat("  Load multi-site data (Case 2 API) before running occu().\n\n")

# ── Step 5: Export all 3 occupancy CSVs ──────────────────────────────────────

# y matrix — add site_id column
y_export <- y_matrix |> rownames_to_column("site_id")
write_csv(y_export, out("occupancy_y_matrix.csv"))
message("✓ Exported: occupancy_y_matrix.csv (",
        nrow(y_matrix), " sites × ", ncol(y_matrix), " occasions)")

write_csv(site_covs, out("occupancy_site_covs.csv"))
message("✓ Exported: occupancy_site_covs.csv (", nrow(site_covs), " sites)")

effort_export <- effort_matrix |> rownames_to_column("site_id")
write_csv(effort_export, out("occupancy_obs_covs_effort.csv"))
message("✓ Exported: occupancy_obs_covs_effort.csv (effort per occasion)")

# ── Step 6: Show exactly how to feed these into EcoCommons notebook ──────────
cat("\n── HOW TO USE IN EC_Occupancy_practical NOTEBOOK ─────────────\n")
cat("
# Read your prepared files:
y_mat  <- read_csv(file.path('outputs','occupancy_y_matrix.csv')) |>
  column_to_rownames('site_id') |> as.matrix()

s_covs <- read_csv(file.path('outputs','occupancy_site_covs.csv')) |>
  column_to_rownames('site_id') |> as.data.frame()

eff    <- read_csv(file.path('outputs','occupancy_obs_covs_effort.csv')) |>
  column_to_rownames('site_id') |> as.data.frame()

# Build unmarkedFrameOccu (matches EC_Occupancy_practical notebook):
umf <- unmarkedFrameOccu(
  y        = y_mat,
  siteCovs = s_covs,           # latitude, longitude, survey_days
  obsCovs  = list(effort = eff) # effort in days per window
)

# Null model (intercept only — start here):
fit_null <- occu(~ 1 ~ 1, data = umf)
summary(fit_null)

# With site covariates (e.g. survey_days as proxy for effort):
fit_cov <- occu(~ effort ~ survey_days, data = umf)

# Model comparison with AIC:
library(MuMIn)
model.sel(fit_null, fit_cov)
")
cat("── Reference: MacKenzie et al. (2002) Ecology 83(8):2248 ──────\n\n")


# ── 6. ALL-CASES SUMMARY PLOT ─────────────────────────────────────────────────
cat("── Verification Plot 6: All cases summary ────────────────────\n")

all_cases <- bind_rows(
  read_csv(out("occurrences_case1.csv"), show_col_types = FALSE) |>
    mutate(case = "Case 1\nWildObs + ALA"),
  read_csv(out("occurrences_case2.csv"), show_col_types = FALSE) |>
    select(species, decimalLatitude, decimalLongitude) |>
    mutate(case = "Case 2\nWildObs API"),
  read_csv(out("occurrences_case3.csv"), show_col_types = FALSE) |>
    mutate(case = "Case 3\nJitter (500m)")
) |> filter(!is.na(decimalLatitude))

p6 <- ggplot(all_cases, aes(decimalLongitude, decimalLatitude)) +
  geom_point(colour = "#E84B23", alpha = 0.55, size = 1.8) +
  facet_wrap(~case, scales = "free") +
  labs(title = paste0("Verification 6: All cases — ", TARGET_SPECIES),
       subtitle = "Check: spatial spread appropriate for intended SDM?",
       x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 10) +
  theme(strip.text    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40", size = 9))
print(p6)
ggsave(out("plot06_all_cases_summary.png"), p6, width=12, height=5, dpi=150)


# ── 7. FINAL SUMMARY ──────────────────────────────────────────────────────────
cat("\n══ EXPORT SUMMARY ════════════════════════════════════════════\n")
cat(sprintf("%-42s %5s  %s\n", "File", "Recs", "Unique coords / notes"))
cat(strrep("─", 70), "\n")

for (f in c("occurrences_case1.csv", "occurrences_case2.csv",
             out("occurrences_case3.csv"))) {
  d <- read_csv(f, show_col_types = FALSE) |> filter(!is.na(decimalLatitude))
  cat(sprintf("%-42s %5d  %d unique coords\n", f, nrow(d),
              n_distinct(paste(d$decimalLatitude, d$decimalLongitude))))
}

cat(sprintf("%-42s %5d  sites × %d occasions\n",
            "occupancy_y_matrix.csv", nrow(y_matrix), ncol(y_matrix)))
cat(sprintf("%-42s %5d  rows (1 per site)\n",
            "occupancy_site_covs.csv", nrow(site_covs)))
cat(sprintf("%-42s %5d  sites × %d occasions\n",
            "occupancy_obs_covs_effort.csv", nrow(effort_matrix), ncol(effort_matrix)))

cat(strrep("─", 70), "\n")
cat("\nPlots saved to:", OUT_DIR, "\n")
pngs <- list.files(OUT_DIR, pattern="\\.png$", full.names=FALSE)
if (length(pngs) > 0) {
  for (p in pngs) cat(sprintf("  %-45s\n", p))
} else {
  cat("  (no PNGs found — check OUT_DIR)\n")
}
cat(strrep("─", 70), "\n")
cat("\nSDM format:   species | decimalLatitude | decimalLongitude\n")
cat("Occupancy:    y_matrix + site_covs + obs_covs → unmarkedFrameOccu()\n")
cat("\n⚠ KEY LIMITATIONS:\n")
cat("  Case 2: placeholder if WildObs API returns no additional sites\n")
cat("  Case 3: jittered points = modelled estimates — declare in methods\n")
cat("  Case 4: current data = 1 lat/lon; occupancy model needs multi-site\n")
cat("          Load multi-site data (Case 2) to meet MacKenzie et al. assumptions\n")
cat("══════════════════════════════════════════════════════════════\n")