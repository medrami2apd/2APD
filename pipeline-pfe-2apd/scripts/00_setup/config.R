CONFIG <- list(
  P_LAGS        = 1L,
  HORIZON_GFEVD = 10L,
  HORIZON_GIRF  = 20L,
  RESIST_HORIZON = 4L,

  SECTORS = c("AGR","MIN","MFG","UTL","CON","TRD","HOS",
              "TRN","ICT","FIN","PUB","EDH","REA","OTS"),

  SECTOR_LABELS = c(
    AGR="Agriculture", MIN="Mining", MFG="Manufacturing", UTL="Utilities",
    CON="Construction", TRD="Trade", HOS="Hospitality", TRN="Transport",
    ICT="ICT", FIN="Finance", PUB="Public Admin", EDH="Education & Health",
    REA="Real Estate", OTS="Other Services"),

  N_CHAINS      = 8L,
  N_ITER        = 12000L,
  N_BURN        = 6000L,
  THIN          = 6L,
  N_TRAIN       = 36L,
  SEED          = 12345L,

  KAPPA_Q       = 0.01,
  KAPPA_Q_GRID  = c(0.003, 0.01, 0.03, 0.1),
  P0_SCALE      = 10,
  SV_MU_PRIOR_V = 10,
  SV_PHI_A      = 20,
  SV_PHI_B      = 1.5,
  SV_S2H_C0     = 2.5,
  SV_S2H_D0     = 0.025,

  UNCERTAINTY_MODE = "draws",
  N_DRAWS_PROP     = 500L,
  HAC_LAG          = NULL,
  PANEL_CLUSTER    = "group",
  MULTIPLE_TEST    = "BH",
  FDR_ALPHA        = 0.05,
  H3_SMALL_G_CORRECTION = TRUE,
  H3_WCB_REPS      = 999L,

  NETWORK_THRESHOLD_TYPE = "proportional",
  NETWORK_THRESHOLD      = 0.30,
  NETWORK_THRESHOLD_GRID = c(0.10, 0.20, 0.30, 0.50),

  ROLE_NEUTRAL_BAND = 0.0,
  BOOT_B            = 1000L,
  BLOCK_LEN         = NULL,

  RESIST_EPS        = 1e-2,
  RESIST_EPS_GRID   = c(1e-3, 1e-2, 1e-1),

  BREAK_MAX        = 5L,
  BREAK_MIN_SIZE   = 0.10,
  BREAK_ALIGN_TOL  = 1L,

  ROBUSTNESS = list(
    ROLLING_WINDOW  = 40L,
    WAIC_MAX_DRAWS  = 200L,
    GEWEKE_FRAC1    = 0.1,
    GEWEKE_FRAC2    = 0.5,
    RHAT_THRESHOLD  = 1.01,
    ESS_THRESHOLD   = 400L,
    PIT_BINS        = 10L,
    DERIVED_TCI_DRAWS_PER_CHAIN = 40L,
    DERIVED_BETA_TIMEPOINTS     = 3L,
    TRACE_N_EXTREME_SV = 3L
  ),

  CREDIBLE_BANDS = list(
    N_DRAWS        = 300L,
    QUANTILE_PROBS = c(0.025, 0.5, 0.975)
  ),

  EPS = 1e-8,

  SUBSET_START = NULL,
  SUBSET_END   = NULL,

  SENSITIVITY = list(
    N_ITER  = 5000L,
    N_BURN  = 2500L,
    THIN    = 4L,
    N_DRAWS = 200L
  ),

  PATHS = list(
    raw_data   = "data/raw/sectoral_growth.csv",
    processed  = "outputs/processed/Y_processed.rds",
    posterior       = "outputs/posterior/posterior_mean_arrays.rds",
    posterior_draws = "outputs/posterior/posterior_draws.rds",
    fit             = "outputs/posterior/tvpvar_fit.rds",
    gfevd      = "outputs/connectedness/gfevd_results.rds",
    network    = "outputs/network/network_results.rds",
    resilience = "outputs/resilience/resilience_results.rds",
    breaks     = "outputs/diagnostics/structural_breaks.rds",
    draws_inference = "outputs/hypotheses/draws_inference.rds",
    credible_bands  = "outputs/credible_bands/credible_bands.rds",
    robustness = "outputs/diagnostics/robustness_evaluation.rds",
    hypotheses = "outputs/hypotheses",
    sensitivity= "outputs/sensitivity",
    tables     = "outputs/tables",
    figures    = "outputs/figures",
    logs       = "outputs/logs"
  ),

  HISTORICAL_EPISODES = c(
    "1983"="Structural Adjustment Program",
    "1996"="EU Association Agreement",
    "2004"="US-Morocco Free Trade Agreement",
    "2005"="Plan Emergence (industrial policy)",
    "2008"="Global Financial Crisis",
    "2014"="Industrial Acceleration Plan",
    "2015"="Samir Refinery Closure",
    "2020"="COVID-19 pandemic",
    "2022"="Energy & food price shock"),

  PREREGISTERED_EPISODES = c(1996L, 2005L, 2008L, 2015L, 2020L)
)

CONFIG$K  <- length(CONFIG$SECTORS)
CONFIG$NX <- CONFIG$K * CONFIG$P_LAGS

if (exists(".CONFIG_OVERRIDE", inherits = TRUE) && is.list(.CONFIG_OVERRIDE)) {
  for (.nm in names(.CONFIG_OVERRIDE)) CONFIG[[.nm]] <- .CONFIG_OVERRIDE[[.nm]]
  CONFIG$K  <- length(CONFIG$SECTORS)
  CONFIG$NX <- CONFIG$K * CONFIG$P_LAGS
}

invisible(CONFIG)
