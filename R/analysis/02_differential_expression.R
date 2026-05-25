# =============================================================================
# 02_differential_expression.R
# Proteomics Analysis — 3xTg-AD Mouse Model
#
# Description:
#   Performs statistical comparison (3xTg vs WT) protein-by-protein.
#   Selects Welch's t-test or Wilcoxon test based on Shapiro-Wilk normality.
#   Computes log2 fold change and -log10 p-value for volcano plot.
#   Classifies proteins as up-regulated, down-regulated, or unchanged.
#
# Input:
#   - data/processed/ajustado.xlsx       (normalized abundance matrix)
#   - data/processed/entera_3xtg.xlsx    (filtered full protein table)
#   - data/raw/sample-metadata.xlsx      (sample annotations)
#
# Output:
#   - results/tables/vulcanot_3xtg.xlsx          (all proteins with stats)
#   - results/tables/df_desreguladas_3xtg.xlsx   (significant DE proteins)
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(readxl)
library(openxlsx)
library(dplyr)
library(tibble)
library(effsize)

options(scipen = 999)

# --- Load data ----------------------------------------------------------------
entera_experimento_3xtg <- read_excel("data/processed/entera_3xtg.xlsx")
sample_metadata          <- read_xlsx("data/raw/sample-metadata.xlsx")
ajustado_raw             <- read_xlsx("data/processed/ajustado.xlsx")

# Reconstruct row names from Accession column
accesion_col    <- dplyr::select(entera_experimento_3xtg, Accession)
gene_symbol_col <- dplyr::select(entera_experimento_3xtg, `Gene Symbol`)
indice_proteina <- cbind(accesion_col, gene_symbol_col)

ajustado <- as.data.frame(ajustado_raw)
rownames(ajustado) <- ajustado_raw$...1   # first column holds row names after read
ajustado$...1 <- NULL

# --- Separate samples by genotype ---------------------------------------------
prot_control <- sample_metadata %>% filter(Genotype == "WT")    %>% pull(`Sample-id`)
prot_3xtg    <- sample_metadata %>% filter(Genotype == "3xTg")  %>% pull(`Sample-id`)

tabla_control <- ajustado %>% select(all_of(prot_control))
tabla_3xtg    <- ajustado %>% select(all_of(prot_3xtg))

# Group means
tabla_control_datos       <- tabla_control
tabla_control_datos$media <- rowSums(tabla_control, na.rm = TRUE) / ncol(tabla_control)

tabla_3xtg_datos       <- tabla_3xtg
tabla_3xtg_datos$media <- rowSums(tabla_3xtg, na.rm = TRUE) / ncol(tabla_3xtg)

nprot <- nrow(tabla_control)

# --- Statistical testing (protein-by-protein) ---------------------------------
p_3xtg          <- numeric(nprot)
normalidad_3xtg <- numeric(nprot)
normalidad_ctrl <- numeric(nprot)

for (n in seq_len(nprot)) {
  c_3xtg   <- as.numeric(tabla_3xtg[n, ])
  c_control <- as.numeric(tabla_control[n, ])

  # Shapiro-Wilk normality test (handle constant vectors)
  p_norm_3xtg  <- if (length(unique(c_3xtg))   == 1) 1 else shapiro.test(c_3xtg)$p.value
  p_norm_ctrl  <- if (length(unique(c_control)) == 1) 1 else shapiro.test(c_control)$p.value

  normalidad_3xtg[n] <- p_norm_3xtg
  normalidad_ctrl[n] <- p_norm_ctrl

  # Choose test: t-test if both normal, Wilcoxon otherwise
  if (all(c_3xtg == c_control)) {
    p_3xtg[n] <- 1
  } else if (p_norm_3xtg > 0.05 && p_norm_ctrl > 0.05) {
    p_3xtg[n] <- t.test(c_3xtg, c_control, var.equal = FALSE)$p.value
  } else {
    p_3xtg[n] <- wilcox.test(c_3xtg, c_control)$p.value
  }
}

# --- Assemble results table ---------------------------------------------------
gene_symbols <- indice_proteina$`Gene Symbol`[
  match(rownames(tabla_3xtg), indice_proteina$Accession)
]

result_3xtg <- data.frame(
  Protein           = rownames(tabla_3xtg),
  Gene              = gene_symbols,
  P_Value           = p_3xtg,
  mean_control      = tabla_control_datos$media,
  mean_3xtg         = tabla_3xtg_datos$media,
  normalidad_3xtg   = normalidad_3xtg,
  normalidad_control = normalidad_ctrl
)

# --- Compute fold change and -log10 p-value -----------------------------------
log2FC   <- log2(tabla_3xtg_datos$media / tabla_control_datos$media)
neg_log10_pval <- -log10(result_3xtg$P_Value)

descriptions_selected <- entera_experimento_3xtg$Description[
  match(result_3xtg$Protein, entera_experimento_3xtg$Accession)
]

vulcanot_3xtg <- data.frame(
  Accession          = result_3xtg$Protein,
  Gene.Symbol        = gene_symbols,
  log2cociente_3xtg  = log2FC,
  log10pvalor_3xtg   = neg_log10_pval,
  P_Value            = result_3xtg$P_Value,
  Normality_3xtg     = normalidad_3xtg,
  Normality_control  = normalidad_ctrl,
  mean_3xtg          = tabla_3xtg_datos$media,
  mean_control       = tabla_control_datos$media,
  Description        = descriptions_selected
)

# --- Classification thresholds ------------------------------------------------
p_umbral <- -log10(0.05)   # ≈ 1.301
FC       <- 0.25
LFCP     <- log2(1 + FC)   # upper FC threshold  (~0.322)
LFCN     <- log2(1 - FC)   # lower FC threshold  (~-0.415)

up_regulated   <- filter(vulcanot_3xtg, log2cociente_3xtg >  LFCP & log10pvalor_3xtg > p_umbral)
down_regulated <- filter(vulcanot_3xtg, log2cociente_3xtg <  LFCN & log10pvalor_3xtg > p_umbral)

up_regulated$Status   <- "up-regulated"
down_regulated$Status <- "down-regulated"

df_desreguladas_3xtg <- rbind(down_regulated, up_regulated)

# --- Save outputs -------------------------------------------------------------
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
write.xlsx(vulcanot_3xtg,        "results/tables/vulcanot_3xtg.xlsx",        sheetName = "todas")
write.xlsx(df_desreguladas_3xtg, "results/tables/df_desreguladas_3xtg.xlsx", sheetName = "Desreguladas")

message(sprintf(
  "Differential expression complete.\n  Up-regulated:   %d\n  Down-regulated: %d",
  nrow(up_regulated), nrow(down_regulated)
))
