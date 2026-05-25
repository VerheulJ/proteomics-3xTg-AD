# =============================================================================
# 01_preprocessing.R
# Proteomics Analysis — 3xTg-AD Mouse Model
#
# Description:
#   Loads raw LFQ abundance data, cleans gene symbols, removes outliers,
#   filters low-coverage proteins, imputes missing values, and normalizes
#   abundances by a sample correction factor.
#
# Input:
#   - data/raw/230414-LFQ_Filtradas.xlsx   (Proteome Discoverer output)
#   - data/raw/sample-metadata.xlsx        (sample annotations)
#
# Output:
#   - data/processed/ajustado.xlsx         (normalized abundance matrix)
#   - data/processed/entera_3xtg.xlsx      (filtered full protein table)
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(readxl)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(factoextra)
library(FactoMineR)
library(missMDA)
library(VIM)

options(scipen = 999)

# --- Load data ----------------------------------------------------------------
entera_3xtg    <- read_excel("data/raw/230414-LFQ_Filtradas.xlsx")
sample_metadata <- read_xlsx("data/raw/sample-metadata.xlsx")

# --- Fix missing Gene Symbols -------------------------------------------------
# Manual dictionary for proteins without an annotated gene symbol
gene_symbols <- c(
  "P47964" = "Rpl36",
  "P97461" = "Rps5",
  "Q8C3W1" = "P_Uncharacterized_protein_C1orf198_homolog",
  "Q91V76" = "P_Ester_hydrolase_C11orf54_homolog",
  "Q9CPR4" = "Rpl17",
  "P03987" = "Ig_gamma_3_chain_C_region"
)

entera_3xtg <- entera_3xtg %>%
  mutate(
    `Gene Symbol` = ifelse(
      is.na(`Gene Symbol`) & Accession %in% names(gene_symbols),
      gene_symbols[Accession],
      `Gene Symbol`
    )
  )

# Clean species tag from descriptions
entera_3xtg$Description <- gsub("\\[OS=Mus musculus\\]", "", entera_3xtg$Description)

# Sort by Accession for consistency
entera_3xtg <- entera_3xtg %>% arrange(Accession)

# --- Extract abundance columns ------------------------------------------------
df_experimento_3xtg <- entera_3xtg

accesion_col    <- dplyr::select(df_experimento_3xtg, Accession)
gene_symbol_col <- dplyr::select(df_experimento_3xtg, `Gene Symbol`)
indice_proteina <- cbind(accesion_col, gene_symbol_col)

# Select columns containing "Abundance: "
colnam  <- colnames(entera_3xtg)
indices <- grep("Abundance: ", colnam)
abundancias_experimento_3xtg <- df_experimento_3xtg[, indices]

# Standardize column names to "mouseXX" format
colnam <- gsub(".*mouse(\\d+).*", "\\1", colnames(abundancias_experimento_3xtg))
colnam <- paste("mouse", gsub("\\D+", "", colnam), sep = "")
colnames(abundancias_experimento_3xtg) <- colnam

id <- indice_proteina[, 1]
df_experimento_3xtg <- cbind(accesion_col, gene_symbol_col, abundancias_experimento_3xtg)
rownames(df_experimento_3xtg) <- id$Accession

# Abundance-only matrix (no Accession / Gene Symbol columns)
df_experimento_3xtg_rowname <- df_experimento_3xtg[, -(1:2)]

# --- QC: Missing values per sample --------------------------------------------
missing_by_sample <- colSums(is.na(df_experimento_3xtg_rowname))
df_missing_sample <- data.frame(
  Muestra = names(missing_by_sample),
  NAs     = missing_by_sample
)

ggplot(df_missing_sample, aes(x = Muestra, y = NAs)) +
  geom_bar(stat = "identity", fill = "red", alpha = 0.7) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    title = "Missing values per sample",
    x     = "Sample",
    y     = "Number of NAs"
  )

# --- QC: Preliminary PCA (complete cases only) --------------------------------
tabla_sin_na   <- df_experimento_3xtg_rowname[complete.cases(df_experimento_3xtg_rowname), ]
pca_preliminar <- prcomp(tabla_sin_na, scale. = TRUE, center = TRUE)

pca_df_pre <- as.data.frame(pca_preliminar$x)
ggplot(pca_df_pre, aes(x = PC1, y = PC2)) +
  geom_point() +
  ggtitle("Preliminary PCA (no imputation, complete cases)")

# --- Separate samples by group ------------------------------------------------
prot_control_hembras <- sample_metadata %>%
  filter(experiment == "WT",   sex == "female") %>% pull(`Sample-id`)

prot_control_machos  <- sample_metadata %>%
  filter(experiment == "WT",   sex == "male")   %>% pull(`Sample-id`)

prot_3xtg_machos     <- sample_metadata %>%
  filter(experiment == "3xTg", sex == "male")   %>% pull(`Sample-id`)

prot_3xtg_hembras    <- sample_metadata %>%
  filter(experiment == "3xTg", sex == "female") %>% pull(`Sample-id`)

tabla_control_hembras <- df_experimento_3xtg_rowname %>% select(all_of(prot_control_hembras))
tabla_control_machos  <- df_experimento_3xtg_rowname %>% select(all_of(prot_control_machos))
tabla_3xtg_machos     <- df_experimento_3xtg_rowname %>% select(all_of(prot_3xtg_machos))
tabla_3xtg_hembras    <- df_experimento_3xtg_rowname %>% select(all_of(prot_3xtg_hembras))

# --- Outlier removal (3×IQR per protein per group) ----------------------------
handle_outliers <- function(x, factor = 3) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  x[x < (Q1 - factor * IQR_val) | x > (Q3 + factor * IQR_val)] <- NA
  return(x)
}

tabla_control_hembras <- as.data.frame(t(apply(tabla_control_hembras, 1, handle_outliers)))
tabla_control_machos  <- as.data.frame(t(apply(tabla_control_machos,  1, handle_outliers)))
tabla_3xtg_machos     <- as.data.frame(t(apply(tabla_3xtg_machos,     1, handle_outliers)))
tabla_3xtg_hembras    <- as.data.frame(t(apply(tabla_3xtg_hembras,    1, handle_outliers)))

# --- Filter: require ≥3 valid values per group --------------------------------
tabla_control_hembras <- tabla_control_hembras[rowSums(!is.na(tabla_control_hembras)) >= 3, ]
tabla_control_machos  <- tabla_control_machos[ rowSums(!is.na(tabla_control_machos))  >= 3, ]
tabla_3xtg_machos     <- tabla_3xtg_machos[    rowSums(!is.na(tabla_3xtg_machos))     >= 3, ]
tabla_3xtg_hembras    <- tabla_3xtg_hembras[   rowSums(!is.na(tabla_3xtg_hembras))    >= 3, ]

# --- Keep only proteins detected in all four groups ---------------------------
proteinas_comunes <- Reduce(
  intersect,
  list(
    rownames(tabla_control_hembras),
    rownames(tabla_control_machos),
    rownames(tabla_3xtg_machos),
    rownames(tabla_3xtg_hembras)
  )
)

tabla_control_hembras <- tabla_control_hembras[proteinas_comunes, ]
tabla_control_machos  <- tabla_control_machos[ proteinas_comunes, ]
tabla_3xtg_machos     <- tabla_3xtg_machos[    proteinas_comunes, ]
tabla_3xtg_hembras    <- tabla_3xtg_hembras[   proteinas_comunes, ]

# --- Imputation (within-group median) -----------------------------------------
impute_median <- function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  return(x)
}

tabla_control_hembras <- as.data.frame(t(apply(tabla_control_hembras, 1, impute_median)))
tabla_control_machos  <- as.data.frame(t(apply(tabla_control_machos,  1, impute_median)))
tabla_3xtg_machos     <- as.data.frame(t(apply(tabla_3xtg_machos,     1, impute_median)))
tabla_3xtg_hembras    <- as.data.frame(t(apply(tabla_3xtg_hembras,    1, impute_median)))

# Restore column names
colnames(tabla_control_hembras) <- prot_control_hembras
colnames(tabla_control_machos)  <- prot_control_machos
colnames(tabla_3xtg_machos)     <- prot_3xtg_machos
colnames(tabla_3xtg_hembras)    <- prot_3xtg_hembras

# Combine all groups
imputados <- cbind(tabla_control_hembras, tabla_control_machos,
                   tabla_3xtg_machos, tabla_3xtg_hembras)

# --- Normalization by correction factor ---------------------------------------
# Factor = column mean / minimum column mean  →  divides out loading differences
media                 <- colSums(imputados, na.rm = TRUE) / nrow(imputados)
valor_medio_minimo    <- min(media)
factor_correlaccion   <- media / valor_medio_minimo
factor_correlaccion_rep <- as.data.frame(
  do.call(cbind, lapply(factor_correlaccion, rep, times = nrow(imputados)))
)

ajustado <- imputados / factor_correlaccion_rep
rownames(ajustado) <- rownames(imputados)

# --- Save outputs -------------------------------------------------------------
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write.xlsx(ajustado, "data/processed/ajustado.xlsx", rowNames = TRUE)

# Filter full table to retained proteins
entera_3xtg_filtrada <- entera_3xtg %>% filter(Accession %in% rownames(ajustado))
write.xlsx(entera_3xtg_filtrada, "data/processed/entera_3xtg.xlsx", rowNames = TRUE)

message("Preprocessing complete. Files saved to data/processed/")
