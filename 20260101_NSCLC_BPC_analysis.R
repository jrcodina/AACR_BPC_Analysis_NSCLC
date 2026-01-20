library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(tibble)
library(survival)
library(survminer)
library(ggplot2)

# Paths
clinical_sample = "./data/data_clinical_sample.txt"
mutations       = "./data/data_mutations_extended.txt"
clin_index      = "./data/cancer_level_dataset_index.csv"
pathology       = "./data/pathology_report_level_dataset.csv"
regimens        = "./data/regimen_cancer_level_dataset.csv"
out_dir         = "./plots"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 1) Genomics (MAF -> wide flags)
patient_sample <- read.table(
  clinical_sample,
  sep = "\t", header = TRUE,
  quote = "", comment.char = "", fill = TRUE,
  stringsAsFactors = FALSE,
  skip = 4
) %>%
  select(PATIENT_ID, SAMPLE_ID)

maf_df <- read.table(
  mutations,
  sep = "\t", header = TRUE,
  quote = "", comment.char = "", fill = TRUE,
  stringsAsFactors = FALSE
) %>%
  left_join(patient_sample, by = c("Tumor_Sample_Barcode" = "SAMPLE_ID")) %>%
  mutate(record_id = PATIENT_ID) %>%
  select(-PATIENT_ID)

genes <- c("EGFR","KRAS","TP53","STK11","KEAP1","MET","BRAF","PIK3CA","ALK","ROS1","RET","ERBB2")

m <- maf_df %>%
  filter(Hugo_Symbol %in% genes, !is.na(record_id)) %>%
  distinct(record_id, Hugo_Symbol) %>%
  mutate(altered = 1L) %>%
  pivot_wider(
    names_from = Hugo_Symbol,
    values_from = altered,
    values_fill = 0L
  )

mut4 <- m %>%
  select(any_of(c("record_id", "EGFR", "KRAS", "STK11", "KEAP1"))) %>%
  mutate(across(any_of(c("EGFR", "KRAS", "STK11", "KEAP1")), ~ replace_na(as.integer(.), 0L)))

# 2) Clinical
clin <- read.csv(clin_index, na.strings = c("", "NA", "Unknown")) %>%
  filter(cohort == "NSCLC") %>%
  select(record_id, ca_seq, stage_dx, stage_dx_iv, ca_lung_cigarette)

# 3) PD-L1 (max tumor-cell percent per patient)
path <- read.csv(pathology, na.strings = c("", "NA", "Unknown")) %>%
  filter(cohort == "NSCLC")

pdl1 <- path %>%
  mutate(
    across(c(pdl1_perc, pdl1_perc_2, pdl1_perc_3), ~ parse_number(as.character(.)))
  ) %>%
  mutate(
    pdl1_row = pmax(pdl1_perc, pdl1_perc_2, pdl1_perc_3, na.rm = TRUE),
    pdl1_row = ifelse(is.infinite(pdl1_row), NA_real_, pdl1_row)
  ) %>%
  group_by(record_id) %>%
  summarise(
    pdl1_perc = if (all(is.na(pdl1_row))) NA_real_ else max(pdl1_row, na.rm = TRUE),
    .groups = "drop"
  )

# 4) Regimens (IO only; first IO per patient)
io_regex <- "pembrolizumab|nivolumab|atezolizumab|durvalumab|cemiplimab"

reg <- read.csv(regimens, na.strings = c("", "NA", "Unknown")) %>%
  filter(cohort == "NSCLC") %>%
  mutate(is_io = str_detect(tolower(regimen_drugs), io_regex)) %>%
  filter(is_io)

reg1 <- reg %>%
  group_by(record_id) %>%
  slice_min(order_by = drugs_startdt_int_1, with_ties = FALSE) %>%
  ungroup()

stopifnot(all(c("ttnt_ca_seq_yrs","ttnt_ca_seq_status","tt_os_g_yrs","os_g_status") %in% names(reg1)))

# 5) Analysis table
dat_all <- clin %>%
  left_join(pdl1, by = "record_id") %>%
  left_join(mut4, by = "record_id") %>%
  left_join(
    reg1 %>% select(record_id, ttnt_ca_seq_yrs, ttnt_ca_seq_status, tt_os_g_yrs, os_g_status),
    by = "record_id"
  ) %>%
  mutate(
    EGFR  = replace_na(EGFR, 0L),
    KRAS  = replace_na(KRAS, 0L),
    STK11 = replace_na(STK11, 0L),
    KEAP1 = replace_na(KEAP1, 0L),
    pdl1_grp = case_when(
      is.na(pdl1_perc) ~ NA_character_,
      pdl1_perc >= 50 ~ "High (>=50%)",
      TRUE ~ "Low (<50%)"
    ),
    pdl1_grp = factor(pdl1_grp, levels = c("Low (<50%)", "High (>=50%)")),
    smoking_ever = case_when(
      ca_lung_cigarette == "Never used" ~ "Never",
      is.na(ca_lung_cigarette) ~ NA_character_,
      TRUE ~ "Ever"
    ),
    smoking_ever = factor(smoking_ever, levels = c("Never","Ever")),
    stageIV = ifelse(stage_dx_iv == "Stage IV", "Stage IV", "Not stage IV"),
    stageIV = factor(stageIV, levels = c("Not stage IV","Stage IV"))
  )

# 6) PD-L1 survival analyses (TTNT and OS) + save combined plot + risk table
dat_ttnt <- dat_all %>%
  filter(!is.na(pdl1_grp), !is.na(ttnt_ca_seq_yrs), !is.na(ttnt_ca_seq_status))

fit_ttnt <- survfit(Surv(ttnt_ca_seq_yrs, ttnt_ca_seq_status) ~ pdl1_grp, data = dat_ttnt)

p_ttnt_km <- ggsurvplot(
  fit_ttnt,
  data = dat_ttnt,
  risk.table = TRUE,
  pval = TRUE,
  legend.title = "PD-L1 (TC %)",
  xlab = "Years since first IO regimen start",
  ylab = "Probability of no next treatment or death (TTNT)"
)

png(file.path(out_dir, "KM_TTNT_by_PDL1.png"), width = 15, height = 8, units = "in", res = 300)
print(p_ttnt_km)
dev.off()

cox_ttnt_unadj <- coxph(Surv(ttnt_ca_seq_yrs, ttnt_ca_seq_status) ~ pdl1_grp, data = dat_ttnt)
cox_ttnt_adj   <- coxph(Surv(ttnt_ca_seq_yrs, ttnt_ca_seq_status) ~ pdl1_grp + smoking_ever + stageIV, data = dat_ttnt)

p_ttnt_forest <- ggforest(
  cox_ttnt_adj,
  data = dat_ttnt,
  main = "TTNT Cox model (adjusted)",
  cpositions = c(0.02, 0.30, 0.45),
  noDigits = 2
)

ggsave(
  filename = file.path(out_dir, "Forest_TTNT_Cox_Adjusted.png"),
  plot = p_ttnt_forest,
  width = 15, height = 6, dpi = 300
)

dat_os <- dat_all %>%
  filter(!is.na(pdl1_grp), !is.na(tt_os_g_yrs), !is.na(os_g_status))

fit_os <- survfit(Surv(tt_os_g_yrs, os_g_status) ~ pdl1_grp, data = dat_os)

p_os_km <- ggsurvplot(
  fit_os,
  data = dat_os,
  risk.table = TRUE,
  pval = TRUE,
  legend.title = "PD-L1 (TC %)",
  xlab = "Years since first IO regimen start",
  ylab = "Overall survival probability"
)

png(file.path(out_dir, "KM_OS_by_PDL1.png"), width = 15, height = 8, units = "in", res = 300)
print(p_os_km)
dev.off()

cox_os_unadj <- coxph(Surv(tt_os_g_yrs, os_g_status) ~ pdl1_grp, data = dat_os)
cox_os_adj   <- coxph(Surv(tt_os_g_yrs, os_g_status) ~ pdl1_grp + smoking_ever + stageIV, data = dat_os)

p_os_forest <- ggforest(
  cox_os_adj,
  data = dat_os,
  main = "OS Cox model (adjusted)",
  cpositions = c(0.02, 0.30, 0.45),
  noDigits = 2
)

ggsave(
  filename = file.path(out_dir, "Forest_OS_Cox_Adjusted.png"),
  plot = p_os_forest,
  width = 15, height = 6, dpi = 300
)

# 7) Driver groups and KRAS co-mutations
dat_g <- dat_all %>%
  mutate(
    driver_grp = case_when(
      EGFR == 1L ~ "EGFR",
      KRAS == 1L ~ "KRAS",
      TRUE ~ "Other"
    ),
    driver_grp = factor(driver_grp, levels = c("EGFR", "KRAS", "Other")),
    kras_comut = case_when(
      KRAS == 1L & (STK11 == 1L | KEAP1 == 1L) ~ "KRAS + (STK11 or KEAP1)",
      KRAS == 1L ~ "KRAS only",
      TRUE ~ NA_character_
    ),
    kras_comut = factor(kras_comut, levels = c("KRAS only", "KRAS + (STK11 or KEAP1)"))
  )

# 8) KM by driver group + save combined plot + risk table
dat_os_drv <- dat_g %>%
  filter(!is.na(driver_grp), !is.na(tt_os_g_yrs), !is.na(os_g_status))

fit_os_drv <- survfit(Surv(tt_os_g_yrs, os_g_status) ~ driver_grp, data = dat_os_drv)

p_os_drv <- ggsurvplot(
  fit_os_drv,
  data = dat_os_drv,
  risk.table = TRUE,
  pval = TRUE,
  xlab = "Years since first IO regimen start",
  ylab = "Overall survival probability",
  legend.title = "Driver group"
)

png(file.path(out_dir, "KM_OS_by_DriverGroup.png"), width = 15, height = 8, units = "in", res = 300)
print(p_os_drv)
dev.off()

dat_ttnt_drv <- dat_g %>%
  filter(!is.na(driver_grp), !is.na(ttnt_ca_seq_yrs), !is.na(ttnt_ca_seq_status))

fit_ttnt_drv <- survfit(Surv(ttnt_ca_seq_yrs, ttnt_ca_seq_status) ~ driver_grp, data = dat_ttnt_drv)

p_ttnt_drv <- ggsurvplot(
  fit_ttnt_drv,
  data = dat_ttnt_drv,
  risk.table = TRUE,
  pval = TRUE,
  xlab = "Years since first IO regimen start",
  ylab = "Probability of no next treatment or death (TTNT)",
  legend.title = "Driver group"
)

png(file.path(out_dir, "KM_TTNT_by_DriverGroup.png"), width = 15, height = 8, units = "in", res = 300)
print(p_ttnt_drv)
dev.off()

# 9) KM within KRAS: KRAS only vs KRAS + (STK11 or KEAP1) + save combined plot + risk table
dat_os_kras <- dat_g %>%
  filter(KRAS == 1L, !is.na(kras_comut), !is.na(tt_os_g_yrs), !is.na(os_g_status))

fit_os_kras <- survfit(Surv(tt_os_g_yrs, os_g_status) ~ kras_comut, data = dat_os_kras)

p_os_kras <- ggsurvplot(
  fit_os_kras,
  data = dat_os_kras,
  risk.table = TRUE,
  pval = TRUE,
  xlab = "Years since first IO regimen start",
  ylab = "Overall survival probability",
  legend.title = "KRAS subgroup"
)

png(file.path(out_dir, "KM_OS_KRAS_Comut.png"), width = 15, height = 8, units = "in", res = 300)
print(p_os_kras)
dev.off()

dat_ttnt_kras <- dat_g %>%
  filter(KRAS == 1L, !is.na(kras_comut), !is.na(ttnt_ca_seq_yrs), !is.na(ttnt_ca_seq_status))

fit_ttnt_kras <- survfit(Surv(ttnt_ca_seq_yrs, ttnt_ca_seq_status) ~ kras_comut, data = dat_ttnt_kras)

p_ttnt_kras <- ggsurvplot(
  fit_ttnt_kras,
  data = dat_ttnt_kras,
  risk.table = TRUE,
  pval = TRUE,
  xlab = "Years since first IO regimen start",
  ylab = "Probability of no next treatment or death (TTNT)",
  legend.title = "KRAS subgroup"
)

png(file.path(out_dir, "KM_TTNT_KRAS_Comut.png"), width = 15, height = 8, units = "in", res = 300)
print(p_ttnt_kras)
dev.off()

# 10) PD-L1 low (<50%) within KRAS + save combined plot + risk table
dat_os_kras_pdl1low <- dat_g %>%
  filter(
    KRAS == 1L,
    !is.na(kras_comut),
    !is.na(pdl1_perc), pdl1_perc < 50,
    !is.na(tt_os_g_yrs), !is.na(os_g_status)
  )

dat_ttnt_kras_pdl1low <- dat_g %>%
  filter(
    KRAS == 1L,
    !is.na(kras_comut),
    !is.na(pdl1_perc), pdl1_perc < 50,
    !is.na(ttnt_ca_seq_yrs), !is.na(ttnt_ca_seq_status)
  )

fit_os_low <- survfit(Surv(tt_os_g_yrs, os_g_status) ~ kras_comut, data = dat_os_kras_pdl1low)

p_os_low <- ggsurvplot(
  fit_os_low,
  data = dat_os_kras_pdl1low,
  risk.table = TRUE,
  pval = TRUE,
  xlab = "Years since first IO regimen start",
  ylab = "Overall survival probability",
  legend.title = "KRAS subgroup",
  title = "PD-L1 Low (<50%): OS by KRAS co-mutation status"
)

png(file.path(out_dir, "KM_OS_KRAS_Comut_PDL1Low.png"), width = 15, height = 8, units = "in", res = 300)
print(p_os_low)
dev.off()

fit_ttnt_low <- survfit(Surv(ttnt_ca_seq_yrs, ttnt_ca_seq_status) ~ kras_comut, data = dat_ttnt_kras_pdl1low)

p_ttnt_low <- ggsurvplot(
  fit_ttnt_low,
  data = dat_ttnt_kras_pdl1low,
  risk.table = TRUE,
  pval = TRUE,
  xlab = "Years since first IO regimen start",
  ylab = "Probability of no next treatment or death (TTNT)",
  legend.title = "KRAS subgroup",
  title = "PD-L1 Low (<50%): TTNT by KRAS co-mutation status"
)

png(file.path(out_dir, "KM_TTNT_KRAS_Comut_PDL1Low.png"), width = 15, height = 8, units = "in", res = 300)
print(p_ttnt_low)
dev.off()

# 11) PD-L1 high (>50%) within KRAS + save combined plot + risk table
dat_os_kras_pdl1high <- dat_g %>%
  filter(
    KRAS == 1L,
    !is.na(kras_comut),
    !is.na(pdl1_perc), pdl1_perc >= 50,
    !is.na(tt_os_g_yrs), !is.na(os_g_status)
  )

#dat_os_kras_pdl1high$kras_comut <- ifelse(is.na(dat_os_kras_pdl1high$kras_comut), "WT", dat_os_kras_pdl1high$kras_comut)

dat_ttnt_kras_pdl1high <- dat_g %>%
  filter(
    KRAS == 1L,
    !is.na(kras_comut),
    !is.na(pdl1_perc), pdl1_perc >= 50,
    !is.na(ttnt_ca_seq_yrs), !is.na(ttnt_ca_seq_status)
  )

fit_os_high <- survfit(Surv(tt_os_g_yrs, os_g_status) ~ kras_comut, data = dat_os_kras_pdl1high)

p_os_high <- ggsurvplot(
  fit_os_high,
  data = dat_os_kras_pdl1high,
  risk.table = TRUE,
  pval = TRUE,
  xlab = "Years since first IO regimen start",
  ylab = "Overall survival probability",
  legend.title = "KRAS subgroup",
  title = "PD-L1 high (>50%): OS by KRAS co-mutation status"
)

png(file.path(out_dir, "KM_OS_KRAS_Comut_PDL1High.png"), width = 15, height = 8, units = "in", res = 300)
print(p_os_high)
dev.off()

fit_ttnt_high <- survfit(Surv(ttnt_ca_seq_yrs, ttnt_ca_seq_status) ~ kras_comut, data = dat_ttnt_kras_pdl1high)

p_ttnt_high <- ggsurvplot(
  fit_ttnt_high,
  data = dat_ttnt_kras_pdl1high,
  risk.table = TRUE,
  pval = TRUE,
  xlab = "Years since first IO regimen start",
  ylab = "Probability of no next treatment or death (TTNT)",
  legend.title = "KRAS subgroup",
  title = "PD-L1 high (>50%): TTNT by KRAS co-mutation status"
)

png(file.path(out_dir, "KM_TTNT_KRAS_Comut_PDL1high.png"), width = 15, height = 8, units = "in", res = 300)
print(p_ttnt_high)
dev.off()


cox_os_low_adj <- coxph(Surv(tt_os_g_yrs, os_g_status) ~ kras_comut + smoking_ever + stageIV, data = dat_os_kras_pdl1low)
cox_ttnt_low_adj <- coxph(Surv(ttnt_ca_seq_yrs, ttnt_ca_seq_status) ~ kras_comut + smoking_ever + stageIV, data = dat_ttnt_kras_pdl1low)

summary(cox_ttnt_unadj)
summary(cox_ttnt_adj)
summary(cox_os_unadj)
summary(cox_os_adj)
summary(cox_os_low_adj)
summary(cox_ttnt_low_adj)

# 12) Denominator summaries
ids_all <- clin %>% distinct(record_id)

ids_maf <- maf_df %>%
  filter(!is.na(record_id)) %>%
  distinct(record_id)

ids_pdl1 <- pdl1 %>%
  mutate(pdl1_perc = as.numeric(pdl1_perc)) %>%
  filter(!is.na(pdl1_perc)) %>%
  distinct(record_id)

ids_io_any <- reg %>% distinct(record_id)

ids_io_os <- reg1 %>%
  filter(!is.na(tt_os_g_yrs), !is.na(os_g_status)) %>%
  distinct(record_id)

ids_io_ttnt <- reg1 %>%
  filter(!is.na(ttnt_ca_seq_yrs), !is.na(ttnt_ca_seq_status)) %>%
  distinct(record_id)

ids_for_os_km <- ids_all %>%
  inner_join(ids_pdl1, by = "record_id") %>%
  inner_join(ids_io_os, by = "record_id")

ids_for_ttnt_km <- ids_all %>%
  inner_join(ids_pdl1, by = "record_id") %>%
  inner_join(ids_io_ttnt, by = "record_id")

denom_tbl <- tibble(
  Metric = c(
    "Total NSCLC",
    "With MAF",
    "With PD-L1",
    "With both PD-L1 + MAF",
    "Received IO (any regimen)",
    "IO + OS fields available (from IO start)",
    "IO + TTNT fields available (from IO start)",
    "PD-L1 + IO + OS (OS KM denominator)",
    "PD-L1 + IO + TTNT (TTNT KM denominator)"
  ),
  N = c(
    n_distinct(ids_all$record_id),
    n_distinct(ids_maf$record_id),
    n_distinct(ids_pdl1$record_id),
    n_distinct(inner_join(ids_pdl1, ids_maf, by = "record_id")$record_id),
    n_distinct(ids_io_any$record_id),
    n_distinct(ids_io_os$record_id),
    n_distinct(ids_io_ttnt$record_id),
    n_distinct(ids_for_os_km$record_id),
    n_distinct(ids_for_ttnt_km$record_id)
  )
)

mut_tbl <- tibble(
  Metric = c(
    "N EGFR mutated",
    "N KRAS only mutated",
    "N KRAS + (STK11 or KEAP1) co-mutated"
  ),
  N = c(
    sum(dat_g$EGFR == 1L, na.rm = TRUE),
    sum(dat_g$KRAS == 1L & dat_g$STK11 == 0L & dat_g$KEAP1 == 0L, na.rm = TRUE),
    sum(dat_g$KRAS == 1L & (dat_g$STK11 == 1L | dat_g$KEAP1 == 1L), na.rm = TRUE)
  )
)

denom_tbl
mut_tbl
