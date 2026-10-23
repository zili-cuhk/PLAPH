# Folder holding luad_pilot_clinical.csv, Z_RNA_var_top2000_os.csv and X_log2CNA_var_top1000_os.csv
root <- Sys.getenv("PLAPH_LUAD_ROOT", "plaph_luad_noimage_pilot_screening")
rank_dir <- file.path(root, "pilot_results")
out_dir <- file.path(root, "working_sets")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

z <- read.csv(file.path(root, "Z_RNA_var_top2000_os.csv"), check.names = FALSE)
x <- read.csv(file.path(root, "X_log2CNA_var_top1000_os.csv"), check.names = FALSE)
clinical <- read.csv(file.path(root, "luad_pilot_clinical.csv"), check.names = FALSE)
rna_rank <- read.csv(file.path(rank_dir, "rna_adjusted_combined_rank_top2000.csv"), check.names = FALSE)
cna_rank <- read.csv(file.path(rank_dir, "cna_adjusted_combined_rank_top1000.csv"), check.names = FALSE)

write_working_set <- function(p, q, label) {
  z_features <- head(rna_rank$feature, p)
  x_features <- head(cna_rank$feature, q)
  z_out <- z[, c("patient_id", z_features), drop = FALSE]
  x_out <- x[, c("patient_id", x_features), drop = FALSE]

  write.csv(clinical, file.path(out_dir, paste0(label, "_clinical.csv")), row.names = FALSE)
  write.csv(z_out, file.path(out_dir, paste0(label, "_Z_RNA_p", p, ".csv")), row.names = FALSE)
  write.csv(x_out, file.path(out_dir, paste0(label, "_X_CNA_q", q, ".csv")), row.names = FALSE)

  z_tab <- rna_rank[match(z_features, rna_rank$feature), ]
  x_tab <- cna_rank[match(x_features, cna_rank$feature), ]
  write.csv(z_tab, file.path(out_dir, paste0(label, "_Z_RNA_p", p, "_features.csv")), row.names = FALSE)
  write.csv(x_tab, file.path(out_dir, paste0(label, "_X_CNA_q", q, "_features.csv")), row.names = FALSE)

  data.frame(
    label = label,
    n = nrow(clinical),
    p = p,
    q = q,
    z_file = paste0(label, "_Z_RNA_p", p, ".csv"),
    x_file = paste0(label, "_X_CNA_q", q, ".csv")
  )
}

summary <- rbind(
  write_working_set(500, 200, "primary_p500_q200"),
  write_working_set(500, 100, "fallback_p500_q100"),
  write_working_set(1000, 100, "sensitivity_p1000_q100")
)
write.csv(summary, file.path(out_dir, "working_sets_summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
