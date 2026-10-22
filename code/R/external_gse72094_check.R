# Optional cross-check of external_gse72094.py with survival::coxph (Efron ties).
# Usage: Rscript external_gse72094_check.R gse72094_cyp17a1_clinical.csv
library(survival)
D <- read.csv(commandArgs(TRUE)[1], check.names = FALSE)
D$time <- suppressWarnings(as.numeric(D$survival_time_in_days)) / 30.4375
D$event <- ifelse(D$vital_status == "Dead", 1, ifelse(D$vital_status == "Alive", 0, NA))
x <- D[["probe_merck-NM_000102_at"]]
D$z <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
D$score <- -0.395 * D$z
E <- subset(D, is.finite(time) & !is.na(event) & is.finite(score) & time > 0)
f <- coxph(Surv(time, event) ~ score, data = E, ties = "efron")
print(summary(f)$conf.int); print(summary(f)$coefficients)
print(concordance(f, reverse = FALSE))
st <- toupper(E$Stage)
E$age <- suppressWarnings(as.numeric(E$age_at_diagnosis)); E$male <- as.numeric(E$gender == "M")
E$s2 <- as.numeric(substr(st, 1, 1) == "2"); E$s3 <- as.numeric(substr(st, 1, 1) == "3"); E$s4 <- as.numeric(substr(st, 1, 1) == "4")
E$s2[!grepl("^[1-4]", st)] <- NA
print(summary(coxph(Surv(time, event) ~ score + age + male + s2 + s3 + s4, data = E, ties = "efron"))$conf.int[1, ])
