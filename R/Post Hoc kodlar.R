##############################
## REQUIRED LIBRARIES
##############################
library(PMCMRplus)
library(stats)
##############################
## REPRODUCIBILITY
##############################
set.seed(12345)
##############################
## PATHS
##############################
checkpoint_dir <- "checkpoints"
dir.create(checkpoint_dir, showWarnings = FALSE)
final_rds <- "results_final.rds"
final_csv <- "results_final.csv"
##############################
## SIMULATION PARAMETERS
##############################
num_replications <- 10000
checkpoint_interval <- 1000
alpha <- 0.05
k_values <- c(3, 4, 5)
n_values <- c(6, 8, 10)
var_conditions <- list(
  Homogeneous = list(k3 = c(1,1,1), k4 = c(1,1,1,1), k5 = c(1,1,1,1,1)),
  Moderate    = list(k3 = c(1,2,3), k4 = c(1,2,3,4), k5 = c(1,2,3,4,5)),
  High        = list(k3 = c(1,4,7), k4 = c(1,4,7,10), k5 = c(1,4,7,10,13))
)
cohens_f <- list(
  k3n6 = 0.4, k3n8 = 0.4, k3n10 = 0.4,
  k4n6 = 0.4, k4n8 = 0.4, k4n10 = 0.4,
  k5n6 = 0.4, k5n8 = 0.4, k5n10 = 0.4
)
##############################
## MEAN GENERATION
##############################
generate_means <- function(k, f, variances) {
  step <- f * sqrt(mean(variances) * k)
  m <- seq(-(k - 1)/2, (k - 1)/2, length.out = k)
  m <- m * step
  m - mean(m)
}
##############################
## SAFE DETECTION FUNCTION
##############################
detect <- function(pvals) {
  if (length(pvals) == 0 || all(is.na(pvals))) return(FALSE)
  any(pvals < alpha, na.rm = TRUE)
}
##############################
## SIMULATION FUNCTION (NORMAL ONLY)
##############################
run_simulation <- function(k, n, variances, means, label) {
  cp_file <- file.path(checkpoint_dir, paste0("checkpoint_", label, ".rds"))
  if (file.exists(cp_file)) {
    cp <- readRDS(cp_file)
    detection  <- cp$detection
    valid_runs <- cp$valid_runs
    start_i    <- cp$i + 1
  } else {
    detection <- c(
      TukeyHSD=0, Scheffe=0, FisherLSD=0,
      GamesHowell=0, DunnettT3=0, TamhaneT2=0, Sidak=0
    )
    valid_runs <- 0
    start_i <- 1
  }
  for (i in start_i:num_replications) {
    value <- unlist(lapply(1:k, function(j)
      rnorm(n, mean = means[j], sd = sqrt(variances[j]))))
    group <- factor(rep(seq_len(k), each = n))
    data  <- data.frame(value, group)
    aov_model <- aov(value ~ group, data = data)
    valid_runs <- valid_runs + 1
    ## Tukey
    detection["TukeyHSD"] <- detection["TukeyHSD"] +
      detect(TukeyHSD(aov_model)$group[, "p adj"])
    ## Scheffe
    detection["Scheffe"] <- detection["Scheffe"] +
      detect(PMCMRplus::scheffeTest(value ~ group, data)$p.value)
    ## Fisher LSD (protected)
    if (summary(aov_model)[[1]][["Pr(>F)"]][1] < alpha) {
      lsd <- pairwise.t.test(value, group, p.adj = "none")$p.value
      detection["FisherLSD"] <- detection["FisherLSD"] + detect(lsd)
    }
    ## Games-Howell
    detection["GamesHowell"] <- detection["GamesHowell"] +
      detect(PMCMRplus::gamesHowellTest(value ~ group, data)$p.value)
    ## Dunnett T3
    detection["DunnettT3"] <- detection["DunnettT3"] +
      detect(PMCMRplus::dunnettT3Test(value ~ group, data)$p.value)
    ## Tamhane T2
    detection["TamhaneT2"] <- detection["TamhaneT2"] +
      detect(PMCMRplus::tamhaneT2Test(value ~ group, data)$p.value)
    ## Sidak
    raw_p <- na.omit(as.vector(
      pairwise.t.test(value, group, p.adj = "none")$p.value
    ))
    sidak_p <- 1 - (1 - raw_p)^length(raw_p)
    detection["Sidak"] <- detection["Sidak"] + detect(sidak_p)
    ## CHECKPOINT
    if (valid_runs %% checkpoint_interval == 0) {
      saveRDS(
        list(detection = detection,
             valid_runs = valid_runs,
             i = i),
        cp_file
      )
      message("Checkpoint saved: ", label,
              " | valid_runs = ", valid_runs)
    }
  }
  detection / valid_runs
}
##############################
## MAIN LOOP
##############################
if (file.exists(final_rds)) {
  results <- readRDS(final_rds)
  message("Final results loaded.")
} else {
  results <- data.frame()
  for (k in k_values) {
    for (n in n_values) {
      for (vc in names(var_conditions)) {
        label <- paste0("k",k,"_n",n,"_",vc)
        message("Running: ", label)
        variances <- var_conditions[[vc]][[paste0("k",k)]]
        f_val     <- cohens_f[[paste0("k",k,"n",n)]]
        power <- run_simulation(
          k, n, variances,
          generate_means(k, f_val, variances),
          paste0(label, "_power")
        )
        type1 <- run_simulation(
          k, n, variances,
          rep(0, k),
          paste0(label, "_type1")
        )
        results <- rbind(
          results,
          data.frame(
            k = k,
            n = n,
            Variance = vc,
            Test = names(power),
            Power = as.numeric(power),
            TypeI_Error = as.numeric(type1)
          )
        )
      }
    }
  }
  saveRDS(results, final_rds)
  write.csv(results, final_csv, row.names = FALSE)
} 
