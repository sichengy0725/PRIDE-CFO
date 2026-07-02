setwd('/rsrch8/home/biostatistics/syang10/PRIDE/results/PRIDE')

ntrial <- 4000L
ndose  <- 5L
nsel   <- ndose + 1L   # 6th = STOP
target <- 0.33

sce_list <- list(
  sce1 = c(0.33, 0.45, 0.58, 0.70, 0.80),
  sce2 = c(0.18, 0.33, 0.52, 0.60, 0.70),
  sce3 = c(0.11, 0.22, 0.33, 0.40, 0.50),
  sce4 = c(0.01, 0.02, 0.03, 0.33, 0.50),
  sce5 = c(0.00, 0.00, 0.05, 0.10, 0.33),
  sce6 = c(0.45, 0.55, 0.65, 0.75, 0.85)
)

summarize_one_scenario <- function(res_dir,
                                   ntrial,
                                   truth,
                                   target,
                                   ndose = length(truth),
                                   stop_code = 99L) {
  true_mtd <- which.min(abs(truth - target))
  sel <- integer(ndose + 1L)  # 1..ndose = dose selected, (ndose+1)=STOP
  pat <- integer(ndose)       # patient counts by highest dose received
  sel_correct <- 0L
  sel_overdose <- 0L
  pat_mtd <- 0L
  pat_overdose <- 0L
  pat_total <- 0L
  
  for (i in seq_len(ntrial)) {
    if(!file.exists(file.path(res_dir, paste0("trial-", i)))){
      cat(paste0('trial-', i), '\n')
      next
    }
    obj <- readRDS(file.path(res_dir, paste0("trial-", i)))
    
    mtd <- obj$final$MTD
    if (!is.na(mtd)) {
      if (mtd == stop_code) {
        sel[ndose + 1L] <- sel[ndose + 1L] + 1L
      } else if (mtd >= 1L && mtd <= ndose) {
        sel[mtd] <- sel[mtd] + 1L
      }
      if (mtd == true_mtd) sel_correct <- sel_correct + 1L
      if (mtd >= 1L && mtd <= ndose && mtd > true_mtd) {
        sel_overdose <- sel_overdose + 1L
      }
    }
    
    admin <- obj$admin
    if (!is.null(admin) && nrow(admin) > 0L) {
      admin2 <- admin[!is.na(admin$id) & !is.na(admin$dose), c("id", "dose"), drop = FALSE]
      
      if (nrow(admin2) > 0L) {
        max_dose_by_pt <- tapply(admin2$dose, admin2$id, max)
        max_dose_by_pt <- as.integer(max_dose_by_pt)
        max_dose_by_pt <- max_dose_by_pt[max_dose_by_pt >= 1L & max_dose_by_pt <= ndose]
        
        if (length(max_dose_by_pt) > 0L) {
          pat <- pat + tabulate(max_dose_by_pt, nbins = ndose)
          pat_mtd <- pat_mtd + sum(max_dose_by_pt == true_mtd)
          pat_overdose <- pat_overdose + sum(max_dose_by_pt > true_mtd)
          pat_total <- pat_total + length(max_dose_by_pt)
        }
      }
    }
  }
  
  tot_assign <- sum(pat)
  
  list(
    true_mtd   = true_mtd,
    sel_count  = sel,
    sel_rate   = sel / ntrial,
    mtd_selection = sel_correct / ntrial,
    overdose_selection = sel_overdose / ntrial,
    pat_count  = pat / ntrial,
    pat_rate   = if (tot_assign > 0) pat / tot_assign else rep(NA_real_, ndose),
    mtd_allocation = if (pat_total > 0) pat_mtd / pat_total else NA_real_,
    overdose_allocation = if (pat_total > 0) pat_overdose / pat_total else NA_real_,
    tot_assign = tot_assign / ntrial
  )
}

# Build output table: one row per scenario
nscen <- length(sce_list)
out <- data.frame(scenario = seq_len(nscen), stringsAsFactors = FALSE)

# truth
for (d in 1:ndose) out[[paste0("true_p", d)]] <- NA_real_
out$true_mtd <- NA_integer_

# selection rates
for (d in 1:ndose) out[[paste0("sel_rate_", d)]] <- NA_real_
out$sel_rate_stop <- NA_real_
out$mtd_selection <- NA_real_
out$overdose_selection <- NA_real_

# allocation rates & counts
for (d in 1:ndose) out[[paste0("pat_rate_", d)]] <- NA_real_
for (d in 1:ndose) out[[paste0("pat_n_", d)]]    <- NA_real_
out$mtd_allocation <- NA_real_
out$overdose_allocation <- NA_real_

out$tot_assign <- NA_real_
CFO <- 0
for (s in seq_len(nscen)) {
  truth <- sce_list[[paste0("sce", s)]]
  for (d in 1:ndose) out[s, paste0("true_p", d)] <- truth[d]
  dose_cap = 12
  Nmax = 24
  summ <- summarize_one_scenario(
    res_dir = paste0("res", s, "-CFO-", CFO, '-dose_cap-', dose_cap, 
                     '-Nmax-', Nmax),
    ntrial = ntrial,
    truth = truth,
    target = target,
    ndose = ndose,
    stop_code = 99L
  )
  
  out$true_mtd[s] <- summ$true_mtd
  for (d in 1:ndose) out[s, paste0("sel_rate_", d)] <- summ$sel_rate[d]
  out$sel_rate_stop[s] <- summ$sel_rate[ndose + 1L]
  out$mtd_selection[s] <- summ$mtd_selection
  out$overdose_selection[s] <- summ$overdose_selection
  
  for (d in 1:ndose) out[s, paste0("pat_rate_", d)] <- summ$pat_rate[d]
  for (d in 1:ndose) out[s, paste0("pat_n_", d)]    <- summ$pat_count[d]
  out$mtd_allocation[s] <- summ$mtd_allocation
  out$overdose_allocation[s] <- summ$overdose_allocation
  
  out$tot_assign[s] <- summ$tot_assign
}
filename <- paste0("PRIDE_summary_all_scenarios-CFO-", CFO, '-dose_cap-', dose_cap, 
                   '-Nmax-', Nmax, '.csv')
write.csv(out, file = filename, row.names = FALSE)
