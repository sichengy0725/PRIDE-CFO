library(coda)
library(rjags)
library(mvtnorm)

make_pride_jags_data <- function(tmp,
                                 K,
                                 mu,              # length K: logit(c_k)
                                 sigma2_beta = 10,
                                 eta = 1) {
  
  stopifnot(is.data.frame(tmp))
  if (!all(c("id", "dose") %in% names(tmp))) stop("tmp must contain id and dose.")
  if (!("y_obs" %in% names(tmp)) && !("y" %in% names(tmp))) stop("tmp must contain y_obs or y.")
  
  y <- if ("y_obs" %in% names(tmp)) tmp$y_obs else tmp$y
  if (any(!y %in% c(0, 1))) stop("y must be 0/1.")
  
  uid <- sort(unique(tmp$id))
  id_map <- setNames(seq_along(uid), uid)
  id <- as.integer(id_map[as.character(tmp$id)])
  
  dose <- as.integer(tmp$dose)
  if (any(dose < 1 | dose > K)) stop("dose must be in 1..K (dose level index).")
  
  list(
    Nobs = nrow(tmp),
    Nid  = length(uid),
    K    = as.integer(K),
    y    = as.integer(y),
    id   = id,
    dose = dose,
    mu   = as.numeric(mu),
    sigma2_beta = as.numeric(sigma2_beta),
    eta  = as.numeric(eta)
  )
}

estimate_PRIDE_JAGS <- function(tmp,
                                K,
                                mu,
                                TARGET = 0.30,
                                cutoff = 0.95,
                                model_file = "PRIDE.bug",
                                sigma2_beta = 10,
                                eta = 1,
                                # --- speed controls ---
                                m_use = Inf,                 # e.g., 1000 to subsample posterior draws
                                pk_method = c("approx","mc"),
                                n_mc_w = 200,                # used only if pk_method="mc"
                                seed_subsample = NULL,
                                # --- MCMC controls ---
                                n.chains = 3,
                                n.adapt  = 1000,
                                n.burn   = 2000,
                                n.iter   = 5000,
                                thin     = 2,
                                monitor = c("beta", "sigma_w")) {
  
  pk_method <- match.arg(pk_method)
  
  # -------- helper: compute pk_draws from beta_draws + sigma_w
  pk_from_beta <- function(beta_draws, sigma_w,
                           method = c("approx","mc"),
                           n_mc_w = 200) {
    method <- match.arg(method)
    beta_draws <- as.matrix(beta_draws)     # M x K
    sigma_w <- as.numeric(sigma_w)          # length M
    M <- nrow(beta_draws); K <- ncol(beta_draws)
    
    if (method == "approx") {
      # E[logit^{-1}(beta + W)] approx logit^{-1}( beta / sqrt(1 + (pi^2/3)*sigma^2) )
      denom <- sqrt(1 + (pi^2 / 3) * (sigma_w^2))   # length M
      return(plogis(beta_draws / denom))            # broadcasts denom across columns
    }
    
    # method == "mc": vectorized MC over W, no per-dose loops
    pk <- matrix(NA_real_, nrow = M, ncol = K)
    for (m in 1:M) {
      Wm <- rnorm(n_mc_w, 0, sigma_w[m])           # length n_mc_w
      # n_mc_w x K matrix: each row is beta[m,] + Wm[j]
      # replicate beta row n_mc_w times and add Wm
      mat <- sweep(matrix(beta_draws[m, ], nrow = n_mc_w, ncol = K, byrow = TRUE),
                   1, Wm, "+")
      pk[m, ] <- colMeans(plogis(mat))
    }
    pk
  }
  
  # -------- prior-only case
  if (is.null(tmp) || nrow(tmp) == 0L) {
    M <- n.iter
    beta_draws <- matrix(rnorm(M * K, mean = rep(mu, each = M), sd = sqrt(sigma2_beta)),
                         nrow = M, ncol = K, byrow = FALSE)
    
    sigma2_w <- 1 / rgamma(M, shape = eta, rate = eta)
    sigma_w  <- sqrt(sigma2_w)
    
    # optionally subsample prior draws too
    if (is.finite(m_use) && m_use < M) {
      if (!is.null(seed_subsample)) set.seed(seed_subsample)
      idx <- sample.int(M, size = m_use)
      beta_draws <- beta_draws[idx, , drop = FALSE]
      sigma_w <- sigma_w[idx]
      M <- nrow(beta_draws)
    }
    
    pk_draws <- pk_from_beta(beta_draws, sigma_w, method = pk_method, n_mc_w = n_mc_w)
    
    prob_overtox <- mean(pk_draws[, 1] > TARGET)
    stop_flag <- as.integer(prob_overtox > cutoff)
    
    posttox <- colMeans(pk_draws)
    mtd <- which.min(abs(posttox - TARGET))
    
    return(list(
      MTD = mtd,
      posttox = posttox,
      pk_draws = pk_draws,
      beta_draws = beta_draws,
      sigma_w_draws = sigma_w,
      prob_overtox = prob_overtox,
      stop = stop_flag
    ))
  }
  
  # -------- posterior case
  data_jags <- make_pride_jags_data(tmp, K = K, mu = mu,
                                    sigma2_beta = sigma2_beta, eta = eta)
  
  jags <- rjags::jags.model(
    file     = model_file,
    data     = data_jags,
    n.chains = n.chains,
    n.adapt  = n.adapt,
    quiet    = TRUE
  )
  update(jags, n.burn, progress.bar = "none")
  
  smp <- rjags::coda.samples(
    model          = jags,
    variable.names = monitor,
    n.iter         = n.iter,
    thin           = thin,
    progress.bar   = "none"
  )
  draws <- as.matrix(smp)
  
  beta_cols <- paste0("beta[", 1:K, "]")
  if (!all(beta_cols %in% colnames(draws))) {
    missing <- beta_cols[!beta_cols %in% colnames(draws)]
    stop("Missing beta columns: ", paste(missing, collapse = ", "))
  }
  beta_draws <- draws[, beta_cols, drop = FALSE]
  
  if (!("sigma_w" %in% colnames(draws))) stop("sigma_w not found in posterior draws.")
  sigma_w <- draws[, "sigma_w"]
  
  # subsample posterior draws if requested
  M <- nrow(beta_draws)
  if (is.finite(m_use) && m_use < M) {
    if (!is.null(seed_subsample)) set.seed(seed_subsample)
    idx <- sample.int(M, size = m_use)
    beta_draws <- beta_draws[idx, , drop = FALSE]
    sigma_w <- sigma_w[idx]
  }
  
  pk_draws <- pk_from_beta(beta_draws, sigma_w, method = pk_method, n_mc_w = n_mc_w)
  
  prob_overtox <- mean(pk_draws[, 1] > TARGET)
  stop_flag <- as.integer(prob_overtox > cutoff)
  
  posttox <- colMeans(pk_draws)
  dose.best <- which.min(abs(posttox - TARGET))
  
  list(
    MTD = dose.best,
    posttox = posttox,
    pk_draws = pk_draws,
    beta_draws = beta_draws,
    sigma_w_draws = sigma_w,
    prob_overtox = prob_overtox,
    stop = stop_flag
  )
}
cfo_move <- function(cur_dose,
                     pk_draws,     # matrix: M x K, columns are p_k draws
                     TARGET,
                     gammaL, gammaR,
                     eps = 1e-12) {
  
  # draws_are <- match.arg(draws_are)
  
  if (!is.matrix(pk_draws)) pk_draws <- as.matrix(pk_draws)
  K <- ncol(pk_draws)
  if (K < 2) stop("pk_draws must have at least 2 dose columns.")
  if (cur_dose < 1 || cur_dose > K) stop("cur_dose out of range.")
  
  # Identify neighbors
  L <- cur_dose - 1L
  R <- cur_dose + 1L
  
  # If at boundary, only one direction is possible
  can_left  <- (L >= 1L)
  can_right <- (R <= K)
  
  
  pC <- pk_draws[, cur_dose]
  pL <- if (can_left)  pk_draws[, L] else NULL
  pR <- if (can_right) pk_draws[, R] else NULL
  
  # Helper: posterior odds Ok = Pr(p_k > TARGET) / Pr(p_k <= TARGET)
  post_odds_over <- function(p_draw, target, eps) {
    pr_over <- mean(p_draw > target)
    pr_over <- pmin(1 - eps, pmax(eps, pr_over))
    pr_over / (1 - pr_over)
  }
  
  OC <- post_odds_over(pC, TARGET, eps)
  OL <- if (can_left)  post_odds_over(pL, TARGET, eps) else NA_real_
  OR <- if (can_right) post_odds_over(pR, TARGET, eps) else NA_real_
  
  # Odds ratios used by CFO (paper):
  #   OC/OL supports moving left (de-escalation) when large
  #   OC/OR supports moving right (escalation) when large
  go_left  <- can_left  && ((OC / (1-OL)) > gammaL)
  go_right <- can_right && (((1-OC) / OR) > gammaR)
  # browser()
  # Table 1 logic (competing): if both push or neither push => stay
  if (go_left && !go_right) {
    return(max(cur_dose - 1L,1))
  } else if (!go_left && go_right) {
    return(min(cur_dose + 1L,ncol(pk_draws)))
  } else {
    return(cur_dose)
  }
}


simulate_PRIDE_design <- function(
    # --- trial size / timing ---
  N_pat       = 30L,      # number of unique patients with arrival times generated
  Nmax_eff    = 24L,      # max effective observations (rows in admin)
  C           = 3L,       # cohort size
  T_assess    = 3,        # assessment window length
  cycle_max   = 3L,       # max cycles per patient
  
  # --- arrival process ---
  arrival_rate = 1/2,     # Exponential interarrival rate (mean=2)
  t0 = 0,
  
  # --- truth (DGM) ---
  p_true,                 # length K
  sigma_true_w = 0.0,     # optional within-patient RE in truth
  
  # --- PRIDE model inputs ---
  K,
  mu,                     # length K: logit(c_k)
  TARGET = 0.30,
  cutoff = 0.95,
  model_file = "PRIDE.bug",
  sigma2_beta = 10,
  eta = 1,
  pk_method = "mc",
  n_mc_w = 50,
  m_use = 1000,
  

  # --- misc ---
  seed = NULL,
  verbose = FALSE
) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(length(p_true) == K, length(mu) == K)
  
  # ---------------- helpers ----------------
  # latest cycle record per patient (by t_eval then row_id)
  get_patient_state <- function(admin) {
    if (nrow(admin) == 0) return(admin)
    a <- admin[order(admin$id, admin$t_eval, admin$row_id), , drop = FALSE]
    a[!duplicated(a$id, fromLast = TRUE), , drop = FALSE]
  }
  
  # eligible for IPDE retreat: y==0, ncycle<cycle_max, dose<next_dose
  eligible_ipde_ids <- function(admin, next_dose, cycle_max) {
    st <- get_patient_state(admin)
    if (nrow(st) == 0) return(integer(0))
    ok <- (st$y == 0L) & (st$ncycle < cycle_max) & (st$dose < next_dose)
    st <- st[ok, , drop = FALSE]
    # stable ordering: earlier finished first
    st <- st[order(st$t_eval, st$id), , drop = FALSE]
    st$id
  }
  
  # patient-specific truth random effect (optional)
  w_i <- rnorm(N_pat, 0, sigma_true_w)
  gen_y <- function(pid, dose) {
    lp <- qlogis(p_true[dose]) + w_i[pid]
    rbinom(1, 1, plogis(lp))
  }
  
  # ---------------- arrivals ----------------
  inter <- rexp(N_pat, rate = arrival_rate)
  t_arr <- t0 + cumsum(inter)
  patients <- data.frame(id = 1L:N_pat, t_arrival = t_arr)
  
  next_new_idx <- 1L
  
  # ---------------- storage ----------------
  admin <- data.frame(
    row_id = integer(0),
    id = integer(0),
    t_arrival = numeric(0),
    t_start = numeric(0),
    t_eval = numeric(0),
    dose = integer(0),
    y = integer(0),
    ncycle = integer(0),
    cohort = integer(0),
    type = character(0), # "new" or "retreat"
    stringsAsFactors = FALSE
  )
  
  waiting <- data.frame(
    id = integer(0),
    t_arrival = numeric(0),
    stringsAsFactors = FALSE
  )
  
  decisions <- data.frame(
    cohort = integer(0),
    t_decision = numeric(0),
    cur_dose = integer(0),
    next_dose = integer(0),
    stop = integer(0),
    stringsAsFactors = FALSE
  )
  
  # ---------------- 1) first cohort ----------------
  cur_dose <- 1L
  cohort_id <- 1L
  
  if (N_pat < C) stop("Need N_pat >= C.")
  
  first_ids <- patients$id[1:C]
  t_start1 <- max(t0, max(patients$t_arrival[1:C]))
  
  for (pid in first_ids) {
    yij <- gen_y(pid, cur_dose)
    admin <- rbind(admin, data.frame(
      row_id = nrow(admin) + 1L,
      id = pid,
      t_arrival = patients$t_arrival[pid],
      t_start = t_start1,
      t_eval = t_start1 + T_assess,
      dose = cur_dose,
      y = as.integer(yij),
      ncycle = 1L,
      cohort = cohort_id,
      type = "new",
      stringsAsFactors = FALSE
    ))
  }
  
  next_new_idx <- C + 1L
  t_decision <- max(admin$t_eval[admin$cohort == cohort_id])
  
  if (verbose) message("Cohort 1: dose=", cur_dose, " decision_time=", t_decision)
  
  stop_trial <- FALSE
  
  # ---------------- main loop ----------------
  repeat {
    if (stop_trial) break
    if (nrow(admin) >= Nmax_eff) break
    
    # (2) generate arrivals for the next cohort window:
    # any arrivals before decision time get queued
    while (next_new_idx <= N_pat && patients$t_arrival[next_new_idx] < t_decision) {
      waiting <- rbind(waiting, patients[next_new_idx, c("id","t_arrival"), drop = FALSE])
      next_new_idx <- next_new_idx + 1L
    }
    
    # decision using all data available at decision time
    dat_dec <- admin[admin$t_eval <= t_decision, c("id","dose","y"), drop = FALSE]
    
    post_dec <- estimate_PRIDE_JAGS(
      tmp = dat_dec,
      K = K,
      mu = mu,
      TARGET = TARGET,
      cutoff = cutoff,
      model_file = model_file,
      sigma2_beta = sigma2_beta,
      eta = eta,
      pk_method = pk_method,
      n_mc_w = n_mc_w,
      m_use = m_use
    )
    
    if (isTRUE(post_dec$stop == 1L)) {
      decisions <- rbind(decisions, data.frame(
        cohort = cohort_id + 1L,
        t_decision = t_decision,
        cur_dose = cur_dose,
        next_dose = cur_dose,
        stop = 1L
      ))
      stop_trial <- TRUE
      if (verbose) message("STOP triggered at decision time ", t_decision)
      break
    }
    # browser()
    # (iv) CFO move: decide next cohort dose
      # gammatable <- CFO::gammatable(npatient = nrow(dat_dec), target = TARGET)
     
    n_left  <- sum(dat_dec$dose == (cur_dose - 1L))
    n_curr  <- sum(dat_dec$dose ==  cur_dose)
    n_right <- sum(dat_dec$dose == (cur_dose + 1L))
    alp.prior <- TARGET
    bet.prior <- 1 - TARGET
    if(cur_dose == 1){
      gammaL <- Inf
      gammaR <- optim.gamma.fn(n1 = n_curr, n2 = n_right,
                               phi = TARGET, type = "R",
                               alp.prior = alp.prior, bet.prior = bet.prior)$gamma
    } 
    if(cur_dose == length(p_true)){
      gammaR <- Inf
      gammaL <- optim.gamma.fn(n1 = n_left, n2 = n_curr,
                               phi = TARGET, type = "L",
                               alp.prior = alp.prior, bet.prior = bet.prior)$gamma
    }
    else{
      gammaL <- optim.gamma.fn(n1 = n_left, n2 = n_curr,
                               phi = TARGET, type = "L",
                               alp.prior = alp.prior, bet.prior = bet.prior)$gamma
      
      gammaR <- optim.gamma.fn(n1 = n_curr, n2 = n_right,
                               phi = TARGET, type = "R",
                               alp.prior = alp.prior, bet.prior = bet.prior)$gamma
    }
    
    
      
    next_dose <- cfo_move(
      cur_dose = cur_dose,
      pk_draws = post_dec$pk_draws,
      TARGET = TARGET,
      gammaL = gammaL,
      gammaR = gammaR
    )
    
    decisions <- rbind(decisions, data.frame(
      cohort = cohort_id + 1L,
      t_decision = t_decision,
      cur_dose = cur_dose,
      next_dose = next_dose,
      stop = 0L
    ))
    
    # ---------------- form new cohort ----------------
    cohort_id <- cohort_id + 1L
    
    # cohort start time: at least decision time
    t_start <- t_decision
    
    # A) pull from waiting queue first
    new_ids <- integer(0)
    if (nrow(waiting) > 0) {
      take_n <- min(C, nrow(waiting))
      new_ids <- waiting$id[1:take_n]
      waiting <- waiting[-seq_len(take_n), , drop = FALSE]
    }
    
    # B) if not full, add IPDE retreat patients
    need <- C - length(new_ids)
    ret_ids <- integer(0)
    if (need > 0) {
      cand <- eligible_ipde_ids(admin, next_dose = next_dose, cycle_max = cycle_max)
      if (length(cand) > 0) ret_ids <- head(cand, need)
    }
    
    # C) if still not full, wait for new arrivals after decision time
    need2 <- C - length(new_ids) - length(ret_ids)
    if (need2 > 0) {
      while (need2 > 0 && next_new_idx <= N_pat) {
        pid <- patients$id[next_new_idx]
        tA  <- patients$t_arrival[next_new_idx]
        # if this arrival is after decision time, cohort starts after last needed arrival
        t_start <- max(t_start, tA)
        new_ids <- c(new_ids, pid)
        next_new_idx <- next_new_idx + 1L
        need2 <- need2 - 1L
      }
      if (need2 > 0) {
        if (verbose) message("Ran out of patients: cannot form a full cohort.")
        break
      }
    }
    
    # after choosing cohort members, move any arrivals with t_arrival < t_start into waiting
    while (next_new_idx <= N_pat && patients$t_arrival[next_new_idx] < t_start) {
      waiting <- rbind(waiting, patients[next_new_idx, c("id","t_arrival"), drop = FALSE])
      next_new_idx <- next_new_idx + 1L
    }
    
    # ---------------- generate outcomes & append to admin ----------------
    # new patients (cycle 1)
    for (pid in new_ids) {
      if (nrow(admin) >= Nmax_eff) break
      yij <- gen_y(pid, next_dose)
      admin <- rbind(admin, data.frame(
        row_id = nrow(admin) + 1L,
        id = pid,
        t_arrival = patients$t_arrival[pid],
        t_start = t_start,
        t_eval = t_start + T_assess,
        dose = next_dose,
        y = as.integer(yij),
        ncycle = 1L,
        cohort = cohort_id,
        type = "new",
        stringsAsFactors = FALSE
      ))
    }
    
    # retreated patients (cycle increment)
    if (length(ret_ids) > 0) {
      st <- get_patient_state(admin)
      for (pid in ret_ids) {
        if (nrow(admin) >= Nmax_eff) break
        last <- st[st$id == pid, , drop = FALSE]
        if (nrow(last) != 1) next
        yij <- gen_y(pid, next_dose)
        admin <- rbind(admin, data.frame(
          row_id = nrow(admin) + 1L,
          id = pid,
          t_arrival = last$t_arrival,
          t_start = t_start,
          t_eval = t_start + T_assess,
          dose = next_dose,
          y = as.integer(yij),
          ncycle = as.integer(last$ncycle + 1L),
          cohort = cohort_id,
          type = "retreat",
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # update decision time as eval time of the last patient in this cohort
    t_decision <- max(admin$t_eval[admin$cohort == cohort_id])
    cur_dose <- next_dose
    
    if (verbose) {
      message("Cohort ", cohort_id,
              ": dose=", next_dose,
              " start=", t_start,
              " decision_time=", t_decision,
              " n_eff=", nrow(admin),
              " waiting=", nrow(waiting))
    }
    
    # hard stop if no more patients possible
    if (next_new_idx > N_pat && nrow(waiting) == 0) break
  }
  
  # ---------------- final MTD selection (your exact block) ----------------
  t_end <- max(admin$t_eval, na.rm = TRUE)
  dat_final <- admin[admin$t_eval <= t_end, , drop = FALSE]
  
  post_final <- estimate_PRIDE_JAGS(
    tmp = dat_final,
    K = K,
    mu = mu,
    pk_method = "mc",
    n_mc_w = 50,
    m_use = 1000,
    TARGET = TARGET,
    cutoff = cutoff,
    model_file = model_file,
    sigma2_beta = sigma2_beta,
    eta = eta
  )
  browser()
  posttox_iso <- Iso::pava(post_final$posttox)
  d <- abs(posttox_iso - TARGET)
  final_dose <- which(d == min(d))[1]   # tie-break: lower dose
  final_dose <- as.integer(final_dose)
  
  list(
    admin = admin,
    waiting = waiting,
    decisions = decisions,
    final = list(
      t_end = t_end,
      post = post_final,
      posttox_iso = posttox_iso,
      MTD = final_dose
    )
  )
}