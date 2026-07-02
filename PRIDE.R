library(coda)
source('CFO_tox_utils.R')
select.mtd <- function(target, y, n, cutoff.eli=0.95)
{
  ## isotonic transformation using the pool adjacent violator algorithm (PAVA)
  pava <- function (x, wt = rep(1, length(x))) 
  {
    n <- length(x)
    if (n <= 1) 
      return(x)
    if (any(is.na(x)) || any(is.na(wt))) {
      stop("Missing values in 'x' or 'wt' not allowed")
    }
    lvlsets <- (1:n)
    repeat {
      viol <- (as.vector(diff(x)) < 0)
      if (!(any(viol))) 
        break
      i <- min((1:(n - 1))[viol])
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i + 1]
      ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
      x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
      lvlsets[ilvl] <- lvl1
    }
    x
  }
  ## determine whether the dose has been eliminated during the trial
  ndose=length(n);
  elimi=rep(0, ndose);
  for(i in 1:ndose)
  {
    if(n[i]>2) {if(1-pbeta(target, y[i]+1, n[i]-y[i]+1)>cutoff.eli) {elimi[i:ndose]=1; break;}}
  }
  
  if(elimi[1]==1) { selectdose=99; } ## no dose should be selected if the first dose is already very toxic
  else
  {
    nadmis = min(max(which(elimi==0)), max(which(n!=0))); ## the highest admissble (or un-eliminated) dose level
    ## poster mean and variance of toxicity probabilities using beta(0.005, 0.005) as the prior 
    phat = (y[1:nadmis]+0.005)/(n[1:nadmis]+0.01); 
    phat.var = (y[1:nadmis]+0.005)*(n[1:nadmis]-y[1:nadmis]+0.005)/((n[1:nadmis]+0.01)^2*(n[1:nadmis]+0.01+1))
    
    ## perform the isotonic transformation using PAVA
    phat = pava(phat, wt=1/phat.var) 
    phat = phat + (1:nadmis)*1E-10 ## break ties by adding an increasingly small number 
    selectdose = sort(abs(phat-target), index.return=T)$ix[1]  ## select dose closest to the target as the MTD
  }
  return(selectdose);  
}
make_pride_jags_data <- function(tmp,
                                 K,
                                 mu,              # length K: logit(c_k)
                                 sigma2_beta = 10,
                                 eta = 1) {
  stopifnot(is.data.frame(tmp))
  if (!all(c("id", "dose") %in% names(tmp))) stop("tmp must contain id and dose.")
  if (!("y_obs" %in% names(tmp)) && !("y" %in% names(tmp))) {
    stop("tmp must contain y_obs or y.")
  }
  
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


get_pride_posterior <- function(tmp,
                                K,
                                mu,
                                TARGET = 0.30,
                                model_file = "PRIDE.bug",
                                sigma2_beta = 10,
                                eta = 1,
                                m_use = Inf,
                                pk_method = c("approx", "mc"),
                                n_mc_w = 200,
                                seed_subsample = NULL,
                                n.chains = 3,
                                n.adapt  = 1000,
                                n.burn   = 2000,
                                n.iter   = 5000,
                                thin     = 2,
                                monitor  = c("beta", "sigma_w")) {
  
  pk_method <- match.arg(pk_method)
  
  pk_from_beta <- function(beta_draws, sigma_w,
                           method = c("approx", "mc"),
                           n_mc_w = 200) {
    method <- match.arg(method)
    beta_draws <- as.matrix(beta_draws)   # M x K
    sigma_w <- as.numeric(sigma_w)        # length M
    M <- nrow(beta_draws)
    K <- ncol(beta_draws)
    
    if (method == "approx") {
      # logistic-normal moment approximation
      denom <- sqrt(1 + (pi^2 / 3) * (sigma_w^2))
      return(plogis(beta_draws / denom))
    }
    
    pk <- matrix(NA_real_, nrow = M, ncol = K)
    for (m in seq_len(M)) {
      Wm <- rnorm(n_mc_w, 0, sigma_w[m])
      mat <- sweep(
        matrix(beta_draws[m, ], nrow = n_mc_w, ncol = K, byrow = TRUE),
        1, Wm, "+"
      )
      pk[m, ] <- colMeans(plogis(mat))
    }
    pk
  }
  
  # prior-only case
  if (is.null(tmp) || nrow(tmp) == 0L) {
    M <- n.iter
    beta_draws <- matrix(
      rnorm(M * K, mean = rep(mu, each = M), sd = sqrt(sigma2_beta)),
      nrow = M, ncol = K, byrow = FALSE
    )
    
    sigma2_w <- 1 / rgamma(M, shape = eta, rate = eta)
    sigma_w  <- sqrt(sigma2_w)
    
    if (is.finite(m_use) && m_use < M) {
      if (!is.null(seed_subsample)) set.seed(seed_subsample)
      idx <- sample.int(M, size = m_use)
      beta_draws <- beta_draws[idx, , drop = FALSE]
      sigma_w <- sigma_w[idx]
    }
    
    pk_draws <- pk_from_beta(beta_draws, sigma_w, method = pk_method, n_mc_w = n_mc_w)
    posttox <- colMeans(pk_draws)
    prob_overtox <- colMeans(pk_draws > TARGET)
    
    return(list(
      pk_draws = pk_draws,
      posttox = posttox,
      beta_draws = beta_draws,
      sigma_w_draws = sigma_w,
      prob_overtox = prob_overtox
    ))
  }
  
  data_jags <- make_pride_jags_data(
    tmp = tmp, K = K, mu = mu,
    sigma2_beta = sigma2_beta, eta = eta
  )
  
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
  
  beta_cols <- paste0("beta[", seq_len(K), "]")
  if (!all(beta_cols %in% colnames(draws))) {
    missing <- beta_cols[!beta_cols %in% colnames(draws)]
    stop("Missing beta columns: ", paste(missing, collapse = ", "))
  }
  if (!("sigma_w" %in% colnames(draws))) {
    stop("sigma_w not found in posterior draws.")
  }
  
  beta_draws <- draws[, beta_cols, drop = FALSE]
  sigma_w <- draws[, "sigma_w"]
  
  M <- nrow(beta_draws)
  if (is.finite(m_use) && m_use < M) {
    if (!is.null(seed_subsample)) set.seed(seed_subsample)
    idx <- sample.int(M, size = m_use)
    beta_draws <- beta_draws[idx, , drop = FALSE]
    sigma_w <- sigma_w[idx]
  }
  
  pk_draws <- pk_from_beta(beta_draws, sigma_w, method = pk_method, n_mc_w = n_mc_w)
  posttox <- colMeans(pk_draws)
  prob_overtox <- colMeans(pk_draws > TARGET)
  
  list(
    pk_draws = pk_draws,
    posttox = posttox,
    beta_draws = beta_draws,
    sigma_w_draws = sigma_w,
    prob_overtox = prob_overtox
  )
}


cfo_move_pride <- function(cur_dose,
                           pk_draws,
                           TARGET,
                           gammaL,
                           gammaR,
                           elim = NULL,
                           use_monotone_pair = TRUE,
                           eps = 1e-12,
                           CFO = 0) {
  if (!is.matrix(pk_draws)) pk_draws <- as.matrix(pk_draws)
  K <- ncol(pk_draws)
  if (K < 2) stop("pk_draws must have at least 2 dose columns.")
  if (cur_dose < 1 || cur_dose > K) stop("cur_dose out of range.")
  if (is.null(elim)) elim <- rep(0L, K)
  
  L <- cur_dose - 1L
  R <- cur_dose + 1L
  
  can_left  <- (L >= 1L)
  can_right <- (R <= K)
  
  post_odds_over <- function(p_draw) {
    pr_over <- mean(p_draw > TARGET)
    pr_over <- min(1 - eps, max(eps, pr_over))
    pr_over / (1 - pr_over)
  }
  
  left_stat <- NA_real_
  if (can_left) {
    pL <- pk_draws[, L]
    pC <- pk_draws[, cur_dose]
    keep <- if (use_monotone_pair) (pL < pC) else rep(TRUE, length(pL))
    if (any(keep)) {
      OL <- post_odds_over(pL[keep])
      OC_L <- post_odds_over(pC[keep])
      left_stat <- OC_L * OL
    }
  }
  right_stat <- NA_real_
  if (can_right) {
    pC <- pk_draws[, cur_dose]
    pR <- pk_draws[, R]
    keep <- if (use_monotone_pair) (pC < pR) else rep(TRUE, length(pC))
    if (any(keep)) {
      OC_R <- post_odds_over(pC[keep])
      OR <- post_odds_over(pR[keep])
      right_stat <- (1/OC_R)/ OR
    }
  }
  # browser()
  # If current dose has been eliminated, force de-escalation
  if (elim[cur_dose] == 1L) {
    if (cur_dose == 1L) return(0L)
    lower_ok <- which(elim[seq_len(cur_dose - 1L)] == 0L)
    if (length(lower_ok) == 0L) return(0L)
    return(max(lower_ok))
  }
  
  go_left  <- can_left  && is.finite(left_stat)  && (left_stat  > gammaL)
  go_right <- can_right && is.finite(right_stat) && (right_stat > gammaR) && (elim[R] == 0L)
  
  if (go_left && !go_right) {
    return(max(cur_dose - 1L, 1L))
  } else if (!go_left && go_right) {
    return(min(cur_dose + 1L, K))
  } else {
    return(cur_dose)
  }
}
cfo_move <- function(d, ndose, target, y, n, elimi, gammaL, gammaR, eps = 1e-12) {
  if (d < 1 || d > ndose) stop("d out of range.")
  L <- d - 1L
  R <- d + 1L
  
  can_left  <- (L >= 1L)
  can_right <- (R <= ndose)
  
  alp.prior <- target
  bet.prior <- 1 - target
  
  OL <- if (can_left) {
    OR.values(target, y[L], n[L], y[d], n[d], alp.prior, bet.prior, type = "L")
  } else {
    NA_real_
  }
  
  OR <- if (can_right) {
    OR.values(target, y[d], n[d], y[R], n[R], alp.prior, bet.prior, type = "R")
  } else {
    NA_real_
  }
  
  go_left  <- can_left  && (OL > gammaL)
  go_right <- can_right && (OR > gammaR)
  
  if (go_left && !go_right) {
    return(max(d - 1L, 1L))
  } else if (!go_left && go_right) {
    if (elimi[R] == 0) return(R)
    return(d)
  } else {
    return(d)
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
  
  # --- MCMC controls ---
  n.chains = 3,
  n.adapt  = 1000,
  n.burn   = 2000,
  n.iter   = 5000,
  thin     = 2,
  
  # --- misc ---
  seed = NULL,
  verbose = FALSE,
  CFO = TRUE
) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(length(p_true) == K, length(mu) == K)
  
  # ---------------- helpers ----------------
  get_patient_state <- function(admin) {
    if (nrow(admin) == 0) return(admin)
    a <- admin[order(admin$id, admin$t_eval, admin$row_id), , drop = FALSE]
    a[!duplicated(a$id, fromLast = TRUE), , drop = FALSE]
  }
  
  eligible_ipde_ids <- function(admin, next_dose, cycle_max) {
    st <- get_patient_state(admin)
    if (nrow(st) == 0) return(integer(0))
    ok <- (st$y == 0L) & (st$ncycle < cycle_max) & (st$dose < next_dose)
    st <- st[ok, , drop = FALSE]
    st <- st[order(st$t_eval, st$id), , drop = FALSE]
    st$id
  }
  
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
  
  elimi <- rep(0L, K)
  
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
  stop_trial <- FALSE
  
  if (verbose) {
    message("Cohort 1: dose=", cur_dose, " decision_time=", t_decision)
  }
  
  # ---------------- main loop ----------------
  repeat {
    if (stop_trial) break
    if (nrow(admin) >= Nmax_eff) break
    while (next_new_idx <= N_pat && patients$t_arrival[next_new_idx] < t_decision) {
      waiting <- rbind(waiting, patients[next_new_idx, c("id", "t_arrival"), drop = FALSE])
      next_new_idx <- next_new_idx + 1L
    }
    
    dat_dec <- admin[admin$t_eval <= t_decision, c("id", "dose", "y"), drop = FALSE]
    
    n_left  <- if (cur_dose > 1L) sum(dat_dec$dose == (cur_dose - 1L)) else 0L
    n_curr  <- sum(dat_dec$dose == cur_dose)
    n_right <- if (cur_dose < K) sum(dat_dec$dose == (cur_dose + 1L)) else 0L
    n_by_dose <- tabulate(dat_dec$dose, nbins = K)
    y_by_dose <- tabulate(dat_dec$dose[dat_dec$y == 1], nbins = K)
    # model-based elimination rule
    if(CFO == TRUE){
      
      ## determine if the current dose should be eliminated
      post_overtox <- 1 - pbeta(TARGET, 
                                TARGET + y_by_dose[cur_dose], 
                                1 - TARGET + n_by_dose[cur_dose] - y_by_dose[cur_dose])
      if (n_by_dose[cur_dose] >= 3L && post_overtox > cutoff) {
        elimi[cur_dose:K] <- 1L
        if (cur_dose == 1L) {
          decisions <- rbind(decisions, data.frame(
            cohort = cohort_id + 1L,
            t_decision = t_decision,
            cur_dose = cur_dose,
            next_dose = cur_dose,
            stop = 1L
          ))
          stop_trial <- TRUE
          if (verbose) {
            message("STOP: dose 1 eliminated under CFO rule. Pr(overtox)=",
                    round(post_overtox, 4))
          }
          break
        }
      }
      
    } else {
      post_dec <- get_pride_posterior(
        tmp = dat_dec,
        K = K,
        mu = mu,
        TARGET = TARGET,
        model_file = model_file,
        sigma2_beta = sigma2_beta,
        eta = eta,
        pk_method = pk_method,
        n_mc_w = n_mc_w,
        m_use = m_use,
        n.chains = n.chains,
        n.adapt = n.adapt,
        n.burn = n.burn,
        n.iter = n.iter,
        thin = thin
      )

      post_overtox_curr <- post_dec$prob_overtox[cur_dose]
      if (post_overtox_curr > cutoff) {
        elimi[cur_dose:K] <- 1L
        if (cur_dose == 1L) {
          decisions <- rbind(decisions, data.frame(
            cohort = cohort_id + 1L,
            t_decision = t_decision,
            cur_dose = cur_dose,
            next_dose = cur_dose,
            stop = 1L
          ))
          stop_trial <- TRUE
          if (verbose) {
            message("STOP: dose 1 eliminated. Pr(overtox)=", round(post_overtox_curr, 4))
          }
          break
        }
      }
    }
   
    alp.prior <- TARGET
    bet.prior <- 1 - TARGET
    
    if (cur_dose == 1L) {
      gammaL <- Inf
      gammaR <- optim.gamma.fn(
        n1 = n_curr, n2 = n_right,
        phi = TARGET, type = "R",
        alp.prior = alp.prior, bet.prior = bet.prior
      )$gamma
    } else if (cur_dose == K) {
      gammaR <- Inf
      gammaL <- optim.gamma.fn(
        n1 = n_left, n2 = n_curr,
        phi = TARGET, type = "L",
        alp.prior = alp.prior, bet.prior = bet.prior
      )$gamma
    } else {
      gammaL <- optim.gamma.fn(
        n1 = n_left, n2 = n_curr,
        phi = TARGET, type = "L",
        alp.prior = alp.prior, bet.prior = bet.prior
      )$gamma
      
      gammaR <- optim.gamma.fn(
        n1 = n_curr, n2 = n_right,
        phi = TARGET, type = "R",
        alp.prior = alp.prior, bet.prior = bet.prior
      )$gamma
    }
    if(CFO == TRUE){
      next_dose <- cfo_move(
        d = cur_dose,
        ndose = K,
        target = TARGET,
        y = y_by_dose , n = n_by_dose,
        elimi = elimi,
        gammaL = gammaL,
        gammaR = gammaR
      )
    
    } else {
      next_dose <- cfo_move_pride(
        cur_dose = cur_dose,
        pk_draws = post_dec$pk_draws,
        TARGET = TARGET,
        gammaL = gammaL,
        gammaR = gammaR,
        elim = elimi,
        use_monotone_pair = FALSE
      )
    }
    
    
    if (next_dose == 0L) {
      decisions <- rbind(decisions, data.frame(
        cohort = cohort_id + 1L,
        t_decision = t_decision,
        cur_dose = cur_dose,
        next_dose = cur_dose,
        stop = 1L
      ))
      stop_trial <- TRUE
      if (verbose) message("STOP: no admissible lower dose remains.")
      break
    }
    
    decisions <- rbind(decisions, data.frame(
      cohort = cohort_id + 1L,
      t_decision = t_decision,
      cur_dose = cur_dose,
      next_dose = next_dose,
      stop = 0L
    ))
    
    # ---------------- form new cohort ----------------
    cohort_id <- cohort_id + 1L
    t_start <- t_decision
    
    # A) pull from waiting queue first
    new_ids <- integer(0)
    if (nrow(waiting) > 0) {
      take_n <- min(C, nrow(waiting))
      new_ids <- waiting$id[1:take_n]
      waiting <- waiting[-seq_len(take_n), , drop = FALSE]
    }
    
    # B) fill with retreat patients
    need <- C - length(new_ids)
    ret_ids <- integer(0)
    if (need > 0) {
      cand <- eligible_ipde_ids(admin, next_dose = next_dose, cycle_max = cycle_max)
      if (length(cand) > 0) ret_ids <- head(cand, need)
    }
    
    # C) if still not full, wait for new arrivals
    need2 <- C - length(new_ids) - length(ret_ids)
    if (need2 > 0) {
      while (need2 > 0 && next_new_idx <= N_pat) {
        pid <- patients$id[next_new_idx]
        tA  <- patients$t_arrival[next_new_idx]
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
    
    while (next_new_idx <= N_pat && patients$t_arrival[next_new_idx] < t_start) {
      waiting <- rbind(waiting, patients[next_new_idx, c("id", "t_arrival"), drop = FALSE])
      next_new_idx <- next_new_idx + 1L
    }
    
    # ---------------- append cohort outcomes ----------------
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
    
    cohort_rows <- admin$cohort == cohort_id
    if (!any(cohort_rows)) {
      if (verbose) message("No observations added for cohort ", cohort_id, "; stopping trial.")
      break
    }
    t_decision <- max(admin$t_eval[cohort_rows])
    cur_dose <- next_dose
    
    if (verbose) {
      message("Cohort ", cohort_id,
              ": dose=", next_dose,
              " start=", t_start,
              " decision_time=", t_decision,
              " n_eff=", nrow(admin),
              " waiting=", nrow(waiting))
    }
    
    if (next_new_idx > N_pat && nrow(waiting) == 0) break
  }
  
  # ---------------- final MTD selection ----------------
  if (nrow(admin) == 0L) {
    return(list(
      admin = admin,
      waiting = waiting,
      decisions = decisions,
      final = list(
        t_end = NA_real_,
        post = NULL,
        posttox_iso = rep(NA_real_, K),
        MTD = 99L,
        eliminated = elimi
      )
    ))
  }
  
  t_end <- max(admin$t_eval, na.rm = TRUE)
  dat_final <- admin[admin$t_eval <= t_end, c("id", "dose", "y"), drop = FALSE]
  
  if(CFO == TRUE){
    post_final <- NULL
    posttox_iso <- rep(NA_real_, K)
    elim_final <- elimi
    
    n_by_dose <- tabulate(dat_final$dose, nbins = K)
    y_by_dose <- tabulate(dat_final$dose[dat_final$y == 1], nbins = K)
    final_dose = select.mtd(TARGET, y_by_dose, n_by_dose, cutoff)
    
    if(final_dose == 99){elim_final <- rep(1,K)}
  } else {
    post_final <- get_pride_posterior(
      tmp = dat_final,
      K = K,
      mu = mu,
      TARGET = TARGET,
      model_file = model_file,
      sigma2_beta = sigma2_beta,
      eta = eta,
      pk_method = "mc",
      n_mc_w = 50,
      m_use = 1000,
      n.chains = n.chains,
      n.adapt = n.adapt,
      n.burn = n.burn,
      n.iter = n.iter,
      thin = thin
    )
    
    posttox_iso <- Iso::pava(post_final$posttox)
    
    # apply final admissibility using the same model-based elimination idea
    n_by_dose <- tabulate(dat_final$dose, nbins = K)
    elim_final <- elimi
    for (d in seq_len(K)) {
      if (post_final$prob_overtox[d] > cutoff) {
        elim_final[d:K] <- 1L
        break
      }
    }
    
    dvec <- abs(posttox_iso - TARGET)
    admissible <- which((elim_final == 0L) & (n_by_dose > 0L))
    
    if (length(admissible) == 0L) {
      final_dose <- 99L
    } else {
      dvec2 <- rep(Inf, K)
      dvec2[admissible] <- dvec[admissible]
      final_dose <- which(dvec2 == min(dvec2))[1]   # tie-break: lower dose
      final_dose <- as.integer(final_dose)
    }
  }
  
  
  
  list(
    admin = admin,
    waiting = waiting,
    decisions = decisions,
    final = list(
      t_end = t_end,
      post = post_final,
      posttox_iso = posttox_iso,
      MTD = final_dose,
      eliminated = elim_final
    )
  )
}
