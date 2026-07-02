decide_next_step_PRIDE <- function(
    admin,
    patients = NULL,             # optional data.frame(id, t_arrival) for future/new patients
    K,
    mu,
    TARGET = 0.30,
    cutoff = 0.95,
    C = 3L,
    cycle_max = 3L,
    T_assess = 3,
    model_file = "PRIDE.bug",
    sigma2_beta = 10,
    eta = 1,
    pk_method = "mc",
    n_mc_w = 50,
    m_use = 2000,
    n.chains = 3,
    n.adapt  = 1000,
    n.burn   = 2000,
    n.iter   = 5000,
    thin     = 2,
    CFO = TRUE,
    waiting = NULL,              # optional current waiting queue: data.frame(id, t_arrival)
    verbose = FALSE
) {
  if (is.null(admin) || nrow(admin) == 0L) {
    stop("admin must contain at least one observed administration.")
  }
  
  req_cols <- c("row_id", "id", "t_arrival", "t_start", "t_eval",
                "dose", "y", "ncycle", "cohort", "type")
  miss <- setdiff(req_cols, names(admin))
  if (length(miss) > 0L) {
    stop("admin is missing required columns: ", paste(miss, collapse = ", "))
  }
  
  admin <- admin[order(admin$t_eval, admin$row_id), , drop = FALSE]
  
  if (is.null(waiting)) {
    waiting <- data.frame(
      id = integer(0),
      t_arrival = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  
  get_patient_state <- function(admin) {
    a <- admin[order(admin$id, admin$t_eval, admin$row_id), , drop = FALSE]
    a[!duplicated(a$id, fromLast = TRUE), , drop = FALSE]
  }
  
  eligible_ipde_ids <- function(admin, next_dose, cycle_max, decision_time) {
    st <- get_patient_state(admin)
    if (nrow(st) == 0L) return(integer(0))
    
    ok <- (st$y == 0L) &
      (st$ncycle < cycle_max) &
      (st$dose < next_dose) &
      (st$t_eval <= decision_time)
    
    st <- st[ok, , drop = FALSE]
    if (nrow(st) == 0L) return(integer(0))
    
    # priority:
    # 1) fewer cycles first
    # 2) earlier enrollment first
    # 3) earlier last evaluation first
    # 4) smaller id
    st <- st[order(st$ncycle, st$t_arrival, st$t_eval, st$id), , drop = FALSE]
    st$id
  }
  
  # current decision time: last completed cohort
  last_cohort <- max(admin$cohort, na.rm = TRUE)
  t_decision <- max(admin$t_eval[admin$cohort == last_cohort], na.rm = TRUE)
  cur_dose <- admin$dose[which.max(admin$row_id[admin$cohort == last_cohort])]
  cur_dose <- admin$dose[admin$cohort == last_cohort][1]
  
  # use only completed observations by decision time
  dat_dec <- admin[admin$t_eval <= t_decision, c("id", "dose", "y"), drop = FALSE]
  
  # update waiting queue using patient arrival list if supplied
  if (!is.null(patients)) {
    if (!all(c("id", "t_arrival") %in% names(patients))) {
      stop("patients must contain columns id and t_arrival.")
    }
    
    used_ids <- unique(c(admin$id, waiting$id))
    newly_available <- patients[
      patients$t_arrival <= t_decision & !(patients$id %in% used_ids),
      c("id", "t_arrival"),
      drop = FALSE
    ]
    
    if (nrow(newly_available) > 0L) {
      waiting <- rbind(waiting, newly_available)
      waiting <- waiting[order(waiting$t_arrival, waiting$id), , drop = FALSE]
    }
  }
  
  n_left  <- if (cur_dose > 1L) sum(dat_dec$dose == (cur_dose - 1L)) else 0L
  n_curr  <- sum(dat_dec$dose == cur_dose)
  n_right <- if (cur_dose < K) sum(dat_dec$dose == (cur_dose + 1L)) else 0L
  
  n_by_dose <- tabulate(dat_dec$dose, nbins = K)
  y_by_dose <- tabulate(dat_dec$dose[dat_dec$y == 1], nbins = K)
  
  elimi <- rep(0L, K)
  
  # reconstruct elimination status from all completed data up to current decision
  if (CFO) {
    for (d in seq_len(K)) {
      if (n_by_dose[d] >= 3L) {
        post_overtox_d <- 1 - pbeta(
          TARGET,
          TARGET + y_by_dose[d],
          1 - TARGET + n_by_dose[d] - y_by_dose[d]
        )
        if (post_overtox_d > cutoff) {
          elimi[d:K] <- 1L
          break
        }
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
    
    for (d in seq_len(K)) {
      if (post_dec$prob_overtox[d] > cutoff) {
        elimi[d:K] <- 1L
        break
      }
    }
  }
  
  # if current dose eliminated, handle stop / forced de-escalation
  if (elimi[cur_dose] == 1L && cur_dose == 1L) {
    return(list(
      stop = TRUE,
      reason = "Current dose eliminated and no lower dose exists",
      t_decision = t_decision,
      cur_dose = cur_dose,
      next_dose = 0L,
      next_start_time = NA_real_,
      next_decision_time = NA_real_,
      waiting = waiting,
      new_ids = integer(0),
      retreat_ids = integer(0),
      elimi = elimi,
      n_by_dose = n_by_dose,
      y_by_dose = y_by_dose
    ))
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
  
  if (CFO) {
    
    next_dose <- cfo_move(
      d = cur_dose,
      ndose = K,
      target = TARGET,
      y = y_by_dose,
      n = n_by_dose,
      elimi = elimi,
      gammaL = gammaL,
      gammaR = gammaR
    )
  } else {
    if (!exists("post_dec")) {
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
    }
    
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
    return(list(
      stop = TRUE,
      reason = "No admissible lower dose remains",
      t_decision = t_decision,
      cur_dose = cur_dose,
      next_dose = 0L,
      next_start_time = NA_real_,
      next_decision_time = NA_real_,
      waiting = waiting,
      new_ids = integer(0),
      retreat_ids = integer(0),
      elimi = elimi,
      n_by_dose = n_by_dose,
      y_by_dose = y_by_dose
    ))
  }
  
  # form the next cohort, but do not generate outcomes
  t_start <- t_decision
  
  new_ids <- integer(0)
  if (nrow(waiting) > 0L) {
    take_n <- min(C, nrow(waiting))
    new_ids <- waiting$id[1:take_n]
    waiting_after_new <- waiting[-seq_len(take_n), , drop = FALSE]
  } else {
    waiting_after_new <- waiting
  }
  
  need <- C - length(new_ids)
  retreat_ids <- integer(0)
  if (need > 0L) {
    cand <- eligible_ipde_ids(
      admin = admin,
      next_dose = next_dose,
      cycle_max = cycle_max,
      decision_time = t_decision
    )
    if (length(cand) > 0L) {
      retreat_ids <- head(cand, need)
    }
  }
  
  need2 <- C - length(new_ids) - length(retreat_ids)
  
  # if future patients list supplied, determine whether later new arrivals fill remaining slots
  future_new_ids <- integer(0)
  if (need2 > 0L && !is.null(patients)) {
    used_ids2 <- unique(c(admin$id, waiting$id, new_ids))
    future_pool <- patients[
      !(patients$id %in% used_ids2) & patients$t_arrival >= t_decision,
      c("id", "t_arrival"),
      drop = FALSE
    ]
    future_pool <- future_pool[order(future_pool$t_arrival, future_pool$id), , drop = FALSE]
    
    if (nrow(future_pool) > 0L) {
      take2 <- min(need2, nrow(future_pool))
      future_new_ids <- future_pool$id[seq_len(take2)]
      t_start <- max(t_start, future_pool$t_arrival[take2])
    }
  }
  
  next_cohort_ids <- c(new_ids, retreat_ids, future_new_ids)
  full_cohort <- length(next_cohort_ids) == C
  
  next_decision_time <- if (full_cohort) t_start + T_assess else NA_real_
  
  if (verbose) {
    message(
      "Decision at t=", t_decision,
      "; current dose=", cur_dose,
      "; next dose=", next_dose,
      "; next cohort size=", length(next_cohort_ids)
    )
  }
  
  list(
    stop = FALSE,
    reason = NA_character_,
    t_decision = t_decision,
    cur_dose = cur_dose,
    next_dose = next_dose,
    gammaL = gammaL,
    gammaR = gammaR,
    next_start_time = if (full_cohort) t_start else NA_real_,
    next_decision_time = next_decision_time,
    next_cohort_full = full_cohort,
    next_cohort_ids = next_cohort_ids,
    new_ids = c(new_ids, future_new_ids),
    retreat_ids = retreat_ids,
    waiting_before = waiting,
    waiting_after_new_pull = waiting_after_new,
    elimi = elimi,
    n_by_dose = n_by_dose,
    y_by_dose = y_by_dose
  )
}

patients <- data.frame(
  id = 1:30,
  t_arrival = (0:29) * 2   # fixed 2-week spacing
)
# my_admin <- data.frame(
#   row_id     = c(1, 2, 3),
#   id         = c(1, 2, 3),
#   t_arrival  = c(0, 2, 4),
#   t_start    = c(4, 4, 4),
#   t_eval     = c(7, 7, 7),
#   dose       = c(1, 1, 1),
#   y          = c(0, 0, 0),
#   ncycle     = c(1, 1, 1),
#   cohort     = c(1, 1, 1),
#   type       = c("new", "new", "new"),
#   stringsAsFactors = FALSE
# )
set.seed(1)
my_admin <- data.frame(
  row_id     = c(1, 2, 3, 4, 5, 6, 7, 8 ,9, 10, 11, 12, 13, 14, 15, 16, 17, 18),
  id         = c(1, 2, 3, 4, 1, 2, 5, 6, 3, 7, 1, 2, 8, 9, 4, 10, 5 ,6),
  t_arrival  = c(0, 2, 4, 6, 0, 2, 8, 10, 4, 12, 0, 2, 14, 16, 6, 18, 8, 10),
  t_start    = c(4, 4, 4, 7, 7, 7, 10, 10, 10, 13, 13, 13, 16, 16, 16, 19, 19, 19),
  t_eval     = c(7, 7, 7, 10, 10 ,10, 13, 13, 13, 16, 16, 16, 19, 19, 19, 22, 22, 22),
  dose       = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5),
  y          = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1),
  ncycle     = c(1, 1, 1, 1, 2, 2, 1, 1, 2, 1, 3, 3, 1, 1, 2, 1, 2, 2),
  cohort     = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
  type       = c("new", "new", "new", "new", "retreat", "retreat", "new",
                 "new", "retreat", "new", "retreat", "retreat", "new", "new", "retreat",
                 "new", "retreat", "retreat"),
  stringsAsFactors = FALSE
)
mu= c(0.1,0.2,0.3,0.4,0.5)
step1 <- decide_next_step_PRIDE(
  admin = my_admin,
  patients = patients,
  K = 5,
  mu = mu,
  TARGET = 0.33,
  cutoff = 0.95,
  C = 3,
  T_assess = 3,
  cycle_max = 3,
  CFO = FALSE
)

step1$next_cohort_ids
step1$new_ids
step1$retreat_ids
step1$next_start_time
step1$next_decision_time
step1$cur_dose
step1$next_dose

# res <- replicate(20, {
#   step1 <- decide_next_step_PRIDE(
#     admin = my_admin,
#     patients = patients,
#     K = 5,
#     mu = mu,
#     TARGET = 0.33,
#     cutoff = 0.95,
#     C = 3,
#     T_assess = 3,
#     cycle_max = 3,
#     CFO = FALSE
#   )
#   step1$next_dose
# })
# 
# table(res)

# patients <- data.frame(
#   id = 1:30,
#   t_arrival = (0:29) * 2
# )
# 
# build_admin_fixed_y <- function(
#     n_cohort,
#     patients,
#     K = 5,
#     mu,
#     TARGET = 0.33,
#     cutoff = 0.95,
#     C = 3,
#     T_assess = 3,
#     cycle_max = 3,
#     CFO = FALSE,
#     cohort_dlt = c(0, 0, 0, 1, 0, 1, 1, 0, 1)
# ) {
#   stopifnot(n_cohort >= 1)
#   stopifnot(length(cohort_dlt) >= n_cohort)
#   
#   # first cohort fixed at dose 1
#   y1 <- c(rep(1L, cohort_dlt[1]), rep(0L, C - cohort_dlt[1]))
#   
#   my_admin <- data.frame(
#     row_id     = 1:C,
#     id         = c(1, 2, 3),
#     t_arrival  = c(0, 2, 4),
#     t_start    = c(4, 4, 4),
#     t_eval     = c(7, 7, 7),
#     dose       = rep(1L, C),
#     y          = y1,
#     ncycle     = c(1, 1, 1),
#     cohort     = rep(1L, C),
#     type       = c("new", "new", "new"),
#     stringsAsFactors = FALSE
#   )
#   
#   if (n_cohort == 1) return(my_admin)
#   
#   for (cc in 2:n_cohort) {
#     step <- decide_next_step_PRIDE(
#       admin = my_admin,
#       patients = patients,
#       K = K,
#       mu = mu,
#       TARGET = TARGET,
#       cutoff = cutoff,
#       C = C,
#       T_assess = T_assess,
#       cycle_max = cycle_max,
#       CFO = CFO
#     )
#     
#     if (isTRUE(step$stop)) {
#       message("Stopped before cohort ", cc)
#       break
#     }
#     
#     ids_this <- step$next_cohort_ids
#     if (length(ids_this) != C) {
#       stop("Cohort ", cc, " not full. ids = ", paste(ids_this, collapse = ", "))
#     }
#     
#     type_this <- ifelse(ids_this %in% step$retreat_ids, "retreat", "new")
#     
#     ncycle_this <- integer(C)
#     t_arrival_this <- numeric(C)
#     
#     for (j in seq_along(ids_this)) {
#       pid <- ids_this[j]
#       old_rows <- my_admin[my_admin$id == pid, , drop = FALSE]
#       
#       if (nrow(old_rows) == 0) {
#         ncycle_this[j] <- 1L
#         t_arrival_this[j] <- patients$t_arrival[patients$id == pid]
#       } else {
#         ncycle_this[j] <- max(old_rows$ncycle) + 1L
#         t_arrival_this[j] <- old_rows$t_arrival[1]
#       }
#     }
#     
#     # fixed cohort-level DLT pattern
#     y_this <- c(rep(1L, cohort_dlt[cc]), rep(0L, C - cohort_dlt[cc]))
#     
#     new_block <- data.frame(
#       row_id     = (nrow(my_admin) + 1):(nrow(my_admin) + C),
#       id         = ids_this,
#       t_arrival  = t_arrival_this,
#       t_start    = rep(step$next_start_time, C),
#       t_eval     = rep(step$next_decision_time, C),
#       dose       = rep(step$next_dose, C),
#       y          = y_this,
#       ncycle     = ncycle_this,
#       cohort     = rep(cc, C),
#       type       = type_this,
#       stringsAsFactors = FALSE
#     )
#     
#     my_admin <- rbind(my_admin, new_block)
#   }
#   
#   my_admin
# }
# test_fixed_y_path <- function(
#     max_cohort = 9,
#     patients,
#     K = 5,
#     mu,
#     TARGET = 0.33,
#     cutoff = 0.95,
#     C = 3,
#     T_assess = 3,
#     cycle_max = 3,
#     CFO = FALSE,
#     cohort_dlt = c(0, 0, 0, 1, 0, 1, 1, 0, 1)
# ) {
#   out_list <- vector("list", max_cohort)
#   
#   for (cc in 1:max_cohort) {
#     admin_cc <- build_admin_fixed_y(
#       n_cohort = cc,
#       patients = patients,
#       K = K,
#       mu = mu,
#       TARGET = TARGET,
#       cutoff = cutoff,
#       C = C,
#       T_assess = T_assess,
#       cycle_max = cycle_max,
#       CFO = CFO,
#       cohort_dlt = cohort_dlt
#     )
#     
#     out_list[[cc]] <- admin_cc
#     
#     cat("\n====================\n")
#     cat("Completed through cohort", cc, "\n")
#     print(admin_cc)
#     
#     step <- decide_next_step_PRIDE(
#       admin = admin_cc,
#       patients = patients,
#       K = K,
#       mu = mu,
#       TARGET = TARGET,
#       cutoff = cutoff,
#       C = C,
#       T_assess = T_assess,
#       cycle_max = cycle_max,
#       CFO = CFO
#     )
#     
#     cat("\nInterim decision after cohort", cc, "\n")
#     cat("cur_dose           =", step$cur_dose, "\n")
#     cat("next_dose          =", step$next_dose, "\n")
#     cat("next_cohort_ids    =", paste(step$next_cohort_ids, collapse = ", "), "\n")
#     cat("new_ids            =", paste(step$new_ids, collapse = ", "), "\n")
#     cat("retreat_ids        =", paste(step$retreat_ids, collapse = ", "), "\n")
#     cat("next_start_time    =", step$next_start_time, "\n")
#     cat("next_decision_time =", step$next_decision_time, "\n")
#     
#     if (isTRUE(step$stop)) {
#       cat("Trial would stop here.\n")
#       break
#     }
#   }
#   
#   invisible(out_list)
# }
# c_k <- c(0.05, 0.10, 0.15, 0.20, 0.25)
# mu  <- qlogis(c_k)
# 
# patients <- data.frame(
#   id = 1:30,
#   t_arrival = (0:29) * 2
# )
# 
# res_path <- test_fixed_y_path(
#   max_cohort = 9,
#   patients = patients,
#   K = 5,
#   mu = mu,
#   TARGET = 0.33,
#   cutoff = 0.95,
#   C = 3,
#   T_assess = 3,
#   cycle_max = 3,
#   CFO = TRUE,
#   cohort_dlt = c(0, 0, 0, 1, 0, 1, 1, 0, 1)
# )
