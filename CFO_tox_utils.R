# source("utilities.R")
# library(magrittr)
# library(BOIN)

# posterior probability of pj >= phi given data
post.prob.fn <- function(phi, y, n, alp.prior=0.1, bet.prior=0.1){
    alp <- alp.prior + y 
    bet <- bet.prior + n - y
    1 - pbeta(phi, alp, bet)
}



Odds.samples <- function(y1, n1, y2, n2, alp.prior, bet.prior){
    alp1 <- alp.prior + y1
    alp2 <- alp.prior + y2
    bet1 <- bet.prior + n1 - y1
    bet2 <- bet.prior + n2 - y2
    sps1 <- c()
    sps2 <- c()
    while (length(sps1)<10000){
        sp1 <- rbeta(1, alp1, bet1)
        sp2 <- rbeta(1, alp2, bet2)
        if (sp1 <= sp2){
            sps1 <- c(sp1, sps1)
            sps2 <- c(sp2, sps2)
        }
    }
    
    list(sps1=sps1, sps2=sps2)
}

prob.int <- function(phi, y1, n1, y2, n2, alp.prior, bet.prior){
    alp1 <- alp.prior + y1
    alp2 <- alp.prior + y2
    bet1 <- bet.prior + n1 - y1
    bet2 <- bet.prior + n2 - y2
    fn.min <- function(x){
        dbeta(x, alp1, bet1)*(1-pbeta(x, alp2, bet2))
    }
    fn.max <- function(x){
        pbeta(x, alp1, bet1)*dbeta(x, alp2, bet2)
    }
    eps <- 1e-5
    const.min <- integrate(fn.min, lower=0+eps, upper=1-eps)$value
    const.max <- integrate(fn.max, lower=0+eps, upper=1-eps)$value
    p1 <- integrate(fn.min, lower=0+eps, upper=phi)$value/const.min
    p2 <- integrate(fn.max, lower=0+eps, upper=phi)$value/const.max
    
    list(p1=p1, p2=p2)
}


OR.values <- function(phi, y1, n1, y2, n2, alp.prior, bet.prior, type){
    ps <- prob.int(phi, y1, n1, y2, n2, alp.prior, bet.prior)
    if (type=="L"){
        pC <- 1 - ps$p2
        pL <- 1 - ps$p1
        oddsC <- pC/(1-pC)
        oddsL <- pL/(1-pL)
        OR <- oddsC*oddsL
        
    }else if (type=="R"){
        pC <- 1 - ps$p1
        pR <- 1 - ps$p2
        oddsC <- pC/(1-pC)
        oddsR <- pR/(1-pR)
        OR <- (1/oddsC)/oddsR
    }
    return(OR)
}

All.OR.table <- function(phi, n1, n2, type, alp.prior, bet.prior){
   ret.mat <- matrix(rep(0, (n1+1)*(n2+1)), nrow=n1+1)
   for (y1cur in 0:n1){
       for (y2cur in 0:n2){
           ret.mat[y1cur+1, y2cur+1] <- OR.values(phi, y1cur, n1, y2cur, n2, alp.prior, bet.prior, type)
       }
   }
   ret.mat
}

# compute the marginal prob when lower < phiL/phiC/phiR < upper
# i.e., Pr(Y=y|lower<phi<upper)
margin.phi <- function(y, n, lower, upper){
    C <- 1/(upper-lower)
    fn <- function(phi) {
       dbinom(y, n, phi)*C
    }
    integrate(fn, lower=lower, upper=upper)$value
}

# Obtain the table of marginal distribution of (y1, y2) 
# after intergrate out (phi1, phi2)
# under H0 and H1
# H0: phi1=phi, phi < phi2 < 2phi
# H1: phi2=phi, 0   < phi1 < phi
margin.ys.table <- function(n1, n2, phi, hyperthesis){
    if (hyperthesis=="H0"){
        p.y1s <- dbinom(0:n1, n1, phi)
        p.y2s <- sapply(0:n2, margin.phi, n=n2, lower=phi, upper=2*phi)
    }else if (hyperthesis=="H1"){
        p.y1s <- sapply(0:n1, margin.phi, n=n1, lower=0, upper=phi)
        p.y2s <- dbinom(0:n2, n2, phi)
    }
    p.y1s.mat <- matrix(rep(p.y1s, n2+1), nrow=n1+1)
    p.y2s.mat <- matrix(rep(p.y2s, n1+1), nrow=n1+1, byrow=TRUE)
    margin.ys <- p.y1s.mat * p.y2s.mat
    margin.ys
}

# Obtain the optimal gamma for the hypothesis test
optim.gamma.fn <- function(n1, n2, phi, type, alp.prior, bet.prior){
    OR.table <- All.OR.table(phi, n1, n2, type, alp.prior, bet.prior) 
    ys.table.H0 <- margin.ys.table(n1, n2, phi, "H0")
    ys.table.H1 <- margin.ys.table(n1, n2, phi, "H1")
    
    argidx <- order(OR.table)
    sort.OR.table <- OR.table[argidx]
    sort.ys.table.H0 <- ys.table.H0[argidx]
    sort.ys.table.H1 <- ys.table.H1[argidx]
    n.tol <- length(sort.OR.table)
    
    if (type=="L"){
        errs <- rep(0, n.tol-1)
        for (i in 1:(n.tol-1)){
            err1 <- sum(sort.ys.table.H0[1:i])
            err2 <- sum(sort.ys.table.H1[(i+1):n.tol])
            err <- err1 + err2
            errs[i] <- err
        }
        min.err <- min(errs)
        if (min.err > 1){
            gam <- 0
            min.err <- 1
        }else {
            minidx <- which.min(errs)
            gam <- sort.OR.table[minidx]
        }
    }else if (type=='R'){
        errs <- rep(0, n.tol-1)
        for (i in 1:(n.tol-1)){
            err1 <- sum(sort.ys.table.H1[1:i])
            err2 <- sum(sort.ys.table.H0[(i+1):n.tol])
            err <- err1 + err2
            errs[i] <- err
        }
        min.err <- min(errs)
        if (min.err > 1){
            gam <- 0
            min.err <- 1
        }else {
            minidx <- which.min(errs)
            gam <- sort.OR.table[minidx]
        }
        
    }
    list(gamma=gam, min.err=min.err)
}

make.decision.CFO.fn <- function(phi, cys, cns, alp.prior, bet.prior, cover.doses, diag=FALSE){
    if (cover.doses[2] == 1){
        return(1)
    }else{
        if (is.na(cys[1]) & (cover.doses[3]==1)){
            return(2)
        }else  if (is.na(cys[1]) & (!(cover.doses[3]==1))){
           gam2 <- optim.gamma.fn(cns[2], cns[3], phi, "R", alp.prior, bet.prior)$gamma 
           OR.v2 <- OR.values(phi, cys[2], cns[2], cys[3], cns[3], alp.prior, bet.prior, type="R")
           if (OR.v2>gam2){
               return(3)
           }else{
               return(2)
           }
        }else  if (is.na(cys[3]) | (cover.doses[3]==1)){
           gam1 <- optim.gamma.fn(cns[1], cns[2], phi, "L", alp.prior, bet.prior)$gamma 
           OR.v1 <- OR.values(phi, cys[1], cns[1], cys[2], cns[2], alp.prior, bet.prior, type="L")
           if (OR.v1>gam1){
               return(1)
           }else{
               return(2)
           }
            
        }else  if (!(is.na(cys[1]) | is.na(cys[3]) | cover.doses[3]==1)){
           gam1 <- optim.gamma.fn(cns[1], cns[2], phi, "L", alp.prior, bet.prior)$gamma 
           gam2 <- optim.gamma.fn(cns[2], cns[3], phi, "R", alp.prior, bet.prior)$gamma 
           OR.v1 <- OR.values(phi, cys[1], cns[1], cys[2], cns[2], alp.prior, bet.prior, type="L")
           OR.v2 <- OR.values(phi, cys[2], cns[2], cys[3], cns[3], alp.prior, bet.prior, type="R")
           v1 <- OR.v1 > gam1
           v2 <- OR.v2 > gam2
           if (v1 & !v2){
               return(1)
           }else if (!v1 & v2){
               return(3)
           }else{
               return(2)
           }
        }
    }
}

overdose.fn <- function(phi, add.args=list()){
    CV <- add.args$CV
    if (is.null(CV)){
        CV <- 0.95
    }
    y <- add.args$y
    n <- add.args$n
    alp.prior <- add.args$alp.prior
    bet.prior <- add.args$bet.prior
    pp <- post.prob.fn(phi, y, n, alp.prior, bet.prior)
    #print(c(phi, y, n, alp.prior, bet.prior))
    if ((pp >= CV) & (add.args$n>=3)){
        return(TRUE)
    }else{
        return(FALSE)
    }
}

# Simulation function for CFO
CFO.simu.fn <- function(phi, p.true, ncohort=12, init.level=1, 
                              cohortsize=1, add.args=list()){
    # phi: Target DIL rate
    # p.true: True DIL rates under the different dose levels
    # ncohort: The number of cohorts
    # cohortsize: The sample size in each cohort
    # alp.prior, bet.prior: prior parameters
    earlystop <- 0
    ndose <- length(p.true)
    cidx <- init.level
    
    tys <- rep(0, ndose) # number of responses for different doses.
    tns <- rep(0, ndose) # number of subject for different doses.
    tover.doses <- rep(0, ndose) # Whether each dose is overdosed or not, 1 yes
    
    
    
    
    for (i in 1:ncohort){
        pc <- p.true[cidx] 
        
        # sample from current dose
        cres <- rbinom(cohortsize, 1, pc)
        
        # update results
        tys[cidx] <- tys[cidx] + sum(cres)
        tns[cidx] <- tns[cidx] + cohortsize
        
        
        
        cy <- tys[cidx]
        cn <- tns[cidx]
        
        add.args <- c(list(y=cy, n=cn, tys=tys, tns=tns, cidx=cidx), add.args)
        
        if (overdose.fn(phi, add.args)){
            tover.doses[cidx:ndose] <- 1
        }
        
        if (tover.doses[1] == 1){
            earlystop <- 1
            break()
        }
        
        
        # the results for current 3 dose levels
        if (cidx!=1){
            cys <- tys[(cidx-1):(cidx+1)]
            cns <- tns[(cidx-1):(cidx+1)]
            cover.doses <- tover.doses[(cidx-1):(cidx+1)]
            #cover.doses <- c(0, 0, 0) # No elimination rule
        }else{
            cys <- c(NA, tys[1:(cidx+1)])
            cns <- c(NA, tns[1:(cidx+1)])
            cover.doses <- c(NA, tover.doses[1:(cidx+1)])
            #cover.doses <- c(NA, 0, 0) # No elimination rule
        }
        
        idx.chg <- make.decision.CFO.fn(phi, cys, cns, add.args$alp.prior, add.args$bet.prior, cover.doses) - 2
            
        
        cidx <- idx.chg + cidx
        
    }
    
    
    if (earlystop==0){
        MTD <- select.mtd(phi, tns, tys)$MTD
    }else{
        MTD <- 99
    }
    list(MTD=MTD, dose.ns=tns, DLT.ns=tys, ps=p.true, target=phi, over.doses=tover.doses, total.time=-1)
    #list(MTD=MTD, dose.ns=tns, DLT.ns=tys, p.true=p.true, target=phi, over.doses=tover.doses)
    # I change it to make it compatible with TITE-CFO simulation
}

