#R library 
list.of.packages <- c("tidyr", "ismev", "fExtremes","lubridate", "tidyverse", "nloptr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(tidyr)
library(ismev)
library(fExtremes)
library(lubridate)
library(tidyverse)
library(nloptr)

#                                     MRC_model
# Inputs:
# ObservedTs: the observed time series with one column 'Precipitation' provided by the user in Step 1;
# n_par: »n» parameter value provided by the user in Step 1

#Outputs:
#Dissagr_Ts: the disaggregated time serie

MRC_model <- function(ObservedTs, n_par = 0.344){
  
  #Routine 1
  ObservedTs = data.frame(ObsTs = ObservedTs[,1])
  starting_day = '2000-01-01'
  start_date <- as.Date(starting_day)
  
  # Generiamo la sequenza anno per anno, evitando il 29 febbraio
  date_series <- unlist(lapply(0:(nrow(ObservedTs)/365-1), function(i) {
    year_days <- seq.Date(as.Date(paste0((year(start_date)+1) + i, "-01-01")), 
                          as.Date(paste0((year(start_date)+1) + i, "-12-31")), 
                          by = "day")
    year_days[!(format(year_days, "%m-%d") == "02-29")]  # Escludiamo il 29 febbraio
  }))
  
  # Convertiamo in vettore di date
  ObservedTs = data.frame(Time=as.Date(date_series), ObsTs=ObservedTs$ObsTs)
  n_anni = length(unique(year(ObservedTs$Time)))
  ObservedTs$Year <- factor(year(ObservedTs$Time))
  
  # Maximum Annual Values & GEV
  Max <-pivot_longer(ObservedTs, 2,values_drop_na=T) %>%
    dplyr::group_by(name, Year) %>%
    dplyr::summarise(Max = max(value, na.rm=TRUE), .groups = "drop")%>%
    ungroup()
  
  ObservedTs = ObservedTs[,1:2] 
  a_par_values <- ((Max$Max / 24) / (24 ^ (n_par - 1)))*1.12 
  
  #2) fit the GEV distribution on the «a» values, using three fitting methods: 'mle', 'pwm', 'L-moments'
  gev_pars=data.frame(matrix(0, ncol = 3 , nrow = 1))
  
  #mle method
  par_gev <- gev.fit(a_par_values, show = FALSE)
  gev_pars[1,] <- par_gev$mle
  
  #probability weighted moment method
  par_gev <- gevFit(a_par_values, type = "pwm")
  gev_pars[2,1]  <- as.numeric(par_gev@fit$par.ests[2]) 
  gev_pars[2,2] <- as.numeric(par_gev@fit$par.ests[3]) 
  gev_pars[2,3] <- as.numeric(par_gev@fit$par.ests[1]) 
  
  names(gev_pars) <- c("loc", "scale", "shape")
  #4) estimate the Mean Annual Rainfall Intensity
  Mean_Intesity <- sum(ObservedTs[,2])/(nrow(ObservedTs)*24)
  
  #DDF «a» parameters
  Return_periods = c(2,4,6,8,10)
  
  A_Tr_mle <- list()
  for(i in 1:length(Return_periods)){
    A_Tr_mle[[i]] <- qgev(1-(1/Return_periods[i]), as.numeric(gev_pars[1,3]), as.numeric(gev_pars[1,1]), as.numeric(gev_pars[1,2]))}

  a_Tr_mle <- round(as.vector(t(do.call(rbind, A_Tr_mle))), 2)

  a_Tr <- data.frame(rbind(a_Tr_mle))
  colnames(a_Tr) =  paste0('Tr', Return_periods)
  
  
  #Routine 2
  n_RP <- length(Return_periods)
  Num_Iterations <- 50000
  lower_bounds = c(0.1, 0.01, 24)
  upper_bounds = c(0.7,0.1,2400)
  
  h_mm <- list()
  for(ii in 1:length(Return_periods)){
    h_mm[[ii]] <- a_Tr[1,ii]*((24)^(n_par-1))*(24)}
  h_mm <- round(as.vector(t(do.call(rbind, h_mm))), 2)
  names(h_mm) <- paste0('Tr', Return_periods)
  
  Kq <- function(q, alpha, c_beta, c_ln) {
    if (alpha < 0) {
      kqc <- c_beta* (q - 1) + c_ln * (q * q - q)
    } else if (alpha > 2) {
      kqc <- log2(c_beta / (c_beta - q)) - q * log2(c_beta / (c_beta - 1))
      kqc[q >= c_beta] <- Inf
    } else {
      kqc <- (c_beta / (alpha - 1)) * (q * (q ^ (alpha - 1) - 1))
    }
    kqlim <- q - 1
    valkq <- abs(kqc) < abs(kqlim)
    kqc[!valkq] <- kqlim[!valkq]
    if (any(q == 1)) {
      qm1 <- (q < 1)
      if (any(qm1)) {
        if (all(valkq[qm1])) {
          kqc[q == 1] <- 0
        }
      } else {
        kqc[q == 1] <- 0
      }
    }
    return(kqc)
  }
  momz <- function(qmax, alpha, c_beta, c_ln) {
    eta <- 1
    qc <- 0:qmax
    if (missing(c_ln)) {
      kqc <- Kq(qc, alpha, c_beta)
    } else {
      kqc <- Kq(qc, alpha, c_beta, c_ln)
    }
    ez <- rep(1, qmax)
    bcoef <- c(1, 1)
    for (k in 2:qmax) {
      ez1 <- ez[2:k]
      ez2 <- rev(ez1)
      bcoef <- c(1, bcoef[-length(bcoef)] + bcoef[-1], 1)
      bc <- bcoef[-c(1, length(bcoef))]
      ezaux <- ez1 * ez2 * bc * (2^(kqc[2:k] + rev(kqc[2:k])))
      if (k > 2) {
        eza <- sum(ezaux)
      } else {
        eza <- sum(sapply(ezaux, function(x) sum(x)))
      }
      ez[k + 1] <- (2^(-k)) / (1 - 2^(kqc[k + 1] - k + 1)) * eza
    }
    ez[!is.finite(ez)] <- Inf
    return(ez)
  }
  
  Equaz15 <- function(Cb,Cln,D,Return_periods, h_mm, MeanIntesity, n_par, Num_Iterations){
    Tr <- coeff_a <- list()
    for(j in 1:length(Return_periods)){
      Tr[[j]] <-  Return_periods[[j]]*(365)*24
      coeff_a[[j]] <- (h_mm[[j]]/(24))/((24)^(n_par-1))  
    }
    res_Tr <- round(as.vector(t(do.call(rbind, Tr))), 2)
    res_a <- round(as.vector(t(do.call(rbind, coeff_a))), 2)
    names(res_Tr) = names(res_a) <- paste0('Tr', Return_periods)
    
    durate <- c(1/3,1,2,6,24)
    erre <- D/durate
    nr <- length(erre)
    
    df_i <- list()
    for (ii in 1:length(Return_periods)) {
      df_i[[paste0("i", ii, "f")]] <- rep(0, nr)
    }
    
    #Calcolo rz inserendo i parametri Cb e Cln
    qstar <- (1 - Cb)/Cln
    alpha <- -1 #alpha>2 for log-exp MF, alpha<0 for Beta-LN MF
    MomentOrder1 <- round(qstar/2,0)
    Ez <- momz(MomentOrder1, alpha, Cb, Cln)
    rz <-Ez[MomentOrder1+1]^(1/(Cb*(MomentOrder1-1)+Cln*(MomentOrder1^2-MomentOrder1)))
    
    #Calcolo i* e T* eq. 16b### T_star non ha senso perchè mi vengono valori altissimi
    i_star<-  (erre*rz)^(2-Cb-Cln)
    rapporto <- ((1-Cb)^2)/Cln
    T_star <- durate*((2*pi*2*log(erre*rz)*rapporto)^0.5)*(erre*rz)^(rapporto+Cb)
    
    # Calcolo per ciascuna Tr
    for (k in 1:length(Return_periods)) {
      # Calcolo per T_k
      argo <- 1 - (((erre * rz) ^ Cb) * D) / (erre * as.numeric(res_Tr[k]))
      i0 <- ((erre * rz) ^ (Cb - Cln)) * exp(((2 * Cln * log(erre * rz)) ^ 0.5) * qnorm(argo))
      i0 <- i0 * MeanIntesity
      istar <- i_star * (res_Tr[k] / T_star) ^ (Cln / (1 - Cb))
      
      # Assegna i valori calcolati a i_kf in base alla condizione Tr <= T_star
      for (jj in 1:nr) {
        if (res_Tr[k] <= T_star[jj]) df_i[[paste0("i", k, "f")]][jj] <- i0[jj] else df_i[[paste0("i", k, "f")]][jj] <- istar[jj]
      }
    }
    
    df_idf <- list()
    for (p in 1:length(Return_periods)) {
      df_idf[[paste0("idf", p)]] <- rep(0, nr)
      df_idf[[p]] <- res_a[p]*(D/erre)^(n_par-1)
    }
    i_f <- data.frame(round(do.call(rbind, df_i), 2))
    idf <- data.frame(round(do.call(rbind, df_idf), 2))
    
    results <- list()
    results$erre <- erre
    results$i_f <- i_f
    results$idf <- idf
    return(results)
  }
  
  obj_fun <- function(par) {
    Cb <- par[1]; Cln <- par[2]; D <- par[3]
    res <- Equaz15(Cb,Cln,D,Return_periods,h_mm,MeanIntesity,n_par, Num_Iterations)
    
    i_f <- res$i_f
    idf <- res$idf
    out <- 0
    for (w in 1:length(Return_periods)) {
      out <- out + sum((i_f[w,] - idf[w,])^2)
    }
    return(out)
  }
  
  opt <- nloptr(x0 = c(lower_bounds[1], lower_bounds[2], lower_bounds[3]), obj_fun, 
                lb = c(lower_bounds[1], lower_bounds[2], lower_bounds[3]), ub = c(upper_bounds[1], upper_bounds[2], upper_bounds[3]), 
                opts = list(algorithm = "NLOPT_GN_ISRES",  maxeval=Num_Iterations, xtol_abs=1e-8, xtol_rel = 1e-8)) #

  #Routine 3
  ending_time_resolution = 15
  N_simulations = 1
  years <- length(unique(year(ObservedTs$Time)))
  
  Num_daysinTS = years*365
  n_minutes = 60
  steps_in_1hour = n_minutes/ending_time_resolution
  num_steps = 24*steps_in_1hour
  
  b = 2
  L <- 1:17
  value <- b^L
  
  res0 <- data.frame(cbind(L, value))
  
  for (i in 1:(nrow(res0))) {
    if (num_steps >= res0[i, "value"] && num_steps <= res0[i + 1, "value"]) {
      L <- res0[i+1, "L"]
      break
    }
  }
  
  start_datetime <- as.POSIXct(paste0(ObservedTs[1,1],"00:00:00"), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  
  # Generiamo la sequenza anno per anno, evitando il 29 febbraio
  time_series <- unlist(lapply(0:(years-1), function(i) {
    year_start <- as.POSIXct(paste0((year(ObservedTs[1,1])) + i, "-01-01 00:00:00"), tz = "UTC")
    year_end <- as.POSIXct(paste0((year(ObservedTs[1,1])) + i, "-12-31 23:45:00"), tz = "UTC")
    
    # Generiamo la sequenza di 15 minuti per l'anno
    year_intervals <- seq.POSIXt(year_start, year_end, by =  paste0(ending_time_resolution, " mins"))
    
    # Escludiamo il 29 febbraio
    year_intervals = year_intervals[!(format(year_intervals, "%m-%d") == "02-29")]
  }))
  
  # Convertiamo in formato POSIXct
  time_series <- as.POSIXct(time_series, tz = "UTC")
  
  Dissagr_Ts <- data.frame(matrix(0, ncol = N_simulations+1 , nrow = length(time_series)))
  colnames(Dissagr_Ts) <- c('Time',paste0('Sim',1:N_simulations))
  
  Dissagr_Ts[,1] <- time_series 
  
  b_eta = as.numeric(round(opt$solution[1],2))
  sigma = sqrt((2*(as.numeric(round(opt$solution[2],2))))/log(b))
  
  SerGdis <- numeric(Num_daysinTS * b^L)
  ind <- which(ObservedTs[, 2] > 0)
  
  set.seed(1234)
  seed0 = sample(1000,N_simulations)
  
  for(sim in 1:N_simulations){
    set.seed(seed0[sim])
    for (z in seq_along(ind)) { 
      l <- 1
      W_up <- b^b_eta * rbinom(b^l, 1, b^(-b_eta)) * b^(-sigma^2 * (log(b) / 2) + sigma * rnorm(b^l))
      
      if (sum(W_up) == 0) {
        repeat {
          l <- 1
          W_up <- b^b_eta * rbinom(b^l, 1, b^(-b_eta)) * b^(-sigma^2 * (log(b) / 2) + sigma * rnorm(b^l))
          if (sum(W_up) != 0) break
        }
      }
      
      for (l in 2:L) {
        W <- b^b_eta * rbinom(b^l, 1, b^(-b_eta)) * b^(-sigma^2 * (log(b) / 2) + sigma * rnorm(b^l))
        
        for (c in seq(1, length(W), 2)) {
          if (sum(W[c] + W[c + 1]) == 0) {
            repeat {
              random_numbers <- b^b_eta * rbinom(b^l, 1, b^(-b_eta)) * b^(-sigma^2 * (log(b) / 2) + sigma * rnorm(b^l))
              W[c] <- random_numbers[1]
              W[c + 1] <- random_numbers[2]
              
              if (sum(W[c] + W[c + 1])  != 0) break
            }
          }
        }
        
        for (j in 1:b^l) {
          W[j] <- W[j] * W_up[ceiling(j / b)]
        }
        W_up <- W
      }
      
      W = W/sum(W)
      SerGdis[(b^L * (ind[z] - 1) + 1): (b^L * ind[z])] <- ObservedTs[ind[z], 2] * W 
    }
    Dissagr_Ts[, sim+1] <- diff(approx((0:(Num_daysinTS * b^L)) / b^L, c(0, cumsum(SerGdis)), (0:nrow(Dissagr_Ts)) / (num_steps))$y)
  }
  
  res <- list()
  res$Dissagr_Ts <- Dissagr_Ts

return(res)
}

daily_data = read.csv('Station1.cosmos2s.daily.csv')
res = MRC_model(ObservedTs,  n_par = 0.344)
res$Dissagr_Ts
