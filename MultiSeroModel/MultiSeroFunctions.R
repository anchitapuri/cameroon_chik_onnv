
# ---- Generate infection status matrix
inf_matrix <- function(N_pathogen, pres=rep(1,N_pathogen)){
  
  # list of possible outcomes for each pathogen
  combos <- list()
  for(c in 1:N_pathogen) combos[[c]] <- c(0,1)
  
  # matrix of all possible infection status combinations
  m <- expand.grid(combos)
  colnames(m) <- letters[1:N_pathogen]
  
  # remove positives of absent pathogens
  if(sum(pres)<N_pathogen){
    for(abs in which(pres==0)) m <- m[m[,abs]==0, ]
  }
  
  return(m)
}

# ---- Build data for inputting into stan model 
prepare_multiplex_sero_data <- function(
    data,
    pathogens,
    present_pathogens) {
  
  # Validate inputs
  if (!all(present_pathogens %in% pathogens)) {
    stop("All present_pathogens must be in the pathogens list")
  }
  
  # Reorder pathogens (present first, then non-present)
  nonpres <- pathogens[!pathogens %in% present_pathogens]
  pathogens_ordered <- c(present_pathogens, nonpres)
  
  # Compile data for model fitting
  stan_data <- list()
  
  # Antibody titer data (log-transformed)
  stan_data$y <- cbind(data[, c(present_pathogens, nonpres)])
  stan_data$y <- stan_data$y
  
  # Pathogen presence indicators
  stan_data$pres <- c(rep(1, length(present_pathogens)), rep(0, length(nonpres)))
  
  # Basic dimensions
  stan_data$N <- nrow(stan_data$y)
  stan_data$nP <- ncol(stan_data$y)
  stan_data$nPp <- sum(stan_data$pres)
  
  # Infection matrix and related variables
  stan_data$infM <- inf_matrix(stan_data$nP, pres = stan_data$pres)
  stan_data$nC <- nrow(stan_data$infM)
  
  # Create indexes for model fitting
  npos <- rowSums(stan_data$infM)
  wpos <- wneg <- matrix(0, ncol = stan_data$nP, nrow = stan_data$nC)
  
  for (c in 1:nrow(stan_data$infM)) {
    for (p in 1:stan_data$nP) {
      if (npos[c] > 0) {
        wpos[c, 1:npos[c]] <- which(stan_data$infM[c, ] == 1)
      }
      if (npos[c] < stan_data$nP) {
        wneg[c, 1:(stan_data$nP - npos[c])] <- which(stan_data$infM[c, ] == 0)
      }
    }
  }
  
  stan_data$npos <- npos  # N pos pathogens per infection status
  stan_data$wpos <- wpos  # index pos pathogens per infection status
  stan_data$wneg <- wneg  # index neg pathogens per infection status
  
  # Return list with data and model
  return(list(
    data = stan_data,
    pathogens = pathogens_ordered
  ))
}



#--- chain starting values
init <- function(data, nChains){
  ii <- init <- list()
  for(i in 1:nChains){
    init$sero <- array(runif(data$nPp, 0.2, 0.8))
    init$sd0 <- runif(data$nP, 0.4, 0.8)
    init$sd1 <- runif(1, 0.3, 0.6)
    init$mu0 <- runif(data$nP, 3, 5)
    init$mu1 <- runif(data$nPp, 1.5, 4.5)
    init$phi <- runif((data$nP * data$nPp - (data$nPp)), 0.01, 0.5)
    init$rho00 <- runif(1, 0.4, 0.7)
    ii[[i]] <- init
  } 
  
  return(ii)
}

init_diffSds <- function(data, nChains){
  ii <- vector("list", nChains)
  for(i in 1:nChains){
    ii[[i]] <- list(
      sero = c(
        runif(1, 0.22, 0.28),
        runif(1, 0.12, 0.18)
      ),
      mu0 = c(
        runif(1, 4.3, 5.0),
        runif(1, 6.0, 6.8),
        runif(1, 5.7, 6.4)
      ),
      mu1 = c(
        runif(1, 1.8, 2.2),
        runif(1, 2.7, 3.3)
      ),
      sd0 = c(
        runif(1, 0.75, 1.00),  # ONNV - narrower
        runif(1, 0.50, 0.60),  # CHIK - narrower
        runif(1, 1.10, 1.30)   # MAYV - wider
      ),
      sd1 = runif(1, 0.25, 0.45),
      phi = runif((data$nP * data$nPp - data$nPp), 0.3, 0.7),
      rho00 = runif(1, 0.25, 0.35)
    )
  }
  ii
}


#----- Extract prevalence estimates
extract_sero <- function(chains, data, pathogens){
  
  sero <- data.frame(pathogen=pathogens, med=NA, ciL=NA, ciU=NA)
  
  for(p in 1:data$nP) sero[p,2:4] <- quantile(chains[,paste('seroAll[', paste(p,']', sep=''), sep='')], c(0.5,0.025,0.975))
  sero <- sero[!sero$med==0, ]
  
  return(sero)
}


#----- Extract covariance matrices per iteration
extract_covM <- function(chains, data){
  
  iter <- length(chains$lp__)
  covM <- list()
  for(i in 1:iter){
    covM[[i]] <- list()
    for(c in 1:data$nC){
      x <- matrix(NA, ncol=data$nP, nrow=data$nP)
      for(p in 1:data$nP) for(p2 in 1:data$nP){
        
        y <- paste(paste(c,p,sep=','),p2,sep=',')
        x[p,p2] <- chains[i,paste(paste('covM[',y,sep=''),']',sep='')]
        
      }
      covM[[i]][[c]] <- x
    }
  }
  return(covM)
}


#----- Plot gaussian distribution fits
plot_fits <- function(chains, data, pathogens, show_crossreactive_for = NULL){
  
  iter <- length(chains$lp__)
  covM <- extract_covM(chains, data)
  
  # simulate multivariate gaussians per combination
  yy <- yp <- list()
  for(p in 1:data$nP) yp[[p]] <- matrix(NA, nrow=512, ncol=iter)
  ypN <- ypP <- ypNcross <- yp
  for(i in 1:iter){
    for(c in 1:data$nC){
      
      # simulate gaussion for combination c, iteration i
      g <- paste('mu[',c, sep='')
      nn <- ceiling(sum(data$N*chains[,paste('theta[',paste(c,']',sep=''), sep='')][i]))
      muu <- vector()
      for(p in 1:data$nP) muu[p] <- chains[i,paste(paste(g,p,sep=','),']',sep='')]
      yy[[c]] <- as.data.frame(rmvnorm(nn, mean=muu,sigma=covM[[i]][[c]]))
      yy[[c]]$C <- c
    }
    
    # density distributions per pathogen
    yc <- do.call('rbind',yy)
    for(p in 1:data$nP){
      
      yp[[p]][,i] <- density(yc[,p], bw=0.01, from=-2.5, to=12)$y
      
      # Check if we should show cross-reactive for this pathogen
      show_cr <- !is.null(show_crossreactive_for) && p %in% show_crossreactive_for
      
      if(data$pres[p]==1){
        z <- which(data$infM[,p]==1)
        pw <- vector()
        for(s in 1:length(z)) pw[s] <- chains[i,paste(paste('theta[',paste(z[s]),sep=''),']',sep='')] 
        propP <- sum(data$N*pw)/data$N
        
        # Positive for this pathogen
        ypP[[p]][,i] <- density(yc[yc$C %in% which(data$infM[,p]==1),p], bw=0.01, from=-2.5, to=12)$y * propP
        
        if(show_cr){
          # Split negatives into true negatives and cross-reactive negatives
          neg_idx <- which(data$infM[,p]==0)
          
          # Cross-reactive negatives: negative for pathogen p but positive for at least one other pathogen
          cross_reactive_idx <- neg_idx[rowSums(data$infM[neg_idx, , drop=FALSE]) > 0]
          # True negatives: negative for all pathogens
          true_neg_idx <- neg_idx[rowSums(data$infM[neg_idx, , drop=FALSE]) == 0]
          
          # Calculate proportions
          pw_cross <- vector()
          pw_true <- vector()
          if(length(cross_reactive_idx) > 0){
            for(s in 1:length(cross_reactive_idx)) pw_cross[s] <- chains[i,paste(paste('theta[',paste(cross_reactive_idx[s]),sep=''),']',sep='')]
            propNcross <- sum(data$N*pw_cross)/data$N
          } else {
            propNcross <- 0
          }
          if(length(true_neg_idx) > 0){
            for(s in 1:length(true_neg_idx)) pw_true[s] <- chains[i,paste(paste('theta[',paste(true_neg_idx[s]),sep=''),']',sep='')]
            propNtrue <- sum(data$N*pw_true)/data$N
          } else {
            propNtrue <- 0
          }
          
          # Density for cross-reactive negatives
          if(length(cross_reactive_idx) > 0){
            ypNcross[[p]][,i] <- density(yc[yc$C %in% cross_reactive_idx,p], bw=0.01, from=-2.5, to=12)$y * propNcross
          } else {
            ypNcross[[p]][,i] <- rep(0, 512)
          }
          
          # Density for true negatives
          if(length(true_neg_idx) > 0){
            ypN[[p]][,i] <- density(yc[yc$C %in% true_neg_idx,p], bw=0.01, from=-2.5, to=12)$y * propNtrue
          } else {
            ypN[[p]][,i] <- rep(0, 512)
          }
        } else {
          # No cross-reactive split, lump all negatives together
          ypN[[p]][,i] <- density(yc[yc$C %in% which(data$infM[,p]==0),p], bw=0.01, from=-2.5, to=12)$y *(1-propP)
          ypNcross[[p]][,i] <- rep(0, 512)
        }
        
      }else{
        ypN[[p]][,i] <- density(yc[yc$C %in% which(data$infM[,p]==0),p], bw=0.01, from=-2.5, to=12)$y
        ypNcross[[p]][,i] <- rep(0, 512)
      }
    }
  }
  
  # quantiles of density distributions
  dpq <- dpqP <- dpqN <- dpqNcross <- list()
  titer <- density(yc[,1], bw=0.01, from=-2.5, to=12)$x
  for(p in 1:data$nP){
    dpq[[p]] <- as.data.frame(rowQuantiles(yp[[p]], probs=c(0.5,0.025,0.975)))
    dpqN[[p]] <- as.data.frame(rowQuantiles(ypN[[p]], probs=c(0.5,0.025,0.975)))
    dpqNcross[[p]] <- as.data.frame(rowQuantiles(ypNcross[[p]], probs=c(0.5,0.025,0.975)))
    dpqP[[p]] <- as.data.frame(rowQuantiles(ypP[[p]], probs=c(0.5,0.025,0.975)))
    dpq[[p]]$titer <- dpqN[[p]]$titer <- dpqNcross[[p]]$titer <- dpqP[[p]]$titer <- titer
    dpq[[p]]$pathogen <- dpqN[[p]]$pathogen <- dpqNcross[[p]]$pathogen <- dpqP[[p]]$pathogen <- pathogens[p]
  }
  dpq <- do.call('rbind',dpq)
  dpqN <- do.call('rbind',dpqN)
  dpqNcross <- do.call('rbind',dpqNcross)
  dpqP <- do.call('rbind',dpqP)
  colnames(dpq)[1:3] <- colnames(dpqN)[1:3] <- colnames(dpqNcross)[1:3] <- colnames(dpqP)[1:3] <- c('med','ciL','ciU')
  
  # compile data for plotting 
  dta <- as.data.frame(data$y)
  colnames(dta) <- pathogens
  dta <- tidyr::gather(dta, key='pathogen', value='t')
  
  # Set factor levels to match input order
  dta$pathogen <- factor(dta$pathogen, levels = pathogens)
  dpq$pathogen <- factor(dpq$pathogen, levels = pathogens)
  dpqN$pathogen <- factor(dpqN$pathogen, levels = pathogens)
  dpqNcross$pathogen <- factor(dpqNcross$pathogen, levels = pathogens)
  dpqP$pathogen <- factor(dpqP$pathogen, levels = pathogens)
  
  # overall fit
  fitD <- ggplot()+ geom_histogram(data=dta, aes(t,y=..density..), bins=150, fill='#dedede', col='#dedede')+
    theme_minimal() + 
    theme(text=element_text(size=20))+
    geom_line(data=dpq, aes(titer,med),col='#008b7f')+ facet_wrap(~pathogen, scales='free_y')+ xlab('titer')+
    geom_ribbon(data=dpq, aes(x=titer,y=med,ymin=ciL,ymax=ciU), fill='#02d0bf', alpha=0.4)
  
  # pos-neg fit with cross-reactive negatives
  fitDPN <- ggplot()+ geom_histogram(data=dta, aes(t,y=..density..), bins=150, fill='#dedede', col='#dedede')+
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.7),
      axis.ticks.x = element_line(color = "black", size = 0.5),
      axis.ticks.y = element_line(color = "black", size = 0.5),
      axis.text = element_text(size = 20),
      axis.text.x = element_text(size = 20),
      axis.title = element_text(size = 24),
      strip.text = element_text(size = 22),
      aspect.ratio = 1,
      legend.title = element_blank(),
      legend.text = element_text(size = 18)
    )+
      # Lines
      geom_line(data = dpqN,
                aes(titer, med, colour = "True negative"),
                linewidth = 1) +
      geom_line(data = dpqNcross,
                aes(titer, med, colour = "Cross-reactive negative"),
                linewidth = 1) +
      geom_line(data = dpqP,
                aes(titer, med, colour = "Positive"),
                linewidth = 1) +

      # Ribbons 
      geom_ribbon(data = dpqN,
                  aes(x = titer, y = med,
                      ymin = ciL, ymax = ciU,
                      fill = "True negative"),
                  alpha = 0.3) +
      geom_ribbon(data = dpqNcross,
                  aes(x = titer, y = med,
                      ymin = ciL, ymax = ciU,
                      fill = "Cross-reactive negative"),
                  alpha = 0.3) +
      geom_ribbon(data = dpqP,
                  aes(x = titer, y = med,
                      ymin = ciL, ymax = ciU,
                      fill = "Positive"),
                  alpha = 0.5) +
    facet_wrap(~pathogen,
        labeller = labeller(
          pathogen = c(
            "ONNV_VLP_log" = "ONNV",
            "CHIKV_sE2_log" = "CHIKV",
            "MAYV_E2_log"  = "MAYV"
          )
        )
      ) +
    labs(x = "Log(titer)", y = "Density") +
    # Manual colours
    scale_colour_manual(
      values = c(
        "True negative" = "#021d37",
        "Cross-reactive negative" = "#035a6f",
        "Positive" = "#530227"
      )
    ) +
    scale_fill_manual(
      values = c(
        "True negative" = "#043565",
        "Cross-reactive negative" = "#028eb1",
        "Positive" = "#c7035b"
      )
    )+ guides(fill = "none",  colour = "none")
  
  # return plots
  return(list(fit=fitD,fitPN=fitDPN)) 
}


#----- Extract gaussian means
extract_mu <- function(chains, data, pathogens){
  
  
  # label combination positives
  pos <- rep('neg', data$nC)
  for(c in 1:data$nC){
    np <- sum(data$infM[c,])
    if(np==1) pos[c] <- pathogens[which(data$infM[c,]==1)]
    else if(np>1) pos[c] <- paste(pathogens[which(data$infM[c,]==1)], sep='&', collapse="&")
  }
  
  # all gaussian means
  mus0 <- data.frame(pg=rep(NA,length(which(data$infM==0))), pos=NA, med=NA, ciL=NA, ciU=NA)
  mus1 <- data.frame(pg=rep(NA,length(which(data$infM==1))), pos=NA, med=NA, ciL=NA, ciU=NA)
  ix0 <- 1
  ix1 <- 1
  for(c in 1:data$nC) for(p in 1:data$nP){
    if(data$infM[c,p]==0){
      mus0$pg[ix0] <- pathogens[p]
      mus0$pos[ix0] <- pos[c] 
      y <- paste(c,p, sep=',')
      mus0[ix0,3:5] <- quantile(chains[,paste(paste('mu[',y,sep=''),']',sep='')], c(0.5,0.025,0.975))
      ix0 <- ix0+1
    }else{
      mus1$pg[ix1] <- pathogens[p]
      mus1$pos[ix1] <- pos[c] 
      y <- paste(c,p, sep=',')
      mus1[ix1,3:5] <- quantile(chains[,paste(paste('mu[',y,sep=''),']',sep='')], c(0.5,0.025,0.975))
      ix1 <- ix1+1
    }
  }
  colnames(mus0)[1] <- 'antigen'
  colnames(mus1)[1] <- 'antigen'
  
  return(list(mus0=mus0, mus1=mus1))
}

extract_phi <- function(chains, data, pathogens){
  
  phi <- data.frame(pos=NA, neg=NA, med=NA, ciL=NA, ciU=NA)
  ind <- 1
  for(p in 1:data$nPp) for(p2 in 1:data$nP){
    if(!p==p2){
      phi[ind,1:2] <- c(pathogens[p], pathogens[p2])
      y <- str_replace_all(toString(c(p,p2))," ","")
      phi[ind,3:5] <- quantile(chains[,paste('CR[', paste(y, ']', sep=''), sep='')], c(0.5,0.025,0.975))
      ind <- ind+1
    }
  }
  
  rho <- data.frame(pars=c('rho00'), med=NA, ciL=NA, ciU=NA)
  rho[1,2:4] <- quantile(chains[,paste('rho00')], c(0.5,0.025,0.975))
  
  return(list(phi=phi,rho=rho))
  
}



# --- Plot cross reactivity - ONNV vs CHIK
plot_titer_increases_comparison <- function(phi_df, mu_mus1) {
  
  # Extract values for each combination
  # ONNV homologous
  onnv_mu1_med <- mu_mus1 %>% filter(pos == "ONNV_VLP_log") %>% pull(med)
  onnv_mu1_ciL <- mu_mus1 %>% filter(pos == "ONNV_VLP_log") %>% pull(ciL)
  onnv_mu1_ciU <- mu_mus1 %>% filter(pos == "ONNV_VLP_log") %>% pull(ciU)
  
  # CHIK homologous
  chik_mu1_med <- mu_mus1 %>% filter(pos == "CHIKV_sE2_log") %>% pull(med)
  chik_mu1_ciL <- mu_mus1 %>% filter(pos == "CHIKV_sE2_log") %>% pull(ciL)
  chik_mu1_ciU <- mu_mus1 %>% filter(pos == "CHIKV_sE2_log") %>% pull(ciU)

  # ONNV → CHIK cross-reactive
  onnv_to_chik_phi_med <- phi_df %>% filter(pos == "ONNV_VLP_log", neg == "CHIKV_sE2_log") %>% pull(med)
  onnv_to_chik_phi_ciL <- phi_df %>% filter(pos == "ONNV_VLP_log", neg == "CHIKV_sE2_log") %>% pull(ciL)
  onnv_to_chik_phi_ciU <- phi_df %>% filter(pos == "ONNV_VLP_log", neg == "CHIKV_sE2_log") %>% pull(ciU)
  
  # CHIK → ONNV cross-reactive
  chik_to_onnv_phi_med <- phi_df %>% filter(pos == "CHIKV_sE2_log", neg == "ONNV_VLP_log") %>% pull(med)
  chik_to_onnv_phi_ciL <- phi_df %>% filter(pos == "CHIKV_sE2_log", neg == "ONNV_VLP_log") %>% pull(ciL)
  chik_to_onnv_phi_ciU <- phi_df %>% filter(pos == "CHIKV_sE2_log", neg == "ONNV_VLP_log") %>% pull(ciU)
  
  # Calculate cross-reactive increases (phi × mu1)
  onnv_to_chik_increase_med <- onnv_to_chik_phi_med * onnv_mu1_med
  onnv_to_chik_increase_ciL <- onnv_to_chik_phi_ciL * onnv_mu1_ciL
  onnv_to_chik_increase_ciU <- onnv_to_chik_phi_ciU * onnv_mu1_ciU
  
  chik_to_onnv_increase_med <- chik_to_onnv_phi_med * chik_mu1_med
  chik_to_onnv_increase_ciL <- chik_to_onnv_phi_ciL * chik_mu1_ciL
  chik_to_onnv_increase_ciU <- chik_to_onnv_phi_ciU * chik_mu1_ciU
  
  # Create plotting data - 4 rows total
  plot_data <- tibble(
    antigen = c("ONNV_VLP_log", "ONNV_VLP_log", "CHIKV_sE2_log", "CHIKV_sE2_log"),
    response_type = c("Homologous", "Cross-reactive \n(φ × μ1)", 
                      "Homologous", "Cross-reactive \n(φ × μ1)"),
    increase = c(
      onnv_mu1_med,                  # ONNV homologous
      chik_to_onnv_increase_med,     # CHIK → ONNV cross-reactive
      chik_mu1_med,                  # CHIK homologous
      onnv_to_chik_increase_med      # ONNV → CHIK cross-reactive
    ),
    ciL = c(
      onnv_mu1_ciL,
      chik_to_onnv_increase_ciL,
      chik_mu1_ciL,
      onnv_to_chik_increase_ciL
    ),
    ciU = c(
      onnv_mu1_ciU,
      chik_to_onnv_increase_ciU,
      chik_mu1_ciU,
      onnv_to_chik_increase_ciU
    ),
    infecting_virus = c("ONNV", "CHIK", "CHIK", "ONNV")
  ) %>%
    mutate(
      response_type = factor(response_type, levels = c("Homologous", "Cross-reactive \n(φ × μ1)")),
      antigen_label = ifelse(antigen == "ONNV_VLP_log", "ONNV", "CHIK")
    )
  
  p <- ggplot(plot_data, aes(x = response_type, y = increase, fill = infecting_virus)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.8) +
    geom_errorbar(aes(ymin = ciL, ymax = ciU), 
                  position = position_dodge(width = 0.8), 
                  width = 0.25, linewidth = 0.8) +
    facet_wrap(~ antigen_label) +
    scale_fill_manual(
      values = c("ONNV" = "#c7035b", "CHIK" = "#028eb1"),
      name = "Infecting pathogen"
    )  +
    labs(
      x = "",
      y = "Antibody titer increase (log units)"
    ) +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.text.x = element_text(size = 20, angle = 20, hjust = 0.8),
      axis.text.y = element_text(size = 20),
      axis.title = element_text(size = 20),
      strip.text = element_text(size = 20),
      strip.background = element_rect(fill = "#ffffff"),
      legend.position = "bottom",
      legend.box.margin = margin(t = -20, r = 0, b = 0, l = 0),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      panel.grid.major.x = element_blank(),
    )
  
  return(p)
}

# --- Plot proportion positive 
plot_seroprevalence <- function(chains_df) {
  
  # Extract sero parameters
  sero_onnv <- chains_df$`sero[1]`
  sero_chik <- chains_df$`sero[2]`
  
  # Create summary data
  sero_data <- tibble(
    pathogen = c("ONNV", "CHIK"),
    med = c(median(sero_onnv), median(sero_chik)),
    ciL = c(quantile(sero_onnv, 0.025), quantile(sero_chik, 0.025)),
    ciU = c(quantile(sero_onnv, 0.975), quantile(sero_chik, 0.975))
  ) %>%
    mutate(
      pathogen = factor(pathogen, levels = c("ONNV", "CHIK"))
    )
  
  p <- ggplot(sero_data, aes(x = pathogen, y = med, fill = pathogen)) +
    geom_col(width = 0.6, alpha = 0.8) +
    geom_errorbar(aes(ymin = ciL, ymax = ciU), width = 0.25, linewidth = 0.8) +
    scale_fill_manual(
      values = c("ONNV" = "#c7035b", "CHIK" = "#028eb1"),
      guide = "none"
    ) +
    scale_y_continuous(limits = c(0, NA), labels = scales::percent) +
    labs(
      x = "",
      y = "Seroprevalence") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 24),
      axis.text.y = element_text(size = 24),
      axis.title = element_text(size = 20),
      strip.text = element_text(size = 20),
      strip.background = element_rect(fill = "#ffffff"),
      legend.position = "right",
      panel.grid.major.x = element_blank(),
      aspect.ratio = 0.6, 
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
  
  # Print the values
  cat("\nSeroprevalence estimates:\n")
  print(sero_data)
  
  return(p)
}

