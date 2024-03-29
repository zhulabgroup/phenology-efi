for(site in 1:nrow(coord_df)) {
  path_train<-paste0(path, focal_var,"/train/",site)
  path_results<-paste0(path,focal_var, "/results/", site)
  dir.create(path_results, recursive = T)
  
  ########################
  # Get new data
  ########################
  X_all<-read_csv(paste0(path_train,"/X.csv")) %>% as.matrix()
  Y_all<-read_csv(paste0(path_train,"/Y.csv"))%>% as.matrix()
  D_all<-read_csv(paste0(path_train,"/D.csv")) %>% as.matrix()
  P_all<-read_csv(paste0(path_train,"/P.csv")) %>% as.matrix()
  
  missing_id<-which(rowSums(is.na(cbind(Y_all,X_all)))!=0)
  if (length(missing_id)>0) {
    X_all<-X_all[-missing_id,,drop=F]
    Y_all<-Y_all[-missing_id,,drop=F]
    P_all<-P_all[-missing_id,,drop=F]
    D_all<-D_all[-missing_id,,drop=F]
  }
  
  train_id <- sample(1:nrow(X_all), round(nrow(X_all) * 8 / 10))
  valid_id <- setdiff(1:nrow(X_all), train_id)
  
  X_train <- X_all[train_id, , drop = F]
  Y_train <- Y_all[train_id, , drop = F]
  P_train <- P_all[train_id, , drop = F]
  D_train <- D_all[train_id, , drop = F]
  
  J_train<-date2doy(D_train)
  window_id<-which (tdistMat[J_train,date2doy(matrix(as.character(today)))]<=7)
  X_train <- X_train[window_id, , drop = F]
  Y_train <- Y_train[window_id, , drop = F]
  P_train <- P_train[window_id, , drop = F]
  D_train <- D_train[window_id, , drop = F]
  
  X_valid <- X_all[valid_id, , drop = F]
  Y_valid <- Y_all[valid_id, , drop = F]
  P_valid <- P_all[valid_id, , drop = F]
  D_valid <- D_all[valid_id, , drop = F]
  
  basisnumber<-min(basisnumber, nrow(X_all))
  if (basisnumber==nrow(X_all)) {
    X_basis <- X_all
    Y_basis <- Y_all
    P_basis <- P_all
    D_basis <- D_all
  } else {
    J_all<-date2doy(D_all)
    set.seed(42)
    cluster <- kmeans(cbind(P_all,J_all,X_all), basisnumber )
    cluster_id_sort <- seq(basisnumber )[order(cluster$size, decreasing = T)]
    basis_id <- rep(NA, basisnumber)
    for (i in 1:basisnumber) {
      # basis_id<-which(cl$cluster==cluster_id_sort[i])[length(which(cl$cluster==cluster_id_sort[i]))]
      basis_id[i] <- sample(which(cluster$cluster == cluster_id_sort[i]), 1)
    }
    X_basis <- X_all[basis_id, , drop = F]
    Y_basis <- Y_all[basis_id, , drop = F]
    P_basis <- P_all[basis_id, , drop = F]
    D_basis <- D_all[basis_id, , drop = F]
  }
  
  write_csv(as.data.frame(X_basis), paste0(path_results,"/X_basis", ".csv"))
  write_csv(as.data.frame(Y_basis), paste0(path_results,"/Y_basis", ".csv"))
  write_csv(as.data.frame(P_basis), paste0(path_results,"/P_basis", ".csv"))
  write_csv(as.data.frame(D_basis), paste0(path_results,"/D_basis", ".csv"))
  
  if (epoch == 1) {
    pars_id_prev <- numeric()
    pars_prev_updated <- matrix(NA, nrow = ndim + 4, ncol = 0)
    pars_var_prev_updated <- vector(mode = "list", length = 0)
    log_p_prev_updated <- numeric(0)
  }
  if (epoch > 1) {
    ########################
    # Get previous hyperparameters
    ########################
    # pars_results_prev<-paste0(path, "/results/", "LC", lc_code, "/", "DOY", doy_oi, "/", "EPOCH",epoch-1)
    # change this
    
    pars_prev <- as.matrix(read_csv(list.files(pars_results_prev, pattern="pars", recursive = T)))
    pars_id_prev <- colnames(pars_prev)
    pars_var_prev <- vector(mode = "list", length = num_part)
    for (i in 1:num_part) {
      pars_var_prev[[i]] <- as.matrix(nearPD(as.matrix(read_csv(list.files(pars_results_prev, pattern=paste0("pars_var_", i), recursive = T))))$mat)
    }
    
    ########################
    # Learn hyperparameters
    ########################
    Rprop_res <-
      foreach(
        i = 1:num_part,
        # .export = ls(globalenv()),
        .export = c(
          "phimax", "phimin", "vemax", "vemin", "taumax", "taumin", "rhomax", "rhomin",
          "pars_prev", "pars_var_prev", "distMat", "priors",
          "X_train", "P_train", "Y_train",
          "X_basis", "P_basis", "Y_basis",
          "GPSDM", "fmingrad_Rprop"
        ),
        .packages = c("Matrix", "LaplacesDemon", "RhpcBLASctl")
      ) %dopar% {
        blas_set_num_threads(1)
        omp_set_num_threads(1)
        
        # subset<-sample(1:nrow(X_new), nrow(X_new))
        # X_new<-X_new[subset,,drop=F]
        # P_new<-P_new[subset,,drop=F]
        # Y_new<-Y_new[subset,,drop=F]
        
        lpost <- function(p,sample_n) {
          minibatch <- sample(1:nrow(X_train), sample_n)
          return(GPSDM(pars = p, distMat = distMat, basisX = X_basis, basisP = P_basis, basisY = Y_basis, newX = X_train[minibatch, , drop = F], newP = P_train[minibatch, , drop = F], newY = Y_train[minibatch, , drop = F], priors_tr = list(E = pars_prev[, i, drop = F], V = pars_var_prev[[i]]), mode = c("optimize")))
        }
        print("start")
        
        fmingrad_Rprop(lpost, pars_prev[, i, drop = F] + matrix(rnorm(ndim + 3, 0, 0.1)), sample_n = min(nrow(X_train),100), maxcount = maxcount)
      }
    
    ########################
    # Estimate hyperparameter uncertainty
    ########################
    pars_prev_updated <- matrix(NA, nrow = ndim + 3, ncol = num_part)
    log_p_prev_updated <- rep(NA, num_part)
    for (i in 1:num_part) {
      pars_prev_updated[, i] <- Rprop_res[[i]]$xopt
      log_p_prev_updated[i] <- -Rprop_res[[i]]$fopt
    }
    pars_var_prev_updated <-
      foreach(
        i = 1:num_part,
        .export = c(
          "phimax", "phimin", "vemax", "vemin", "taumax", "taumin", "gammamax", "gammamin",
          "pars_prev", "pars_var_prev", "pars_prev_updated", "distMat", "priors",
          "X_train", "P_train", "Y_train",
          "X_basis", "P_basis", "Y_basis",
          "GPSDM", "fmingrad_Rprop"
        ),
        .packages = c("Matrix", "LaplacesDemon", "RhpcBLASctl")
      ) %dopar% {
        blas_set_num_threads(1)
        omp_set_num_threads(1)
        
        num_sample <- 50
        pars_sample <- pars_prev_updated[, i, drop = F] %*% matrix(1, nrow = 1, ncol = num_sample) + t(rmvn(num_sample, mu = rep(0, ndim + 3), Sigma = diag(1, nrow = ndim + 3, ncol = ndim + 3)))
        loglik_res <- rep(NA, num_sample)
        
        for (k in 1:num_sample) {
          lpost <- function(p, frac) {
            minibatch <- sample(1:nrow(X_train), round(nrow(X_train) * frac))
            return(GPSDM(pars = p, distMat = distMat, basisX = X_basis, basisP = P_basis, basisY = Y_basis, newX = X_train[minibatch, , drop = F], newP = P_train[minibatch, , drop = F], newY = Y_train[minibatch, , drop = F], priors_tr = list(E = pars_prev[, i, drop = F], V = pars_var_prev[[i]]), mode = c("optimize")))
          }
          
          res <- lpost(pars_sample[, k, drop = F], frac = 1)
          loglik_res[k] <- -res$neglpost
          print(k)
        }
        loglik_res <- loglik_res - mean(loglik_res)
        loglik_res[loglik_res < -50] <- -50
        loglik_res[loglik_res > 50] <- 50
        lik_res <- exp(loglik_res)
        pars_var_prev_updated_i <- cov.wt(t(pars_sample),
                                          wt = as.numeric(exp(loglik_res)),
                                          cor = FALSE, center = t(pars_prev_updated[, i, drop = F])
        )$cov
        pars_var_prev_updated_i
      }
  }
  
  ########################
  # Add new sets of hyperparameters
  ########################
  particles <- matrix(NA, nrow = ndim + 4, ncol = num_part)
  for (i in 1:ndim) {
    particles[i, ] <- rtrunc(n = num_part, spec = "norm", a = phimin, b = phimax, mean = priors$E_phi[i], sd = sqrt(priors$V_phi[i]))
    particles[i, ] <- -log((phimax - phimin) / (particles[i, ] - phimin) - 1) # transform
  }
  
  particles[ndim + 1, ] <- rtrunc(n = num_part, spec = "norm", a = vemin, b = vemax, mean = priors$E_ve, sd = sqrt(priors$V_ve))
  particles[ndim + 1, ] <- -log((vemax - vemin) / (particles[ndim + 1, ] - vemin) - 1) # transform
  
  particles[ndim + 2, ] <- rtrunc(n = num_part, spec = "norm", a = taumin, b = taumax, mean = priors$E_tau, sd = sqrt(priors$V_tau))
  particles[ndim + 2, ] <- -log((taumax - taumin) / (particles[ndim + 2, ] - taumin) - 1) # transform
  
  particles[ndim + 3, ] <- rtrunc(n = num_part, spec = "norm", a = rhomin, b = rhomax, mean = priors$E_rho, sd = sqrt(priors$V_rho))
  particles[ndim + 3, ] <- -log((rhomax - rhomin) / (particles[ndim + 3, ] - rhomin) - 1) # transform
  
  particles[ndim + 4, ] <- rtrunc(n = num_part, spec = "norm", a = gammamin, b = gammamax, mean = priors$E_gamma, sd = sqrt(priors$V_gamma))
  particles[ndim + 4, ] <- -log((gammamax - gammamin) / (particles[ndim + 4, ] - gammamin) - 1) # transform
  
  pars_add <- particles
  pars_id_add <- ((epoch - 1) * num_part + 1):(epoch * num_part)
  
  Rprop_res_add <-
    foreach(
      i = 1:num_part,
      # .export = ls(globalenv()),
      .export = c(
        "phimax", "phimin", "vemax", "vemin", "taumax", "taumin", "gammamax", "gammamin",
        "pars_add", "distMat", "priors",
        "X_train", "P_train", "D_train","Y_train",
        "X_basis", "P_basis", "D_basis","Y_basis",
        "GPSDM", "fmingrad_Rprop"
      ),
      .packages = c("Matrix", "LaplacesDemon", "RhpcBLASctl")
    ) %dopar% {
      blas_set_num_threads(1)
      omp_set_num_threads(1)
      
      # subset<-sample(1:nrow(X_new), nrow(X_new))
      # X_new<-X_new[subset,,drop=F]
      # P_new<-P_new[subset,,drop=F]
      # Y_new<-Y_new[subset,,drop=F]
      
      lpost <- function(p, sample_n, seed) {
        set.seed(seed)
        minibatch <- sample(1:nrow(X_train), sample_n)
        return(GPSDM(pars = p, distMat = distMat, basisX = X_basis, basisP = P_basis, basisD=D_basis,basisY = Y_basis, newX = X_train[minibatch, , drop = F], newP = P_train[minibatch, , drop = F],newD = D_train[minibatch, , drop = F], newY = Y_train[minibatch, , drop = F], priors = priors, mode = c("optimize")))
      }
      print("start")
      
      fmingrad_Rprop(lpost, pars_add[, i, drop = F] + matrix(rnorm(ndim + 4, 0, 0.1)), sample_n = min(nrow(X_train),50), maxcount = maxcount)
    }
  
  pars_add_updated <- matrix(NA, nrow = ndim + 4, ncol = num_part)
  log_p_add_updated <- rep(NA, num_part)
  for (i in 1:num_part) {
    pars_add_updated[, i] <- Rprop_res_add[[i]]$xopt
    log_p_add_updated[i] <- -Rprop_res_add[[i]]$fopt
  }
  
  pars_var_add_updated <-
    foreach(
      i = 1:num_part,
      .export = c( 
        "phimax", "phimin", "vemax", "vemin", "taumax", "taumin", "gammamax", "gammamin",
        "pars_add_updated", "distMat", "priors",
        "X_train", "P_train", "Y_train",
        "X_basis", "P_basis", "Y_basis",
        "GPSDM", "fmingrad_Rprop"
      ),
      .packages = c("Matrix", "LaplacesDemon", "RhpcBLASctl")
    ) %dopar% {
      blas_set_num_threads(1)
      omp_set_num_threads(1)
      
      num_sample <- 50
      pars_sample <- pars_add_updated[, i, drop = F] %*% matrix(1, nrow = 1, ncol = num_sample) + t(rmvn(num_sample, mu = rep(0, ndim + 4), Sigma = diag(1, nrow = ndim + 4, ncol = ndim + 4)))
      loglik_res <- rep(NA, num_sample)
      
      for (k in 1:num_sample) {
        lpost <- function(p, sample_n) {
          minibatch <- sample(1:nrow(X_train), sample_n)
          return(GPSDM(pars = p, distMat = distMat, basisX = X_basis, basisP = P_basis,basisD=D_basis, basisY = Y_basis, newX = X_train[minibatch, , drop = F], newP = P_train[minibatch, , drop = F],newD = D_train[minibatch, , drop = F], newY = Y_train[minibatch, , drop = F], priors = priors, mode = c("optimize")))
        }
        
        res <- lpost(pars_sample[, k, drop = F], sample_n = min(nrow(X_train),100))
        loglik_res[k] <- -res$neglpost
        print(k)
      }
      loglik_res <- loglik_res - mean(loglik_res)
      loglik_res[loglik_res < -50] <- -50
      loglik_res[loglik_res > 50] <- 50
      lik_res <- exp(loglik_res)
      pars_var_add_updated_i <- cov.wt(t(pars_sample),
                                       wt = as.numeric(exp(loglik_res)),
                                       cor = FALSE, center = t(pars_add_updated[, i, drop = F])
      )$cov
      pars_var_add_updated_i
    }
  
  pars_new <- cbind(pars_prev_updated, pars_add_updated)
  pars_var_new <- c(pars_var_prev_updated, pars_var_add_updated)
  log_p_new <- c(log_p_prev_updated, log_p_add_updated)
  pars_id_new <- c(pars_id_prev, pars_id_add)
  
  ########################
  # Choose sets of hyperparameters
  ########################
  # valid_id_subset<-vector(mode="list")
  # for (n in 1:10) {
  #   valid_id_subset[[n]]<-sample(1:nrow(X_valid),100)
  # }
  
  log_p_valid_new <-
    foreach(
      i = 1:ncol(pars_new),
      # .export = ls(globalenv()),
      .export = c(
        "phimax", "phimin", "vemax", "vemin", "taumax", "taumin", "gammamax", "gammamin",
        "pars_new", "distMat", "priors",
        "X_valid", "P_valid", "Y_valid",
        "X_basis", "P_basis", "Y_basis",
        "GPSDM", "fmingrad_Rprop"
      ),
      .combine = cbind
    ) %dopar% {
      blas_set_num_threads(1)
      omp_set_num_threads(1)


      num_sample <- 50
      loglik_res <- rep(NA, num_sample)

      for (k in 1:num_sample) {
        minibatch <- sample(1:nrow(X_valid), min(nrow(X_valid),100))
        res <- GPSDM(pars = pars_new[, i, drop = F], distMat = distMat, basisX = X_basis, basisP = P_basis, basisY = Y_basis, newX = X_valid[minibatch, , drop = F], newP = P_valid[minibatch, , drop = F], newY = Y_valid[minibatch, , drop = F],priors=priors, mode = c("optimize"))
        loglik_res[k] <- -res$neglpost
        print(k)
      }

      median(loglik_res)
    }
  
  sort_id <- rev(order(log_p_valid_new))[1:num_part]
  pars <- pars_new[, sort_id, drop = F]
  pars_var <- pars_var_new[sort_id]
  log_p <- log_p_new[sort_id]
  pars_id <- pars_id_new[sort_id]
  log_p_valid <- log_p_valid_new[sort_id]
  colnames(pars) <- pars_id
  
  write_csv(as.data.frame(pars), paste0(path_results,"/pars", ".csv"))
  for (i in 1:num_part) {
    write_csv(as.data.frame(pars_var[[i]]), paste0(path_results,"/pars_var_", i, ".csv"))
  }
  write_csv(as.data.frame(log_p), paste0(path_results,"/log_p_train", ".csv"))
  write_csv(as.data.frame(log_p_valid), paste0(path_results,"/log_p_valid", ".csv"))
  
  
  
  
  ######## Draw phi grid
  C <- pars[1:ndim, , drop = F]
  particle_wts <- exp(log_p_valid - mean(log_p_valid))
  
  C_df <- C %>%
    as_tibble() %>%
    mutate(embedding = colnames(read_csv(paste0(path_results, "/X_basis.csv")))) %>%
    mutate(variable = unlist(purrr::map(strsplit(embedding, "_"), 1))) %>%
    mutate(neighbor = unlist(purrr::map(strsplit(embedding, "_"), 2))) %>%
    mutate(lag = unlist(purrr::map(strsplit(embedding, "_"), 3))) %>%
    rowwise() %>% 
    mutate(lag=lags[[which(var_list==variable)]][[as.integer(lag)]]  %>% max()) %>% 
    gather(key = "combination", value = "phi", -embedding, -variable, -neighbor, -lag) %>%
    mutate(phi = (phimax - phimin) / (1 + exp(-phi)) + phimin) %>%
    mutate(combination = as.factor(combination)) %>%
    mutate(combination = forcats::fct_relevel(combination, rev(colnames(C)))) %>%
    # mutate(embedding = as.factor(embedding)) %>%
    mutate(variable = as.factor(variable)) %>%
    mutate(variable = forcats::fct_relevel(variable, as.character(var_list[vars]))) %>%
    mutate(neighbor = as.factor(neighbor)) %>%
    mutate(neighbor = forcats::fct_relevel(neighbor, as.character(seq(max(unlist(neighbors)))))) %>%
    mutate(lag = as.factor(lag))
  
  
  phi_weighted_mean<-C_df %>% 
    group_by(embedding,variable, neighbor, lag) %>% 
    summarize (phi=stats::weighted.mean(phi, particle_wts)) %>% 
    arrange(variable, neighbor, lag) 
  
  cairo_pdf(paste0(path_results,"/phi_grid.pdf"))
  p<-ggplot() +
    geom_raster(data=phi_weighted_mean %>% filter(phi>0),aes(x = lag, y = variable, fill = log(phi))) +
    geom_raster(data=phi_weighted_mean %>% filter(phi==0),aes(x = lag, y = variable), fill="black") +
    coord_equal()+
    theme_light(base_size = 15)+
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text.x  = element_text(size=8)
          # legend.title = element_blank()
    )+
    viridis::scale_fill_viridis()
  print(p)
  dev.off()
}

