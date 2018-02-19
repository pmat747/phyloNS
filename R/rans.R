#' @title Re-analysis of NS output(s)
#' @description This ...
#' @param data An alignment, object of class phyDat
#' @return list
#' @export rans


rans = function(R, S = 100, n_ns_outputs = 1,
                act_plot = TRUE){

  if(n_ns_outputs == 1){ # case of one NS run

    logL  = R$logLd;
    n     = length(logL);
    orLs  = order(logL, decreasing = FALSE); # positions
    Theta = R$dTheta;
    Theta = Theta[orLs, ]
    logL  = sort(logL); # sorted likelihoods
    N     = R$info$N;

  }else{ # case of multiple NS runs

    n_ns_outputs = length(objects(R));
    logL  = NULL;
    Theta = NULL;
    N     = 0;

    for(i in 1:n_ns_outputs){

      logL  = c(logL, R[[i]]$logLd);
      Theta = rbind(Theta, R[[i]]$dTheta);
      N     = N + R[[i]]$info$N;

    }
    n     = length(logL);
    orLs  = order(logL, decreasing = FALSE); # positions
    Theta = Theta[orLs, ]
    logL  = sort(logL); # sorted likelihoods

  }

  sim_mu = NULL;
  sim_sd = NULL;
  sim_med= NULL;
  vecZ   = NULL;
  vecH   = NULL;

  ### Trapezoidal rule ###

  for(j in 1:S){

    print(j);

    us        = rbeta(n, shape1 = N, shape2 = 1);

    logX      = cumsum(log(us));

    aux       = us * c(us[-1], rbeta(1,shape1 = N, shape2 = 1));

    lw        = c(0, logX[-n]) + log(1 - aux) - log(2); # Trapezoidal rule

    lw_lL     = lw + logL;

    logZ      = logplusvec(lw_lL);

    H         = sum(exp(lw_lL - logZ) * logL) - logZ;

    sd        = sqrt(H/N);

    ### Sampling ###

    logPostW = lw_lL - logZ; # log-posterior weigths

    postW    = exp(logPostW); # posterior weigths

    M        = round( exp(-sum(postW * logPostW)) );

    sam      = sample(1:n, size = M, replace = TRUE, prob = postW / max(postW));
    #sam      = sample.lp(1:n, size = M, logpbb = logPostW);

    sim      = as.matrix(Theta[sam, ]);

    sim_mu  = rbind(sim_mu, apply(sim, 2, mean));

    sim_sd  = rbind(sim_sd, apply(sim, 2, sd));

    sim_med = rbind(sim_med, apply(sim, 2, median));

    vecZ[j] = logZ;
    vecH[j] = H;

  }

  rownames(sim_mu)  = NULL;
  rownames(sim_sd)  = NULL;
  rownames(sim_med) = NULL;

  ### mean ###
  aux1 = apply(sim_mu,2, mean);
  aux2 = apply(sim_mu,2, sd);
  sim_mu_est  = cbind(aux1, aux2, aux1 - 1.96 * aux2,  aux1 + 1.96 * aux2);
  sim_mu_est  = round(sim_mu_est,4);
  colnames(sim_mu_est) = c('mu', 'error', 'lower', 'upper');

  ### standard deviation ###
  aux1 = apply(sim_sd,2, mean);
  aux2 = apply(sim_sd,2, sd);
  sim_sd_est  = cbind(aux1, aux2, aux1 - 1.96 * aux2,  aux1 + 1.96 * aux2);
  sim_sd_est  = round(sim_sd_est,4);
  colnames(sim_sd_est) = c('sd', 'error', 'lower', 'upper');

  ### median ###
  aux1 = apply(sim_med,2, mean);
  aux2 = apply(sim_med,2, sd);
  sim_med_est  = cbind(aux1, aux2, aux1 - 1.96 * aux2,  aux1 + 1.96 * aux2);
  sim_med_est  = round(sim_med_est, 4);
  colnames(sim_med_est) = c('median', 'error', 'lower', 'upper');

  vecZ_est = c(mean(vecZ), mean(vecZ) - 1.96 * sd(vecZ),
               mean(vecZ) + 1.96 * sd(vecZ));
  vecZ_est = matrix(vecZ_est, nrow = 1);
  colnames(vecZ_est) = c('logZ', 'lower', 'upper');

  vecH_est = c(mean(vecH), mean(vecH) - 1.96 * sd(vecH),
               mean(vecH) + 1.96 * sd(vecH));
  vecH_est = matrix(vecH_est, nrow = 1);
  colnames(vecH_est) = c('H', 'lower', 'upper');

  if(act_plot == TRUE){
    par(mar = c(4,4,4,4))
    hist(vecZ,
         main = paste("Marginal likelihood distribution for N=", N, sep = ""))

  }

  return( list(sim_mu_est = sim_mu_est, sim_sd_est = sim_sd_est,
               sim_med_est = sim_med_est, vecZ_est = vecZ_est,
               vecH_est = vecH_est, vecZ = vecZ, vecH = vecH,
               sim_mu = sim_mu, sim_sd = sim_sd));

}
