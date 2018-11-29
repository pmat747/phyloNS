#' @title Nested Sampling in phylogenetics
#' @description It allows the use of Nested Sampling (NS) in phylogenetics for the fixed tree topology case.  NS is basically used to estimate the marginal likelihood, quantity used for model selection under a Bayesian framework.
#' @param data An alignment, object of class phyDat
#' @return list
#' @export ns

#library(phangorn)
#library(MCMCpack)
#library(phytools) # writenexus

##################################
##### log(x+y)=log(x)+log(y) #####
##################################

# log.plus <- function(x,y)
# {
#   if(x>y) x + log(1+exp(y-x))
#   else    y + log(1+exp(x-y))
# }

##############
### Priors ###
##############

### JC ###

log_prior_JC = function(vec, a , b, al, bl, n_brs){

  if( any(vec < 0) ){

    return( c(-Inf) )

  }

  # mu = vec[13];

  # branches + mean branch

  p1 = -(a + n_brs + 1) * log(vec[13]) - (sum(vec[-c(1:13)]) + b) / vec[13];

  return( p1 );

}

### HKY ###

log_prior_HKY = function(vec, a , b, al, bl, n_brs){

  if( any(vec < 0) ){

    return( c(-Inf) )

  }

  # mu = vec[13];

  # branches + mean branch

  p1 = -(a + n_brs + 1) * log(vec[13]) - (sum(vec[-c(1:13)]) + b) / vec[13];

  # Rates: q = vec[7] (HKY) = vec[6:10] (GTR); phi = vec[12];

  p  = log(vec[12]) - vec[12] * (vec[7] + 1);

  return( p1 + p );

}

### GTR ###

log_prior_GTR = function(vec, a , b, al, bl, n_brs){

  if( any(vec < 0) ){

    return( c(-Inf) )
  }

  # mu = vec[13];

  # branches + mean branch

  p1 = -(a + n_brs + 1) * log(vec[13]) - (sum(vec[-c(1:13)]) + b) / vec[13];

  # Rates: q = vec[7] (HKY) = vec[6:10] (GTR); phi = vec[12];

  p  = 5 * log(vec[12]) - vec[12] * ( sum(vec[6:10]) + 1 );

  return( p1 + p );

}

### JC_G ###

log_prior_JC_G = function(vec, a , b, al, bl, n_brs){

  if( any(vec < 0) ){

    return( c(-Inf) )
  }

  # mu = vec[13];

  # branches + mean branch

  p1 = -(a + n_brs + 1) * log(vec[13]) - (sum(vec[-c(1:13)]) + b) / vec[13];

  # Rates: q = vec[7] (HKY) = vec[6:10] (GTR); phi = vec[12];

  # Gamma parameter: lambda = vec[1];

  p2     = (al - 1.0) * log(vec[1]) - vec[1] / bl;

  return( p1 + p2);

}

### HKY+G ###

log_prior_HKY_G = function(vec, a , b, al, bl, n_brs){

  if( any(vec < 0) ){

    return( c(-Inf) )
  }

  # mu = vec[13];

  # branches + mean branch

  p1 = -(a + n_brs + 1) * log(vec[13]) - (sum(vec[-c(1:13)]) + b) / vec[13];

  # Rates: q = vec[7] (HKY) = vec[6:10] (GTR); phi = vec[12];

  # Gamma parameter: lambda = vec[1];

  p2     = (al - 1.0) * log(vec[1]) - vec[1] / bl;

  p  = log(vec[12]) - vec[12] * (vec[7] + 1);

  return( p1 + p2 + p );

}

### GTR+G ###

log_prior_GTR_G = function(vec, a , b, al, bl, n_brs){

  if( any(vec < 0) ){

    return( c(-Inf) )
  }

  # mu = vec[13];

  # branches + mean branch

  p1 = -(a + n_brs + 1) * log(vec[13]) - (sum(vec[-c(1:13)]) + b) / vec[13];

  # Rates: q = vec[7] (HKY) = vec[6:10] (GTR); phi = vec[12];

  # Gamma parameter: lambda = vec[1];

  p2     = (al - 1.0) * log(vec[1]) - vec[1] / bl;

  p  = 5 * log(vec[12]) - vec[12] * ( sum(vec[6:10]) + 1 );

  return( p1 + p2 + p );

}

##############################
### Positions for sampling ###
##############################

# It yields the positions and its weigths
#   for sampling according to the model.

pos_sam = function(model, n_brs){

  if(model == "JC"){

    pos  = c(13, 14:(13 + n_brs));
    prob = rep(1, 1 + n_brs);

  }else if(model == "JC_G"){

    pos  = c(1, 13, 14:(13 + n_brs));
    prob = rep(1, 2 + n_brs);

  }else if(model == "HKY"){

    pos  = c(2:5, 7, 12, 13, 14:(13 + n_brs));
    prob = rep(1, 7 + n_brs);

  }else if(model == "HKY_G"){

    pos  = c(1, 2:5, 7, 12, 13, 14:(13 + n_brs));
    prob = rep(1, 8 + n_brs);

  }else if(model == "GTR"){

    pos  = c(2:5, 6:10, 12, 13, 14:(13 + n_brs));
    prob = rep(1, 11 + n_brs);

  }else if(model == "GTR_G"){

    pos = c(1, 2:5, 6:10, 12, 13, 14:(13 + n_brs));
    prob = rep(1,12 + n_brs);

  }

  return( list(pos = pos, prob = prob) );

}

##########################
### Generator function ###
##########################

proposal <- function(vec, pos, model, prior_width){

  # abs acts as reflection for negative values
  out = abs(vec[pos] + prior_width * 10.0^(1.5 -
                                             3.0 * abs(stats::rt(1, df = 2))) * stats::rnorm(1));

  vec[pos] = out;

  if(pos >= 12){

    NULL;

  }else if( (pos >= 2) && ( pos <= 5) ){ # frequencies

    vec[2:5] = vec[2:5] / sum(vec[2:5]);

  }else if( (pos >=6) && ( pos <= 11)){

    # rates

    if((model == "HKY") || (model == "HKY_G")){

      vec[10] = out;

    }

  }else if( pos == 1 ){

    while( out < 0.10 ){ # Gamma parameter (0.14)

      out = abs(vec[pos] + prior_width * 10.0^(1.5 -
            3.0 * abs(stats::rt(1, df = 2))) * stats::rnorm(1));
    }

    vec[pos] = out;

  }

  return(vec);

}

genValue <- function(data, tree, vec0, l0, worstLike,
                     a, b, al, bl, k, step, spe){

  # Random numbers

  u = stats::runif(step);

  pos = sample(x = spe$pos, size = step, replace = TRUE, prob = spe$prob);

  # Starting Value

  p0  = spe$log_prior(vec0, a, b, al, bl, spe$n_brs); # prior

  # loop

  i = 0;

  for(ipos in pos){

    if(ipos > 13){

      prior_width = (log(4) + 1.5*log(3)) * vec0[13]; # Tukey criteria (rate = 1/vec0[13])

    }else if(ipos == 13){ prior_width = spe$mu_w;

    }else if(ipos == 1){ prior_width = spe$lambda_w ;

    }else if((ipos >= 2) && (ipos <= 5)){ prior_width = 0.77; # 4*sqrt(3/(16*5))

    }else if((ipos >= 6) && (ipos <= 11)){

      prior_width = (log(4) + 1.5*log(3)) / vec0[12];

    }else if(ipos == 12){ prior_width = 4;} # phi~Exp(1)=>dexp(0.99)-dexp(.01)>6

    # proposal value

    vec1 = proposal(vec0, ipos, spe$model, prior_width);

    if( (ipos == 12) || (ipos == 13) ){

      l1 = l0; # phi & mu do not alter the likelihood

    }else{

      if( ipos > 13 ){

        tree$edge.length[ipos - 13] = vec1[ipos]; # branches

      }

      l1 = phangorn::pml(tree, data, Q = vec1[6:11], bf = vec1[2:5],
                         shape = vec1[1], k = k)$logLik;
    }

    p1   = spe$log_prior(vec1, a, b, al, bl, spe$n_brs);

    i    = i + 1;

    if( ( u[i] < exp(p1 - p0) ) &&
        ( worstLike < l1 ) ){

      # this avoids copying too many unneccesary elements

      if( ipos >= 14 ){ # branch

        vec0[ipos] = vec1[ipos];

      }else if( (ipos >= 2) && (ipos <= 5) ){ # frequency

        vec0[2:5] = vec1[2:5];

      }else if((ipos == 1) || (ipos == 12) ||
               (ipos == 13) ){ # Gamma par, mu & phi

        vec0[ipos] = vec1[ipos];

      }else if( (ipos >= 6) && (ipos <= 11)){ # Rates for HKY & GTR versions

        if( (spe$model == "HKY") || (spe$model == "HKY_G") ){

          vec0[c(7, 10)] = vec1[ipos];

        }else{ # GTR model

          vec0[ipos] = vec1[ipos];

        }
      }

      l0   = l1;
      p0   = p1;

    }else{

      if( ipos > 13 ){ # reseting proposal value in "tree" branch length

        tree$edge.length[ipos - 13] = vec0[ipos];

      }
    }

  } # ending step

  return( list( theta = vec0, like = l0 ) );

}

#######################
### Nested Sampling ###
#######################

# "vns" allows variable topology

##############
### Priors ###
##############

# Vector of prior parameters
# theta = c(l[1], freq[2:5], q[6:11], phi[12], mu[13], t[14:])

# Branch length - t_i|mu~Exp(1/mu)
#   f(t_i|mu) = (1/mu) * exp(-t_i/mu)

# Branch length mean - mu~Inverse Gamma
#   f(mu) = (b^a/G(a)) * mu^(-a-1) * exp(-b/mu)
#   a = shape
#   b = scale

# Transition rates - q_i|phi~Exp(phi)
#   f(q_i|phi) = phi * exp(-q_i*phi)

# Rate mean - phi~Exp(1)
#   f(phi) = exp(-phi)

# Gamma parameter l~Gamma(a, b)
#   f(l) = 1/(G(a)*b^a) * l^(a-1) * exp(-l/b)
#   al = shape
#   bl = scale

# Frequencies (Dirichlet(1,1,1,1))

ns = function( data, tree, N = 5, k = 4,
               a = 3, b = 0.2, al = 1, bl = 1,
               step = 50, max.iter = 2000, end = 2,
               act_plot = TRUE, model){

  # checking if the model is defined

  models = c("JC", "JC_G", "HKY", "HKY_G", "GTR", "GTR_G");

  if( !any(model == models) )stop("The model is not defined");

  if( any(c(a,b,al,bl) < 0) ) stop("Prior parameters must be positive");

  if( length(tree$tip.label) != length(data) ) stop("Data and tree do not match");

  t0     = proc.time(); # starting time

  ### Priors definitions ###

  if(model == "JC"){          log_prior = log_prior_JC;

  }else if(model == "HKY"){   log_prior = log_prior_HKY;

  }else if(model == "GTR"){   log_prior = log_prior_GTR;

  }else if(model == "JC_G"){  log_prior = log_prior_JC_G;

  }else if(model == "HKY_G"){ log_prior = log_prior_HKY_G;

  }else if(model == "GTR_G"){ log_prior = log_prior_GTR_G;

  }

  ###

  n_taxa = length(data); # number of taxa/leaves

  if(ape::is.rooted(tree) == TRUE){ # number of branches

    n_brs = 2 * n_taxa - 2; # rooted case

  }else{

    n_brs = 2 * n_taxa - 3; # unrooted case

  }

  spe   = list();

  ### Sampling the prior according to the model ###

  # Branch lengths

  mu = as.matrix(MCMCpack::rinvgamma(N, shape = a, scale = b), ncol = 1);

  t  = t(apply(mu, 1, function(x)stats::rexp(n_brs, rate = 1 / x) ));

  m  = unlist(strsplit(model, ""));

  ### Gamma parameter ###

  if(m[length(m)] != "G"){

    lambda = rep(1, N); # JC - HKY

    k = 1;

  }else{

    if( is.null(al) || is.null(bl) )
      stop('Provide full specified prior for the gamma shape parameter');

    lambda_w = stats::qgamma(0.999, shape = al, scale = bl, lower.tail = TRUE) -
               stats::qgamma(0.001, shape = al, scale = bl, lower.tail = TRUE);

    spe$lambda_w = lambda_w; # prior width

    lambda = stats::rgamma(N, shape = al, scale = bl); # JC_G HKY_G GTR_G

    for(i in 1:N){ # It prevents 'phangorn' from crushing

      while( lambda[i] < 0.07 ){ # It could even crush with 0.07464181

        lambda[i] = stats::rgamma(1, shape = al, scale = bl);

      }
    }
  }

  ### Frequency ###

  if( m[1] != "J"){

    freq  = matrix( stats::rexp( N * 4, 1), ncol = 4);

    freq  = freq / apply(freq, 1, sum); # Dirichlet(1,1,1,1)

  }else{

    freq   = matrix( 0.25, nrow = N, ncol = 4); # JC JC_G

  }

  ### Rates ###

  if( m[1] == "J" ){

    q   = matrix(1, nrow = N, ncol = 6);
    phi = rep(1, N);

  }else if( m[1] == "H" ){

    phi =  as.matrix(stats::rexp(N), ncol = 1);

    q   = c(apply(phi, 1, function(x)stats::rexp(n = 1, rate = x)));

    q   = cbind(rep(1, N), q, rep(1, N), rep(1, N), q, rep(1, N));

  }else if( m[1] == "G" ){

    phi =  as.matrix(stats::rexp(N), ncol = 1);

    q   = t(apply(phi, 1, function(x)stats::rexp(n = 5, rate = x)));

    q   = cbind(q, rep(1,N));

  }

  ### end sampling from the prior ###

  # specifications

  v   = pos_sam(model, n_brs);

  # prior width

  mu_w = MCMCpack::rinvgamma(2000, shape = a, scale = b)

  mu_w = max(mu_w) - min(mu_w); # mu

  spe$n_taxa = n_taxa;
  spe$n_brs  = n_brs;
  spe$model  = model;
  spe$pos    = v$pos;
  spe$prob   = v$prob;
  spe$log_prior = log_prior;
  spe$mu_w = mu_w;

  ### Objects: parameters & trees ###

  obj_par = cbind(lambda, freq, q, phi, mu, t); # parameter

  logL    = apply(obj_par, 1, function(x){
                                     tree$edge.length = x[-c(1:13)];
                                     phangorn::pml(tree, data, Q = x[6:11],
                                     bf = x[2:5], shape = x[1], k = k)$logLik});
  # definitions

  iter      = 0;
  logwidth1 = log( 1.0 - exp( - 1.0 / N ) ); # simple
  #logwidth1 = log(1 - exp( - 2 / N )) - log(2); # Trapezoidal
  H         = 0.0;
  dTheta    = NULL; # discarded values
  lw        = NULL; # like * prior
  logX      = NULL; # prior mass
  logLd     = NULL; # discarded log-likelihoods
  s_z       = NULL; # sequence of z-evolution
  s_z[1]    = -1.79769e+308;

  repeat{

    iter        <- iter + 1;

    worst       <- which(logL == min(logL))[1];

    logX[iter]  <- - iter / N;

    logLd[iter] <- logL[worst];

    logwidth    <- logwidth1 - (iter - 1) / N;

    lw[iter]    <- logwidth + logL[worst];

    logZ        <- log.plus(s_z[iter], lw[iter]);

    H           <- exp(lw[iter]-logZ)*logL[worst] - logZ +
                         exp(s_z[iter]-logZ)*(H+s_z[iter]);

    s_z[iter+1] <- logZ; # Storing z history

    cat(c(iter, logZ), sep = "  ", "\n"); # printing iter & logZ

    # Storing discarded point

    dTheta = rbind(dTheta, obj_par[worst,] ); # parameters

    if( (iter >= end * N * H) || (iter >= max.iter) ) break;

      # starting value for new active point #

      #init = sample(c(1:N)[-worst], size = 1);
      init = ifelse(N == 1, 1, sample(c(1:N)[-worst], size = 1));

      tree$edge.length = obj_par[init, -c(1:13)];

      values = genValue(data, tree, vec0 = obj_par[init,],
                        l0 = logL[init], worstLike = logL[worst],
                        a, b, al, bl, k, step, spe);

      obj_par[worst, ]  = values$theta;

      logL[worst]       = values$like;
  }

  logZc   <- log.plus(logZ, log(mean(exp(logL))) - iter/N )

  Lw.z    <- exp(lw - logZc)

  max.post<- exp(- sum(Lw.z*log(Lw.z), na.rm = TRUE)) # number of posterior samples

  Lw.z    <- Lw.z / max(Lw.z)

  d       <- sample(iter, replace = TRUE, size = max.post, prob = Lw.z )

  ### Posterior sampling ###

  sampled_par <- matrix(NA, nrow = max.post, ncol = dim(obj_par)[2]);

  for( i in 1:max.post){

    sampled_par[i, ] = dTheta[d[i],];

  }

  # Naming posterior samples & discarded points
  colnames(sampled_par) = c('lambda', 'freqA', 'freqC', 'freqG', 'freqT',
                            'qCA', 'qGA', 'qTA', 'qGC', 'qTC', 'qTG',
                            'phi', 'mu',
                            paste(rep('t', n_brs), 1:n_brs, sep = ""));

  colnames(dTheta) = c('lambda', 'freqA', 'freqC', 'freqG', 'freqT',
                       'qCA', 'qGA', 'qTA', 'qGC', 'qTC', 'qTG',
                       'phi', 'mu',
                        paste(rep('t', n_brs), 1:n_brs, sep = ""));

  ###

  t1   = proc.time(); # ending time

  time = (t1 - t0)[1] # calculating time

  ### Plots ###

  if(act_plot){

    graphics::par(mar = c(4,4,4,4));

    graphics::plot(logX, logLd, xlab = expression("log"* xi), ylab='log(L)', pch='.');

    graphics::plot(logX, Lw.z, xlab = expression("log"* xi), ylab = 'Weight for posterior');

    meanBranch = apply(sampled_par[, - c(1:13)], 2, mean);

    tree$edge.length = meanBranch;

    phytools::plotTree(ape::ladderize(tree)); # phytool package ape respectively
    ape::add.scale.bar(cex = 0.7, font = 1) # ape package

  }

  # Storing process info
  info = list(N = N, k = k, a = a, b = b, al = al, bl = bl,
              step = step, model = model, tree = tree, ns = "fixed_tree");

  return(list(logZ = logZ, logZc = logZc,
              sd_logZ = sqrt(H / N),
              H = H, iter = iter, time = time,
              logX = logX, logLd = logLd,
              Lw.z = Lw.z, seq_lz = s_z[-1],
              dTheta = dTheta, sampled_par = sampled_par,
              info = info))
}

###########
### Run ###
###########

# N = 15; k = 4; a = 3; b = 0.2; al = 10; bl = 0.026; step = 100;
# model ='HKY_G'; max.iter = 100; end = 2; rooted = FALSE

#R = ns(rbcL, tree, N = 100, k = 4, a = 3, b = 0.2, al = 10, bl = 0.026,
#               step = 50, max.iter = 2000, end = 2, act_plot = TRUE,
#               model = "GTR_G")

# Output
# (lambda, freq, q, phi, mu, t); # parameters
