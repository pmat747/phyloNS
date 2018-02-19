#' @title Nested importance Sampling in phylogenetics
#' @description It allows the use of Nested Sampling (NS) in phylogenetics for the variable tree topology case.  NS is basically used to estimate the marginal likelihood, quantity used for model selection under a Bayesian framework.
#' @param data An alignment, object of class phyDat
#' @return list
#' @export nis

#source("prior&ref(fixed_top).R")

#########################
### Pseudo likelihood ###
#########################

dpselike = function(data, tree, vec, ref_par, pr_par,
                    k, dprior, dpseprior, rec_val){

  if( any(vec < 0) ){

    return( c(-Inf) );

  }

  if(is.null(rec_val)){

    #tree$edge.length = vec[-c(1:13)];

    like = pml(tree, data, Q = vec[6:11], bf = vec[2:5],
               shape = vec[1], k = k)$logLik;

  }else{

    like = rec_val$like; # not affected (when mu & phi are updated)

  }

  prior    = dprior(vec, pr_par);

  pseprior = dpseprior(vec, ref_par);

  return( list(pselike = like + prior - pseprior,
               like = like, pseprior = pseprior));
}

##############################
### Positions for sampling ###
##############################

# It yields the positions and its weigths
#   for sampling according to the model.

# 'pos_sam' uses this function from ns.R

##########################
### Generator function ###
##########################

# it generates one value

# proposal uses this function from ns.R

### generator

ref_gen = function(data, tree, vec0, inf0, worstLike,
                   step, spe, ref_par, pr_par, k){

  # Random numbers

  u = runif(step);

  pos = sample(x = spe$pos, size = step, replace = TRUE, prob = spe$prob);

  # loop

  i = 0;

  for(ipos in pos){

    if(ipos > 13){ # branches

      inf1 = NULL;
      prior_width = spe$t_w[ipos - 13];

    }else if(ipos == 13){ # mu

      inf1 = inf0;
      prior_width = spe$mu_w;

    }else if(ipos == 1){# gamma parameter

      inf1 = NULL;
      prior_width = spe$lambda_w ;

    }else if((ipos >= 2) && (ipos <= 5)){ # frequencies

      inf1 = NULL;
      prior_width = spe$f_w[ipos - 1];

    }else if((ipos >= 6) && (ipos <= 11)){ # rates

      inf1 = NULL;
      prior_width = spe$q_w[ipos - 5];

    }else if(ipos == 12){ # phi

      inf1 = inf0;
      prior_width = spe$phi_w;

    }

    vec1 = proposal(vec0, ipos, spe$model, prior_width);

    if( ipos > 13 ){

      tree$edge.length[ipos - 13] = vec1[ipos];

    }

    inf1 = dpselike(data, tree, vec1, ref_par, pr_par,
                    k, spe$dprior, spe$dpseprior, inf1);

    i   = i + 1;

    if( ( u[i] < exp(inf1$pseprior - inf0$pseprior) ) &&
        ( worstLike < inf1$pselike ) ){

      # this avoids copying too many unneccesary elements

      if( (ipos >= 12) || (ipos == 1) ){ # lambda + phi + mu + t

        vec0[ipos] = vec1[ipos];

      }else if( (ipos >= 2) && (ipos <= 5) ){ # frequency

        vec0[2:5] = vec1[2:5];

      }else if( (ipos >= 6) && (ipos <= 11)){ # Rates for HKY & GTR versions

        if( (spe$model == "HKY") || (spe$model == "HKY_G") ){

          vec0[c(7, 10)] = vec1[ipos];

        }else{ # GTR model

          vec0[ipos] = vec1[ipos];

        }
      }

      inf0 = inf1;

    }else{

      if( ipos > 13 ){ # reseting proposal value in "tree" branch length

        tree$edge.length[ipos - 13] = vec0[ipos];

      }
    }

  } # ending step

  return( list( theta = vec0, pselike = inf0$pselike,
                pseprior = inf0$pseprior, like = inf0$like) );
}

##################################
### Nested Importance Sampling ###
##################################

nis = function( data, tree, N = 1, k = 4,
                a = 3, b = 0.2, al = NULL, bl = NULL,
                step = 50, max.iter = 2000, tol = 1e-07,
                model, post_sample = NULL){

  #############
  ### Check ###
  #############

  # checking if the model is defined

  models = c("JC", "JC_G", "HKY", "HKY_G", "GTR", "GTR_G");

  if( !any(model == models) )stop("The model is not defined");

  if( any(c(a,b,al,bl) < 0) )stop("Prior parameters must be positive");

  if( length(tree$tip.label) != length(data) ) stop("Data and tree do not match");

  t0     = proc.time(); # starting time

  #########################
  ### Prior & reference ###
  ###   distributions   ###
  #########################

  if(model == "JC"){          dprior = dprior_JC;

                              dpseprior = dpseprior_JC;

  }else if(model == "HKY"){   dprior = dprior_HKY;

                              dpseprior = dpseprior_HKY;

  }else if(model == "GTR"){   dprior = dprior_GTR;

                              dpseprior = dpseprior_GTR;

  }else if( model == "JC_G"){ dprior = dprior_JC_G;

                              dpseprior = dpseprior_JC_G;

  }else if(model == "HKY_G"){ dprior = dprior_HKY_G;

                              dpseprior = dpseprior_HKY_G;

  }else if( model == "GTR_G"){dprior = dprior_GTR_G;

                              dpseprior = dpseprior_GTR_G;
  }

  # definitions

  spe    = list(); # to store prior widths + model specifications + ...

  pr_par = list(a = a, b = b); # compulsory

  logL = NULL; # pseudo likelihood

  logP = NULL; # pseudo prior

  like = NULL; # likelihood

  n_taxa = length(data); # number of taxa/leaves

  n_brs = dim(tree$edge)[1]; # number of branches

  obj_par = matrix(NA, nrow = N, ncol = (13 + n_brs));

  #########################
  ### Posterior samples ###
  #########################

  if( is.null(post_sample) ){

    stop('Sample the posterior first!')

    # sample posterior

  }else{

    post_mean = apply(post_sample, 2, mean);

    post_sd   = apply(post_sample, 2, sd);

  }

  ############################
  ### Reference parameters ###
  ############################

  m  = unlist(strsplit(model, ""));

  ##########
  ### mu ###
  ##########

  a_r  = (post_mean[13] / post_sd[13] )^2 + 2;
  b_r  = post_mean[13] * ((post_mean[13] / post_sd[13] )^2 + 1);

  obj_par[, 13] = rinvgamma(N, shape = a_r, scale = b_r); # active points

  mu_w = rinvgamma(2000, shape = a_r, scale = b_r);

  mu_w = max(mu_w) - min(mu_w); # sampling prior width

  ################
  ### Branches ###
  ################

  at_r = (post_mean[-c(1:13)] / post_sd[-c(1:13)])^2;
  bt_r = (post_sd[-c(1:13)])^2 / post_mean[-c(1:13)];

  t = NULL;

  for(i in 1:n_brs){

    t = cbind(t, rgamma(N, shape = at_r[i], scale = bt_r[i]));

  }

  obj_par[, 14:(13+n_brs)] = t; # active points

  t_w = qgamma(0.999, shape = at_r, scale = bt_r, lower.tail = TRUE) -

        qgamma(0.001, shape = at_r, scale = bt_r, lower.tail = TRUE); # width

  ref_par = list(a_r = a_r, b_r = b_r, at_r = at_r, bt_r = bt_r); # ref par

  #######################
  ### Gamma parameter ###
  #######################

  if(m[length(m)] == "G"){

    if( is.null(al) || is.null(bl) )
      stop('Provide full specified prior for the gamma parameter');

    pr_par$al = al; # prior parameters
    pr_par$bl = bl; # prior parameters

    al_r = (post_mean[1] / post_sd[1])^2;
    bl_r = (post_sd[1])^2 / post_mean[1];

    ref_par$al_r = al_r; # reference paramaters
    ref_par$bl_r = bl_r; # reference paramaters

    lambda_w = qgamma(0.999, shape = al_r, scale = bl_r, lower.tail = TRUE) -
               qgamma(0.001, shape = al_r, scale = bl_r, lower.tail = TRUE);

    spe$lambda_w = lambda_w; # prior sampling width

    obj_par[, 1] = rgamma(N, shape = al_r, scale = bl_r); # active points

    # preventing from crashing

    for(i in 1:N){ # It prevents phangorn from crashing

      while( obj_par[i, 1] < 0.07 ){ # It could even crush with 0.07464181

        obj_par[i, 1] = rgamma(1, shape = al_r, scale = bl_r);

      }
    }

  }else{

    obj_par[, 1] = rep(1, N);

    k = 1;

  }

  #########################
  ### Frequency & rates ###
  #########################

  if( m[1] == "J"){

    obj_par[, 2:5] = 0.25; # frequencies

    obj_par[, 6:12] = 1.0; # rates and mean rate

  }else{

    ### Frequency ###

    mf  = sum(post_mean[2:5]^2 * (1.0 - post_mean[2:5])^2) /
          sum(post_sd[2:5]^2 * post_mean[2:5] * (1.0 - post_mean[2:5])) - 1.0;

    f_r = post_mean[2:5] * mf;

    ref_par$f_r = f_r; # reference parameter

    f_w = rdirichlet(1000, f_r);

    f_w = apply(f_w, 2, max) - apply(f_w, 2, min);

    spe$f_w = f_w; # width

    obj_par[, 2:5] = rdirichlet(N, f_r); # active points

    ### Rates ###

    aphi_r = (post_mean[12] / post_sd[12])^2;
    bphi_r = (post_sd[12])^2 / post_mean[12];

    ref_par$aphi_r = aphi_r; # reference parameter
    ref_par$bphi_r = bphi_r; # reference parameter

    phi_w = qgamma(0.999, shape = aphi_r, scale = bphi_r, lower.tail = TRUE) -
            qgamma(0.001, shape = aphi_r, scale = bphi_r, lower.tail = TRUE);

    spe$phi_w = phi_w; # width

    obj_par[, 12] = rgamma(N, shape = aphi_r, scale = bphi_r); # active point

  if( m[1] == "H"){

    aq_r = (post_mean[7] / post_sd[7])^2;
    bq_r = (post_sd[7])^2 / post_mean[7];

    q_w = qgamma(0.999, shape = aq_r, scale = bq_r, lower.tail = TRUE) -
          qgamma(0.001, shape = aq_r, scale = bq_r, lower.tail = TRUE);

    spe$q_w = c(NA, q_w); # width  (match position with HKY+G)

    q = rgamma(N, shape = aq_r, scale = bq_r);

    obj_par[,6:11] = cbind(rep(1,N), q, rep(1,N), rep(1,N), q, rep(1,N));#active points

  }else if( m[1] == "G"){

    aq_r = (post_mean[6:10] / post_sd[6:10])^2;
    bq_r = (post_sd[6:10])^2 / post_mean[6:10];

    q_w = qgamma(0.999, shape = aq_r, scale = bq_r, lower.tail = TRUE) -
          qgamma(0.001, shape = aq_r, scale = bq_r, lower.tail = TRUE);

    spe$q_w = q_w; # width

    q = NULL

    for(i in 1:5){

      q = cbind(q, rgamma(N, shape = aq_r[i], scale = bq_r[i]));

    }

    obj_par[, 6:11] = cbind(q, rep(1,N)); # active points

  }

  ref_par$aq_r = aq_r; # reference parameter
  ref_par$bq_r = bq_r; # reference parameter

  }

  v = pos_sam(model, n_brs); # position to sample according to the model

  # storing specifications in "spe"

  spe$dprior    = dprior;    # prior function
  spe$dpseprior = dpseprior; # pseudo-prior function
  spe$n_taxa = n_taxa;
  spe$n_brs = n_brs;
  spe$model  = model;
  spe$pos    = v$pos;  # positions to sample
  spe$prob   = v$prob; # ppb to sample each parameter
  spe$mu_w   = mu_w;
  spe$t_w    = t_w;

  ### likelihood of active points ###

  for(i in 1:N){

    tree$edge.length = obj_par[i, -c(1:13)];

    v = dpselike(data, tree, obj_par[i,], ref_par, pr_par, k,
                 dprior, dpseprior, rec_val = NULL);

    logL[i] = v$pselike;
    logP[i] = v$pseprior;
    like[i] = v$like;

  }

  # definitions

  iter      = 0;
  logwidth1 = log( 1.0 - exp( - 1.0 / N ) );
  H         = 0;
  dTheta    = NULL; # discarded values
  lw        = NULL; # like * prior
  logX      = NULL; # prior mass
  logPsLd   = NULL; # discarded log-pseudo-likelihoods
  logPsPd     = NULL; # discarded log-pseudo-priors
  logLd     = NULL; # discarded log-likelihoods
  s_z       = NULL; # sequence of z-evolution
  s_z[1]    = -1.79769e+308;

  print("Nested Sampling loop")

  repeat{

    iter          <- iter + 1;

    worst         <- which(logL == min(logL))[1];

    logX[iter]    <- - iter / N;

    logPsLd[iter] <- logL[worst];
    logPsPd[iter] <- logP[worst];
    logLd[iter]   <- like[worst];

    logwidth      <- logwidth1 - (iter - 1) / N;

    lw[iter]      <- logwidth + logL[worst];

    logZ          <- log.plus(s_z[iter], lw[iter]);

    H             <- exp(lw[iter]-logZ)*logL[worst] - logZ +
                     exp(s_z[iter]-logZ)*(H+s_z[iter]);

    s_z[iter+1] <- logZ; # Storing z history

    cat(c(iter, logZ), sep = "  ", "\n"); # printing iter & logZ

    # Storing discarded point

    dTheta = rbind(dTheta, obj_par[worst,] ); # parameters

    ###

    if( (abs((logZ - s_z[iter])/logZ) <= tol) || (iter >= max.iter) ) break;

    # starting value to generate a new one

    #init      = sample(c(1:N)[-worst], size = 1);
    init = ifelse(N == 1, 1, sample(c(1:N)[-worst], size = 1));

    vec0      = obj_par[init,];

    inf0      = list(pselike = logL[init], pseprior = logP[init],
                     like = like[init]);

    tree$edge.length = obj_par[init, -c(1:13)];

    newvalue  = ref_gen(data, tree, vec0, inf0, worstLike = logL[worst],
                        step, spe, ref_par, pr_par, k);

    obj_par[worst,] = newvalue$theta;

    logL[worst]     = newvalue$pselike;

    logP[worst]     = newvalue$pseprior;

    like[worst]     = newvalue$like;

  }

  t1   = proc.time(); # ending time

  time = (t1 - t0)[1] # calculating time

  info = list(tree = tree, N = N, k = k, a = a, b = b,
              al = al, bl = bl, step = step, tol = tol,
              model = model, nis = "fixed_tree");

  return( list(logZ = logZ, H = H, logPsLd = logPsLd,
               logPsPd = logPsPd, logLd = logLd,
               dTheta = dTheta, seq_lz = s_z[-1],
               time = time, info = info) );

} # end function
