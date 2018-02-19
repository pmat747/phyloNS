#' @title Variable tree topology Nested Sampling in phylogenetics
#' @description It allows the use of Nested Sampling (NS) in phylogenetics for the variable tree topology case.  NS is basically used to estimate the marginal likelihood, quantity used for model selection under a Bayesian framework.
#' @param data An alignment, object of class phyDat
#' @return list
#' @export vns

##################################
##### log(x+y)=log(x)+log(y) #####
##################################
#
# log.plus <- function(x,y)
# {
#   if(x>y) x + log(1+exp(y-x))
#   else    y + log(1+exp(x-y))
# }

# logplusvec = function(x){
#
#   r = -Inf;
#
#   for(i in x){
#
#     r = log.plus(r, i);
#
#   }
#
#   return(r);
#
# }

##############################
### Positions for sampling ###
##############################

# It yields the positions and its weigths
#   for sampling according to the model.

pos_sam_vt = function(model, n_brs, t_pos, m = 1){

  # m  = 5 effort on proposal trees respect to other parameters

  if(model == "JC"){

    pos  = c(13, 14:(13 + n_brs), t_pos);
    prob = c(rep(1, 1 + n_brs), m);

  }else if(model == "JC_G"){

    pos  = c(1, 13, 14:(13 + n_brs), t_pos);
    prob = c(rep(1, 2 + n_brs), m);

  }else if(model == "HKY"){

    pos  = c(2:5, 7, 12, 13, 14:(13 + n_brs), t_pos);
    prob = c(rep(1, 7 + n_brs), m);

  }else if(model == "HKY_G"){

    pos  = c(1, 2:5, 7, 12, 13, 14:(13 + n_brs), t_pos);
    prob = c(rep(1, 8 + n_brs), m);

  }else if(model == "GTR"){

    pos  = c(2:5, 6:10, 12, 13, 14:(13 + n_brs), t_pos);
    prob = c(rep(1, 11 + n_brs), m);

  }else if(model == "GTR_G"){

    pos = c(1, 2:5, 6:10, 12, 13, 14:(13 + n_brs), t_pos);
    prob = c(rep(1,12 + n_brs), m);

  }

  return( list(pos = pos, prob = prob / sum(prob)) );

}

##########################
### Frequency of Trees ###
##########################

# library(ape)
# It yieds the the frequency for a set of trees

freq_tree = function(trees){

  tree = list();

  class(tree) = "multiPhylo";

  freq = NULL;

  i = 0;

  repeat{

    i       = i + 1;

    tree1   = trees[[1]];

    fun     = function (x) all.equal(x, tree1, use.edge.length = FALSE);

    x       = sapply(trees, fun);

    pos     = which(x == TRUE);

    freq[i] = length(pos);

    tree[[i]]  = tree1;

    #trees[pos] = NULL; # eliminating counted elements
    trees = trees[-pos]; # eliminating counted elements

    if( length(trees) == 0) break;

  }

  n = sum(freq); # total

  # reordering trees according to the maximum posterior tree

  maxTree = tree[[ which(freq == max(freq)) ]];

  maxTree = read.tree(text = write.tree(ladderize(maxTree)));

  tree = lapply(tree, function(x)
    read.tree(text =
                write.tree(rotateConstr(x, maxTree$tip.label))));

  class(tree) = "multiPhylo";

  return( list( tree = tree, freq = freq, pbb = freq/n, n = n))

}

##############
### Priors ###
##############

### JC ###

log_prior_JC_vt = function(tree, vec, a , b, al, bl, n_brs){

  if( any(vec < 0) || any(tree$edge.length < 0)){

    return( c(-Inf) )

  }

  # mu = vec[13];

  # branches + mean branch

  p1 = -(a + n_brs + 1) * log(vec[13]) - (sum(tree$edge.length) + b) / vec[13];

  return( p1 );

}

### HKY ###

log_prior_HKY_vt = function(tree, vec, a , b, al, bl, n_brs){

  if( any(vec < 0) || any(tree$edge.length < 0)){

    return( c(-Inf) )

  }

  # mu = vec[13];

  # branches + mean branch

  p1 = -(a + n_brs + 1) * log(vec[13]) - (sum(tree$edge.length) + b) / vec[13];

  # Rates: q = vec[7] (HKY) = vec[6:10] (GTR); phi = vec[12];

  p  = log(vec[12]) - vec[12] * (vec[7] + 1);

  return( p1 + p );

}

### GTR ###

log_prior_GTR_vt = function(tree, vec, a , b, al, bl, n_brs){

  if( any(vec < 0) || any(tree$edge.length < 0)){

    return( c(-Inf) )
  }

  # mu = vec[13];

  # branches + mean branch

  p1 = -(a + n_brs + 1) * log(vec[13]) - (sum(tree$edge.length) + b) / vec[13];

  # Rates: q = vec[7] (HKY) = vec[6:10] (GTR); phi = vec[12];

  p  = 5 * log(vec[12]) - vec[12] * ( sum(vec[6:10]) + 1 );

  return( p1 + p );

}

### JC_G ###

log_prior_JC_G_vt = function(tree, vec, a , b, al, bl, n_brs){

  if( any(vec < 0) || any(tree$edge.length < 0)){

    return( c(-Inf) )
  }

  # mu = vec[13];

  # branches + mean branch

  p1 = -(a + n_brs + 1) * log(vec[13]) - (sum(tree$edge.length) + b) / vec[13];

  # Rates: q = vec[7] (HKY) = vec[6:10] (GTR); phi = vec[12];

  # Gamma parameter: lambda = vec[1];

  p2     = (al - 1.0) * log(vec[1]) - vec[1] / bl;

  return( p1 + p2);

}

### HKY+G ###

log_prior_HKY_G_vt = function(tree, vec, a , b, al, bl, n_brs){

  if( any(vec < 0) || any(tree$edge.length < 0)){

    return( c(-Inf) )
  }

  # mu = vec[13];

  # branches + mean branch

  p1 = -(a + n_brs + 1) * log(vec[13]) - (sum(tree$edge.length) + b) / vec[13];

  # Rates: q = vec[7] (HKY) = vec[6:10] (GTR); phi = vec[12];

  # Gamma parameter: lambda = vec[1];

  p2     = (al - 1.0) * log(vec[1]) - vec[1] / bl;

  p  = log(vec[12]) - vec[12] * (vec[7] + 1);

  return( p1 + p2 + p );

}

### GTR+G ###

log_prior_GTR_G_vt = function(tree, vec, a , b, al, bl, n_brs){

  if( any(vec < 0) || any(tree$edge.length < 0)){

    return( c(-Inf) )
  }

  # mu = vec[13];

  # branches + mean branch

  p1 = -(a + n_brs + 1) * log(vec[13]) - (sum(tree$edge.length) + b) / vec[13];

  # Rates: q = vec[7] (HKY) = vec[6:10] (GTR); phi = vec[12];

  # Gamma parameter: lambda = vec[1];

  p2     = (al - 1.0) * log(vec[1]) - vec[1] / bl;

  p  = 5 * log(vec[12]) - vec[12] * ( sum(vec[6:10]) + 1 );

  return( p1 + p2 + p );

}

##########################
### Generator function ###
##########################

proposal_vt <- function(vec, tree, pos, model, t_pos,
                     pbb_rSPR, rooted, prior_width){

  if( (pos >= 14) && (pos < t_pos) ){

    out = abs(tree$edge.length[pos-13] + prior_width * 10.0^(1.5 -
                                                               3.0 * abs(rt(1, df = 2))) * rnorm(1));

    tree$edge.length[pos-13] = out; # Branches

  }else if( pos == t_pos){

    if( runif(1) < pbb_rSPR ){

      mov_k = NULL; # rSPR

    }else{

      mov_k = 1; # rNNI

    }

    if(rooted){

      repeat{

        ptree = rSPR(tree, move = 1, k = mov_k);

        if(is.rooted(ptree)){ # selecting only a rooted tree

          tree = ptree;
          break

        }
      }

    }else{

      tree = rSPR(tree, move = 1, k = mov_k); # tree

    }
  }else{

    # abs acts as reflection for negative values
    out = abs(vec[pos] + prior_width * 10.0^(1.5 -
                                               3.0 * abs(rt(1, df = 2))) * rnorm(1));

    if( (pos >= 2) && ( pos <= 5) ){

      vec[pos] = out; # frequencies

      vec[2:5] = vec[2:5] / sum(vec[2:5]);

    }else if( (pos >=6) && ( pos <= 11)){

      vec[pos] = out; # rates

      if((model == "HKY") || (model == "HKY_G")){

        vec[10] = out;

      }

    }else if( (pos == 12) || (pos == 13) ){

      vec[pos] = out; # mu & phi

    }else if( pos == 1 ){

      while( out < 0.10 ){ # Gamma parameter (0.14)

        out = abs(vec[pos] + prior_width * 10.0^(1.5 -
                                                   3.0 * abs(rt(1, df = 2))) * rnorm(1))

      }

      vec[pos] = out;

    }
  }

  return(list( theta = vec, tree = tree));

}

genValue_vt <- function(data, vec0, l0, worstLike,
                        a, b, al, bl, k, step, spe){

  # Random numbers

  u = runif(step);

  pos = sample(x = spe$pos, size = step, replace = TRUE, prob = spe$prob);

  # Starting Value

  p0   = spe$log_prior_vt(vec0$tree, vec0$theta, a, b, al, bl, spe$n_brs); # prior

  # loop

  i = 0;

  for(ipos in pos){

    if( ipos == spe$t_pos){ # NULL

    }else if(ipos > 13){

      r = 1/vec0$theta[13];

      prior_width = qexp(0.99, rate = r ) - qexp(0.01, rate = r);

    }else if(ipos == 13){ prior_width = spe$mu_w;

    }else if(ipos == 1){ prior_width = spe$lambda_w ;

    }else if((ipos >= 2) && (ipos <= 5)){ prior_width = 0.77; # 4*sqrt(3/(16*5))

    }else if((ipos >= 6) && (ipos <= 11)){

      r = vec0$theta[12];

      prior_width = qexp(0.99, rate = r) - qexp(0.01, rate = r);

    }else if(ipos == 12){ prior_width = 4;} # phi~Exp(1)=>dexp(0.99)-dexp(.01)>6

    vec1 = proposal_vt(vec0$theta, vec0$tree, ipos,
                       spe$model, spe$t_pos, spe$pbb_rSPR,
                       spe$rooted, prior_width);

    if( (ipos == 12) || (ipos == 13) ){

      l1 = l0; # phi & mu do not alter the likelihood

    }else{

      l1   = pml(vec1$tree, data, Q = vec1$theta[6:11],
                 bf = vec1$theta[2:5], shape = vec1$theta[1],
                 k = k)$logLik;
    }

    p1   = spe$log_prior_vt(vec1$tree, vec1$theta, a, b, al, bl, spe$n_brs);

    i    = i + 1;

    if( ( u[i] < exp(p1 - p0) ) &&
        ( worstLike < l1 ) ){

      # this avoids copying too many unneccesary elements

      if( ipos == spe$t_pos ){ # tree

        vec0$tree = vec1$tree;

      }else if( (ipos >= 14) && (ipos < spe$t_pos) ){ # branch

        vec0$tree$edge.length[ipos-13] = vec1$tree$edge.length[ipos-13];

      }else if( (ipos >= 2) && (ipos <= 5) ){ # frequency

        vec0$theta[2:5] = vec1$theta[2:5];

      }else if((ipos == 1) || (ipos == 12) ||
               (ipos == 13) ){ # Gamma par, mu & phi

        vec0$theta[ipos] = vec1$theta[ipos];

      }else if( (ipos >= 6) && (ipos <= 11)){ # Rates for HKY & GTR versions

        if( (spe$model == "HKY") || (spe$model == "HKY_G") ){

          vec0$theta[c(7, 10)] = vec1$theta[ipos];

        }else{ # GTR model

          vec0$theta[ipos] = vec1$theta[ipos];

        }

      }

      l0   = l1;
      p0   = p1;

    }

  } # ending step

  return( list( theta = vec0$theta, tree = vec0$tree, like = l0 ) );

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

vns = function( data, N = 100, k = 4, a = 3, b = 0.2, al = NULL, bl = NULL,
                step = 50, max.iter = 2000, end = 2, act_plot = TRUE, model,
                pbb_rSPR = 0.5){

  # checking if the model is defined

  models = c("JC", "JC_G", "HKY", "HKY_G", "GTR", "GTR_G");

  if( !any(model == models) )stop("The model is not defined");

  rooted = FALSE; # This option could be useful in the future

  if( rooted == FALSE){

    if( length(data) <= 3 )
      stop("Number of taxa must be greater than 3 for unrooteed trees");
  }

  if( any(c(a,b,al,bl) < 0) ) stop("Prior parameters must be positive");

  t0     = proc.time(); # starting time

  ### Priors definitions ###

  if(model == "JC"){          log_prior_vt = log_prior_JC_vt;

  }else if(model == "HKY"){   log_prior_vt = log_prior_HKY_vt;

  }else if(model == "GTR"){   log_prior_vt = log_prior_GTR_vt;

  }else if(model == "JC_G"){  log_prior_vt = log_prior_JC_G_vt;

  }else if(model == "HKY_G"){ log_prior_vt = log_prior_HKY_G_vt;

  }else if(model == "GTR_G"){ log_prior_vt = log_prior_GTR_G_vt;

  }

  ###

  n_taxa = length(data); # number of taxa/leaves

  n_brs  = ifelse(rooted, 2*n_taxa - 2, 2*n_taxa - 3); # number of branches (rooted-unrooted)

  tips   = names(data); # names of species

  spe    = list();

  ### Sampling the prior according to the model ###

  # Branch lengths

  mu = as.matrix(rinvgamma(N, shape = a, scale = b), ncol = 1);

  t  = t(apply(mu, 1, function(x)rexp(n_brs, rate = 1 / x) ));

  m  = unlist(strsplit(model, ""));

  ### Gamma parameter ###

  if(m[length(m)] != "G"){

    lambda = rep(1, N); # JC - HKY

    k = 1;

  }else{

    if( is.null(al) || is.null(bl) )
      stop('Provide full specified prior for the gamma parameter');

    lambda = rgamma(N, shape = al, scale = bl); # JC_G HKY_G GTR_G

    for(i in 1:N){ # It prevents phangorn from crushing

      while( lambda[i] < 0.07 ){ # It could even crush with 0.07464181

        lambda[i] = rgamma(1, shape = al, scale = bl)

      }
    }

    lambda_w = qgamma(0.999, shape = al, scale = bl, lower.tail = TRUE) -
      qgamma(0.001, shape = al, scale = bl, lower.tail = TRUE); # Gamma par

    spe$lambda_w = lambda_w;

  }

  ### Frequency ###

  if( m[1] != "J"){

    freq  = matrix( rexp( N * 4, 1), ncol = 4);

    freq  = freq / apply(freq, 1, sum); # Dirichlet(1,1,1,1)

  }else{

    freq   = matrix( 0.25, nrow = N, ncol = 4); # JC JC_G

  }

  ### Rates ###

  if( m[1] == "J" ){

    q   = matrix(1, nrow = N, ncol = 6);
    phi = rep(1, N)

  }else if( m[1] == "H" ){

    phi =  as.matrix(rexp(N), ncol = 1);

    q   = c(apply(phi, 1, function(x)rexp(n = 1, rate = x)));

    q   = cbind(rep(1, N), q, rep(1, N), rep(1, N), q, rep(1, N));

  }else if( m[1] == "G" ){

    phi =  as.matrix(rexp(N), ncol = 1);

    q   = t(apply(phi, 1, function(x)rexp(n = 5, rate = x)));

    q   = cbind(q, rep(1,N));

  }

  ### end sampling from the prior ###

  # specifications

  t_pos = 14 + n_brs; # position ocupied by the tree (13+1+int_branches)

  v   = pos_sam_vt(model, n_brs, t_pos, m = 3); # m = effort in tree sampling

  # prior width

  mu_w = rinvgamma(2000, shape = a, scale = b)

  mu_w = max(mu_w) - min(mu_w); # mu

  spe$n_taxa = n_taxa;
  spe$n_brs  = n_brs;
  spe$t_pos  = t_pos;
  spe$mu_w   = mu_w;
  spe$pbb_rSPR = pbb_rSPR;
  spe$log_prior_vt = log_prior_vt;
  spe$rooted = rooted;
  spe$model  = model;
  spe$pos  = v$pos;
  spe$prob = v$prob;

  ### Objects: parameters & trees ###

  Obj_par   = cbind(lambda, freq, q, phi, mu); # parameter

  Obj_trees = list();

  class(Obj_trees) = "multiPhylo";

  logL    = NULL;

  for(i in 1:N){

    ### prior sample trees ###

    ptree = rtree(n_taxa, rooted = rooted);
    ptree$tip.label = sample(tips, n_taxa); # branch names

    ptree$edge.length = t[i, ]; # branch length

    logL[i] = pml(ptree, data, Q = Obj_par[i,6:11], bf = Obj_par[i,2:5],
                  shape = Obj_par[i,1], k = k)$logLik;

    Obj_trees[[i]] = ptree;

  }

  # definitions

  iter      = 0;
  #logwidth1 = log( 1.0 - exp( - 1.0 / N ) ); # simple
  logwidth1 = log(1 - exp( - 2 / N )) - log(2); # Trapezoidal
  H         = 0;
  dTheta    = NULL; # discarded values
  dTrees    = list(); # discarded trees
  class(dTrees) = "multiPhylo";
  lw        <- NULL; # like * prior
  logX      <- NULL; # prior mass
  logLd     <- NULL; # discarded log-likelihoods
  s_z       <- NULL; # sequence of z-evolution
  s_z[1]    <- -1.79769e+308;

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

    dTheta = rbind(dTheta, Obj_par[worst,] ); # parameters

    dTrees[[iter]] = Obj_trees[[worst]]; # tree

    if( (iter >= end * N * H) || (iter >= max.iter) ) break;

    # starting value for new active point #

    #init = sample(c(1:N)[-worst], size = 1);
    init = ifelse(N == 1, 1, sample(c(1:N)[-1], 1));

    vec0 = list(theta = Obj_par[init,], tree = Obj_trees[[init]]);

    l0   = logL[init];

    values = genValue_vt(data, vec0, l0, worstLike = logL[worst],
                         a, b, al, bl, k, step, spe);

    Obj_par[worst, ]  = values$theta;

    Obj_trees[[worst]]= values$tree;

    logL[worst]       = values$like;

    }

  logZc   <- log.plus(logZ, log(mean(exp(logL))) - iter/N )

  Lw.z    <- exp(lw - logZc)

  max.post<- exp(- sum(Lw.z*log(Lw.z), na.rm = TRUE)) # number of posterior samples

  Lw.z    <- Lw.z / max(Lw.z)

  d       <- sample(iter, replace = TRUE, size = max.post, prob = Lw.z )

  ### Posterior sampling ###

  sampled_par <- matrix(NA, nrow = max.post, ncol = 13 + n_brs);

  for( i in 1:max.post){

    sampled_par[i, ] = c(dTheta[d[i],], dTrees[[d[i]]]$edge.length);

  }

  # Naming posterior samples & discarded points
  colnames(sampled_par) = c('lambda', 'freqA', 'freqC', 'freqG', 'freqT',
                            'qCA', 'qGA', 'qTA', 'qGC', 'qTC', 'qTG',
                            'phi', 'mu',
                            paste(rep('t', n_brs), 1:n_brs, sep = ""));

  colnames(dTheta) = c('lambda', 'freqA', 'freqC', 'freqG', 'freqT',
                       'qCA', 'qGA', 'qTA', 'qGC', 'qTC', 'qTG',
                       'phi', 'mu');

  # Posterior trees

  sampled_trees    = list();

  samp_tree_length = NULL;

  sampled_likes = logLd[d];

  class(sampled_trees) <- "multiPhylo";

  for(i in 1:max.post){

    sampled_trees[[i]]  = dTrees[[ d[i] ]]; # storing tree

    samp_tree_length[i] = sum(dTrees[[ d[i] ]]$edge.length); # tree length

  }

  # Maximum posterior tree #

  freqTree = freq_tree(sampled_trees)

  maxTree = which(freqTree$freq == max(freqTree$freq) )

  maxTree = freqTree$tree[[maxTree]]

  meanBranch = NULL

  fun = function(x) all.equal(x, maxTree, use.edge.length = FALSE);

  if( rooted == TRUE ){

    for(i in 1:length(sampled_trees)){

      if( fun(sampled_trees[[i]]) ){

        aux = rotateConstr(sampled_trees[[i]], maxTree$tip.label);

        # aux = read.tree(text = write.tree(aux));

        meanBranch = rbind(meanBranch, aux$edge.length);

      }
    }

  }else{ # unrooted tree

    for(i in 1:length(sampled_trees)){

      if( fun(sampled_trees[[i]]) ){

        aux = root(sampled_trees[[i]], maxTree$tip.label[1], resolve.root = F);

        aux = rotateConstr(aux, maxTree$tip.label);

        aux = read.tree(text = write.tree(aux));

        meanBranch = rbind(meanBranch, sampled_trees[[i]]$edge.length);

      }
    }
  }

  meanBranch = apply(meanBranch, 2, mean);

  maxTree$edge.length = meanBranch;

  ###

  t1   = proc.time(); # ending time

  time = (t1 - t0)[1] # calculating time

  ### Plots ###

  if(act_plot){

    par(mar = c(4,4,4,4));

    plot(logX,logLd, xlab = expression("log"* xi), ylab='log(L)', pch='.');

    plot(logX,Lw.z, xlab = expression("log"* xi), ylab = 'Weight for posterior');

    plot(samp_tree_length, type = 'l', xlab = 'Sample',
         ylab = 'Posterior tree length')

    plot(sampled_likes, type = 'l', xlab = 'Sample',
         ylab = 'Log-likelihood')

    if(length(data) > 4 ){

      # consensus net is defined for more than 4 taxa

      cnt  = consensusNet(sampled_trees); # phangorn

      plot(cnt, "2", tip.color = "black", edge.width = 1, font = 1,
           show.nodes = TRUE, type = "2D");

      title("Consensus Network");

    }

    plotTree(ladderize(maxTree), main = "Maximum posterior tree"); # phytool package ape respectively
    add.scale.bar(cex = 0.7, font = 1) # ape package

  }

  # Storing process info
  info = list(N = N, k = k, a = a, b = b, al = al, bl = bl,
              step = step, model = model, ns = "variable_tree")

  return(list(logZ = logZ, logZc = logZc,
              sd_logZ = sqrt(H / N),
              H = H, iter = iter, time = time,
              sampled_par = sampled_par,
              sampled_trees = sampled_trees,
              samp_tree_length = samp_tree_length,
              sampled_likes = sampled_likes,
              maxTree = maxTree,
              logX = logX, logLd = logLd,
              Lw.z = Lw.z, seq_lz = s_z[-1],
              dTheta = dTheta, dTrees = dTrees,
              info = info
              ) )
}
