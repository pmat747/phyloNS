#' @title Prior and reference distributions for fixed tree case
#' @description These functions are used NIS
#' @param
#' @return
#' @export

##########################
### Priors & Ref distb ###
##########################

# Prior and pseudo prior functions do not check parameter space.
# This is checked in the pseudo-likelihood function.

##########
### JC ###
##########

### Prior ###

dprior_JC = function(vec, pr_par){

  # branches: vec[14: ...]

  sum(stats::dexp(vec[-c(1:13)], rate = 1/vec[13], log = TRUE)) +

    # mean branch: mu = vec[13];

    log(MCMCpack::dinvgamma(vec[13], shape = pr_par$a, scale = pr_par$b));

}

### Ref ###

dpseprior_JC = function(vec, ref_par, ref_trees){

  # branches: vec[14: ...]

  sum(stats::dgamma(vec[-c(1:13)], shape = ref_par$at_r,
             scale = ref_par$bt_r, log = TRUE)) +

    #  mean branch: mu = vec[13];

    log(MCMCpack::dinvgamma(vec[13], shape = ref_par$a_r, scale = ref_par$b_r) );

}

###########
### HKY ###
###########

### prior ###

dprior_HKY = function(vec, pr_par){

  # branches: vec[14: ...]

  sum(stats::dexp(vec[-c(1:13)], rate = 1/vec[13], log = TRUE)) +

    # mean branch: mu = vec[13];

    log(MCMCpack::dinvgamma(vec[13], shape = pr_par$a, scale = pr_par$b)) +

    # Rates: q = vec[7] (HKY); phi = vec[12];

    log(vec[12]) - vec[12] * (vec[7] + 1) +

    # Frequencies: vec[2:5]

    lgamma(4); # lgamma(4) log( ddirichlet(vec[2:5], rep(1.0, 4)) );

}

### Ref ###

dpseprior_HKY = function(vec, ref_par){

  # branches: vec[14: ...]

  sum(stats::dgamma(vec[-c(1:13)], shape = ref_par$at_r,
             scale = ref_par$bt_r, log = TRUE)) +

    #  mean branch: mu = vec[13];

    log(MCMCpack::dinvgamma(vec[13], shape = ref_par$a_r, scale = ref_par$b_r) ) +

    # Rates: q = vec[7];

    stats::dgamma(vec[7], shape = ref_par$aq_r, scale = ref_par$bq_r, log = TRUE) +

    # rate mean: phi = vec[12]

    stats::dgamma(vec[12], shape = ref_par$aphi_r, scale = ref_par$bphi_r, log = TRUE) +

    # Frequencies: vec[2:5]

    lgamma(sum(ref_par$f_r)) - sum(lgamma(ref_par$f_r)) +
    sum((ref_par$f_r - 1 ) * log(vec[2:5]));

}

############
### GTR ####
############

### prior ###

dprior_GTR = function(vec, pr_par){

  # branches: vec[14: ...]

  sum(stats::dexp(vec[-c(1:13)], rate = 1/vec[13], log = TRUE)) +

    # mean branch: mu = vec[13];

    log(MCMCpack::dinvgamma(vec[13], shape = pr_par$a, scale = pr_par$b) ) +

    # Rates: q = vec[6:10] (GTR); phi = vec[12];

    5.0 * log(vec[12]) - vec[12] * ( sum(vec[6:10]) + 1.0 ) +

    # Frequencies: vec[2:5]

    lgamma(4); # lgamma(4) log( ddirichlet(vec[2:5], rep(1.0, 4)) );

}

### Ref ###

dpseprior_GTR = function(vec, ref_par){

  # branches: vec[14: ...]

  sum(stats::dgamma(vec[-c(1:13)], shape = ref_par$at_r,
             scale = ref_par$bt_r, log = TRUE)) +

    #  mean branch: mu = vec[13];

    log(MCMCpack::dinvgamma(vec[13], shape = ref_par$a_r, scale = ref_par$b_r) ) +

    # Rates: q = vec[6:10] (GTR)

    sum(stats::dgamma(vec[6:10], shape = ref_par$aq_r, scale = ref_par$bq_r, log = TRUE)) +

    # rate mean: phi = vec[12]

    stats::dgamma(vec[12], shape = ref_par$aphi_r, scale = ref_par$bphi_r, log = TRUE) +

    # Frequencies: vec[2:5]

    lgamma(sum(ref_par$f_r)) - sum(lgamma(ref_par$f_r)) +
    sum((ref_par$f_r - 1 ) * log(vec[2:5]));

}

############
### JC_G ###
############

### Prior ###

dprior_JC_G = function(vec, pr_par){

  # branches: vec[14: ...]

  sum(stats::dexp(vec[-c(1:13)], rate = 1/vec[13], log = TRUE)) +

    # mean branch: mu = vec[13];

    log(MCMCpack::dinvgamma(vec[13], shape = pr_par$a, scale = pr_par$b)) +

    # Gamma parameter: lambda = vec[1];

    stats::dgamma(vec[1], shape = pr_par$al, scale = pr_par$bl, log = TRUE);

}

### Ref ###

dpseprior_JC_G = function(vec, ref_par){

  # branches: vec[14: ...]

  sum(stats::dgamma(vec[-c(1:13)], shape = ref_par$at_r,
             scale = ref_par$bt_r, log = TRUE)) +

    #  mean branch: mu = vec[13];

    log(MCMCpack::dinvgamma(vec[13], shape = ref_par$a_r, scale = ref_par$b_r) ) +

    # Gamma parameter: lambda = vec[1]

    stats::dgamma(vec[1], shape = ref_par$al_r, scale = ref_par$bl_r, log = TRUE);

}

#############
### HKY+G ###
#############

### Prior ###

dprior_HKY_G = function(vec, pr_par){

  # branches: vec[14: ...]

  sum(stats::dexp(vec[-c(1:13)], rate = 1/vec[13], log = TRUE)) +

    # mean branch: mu = vec[13];

    log(MCMCpack::dinvgamma(vec[13], shape = pr_par$a, scale = pr_par$b) ) +

    # Gamma parameter: lambda = vec[1];

    stats::dgamma(vec[1], shape = pr_par$al, scale = pr_par$bl, log = TRUE) +

    # Rates: q = vec[7] (HKY); phi = vec[12];

    log(vec[12]) - vec[12] * ( sum(vec[7]) + 1.0 ) +

    # Frequencies: vec[2:5]

    lgamma(4); # lgamma(4) log( ddirichlet(vec[2:5], rep(1.0, 4)) );

}

### Ref ###

dpseprior_HKY_G = function(vec, ref_par){

  # branches: vec[14: ...]

  sum(stats::dgamma(vec[-c(1:13)], shape = ref_par$at_r,
             scale = ref_par$bt_r, log = TRUE)) +

    #  mean branch: mu = vec[13];

    log(MCMCpack::dinvgamma(vec[13], shape = ref_par$a_r, scale = ref_par$b_r) ) +

    # Gamma parameter: lambda = vec[1]

    stats::dgamma(vec[1], shape = ref_par$al_r, scale = ref_par$bl_r, log = TRUE) +

    # Rates: q = vec[7]

    stats::dgamma(vec[7], shape = ref_par$aq_r, scale = ref_par$bq_r, log = TRUE) +

    # rate mean: phi = vec[12]

    stats::dgamma(vec[12], shape = ref_par$aphi_r, scale = ref_par$bphi_r, log = TRUE) +

    # Frequencies: vec[2:5]

    lgamma(sum(ref_par$f_r)) - sum(lgamma(ref_par$f_r)) +
    sum((ref_par$f_r - 1 ) * log(vec[2:5]));

}

#############
### GTR+G ###
#############

### Prior ###

dprior_GTR_G = function(vec, pr_par){

  # branches: vec[14: ...]

  sum(stats::dexp(vec[-c(1:13)], rate = 1/vec[13], log = TRUE)) +

    # mean branch: mu = vec[13];

    log(MCMCpack::dinvgamma(vec[13], shape = pr_par$a, scale = pr_par$b) ) +

    # Gamma parameter: lambda = vec[1];

    stats::dgamma(vec[1], shape = pr_par$al, scale = pr_par$bl, log = TRUE) +

    # Rates: q = vec[6:10] (GTR); phi = vec[12];

    5.0 * log(vec[12]) - vec[12] * ( sum(vec[6:10]) + 1.0 ) +

    # Frequencies: vec[2:5]

    lgamma(4); # lgamma(4) log( ddirichlet(vec[2:5], rep(1.0, 4)) );

}

### Ref ###

dpseprior_GTR_G = function(vec, ref_par){

  # branches: vec[14: ...]

  sum(stats::dgamma(vec[-c(1:13)], shape = ref_par$at_r,
             scale = ref_par$bt_r, log = TRUE)) +

    #  mean branch: mu = vec[13];

    log(MCMCpack::dinvgamma(vec[13], shape = ref_par$a_r, scale = ref_par$b_r) ) +

    # Gamma parameter: lambda = vec[1]

    stats::dgamma(vec[1], shape = ref_par$al_r, scale = ref_par$bl_r, log = TRUE) +

    # Rates: q = vec[6:10] (GTR)

    sum(stats::dgamma(vec[6:10], shape = ref_par$aq_r, scale = ref_par$bq_r, log = TRUE)) +

    # rate mean: phi = vec[12]

    stats::dgamma(vec[12], shape = ref_par$aphi_r, scale = ref_par$bphi_r, log = TRUE) +

    # Frequencies: vec[2:5]

    lgamma(sum(ref_par$f_r)) - sum(lgamma(ref_par$f_r)) +
    sum((ref_par$f_r - 1 ) * log(vec[2:5]));

}
