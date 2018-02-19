#' @title General functions
#' @description These functions are transversally used by the rest of functions
#' @param
#' @return
#' @export

##################################
##### log(x+y)=log(x)+log(y) #####
##################################

log.plus <- function(x,y)
{
  if(x>y) x + log(1+exp(y-x))
  else    y + log(1+exp(x-y))
}

### Vector version ###

logplusvec = function(x){

  r = -Inf;

  for(i in x){

    r = log.plus(r, i);

  }

  return(r);

}

### ###
