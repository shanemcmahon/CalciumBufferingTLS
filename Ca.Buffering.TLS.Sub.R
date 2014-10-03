
if(interactive.file.chooser){
  print("All data files are assumed to be in comma seperated value (csv) format. The first row is assumed contain column names. For details regarding the format see the documentation.")
  
  print("Select the csv file containing the fluorescence data.")
  fluorescence.filename <- choose.files(caption = "Select the csv file containing the fluorescence data.")
  
  print("Select the csv file containing standard errors for the fluorescence measurements.")
  fluorescence.se.filename <- choose.files(caption = "Select the csv file containing standard errors for the fluorescence measurements.")
  
  print("Select the csv file containing the change in total calcium concentration data.")
  delta.ca.total.filename <- choose.files(caption = "Select the csv file containing the change in total calcium concentration data.")
  print("Select the csv file containing standard errors for the change in total calcium measurements.")
  delta.ca.total.se.filename <- choose.files(caption = "Select the csv file containing standard errors for the change in total calcium measurements.")  
}

#############################################################################################################################
#############################################################################################################################

# End Section 1

#############################################################################################################################
#############################################################################################################################


#############################################################################################################################
#############################################################################################################################

# Begin Section 2: function library

#############################################################################################################################
#############################################################################################################################

require(compiler); require(foreach); require(doParallel); require(minpack.lm); require(mosaic)

String.Length <- function(str.1){
  i <- 1
  repeat{
    if(substring(str.1,i,i+1)=="")
      return(i-1)
    i <- i +1
  }
}

T.L.S.Distance <- function(x.error, beta, x, y, x.weight, y.weight, Fxn, return.vector=FALSE,x.max=NULL,x.min=NULL, return.matrix = FALSE, ...){
  
  if(is.null(dim(x))){
    x <- matrix(x,nrow=1)
    x.error <- matrix(x.error, nrow=1)
    x.weight <- matrix(x.weight, nrow=1)
    y.weight <- matrix(y.weight, nrow=1)
  }
  
  
  y.hat <- Fxn(x+x.error,beta, ...)
  #   y.hat[!is.finite(y.hat)] <- .Machine$double.xmax^.25
  if(any(!is.finite(y.hat)) | any((x + x.error > x.max)) | any((x + x.error < x.min))){
    y.hat <- y.hat + .Machine$double.xmax^.25
  }
  y.error <- y.hat - y
  out <- NULL  
  x.out <- NULL  
  
  
  
  
  if(is.null(ncol(x.weight))){
    for(i in 1:length(x.weight)){
      x.out <- cbind(x.out, x.error[,i]^2*x.weight[i]^2)
    }
  }else{
    x.out <- (x.error^2*x.weight^2)
  }
  
  y.out <- NULL  
  if(is.null(ncol(y.weight))){
    
    for(i in 1:length(y.weight)){
      y.out <- cbind(y.out, y.error[,i]^2*y.weight[i]^2)
    }
  }else{
    y.out <- (y.error^2*y.weight^2)
  }
  
  
  
  
  out <- cbind(y.out,x.out)
  if(return.matrix){
    return(out)
  }
  out <- rowSums(out)
  
  if(!return.vector){
    return(sum(out)^.5)
  }
  return(out^.5)
}
T.L.S.Distance.C <- cmpfun(T.L.S.Distance)

T.L.S <- function(beta, ifixb, ifixx, y, x, x.weight, y.weight, beta.lower=-Inf, beta.upper=Inf, return.vector=FALSE, Fxn, x.error, beta2 = NULL,tls.env , x.max=NULL, x.min=NULL, return.matrix = FALSE, ...){
  
  if(!is.null(beta2)){
    beta2[!ifixb] <- beta
    beta <- beta2
  }
  values <- array(0, nrow(x))
  ifixb <- as.logical(ifixb)
  beta <- unlist(beta)
  if(length(beta.lower)!=length(beta)){
    beta.lower <- rep(-Inf,length(beta))
  }
  if(length(beta.upper)!=length(beta)){
    beta.upper <- rep(Inf,length(beta))
  }
  if(!ifixx){
    if(length(y.weight) < length(y)){
      y.weight <- matrix(replicate(nrow(y),matrix(y.weight,nrow=1)),nrow=nrow(y))
    }
    if(length(x.weight) < length(x)){
      x.weight <- matrix(replicate(nrow(x),matrix(x.weight,nrow=1)),nrow=nrow(x))
    }
    for(data.index in seq(nrow(x))){
      x.err.start = .5 - x[data.index, ]
      
      fit <- optim(par=x.err.start, fn = T.L.S.Distance.C, beta = beta, x = x[data.index, ] , y = y[data.index, ],x.weight = x.weight[data.index,] ,y.weight = y.weight[data.index,] ,Fxn = Fxn, x.max = x.max, x.min = x.min, ...)
      values[data.index] <-  fit$value
      x.error[data.index,] <-  fit$par
    }
    
    tls.env$x.error <- x.error
  }else{ #ifixx
    values <- T.L.S.Distance.C(x.error = x.error,beta = beta,x = x,y = y,x.weight = x.weight,y.weight = y.weight, Fxn = Fxn, return.vector = TRUE,x.max = x.max, x.min = x.min, return.matrix = return.matrix ,...)
  } #ifixb
  
  
  if (any(beta[!ifixb] < beta.lower[!ifixb])){
    values <- values * (1 + sum(abs((beta.lower[!ifixb] - beta[!ifixb])[which(beta[!ifixb] < beta.lower[!ifixb])]) ) ) 
  }
  if (any(beta[!ifixb] > beta.upper[!ifixb])){
    values <- values * (1 + sum(abs((beta[!ifixb] - beta.upper[!ifixb])[which(beta[!ifixb] > beta.upper[!ifixb])]) ) )
  }
  
  
  if(return.vector==TRUE){
    return(values)
  }
  return(sum(values^2))
}
T.L.S.C <- cmpfun(T.L.S)

Delta.Ca.Total.Hat <- function(x,beta){
  #   x[,1] == final fluorescence
  #   x[,2] == final fluorescence
  # beta[1]   == dynamic range
  #beta[2]   == dye kd
  #   beta[3] == dye concentration
  #   beta[4] == endogenous buffer 1 kd
  #   beta[5] == endogenous buffer 1 concentration
  #   beta[6] == endogenous buffer 2 kd
  #   beta[7] == endogenous buffer 2 concentration
  #   beta[8] == nonsaturable endogenous buffer
  #   beta[9] == accessible volume fraction
  
  calcium.final   <- ((x[,2] - (1/beta[1]))/(1 - x[,2])) * beta[2]
  
  calcium.initial <- ((x[,1] - (1/beta[1]))/(1 - x[,1])) * beta[2]
  
  delta.ca.total.hat  <- (1+beta[8])*(calcium.final - calcium.initial) - beta[4] * beta[5]/(beta[4] + calcium.final) + beta[4] * beta[5]/(beta[4] + calcium.initial) - beta[6] * beta[7]/(beta[6] + calcium.final) + beta[6] * beta[7]/(beta[6] + calcium.initial) - beta[2] * beta[3]/(beta[2] + calcium.final) + beta[2] * beta[3]/(beta[2] + calcium.initial)
  
  delta.ca.total.hat <- cbind(delta.ca.total.hat * beta[9])
}
Delta.Ca.Total.Hat.C <- cmpfun(Delta.Ca.Total.Hat)

T.L.S.Gradient <- function(beta, x, y, x.error, x.weight, y.weight, beta2=NULL, beta.lower=-Inf, beta.upper=Inf, return.vector=TRUE, Distance.Function = T.L.S.Distance.C, Fxn, return.matrix=FALSE, F.Grad = NULL, ifixx = TRUE, ifixb, x.max = NULL, x.min = NULL, tls.env, grad.tol = .Machine$double.eps^.5, ...)
{
  gradient <- double(length(beta))
  
  if(!is.null(F.Grad)){
    gradient <-  colSums(2*(y - Fxn(x = x + x.error , beta = beta, ...)) * F.Grad(x = x + x.error, beta = beta, ifixb, ...))
    gradient[ifixb] <- 0
  }else{
    
    gradient <- NULL
    # no user supplied gradient
    for(par.index in 1:length(beta)){
      if(ifixb[par.index]){
        gradient <- c(gradient,0)
        next
      }
      delta <- double(length(beta))
      delta[[par.index]] <- abs(beta[[par.index]] * .Machine$double.eps^.5) + .Machine$double.eps^.5
      delta <- unlist(delta)
      
      grad0 <-  (T.L.S.C(beta = beta + delta, x = x, y = y, x.weight = x.weight, y.weight = y.weight, beta.lower = beta.lower, beta.upper = beta.upper, return.vector = FALSE,Fxn = Fxn, x.error = x.error, x.max = x.max, x.min = x.min, ifixx = ifixx, ifixb = as.logical(ifixb), tls.env = tls.env) - T.L.S.C(beta = beta - delta, x = x, y = y, x.weight = x.weight, y.weight = y.weight, beta.lower = beta.lower, beta.upper = beta.upper, return.vector = FALSE,Fxn = Fxn, x.error = x.error, x.max = x.max, x.min = x.min, ifixx = ifixx, ifixb = as.logical(ifixb), tls.env = tls.env) )/(2*delta[[par.index]])
      
      gradient <- c(gradient, grad0)
      
    }    
    
  }
  
  return(gradient)
}
T.L.S.Gradient.C <- cmpfun(T.L.S.Gradient)

# S.G.D <- function(beta, x, y, x.error=NULL, x.weight, y.weight, x.min = NULL, x.max = NULL, beta.lower, beta.upper, ifixb, Objective.Fxn = T.L.S.C, solution.tolerance=.Machine$double.eps^.5, max.iter = 5000, num.steps = 500, Grad = T.L.S.Gradient.C, decay.constant=.9995, data.env = new.env(), Fxn, parameter.boundary.margin = .01, p.restart = 1, p.restart.decay = .5,   p.restart.grow = 1.1, p.check.obj = .05, p.goto.best.par = .005, p.check.boundary = .01, p.rescale = .005, rescale.vector = NULL, ...){
#   #   Function: S.G.D (Stochastic Gradient Descent)
#   # S.G.D is a implementation of a stochastic gradient descent algorithm, in principle it can be used to minimize arbitrary user specefied functions but it was designed with the total least squares problem in mind
#   # Fxn: a user specefied function, here by default the function is a model of calcium buffering. the calling sequence for Fxn is Fxn(x, beta)
#   
#   # x: matrix of independent variables
#   #   y: matrix of dependent variables
#   # x.error: matrix of delta(x) values
#   #   x.weight: matrix of "weights" for x errors. The correct weight is the inverse of the standard error.
#   #   y.weight: matrix of "weights" for y errors. The correct weight is the inverse of the standard error.
#   #   x.min: lower limit for x values. may be a vector or scalar. during the minimization procedure we will require that x + x.error > x.min
#   #  x.max: analagous to x.min
#   # beta.lower: lower limit for parameters
#   # beta.upper: upper limit for parameters
#   # ifixb: logical vector with ifixb[i] == TRUE indicating the the ith component of beta is fixed
#   # Objective.Fxn: the function to be minized, by default a model of calcium buffering
#   # solution.tolerance: used to check convergence and whether two numbers are "equal"
#   # max.iter: maximum number of iterations to perform
#   # num.steps: sets the initial learning rate eta = (beta.upper - beta.lower)/num.steps. Here eta is the learning rate. That is, the update rule for each iteration is beta = beta + eta * grad(Objective.Fxn)
#   # Grad: function for computing the gradient of Objective.Fxn
#   # decay.constant: controls the rate at which the learning rate eta "decays" for each iteration we have eta = eta * decay.constant
#   # data.env: an R environment that is used to store some working data
#   # Fxn: a user supplied function to be fit
#   # parameter.boundary.margin: controls how close to the boundary the solution is allowed to be before restarting. If at any point we have abs(beta - beta.upper) < abs(beta.upper - beta.lower) * parameter.boundary.margin or abs(beta - beta.upper ) < abs(beta.upper - beta.lower) * parameter.boundary.margin, the beta is reset either randomly or the the best parameters so far
#   # p.restart: if the parameters are "too close to the boundary" as defined above, then the reset behavior described above is executed with probability p.restart
#   #p.restart.decay: if the parameters are reset, then p.restart = p.restart * p.restart.decay, this option is included to prevent the algorithm from restarting at every iteration, which may happen if the parameter limits are poorly chosen or on certain data sets
#   # p.restart.grow: if the parameters are on the boundary and a restart action is not performed, then p.restart = p.restart * p.restart.grow
#   # p.check.obj: evaluating the toal least squares objective function is very computationaly intensive, therefore we perform the evaluation with probability p.check.obj
#   # p.goto.best.par: with probably p.goto.best.par set beta to the best parameters observed so far
#   
#   #   if x.error wasn't initialized, initialize at 0
#   if(is.null(x.error)){
#     x.error <- x*0
#   }
#   #   cast ifixb as logical in case it was passed as an integer
#   ifixb <- as.logical(ifixb)
#   # get some useful measurements
#   n.row <- nrow(x); n.col.x <- ncol(x); n.col.y <- ncol(y)
#   #   call the tls objective function. when the objective function is called, the values for x.error that maximize the probability with beta fixed are stored in data.env$x.error. We will read those values and use them as the starting point for x.error for the minimizations.
#   Objective.Fxn(beta = beta,ifixb = ifixb,ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, Fxn = Fxn, x.error = x.error, tls.env = data.env, x.min = x.min, x.max = x.max  )
#   x.error <- data.env$x.error
#   
#   # minimize total squared error with levenberg-marquadt algorithm
#   nls.fit <- nls.lm(par = beta[!ifixb],lower = beta.lower[!ifixb] ,upper = beta.upper[!ifixb] ,fn = Objective.Fxn, ifixb = ifixb, ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, beta.lower = beta.lower, beta.upper = beta.upper, return.vector = TRUE, Fxn = Fxn, x.error = x.error, beta2 = beta, tls.env = data.env, x.max = x.max, x.min = x.min )
#   # parameters from lm fit are guaranteed to be the best so far, so store beta in best.par and the value in best.val
#   beta[!ifixb] <- nls.fit$par
#   best.par <- beta
#   best.int.par <- beta
#   best.val <- Objective.Fxn(beta = beta,ifixb = ifixb,ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, Fxn = Fxn, x.error = x.error, tls.env = data.env, x.min = x.min, x.max = x.max  )
#   x.error <- data.env$x.error
#   print(best.val)
#   print(best.par)
#   
#   j <- 1
#   k<-1
#   # initiate the learning rate based on user defined parameters
#   eta <- (beta.upper - beta.lower)/num.steps
#   if(any(!is.finite(eta))){
#     eta <- rep(x = 1/num.steps, times = sum(!ifixb))
#   }
#   reshuffle.counter <- 0
#   # shuffling a sequence of integers is a computationally efficient method for doing a balanced bootstrap
#   shuffle.seq <- resample(seq(n.row))
#   
#   while(j <= max.iter){
#     #     decrease learning rate at each iteration
#     eta <- eta * decay.constant
#     
#     #     increment counter for reshuffling, and reshuffle if needed 
#     reshuffle.counter <- reshuffle.counter + 1
#     if(reshuffle.counter > nrow(x)-1){
#       print(j)
#       shuffle.seq <- resample(seq(nrow(x)))
#       reshuffle.counter <- 1 
#     }
#     
#     # check the objective function with probability p.check.obj
#     if(runif(1) < p.check.obj){
#       #       calculate objective at current parameters
#       current.val <- Objective.Fxn(beta = beta,ifixb = ifixb,ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, Fxn = Fxn, x.error = x.error, tls.env = data.env, x.min = x.min, x.max = x.max  )
#       # retrieve x.errors that maximize probability density with beta fixed at current value
#       x.error <- data.env$x.error
#       # if this beta is better than the best parameters found so far, then update best.par
#       if(current.val < best.val){
#         nls.fit <- nls.lm(par = beta[!ifixb],lower = beta.lower[!ifixb] ,upper = beta.upper[!ifixb] ,fn = Objective.Fxn, ifixb = ifixb, ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, beta.lower = beta.lower, beta.upper = beta.upper, return.vector = TRUE, Fxn = Fxn, x.error = x.error, beta2 = beta, tls.env = data.env, x.max = x.max, x.min = x.min )
#         beta[!ifixb] <- nls.fit$par
#         best.par <- beta
#         best.val <- Objective.Fxn(beta = beta,ifixb = ifixb,ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, Fxn = Fxn, x.error = x.error, tls.env = data.env, x.min = x.min, x.max = x.max  )
#         x.error <- data.env$x.error
#         print(best.val)
#         print(best.par)
#       }
#     }
#     
#     # check if parameters are on the boundary with probability p.check.boundary 
#     if(runif(1) < p.check.boundary){
#       #       check for parameters within the margin of the boundary
#       parameter.on.boundary <- any(beta[!ifixb] < beta.lower[!ifixb] + (beta.upper[!ifixb] - beta.lower[!ifixb]) * parameter.boundary.margin) | any(beta[!ifixb] > beta.upper[!ifixb] - (beta.upper[!ifixb] - beta.lower[!ifixb]) * parameter.boundary.margin)
#       if(parameter.on.boundary){
#         #         if parameters are on the boundary, increase the current value of p.restart by a factor of p.restart.grow
#         p.restart <- p.restart * p.restart.grow
#         # with probability p.restart, move the current estimate off of the boundary        
#         if(runif(1) < p.restart){
#           #           if we restart, then decrease the probability that we will restart again
#           p.restart <- p.restart * p.restart.decay
#           #if parameters on boundary, then randomly reset parameters. The method for resetting parameters is randomly chosen between:
#           #1: randomly reset parameters
#           #2: set beta to previous best point on interior of parameter space
#           if(resample(c(TRUE,FALSE), 1)){
#             beta[!ifixb] <-  beta.lower[!ifixb] + (beta.upper[!ifixb] - beta.lower[!ifixb]) * runif(sum(!ifixb), min = parameter.boundary.margin, max = (1-parameter.boundary.margin))
#           }else{
#             beta <- best.int.par 
#           }
#           #           if the parameters were reset, then use lm for local minimization
#           x.error <- matrix(0, nrow = n.row, ncol = n.col.x)
#           Objective.Fxn(beta = beta,ifixb = ifixb,ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, Fxn = Fxn, x.error = x.error, tls.env = data.env, x.min = x.min, x.max = x.max  )
#           x.error <- data.env$x.error
#           nls.fit <- nls.lm(par = beta[!ifixb],lower = beta.lower[!ifixb] ,upper = beta.upper[!ifixb] ,fn = Objective.Fxn, ifixb = ifixb, ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, beta.lower = beta.lower, beta.upper = beta.upper, return.vector = TRUE, Fxn = Fxn, x.error = x.error, beta2 = beta, tls.env = data.env, x.max = x.max, x.min = x.min )
#           beta[!ifixb] <- nls.fit$par
#         }
#         
#       } 
#     }
#     # set beta to the best parameters with probability p.goto.best.par
#     if(runif(1) < p.goto.best.par){
#       beta <- best.par
#       current.val <- Objective.Fxn(beta = beta,ifixb = ifixb,ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, Fxn = Fxn, x.error = x.error, tls.env = data.env, x.min = x.min, x.max = x.max  )
#       x.error <- data.env$x.error
#     }
#     
#     if(runif(1) < p.rescale){
#       if(!is.null(rescale.vector)){
#         
#         if(all((beta[!ifixb] * rescale.vector[!ifixb] < beta.upper[!ifixb]) & (beta[!ifixb] * rescale.vector[!ifixb] > beta.lower[!ifixb]))){
#           beta[!ifixb] <- beta[!ifixb] * rescale.vector[!ifixb]
#           current.val <- Objective.Fxn(beta = beta,ifixb = ifixb,ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, Fxn = Fxn, x.error = x.error, tls.env = data.env, x.min = x.min, x.max = x.max  )
#           x.error <- data.env$x.error
#         }else{
#           if(all((beta[!ifixb]/rescale.vector[!ifixb] < beta.upper[!ifixb]) & (beta[!ifixb]/rescale.vector[!ifixb] > beta.lower[!ifixb]))){
#             beta[!ifixb] <- beta[!ifixb]/rescale.vector[!ifixb]
#             current.val <- Objective.Fxn(beta = beta,ifixb = ifixb,ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, Fxn = Fxn, x.error = x.error, tls.env = data.env, x.min = x.min, x.max = x.max  )
#             x.error <- data.env$x.error
#           } 
#         }
#         
#         
#       }
#     }
#     # calculate gradient
#     delta <- Grad(beta = beta,x = matrix(x[shuffle.seq[reshuffle.counter],] ,ncol = ncol(x)), y = matrix(y[shuffle.seq[reshuffle.counter],], ncol = ncol(y)), x.error = matrix(x.error[shuffle.seq[reshuffle.counter],], ncol = ncol(x)), x.weight = matrix(x.weight[shuffle.seq[reshuffle.counter],], ncol = ncol(x)), y.weight = matrix(y.weight[shuffle.seq[reshuffle.counter],], ncol = ncol(y)), beta.lower = beta.lower, beta.upper = beta.upper, Fxn = Fxn, ifixx = TRUE, ifixb = ifixb,tls.env = data.env, x.min = x.min, x.max = x.max )
#     
#     
#     # if the update would take beta outside of the parameter bounds, then decrease the current step size    
#     while(any(((beta - eta*delta) < beta.lower) | ((beta - eta*delta) > beta.upper))){
#       delta <- delta/10
#     }
#     
#     # update beta
#     beta <- beta - eta * delta
#     
#     # if beta is outside the parameter space, place it on the edge
#     #     beta[beta < beta.lower] <- beta.lower[beta < beta.lower]
#     #     beta[beta > beta.upper] <- beta.upper[beta > beta.upper]
#     
#     
#     j <- j + 1
#   }
#   return(list(value = best.val, par = best.par, fvec = Objective.Fxn(beta = beta,ifixb = ifixb,ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, Fxn = Fxn, x.error = x.error, tls.env = data.env, x.min = x.min, x.max = x.max, return.vector = TRUE )))    
# }
# S.G.D.C <- cmpfun(S.G.D)

S.G.D <- function(beta, x, y, x.error=NULL, x.weight, y.weight, x.min = NULL, x.max = NULL, beta.lower, beta.upper, ifixb, Objective.Fxn = T.L.S.C, solution.tolerance=.Machine$double.eps^.5, max.iter = 5000, num.steps = 500, Grad = T.L.S.Gradient.C, decay.constant=.9995, data.env = new.env(), Fxn, parameter.boundary.margin = .01, p.restart = 1, p.restart.decay = .5,   p.restart.grow = 1.1, p.check.obj = .05, p.goto.best.par = .005, p.check.boundary = .01, p.rescale = .005, rescale.vector = NULL, ...){
  #   Function: S.G.D (Stochastic Gradient Descent)
  # S.G.D is a implementation of a stochastic gradient descent algorithm, in principle it can be used to minimize arbitrary user specefied functions but it was designed with the total least squares problem in mind
  # Fxn: a user specefied function, here by default the function is a model of calcium buffering. the calling sequence for Fxn is Fxn(x, beta)
  
  # x: matrix of independent variables
  #   y: matrix of dependent variables
  # x.error: matrix of delta(x) values
  #   x.weight: matrix of "weights" for x errors. The correct weight is the inverse of the standard error.
  #   y.weight: matrix of "weights" for y errors. The correct weight is the inverse of the standard error.
  #   x.min: lower limit for x values. may be a vector or scalar. during the minimization procedure we will require that x + x.error > x.min
  #  x.max: analagous to x.min
  # beta.lower: lower limit for parameters
  # beta.upper: upper limit for parameters
  # ifixb: logical vector with ifixb[i] == TRUE indicating the the ith component of beta is fixed
  # Objective.Fxn: the function to be minized, by default a model of calcium buffering
  # solution.tolerance: used to check convergence and whether two numbers are "equal"
  # max.iter: maximum number of iterations to perform
  # num.steps: sets the initial learning rate eta = (beta.upper - beta.lower)/num.steps. Here eta is the learning rate. That is, the update rule for each iteration is beta = beta + eta * grad(Objective.Fxn)
  # Grad: function for computing the gradient of Objective.Fxn
  # decay.constant: controls the rate at which the learning rate eta "decays" for each iteration we have eta = eta * decay.constant
  # data.env: an R environment that is used to store some working data
  # Fxn: a user supplied function to be fit
  # parameter.boundary.margin: controls how close to the boundary the solution is allowed to be before restarting. If at any point we have abs(beta - beta.upper) < abs(beta.upper - beta.lower) * parameter.boundary.margin or abs(beta - beta.upper ) < abs(beta.upper - beta.lower) * parameter.boundary.margin, the beta is reset either randomly or the the best parameters so far
  # p.restart: if the parameters are "too close to the boundary" as defined above, then the reset behavior described above is executed with probability p.restart
  #p.restart.decay: if the parameters are reset, then p.restart = p.restart * p.restart.decay, this option is included to prevent the algorithm from restarting at every iteration, which may happen if the parameter limits are poorly chosen or on certain data sets
  # p.restart.grow: if the parameters are on the boundary and a restart action is not performed, then p.restart = p.restart * p.restart.grow
  # p.check.obj: evaluating the toal least squares objective function is very computationaly intensive, therefore we perform the evaluation with probability p.check.obj
  # p.goto.best.par: with probably p.goto.best.par set beta to the best parameters observed so far
  
  #   if x.error wasn't initialized, initialize at 0
  if(is.null(x.error)){
    x.error <- x*0
  }
  #   cast ifixb as logical in case it was passed as an integer
  ifixb.2 <- ifixb
  ifixb.2["kd.dye"] <- TRUE
  ifixb.2 <- as.logical(ifixb.2)
  ifixb <- as.logical(ifixb)
  # get some useful measurements
  n.row <- nrow(x); n.col.x <- ncol(x); n.col.y <- ncol(y)
  #   call the tls objective function. when the objective function is called, the values for x.error that maximize the probability with beta fixed are stored in data.env$x.error. We will read those values and use them as the starting point for x.error for the minimizations.
  Objective.Fxn(beta = beta,ifixb = ifixb,ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, Fxn = Fxn, x.error = x.error, tls.env = data.env, x.min = x.min, x.max = x.max  )
  x.error <- data.env$x.error
  nls.fit <- nls.lm(par = beta[!ifixb.2],lower = beta.lower[!ifixb.2] ,upper = beta.upper[!ifixb.2] ,fn = Objective.Fxn, ifixb = ifixb.2, ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, beta.lower = beta.lower, beta.upper = beta.upper, return.vector = TRUE, Fxn = Fxn, x.error = x.error, beta2 = beta, tls.env = data.env, x.max = x.max, x.min = x.min )
  beta[!ifixb.2] <- nls.fit$par
  Objective.Fxn(beta = beta,ifixb = ifixb,ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, Fxn = Fxn, x.error = x.error, tls.env = data.env, x.min = x.min, x.max = x.max  )
  x.error <- data.env$x.error
  
  # minimize total squared error with levenberg-marquadt algorithm
  nls.fit <- nls.lm(par = beta[!ifixb],lower = beta.lower[!ifixb] ,upper = beta.upper[!ifixb] ,fn = Objective.Fxn, ifixb = ifixb, ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, beta.lower = beta.lower, beta.upper = beta.upper, return.vector = TRUE, Fxn = Fxn, x.error = x.error, beta2 = beta, tls.env = data.env, x.max = x.max, x.min = x.min )
  # parameters from lm fit are guaranteed to be the best so far, so store beta in best.par and the value in best.val
  beta[!ifixb] <- nls.fit$par
  best.par <- beta
  best.int.par <- beta
  best.val <- Objective.Fxn(beta = beta,ifixb = ifixb,ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, Fxn = Fxn, x.error = x.error, tls.env = data.env, x.min = x.min, x.max = x.max  )
  x.error <- data.env$x.error
  print(best.val)
  print(best.par)
  
  j <- 1
  k<-1
  # initiate the learning rate based on user defined parameters
  eta <- (beta.upper - beta.lower)/num.steps
  if(any(!is.finite(eta))){
    eta <- rep(x = 1/num.steps, times = sum(!ifixb))
  }
  reshuffle.counter <- 0
  # shuffling a sequence of integers is a computationally efficient method for doing a balanced bootstrap
  shuffle.seq <- resample(seq(n.row))
  
  while(j <= max.iter){
    #     decrease learning rate at each iteration
    eta <- eta * decay.constant
    
    #     increment counter for reshuffling, and reshuffle if needed 
    reshuffle.counter <- reshuffle.counter + 1
    if(reshuffle.counter > nrow(x)-1){
      print(j)
      shuffle.seq <- resample(seq(nrow(x)))
      reshuffle.counter <- 1 
    }
    
    # check the objective function with probability p.check.obj
    if(runif(1) < p.check.obj){
      #       calculate objective at current parameters
      current.val <- Objective.Fxn(beta = beta,ifixb = ifixb,ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, Fxn = Fxn, x.error = x.error, tls.env = data.env, x.min = x.min, x.max = x.max  )
      # retrieve x.errors that maximize probability density with beta fixed at current value
      x.error <- data.env$x.error
      # if this beta is better than the best parameters found so far, then update best.par
      if(current.val < best.val){
        nls.fit <- nls.lm(par = beta[!ifixb],lower = beta.lower[!ifixb] ,upper = beta.upper[!ifixb] ,fn = Objective.Fxn, ifixb = ifixb, ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, beta.lower = beta.lower, beta.upper = beta.upper, return.vector = TRUE, Fxn = Fxn, x.error = x.error, beta2 = beta, tls.env = data.env, x.max = x.max, x.min = x.min )
        beta[!ifixb] <- nls.fit$par
        best.par <- beta
        best.val <- Objective.Fxn(beta = beta,ifixb = ifixb,ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, Fxn = Fxn, x.error = x.error, tls.env = data.env, x.min = x.min, x.max = x.max  )
        x.error <- data.env$x.error
        print(best.val)
        print(best.par)
      }
    }
    
    # check if parameters are on the boundary with probability p.check.boundary 
    if(runif(1) < p.check.boundary){
      #       check for parameters within the margin of the boundary
      parameter.on.boundary <- any(beta[!ifixb] < beta.lower[!ifixb] + (beta.upper[!ifixb] - beta.lower[!ifixb]) * parameter.boundary.margin) | any(beta[!ifixb] > beta.upper[!ifixb] - (beta.upper[!ifixb] - beta.lower[!ifixb]) * parameter.boundary.margin)
      if(parameter.on.boundary){
        #         if parameters are on the boundary, increase the current value of p.restart by a factor of p.restart.grow
        p.restart <- p.restart * p.restart.grow
        # with probability p.restart, move the current estimate off of the boundary        
        if(runif(1) < p.restart){
          #           if we restart, then decrease the probability that we will restart again
          p.restart <- p.restart * p.restart.decay
          #if parameters on boundary, then randomly reset parameters. The method for resetting parameters is randomly chosen between:
          #1: randomly reset parameters
          #2: set beta to previous best point on interior of parameter space
          if(resample(c(TRUE,FALSE), 1)){
            beta[!ifixb] <-  beta.lower[!ifixb] + (beta.upper[!ifixb] - beta.lower[!ifixb]) * runif(sum(!ifixb), min = parameter.boundary.margin, max = (1-parameter.boundary.margin))
          }else{
            beta <- best.int.par 
          }
          #           if the parameters were reset, then use lm for local minimization
          x.error <- matrix(0, nrow = n.row, ncol = n.col.x)
          Objective.Fxn(beta = beta,ifixb = ifixb,ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, Fxn = Fxn, x.error = x.error, tls.env = data.env, x.min = x.min, x.max = x.max  )
          x.error <- data.env$x.error
          nls.fit <- nls.lm(par = beta[!ifixb],lower = beta.lower[!ifixb] ,upper = beta.upper[!ifixb] ,fn = Objective.Fxn, ifixb = ifixb, ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, beta.lower = beta.lower, beta.upper = beta.upper, return.vector = TRUE, Fxn = Fxn, x.error = x.error, beta2 = beta, tls.env = data.env, x.max = x.max, x.min = x.min )
          beta[!ifixb] <- nls.fit$par
        }
        
      } 
    }
    # set beta to the best parameters with probability p.goto.best.par
    if(runif(1) < p.goto.best.par){
      beta <- best.par
      current.val <- Objective.Fxn(beta = beta,ifixb = ifixb,ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, Fxn = Fxn, x.error = x.error, tls.env = data.env, x.min = x.min, x.max = x.max  )
      x.error <- data.env$x.error
    }
    
    if(runif(1) < p.rescale){
      if(!is.null(rescale.vector)){
        
        if(all((beta[!ifixb] * rescale.vector[!ifixb] < beta.upper[!ifixb]) & (beta[!ifixb] * rescale.vector[!ifixb] > beta.lower[!ifixb]))){
          beta[!ifixb] <- beta[!ifixb] * rescale.vector[!ifixb]
          current.val <- Objective.Fxn(beta = beta,ifixb = ifixb,ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, Fxn = Fxn, x.error = x.error, tls.env = data.env, x.min = x.min, x.max = x.max  )
          x.error <- data.env$x.error
        }else{
          if(all((beta[!ifixb]/rescale.vector[!ifixb] < beta.upper[!ifixb]) & (beta[!ifixb]/rescale.vector[!ifixb] > beta.lower[!ifixb]))){
            beta[!ifixb] <- beta[!ifixb]/rescale.vector[!ifixb]
            current.val <- Objective.Fxn(beta = beta,ifixb = ifixb,ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, Fxn = Fxn, x.error = x.error, tls.env = data.env, x.min = x.min, x.max = x.max  )
            x.error <- data.env$x.error
          } 
        }
        
        
      }
    }
    # calculate gradient
    delta <- Grad(beta = beta,x = matrix(x[shuffle.seq[reshuffle.counter],] ,ncol = ncol(x)), y = matrix(y[shuffle.seq[reshuffle.counter],], ncol = ncol(y)), x.error = matrix(x.error[shuffle.seq[reshuffle.counter],], ncol = ncol(x)), x.weight = matrix(x.weight[shuffle.seq[reshuffle.counter],], ncol = ncol(x)), y.weight = matrix(y.weight[shuffle.seq[reshuffle.counter],], ncol = ncol(y)), beta.lower = beta.lower, beta.upper = beta.upper, Fxn = Fxn, ifixx = TRUE, ifixb = ifixb,tls.env = data.env, x.min = x.min, x.max = x.max )
    
    
    # if the update would take beta outside of the parameter bounds, then decrease the current step size    
    while(any(((beta - eta*delta) < beta.lower) | ((beta - eta*delta) > beta.upper))){
      delta <- delta/10
    }
    
    # update beta
    beta <- beta - eta * delta
    
    # if beta is outside the parameter space, place it on the edge
    #     beta[beta < beta.lower] <- beta.lower[beta < beta.lower]
    #     beta[beta > beta.upper] <- beta.upper[beta > beta.upper]
    
    
    j <- j + 1
  }
  
  
  Objective.Fxn(beta = beta,ifixb = ifixb,ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, Fxn = Fxn, x.error = x.error, tls.env = data.env, x.min = x.min, x.max = x.max  )
  x.error <- data.env$x.error
  nls.fit <- nls.lm(par = beta[!ifixb],lower = beta.lower[!ifixb] ,upper = beta.upper[!ifixb] ,fn = Objective.Fxn, ifixb = ifixb, ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, beta.lower = beta.lower, beta.upper = beta.upper, return.vector = TRUE, Fxn = Fxn, x.error = x.error, beta2 = beta, tls.env = data.env, x.max = x.max, x.min = x.min )
  
  
  return(list(value = best.val, par = best.par, fvec = Objective.Fxn(beta = beta,ifixb = ifixb,ifixx = FALSE, y = y, x = x, x.weight = x.weight, y.weight = y.weight, Fxn = Fxn, x.error = x.error, tls.env = data.env, x.min = x.min, x.max = x.max, return.vector = TRUE ), coef = try(summary(nls.fit)$coef)) )    
}
S.G.D.C <- cmpfun(S.G.D)

#############################################################################################################################
#############################################################################################################################

# End Section 2

#############################################################################################################################
#############################################################################################################################



#############################################################################################################################
#############################################################################################################################

# Begin Section 3: executable code

#############################################################################################################################
#############################################################################################################################

if(!interactive.file.chooser){
  if(substring(data.directory,String.Length(data.directory)) == "/"){
    data.directory <-  substring(data.directory,1, String.Length(data.directory)-1)
  }  
  fluorescence.filename <- file.path(data.directory,fluorescence.filename)
  fluorescence.se.filename <- file.path(data.directory, fluorescence.se.filename)
  delta.ca.total.filename <- file.path(data.directory, delta.ca.total.filename)
  delta.ca.total.se.filename <- file.path(data.directory, delta.ca.total.se.filename)
}


x <- as.matrix(read.csv(fluorescence.filename))
x.weight <- as.matrix(read.csv(fluorescence.se.filename))
y <- as.matrix(read.csv(delta.ca.total.filename))
y.weight <- as.matrix(read.csv(delta.ca.total.se.filename))

if( length(unique(c(nrow(x),nrow(y),nrow(x.weight),nrow(y.weight)))) > 1){
  stop("fluorescence, change in total calcium, and their respective measurement errors must have the same number of rows")
}


x.weight <- 1/x.weight
y.weight <- 1/y.weight
tls.data.env <- new.env()

x.min <- c(F.MIN, F.MIN)
x.max <- c(F.MAX, F.MAX)
ifixb <- parameter.is.fixed


if(!do.bootstrap.estimate){
  sgd.result <- S.G.D.C(beta = beta, x = x, y = y, x.weight = x.weight, y.weight = y.weight, x.min = x.min, x.max = x.max, beta.lower = beta.lower, beta.upper = beta.upper, ifixb = ifixb, max.iter = max.iterations, data.env = tls.data.env, Fxn = Delta.Ca.Total.Hat.C, parameter.boundary.margin = parameter.boundary.margin, p.restart = p.restart, p.restart.decay = p.restart.decay,   p.restart.grow = p.restart.grow, p.check.obj = p.check.obj, p.goto.best.par = p.goto.best.par, rescale.vector = rescale.vector, p.rescale = p.rescale)
  write(paste("Sum of total squared errors at final parameter estimates:",sgd.result$value,sep=" "), file = "results.txt")
  write("",file = "results.txt", append = TRUE)
  write("Paramter estimates:",file = "results.txt", append = TRUE)
  write("",file = "results.txt", append = TRUE)
#   write(parameter.names, file = "results.txt", append = TRUE)
  write(paste(parameter.names, sgd.result$par, sep = ": "), file = "results.txt", append = TRUE)
  write("",file = "results.txt", append = TRUE)

write.table(sgd.result$coef, file="results.txt", append=TRUE)
  

print(paste("Sum of total squared errors at final parameter estimates:",sgd.result$value,sep=" "))
print("Paramter estimates:")
print(paste(parameter.names, sgd.result$par, sep = ": "))
print("Standard errors were not estimated.")
print(paste("A copy of these results can be found in", file.path(getwd(), "results.txt") ))
}else{
  boot.estimates <- NULL
  
  cl <- makeCluster(n.threads)
  registerDoParallel(cl)
  
  boot.estimates <- NULL
  boot.estimates <- foreach(i = seq(bootstrap.replicates), .packages=c("minpack.lm","mosaic"), .combine = rbind ) %dopar% {
    boot.samples <- ((seq(nrow(x)) + resample(size = nrow(x),seq(-floor(nrow(x)/10), floor(nrow(x)/10)))) %% nrow(x)) + 1
    tls.data.env <- new.env()
    tls.data.env$x.error <- matrix(0,nrow = nrow(x), ncol = ncol(x))
    sgd.result <- S.G.D.C(beta = beta, x = x[boot.samples,], y = matrix(y[boot.samples,], nrow = nrow(y), ncol = ncol(y)), x.weight = x.weight[boot.samples,], y.weight = matrix(y.weight[boot.samples,], nrow = nrow(y), ncol = ncol(y)), x.min = x.min, x.max = x.max, beta.lower = beta.lower, beta.upper = beta.upper, ifixb = ifixb, max.iter = max.iterations, data.env = tls.data.env, Fxn = Delta.Ca.Total.Hat.C, parameter.boundary.margin = parameter.boundary.margin, p.restart = p.restart, p.restart.decay = p.restart.decay,   p.restart.grow = p.restart.grow, p.check.obj = p.check.obj, p.goto.best.par = p.goto.best.par, rescale.vector = rescale.vector, p.rescale = p.rescale)
    
    if(replace.on.boundary){
      while(sum((abs(sgd.result$par[!ifixb] - beta.lower[!ifixb]) < boundary.margin[!ifixb]) + (abs(sgd.result$par[!ifixb] - beta.upper[!ifixb]) < boundary.margin[!ifixb])) > 0){
        boot.samples <- ((seq(nrow(x)) + resample(size = nrow(x),seq(-floor(nrow(x)/10), floor(nrow(x)/10)))) %% nrow(x)) + 1
        tls.data.env <- new.env()
        tls.data.env$x.error <- matrix(0,nrow = nrow(x), ncol = ncol(x))
        sgd.result <- S.G.D.C(beta = beta, x = x[boot.samples,], y = matrix(y[boot.samples,], nrow = nrow(y), ncol = ncol(y)), x.weight = x.weight[boot.samples,], y.weight = matrix(y.weight[boot.samples,], nrow = nrow(y), ncol = ncol(y)), x.min = x.min, x.max = x.max, beta.lower = beta.lower, beta.upper = beta.upper, ifixb = ifixb, max.iter = max.iterations, data.env = tls.data.env, Fxn = Delta.Ca.Total.Hat.C, parameter.boundary.margin = parameter.boundary.margin, p.restart = p.restart, p.restart.decay = p.restart.decay,   p.restart.grow = p.restart.grow, p.check.obj = p.check.obj, p.goto.best.par = p.goto.best.par, rescale.vector = rescale.vector, p.rescale = p.rescale)
      }
      
      
    }
    if(sgd.result$par[4] > sgd.result$par[6]){
      beta.4 <- sgd.result$par[4]
      sgd.result$par[4] <- sgd.result$par[6]
      sgd.result$par[6] <- beta.4
      
      beta.5 <- sgd.result$par[5]
      sgd.result$par[5] <- sgd.result$par[7]
      sgd.result$par[7] <- beta.5
    }
    c(sgd.result$par,TSE=sgd.result$value)
  }
  
  stopImplicitCluster()
  stopCluster(cl)
  
  write.csv(boot.estimates, "results.txt", row.names = FALSE)
  print(paste("Finished. Results are contained in", file.path(getwd(), "results.txt") ))
}


