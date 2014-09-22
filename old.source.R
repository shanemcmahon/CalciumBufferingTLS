Delta.Ca.Total.Hat <- function (par.1, data.1=data.1){
  
  calcium.final   <- ((data.1$f.f - (1/par.1["rf"]))/(1 - data.1$f.f)) * par.1["kd.3"]
  calcium.initial <- ((data.1$f.i - (1/par.1["rf"]))/(1 - data.1$f.i)) * par.1["kd.3"]
  delta.ca.total  <- (1+par.1[["kappa.nonsaturable"]])*(calcium.final - calcium.initial) - par.1["kd.1"] * par.1["bt.1"]/(
    par.1["kd.1"] + calcium.final) + par.1["kd.1"] * par.1["bt.1"]/(par.1["kd.1"] + calcium.initial
    ) - par.1["kd.2"] * par.1["bt.2"]/(par.1["kd.2"] + calcium.final) + par.1["kd.2"] * par.1["bt.2"]/(
      par.1["kd.2"] + calcium.initial) - par.1["kd.3"] * par.1["bt.3"]/(par.1["kd.3"] + calcium.final)+
    par.1["kd.3"] * par.1["bt.3"]/(par.1["kd.3"] + calcium.initial)
  
  return(delta.ca.total)
}
Delta.Ca.Total.Hat.compiled <- cmpfun(Delta.Ca.Total.Hat)

TLS.Distance <- function(fluorescence=c(.5, .5), par.1, data.1, return.vector=FALSE){
  
  delta.ca.total.hat <- Delta.Ca.Total.Hat(data.1=data.frame(f.i.=fluorescence[1],f.f=fluorescence[2]), par.1=par.1)  
  if (is.null(par.1["f.noise"])){
    stop("par.1[[\"f.noise\"]] can not be of type NULL ")
  }
  if (is.null(par.1["ca.noise"])){
    stop("par.1[[\"ca.noise\"]] can not be of type NULL ")
  }
  
  
  out <- ((fluorescence[1] - data.1["f.i"])^2)/(par.1["f.noise"]^2) + ((fluorescence[2] - data.1["f.f"])^2)/(par.1["f.noise"]^2) + ((delta.ca.total.hat-data.1["delta.ca.total"] * par.1["ev"])^2)/((par.1["ca.noise"])^2)
  
  if (any(fluorescence>1)){
    out <- Inf
  }
  attributes(out)$f.i <- fluorescence[1]
  attributes(out)$f.i <- fluorescence[2]  
  return(out^.5)
}
TLS.Distance.Compiled <- cmpfun(TLS.Distance)



TLS <- function(par.1, data.1, par.2=NULL, my.lower=-Inf, my.upper=Inf, par.1.names=NULL,par.2.names=NULL, return.vector=FALSE, Distance.Function = TLS.Distance.Compiled,trim=FALSE){
  if(is.null(par.1.names)){
    par.1.names <- attributes(par.1)$names
  }
  if(is.null(par.2.names)){
    par.2.names <-  attributes(par.2)$names
  }
  p.3 <- c(par.1,par.2)
  p.3 <- unlist(p.3)
  attributes(p.3)$names <- c(par.1.names, par.2.names)
  
  
  
  values <- array(0, length(data.1[, 1]))
  #   attributes(values)$f.i <- values
  #   attributes(values)$f.f <- values
  values.f.i <- values
  values.f.f <- values
  system.time({
    for(data.index in 1:length(data.1[, 1])){
      point <- data.1[data.index, ]
      fit <- optim(par=c(point["f.i"],point["f.f"]),fn=Distance.Function, par.1=p.3,data.1=point)
      # attributes(values)$f.i[data.index]<-fit$par[1]
      values.f.i[data.index] <- fit$par[1]
      # attributes(values)$f.f[data.index]<-fit$par[2]
      values.f.f[data.index]<-fit$par[2]
      values[data.index]<-fit$value
    }
  })
  if (any(par.1 < my.lower)){
    values <- values * (1 + sum(abs((my.lower - par.1)[which(par.1 < my.lower)]) ) ) 
  }
  if (any(par.1 > my.upper)){
    values <- values * (1 + sum(abs((par.1 - my.upper)[which(par.1 > my.upper)]) ) )
  }
  
  if (trim){
    values[values^2>median(values^2)+16.27]<-median(values)
  }
  if(length(values) > 10){
    print(par.1)
    print(sum(values^2))
  }
  if(return.vector==TRUE){
    return(values)
  }
  return(sum(values^2))
}
TLS.Compiled <- cmpfun(TLS)

TLS.Gradient <- function(par.1, data.1, par.2=NULL, my.lower=-Inf, my.upper=Inf, par.1.names=NULL, par.2.names=NULL, return.vector=TRUE, Distance.Function = TLS.Distance.Compiled, return.matrix=FALSE){
  gradient <- NULL
  
  for(par.index in 1:length(par.1)){
    delta <- rep(0, length(par.1))
    delta[[par.index]] <- c(1E-6, par.1[par.index] * 1E-6)[which.max( c(1E-6,abs((par.1[par.index])*1E-6)))]
    delta <- unlist(delta)
    gradient <- cbind(gradient, (TLS.Compiled(par.1=(par.1 + delta), data.1=data.1, par.2=par.2, my.lower=my.lower, my.upper = my.upper, par.1.names=par.1.names, par.2.names=par.2.names,  return.vector=TRUE, Distance.Function=Distance.Function) - TLS.Compiled(par.1=(par.1-delta), data.1=data.1, par.2=par.2, my.lower=my.lower, my.upper=my.upper, par.1.names=par.1.names, par.2.names=par.2.names,  return.vector=TRUE, Distance.Function=Distance.Function))/(2*delta[par.index]))
  }
  if (return.vector){
    return(colSums(gradient))
  }
  
  if (return.matrix){
    return(gradient)
  }
  return(colSums(gradient))
}
TLS.Gradient.Compiled <- cmpfun(TLS.Gradient)

SGD.6 <- function(par.1, data.1, par.2, 
                  my.lower=lower.1, my.upper=upper1, par.1.names, par.2.names,
                  Model=TLS.Compiled, eta=.1, solution.tolerance=1E-12, p.1 = 17,  p.2 = 
                    11, p.3 = 13, p.5 = 29, max.iter = 20000, num.steps = 100, Grad=
                    TLS.Gradient.Compiled, verbose=0, decay.constant=.99){
  counter.3 <- 1
  pars <- par.1
  counter.2 <- 1
  counter <- 1
  best.val <- (Model(par.1=par.1, data.1=data.1, par.2=
                       par.2, my.lower=my.lower, my.upper=my.upper, par.1.names=par.1.names, 
                     par.2.names=par.2.names))
  new.val <- best.val
  best.par <- par.1
  print(c(best.par, best.val))
  best.vals <- array(0, p.1)
  j <- 1
  k<-1
  eta <- (my.upper - my.lower)/num.steps  
  
  while(j <= max.iter * length(par.1)){
    eta <- eta * decay.constant
    if(counter.2 > 20){
      par.1 <- best.par
      counter.2 <- 0
    }
    point <- resample(x=data.1, size=3)
    
    delta <- (Grad(par.1=(par.1), data.1=point, par.2=
                     par.2, my.lower=my.lower, my.upper=my.upper, par.1.names=par.1.names, 
                   par.2.names=par.2.names, return.vector=FALSE))
    
    asdf <- runif(n=1, min=0, max=1)
    
    while(any(((par.1 - eta*delta) < my.lower) | ((par.1 - eta*delta) >
                                                    my.upper))){
      delta <- delta/10
    }
    
    counter.2 <- counter.2  +  1
    par.1 <- par.1 - eta * delta
    par.1[par.1 < my.lower] <- my.lower[par.1 < my.lower]
    par.1[par.1 > my.upper] <- my.upper[par.1 > my.upper]
    counter.2 <- 1
    
    pars <- rbind(pars, par.1)
    print(dim(pars))
    if(verbose>=1){
      print(eta*delta)
    }
    if(j %% p.2 * length(par.1)==0){
      counter.3 <- counter.3  +  1
      par.1 <- (apply(X=pars, MARGIN=2, median))
      print(par.1)
      pars <- par.1
      new.val <- (Model(par.1=(par.1), data.1=data.1, par.2
                        =par.2, my.lower=my.lower, my.upper=my.upper, par.1.names=par.1.names, 
                        par.2.names=par.2.names))
      
      if(new.val < best.val){
        fit <- nls.lm(fn=TLS.Compiled, par=par.1, data.1=data.1,
                      par.2=par.2, my.lower=my.lower, lower=my.lower, 
                      my.upper=my.upper, upper=my.upper, rf=50, par.1.names=par.1.names, par.2.names
                      =par.2.names, return.vector=TRUE, control=nls.lm.control(ptol=1E-4, 
                                                                               ftol=1E-4))
        
        par.1 <- fit$par  
        new.val <- (Model(par.1=(par.1), data.1=data.1, par.2
                          =par.2, my.lower=my.lower, my.upper=my.upper, par.1.names=par.1.names,
                          par.2.names=par.2.names))
        best.par <- fit$par
        best.val <- new.val
        counter <- 1
        print(best.par)
        counter.3 <- 1
      }
      best.vals[k %% p.1 + 1] <- best.val
      if  (!(best.vals[p.1]==0)){
        if (sd(best.vals) < best.val * solution.tolerance){
          break
        }
      }
      if(!(best.vals[p.1]==0)){
        if(all(best.vals==min(best.vals)) ){
          break
        }
      }
      if(counter.3 > p.5){
        break
      }
      counter <- counter  +  1
      if(counter > 100){
        break
      }
      k <- k + 1
      if(k %% p.3 * length(par.1)){
        par.1 <- best.par
      }
      print(c(best.par, best.val))
    }    
    j <- j + 1
  }
  if(j >= max.iter){
    print(warning("max.iter exceeded"))
  }
  return(list(value = best.val, par = best.par, fvec = (Model(par.1=
                                                                (best.par), data.1=data.1, par.2=par.2, my.lower=my.lower, 
                                                              my.upper=my.upper, par.1.names=par.1.names, par.2.names=par.2.names, Model=Model, 
                                                              return.vector=TRUE))))    
}
SGD.6.Compiled <- cmpfun(SGD.6)