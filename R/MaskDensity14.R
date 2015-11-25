
########################################
##########################################
#    Version 10 is the same as Version 9
#    Version 6.3: modify MaskDensity6.2
#      Version 6.2 does not work for categorical variable
#       neworder in function"findOrder_Rmask" is also changed
#
#     let a and b decide by data agency as well as a and b determined by alpha
#
#  This version is modified based on MaskDensity11.R
#  The method of estimating the mass function of a categorical variable
#  is changed.
#############################################
###########################################


# In this version of MaskDensity, the value of CORmax is reported.
# The value will be used to show the performance of  the approximant of density function.


encriptNoise <- function(noise, a, b, maxorder, levels, EPS, noisefile)
{
  lvls <- levels
  save(a, b, maxorder, lvls, noise, EPS, file=noisefile)
}

calc_muX <- function (muY, a, b) 
{
  muX <- 0 * muY
  muX[1] <- 1
  for (j in 1:(length(muY) - 1)) {
    k <- 0:j
    muX[j + 1] <- sum(choose(j, k) * 2^k * muY[k + 1] *
                        (-1)^(j - k) * (a + b)^(j - k))/(b - a)^j
  }
  return(muX)
}

P_Rmask <- function (x, k) 
{
  SUM <- 0 * x
  for (i in 0:floor(k/2)) SUM <- SUM + (-1)^i/2^k *
    exp(lgamma(2 * k - 2 * i + 1) - lgamma(i + 1) -
          lgamma(k - i + 1) - lgamma(k - 2 * i + 1)) * x^(k - 2 * i)
  return(SUM)
}

lambda_Rmask <- function (muX, k) 
{
  i <- 0:floor(k/2)
  val <- (2 * k + 1)/2 * sum((-1)^i/2^k * exp(lgamma(2 *
                                                       k - 2 * i + 1) - lgamma(i + 1) - lgamma(k - i + 1) -
                                                lgamma(k - 2 * i + 1)) * muX[k - 2 * i + 1])
  return(val)
}

fY_Rmask <- function(y, muY, a, b) 
{
  x <- (2 * y - (a + b))/(b - a)
  muX <- calc_muX(muY, a, b)
  SUM <- 0 * y
  for (k in 0:(length(muY) - 1)) SUM <- SUM + lambda_Rmask(muX, k) * P_Rmask(x, k)
  SUM[(y < a) || (y > b)] <- 0
  val <- SUM * 2/(b - a)
  val[is.na(val)] <- 0
  return(val)
}


moments <- function (y, order = 8) 
{
  val <- 1
  if (order > 0) for (k in 1:order)
    val <- c(val, mean(y^k))
  return(val)
}


density_Rmask <- function (moments, a, b, n = 512)
{
  span <- 25
  fYseq <- yseq<- seq(a, b, length.out = n)
  for(i in 1:n) {
    jseq <- i + (-span):span
    if(i <= span) jseq <- 1:i
    if(i > n-span) jseq <- i:n
    fYseq[i] <- length(jseq)/(2*span+1) * mean(fY_Rmask(yseq[jseq], moments, a, b))
  }
  scale=sum(fYseq)
  negative <- (1:n)[fYseq <= 0]
  fYseq[negative] <- 0
  output <- list(x = yseq, y = fYseq, call=sys.call())
  #class(output) <- "density"
  return(output)
}


findOrder_Rmask <- function(ystar, noise, a, b, maxorder=100, EPS)
{
  M <- mean(noise)
  ymoments <- moments(ystar/M, order=maxorder) /
    moments(noise/M, order=maxorder)
  
  #############################################
  #   this part is changed
  #
  #neworder<-NULL                    # the begining of changes
  #for(i in 1:maxorder){
  #  if(is.finite(ymoments[i]))
  #     neworder<-i
  #    }
  # maxorder<-neworder-1    # end of changes
  
  #############################
  # this code is new
  ################
  
  neworder<-maxorder                    # the begining of changes
  for(i in 1:maxorder+1){
    if(is.infinite(ymoments[i]))
      neworder<-i-1
    break
  }
  maxorder<-neworder    # end of changes
  
  ###################
  
  C <- sample(noise, size=length(ystar), replace=TRUE)
  
  CORsamp <- NULL
  
  CORmax <- -Inf
  
  for(m in 1:maxorder) {
    cat("Trying", m, "moments ../\n")
    Provost <- density_Rmask(ymoments[1:m], a, b)
    
    ############
    #####
    #
    ###  modify here
    #
    ##################
    
    # if (sum(Provost$y)==0) 
    if (sum(is.nan(Provost$y)) >=1) 
    {
      cat("NaN in probabilty vector../\n")
      break
    } 
    if (sum(Provost$y)==0)  
    {
      cat("NA in probabilty vector../\n")
      break
    } 
    else 
    {
      prob <- Provost$y/sum(Provost$y)
      y1 <- sample(Provost$x, size = length(ystar), prob=prob, replace=TRUE)
      y1star <- y1*C
      COR <- round(10000*cor(sort(ystar), sort(y1star)))/10000
      if(COR > (1+EPS)*CORmax) {
        CORmax <- COR
        mOPT <- m
      }
    }
    CORsamp <- c(CORsamp, COR)
    
    if(COR < 1 - 10*(1-CORmax)) break
    
  }
  
  #cat("unmask:", mOPT, "moments used.\n")
  cat("unmask:", mOPT, "moments used.\n", CORmax, "correlation. \n")
  
  return(mOPT)
}





createNoise <- function (n, mean = rep(NA, 5), sd = rep(1, length(mean)),prob = rep(1, length(mean))/length(mean))  
{
  k <- length(mean)
  
  if (is.na(sum(mean))) 
  {    
    best <- rnorm(k, 500, 200)
    MAX <- min(dist(best))
    
    for (i in 1:1000) {
      mean <- rnorm(k, 500, 200)
      if ((min(dist(mean)) > MAX) && (min(mean) > 50)) {
        MAX <- min(dist(mean))
        best <- mean
      }
    }
    mean <- best
  }
  
  C <- rnorm(n)
  
  grp <- sample(1:length(mean), n, prob = prob, replace = TRUE)
  
  for (i in 1:length(mean)) 
    C[grp == i] = mean[i] + sd[i] * C[grp == i]  
  C[C < 0] <- -C[C < 0]
  
  return(C)
}



######## the following mask function is modified #####

mask <- function (y, noisefile, noise = createNoise(length(y)), a1=min(y), b1=max(y),
                  maxorder = 100, EPS = 1e-06) 
{
  n <- length(y)
  
  if(length(noise) != n) stop("'y' and 'noise' lengths differ")
  
  if(sum(noise < 0) > 0) stop("Negative noise not allowed")
  
  lvls <- NULL
  
  FACTOR <- is.factor(y)
  
  if (FACTOR) 
  {
    k <- length(levels(y))
    lvls <- levels(y)
    y1 <- rep(NA, length(y))
    for (i in 1:k) y1[y == lvls[i]] <- i
    y <- y1
  }
  
  if (sum(noise < 0) > 0) 
    stop("Negative noise not allowed")
  
  ndens <- density(noise)
  prob <- ndens$y/sum(ndens$y)
  ystar = y * noise
  
  C2 <- sample(noise, size = 10*n, replace = TRUE)
  
  if (FACTOR) 
  {
    a <- 0
    b <- k + 1
  }  
  else {
    a<- a1
    b<- b1
    k <- 0
  }
  
  encriptNoise(C2, a, b, maxorder, lvls, EPS, noisefile)
  
  return(list(ystar = ystar, noisefile = noisefile))
}

####the following unmask_Rmask is modified  #######

#########
#   unmask_Rmask1
#########


unmask <- function (ystar, noisefile, noise)
{
  load(noisefile)
  outMeanOfNoise <- mean(noise)
  outMeanOfSquaredNoise <- mean(noise^2)
  unmask_Rmask1 <- function(ystar)  #modified
  {
    if(!exists("noise")) stop("Problem with 'noisefile'")
    if(!exists("a")) stop("Problem with 'noisefile'")
    if(!exists("b")) stop("Problem with 'noisefile'")
    if(!exists("maxorder")) stop("Problem with 'noisefile'")
    if(!exists("EPS")) stop("Problem with 'noisefile'")
    
    order <- findOrder_Rmask(ystar, noise, a, b, maxorder = maxorder, EPS = EPS)
    M <- mean(noise)
    ymoments <- moments(ystar/M, order=order) /moments(noise/M, order=order)
    
    return(list(a=a, b=b, levels=lvls, ymoments=ymoments))
  }
  
  
  ##########
  #  unmask_Rmask2
  #########
  
  unmask_Rmask2 <- function(ystar, alpha)  #modified
  {
    if(!exists("noise")) stop("Problem with 'noisefile'")
    if(!exists("a")) stop("Problem with 'noisefile'")
    if(!exists("b")) stop("Problem with 'noisefile'")
    if(!exists("maxorder")) stop("Problem with 'noisefile'")
    if(!exists("EPS")) stop("Problem with 'noisefile'")
    
    if (is.null(lvls)){
      a1<- mean(ystar)/mean(noise)-sqrt(1/alpha)*sqrt(mean(ystar^2)/mean(noise^2)-(mean(ystar)/mean(noise))^2)
      b1<- mean(ystar)/mean(noise)+sqrt(1/alpha)*sqrt(mean(ystar^2)/mean(noise^2)-(mean(ystar)/mean(noise))^2)
      
      a<-max(a,a1)
      b<-min(b,b1)
    }
    
    order <- findOrder_Rmask(ystar, noise, a, b, maxorder = maxorder, EPS = EPS)
    M <- mean(noise)
    ymoments <- moments(ystar/M, order=order) /moments(noise/M, order=order)
    
    return(list(a=a, b=b, levels=lvls, ymoments=ymoments))
  }
  ######  original code  ####   
  #    info <- unmask_Rmask(ystar, noisefile, alpha)   ##### changed
  #    size <- length(ystar)
  #    a <- info$a
  #    b <- info$b
  #    lvls <- info$levels
  #    ymoments <- info$ymoments
  #    dens <- density_Rmask(ymoments, a, b)
  ######### end ##########
  
  #DensM<-NULL        # modified 
  DensA<-NULL
  DensB<- NULL
  #DensLvls<-NULL
  
  sampDens<-NULL   # use thid to instore  the sample from the unmasked density
  
  corMatrix<-matrix(0,6,2)
  
  ################
  # Add a function for discrete case. For discrete variables, their mass functions
  # are estimated from a different way.
  ##############
  
  k<-length(lvls)
  if (k !=0){
    Ma<-seq(1,k, by=1)
    MA<-Ma
    Mb<-Ma
    for(i in 1:(k-1)){
      Mb<-Mb*Ma
      MA<-cbind(MA,Mb)
    }
    
    MA<-t(MA)
    
    Mnoise<-diag(k)
    Mnoise[1,1]<-mean(noise)
    noisePower<-noise
    for( i in 1:(k-1)){
      noisePower<-noisePower*noise
      Mnoise[i+1,i+1]<-mean(noisePower)
    }
    MstarY<-mean(ystar)
    ystarPower<-ystar
    for(i in 1:(k-1)){
      ystarPower<-ystarPower*ystar
      MstarY<-c(MstarY, mean(ystarPower))
    }
    
    Poutput<-solve(MA)%*%solve(Mnoise)%*%MstarY
    
    if(min(Poutput) > 0){
      sizeOut<-length(ystar)
      y1<-sample(Ma, size=sizeOut, replace=TRUE, prob=Poutput)
      return(list(unmaskedVariable = y1, prob=Poutput, meanOfNoise = outMeanOfNoise, meanOfSquaredNoise = outMeanOfSquaredNoise))
    }else{ # if matrix is singluar
      ############ end of modification
      
      ###############
      #  the method used to estimate the mass function of a categorical variable is changed.
      # the old method is kept below
      #############
      #  (from here)
      #############
      warning("Method of moments failed, using k-means instead",immediate. =T)
      
      info <- unmask_Rmask1(ystar)   ##### changed
      size <- length(ystar)
      a <- info$a
      b <- info$b
      ymoments <- info$ymoments
      
      dens <- density_Rmask(ymoments, a, b)
      densK<-dens     #  this density information is for categorical case
      
      y1 <- sample(densK$x, size = size, replace = TRUE, prob = densK$y/sum(densK$y))
      centers <- 1:k
      OK <- FALSE
      rmv <- NULL
      repeat {
        try({ out <- kmeans(y1, centers=centers[!centers %in% rmv])$cluster;
              OK <- TRUE })
        if(OK) {
          break
        }
        y2 <- round(y1)
        for(i in centers) {
          if(sum(y2 == i) == 0) {
            rmv <- c(rmv, i)
          }
        }
      }
      
      out <- as.factor(lvls[out])      
      levels(out) <- lvls
      Poutput<-summary(out)/length(out)

      #   return(out)
      return(list(unmaskedVariable = out, prob=Poutput, meanOfNoise = outMeanOfNoise, meanOfSquaredNoise = outMeanOfSquaredNoise))
    }
  } 
  else{ # continuous case
    ###########
    #  determine by a and b
    #######
    
    info <- unmask_Rmask1(ystar)   ##### changed
    size <- length(ystar)
    a <- info$a
    b <- info$b
    lvls <- info$levels
    
    ymoments <- info$ymoments
    dens <- density_Rmask(ymoments, a, b)
    yy1 <- sample(dens$x, size = size, replace = TRUE, prob = dens$y/sum(dens$y))
    c1<-sample(noise, size=length(ystar), replace=TRUE)
    yy1star<-yy1*c1
    
    corMatrix[1,1] <- round(10000*cor(sort(ystar), sort(yy1star)))/10000
    
    #DensM<-cbind(DensM, dens)     # modeified
    DensA<-c(DensA,a)
    DensB<-c(DensB,b)
    
    sampDens<-cbind(sampDens, yy1)
    
    #######
    # determined by alpha
    ###########
    
    size <- length(ystar)
    corMatrix[1,2]<-1
    
    for (i in 1:5){
      alpha=0.01*i
      corMatrix[i+1,2]<-i+1
      info <- unmask_Rmask2(ystar, alpha)   ##### changed
      
      a <- info$a
      b <- info$b
      
      lvls <- info$levels
      ymoments <- info$ymoments
      
      dens <- density_Rmask(ymoments, a, b)
      yy1 <- sample(dens$x, size = size, replace = TRUE, prob = dens$y/sum(dens$y))
      
      c1<-sample(noise, size=length(ystar), replace=TRUE)
      
      yy1star<-yy1*c1
      
      corMatrix[i,1] <- round(10000*cor(sort(ystar), sort(yy1star)))/10000
      
      #DensM<-cbind(DensM, dens)     # modeified
      DensA<-c(DensA,a)
      DensB<-c(DensB,b)
      
      sampDens<-cbind(sampDens, yy1)
    }
    
    corMatrix<-corMatrix[order(corMatrix[,1]),]  #increased based on the values of cor
    
    selected<-corMatrix[6,2]
    size <- length(ystar)
    a <- DensA[selected]
    b <- DensB[selected]
    
    y1<- sampDens[, selected]
  
    
    return(list(unmaskedVariable = y1,  meanOfNoise = outMeanOfNoise, meanOfSquaredNoise = outMeanOfSquaredNoise))
  }
}








