################################################################################
#                                                                              #
# Systems analysis R library                                                   #
# ==========================                                                   #
#                                                                              #
# Peter Reichert, Eawag (basic functions)                                      #
# Dario Del Giudice, Eawag - Carnegie (SIP. Bias, multipliers, log-sinh etc)   #
# last modified: Feb 2016                                                      # 
#                                                                              #
################################################################################

# Overview of functions provided by this library
# (see individual descriptions of functions for more details):
#
# Data handling
# -------------
# sysanal.decode                  Decode a result code used by plotting funct.
#
# Plotting:
# ---------
# sysanal.plot.res                Plot results of a simulation
# sysanal.plot.vars               Plot several sets of results
# sysanal.contourpdf              Contour a 2d probability density function
# sysanal.plot.chains             Plot components of a Markov chain
# sysanal.plot.margs              Plot marginals of a sample
#
# Sampling:
# ---------
# sysanal.randsamp                Sample from a mv. normal or lognormal dist.
# sysanal.markovchain.metropolis  Calculate a Metropolis Markov chain of a pdf
# sysanal.jump.dist               Calculate a jump dist. for a Metropolis MC
#
# Model formulation:
# ------------------
# sysanal.calcpdf                 Calculate density of univariate distributions
# sysanal.calcpdf_mv              Calculate density of multivariate distrib.
# sysanal.loglikeli               Calculate log likelihood of simple model
# sysanal.logposterior            Calculate log posterior of simple model
#
# Sensitivity analysis:
# ---------------------
# sysanal.sens.loc                Calculate local sensitivities of det. model
# sysanal.sens.var.samp           Calculate variance-based sens. of det. model
#
# Identifiability analysis:
# -------------------------
# sysanal.comb                    Calculate combinations of subsets of indices
# sysanal.collind                 Calculate collinearity index
# sysanal.ident                   Calculate identifiability measures
#
# Output transformation:
# ----------------------
# sysanal.boxcox                  Box-Cox transformation
# sysanal.boxcox.deriv            derivative of Box-Cox transformation
# sysanal.boxcox.inv              Inverse of Box-Cox transformation
#
# Multivariate regression:
# ------------------------
# sysanal.gnls                    Generalized nonlinear least squares regression
# sysanal.gnls.diag               Diagnostics for gnls
# sysanal.gnls.test               Statistical tests for gnls
# sysanal.gnls.predict            Prediction and confidende intervals for gnls
# sysanal.confint                 Estimation of confidence intervals
#
# Residual diagnostics:
# ---------------------
# sysanal.resid.diag              Residual diagnostics plots
# sysanal.resid.diag.boxcox       Residual diagnostics plots with Box-Cox trans.
#
# Miscellaneous:
# --------------
# sysanal.hessian                 Calculate the Hessian of a real-valued funct.
#
################################################################################

# sysanal.decode
# ==============

# purpose:
# decode a simple model layout code into variable names and value of input 
# variable
# (the function assumes that the code consists of the name of the output 
# variable and the value of the input variable separated by an underline 
# character; the variable name may still contain underline characters)

# arguments:
# rescode:     vector of result codes or vector of results with codes as 
#              names of the elements (the function assumes that the code 
#              consists of the name of the output variable and the value of the 
#              input variable separated by an underline character; the variable 
#              name may still contain underline characters)

# output:
# data frame with variables:
# var:         vector of character strings of variables
# val:         vector of corresponding values of input variable
# the row names of the data frame contain the original codes

sysanal.decode <- function(L)
{
   if ( is.numeric(L) ) codes <- names(L)
   else                 codes <- L
   var <- character(length(codes))
   val <- numeric(length(codes))
   for ( i in 1:length(codes))
   {
      s <- strsplit(as.character(codes[i]),split="_")
      n <- length(s[[1]])
      if ( n < 2 ) stop("error in sysanal.decode: illegal result code")
      var[i] <- paste(s[[1]][1:(n-1)],collapse="_")
      val[i] <- as.numeric(s[[1]][n])
   }
   L.decoded <- data.frame(var=var,val=val)
   L.decoded$var <- as.character(L.decoded$var) #access the variable var
   rownames(L.decoded) <- codes
   return(L.decoded)
}

################################################################################

sysanal.package <- function(package)
{
   if ( is.na(match(package,installed.packages()[,1])) )
   {
      print(paste("Do you agree to install the package \"",
                  package,
                  "\"?",
                  sep=""))
      if ( menu(choices=c("yes","no")) == 1 )
      {
         install.packages(package)
      }
      else
      {
         stop(paste("Package \"",package,"\" not found",sep=""))
      }
   }
   library(package,character.only=TRUE)
}

################################################################################

# sysanal.plot.res
# ================

# purpose:
# plot results provided as a numeric vector with result codes

# arguments:
# res:         vector of results named by result codes 
#              (see sysanal.decode for an explanation of result codes)
# xlim:        optional limits of the x axis
# ylim:        optional limits of the y axis
# markers:     if TRUE plot markers instead of lines
# header:      optional header of the plot
# xlab:        optional label of the x axis
# ylab:        optional label of the y axis
# pos:         position of legend (only if more than one variable encoded)

# output:
# plot of all variables as a function of the independent variable
# (note that variable names and values of the independent variable are
# encoded in the names of the components of the result vector)

sysanal.plot.res <- function(res,xlim=NA,ylim=NA,
                             markers=F,header="",
                             xlab="",ylab="",pos="topright")
{
   codes <- sysanal.decode(res)
   varnames <- unique(codes$var)
   if ( is.na(xlim[1]) )      xlim <- range(codes$val)
   if ( is.na(ylim[1]) )      ylim <- range(res)
   if ( nchar(ylab[1]) == 0 ) ylab <- paste(varnames,collapse=", ")
   plot(numeric(0),numeric(0),type="n",xlim=xlim,ylim=ylim,
        xlab=xlab,ylab=ylab,main=header)
   if ( markers )
   {
      for ( i in 1:length(varnames) )
      {
         ind <- codes$var == varnames[i]
         points(codes$val[ind],res[ind],pch=i)
      }
      if ( length(varnames) > 1 )
      {
         legend(x=pos,legend=varnames,pch=1:length(varnames))
      }
   }
   else
   {
      for ( i in 1:length(varnames) )
      {
         ind <- codes$var == varnames[i]
         lines(codes$val[ind],res[ind],lty=i)
      }
      if ( length(varnames) > 1 )
      {
         legend(x=pos,legend=varnames,lty=1:length(varnames))
      }
   }
}

################################################################################

# sysanal.plot.vars
# =================

# purpose:
# plot variables provided as a data frame or matrix with result codes given 
# by the row names

# arguments:
# vars:        matrix or data frame with variables and result codes as row names 
#              (see sysanal.decode for an explanation of result codes)
# ncol:        optional number of columns of sub-panels of the plot
# mar:         optional specification of margins in the form 
#              c(bottom,left,top,right)
# ylim:        optional named (by variable name) list of limits of the y axes
# markers:     if TRUE plot markers instead of lines
# header:      optional named (by variable name) list of headers of the plots
# xlab:        optional label of the x axis
# ylab:        optional label of the y axis
# pos:         position of legend (only if more than one variable encoded)

# output:
# plot of all variables as a function of the independent variable
# (note that variable names and values of the independent variable are
# encoded in the names of the components of the result vector)

sysanal.plot.vars <- function(vars,ncol=NA,mar=NA,
                              ylim=list(),markers=F,
                              headers=list(),xlab="",ylab="",pos="topright")
{
   nvar <- ncol(vars)
   if ( is.na(ncol) ) nc <- floor(sqrt(nvar))
   nr <- ceiling(nvar/nc)
   marg <- mar
   if ( is.na(marg[1]) ) marg <- c(4.5,4.0,2.5,1.0) # c(bottom, left, top, right)
   par.def <- par(no.readonly=T)
   par(mfrow=c(nr,nc),xaxs="i",yaxs="i",mar=marg) 
   for ( i in 1:nvar )
   {
      name <- colnames(vars)[i]
      ylim.i <- c(min(vars[,i]),max(vars[,i]))
      if ( ylim.i[1] == ylim.i[2] ) 
      {
         ylim.i[1] <- 0.5*ylim.i[1]
         ylim.i[2] <- 1.5*ylim.i[2]
      }
      if ( length(ylim) > 0 )
      {
         ind <- match(name,names(ylim))
         if ( !is.na(ind) )
         {
            ylim.i <- ylim[[ind]]
         }
      }
      header.i <- name
      if ( length(headers) > 0 )
      {
         ind <- match(name,names(headers))
         {
            if ( !is.na(ind) )
            {
               header.i <- headers[[ind]]
            }
         }
      }
      res <- as.numeric(vars[,i])
      names(res) <- rownames(vars)
      sysanal.plot.res(res,
                       header=header.i,markers=markers,
                       xlab=xlab,ylab=ylab,pos=pos,
                       xlim=NA,ylim=ylim.i)
   }
   par(par.def)
}

################################################################################

# Contour the probability density function of a bivariate normal or 
# lognormal distribution

sysanal.contourpdf <- function(calcpdf_mv,norm=T,xlim=c(-3,3),ylim=c(-3,3),
                               levels=c(0.05,0.5,0.95),res=20,lty="solid",
                               ...)
{
   # -----------------------------------------------------------------------
   # This function plots contour lines of a normalized or unnormalized
   # bivariate probability density function. If the function is not 
   # normalized, the integral over the domain specified by xlim and ylim is  
   # used for normalization (this leads to incorrect results if this domain 
   # does not contain most of the distribution).
   #
   # Arguments:
   # calcpdf_mv: function to calculate the log pdf values at a set of
   #             locations specified by its first argument. Further 
   #             arguments specified under ... will be passed to this 
   #             function
   # norm:       TRUE if the probability density is normalized,
   #             FALSE if integration over the domain given by xlim and ylim
   #             should be used for normalization
   # xlim:       bounds of the integration range for the first variable
   # ylim:       bounds of the integration range for the second variable 
   # levels:     vector of probabilities to be contained in the contour line
   # res:        resolution of grid used to countour
   #             (number of points in each dimension)
   # lty:        line type of contour lines
   #
   # Return Value:
   # integral of the probability density over the given range at the given
   # resolution
   #
   #                                        Peter Reichert    Feb.  02, 2005
   #                                        last modification March 31, 2008
   # -----------------------------------------------------------------------

   dx.grid <- (xlim[2]-xlim[1])/res
   dy.grid <- (ylim[2]-ylim[1])/res
   x.grid <- seq(xlim[1]+dx.grid/2,xlim[2]-dx.grid/2,len=res)
   y.grid <- seq(ylim[1]+dy.grid/2,ylim[2]-dy.grid/2,len=res)
   xarray <- vector()
   for ( i in 1:res ) xarray = c(xarray,rep(x.grid[i],res))
   yarray <- rep(y.grid,res)
   z.sample <- cbind(xarray,yarray)
   logpdf <- calcpdf_mv(z.sample,...)
   if ( norm == F) logpdf <- logpdf-max(logpdf,na.rm=T)
   pdf <- ifelse(is.na(logpdf),0,exp(logpdf))
   pdf.sort <- sort(pdf,decreasing=T)
   integral.cum <- pdf.sort*dx.grid*dy.grid
   for ( i in 2:(res*res) ) 
   {
      integral.cum[i] <- integral.cum[i-1]+integral.cum[i]
   }
   integral <- integral.cum[res*res]
   if ( norm == F) integral.cum <- integral.cum/integral
   logpdflevels <- vector()
   index <- 1
   levelssort <- sort(levels)
   for ( i in 1:(res*res) )
   {
      if ( integral.cum[i] > levelssort[index] )
      {
         logpdflevels[index] <- log(pdf.sort[i])
         index <- index + 1
         if ( index > length(levelssort) ) break
      }
   }
   logpdf <- matrix(logpdf,nrow=res,ncol=res,byrow=TRUE)
   contour(x=x.grid,y=y.grid,z=logpdf,levels=logpdflevels,
           add=TRUE,drawlabels=FALSE,lty=lty)
   return(integral)
}

################################################################################

# function to plot markov chains:
# -------------------------------

sysanal.plot.chains <- function(postsamp,ncol=NA,mar=NA,
                                ylim=list(),
                                titles=list(),xlab="chain index",ylab=list())
{
   nvar <- ncol(postsamp)
   if ( is.na(ncol) ) nc <- floor(sqrt(nvar))
   nr <- ceiling(nvar/nc)
   marg <- mar
   if ( is.na(marg[1]) ) marg <- c(2.0,2.0,2.5,0.5) # c(bottom, left, top, right)
   par.def <- par(no.readonly=T)
   par(mfrow=c(nr,nc),xaxs="i",yaxs="i",mar=marg) 
   for ( i in 1:nvar )
   {
      name <- colnames(postsamp)[i]
      data <- postsamp[,i]
      ylim.i <- c(min(data),max(data))
      if ( ylim.i[1] == ylim.i[2] ) 
      {
         ylim.i[1] <- 0.5*ylim.i[1]
         ylim.i[2] <- 1.5*ylim.i[2]
      }
      if ( length(ylim) > 0 )
      {
         ind <- match(name,names(ylim))
         if ( !is.na(ind) )
         {
            ylim.i <- ylim[[ind]]
         }
      }
      title.i <- name
      if ( length(titles) > 0 )
      {
         ind <- match(name,names(titles))
         {
            if ( !is.na(ind) )
            {
               title.i <- titles[[ind]]
            }
         }
      }
      ylab.i <- name
      if ( length(ylab) > 0 )
      {
         ind <- match(name,names(ylab))
         if ( !is.na(ind) )
         {
            ylab.i <- ylab[[ind]]
         }
      }
      plot(data,ylim=ylim.i,type="l",main=title.i,xlab=xlab,ylab=ylab.i)
   }
   par(par.def)
}

################################################################################

# function to plot prior and posterior marginals:
# -----------------------------------------------

sysanal.plot.margs <- function(postsamp,pridist=list(),ncol=NA,mar=NA,
                               xlim=list(),ymax=list(),
                               titles=list(),xlab=list(),ylab=list(),
                               lty=NA)
{
   # transform samples to a list of all samples to plot the marginals of:
   
   postsamp.list <- list()
   if ( is.data.frame(postsamp) )
   {
      postsamp.list[[1]] <- postsamp
   }
   else
   {
      if ( is.matrix(postsamp) )
      {
         postsamp.list[[1]] <- postsamp
      }
      else
      {
         if ( is.list(postsamp) )  # be careful: a data frame is a list!
         {
            postsamp.list <- postsamp
         }
         else
         {
            stop("sysanal.plot.margs: postsamp is of illegal type")
         }
      }
   }
   nsamp <- length(postsamp.list)
   
   # transform prior definitions to a list of priors to plot:
   
   pridist.list <- list()
   if ( length(pridist) > 0 )
   {
      if ( length(names(pridist)) != length(pridist) )
      {
         pridist.list <- pridist
      }
      else
      {
         pridist.list[[1]] <- pridist
      }
   }
   npri <- length(pridist.list)
   
   # get all variable names of the sample: 
   
   var <- colnames(postsamp.list[[1]])
   if ( nsamp > 1 )
   {
      for ( j in 2:nsamp ) 
      {
         var <- c(var,colnames(postsamp.list[[j]]))
      }
      var <- unique(var)
   }
   nvar <- length(var)
   
   # define layout of plot panels:

   if ( !is.na(ncol) ) nc <- ncol
   else                nc <- floor(sqrt(nvar))
   nr <- ceiling(nvar/nc)
   marg <- mar
   if ( is.na(marg[1]) ) marg <- c(2.0,2.0,2.5,0.5) # c(bottom, left, top, right)

   # define line types:
   
   lty.loc <- lty
   if ( is.na(lty.loc[1]) ) lty.loc <- 1
   if ( length(lty.loc) < nsamp+npri ) lty.loc <- 1:(nsamp+npri)
   
   # plot marginals of samples and priors:
   
   par.def <- par(no.readonly=T)
   par(mfrow=c(nr,nc),xaxs="i",yaxs="i",mar=marg) # c(bottom, left, top, right)
   for ( i in 1:nvar )
   {
      name <- var[i]
      data.min    <- NA
      data.max    <- NA
      density.max <- NA
      d <- as.list(rep(NA,nsamp))
      for ( j in 1:nsamp )
      {
         ind <- match(var[i],colnames(postsamp.list[[j]]))
         if ( !is.na(ind) )
         {
            data.min    <- min(c(data.min,postsamp.list[[j]][,ind]),na.rm=TRUE)
            data.max    <- max(c(data.max,postsamp.list[[j]][,ind]),na.rm=TRUE)
            bw <- .5*sd(postsamp.list[[j]][,ind])
            d[[j]]      <- density(postsamp.list[[j]][,ind],bw=bw) ### changed!###
#             d[[j]]      <- density(postsamp.list[[j]][,ind],adjust=1.0)
            density.max <- max(c(density.max,d[[j]]$y),na.rm=TRUE)
         }
      }
      
      xlim.i <- c(data.min,data.max)
      if ( xlim.i[1] == xlim.i[2] ) 
      {
         xlim.i[1] <- 0.5*xlim.i[1]
         xlim.i[2] <- 1.5*xlim.i[2]
      }
      if ( length(xlim) > 0 )
      {
         ind <- match(name,names(xlim))
         if ( !is.na(ind) )
         {
            xlim.i <- xlim[[ind]]
         }
      }

      if ( is.na(density.max) ) ylim.i <- c(0,1)
      else                      ylim.i <- c(0,1.1*density.max)
      if ( length(ymax) > 0 )
      {
         ind <- match(name,names(ymax))
         if ( !is.na(ind) )
         {
            ylim.i <- c(0,ymax[[ind]])
         }
      }

      title.i <- name
      if ( length(titles) > 0 )
      {
         ind <- match(name,names(titles))
         {
            if ( !is.na(ind) )
            {
               title.i <- titles[[ind]]
            }
         }
      }
      
      xlab.i <- name
      if ( length(xlab) > 0 )
      {
         ind <- match(name,names(xlab))
         if ( !is.na(ind) )
         {
            xlab.i <- xlab[[ind]]
         }
      }

      ylab.i <- "f"
      if ( length(ylab) > 0 )
      {
         ind <- match(name,names(ylab))
         if ( !is.na(ind) )
         {
            ylab.i <- ylab[[ind]]
         }
      }

      plot(numeric(0),numeric(0),main=title.i,xlab=xlab.i,ylab=ylab.i,
           xlim=xlim.i,ylim=ylim.i,type="n")

      for ( j in 1:nsamp )
      {
         if ( !is.na(d[[j]])[[1]] )
         {
#             lines(d[[j]],lty=lty.loc[j])
            polygon(x=c(d[[j]]$x, rev(d[[j]]$x)),
                    y=c(d[[j]]$y, rev(rep(0, length(d[[j]]$y)))),### changed!###
                    col = "lightgrey", border=NA)
         }
      }

      if ( npri > 0 )
      {
         x <- seq(xlim.i[1],xlim.i[2],length=101)
         y <- numeric(0)
         for ( j in 1:npri )
         {
            if ( length(pridist.list[[j]]) > 0 )
            {
               ind <- sysanal.multmatch(name,names(pridist.list[[j]]))
               if ( !is.na(ind[1]) )
               {
                  for ( k in 1:length(ind) )
                  {
                     y <- sysanal.calcpdf(x,pridist.list[[j]][[ind[k]]])
                     lines(x,y,lty=lty.loc[nsamp+j])
                  }
               }
            }
         }
      }
   }
   par(par.def)
}


sysanal.multmatch <- function(x,table,incomparables=NULL)
{
   table.loc <- table
   inds <- numeric(0)
   while (TRUE)
   {
      ind <- match(x,table.loc,incomparables)
      if ( is.na(ind) )
      {
         if ( length(inds) == 0 ) return(ind)
         else                     return(inds)
      }
      inds <- c(inds,ind)
      table.loc[ind] <- paste(x,"_",sep="")
   }  
}

################################################################################

# generate a random sample from a multivariate normal or lognormal 
# distribution

sysanal.randsamp <- function(sampsize=1,dist="normal",mean=0,sd=1,cor=NA,
                             file=NA)
{
   # -----------------------------------------------------------------------
   # This function generates a random sample from a multivariate
   # normal or lognormal distribution.
   # This function is a simplified version of the program "randsamp"
   # available as part of the package UNCSIM at http://www.uncsim.eawag.ch
   #
   # Arguments:
   # sampsize:   sample size
   # mean:       vector of means
   # sd:         vector of standard deviations
   # cor:        correlation matrix of the distribution
   # dist:       distribution type: "normal" or lognormal".
   #
   # Return Value:
   # List of:
   # mean:       vector of means
   # sd:         vector of standard deviations
   # corr:       correlation matrix of the distribution
   # sampsize:   sample size
   # sample:     matrix of parameter samples (each row corresponds 
   #             to a draw)
   # pdf:        vector of values of probability density function at the
   #             sample points
   #
   #                                        Peter Reichert    Dec.  29, 2003
   #                                        last modification March 30, 2008
   # -----------------------------------------------------------------------

   # consistency checks and initializations:
   numpar <- length(mean)
   if ( length(sd) != numpar )
   {
      stop("sysanal.randsamp: mean and sd do not have the same length")
   }
   R <- diag(rep(1,numpar))
   if ( is.matrix(cor) ) R <- cor
   if ( nrow(R) != numpar || ncol(R) != numpar )
   {
      stop("sysanal.randsamp: illegal dimension of correlation matrix")
   }

   # calculate sample from multivariate uniform distribution:
   samp <- runif(sampsize*numpar)
   dim(samp) <- c(sampsize,numpar)

   # transform sample to multivariate normal or lognormal and calculate
   # logarithm of probability density function:
   logpdf <- numeric(sampsize)
   if ( dist == "normal" )
   {
      # calculate transformation matrix and transform the sample:
      sigma <- diag(sd) %*% R %*% diag(sd)
      sigma.inv = solve(sigma)
      det.sigma = det(sigma)
      A <- t(chol(sigma))
      for ( i in 1:sampsize )
      {
         samp[i,]  <- A %*% qnorm(samp[i,]) + mean
         v <- samp[i,]-mean
         v <- matrix(v,nrow=numpar,ncol=1) # ensure that v is a column vector
         logpdf[i] <- -numpar/2*log(2*pi) - 0.5*log(det.sigma) -
                      0.5 * t(v) %*% sigma.inv %*% (v)
      }
   }
   else
   {
      if ( dist == "lognormal" | dist == "Lognormal" )
      {
         # parameters of the log of the variable, calculate transformation 
         # matrix and transform the sample:
         sdlog    <- sqrt(log(1+sd*sd/(mean*mean)))
         meanlog  <- log(mean) - sdlog*sdlog/2
         if ( numpar > 1 )
         {
            ln.sigma <- log( 1 + diag(sqrt(exp(sdlog*sdlog)-1)) %*%
                                 R %*% diag(sqrt(exp(sdlog*sdlog)-1)) )
         }
         else
         {
            ln.sigma <- as.matrix( log( 1 + sqrt(exp(sdlog*sdlog)-1)^2 ) )
         }
         ln.sigma.inv = solve(ln.sigma)
         det.ln.sigma = det(ln.sigma)
         ln.A <- t(chol(ln.sigma))
         for ( i in 1:sampsize )
         {
            log.samp.i <- ln.A %*% qnorm(samp[i,]) + meanlog
            samp[i,]  <- exp(log.samp.i)
            v <- log.samp.i-meanlog
            v <- matrix(v,nrow=numpar,ncol=1) # ensure that v is a column vector
            logpdf[i] <- -numpar/2*log(2*pi) - 0.5*log(det.ln.sigma) - 
                         log(prod(samp[i,])) - 0.5 * t(v) %*% ln.sigma.inv %*% v
         }
      }
      else
      {
         stop("sysanal.randsamp: unknown distribution type")
      }
   }

   # collect results:
   colnames(samp) <- names(mean)
   res <- list(
               mean       = mean,
               sd         = sd,
               cor        = R,
               sampsize   = sampsize,
               sample     = samp,
               logsamppdf = logpdf
              )
              
   # write results:
   if ( !is.na(file) )
   {
      write.table(data.frame(samp,logsamppdf=logpdf),file=file,
                  col.names=TRUE,row.names=FALSE,sep="\t")
   }

   # return results:
   return(res)
}

################################################################################

# generate a random sample from a uniform distribution in a ball of given radius

# see introduction of Harman, R. and Lacko, V., On decompositional algorithms
# for uniform sampling from n-spheres and n-balls, Journal of Multivariate 
# Analysis 101, 2297-2304, 2010 for details and references 

sysanal.randsamp.ball <- function(sampsize=1,dim=1,radius=1)
{
   res <- matrix(data=rnorm(sampsize*dim),nrow=sampsize)
   norm <- function(x) { return(sqrt(sum(x^2))) }
   rnorms <- apply(res,1,norm)
   res <- diag(radius/rnorms) %*% res
   res <- diag(runif(sampsize)^(1/dim)) %*% res
   return(res) 
}

################################################################################

# Calculate a Metropolis Markov Chain sample of a distribution

sysanal.markovchain.metropolis <- function(log.pdf,z.ini,
                                           jump.sd,jump.cor=0,
                                           sampsize,thin=1,...)
{
   # -----------------------------------------------------------------------
   # This function calculates a Markov Chain of a probability distribution
   # of a vector of continuous random variables using the Metropolis
   # algorithm with a normal jump distribution with given standard
   # deviations and correlation structure.
   # The log of the probability density of the distribution must be 
   # specified as a function log.pdf(z,...) where z is the vector of values
   # for which the density has to be evaluated.
   #
   # Arguments:
   # log.pdf      function "log.pdf(z,...)" that calculates the log of the
   #              probability density function of the vector of random 
   #              variables that are to be sampled.
   #              log.pdf must be given a single vector z at which the
   #              density has to be evaluated. Additional arguments from
   #              the call to calc.markovchain.metropolis will be passed on.
   #              If the probability density is zero, NA must be returned.
   # z.ini        vector of values at which the chain is to be started.
   # jump.sd      vector of standard deviations of the jump distribution.
   # jump.cor     correlation matrix of jump distribution or NA if all
   #              correlations are zero.
   # sampsize     sample size (length of the chain)
   # thin         factor with which to thin storage of results (thin=n: 
   #              only each nth result is returned; this saves memory)
   # ...          additional parameters passed to "log.pdf"
   #
   # Return Value:
   # List with the following elements:
   # z            sample as a matrix with sample points in its rows.
   # log.pdf      vector with log pdf values of the sample.
   # reject.freq  rejection frequency of the jumps.
   # error        error message (empty string if no error occurred).
   #
   #                                first version:         Dec.  08, 2007 PR
   #                                add parameter "thin":  March 26, 2008 PR
   #                                minor modification:    March 28, 2009 PR
   # -----------------------------------------------------------------------

   # set up and initialize arrays:

   returnsize         <- floor(sampsize/thin)+1
   z.sample           <- matrix(data=NA,nrow=returnsize,ncol=length(z.ini))
   colnames(z.sample) <- names(z.ini)
   log.pdf.sample     <- rep(NA,returnsize)
   reject.freq        <- 1
   error              <- ""

   # calculate Cholesky decomposition of variance-covariance matrix:

   R <- diag(rep(1,length(jump.sd)))
   if ( is.matrix(jump.cor) ) R <- jump.cor
   if ( (nrow(R)!=length(jump.sd)) | (ncol(R)!=length(jump.sd)) )
   {
      error <- paste("sysanal.markovchain.metropolis:",
                     "illegal dimension of correlation matrix")
   }
   if ( nchar(error) == 0 )
   {
      sigma <- diag(jump.sd) %*% R %*% diag(jump.sd)
      A.chol <- try(t(chol(sigma)),silent=FALSE)
      if( inherits(A.chol,"try-error") )
      {
         error <- paste("sysanal.markovchain.metropolis:",
                        "unable to calculate Cholesky decomposition of variance-covariance matrix")
      }
   }
   
   # initialize Markov chain:

   if ( nchar(error) == 0 )
   {
      z.current         <- z.ini
      log.pdf.current   <- log.pdf(z.current,...)
      z.sample[1,]      <- z.current
      log.pdf.sample[1] <- log.pdf.current
      if ( is.na(log.pdf.sample[1]) )
      {
         error <- paste("sysanal.markovchain.metropolis:",
                        "probability density is zero at initial point of chain")
      }
   }
   if ( nchar(error) == 0 )
   {
      num.accept <- 0
      for ( i in 2:returnsize )
      {
         for ( j in 1:thin )
         {
            # calculate suggested new sample point:

            jump.unif <- runif(length(z.ini),min=0,max=1)
            jump <- A.chol %*% qnorm(jump.unif)
            z.suggested <- z.current + as.vector(jump)

            # calculate log pdf at suggested new sample point:

            log.pdf.suggested <- log.pdf(z.suggested,...)

            # accept new point with probability r=pdf.suggested/pdf.prev

            accept <- FALSE
            if ( is.finite(log.pdf.suggested) )
            {
               if ( log.pdf.suggested > log.pdf.current )
               {
                  accept <- TRUE
               }
               else
               {
                  r <- exp(log.pdf.suggested-log.pdf.sample[i-1])
                  if ( runif(n=1,min=0,max=1) <= r )
                  {
                     accept <- TRUE
                  }
               }
            }
            if ( accept == TRUE )
            {
               z.current       <- z.suggested
               log.pdf.current <- log.pdf.suggested
               num.accept      <- num.accept+1
            }
            reject.freq <- ((i-2)*thin+j-1-num.accept)/((i-2)*thin+j-1)
         }
         z.sample[i,]      <- z.current
         log.pdf.sample[i] <- log.pdf.current
      }
   }
   
   # collect and return results:

   res <- list(z           = z.sample,
               log.pdf     = log.pdf.sample,
               reject.freq = reject.freq,
               error       = error)
   return(res)
}

################################################################################

# Calculate an importance sample of a distribution

sysanal.impsamp <- function(log.pdf,z,z.log.pdf,...)
{
   w <- rep(NA,nrow(z))
   log.pdf.values <- w
   for ( i in 1:nrow(z) )
   {
      par <- z[i,]
      names(par) <- colnames(z)
      log.pdf.values[i] <- log.pdf(par,...)
   }
   log.pdf.max <- max(log.pdf.values,na.rm=TRUE)
   w <- exp(log.pdf.values-log.pdf.max-z.log.pdf)
   w <- ifelse ( is.na(w),0,w )
   w <- w/sum(w)
   ess <- sum(w)^2/sum(w^2)
   return(list(z=z,w=w,ess=ess,
               log.pdf.dist=log.pdf.values,log.pdf.samp=z.log.pdf))
}


################################################################################

# function to improve jump distribution;

sysanal.jump.dist <- function(postsamp,fact.sd,fract.burnin=0,fact.cor=1,
                              plot=F)
{
   ind.end   <- nrow(postsamp)
   ind.start <- as.integer(fract.burnin*(ind.end-1))+1
   postsamp.local <- postsamp[ind.start:ind.end,sd(postsamp)!=0]
   if ( is.vector(postsamp.local) )
   {
      postsamp.local <- as.matrix(postsamp.local,nrow=length(postsamp.local))
      colnames(postsamp.local) <- colnames(postsamp)[sd(postsamp)!=0]
      postsamp.local <- data.frame(postsamp.local)
   }
   sd  <- fact.sd*sd(postsamp.local)
   cor <- NA
   if ( ncol(postsamp.local) > 1 )
   {
      corr <- cor(postsamp.local)
      if ( plot )
      {
         image(x=1:ncol(corr),y=1:nrow(corr),z=abs(corr),zlim=c(0,1),
               col=grey((100:0)/100),xlab="variable index",ylab="variable index",
               main="structure of correlation matrix")
         abline(v=0.5)
         abline(h=0.5)
         abline(v=ncol(corr)+0.5)
         abline(h=nrow(corr)+0.5)
      }
      corr <- fact.cor*corr
      diag(corr) <- rep(1,nrow(corr))
      try(chol(corr)) # test if Cholesky factorization works 
                      # (important for subsequent sampling)
   }
   return(list(sd=sd,cor=corr))
}

################################################################################

# function to calculate probability densities of univariate distributions:
# ------------------------------------------------------------------------

sysanal.calcpdf <- function(x,distpar,log=FALSE)
{
   if ( distpar[1] == "Uniform" )
   {
      # uniform distribution; parameters are min and max
      min <- as.numeric(distpar[2])
      max <- as.numeric(distpar[3])
      return(dunif(x,min=min,max=max,log=log))
   }
   if ( distpar[1] == "Normal" )
   {
      # normal distribution; parameters are mean and sd:
      mean <- as.numeric(distpar[2])
      sd   <- as.numeric(distpar[3])
      return(dnorm(x,mean=mean,sd=sd,log=log))
   }
   if ( distpar[1] == "NormalTrunc" )
   {
      # truncated normal distribution; parameters are mean, sd, min and max
      mean <- as.numeric(distpar[2])
      sd   <- as.numeric(distpar[3])
      min  <- as.numeric(distpar[4])
      max  <- as.numeric(distpar[5])
      fact <- 1/(pnorm(q=max,mean=mean,sd=sd)-pnorm(q=min,mean=mean,sd=sd))
      if ( !log )
      {
         return(ifelse(x<min|x>max,0,fact*dnorm(x,mean=mean,sd=sd)))
      }
      else
      {
         return(ifelse(x<min|x>max,-Inf,
                       log(fact)+dnorm(x,mean=mean,sd=sd,log=TRUE)))
      }
   }
   if ( distpar[1] == "Lognormal" )
   {
      # lognormal distribution; parameters are mean and sd;
      # R parameters are mean and sd of the log of the random variable
      mean    <- as.numeric(distpar[2])
      sd      <- as.numeric(distpar[3])
      sdlog   <- sqrt(log(1+sd^2/mean^2))
      meanlog <- log(mean) - 0.5*sdlog^2
      return(dlnorm(x,meanlog=meanlog,sdlog=sdlog,log=log))
   }
   if ( distpar[1] == "LognormalTrunc" )
   {
      # truncated lognormal distribution; parameters are mean, sd, min and max;
      # R parameters are mean and sd of the log of the random variable
      mean    <- as.numeric(distpar[2])
      sd      <- as.numeric(distpar[3])
      sdlog   <- sqrt(log(1+sd^2/mean^2))
      meanlog <- log(mean) - 0.5*sdlog^2
      min     <- as.numeric(distpar[4])
      max     <- as.numeric(distpar[5])
      fact <- 1/(plnorm(q=max,meanlog=meanlog,sdlog=sdlog)-plnorm(q=min,meanlog=meanlog,sdlog=sdlog))
      if ( !log )
      {
         return(ifelse(x<min|x>max,0,fact*dlnorm(x,meanlog=meanlog,sdlog=sdlog)))
      }
      else
      {
         return(ifelse(x<min|x>max,NA,
                       log(fact)+dlnorm(x,meanlog=meanlog,sdlog=sdlog,log=TRUE)))
      }
   }
   if ( distpar[1] == "Inv" )
   {
      # inverse distribution; parameters are min and max:
      min <- as.numeric(distpar[2])
      max <- as.numeric(distpar[3])
      if ( !log )
      {
         return(ifelse(x<min|x>max,0,1/(log(max/min)*x)))   
      }
      else
      {
         return(ifelse(x<min|x>max,NA,-log(log(max/min)) - log(x)))   
      }
   }
   if ( distpar[1] == "Exponential" )
   {
      # exponential distribution; parameter is mean:
      mean <- as.numeric(distpar[2])
      if ( !log )
      {
         return(ifelse(x<0,0,1/mean*exp(-x/mean)))   
      }
      else
      {
         return(ifelse(x<0,-Inf,-log(mean)-x/mean)) #!!!
#          return(ifelse(x<0,NA,-log(mean)-x/mean))
      }
      
#     if ( distpar[1] == "Delta" )
#       {
#         # exponential distribution; parameter is mean:
#         mean <- as.numeric(distpar[2])
# 
#           return(ifelse(x!= 0,-Inf,1) #!!!
#           #          return(ifelse(x<0,NA,-log(mean)-x/mean))
#         }   
   }
   stop(paste("Distribution",dist,"not yet implemented"))
}

################################################################################

# Calculate the logarithm of the probability density function of a 
# multivariate normal or lognormal distribution or of a product of independent
# marginals

sysanal.calcpdf_mv <- function(z,dist="normal",mean=0,sd=1,cor=0,
                               cor.inv=NA,log=TRUE,distdef=NA,file=NA)
{
   # -----------------------------------------------------------------------
   # This function calculates the logarithm of the probability density 
   # function of a multivariate normal or lognormal distribution.
   #
   # Arguments:
   # z:          vector, matrix or data frame at which the logarithm of the 
   #             probability density function has to be evaluated
   # dist:       distribution type: "normal" or lognormal".
   # mean:       vector of means
   # sd:         vector of standard deviations
   # cor:        correlation matrix of the distribution
   #
   # Return Value:
   # logarithm of the probability density function for all samples in x
   #
   #                                        Peter Reichert    Jan.  01, 2004
   #                                        last modification June  18, 2010
   # -----------------------------------------------------------------------

   # consistency checks and initializations:
   mean <- as.vector(mean)
   sd <- as.vector(sd)
   if ( length(sd) != length(mean) )
   {
      stop("sysanal.calcpdf_mv: illegal dimension of standard deviations")
   }
   if ( is.vector(z) )
   {
      len <- length(z)
      names <- names(z)
      z <- as.matrix(z)
      dim(z) <- c(1,len)
      if ( length(names) == len ) colnames(z) <- names
   }
   numpar <- ncol(z)
   R <- diag(rep(1,numpar))
   if ( is.matrix(cor) ) R <- cor
   if ( nrow(R) != numpar || ncol(R) != numpar )
   {
      stop("sysanal.calcpdf_mv: illegal dimension of correlation matrix")
   }

   # calculate logarithm of probability density function:
   sampsize <- nrow(z)
   logpdf <- numeric(sampsize)
   if ( dist == "normal" )
   {
      # multivariate normal distribution:
      n <- length(sd)
      sigma <- diag(sd,nrow=n,ncol=n) %*% R %*% diag(sd,nrow=n,ncol=n)
      if ( is.matrix(cor.inv) )
      {
         R.inv <- cor.inv
      }
      else
      {
         R.inv <- solve(R)
      }
      det.R.inv <- det(R.inv)
      if ( det.R.inv > 0 )
      {
         for ( i in 1:sampsize )
         {
            v <- as.matrix(z[i,]-mean,nrow=numpar,ncol=1)/
                 as.matrix(sd,nrow=numpar,ncol=1)
            logpdf[i] <- -numpar/2*log(2*pi) + 0.5*log(det.R.inv) - log(prod(sd)) -
                         0.5 * t(v) %*% R.inv %*% v
         }
      }
      else
      {
         logpdf <- rep(NA,sampsize)
      }
   }
   else
   {
      if ( dist == "lognormal" )
      {
         # multivariate lognormal distribution:
         sdlog    <- sqrt(log(1+sd*sd/(mean*mean)))
         meanlog  <- log(mean) - sdlog*sdlog/2
         if ( numpar > 1 )
         {
            n <- length(sdlog)
            ln.sigma <- log( 1 + diag(sqrt(exp(sdlog*sdlog)-1),nrow=n,ncol=n) %*% 
                                 R %*% diag(sqrt(exp(sdlog*sdlog)-1),nrow=n,ncol=n) )
         }
         else
         {
            ln.sigma <- as.matrix( log( 1 + sqrt(exp(sdlog*sdlog)-1)^2 ) )
         }
         ln.sigma.inv = solve(ln.sigma)
         det.ln.sigma = det(ln.sigma)
         if ( det.ln.sigma > 0 )
         {
            for ( i in 1:sampsize )
            {
               if ( min(z[i,]) <= 0 )
               {
                  logpdf[i] <- NA 
               }
               else
               {
                  v <- log(z[i,])-meanlog
                  v <- matrix(v,nrow=numpar,ncol=1) # ensure that v is a column vector
                  logpdf[i] <- -numpar/2*log(2*pi) - 0.5*log(det.ln.sigma) - 
                               log(prod(z[i,])) - 0.5 * t(v) %*% ln.sigma.inv %*% v
               }
            }
         }
         else
         {
            logpdf <- rep(NA,sampsize)
         }
      }
      else
      {
         if ( dist == "indep" ) 
         { 
            for ( i in 1:sampsize )
            {
               logpdf[i] <- 0
               for ( j in 1:ncol(z) )
               {
                  ind <- j
                  if ( length(colnames(z)) == ncol(z) )
                  {
                     ind <- match(colnames(z)[j],names(distdef))
                     if ( is.na(ind) )
                     {
                        stop(paste("error in calcpdf_mv:",
                                   "variable",colnames(z)[j],"not found"))
                     }
                  }
                  # ****** new for hierarchical models ********** 
                  # when j is found so that th1 is a latent variable...
                  if (  !is.na(any(match(distdef[[ind]],names(distdef)))) ) 
                 
                  # compute p(th1|th2)
                  {
                                
                   pos.th2  <-  which(!is.na(match(distdef[[ind]],names(distdef))))    #position of the hierarch param
                   
                   distdef[[ind]][pos.th2] = z[i,distdef[[ind]][pos.th2]] # instead of the name put the true sampled hyperparameter
                    
                  }
                  
                  # *********************************************
                  
                  logpdf[i] <- logpdf[i] + sysanal.calcpdf(z[i,j],
                                                           distdef[[ind]],
                                                           log=TRUE)
#           
               }
            }
         }
         else
         {
            stop(paste("sysanal.calcpdf_mv: unknown distribution type:",
                       dist))
         }
      }
   }
   
   # write results:
   if ( !is.na(file) )
   {
      write.table(data.frame(z,logpdf=logpdf),file=file,
                  col.names=TRUE,row.names=FALSE,sep="\t")
   }

   # return result:
   if ( log )
   {
      return(logpdf)
   }
   else
   {
      return(exp(logpdf))
   }
}

################################################################################

# Function implementing a simple standard likelihood function for a model
# provided by a deterministic function. This function provides the density of
# a multivariate normal distribution centered at the results of the 
# deterministic model.
# If there exists a component of the parameter vector labelled "sd_rel" 
# instead of the provided standard deviations, "error.sd" the standard 
# deviations "error.sd * par["sd_rel"]" are used. This makes it possible to
# estimate a common factor of a given variance-covariance structure.
# This likelihood implementation serves as a template for implementing more
# specific likelihood functions. It is called by "sysanal.logposterior".

sysanal.loglikeli <- function(par,model,L,y,
                              error.sd=NA,error.cor.inv=NA,...)
{
   sd.rel <- 1; if ( !is.na(par["sd_rel"]) ) sd.rel <- par["sd_rel"]
   res <- numeric(0)
   if ( !is.matrix(par) )  # single model evaluation at parameter vector par
   {
      mean <- model(par,L=L,...)
      sd <- error.sd; if ( is.na(sd[1]) ) sd <- rep(1,length(mean))
      if ( is.na(error.cor.inv[1]) )
      {
         res <- 0;
         for ( i in 1:length(mean) )
         {
            res <- res + dnorm(y[i],mean[i],sd[i]*sd.rel,log=TRUE)
         }
      }
      else 
      {
         res <- sysanal.calcpdf_mv(z=y,dist="normal",
                                   mean=mean,sd=sd*sd.rel,
                                   cor.inv=error.cor.inv)
      }
   }
   else                    # multiple model evaluations for parameter sample
   {
      for ( i in 1:nrow(par) )
      {
         mean <- model(par[i,],L=L,...)
         sd <- error.sd; if ( is.na(sd[1]) ) sd <- rep(1,length(mean))
         if ( is.na(error.cor.inv[1]) )
         {
            res[i] <- 0;
            for ( i in 1:length(mean) )
            {
               res[i] <- res[i] + dnorm(y[i],mean[i],sd[i]*sd.rel,log=TRUE)
            }
         }
         else 
         {
            res[i] <- sysanal.calcpdf_mv(z=y,dist="normal",
                                         mean=mean,sd=sd*sd.rel,
                                         cor.inv=error.cor.inv)
         }
      }
   }
   return(res)
}

################################################################################

# Function implementing a simple standard posterior probability density (up to
# a multiplicative normalizing constant). It is based on a likelihood 
# definition provided by the user and a multivariate normal or
# lognormal distribution at the prior.
# This log posterior implementation serves as a template for implementing more
# specific posteriors.

sysanal.logposterior <- function(par,model,L,y,
                                 prior.dist="lognormal",prior.mean=1,prior.sd=1,
                                 prior.cor=NA,prior.def=NA,
                                 loglikeli=sysanal.loglikeli,
                                 ...)
{
#    print(par)##
   logprior  <- sysanal.calcpdf_mv(z=par,dist=prior.dist,mean=prior.mean,
                                   sd=prior.sd,cor=prior.cor,distdef=prior.def)
#    print(paste("log prior: ", format(logprior, digits=2,scientific=T)))

   if ( is.na(logprior) ) return(NA)
    loglikeli <- loglikeli(par=par,model=model,L=L,y=y,...)


#   print(paste("log likeli: ", format(loglikeli, digits=2,scientific=T)))
   return(logprior+loglikeli)
}


################################################################################

# sysanal.sens.loc
# ================

# purpose:
# calculate the local sensitivity matrix of a model
# (matrix of partial derivatives of model results with respect to parameters)

# arguments:
# par:         named parameter vector; passed to the model as separate arguments
# model        function representing the model; typically returns a vector,
#              but could also return a matrix
#              (cf. sysanal.model)
# par.inc:     increments of parameters used to approximate the derivatives
# ...          further arguments are passed to model

# output:
# matrix of partial derivatives of model results (rows) 
# with respect to parameters (columns)
# or 3 dimensional array of partial derivatives if model output is a matrix

sysanal.sens.loc <- function(par,model,par.inc=0.01*par,...)
{
   if ( length(par) != length(par.inc) ) 
      stop("*** error in sysanal.sens.loc: unequal length of par and par.inc")
   if ( min(abs(par.inc)) == 0 )
      stop("*** error in sysanal.sens.loc: elements of par.inc must not be zero")
   res.par <- model(par,...)
   V = NA
   if ( is.vector(res.par) )
   {
      V <- matrix(NA,nrow=length(res.par),ncol=length(par))
      colnames(V) <- names(par)
      rownames(V) <- names(res.par)
      for( j in 1:length(par) )
      {
         par.j <- par
         par.j[j] <- par.j[j] + par.inc[j]
         V[,j] <- ( model(par.j,...) - res.par ) / par.inc[j]
      }
   }
   else
   {
      if ( is.matrix(res.par) )
      {
         V <- array(NA,dim=c(dim(res.par),length(par)),
                    dimnames=list(rownames(res.par),
                                  colnames(res.par),
                                  names(res.par)))
         for ( j in 1:length(par) )
         {
            par.j <- par
            par.j[j] <- par.j[j] + par.inc[j]
            V[,,j] <- ( model(par.j,...) - res.par ) / par.inc[j]
         }
      }
   }
   return(V)
}

# version of the same function that calls a model with explicitly listed 
# arguments (this is primarily for internal use when internally calling nls)
# model.name is the the name of the function with explicit parameter arguments

sysanal.sens.loc.explicitpars <- function(par,model.name,par.inc=0.01*par,x)
{
   if ( length(par) != length(par.inc) ) 
      stop("*** error in sysanal.sensfun: unequal length of par and par.inc")
   if ( min(abs(par.inc)) == 0 )
      stop("*** error in sysanal.sensfun: elements of par.inc must not be zero")
   args <- as.list(par)
   if ( length(x) > 0 ) args <- c(args,x)
   res.par <- do.call(model.name,args=args)
   V <- matrix(NA,nrow=length(res.par),ncol=length(par))
   colnames(V) <- names(par)
   rownames(V) <- names(res.par)
   for( j in 1:length(par) )
   {
      par.j <- par
      par.j[j] <- par.j[j] + par.inc[j]
      args <- as.list(par.j)
      if ( length(x) > 0 ) args <- c(args,x)
      res.j <- do.call(model.name,args=args)
      V[,j] <- ( res.j - res.par ) / par.inc[j]
   }
   return(V)
}

################################################################################

# sysanal.sens.var.samp
# ---------------------

# function to calculate the first order variance based sensitivity
# coefficients from a parameter sample and the corresponding result sample
													
# Call:      sysanal.sens.var.samp(parsamp,ressamp,nbin=NA,
#                                  method="smooth",order=2,
#                                  bandwidth=1,span=0.75,sd.rel=0.05,
#                                  plot=F)					
# ------------------------------------------------------------------										
#													
# parsamp	   matrix containing the parameter sample (each row corresponds to
#            a sample point)						
#													
# ressamp	   vector (of only one result per parameter sample point) or matrix 
#            of results corresponding to the parameter sample (each row 
#            provides the results corresponding to the parameter values in 
#            the same row of parsamp)									
#													
# nbin 	     number of quantile intervals for which the conditional means	
#		         are calculated, default is the square root of the sample size								
#
# method     "smooth", "loess", "glkerns", "lokerns", "lpepa", "lpridge"	
#            routine to be used for smoothing
#
# order		   order of local regression polynomial or of kernel
#
# bandwidth  method="lpepa" or method="lpdidge" only: bandwidth of the
#            smoothing algorithm
#
# span       method="loess" only: fration of points used for local regression 
#
# sd.rel     method="smooth" only: standard deviation of normal distribution of 
#            smoothing algorithm relative to the 99% quantile interval
#
#	plot       logical variable indicating if a scatter plot of the relationship
#            between parameter and model output should be plotted								
#													
# Output:												
# -------												
#													
# List with two elements:	
#								
# var	       vector with the total variance of each column (each model	 
#            output) of the ressamp matrix						
#													
# var.cond.E matrix with the variance of the conditional expected value	
#			       of each parameter (columns) and each model output (rows)	

sysanal.sens.var.samp <- function(parsamp,ressamp,nbin=NA,
                                  method="smooth",order=2,
                                  bandwidth=1,span=0.75,sd.rel=0.1,
                                  plot=F)
{
   sysanal.package("lpridge")
   sysanal.package("lokern")

   if ( is.vector(ressamp) ) ressamp <- as.matrix(ressamp,ncol=1)
   npar  <- ncol(parsamp)
   nsamp <- nrow(parsamp)
   nres  <- ncol(ressamp)
   if ( is.na(nbin) ) nbin <- ceiling(sqrt(nsamp))
   if ( nrow(parsamp ) != nrow(ressamp)) 
   {
      stop ("ressamp and parsamp do not have the same row length")
   }

   var_k <- rep(NA,nres)
   names(var_k) <- colnames(ressamp)
   for ( k in 1:nres ) var_k[k] <- var(ressamp[,k])

   var_q <- data.frame(matrix(NA,nrow=nres,ncol=npar))
   colnames(var_q) <- colnames(parsamp)
   rownames(var_q) <- colnames(ressamp)

   for ( i in 1:npar )
   {
      q  <- quantile(parsamp[,i],probs=(((1:nbin)-0.5)/nbin),
                     na.rm=FALSE,names=TRUE,type=7)
      for ( k in 1:nres )
      {
         if ( method == "lpepa" )
         {
            mean_res <- lpepa(parsamp[,i],ressamp[,k],bandwidth=bandwidth,
                              x.out=q,order=order)$est
         }
         else
         {
            if ( method == "lpridge" )
            {
               mean_res <- lpridge(parsamp[,i],ressamp[,k],bandwidth=bandwidth,
                                   x.out=q,order=order)$est
            }
            else
            {
               if ( method == "glkerns" )
               {
                  mean_res <- glkerns(parsamp[,i],ressamp[,k],
                                      x.out=q,korder=order,hetero=T)$est
               }
               else
               {
                  if ( method == "lokerns" )
                  {
                     mean_res <- lokerns(parsamp[,i],ressamp[,k],
                                         x.out=q,korder=order,hetero=T)$est
                  }
                  else
                  {
                     if ( method == "loess" )
                     {
                        #mean_res <- lowess(parsamp[,i],ressamp[,k],f,
                        #                   x.out=q,korder=order,hetero=T)$est
                        res <- 
                           loess(y~x,
                                 data=data.frame(x=parsamp[,i],y=ressamp[,k]),
                                 span=span,degree=order)
                        mean_res <- predict(res,newdata=data.frame(x=q))
                     }
                     else
                     {
                        if ( method == "smooth" )
                        {
                           if ( order == 2 ) m <- "quadratic"
                           else              m <- "linear"
                     	     mean_res <- sysanal.smooth(parsamp[,i],ressamp[,k],
                                                      sigma=((q[nbin]-q[1])*sd.rel),
                                                      newx=q,
                                                      method=m)$y
                        }
                        else
                        {
                           stop(paste("calc.var.sens: method",
                                      method,
                                      "not implemented"))
                        }
                     }
                  }
               }
            }
         }
                                   
  	     var_q[k,i] <- var(mean_res)
         if ( plot )
         {
            plot(parsamp[,i],ressamp[,k],pch=19,cex=0.2,
                 xlab=colnames(parsamp)[i],ylab="Y")
            lines(q,mean_res,lwd=3,col="red")
         }
      }
   }

   return (list(var=var_k,var.cond.E=var_q))
}

################################################################################

# calculation of index combinations (auxiliary function used in ident)
# ====================================================================

sysanal.comb <- function(n,p)
{
   # -----------------------------------------------------------------------
   # This function calculates all combination of subsets of length p out
   # of n indices.
   #
   # Arguments:
   # n:   number of indices.
   # p:   length of subset of indices.
   #
   # Return Value:
   # matrix with subsets of length p as rows.
   #
   #                                         Peter Reichert    Dec. 27, 2002
   # -----------------------------------------------------------------------

   # check input:
   if ( p > n ) stop("comb: illeal arguments (p>n)")

   # initialize array and auxiliary variables:
   num.comb <- choose(n,p)
   comb <- matrix(nrow=num.comb,ncol=p)
   ind <- 1:p
   pointer <- p

   # calculate index combinations:
   for ( i in 1:num.comb )
   {
      comb[i,] <- ind
      ind[pointer] <- ind[pointer]+1
      if ( ind[pointer] > n )
      {
         while ( pointer > 1 )
         {
            pointer <- pointer-1
            if ( ind[pointer] < n-(p-pointer) )
            {
               ind[pointer] <- ind[pointer]+1
               for ( j in (pointer+1):p )
               {
                  ind[j] <- ind[pointer]+j-pointer
               }
               pointer <- p
               break
            }
         }
      }
   }

   # return results:
   return(comb)
}

################################################################################

# calculation of collinearity index (auxiliary function used in ident)
# ====================================================================

sysanal.collind <- function(sen.scaled)
{
   # -----------------------------------------------------------------------
   # This function calculates the collinearity index from a scaled 
   # sensitivity matrix.
   #
   # Arguments:
   # sen:       matrix of model sensitivities (scaled partial derivatives
   #            of model outcomes with respect to model parameters:
   #            delta.par/scale dy/dpar); the columns of sen refer to
   #            different model parameters, the rows to different model
   #            outcomes.
   #
   # Return Value:
   # collinearity index (real number).
   #
   #                                         Peter Reichert    Dec. 27, 2002
   # -----------------------------------------------------------------------

   # normalize sensitivity functions:
   num.par <- ncol(sen.scaled)
   norms <- numeric(num.par)
   for ( i in 1:num.par )
   {
      norms[i] <- sqrt( sen.scaled[,i] %*% sen.scaled[,i] )
   }
   sen.norm <- sen.scaled %*% diag( 1/norms )

   # calculate collinearity index:
   collind <- 1/sqrt(min(eigen( t(sen.norm) %*% sen.norm )$values))

   # return result:
   return(collind)
}

################################################################################

# calculation of identifiability measures
# =======================================

sysanal.ident <- function(sen,delta.par=0,scale=0,max.subset.size=0)
{
   # -----------------------------------------------------------------------
   # This function calculates a parameter sensitivity ranking and 
   # collinearity indices for a series of parameter combinations
   # based on linear sensitivity functions of a model, parameter
   # uncertainty ranges and scale factors of model results.
   # This function is a simplified version of the program "ident"
   # available at http://www.ident.eawag.ch
   #
   # Arguments:
   # sen:       matrix of model sensitivities (partial derivatives
   #            of model outcomes with respect to model parameters:
   #            dy/dpar); the columns of sen refer to different 
   #            model parameters, the rows to different model outcomes.
   # delta.par: model parameter uncertainty ranges (length equal to 
   #            the number of columns of sen); if zero, then all ranges
   #            are assumed to be unity.
   # scale:     scaling factors of model results (if zero, then all 
   #            scaling factors are assumed to be unity).
   #
   # Return Value:
   # List of delta.msqr, collind.
   #
   #                                         Peter Reichert    Dec. 27, 2002
   # -----------------------------------------------------------------------

   # determine number of model parameters:
   num.out <- nrow(sen)
   num.par <- ncol(sen)
   names.par <- colnames(sen)
   if ( length(names.par) != num.par ) names.par <- paste("par",1:num.par,sep="")
   if ( max.subset.size == 0 ) max.subset.size <- min(num.par,4)

   # apply parameter uncertainty ranges and scale factors if available:
   sen.scaled <- sen
   if ( length(delta.par) == num.par ) sen.scaled <- sen.scaled %*% diag(delta.par)
   if ( length(scale)     == num.out ) sen.scaled <- diag(1/scale) %*% sen.scaled

   # calculate sensitivity ranking:
   delta.msqr <- numeric(num.par)
   names(delta.msqr) <- names.par
   for ( i in 1:num.par )
   {
      delta.msqr[i] <- sqrt( t(sen.scaled[,i]) %*% sen.scaled[,i] ) / sqrt(num.out)
   }
   res <- list(delta.msqr=delta.msqr)

   if ( max.subset.size > 1 )
   {
      for ( i in 2:min(max.subset.size,num.par) )
      {
         ind <- sysanal.comb(num.par,i)
         collind <- numeric(nrow(ind))
         par.set <- matrix(nrow=nrow(ind),ncol=i)
         colnames(par.set) <- paste("par.",1:i,sep="")
         for ( j in 1:nrow(ind) )
         {
            collind[j] <- sysanal.collind(sen.scaled[,ind[j,]])
            for ( k in 1:i )
            {
               par.set[j,k] <- names.par[ind[j,k]]
            }
         }
         if ( nrow(par.set) > 1 )
         {
            ind.sorted <- order(collind)
            res[[paste("collind.",i,sep="")]] <- 
               data.frame(par.set[ind.sorted,],collind=collind[ind.sorted])
         }
         else
         {
            res[[paste("collind.",i,sep="")]] <- 
               data.frame(par.set,collind=collind)
         }
      }
   }

   # return results:
   return(res)
}

################################################################################

# sysanal.boxcox
# ==============

# purpose:
# Box-Cox transformation

# arguments:
# data:        data to be transformed
# lambda1:     value of parameter lambda1 (exponent) of Box-Cox transformation
# lambda2:     value of parameter lambda2 (offset) of Box-Cox transformation

# output:
# transformed data

# The prior for the parameter l1 of the Box-Cox transformation was chosen to be uniform
# in the interval [0, 1], and l2 was kept fixed at a value of 0

sysanal.boxcox <- function(data,lambda1=1,lambda2=1)
{
   if ( lambda1 == 0 )
   {
      return(ifelse(data>-lambda2,log(data+lambda2),NA))
   }
   else
   {
      return(ifelse(data>=-lambda2,((data+lambda2)^lambda1 - 1)/lambda1,NA))
   }
}

################################################################################

# sysanal.boxcox.deriv
# ====================

# purpose:
# calculate derivative of Box-Cox transformation

# arguments:
# data:        data to be transformed
# lambda1:     value of parameter lambda1 (exponent) of Box-Cox transformation
# lambda2:     value of parameter lambda2 (offset) of Box-Cox transformation

# output:
# derivative of Box-Cox transformation

sysanal.boxcox.deriv <- function(data,lambda1=1,lambda2=1)
{
   return(ifelse(data>-lambda2,(data+lambda2)^(lambda1 - 1),NA))
}

################################################################################

# sysanal.boxcox.inv
# ==================

# purpose:
# inverse Box-Cox transformation

# arguments:
# data:        data to be back-transformed
# lambda1:     value of parameter lambda1 (exponent) of Box-Cox transformation
# lambda2:     value of parameter lambda2 (offset) of Box-Cox transformation

# output:
# back-transformed data

sysanal.boxcox.inv <- function(data,lambda1=1,lambda2=1)
{
   if ( lambda1 == 0 )
   {
      return(exp(data)-lambda2)
   }
   else
   {
      return(ifelse(lambda1*data>-1,(lambda1*data+1)^(1/lambda1)-lambda2,
                                    -lambda2))
   }
}

################################################################################

# sysanal.trans.to.interval
# =========================

# purpose:
# transforms the real axis to an interval (default: unit interval)

# arguments:
# x:       data to be transformed
# min:     minimum of the interval (default: 0)
# max:     maximum of the interval (default: 1)

# output:
# transformed data

sysanal.trans.to.interval <- function(x,min=0,max=1)
{
   y <- 0.5*(min+max) + (max-min)/pi*atan(x)
   return(y)
}

################################################################################

# sysanal.trans.from.interval
# ===========================

# purpose:
# transforms an interval (default: unit interval) to the real axis

# arguments:
# y:       data to be transformed
# min:     minimum of the interval (default: 0)
# max:     maximum of the interval (default: 1)

# output:
# transformed data

sysanal.trans.from.interval <- function(y,min=0,max=1)
{
   x <- tan(0.5*pi*(2*y-max-min)/(max-min))
   return(x)
}

################################################################################
            
sysanal.trans.par.normal.tovec <- function(mean,sd,cor,trans=T,max.cor=0.5)
{
   n <- length(mean)
   
   par <- rep(NA,n*(n+3)/2)
   par[1:n] <- mean
   par[(n+1):(2*n)] <- sd
   if ( n > 1 )
   {
      k <- 2*n
      for ( i in 1:(n-1) ) 
      {
         par[k+1:(n-i)] <- cor[i,(i+1):n]
         k <- k + n - i
      }
   }
   names(par) <- c(names(mean),names(mean),rep("cor",n*(n-1)/2))
   
   if ( trans )
   {
      par[(n+1):(2*n)] <- log(par[(n+1):(2*n)])
      par[(2*n+1):(n*(n+3)/2)] <- 
         sysanal.trans.from.interval(par[(2*n+1):(n*(n+3)/2)],
                                     min=-max.cor,max=max.cor)
   }
   
   return(par)
}

sysanal.trans.par.normal.fromvec <- function(par,trans=T,max.cor=0.5)
{
   n <- (-3+sqrt(9+8*length(par)))/2
   
   if ( length(par) != n*(n+3)/2 )
   {
      cat("sysanal.trans.par.normal.fromvec:",
          "illegal length of parameter vector:",length(par),"\n")
      mean <- NA
      sd   <- NA
      cor  <- NA
   }
   else
   {
      if ( trans )
      {
         par[(n+1):(2*n)] <- exp(par[(n+1):(2*n)])
         par[(2*n+1):(n*(n+3)/2)] <- 
            sysanal.trans.to.interval(par[(2*n+1):(n*(n+3)/2)],
                                      min=-max.cor,max=max.cor)
      }

      mean <- par[1:n]
      sd   <- par[(n+1):(2*n)]
      cor  <- diag(rep(1,n),nrow=n)
      if ( n > 1 )
      {
         k <- 2*n
         for ( i in 1:(n-1) )
         {
            cor[i,(i+1):n] <- par[k+1:(n-i)]
            cor[(i+1):n,i] <- par[k+1:(n-i)]
            k <- k + n - i
         }
      }
      names(mean)   <- names(par)[1:n]
      names(sd)     <- names(par)[1:n]
      rownames(cor) <- names(par)[1:n]
      colnames(cor) <- names(par)[1:n]
   }

   return(list(mean=mean,sd=sd,cor=cor))
}            

################################################################################

# sysanal.gnls
# ============

# purpose:
# calculate the generalized least squares parameter estimates of a nonlinear model

# arguments:
# model.name:  name of the function representing the model
#              (note that this model requires the parameters to be specified
#              explicitly in the function headers)
# y.obs:       observed data corresponding to model output
# var.inv:     inverse variance-covariance matrix of the multivariate normal distribution
#              that characterizes the error model
# par:         named parameter vector with initial values; 
#              passed to the model as separate arguments
# x:           list of named inputs passed to the model
# ...          further optional arguments passed to nls

# output:
# list of:
# model.name:  name of the function representing the model
#              (cf. sysanal.model)
# par:         named parameter vector with parameter estimates 
# y.obs:       observed data corresponding to model output
# y.det:       deterministic model results corresponding to the estimated parameters
# resid:       residuals
# var.inv:     inverse variance-covariance matrix of the multivariate normal distribution
#              that characterizes the error model
# x:           input object passed as first argument to the model
# res.nls:     results from application of the nls function

sysanal.gnls <- function(model.name,y.obs,var.inv,par,x=list(),...)
{
   A <- chol(var.inv)
   add.args <- ""
   if ( names(x) > 0 )
   {
      add.args <- paste(",",paste(names(x),collapse=","),sep="")
   }
   model.call <- paste(model.name,"(",paste(names(par),collapse=","),
                       add.args,")",sep="")
   data <- x
   data$y.obs.trans <- A%*%y.obs
   res.nls <- nls(as.formula(paste("y.obs.trans ~ A%*%",model.call)),
                  data=data,
                  start=par,...)

   coef  <- coef(res.nls)
   args <- as.list(coef)
   if ( length(x) > 0 ) args <- c(args,x)
   y.det <- do.call(model.name,args=args)  # reevaluate to avoid need for
   resid <- y.obs - y.det                  # back transformation

   return(list(model.name = model.name,
               par        = coef,
               y.obs      = y.obs,
               y.det      = y.det,
               resid      = resid,
               var.inv    = var.inv,
               x          = x,
               res.nls    = res.nls))
}

################################################################################

# sysanal.gnls.diag
# =================

# purpose:
# calculate regression diagnostics for generalized least squares

# arguments:
# res.gnls:    output from sysanal.gnls
# par.inc:     increments of parameters used to approximate the derivatives

# output:
# list of:
# par.inc:     parameter increments used to calculate V
# V:           matrix of partial derivatives of model results with respect to
#              parameters at the estimated parameter values
# sd.par:      standard errors of the parameters assuming known error variances
# var.par:     estimated variance-covariance matrix of the parameters 
#              assuming known error variances
# corr.par:    estimated correlation matrix of the parameters
# sd.rel:      estimated correction factor in standard deviations of the error model
# sd.par.rel:  standard errors of the parameters when estimating a common factor
#              in error variances
# var.par.rel: estimated variance-covariance matrix of the parameters 
#              when estimating a common factor in error variances

sysanal.gnls.diag <- function(res.gnls,par.inc)
{
   # sensitivites:
   # -------------

   V <- sysanal.sens.loc.explicitpars(par        = res.gnls$par,
                                      model.name = res.gnls$model.name,
                                      par.inc    = par.inc,
                                      x          = res.gnls$x)

   # calculate variance-covariance matrix of the estimator for given error variance:
   # -------------------------------------------------------------------------------

   var.par  <- solve( t(V) %*% res.gnls$var %*% V )
   rownames(var.par) = names(res.gnls$par)
   colnames(var.par) = names(res.gnls$par)

   sd.par   <- sqrt(diag(var.par))
   names(sd.par) <- names(res.gnls$par)

   corr.par <- (1/sd.par) %*% t(1/sd.par) * var.par
   rownames(corr.par) = names(res.gnls$par)
   colnames(corr.par) = names(res.gnls$par)

   # variance-covariance matrix of the estimator for estimated error variance:
   # -------------------------------------------------------------------------

   sd.rel <- sqrt(as.numeric(
                     t(res.gnls$y.obs-res.gnls$y.det) %*% res.gnls$var.inv %*% (res.gnls$y.obs-res.gnls$y.det) /
                     (length(res.gnls$y.obs)-length(res.gnls$par))
                 ))
   var.par.rel <- var.par*sd.rel^2
   sd.par.rel  <- sd.par*sd.rel

   return(c(res.gnls,
          list(par.inc     = par.inc,
               V           = V,
               sd.par      = sd.par,
               var.par     = var.par,
               corr.par    = corr.par,
               sd.rel      = sd.rel,
               sd.par.rel  = sd.par.rel,
               var.par.rel = var.par.rel
         )))
}

################################################################################

# sysanal.gnls.test
# =================

# purpose:
# calculate test statistics for generalized least squares

# arguments:
# res.gnls:    output from sysanal.gnls

# output:
# list of:
# V:           matrix of partial derivatives of model results with respect to
#              parameters at the estimated parameter values
# sd.par:      standard errors of the parameters assuming known error variances
# var.par:     estimated variance-covariance matrix of the parameters 
#              assuming known error variances
# corr.par:    estimated correlation matrix of the parameters
# sd.rel:      estimated correction factor in standard deviations of the error model
# sd.par.rel:  standard errors of the parameters
#              when estimating a common factor in error variances
# var.par.rel: estimated variance-covariance matrix of the parameters 
#              when estimating a common factor in error variances

sysanal.gnls.test <- function(res.gnls,V,par,y.det)
{
   # chi2 test; to be compared with 
   # qchisq(1-alpha,length(y.obs)):
   # ------------------------------

   chi2 <- t(res.gnls$y.obs-y.det) %*% res.gnls$var.inv %*% (res.gnls$y.obs-y.det)

   # exact F test; to be compared with 
   # qf(1-alpha,length(par),length(y.obs)-length(par)):
   # --------------------------------------------------

   F.exact <- t(res.gnls$y.obs-y.det) %*% res.gnls$var.inv %*% V %*% 
              solve(t(V)%*%res.gnls$var.inv%*%V) %*% 
              t(V) %*% res.gnls$var.inv %*% (res.gnls$y.obs-y.det) /
                 t(res.gnls$y.obs-y.det) %*% 
                 ( res.gnls$var.inv -
                   res.gnls$var.inv %*% V %*% 
                   solve(t(V)%*%res.gnls$var.inv%*%V) %*% 
                   t(V) %*% res.gnls$var.inv ) %*% 
                 (res.gnls$y.obs-y.det) *
                 (length(res.gnls$y.obs)-length(res.gnls$par))/length(res.gnls$par)

   # F test for linearized model; to be compared with 
   # qf(1-alpha,length(par),length(y.obs)-length(par)):
   # --------------------------------------------------

   F.lin <- t(par-res.gnls$par) %*% t(V) %*% res.gnls$var.inv %*% V  %*% (par-res.gnls$par) /
               ( t(res.gnls$y.obs-res.gnls$y.det) %*% res.gnls$var.inv %*% (res.gnls$y.obs-res.gnls$y.det) ) *
               (length(res.gnls$y.obs)-length(res.gnls$par))/length(res.gnls$par)

   return(list(chi2        = chi2,
               F.exact     = F.exact,
               F.lin       = F.lin))
}

################################################################################

# sysanal.gnls.predict
# ====================

# purpose:
# calculate predictions based on generalized least squares regression results

# arguments:
# diag.gnls:   output from sysanal.gnls.diag
# newx:        new model input x at which confidence intervals are to be calculated
# level:       probability level of confidence intervals (default 0.95)
# var:         error variance-covariance matrix at new model input (only variances used)

# output:
# list of:
# x:           x values at which model results are calculated
# y.pred:      predicted model results
# confint:     confidence intervals of deterministic model results
# confint.rel: confidence intervals of deterministic model results
#              when estimating a common factor in error variances
# predint:     prediction intervals (only if error variance.covariance matrix is provided)
# predint.rel: prediction intervals (only if error variance.covariance matrix is provided)
#              when estimating a common factor in error variances

sysanal.gnls.predict <- function(diag.gnls,newx,level=0.95,var=NA)
{
   # prediction and sensitivites:
   # ----------------------------

   args <- as.list(diag.gnls$par)
   if ( length(newx) > 0 ) args <- c(args,newx)
   res.par <- do.call(diag.gnls$model.name,args=args)
   newV <- sysanal.sens.loc.explicitpars(par        = diag.gnls$par,
                                         model.name = diag.gnls$model.name,
                                         par.inc    = diag.gnls$par.inc,
                                         x          = newx)

   # variance-covariance matrix and standard deviaitons of model results:
   # --------------------------------------------------------------------

   var.y    <- newV %*% diag.gnls$var.par %*% t(newV)
   sd.y     <- sqrt(diag(var.y))
   sd.y.rel <- diag.gnls$sd.rel*sd.y

   # calculate confidence intervals:
   # -------------------------------

   df          <- length(diag.gnls$y.obs)-length(diag.gnls$par)
   confint     <- sysanal.confint(est=res.par,sd=sd.y,df=NA,level=level)
   confint.rel <- sysanal.confint(est=res.par,sd=sd.y.rel,df=df,level=level)

   res <- list(x           = newx,
               y.pred      = res.par,
               confint     = confint,
               confint.rel = confint.rel)

   if ( is.matrix(var) )
   {
      sd.y.err     <- sqrt(sd.y^2+diag(var))
      sd.y.rel.err <- sqrt(sd.y.rel^2+diag(var)*diag.gnls$sd.rel^2)
      predint     <- sysanal.confint(est=res.par,sd=sd.y.err,df=NA,level=level)
      predint.rel <- sysanal.confint(est=res.par,sd=sd.y.rel.err,df=df,level=level)
      res <- c(res,list(predint=predint,predint.rel=predint.rel))
   }

   return(res)
}

################################################################################

# sysanal.confint
# ===============

# purpose:
# calculate confidence intervals based on the Student t distribution
# or on the normal distribution (if df=NA)

# arguments:
# est:         vector of point estimates
# sd:          vector of estimated standard deviations
# df:          degrees of freedom
# level:       probability level of confidence intervals (default 0.95)

# output:
# table of lower and upper bounds of confidence intervals

sysanal.confint <- function(est,sd,df,level)
{
   alpha <- 1- level
   confint <- matrix(NA,nrow=length(est),ncol=2)
   colnames(confint) <- c(paste(100*alpha/2,"%",sep=""),paste(100*(1-alpha/2),"%",sep=""))
   rownames(confint) <- names(est)
   if ( is.na(df) )
   {
      fact <- qnorm(1-alpha/2)
   }
   else
   {
      fact <- qt(1-alpha/2,df)
   }
   confint[,1] <- est - fact*sd
   confint[,2] <- est + fact*sd
   return(confint)
}

################################################################################

# sysanal.resid.diag
# ==================

# purpose:
# plot residual diagnostics plots

# arguments:
# obs:         vector of observations
# calc:        vector of calculated results

# output:
# diagnostic plots

sysanal.resid.diag <- function(obs,calc,header="")
{
   # calculate range of measurements and normalized residuals:
   # ---------------------------------------------------------

   obs.min        <- min(obs)
   obs.max        <- max(obs)
   calc.min       <- min(calc)
   calc.max       <- max(calc)
   resid.norm     <- (obs-calc)/sd(obs-calc)
   resid.norm.abs <- abs(resid.norm)
   resid.max      <- max(abs(resid.norm.abs))
   resid.lim      <- 1.1*resid.max*c(-1,1)
   marker         <- 19

   # divide plot area into four panels:
   # ----------------------------------

   par.def <- par(no.readonly=TRUE)
   par(mfrow=c(2,2),xaxs="i",yaxs="i",mar=c(4.5,4,3,1),oma=c(0,0,2,0))

   # plot sequence of residuals:
   # ---------------------------

   plot(resid.norm,main="Sequence of Residuals",
        ylim=resid.lim,pch=marker,cex=0.8,
        xlab="Data Points",ylab="Normalized Residuals")
   lines(c(0,length(resid.norm)),c(0,0))

   # plot residuals as function of predicted values:
   # -----------------------------------------------

   plot(calc,resid.norm,main="Residuals vs. Predicted",
        xlim=c(calc.min,calc.max),ylim=resid.lim,
        xlab="Predicted",ylab="Normalized Residuals",pch=marker,cex=0.8)
   lines(c(calc.min,calc.max),c(0,0))
   res.lm <- lm(resid.norm.abs ~ calc)
   x.new <- c(calc.min,calc.max)
   y.new <- predict(res.lm,newdata=data.frame(calc=x.new))
   lines(x.new, y.new)
   lines(x.new,-y.new)

   # plot histogram of residuals:
   # ----------------------------

   hist(resid.norm,freq=FALSE,main="Hist. of Residuals",
        xlim=resid.lim,
        xlab="Normalized Residuals",ylab="Density")
   lines(seq(-3,3,by=0.1),dnorm(seq(-3,3,by=0.1)))
   lines(resid.lim,c(0,0))

   # normal quantile plot:
   # ---------------------

   lim <- max(resid.max,qnorm(1-0.5/length(obs))+0.1)
   qqnorm(resid.norm,main="Sample vs. Normal Quant.",
          xlab="Normal Quantiles",ylab="Sample Quantiles",
          pch=marker,cex=0.8,
          xlim=1.1*lim*c(-1,1),ylim=1.1*lim*c(-1,1))
   lines(c(-10,10),c(-10,10))

   # reset plot attributes:
   # ----------------------

   mtext(header,side=3,outer=T,adj=0.5,cex=1.2)

   par(par.def)
}

################################################################################

# sysanal.resid.diag.boxcox
# =========================

# purpose:
# plot residual diagnostics plots for given Box-Cox transformation parameters

# arguments:
# obs:         vector of observations
# calc:        vector of calculated results
# lambda1:     value of parameter lambda1 (exponent) of Box-Cox transformation
# lambda2:     value of parameter lambda2 (offset) of Box-Cox transformation

# output:
# diagnostic plots

sysanal.resid.diag.boxcox <- function(obs,calc,lambda1=1,lambda2=1)
{
   sysanal.resid.diag(sysanal.boxcox(obs,lambda1,lambda2),
                      sysanal.boxcox(calc,lambda1,lambda2))
}

################################################################################

# sysanal.hessian
# ===============

sysanal.hessian <- function(fn,par,par.inc=NA,...)
{
   if ( is.na(par.inc[1]) ) par.inc <- 0.01*par
   n <- length(par)
   h <- matrix(NA,nrow=n,ncol=n)
   colnames(h) <- names(par)
   rownames(h) <- names(par)
   trace <- matrix(NA,nrow=2*n^2+1,ncol=n+1)
   colnames(trace) <- c(names(par),"fn")
   minimum <- TRUE
   maximum <- TRUE
   par.0 <- par
   res.0 <- fn(par.0,...)
   num.eval <- 1
   trace[1,] <- c(par.0,res.0)
   row.trace <- 1
   for ( i in 1:n )
   {
      for ( j in 1:i )
      {
         if ( i==j )
         {
            par.inc.tmp <- par.inc
            counter <- 0
            while (TRUE)
            {
               par.1 <- par.0
               par.1[i] <- par.1[i] + par.inc.tmp[i]
               res.1 <- fn(par.1,...)
               num.eval <- num.eval+1
               if ( !is.na(res.1) )
               { 
                  if ( res.1 > res.0 ) maximum <- FALSE
                  if ( res.1 < res.0 ) minimum <- FALSE
                  par.2 <- par.0
                  par.2[i] <- par.2[i] - par.inc.tmp[i]
                  res.2 <- fn(par.2,...)
                  num.eval <- num.eval+1
                  if ( !is.na(res.2) )
                  {
                     if ( res.2 > res.0 ) maximum <- FALSE
                     if ( res.2 < res.0 ) minimum <- FALSE
                     h[i,i] <- (res.1 - 2*res.0 + res.2)/par.inc.tmp[i]^2
                     trace[row.trace+1,] <- c(par.1,res.1)
                     trace[row.trace+2,] <- c(par.2,res.2)
                     row.trace <- row.trace + 2
                     break
                  }
               }
               counter <- counter + 1
               if ( counter > 5 ) stop("sysanal.hessian: unable to calculate hessian")
               par.inc.tmp <- 0.5*par.inc.tmp
            }
         }
         else
         {
            par.inc.tmp <- par.inc
            counter <- 0
            while (TRUE)
            {
               par.1 <- par.0
               par.1[i] <- par.1[i] + par.inc.tmp[i]
               par.1[j] <- par.1[j] + par.inc.tmp[j]
               res.1 <- fn(par.1,...)
               num.eval <- num.eval + 1
               if ( !is.na(res.1) )
               {
                  if ( res.1 > res.0 ) maximum <- FALSE
                  if ( res.1 < res.0 ) minimum <- FALSE
                  par.2 <- par.0
                  par.2[i] <- par.2[i] + par.inc.tmp[i]
                  par.2[j] <- par.2[j] - par.inc.tmp[j]
                  res.2 <- fn(par.2,...)
                  num.eval <- num.eval + 1
                  if ( !is.na(res.2) )
                  {
                     if ( res.2 > res.0 ) maximum <- FALSE
                     if ( res.2 < res.0 ) minimum <- FALSE
                     par.3 <- par.0
                     par.3[i] <- par.3[i] - par.inc.tmp[i]
                     par.3[j] <- par.3[j] + par.inc.tmp[j]
                     res.3 <- fn(par.3,...)
                     num.eval <- num.eval + 1
                     if ( !is.na(res.3) )
                     {
                        if ( res.3 > res.0 ) maximum <- FALSE
                        if ( res.3 < res.0 ) minimum <- FALSE
                        par.4 <- par.0
                        par.4[i] <- par.4[i] - par.inc.tmp[i]
                        par.4[j] <- par.4[j] - par.inc.tmp[j]
                        res.4 <- fn(par.4,...)
                        num.eval <- num.eval + 1
                        if ( !is.na(res.4) )
                        {
                           if ( res.4 > res.0 ) maximum <- FALSE
                           if ( res.4 < res.0 ) minimum <- FALSE
                           h[i,j] <- (res.1 - res.2 - res.3 + res.4)/
                                     (4*par.inc.tmp[i]*par.inc.tmp[j])
                           h[j,i] <- h[i,j]
                           trace[row.trace+1,] <- c(par.1,res.1)
                           trace[row.trace+2,] <- c(par.2,res.2)
                           trace[row.trace+3,] <- c(par.3,res.3)
                           trace[row.trace+4,] <- c(par.4,res.4)
                           row.trace <- row.trace + 4
                           break
                        }
                     }
                  }
               }
               counter <- counter + 1
               if ( counter > 5 ) stop("sysanal.hessian: unable to calculate hessian")
               par.inc.tmp <- 0.5*par.inc.tmp
            }
         } 
      }
   }
   return(list(h        = h,
               minimum  = minimum,
               maximum  = maximum,
               num.eval = num.eval,
               trace    = trace))
}      

##############################################################################
#                                                                            #
# sysanal.smooth                                                             #
# --------------                                                             #
#                                                                            #
# Function for smoothing data points and estimating the derivative of the    #
# smoothed curve by local quadratic or optionally local linear regression.   #
# Local regression is implemented by using a Gaussian distribution of        #
# of weights centered a the point at which the smoothed curve has to be      #
# evaluated.                                                                 #
# The smoothing parameter is the standard deviation of the Gaussian weights. #
#                                                                            #
# Peter Reichert 05.09.2008 , last modification 01.05.2009                   #
#                                                                            #
##############################################################################

# Call:      sysanal.smooth(x,y,sigma,newx=NA,fac.extrap=1,method="quadratic")
# -----
#
# Input:
# ------
#
# x          vector of x-coordinates of data points to be smoothed
# y          vector of y-coordinates of data points to be smoothed
#            (x and y must be of the same length)
# sigma      standard deviation of Gaussian distribution used as weights for
#            local quadratic or local linear regression
# newx       optional vector of x-coordinates at which smoothed results and
#            derivatives are to be calculated (if not specified, results
#            are provided at the same locations as there is data available)
# fac.extrap calculate smoothed value only if data is available within
#            fac.extrap*sigma (default is 1)
# method     "quadratic" (default) indicates local quadratic regression, 
#            otherwise local linear regression is used
#
# Output:
# -------
#
# Data frame of
# x          x-coordinates at which smoothed results and derivatives are 
#            available
# y          smoothed results at the locations x
# ydot       derivatives of smoothed results at the locations x

sysanal.smooth <- function(x,y,sigma,newx=NA,fac.extrap=1,method="quadratic")
{
   # calculate x and y vectors for available data:
   ind <- !is.na(y)
   x.data <- x[ind]
   y.data <- y[ind]

   # check consistency of input:
   if ( length(x.data) < 1 )
   {
      stop("*** error in sysanal.smooth: no data available ***")
   }
   if ( length(x.data) != length(y.data) ) 
   {
      stop("*** error in sysanal.smooth: length of x and y different ***")
   }
   if ( ! sigma > 0 ) 
   {
      stop("*** error in sysanal.smooth: sigma is not positive ***")
   }
   
   # select x values for output:
   if ( is.na(newx[1]) ) newx <- x
   
   # set up ysmooth and ydot vectors:
   n <- length(newx)
   ysmooth <- rep(NA,n)
   ydot    <- rep(NA,n)
   
   # calclate smoothed values and derivatives:
   for ( i in 1:n )
   {
      # get indices of data within +/- 2*sigma:
      ind.extrap <- x.data >= newx[i]-fac.extrap*sigma & 
                    x.data <= newx[i]+fac.extrap*sigma
      num.extrap <- sum(ifelse(ind.extrap,1,0))
      
      # calculate smoothed value only if data is available within +/- 2*sigma:
      if ( num.extrap > 0 ) 
      {
         # still use data within a 5 times larger interval
         # to calculate the smoothed value:
         fac.use <- 4*max(1,fac.extrap)
         ind.use <- x.data >= newx[i]-fac.use*sigma & 
                    x.data <= newx[i]+fac.use*sigma
         x1  <- x.data[ind.use]-newx[i]
         x2  <- (x.data[ind.use]-newx[i])^2
         num.use <- sum(ifelse(ind.use,1,0))
         if ( num.use == 1 )  # use value
         {
            ysmooth[i] <- y.data[ind.use][1]
            ydot[i]    <- 0
         }
         else
         {
            if ( num.use == 2 )  # use weighted mean
            {
               weights <- dnorm(x.data[ind.use],mean=newx[i],sd=sigma)
               weights <- weights/sum(weights)
               ysmooth[i] <- weights[1]*y.data[ind.use][1] + 
                             weights[2]*y.data[ind.use][2]
               if ( x.data[ind.use][2] != x.data[ind.use][1] )
               {
                  ydot[i]    <- (y.data[ind.use][2]-y.data[ind.use][1])/
                                (x.data[ind.use][2]-x.data[ind.use][1])
               }
            }
            else
            {
               if ( method != "quadratic" | num.use == 3 ) # use local linear
               {                                           # regression
                  res.lm     <- lm(y.data[ind.use] ~ x1,
                                   weights=dnorm(x.data[ind.use],
                                                 mean=newx[i],sd=sigma))
                  ysmooth[i] <- coef(res.lm)[1]
                  ydot[i]    <- coef(res.lm)[2]
               }
               else  # use local quadratic regression
               {
                  res.lm     <- lm(y.data[ind.use] ~ x1 + x2,
                                   weights=dnorm(x.data[ind.use],
                                                 mean=newx[i],sd=sigma))
                  ysmooth[i] <- coef(res.lm)[1]
                  ydot[i]    <- coef(res.lm)[2]
               }
            }
         }
      }
   }
   
   # return data frame:
   return(data.frame(x=newx,y=ysmooth,ydot=ydot))
}


##############################################################################
#                                                                            #
# sysanal.smooth_fun                                                         #
# ------------------                                                         #
#                                                                            #
# Function for smoothing piecewise linear functions.                         #
# The smoothing parameter is the standard deviation of the Gaussian weights. #
#                                                                            #
# Peter Reichert 16.12.2008 , last modification 01.05.2009                   #
#                                                                            #
##############################################################################

# Call:      sysanal.smooth_fun(data,z,sigma,newz=NA,newx=NA,fac.extrap=1,
#                               method="quadratic")
# -----
#
# Input:
# ------
#
# data       list of matrices specifying piecewise linear functions:
#            for each function, the independent variable x must be provided 
#            in the first column, the dependent variable y in the second
# z          vector of z-coordinates corresponding to the functions
# sigma      standard deviation of Gaussian distribution used as weights for
#            local quadratic regression in z-coordinates (see sysanal.smooth)
# newz       optional vector of z-coordinates at which smoothed functions
#            are to be calculated (if not specified, results are provided
#            at the same locations as there is data available)
# newx       optional vector of x-coordinates  at which smoothed function 
#            values are to be calculated
# fac.extrap calculate smoothed value only if data is available within
#            fac.extrap*sigma (default is 1)
# method     "quadratic" (default) indicates local quadratic regression, 
#            otherwise local linear regression is used
#
# Output:
# -------
#
# List of data frames with columns
# x          x-coordinates at which smoothed results are available
# y          smoothed results at the locations x

sysanal.smooth_fun <- function(data,z,sigma,newz=NA,newx=NA,fac.extrap=1,
                               method="quadratic")
{
   # check consistency of input:
   n <- length(data)
   if ( length(z) != n )
   { 
      stop("error in sysanal.smooth_fun: not same number of locations as functions")
   }
   
   # select z and x values for output:
   if ( is.na(newz[1]) ) newz <- z
   if ( is.na(newx[1]) )
   {
      x <- numeric(0)
      for ( j in 1:n ) x <- c(x,data[[j]][,1])
      newx <- sort(unique(x))
   }

   # interpolate input functions to selected x values:    
   y <- matrix(nrow=length(newx),ncol=n)
   for ( j in 1:n )
   {
      y[,j] <- approx(x=data[[j]][,1],y=data[[j]][,2],xout=newx)$y
   }

   # set up result data structure:
   res <- list()
   for ( j in 1:length(newz) )
   {
      res[[j]] <- data.frame(x=newx,y=rep(NA,length(newx)))
   }
   names(res) <- newz

   # calculate results:
   for ( i in 1:length(newx) )
   {
      res.smooth <- sysanal.smooth(x=z,y=y[i,],sigma=sigma,newx=newz,
                                   fac.extrap=fac.extrap,method=method)$y
      for ( j in 1:length(newz) )
      {
         res[[j]]$y[i] <- res.smooth[j] 
      }     
   }
   
   # return results:
   return(res)
}

############################################################################

# ======================================================== #
# Numerical Integration of Ordinary Differential Equations #
# ======================================================== #


# This library contains a didactical implementation of a numerical integrator
# of a system of ordinary differential equations.
# Note that this implementation is intended to demonstrate how such 
# techniques can be implemented in R. It does not represent the state of the 
# art of numerical integration of such differential equations.


# created and maintained by 
# Peter Reichert
# EAWAG
# Duebendorf
# Switzerland
# reichert@eawag.ch


# First version: Dec.  22, 2002
# Last revision: April 09, 2006


# Overview of functions
# =====================

# sysanal.ode:    numerical integration of deterministic ordinary differential 
#                 equations


# =========================================================================== #


# numerical integration of ordinary differential equations
# ========================================================

sysanal.ode <- function(rhs,x.ini,par,t.out,dt.max=NA,algorithm="euler",...)
{
   # -----------------------------------------------------------------------
   # This solver for a system of ordindary differential equations
   # was written for didactical purposes. It does not represent 
   # the state of the art in numerical integration techniques
   # but should rather demonstrate how the simplest integration
   # technique can be implemented for relatively simple use.
   # Check the package "odesolve" available from http://www.r-project.org
   # for professional solvers for ordinary differential equations.
   #
   # Arguments:
   # rhs:       function returning the right hand side of the 
   #            system of differential equations as a function
   #            of the arguments x (state variables), t (time),
   #            and par (model parameters).
   # x.ini:     start values of the state variables.
   # par:       model parameters (transferred to rhs).
   # t.out:     set of points in time at which the solution 
   #            should be calculated; the solution at the first
   #            value in t is set to x.ini, the solution at subsequent
   #            values is calculated by numerical integration.
   # dt.max:    maximum time step; if the difference between points
   #            in time at which output has to be provided, t, is
   #            larger than dt.max, then this output interval is 
   #            divided into step of at most dt.max internally to
   #            improve the accuracy of the solution at the next 
   #            output point.
   # algorithm: right now, the only options are "euler" (explicit first order
   #            Euler algorithm) and "euler.2.order" (explicit second order
   #            Euler algorithm).
   #
   # Return Value:
   # matrix with results for state variables (columns) at all output
   # time steps t.out (rows).
   #
   #                                        Peter Reichert    Dec.  22, 2002
   #                                        last modification April 09, 2006
   # -----------------------------------------------------------------------
 
   # determine number of equations and number of time steps:
   num.eq <- length(x.ini)
   steps <- length(t.out)

   # define and initialize result matrix:
   x <- matrix(nrow=steps,ncol=num.eq)
   colnames(x) <- names(x.ini)
   rownames(x) <- t.out
   x[1,] <- x.ini

   # perform integration:
   if ( algorithm == "euler" )
   {
      for ( i in 2:steps )
      {
         if ( is.na(dt.max) || dt.max >= t.out[i]-t.out[i-1] )
         {
            # output interval <= dt.max (or dt.max not available):
            # perform a single step to the next output time:
            x[i,] <- x[i-1,] + (t.out[i]-t.out[i-1])*rhs(x[i-1,],t.out[i-1],par,...)
         }
         else
         {
            # output interval > dt.max:
            # perform multiple steps of maximum size dt.max:
            x[i,] <- x[i-1,]
            steps.internal <- ceiling((t.out[i]-t.out[i-1])/dt.max)
            dt <- (t.out[i]-t.out[i-1])/steps.internal
            for ( j in 1:steps.internal )
            {
               t.current <- t.out[i-1] + (j-1)*dt
               x[i,] <- x[i,] + dt*rhs(x[i,],t.current,par,...)
            }
         }
      }
   }
   else
   {
      if ( algorithm == "euler.2.order" )
      {
         for ( i in 2:steps )
         {
            if ( is.na(dt.max) || dt.max >= t.out[i]-t.out[i-1] )
            {
               # output interval <= dt.max (or dt.max not available):
               # perform a single step to the next output time:
               t.mid <- 0.5*(t.out[i-1]+t.out[i])
               x.mid <- x[i-1,] + 0.5*(t.out[i]-t.out[i-1])*rhs(x[i-1,],t.out[i-1],par,...)
               x[i,] <- x[i-1,] + (t.out[i]-t.out[i-1])*rhs(x.mid,t.mid,par,...)
            }
            else
            {
               # output interval > dt.max:
               # perform multiple steps of maximum size dt.max:
               x[i,] <- x[i-1,]
               steps.internal <- ceiling((t.out[i]-t.out[i-1])/dt.max)
               dt <- (t.out[i]-t.out[i-1])/steps.internal
               for ( j in 1:steps.internal )
               {
                  t.current <- t.out[i-1] + (j-1)*dt
                  t.mid <- t.current + 0.5*dt
                  x.mid <- x[i,] + 0.5*dt*rhs(x[i,],t.current,par,...)
                  x[i,] <- x[i,] + dt*rhs(x.mid,t.mid,par,...)
               }
            }
         }
      }
      else
      {
         stop(paste("sysanal.ode: algorithm \"",
                    algorithm,
                    "\" not implemented",
                    sep=""))
      }
   }

   # return result matrix: 
   return(x)
}


# =========================================================================== #


sysanal.Var.B.M.L <- function(psi,L,i=NA,j=NA)
{
   L.decoded <- L
   L.encoded <- L
   if ( is.vector(L) )
   {
      L.decoded <- sysanal.decode(L)
   }
   else
   {
      L.encoded <- rownames(L)
   }
   n <- length(L.encoded)

   i.local <- i
   if ( is.na(i.local[1]) ) i.local <- 1:n
   j.local <- j # in predictions it will just be n2+i i.e. a scalar like 956
   if ( is.na(j.local[1]) ) j.local <- 1:n
   Var <- matrix(0,nrow=length(i.local),ncol=length(j.local))
   vars <- unique(L.decoded$var)
   beta <- 1/psi["corrlen"]^2
   if ( is.na(beta) ) beta <- 0
   for ( var in vars )
   {
      ind <- which(L.decoded$var==var)
      ind.i <- match(ind,i.local)
      ind.i <- ind.i[!is.na(ind.i)] # all the positions of L1, e.g. c(1,2,3...
      ind.j <- match(ind,j.local)
      ind.j <- ind.j[!is.na(ind.j)]
      if ( length(ind.i)>0 & length(ind.j)>0 )
      {
         name.sd.B <- paste("sd.B",var,sep="_")
         var.B <- psi[name.sd.B]^2
         if ( is.na(var.B) ) var.B <- 0
         dist <- rep(1,length(ind.i)) %*% t(L.decoded$val[j.local[ind.j]]) -
                 L.decoded$val[i.local[ind.i]] %*% t(rep(1,length(ind.j)))
# rep(1,length(ind.i)): repeats 1 n1 time
# L.decoded$val[j.local[ind.j]]: in predictions it accesses to the n2+i value of L, e.g. 1912
# L.decoded$val[i.local[ind.i]]: accesses to to all the values of L1 e.g. c(2,4,6...
# rep(1,length(ind.j)): in predictions is 1   
         
# dist: difference of a col vec with the n2+i value of L and another with all the values of L1 (in pred)
         Var[ind.i,ind.j] <- var.B*exp(-beta*dist^2) #Covariance structure
# Var: in pred simply indicates the covariance of every value in L2 with all the values in L1
      }
   }
   rownames(Var) <- L.encoded[i.local]
   colnames(Var) <- L.encoded[j.local]
   return(Var)
}


sysanal.sd.Eps.L <- function(xi,L)
{
   L.decoded <- L
   L.encoded <- L
   if ( is.vector(L) )
   {
      L.decoded <- sysanal.decode(L)
   }
   else
   {
      L.encoded <- rownames(L)
   }
   n <- length(L.encoded)

   sd <- rep(NA,n)
   vars <- unique(L.decoded$var)
   for ( var in vars )
   {
      name.sd.Eps <- paste("sd.Eps",var,sep="_")
      sd.Eps <- xi[name.sd.Eps]
      if ( is.na(sd.Eps) ) stop(paste("*** parameter",
                                      name.sd.Eps,
                                      "not found in sd.Eps.L"))
      ind <- which(L.decoded$var==var)
      sd[ind] <- rep(sd.Eps,length(ind))
   }
   names(sd) <- L.encoded   
   return(sd)
}


sysanal.loglikeli.bias <- function(par,model,L,y.obs,Var.B,sd.Eps,
                                   lambda1=1,lambda2=1)
{
#    if (any (par<0))
#    {return(-Inf)}###
     
  
  # decode layout definition: 

   L.decoded <- L
   L.encoded <- L
   if ( is.vector(L) )
   {
      L.decoded <- sysanal.decode(L)
   }
   else
   {
      L.encoded <- rownames(L)
   }
   
   # calculate results of deterministic model:

   y.calc        <- model(par,L.encoded) #changed!! L.decoded
   
   # transform results and observations:
   
   y.calc.trans  <- sysanal.boxcox(y.calc,lambda1,lambda2)
   y.obs.trans   <- sysanal.boxcox(y.obs,lambda1,lambda2)
   boxcox.deriv  <- sysanal.boxcox.deriv(y.obs,lambda1,lambda2)
   
   # evaluate likelihood function:
   
   Sigma.B       <- Var.B(par,L.decoded)
   Sigma.Eps     <- diag(sd.Eps(par,L.decoded)^2)
   
   Sum.Sigma     <- Sigma.B + Sigma.Eps
   
   Sum.Sigma.inv <- solve(Sum.Sigma)

   log.det.Sum.Sigma <- determinant(Sum.Sigma,logarithm=TRUE)
   if ( log.det.Sum.Sigma$sign < 0 ) stop("determinant Sigma.Eps+Sigma.B < 0")
   
   loglikeli <- - 0.5 * length(L.encoded) * log(2*pi) -
                  0.5 * as.numeric(log.det.Sum.Sigma$modulus) -
                  0.5 * t(y.obs.trans-y.calc.trans) %*% (Sum.Sigma.inv ) %*% 
                        (y.obs.trans-y.calc.trans) +
                  sum(log(abs(boxcox.deriv)))

   return(loglikeli)
}
   

sysanal.loglikeli.bias.blockdesign <- function(par,
                                               model,
                                               L,
                                               y.obs,
                                               Var.B,
                                               sd.Eps,
                                               lambda1 = 1,
                                               lambda2 = 1)
{
   # decode layout definition:

   L.decoded <- L
   L.encoded <- L
   if ( is.vector(L) )
   {
      L.decoded <- sysanal.decode(L)
   }
   else
   {
      L.encoded <- rownames(L)
   }
   n <- length(L.encoded)

   # calculate results of deterministic model:

   y.calc <- model(par,L.decoded)
   
   # transform results and observations:
   
   y.calc.trans  <- sysanal.boxcox(y.calc,lambda1,lambda2)
   y.obs.trans   <- sysanal.boxcox(y.obs,lambda1,lambda2)
   boxcox.deriv  <- sysanal.boxcox.deriv(y.obs,lambda1,lambda2)
   
   # evaluate likelihood function:
   
   loglikeli <- 0
   vars <- unique(L.decoded$var)
   for ( var in vars )
   {
      ind <- which(L.decoded$var==var)

      Sigma.B   <- Var.B(par,L.decoded[ind,])
      Sigma.Eps <- diag(sd.Eps(par,L.decoded[ind,])^2)

      Sum.Sigma <- Sigma.Eps + Sigma.B

      Sum.Sigma.inv <- solve(Sum.Sigma)
   
      log.det.Sum.Sigma <- determinant(Sum.Sigma,logarithm=TRUE)
      if ( log.det.Sum.Sigma$sign < 0 ) stop("determinant Sigma.Eps+Sigma.B < 0")

      loglikeli <- loglikeli -
                     0.5 * length(ind) * log(2*pi) -
                     0.5 * as.numeric(log.det.Sum.Sigma$modulus) -
                     0.5 * t(y.obs.trans[ind]-y.calc.trans[ind]) %*%
                           Sum.Sigma.inv %*%
                           (y.obs.trans[ind]-y.calc.trans[ind]) +
                  sum(log(abs(boxcox.deriv[ind])))
   }

   return(loglikeli)
}
   



sysanal.predict.bias.cond.blockdesign <- function(par,model,L1,y.obs,
                                                  Var.B,sd.Eps,
                                                  L2=NA,y.calc=NA,
                                                  lambda1=1,lambda2=1,
                                                  ...)
{
   # decode likelihood definitions:
   
   L1.decoded <- L1
   L1.encoded <- L1
   if ( is.vector(L1) )
   {
      L1.decoded <- sysanal.decode(L1)
   }
   else
   {
      L1.encoded <- rownames(L1)
   }
   n1 <- length(L1.encoded)

   L.decoded <- L1.decoded
   L.encoded <- L1.encoded
   n2 <- 0
   n <- n1
   
   L2.available <- TRUE
   if ( is.vector(L2) ) { if ( is.na(L2[1]) ) L2.available <- FALSE }   
   if ( L2.available )
   { 
      L2.decoded <- L2
      L2.encoded <- L2
      if ( is.vector(L2) )
      {
         L2.decoded <- sysanal.decode(L2)
      }
      else
      {
         L2.encoded <- rownames(L2)
      }
      n2 <- length(L2.encoded)
   
      L.decoded <- rbind(L1.decoded,L2.decoded)
      L.encoded <- c(L1.encoded,L2.encoded)
      n <- n1 + n2
   }

   # calculate results of deterministic model:
   
   if ( is.na(y.calc[1]) )
   {   
      y.calc    <- model(par,L.decoded,...)
   }
   else
   {
      if ( length(y.calc) != n )
      {
         cat("*** y.calc is not of correct length:",length(y.calc),
             "instead of",n,"\n")
         return(NA)
      }
   }
   y.calc.L1 <- y.calc[1:n1]
   y.calc.L2 <- NA
   if ( n2 > 0 )
   {
      y.calc.L2    <- y.calc[n1+(1:n2)]
   }

   # transform results and observations:
   
   y.calc.L1.trans  <- sysanal.boxcox(y.calc.L1,lambda1,lambda2)
   y.calc.L2.trans  <- sysanal.boxcox(y.calc.L2,lambda1,lambda2)
   y.obs.trans      <- sysanal.boxcox(y.obs,lambda1,lambda2)

   # set up arrays and calculate results:

   B.mean.L1 <- rep(NA,n1)
   B.var.L1  <- rep(NA,n1)
   Y.mean.L1 <- rep(NA,n1)
   Y.var.L1  <- rep(NA,n1)
   B.mean.L2 <- rep(NA,max(1,n2))
   B.var.L2  <- rep(NA,max(1,n2))
   Y.mean.L2 <- rep(NA,max(1,n2))
   Y.var.L2  <- rep(NA,max(1,n2))
   list.Sum.Sigma.L1.inv <- list()
   list.Sigma.B.L1       <- list()
   list.Sigma.B.L2L1     <- list()
   vars <- unique(L.decoded$var)
   for ( var in vars )
   {
      ind1 <- which(L1.decoded$var==var)
      
      Sigma.B.L1   <- Var.B(par,L1.decoded[ind1,])
      sigma.Eps.L1 <- sd.Eps(par,L1.decoded[ind1,])
      Sigma.Eps.L1 <- diag(sigma.Eps.L1^2)

      Sum.Sigma.L1 <- Sigma.B.L1 + Sigma.Eps.L1
   
      Sum.Sigma.L1.inv <- solve(Sum.Sigma.L1)
      
      list.Sigma.B.L1[[var]]       <- Sigma.B.L1
      list.Sum.Sigma.L1.inv[[var]] <- Sum.Sigma.L1.inv

      B.var.L1[ind1]  <- diag(Sigma.B.L1 %*% Sum.Sigma.L1.inv %*% Sigma.Eps.L1)
      B.mean.L1[ind1] <- Sigma.B.L1 %*% Sum.Sigma.L1.inv %*% 
                         (y.obs.trans[ind1]-y.calc.L1.trans[ind1])

      # calculate predictions for layout 1:

      Y.mean.L1[ind1] <- y.calc.L1.trans[ind1] + 
                         Sigma.B.L1 %*% Sum.Sigma.L1.inv %*% 
                            ( y.obs.trans[ind1] - y.calc.L1.trans[ind1] )
      Y.var.L1[ind1]  <- diag(Sigma.B.L1 %*% Sum.Sigma.L1.inv %*% Sigma.Eps.L1) +
                         diag(Sigma.Eps.L1)
                         
      # calculate predictions for layout 2: 
     
      if ( n2 > 0 )
      {
         ind2 <- which(L2.decoded$var==var)
         if ( length(ind2) > 0 )
         {
            list.Sigma.B.L2L1[[var]] <- matrix(nrow=length(ind2),ncol=length(ind1))
      
            for ( i in 1:length(ind2) )
            {
               v <- as.vector(Var.B(par,L.decoded,i=ind1,j=n1+ind2[i]))
               list.Sigma.B.L2L1[[var]][i,] <- v
               Sigma.B.L2.i.i   <- as.numeric(Var.B(par,L2.decoded,ind2[i],ind2[i]))

               Y.mean.L2[ind2[i]] <- y.calc.L2.trans[ind2[i]] +  
                                     t(v) %*% Sum.Sigma.L1.inv %*% 
                                     ( y.obs.trans[ind1] - y.calc.L1.trans[ind1] )
               Y.var.L2[ind2[i]]  <- Sigma.B.L2.i.i + 
                                     sd.Eps(par,L2.decoded[ind2[i],])^2 - 
                                     t(v) %*% Sum.Sigma.L1.inv %*% v 
                         
               B.mean.L2[ind2[i]] <- t(v) %*% Sum.Sigma.L1.inv %*% 
                                     ( y.obs.trans[ind1] - y.calc.L1.trans[ind1] )
               B.var.L2[ind2[i]]  <- Sigma.B.L2.i.i - 
                                     t(v) %*% Sum.Sigma.L1.inv %*% v 
            }
         }
      }
   }
   names(B.mean.L1) <- L1.encoded
   names(B.var.L1)  <- L1.encoded
   names(Y.mean.L1) <- L1.encoded
   names(Y.var.L1)  <- L1.encoded
   names(B.mean.L2) <- L2.encoded
   names(B.var.L2)  <- L2.encoded
   names(Y.mean.L2) <- L2.encoded
   names(Y.var.L2)  <- L2.encoded

   return(list(y.calc.L1        = y.calc.L1.trans,
               B.mean.L1        = B.mean.L1,
               B.var.L1         = B.var.L1,
               Y.mean.L1        = Y.mean.L1,
               Y.var.L1         = Y.var.L1,
               y.calc.L2        = y.calc.L2.trans,
               B.mean.L2        = B.mean.L2,
               B.var.L2         = B.var.L2,
               Y.mean.L2        = Y.mean.L2,
               Y.var.L2         = Y.var.L2,
               Sum.Sigma.L1.inv = list.Sum.Sigma.L1.inv,
               Sigma.B.L1       = list.Sigma.B.L1,
               Sigma.B.L2L1     = list.Sigma.B.L2L1))
}


sysanal.predict.bias <- function(parsamp.L1,model,L1,y.obs,
                                 predict.bias.cond,
                                 probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),
                                 L2=NA,y.calc=NA,
                                 lambda1=1,lambda2=1,
                                 ...)
{
   # decode likelihood definitions:
   
   L1.decoded <- L1
   L1.encoded <- L1
   if ( is.vector(L1) )
   {
      L1.decoded <- sysanal.decode(L1)
   }
   else
   {
      L1.encoded <- rownames(L1)
   }
   n1 <- length(L1.encoded)

   L.decoded <- L1.decoded
   L.encoded <- L1.encoded
   n2 <- 0
   n <- n1
   
   L2.available <- TRUE
   if ( is.vector(L2) ) { if ( is.na(L2[1]) ) L2.available <- FALSE }   
   if ( L2.available )
   { 
      L2.decoded <- L2
      L2.encoded <- L2
      if ( is.vector(L2) )
      {
         L2.decoded <- sysanal.decode(L2)
      }
      else
      {
         L2.encoded <- rownames(L2)
      }
      n2 <- length(L2.encoded)
   
      L.decoded <- rbind(L1.decoded,L2.decoded)
      L.encoded <- c(L1.encoded,L2.encoded)
      n <- n1 + n2
   }
   
   # initialize result arrays:

   y.L1.samp <- matrix(nrow=nrow(parsamp.L1),ncol=n1)
   colnames(y.L1.samp) <- L1.encoded
   B.L1.samp <- y.L1.samp
   Y.L1.samp <- y.L1.samp
   neg.var.B.L1 <- rep(0,n1)
   names(neg.var.B.L1) <- L1.encoded
   neg.var.Y.L1 <- neg.var.B.L1
   
   y.L2.samp    <- NA
   B.L2.samp    <- NA
   Y.L2.samp    <- NA
   neg.var.B.L2 <- NA
   neg.var.Y.L2 <- NA
   if ( n2 > 0 )
   {
      y.L2.samp <- matrix(nrow=nrow(parsamp.L1),ncol=n2)
      colnames(y.L2.samp) <- L2.encoded
      B.L2.samp <- y.L2.samp
      Y.L2.samp <- y.L2.samp
      neg.var.B.L2 <- rep(0,n2)
      names(neg.var.B.L2) <- L2.encoded
      neg.var.Y.L2 <- neg.var.B.L2
   }
   
   # calculate samples:
   
   par.old <- rep(NA,ncol(parsamp.L1))
   num.eval <- 0
   for ( j in 1:nrow(parsamp.L1) )
   {
      par <- parsamp.L1[j,]
      
      if ( j==1 | sum((par-par.old)^2) != 0 ) # do not reevaluate if parameter
      {                                       # values stayed the same
         par.old <- par
         num.eval <- num.eval + 1      
         res <- predict.bias.cond(par     = par,
                                  model   = model,
                                  L1      = L1.decoded,
                                  y.obs   = y.obs,
                                  L2      = L2,
                                  y.calc  = y.calc,
                                  lambda1 = lambda1,
                                  lambda2 = lambda2,
                                  ...)
      }
      y.L1.samp[j,] <- res$y.calc.L1;
#       print(j);# print(summary(y.L1.samp[j,]))
      for ( i in 1:n1 )
      {
         var <- res$Y.var.L1[i]
         if ( var < 0 ) 
         { 
            cat("* Warning: negative variance of Y:",var,"at",L1.encoded[i],"\n") 
            var <- 0; neg.var.Y.L1[i] <- neg.var.Y.L1[i]+1
         }
         Y.L1.samp[j,i] <- rnorm(1,res$Y.mean.L1[i],sqrt(var))
      }
      if ( ! is.na(res$B.mean.L1[1]) )
      {
         for ( i in 1:n1 )
         {   
            var <- res$B.var.L1[i]
            if ( var < 0 ) 
            { 
               cat("* Warning: negative variance of B:",var,"at",L1.encoded[i],"\n") 
               var <- 0; neg.var.B.L1[i] <- neg.var.B.L1[i]+1 
            }
            B.L1.samp[j,i] <- rnorm(1,res$B.mean.L1[i],sqrt(var))   
         }
      }
      if ( n2 > 0 )
      {
         y.L2.samp[j,] <- res$y.calc.L2
         for ( i in 1:n2 )
         {
            var <- res$Y.var.L2[i]
            if ( var < 0 ) 
            { 
               cat("* Warning: negative variance of Y:",var,"at",L2.encoded[i],"\n") 
               var <- 0; neg.var.Y.L2[i] <- neg.var.Y.L2[i]+1 
            }
            Y.L2.samp[j,i] <- rnorm(1,res$Y.mean.L2[i],sqrt(var))
         }
         if ( ! is.na(res$B.mean.L2[1]) )
         {
            for ( i in 1:n2 )
            {   
               var <- res$B.var.L2[i]
               if ( var < 0 ) 
               { 
                  cat("* Warning: negative variance of B:",var,"at",L2.encoded[i],"\n") 
                  var <- 0; neg.var.B.L2[i] <- neg.var.B.L2[i]+1 
               }
               B.L2.samp[j,i] <- rnorm(1,res$B.mean.L2[i],sqrt(var))
            }
         }      
      }
   }
   neg.var.B.L1 <- neg.var.B.L1/n1
   neg.var.Y.L1 <- neg.var.Y.L1/n1
   if ( n2 > 0 ) 
   {
      neg.var.B.L2 <- neg.var.B.L2/n2
      neg.var.Y.L2 <- neg.var.Y.L2/n2
   }
   
   # derive quantiles:
   
   y.L1.quant <- matrix(nrow=length(probs),ncol=n1)
   colnames(y.L1.quant) <- L1.encoded
   rownames(y.L1.quant) <- probs
   B.L1.quant <- y.L1.quant
   yplusB.L1.quant <- y.L1.quant
   Y.L1.quant <- y.L1.quant
   for ( i in 1:n1 )
   {
      y.L1.quant[,i]      <- quantile(y.L1.samp[,i],probs=probs,na.rm=TRUE)
      B.L1.quant[,i]      <- quantile(B.L1.samp[,i],probs=probs,na.rm=TRUE)
      yplusB.L1.quant[,i] <- quantile(y.L1.samp[,i]+B.L1.samp[,i],
                                                    probs=probs,na.rm=TRUE)
      Y.L1.quant[,i]      <- quantile(Y.L1.samp[,i],probs=probs,na.rm=TRUE)
   }
   
   y.L2.quant <- NA
   B.L2.quant <- NA
   yplusB.L2.quant <- NA
   Y.L2.quant <- NA
   if ( n2 > 0 )
   {
      y.L2.quant <- matrix(nrow=length(probs),ncol=n2)
      colnames(y.L2.quant) <- L2
      rownames(y.L2.quant) <- probs
      B.L2.quant <- y.L2.quant
      yplusB.L2.quant <- y.L2.quant
      Y.L2.quant <- y.L2.quant
      for ( i in 1:n2 )
      {
         y.L2.quant[,i]      <- quantile(y.L2.samp[,i],probs=probs,na.rm=TRUE)
         B.L2.quant[,i]      <- quantile(B.L2.samp[,i],probs=probs,na.rm=TRUE)
         yplusB.L2.quant[,i] <- quantile(y.L2.samp[,i]+B.L2.samp[,i],
                                                       probs=probs,na.rm=TRUE)
         Y.L2.quant[,i]      <- quantile(Y.L2.samp[,i],probs=probs,na.rm=TRUE)
      }
   }

   return(list(y.L1.samp       = y.L1.samp,
               B.L1.samp       = B.L1.samp,
               Y.L1.samp       = Y.L1.samp,
               y.L2.samp       = y.L2.samp,
               B.L2.samp       = B.L2.samp,
               Y.L2.samp       = Y.L2.samp,
               y.L1.quant      = y.L1.quant,
               B.L1.quant      = B.L1.quant,
               yplusB.L1.quant = yplusB.L1.quant,
               Y.L1.quant      = Y.L1.quant,
               y.L2.quant      = y.L2.quant,
               B.L2.quant      = B.L2.quant,
               yplusB.L2.quant = yplusB.L2.quant,
               Y.L2.quant      = Y.L2.quant,
               neg.var.B.L1    = neg.var.B.L1,
               neg.var.Y.L1    = neg.var.Y.L1,
               neg.var.B.L2    = neg.var.B.L2,
               neg.var.Y.L2    = neg.var.Y.L2))
}



sysanal.predict.bias.lin <- function(par.mean.L1,
                                     par.var.L1,
                                     model,
                                     L1,
                                     y.obs,
                                     predict.bias.cond,
                                     par.inc = 0.01*par.mean.L1,
                                     L2      = NA,
                                     y.calc  = NA,
                                     V.y     = NA,
                                     lambda1 = 1,
                                     lambda2 = 1,
                                     logpar  = F,
                                     ...)
{
   # decode likelihood definitions:

   L1.decoded <- L1
   L1.encoded <- L1
   if ( is.vector(L1) )
   {
      L1.decoded <- sysanal.decode(L1)
   }
   else
   {
      L1.encoded <- rownames(L1)
   }
   n1 <- length(L1.encoded)

   L.decoded <- L1.decoded
   L.encoded <- L1.encoded
   n2 <- 0
   n <- n1
   
   L2.available <- TRUE
   if ( is.vector(L2) ) { if ( is.na(L2[1]) ) L2.available <- FALSE }   
   if ( L2.available )
   { 
      L2.decoded <- L2
      L2.encoded <- L2
      if ( is.vector(L2) )
      {
         L2.decoded <- sysanal.decode(L2)
      }
      else
      {
         L2.encoded <- rownames(L2)
      }
      n2 <- length(L2.encoded)
   
      L.decoded <- rbind(L1.decoded,L2.decoded)
      L.encoded <- c(L1.encoded,L2.encoded)
      n <- n1 + n2
   }
   
   # calculate transformed data:
   
   y.obs.trans <- sysanal.boxcox(y.obs,lambda1,lambda2)
   
   # calculate solution and conditional predictions:

   par <- par.mean.L1
   par.local <- par
   if ( logpar ) par.local <- exp(par)
   res.par <- predict.bias.cond(par     = par.local,
                                model   = model,
                                L1      = L1.decoded,
                                y.obs   = y.obs,
                                L2      = L2.decoded,
                                y.calc  = y.calc,
                                lambda1 = lambda1,
                                lambda2 = lambda2,
                                ...)
   boxcox.deriv <- sysanal.boxcox.deriv(sysanal.boxcox.inv(c(res.par$y.calc.L1,
                                                             res.par$y.calc.L2),
                                                           lambda1,lambda2),
                                        lambda1,lambda2)
                            
   # calculate jacobian of bias and of model plus bias:

   V.yplusB <- matrix(NA,nrow=n,ncol=length(par))
   rownames(V.yplusB) <- L.encoded
   colnames(V.yplusB) <- names(par)
   V.y.test <- V.yplusB
   V.B      <- V.yplusB
   V.Eps    <- matrix(NA,nrow=n1,ncol=length(par))
   rownames(V.Eps) <- L1.encoded
   colnames(V.Eps) <- names(par)
   if ( !is.matrix(V.y) & !is.data.frame(V.y) ) # calculate numerical deriv 
   {                                            # of all expected val.
      for( j in 1:length(par) )
      {
         par.j <- par
         par.j[j] <- par.j[j] + par.inc[j]
         if ( logpar ) par.j <- exp(par.j)
         res.par.j <- predict.bias.cond(par     = par.j,
                                        model   = model,
                                        L1      = L1.decoded,
                                        y.obs   = y.obs,
                                        L2      = L2,
                                        y.calc  = NA,
                                        lambda1 = lambda1,
                                        lambda2 = lambda2,
                                        ...)
         
         V.y.test[,j] <- ( c(res.par.j$y.calc.L1,res.par.j$y.calc.L2) - 
                           c(res.par$y.calc.L1,res.par$y.calc.L2) ) /
                         par.inc[j]
         V.yplusB[,j] <- ( c(res.par.j$y.calc.L1 + res.par.j$B.mean.L1,
                             res.par.j$y.calc.L2 + res.par.j$B.mean.L2) - 
                           c(res.par$y.calc.L1   + res.par$B.mean.L1,
                             res.par$y.calc.L2   + res.par$B.mean.L2) ) /
                         par.inc[j]
         V.B[,j]      <- ( c(res.par.j$B.mean.L1,res.par.j$B.mean.L2) - 
                           c(res.par$B.mean.L1,res.par$B.mean.L2) ) /
                         par.inc[j]
         V.Eps[,j]    <- ( (y.obs.trans - res.par.j$y.calc.L1 - res.par.j$B.mean.L1) -
                           (y.obs.trans - res.par$y.calc.L1 - res.par$B.mean.L1) ) /
                         par.inc[j]
      }
      #write.table(V.y.test,"V.y.calculated.dat",row.names=T,col.names=NA,sep="\t")
   }
   else
   {
      # extend Y.y to other parameters, assume no sensitivity
      
      if ( nrow(V.y) != n ) stop("incorrect number of rows in V.y")
      V.y.local <- matrix(0,nrow=n,ncol=length(par))
      colnames(V.y.local) = names(par)
      rownames(V.y.local) = L.encoded
      for ( j in colnames(V.y) )  # copy entries from V.y
      {
         V.y.local[,j] <- V.y[,j]*boxcox.deriv
      }
      #write.table(V.y.local,"V.y.derived.dat",row.names=T,col.names=NA,sep="\t")
      
      # calculate derivatives of matrix factors:

      if ( ! is.list(res.par$Sum.Sigma.L1.inv) )
      {
         A <- res.par$Sigma.B.L1   %*% res.par$Sum.Sigma.L1.inv
         if ( n2 > 0 ) B <- res.par$Sigma.B.L2L1 %*% res.par$Sum.Sigma.L1.inv
         for( j in 1:length(par) )
         {
            par.j <- par
            par.j[j] <- par.j[j] + par.inc[j]
            if ( logpar ) par.j <- exp(par.j)
            res.par.j <- predict.bias.cond(par     = par.j,
                                           model   = model,
                                           L1      = L1.decoded,
                                           y.obs   = y.obs,
                                           L2      = L2,
                                           y.calc  = c(res.par$y.calc.L1,res.par$y.calc.L2),
                                           lambda1 = lambda1,
                                           lambda2 = lambda2,
                                           ...)
            # note that y.calc is not correct for these parameter values.
            # this makes all results that depend on y.calc incorrect.
            # this is no problem as we will only use results that do not depend on 
            # y.calc
         
            dA <- ( res.par.j$Sigma.B.L1   %*% res.par.j$Sum.Sigma.L1.inv - A ) /
                  par.inc[j] 
            a1 <- dA %*% ( y.obs.trans - res.par$y.calc.L1 )
            b1 <- A %*% ( - V.y.local[1:n1,j] )            
            V.yplusB[1:n1,j] <- V.y.local[1:n1,j] + a1 + b1
            V.B[1:n1,j]      <- a1 + b1
            V.Eps[,j]        <- - V.yplusB[1:n1,j]
            if ( n2 > 0 )
            {      
               dB <- ( res.par.j$Sigma.B.L2L1 %*% res.par.j$Sum.Sigma.L1.inv - B ) /
                     par.inc[j]
               a2 <- dB %*% ( y.obs.trans - res.par$y.calc.L1 )
               b2 <- B %*% ( - V.y.local[1:n1,j] )             
               V.yplusB[n1+(1:n2),j] <- V.y.local[n1+(1:n2),j] + a2 + b2
               V.B[n1+(1:n2),j]      <- a2 + b2
            }
         }
      }
      else   # blockdesign
      {
         for( j in 1:length(par) )
         {
            par.j <- par
            par.j[j] <- par.j[j] + par.inc[j]
            if ( logpar ) par.j <- exp(par.j)
            res.par.j <- predict.bias.cond(par     = par.j,
                                           model   = model,
                                           L1      = L1.decoded,
                                           y.obs   = y.obs,
                                           L2      = L2,
                                           y.calc  = c(res.par$y.calc.L1,res.par$y.calc.L2),
                                           lambda1 = lambda1,
                                           lambda2 = lambda2,
                                           ...)
            # note that y.calc is not correct for these parameter values.
            # this makes all results that depend on y.calc incorrect.
            # this is no problem as we will only use results that do not depend on 
            # y.calc
         
            vars <- unique(L.decoded$var)
            for ( var in vars )
            {
               ind1 <- which(L1.decoded$var==var)
               A <- res.par$Sigma.B.L1[[var]]   %*% res.par$Sum.Sigma.L1.inv[[var]]
               dA <- ( res.par.j$Sigma.B.L1[[var]]   
                          %*% res.par.j$Sum.Sigma.L1.inv[[var]] - A ) /
                     par.inc[j] 
               a1 <- dA %*% ( y.obs.trans[ind1] - res.par$y.calc.L1[ind1] )
               b1 <- A %*% ( - V.y.local[ind1,j] )
               V.yplusB[ind1,j] <- V.y.local[ind1,j] + a1 + b1
               V.B[ind1,j]      <- a1 + b1
               V.Eps[ind1,j]    <- - V.yplusB[ind1,j]
               ind2 <- which(L2.decoded$var==var)
               if ( length(ind2) > 0 )
               {
                  B <- res.par$Sigma.B.L2L1[[var]] %*% res.par$Sum.Sigma.L1.inv[[var]]
                  dB <- ( res.par.j$Sigma.B.L2L1[[var]] 
                             %*% res.par.j$Sum.Sigma.L1.inv[[var]] - B ) /
                        par.inc[j]
                  a2 <- dB %*% ( y.obs.trans[ind1] - res.par$y.calc.L1[ind1] )
                  b2 <- B %*% ( - V.y.local[ind1,j] )
                  V.yplusB[n1+ind2,j]  <- V.y.local[n1+ind2,j] + a2 + b2
                  V.B[n1+ind2,j]       <- a2 + b2
               }
            }
         }
      }
   }

   # calculate results:

   y.cond      <- c(res.par$y.calc.L1,res.par$y.calc.L2)
   B.mean.cond <- c(res.par$B.mean.L1,res.par$B.mean.L2)
   Y.mean.cond <- c(res.par$Y.mean.L1,res.par$Y.mean.L2)
   B.var.cond  <- c(res.par$B.var.L1, res.par$B.var.L2)
   Y.var.cond  <- c(res.par$Y.var.L1, res.par$Y.var.L2)

   y.mean      <- y.cond
   B.mean      <- B.mean.cond
   B.var       <- diag( V.B %*% par.var.L1 %*% t(V.B) ) + B.var.cond
   Eps.mean    <- y.obs.trans - res.par$y.calc.L1 - res.par$B.mean.L1
   Eps.var     <- diag( V.Eps %*% par.var.L1 %*% t(V.Eps) ) + res.par$B.var.L1
   yplusB.mean <- y.cond + B.mean.cond
   yplusB.var  <- diag( V.yplusB %*% par.var.L1 %*% t(V.yplusB) ) + B.var.cond
   Y.mean      <- yplusB.mean
   Y.var       <- diag( V.yplusB %*% par.var.L1 %*% t(V.yplusB) ) + Y.var.cond
   
   names(y.mean)      <- L.encoded
   names(B.mean)      <- L.encoded
   names(B.var)       <- L.encoded
   names(Eps.mean)    <- L1.encoded
   names(Eps.var)     <- L1.encoded
   names(yplusB.mean) <- L.encoded
   names(yplusB.var)  <- L.encoded
   names(Y.mean)      <- L.encoded
   names(Y.var)       <- L.encoded
  
   return(list(y.mean      = y.mean,
               B.mean      = B.mean,
               B.var       = B.var,
               Eps.mean    = Eps.mean,
               Eps.var     = Eps.var,
               yplusB.mean = yplusB.mean,
               yplusB.var  = yplusB.var,
               Y.mean      = Y.mean,
               Y.var       = Y.var))
}

# didactical emulator (preliminary version, added May 30, 2012):
# ====================

evaluate <- function(x, ...) UseMethod("evaluate")

sysanal.emulator.create <- function(inp.design,
                                    res.design,
                                    par,
                                    sd,
                                    lambda,
                                    alpha    = 2,
                                    sd.smooth = 0,
                                    pri.mean = sysanal.emulator.pri.mean.lin,
                                    pri.var  = sysanal.emulator.pri.var.normal)
{
  #sysanal.package("fpc")
  
  n.design <- length(res.design)
  dim.inp  <- length(lambda)
  if ( is.vector(inp.design) ) 
  {
    inp.design <- matrix(inp.design,nrow=n.design)
  }
  if ( nrow(inp.design) != n.design ) stop("number of outputs must be equal to number of inputs")
  if ( ncol(inp.design) != dim.inp )  stop("dimension of input must be equal to length of lambda")
  
  emu <- list()
  emu$inp.design  <- inp.design
  emu$res.design  <- res.design
  emu$par         <- par
  emu$sd          <- sd
  emu$lambda      <- lambda
  emu$alpha       <- alpha
  emu$sd.smooth   <- sd.smooth
  emu$pri.mean    <- pri.mean
  emu$pri.var     <- pri.var
  emu$pri.var.num <- pri.var(inp=inp.design,sd=sd,lambda=lambda,alpha=alpha,sd.smooth=sd.smooth)
  emu$delta       <- res.design - pri.mean(inp=inp.design,par=par)
  emu$inv         <- solve(emu$pri.var.num)
  #emu$inv         <- solvecov(emu$pri.var.num)$inv
  emu$v           <- emu$inv %*% emu$delta
  class(emu) <- "emulator.GASP"
  
  return(emu)
}


sysanal.emulator.pri.mean.lin <- function(inp,par)
{
  if ( is.vector(inp) ) inp <- matrix(inp,ncol=length(par)-1)
  par <- as.numeric(par)
  
  n.inp   <- nrow(inp)
  dim.inp <- ncol(inp)
  dim.par <- length(par)
  
  if ( dim.par != dim.inp+1 ) stop("number of parameters must be dim.inp+1")
  
  res <- par[1] + inp %*% par[-1]
  rownames(res) <- 1:n.inp
  colnames(res) <- "y"
  
  return(res)
}


sysanal.emulator.pri.var.normal <- function(inp,sd,lambda,alpha=2,sd.smooth=0,rows=NA,cols=NA)
{
  if ( is.vector(inp) ) inp <- matrix(inp,ncol=length(lambda))
  n.inp   <- nrow(inp)
  dim.inp <- ncol(inp)
  if ( length(lambda) != dim.inp ) stop("length of lambda must be dim.inp")
  if ( is.na(rows[1]) ) rows <- 1:n.inp
  if ( is.na(cols[1]) ) cols <- 1:n.inp
  var <- matrix(0,nrow=length(rows),ncol=length(cols))
  
  s <- matrix(0,nrow=length(rows),ncol=length(cols))
  for ( k in 1:dim.inp )
  {
    s <- s + (abs(inp[rows,k]%o%rep(1,length(cols))-rep(1,length(rows))%o%inp[cols,k])/
      lambda[k])^alpha
  }
  var <- sd*sd*exp(-s)
  var <- var+diag(sd.smooth*sd.smooth,nrow=n.inp)[rows,cols]
  
  colnames(var) <- cols
  rownames(var) <- rows
  return(var)
}


evaluate.emulator.GASP <- function(emulator,inp)
{
  if ( is.vector(inp) ) inp <- matrix(inp,ncol=length(emulator$lambda))
  n.inp <- nrow(inp)
  n.design <- nrow(emulator$inp.design)
  y   <- rep(NA,n.inp)
  var <- rep(NA,n.inp)
  for ( i in 1:n.inp )
  {
    pri.mean <- emulator$pri.mean(inp      = matrix(inp[i,],nrow=1),par=emulator$par)
    k        <- emulator$pri.var(inp       = rbind(emulator$inp.design,matrix(inp[i,],nrow=1)),
                                 sd        = emulator$sd,
                                 lambda    = emulator$lambda,
                                 alpha     = emulator$alpha,
                                 sd.smooth = emulator$sd.smooth,
                                 rows      = n.design+1,
                                 cols      = 1:n.design)
    K        <- emulator$pri.var(inp       = matrix(inp[i,],nrow=1),
                                 sd        = emulator$sd,
                                 lambda    = emulator$lambda,
                                 alpha     = emulator$alpha,
                                 sd.smooth = emulator$sd.smooth)
    y[i]   <- pri.mean + k %*% emulator$v
    var[i] <- K - k %*% emulator$inv %*% t(k)
  }
  return(list(inp=inp,y=y,var=var))
}

# ---------------------------------------------------------------------------
# Input-dependent bias likelihood (actually also constant or output dependent)
# ---------------------------------------------------------------------------

sysanal.loglikeli.bias.inp <- function(par,model,L,y.obs,Var.Bs, inp=rep(0,length(y.obs)), sd.Eps,
																			 par.tr,par.fix=NULL, ...)
{
  if (any (par[!is.na(par)]<0))
  {return(-Inf)} 
  
  # decode layout definition: 
  
  L.decoded <- L
  L.encoded <- L
  if ( is.vector(L) )
  {
    L.decoded <- sysanal.decode(L)
  }
  else
  {
    L.encoded <- rownames(L)
  }
  
  # calculate results of deterministic model:

  y.calc        <- model(par,L.encoded,...)  
  
  # evaluate likelihood function:
  par.comb      <- c(par,par.fix)
  
  
  if (sum(inp)>0) inp = inp[(par.comb["Del.Max"]+1-round(par.comb["Delta"])):(length(inp)-round(par.comb["Delta"]))]

  # transform results and observations:
  
  if(is.na(par.tr["alpha"])) 
  {
  	y.calc.trans  <- sysanal.boxcox(y.calc,par.tr["l1"],par.tr["l2"])
  	y.obs.trans   <- sysanal.boxcox(y.obs,par.tr["l1"],par.tr["l2"])
  	boxcox.deriv  <- sysanal.boxcox.deriv(y.obs,par.tr["l1"],par.tr["l2"])
  }  else{
  	y.calc.trans  <- sysanal.logsinh(y.calc,par.tr["alpha"],par.tr["beta"])
  	y.obs.trans   <- sysanal.logsinh(y.obs,par.tr["alpha"],par.tr["beta"])
  	boxcox.deriv  <- sysanal.logsinh.deriv(y.obs,par.tr["alpha"],par.tr["beta"])
  }


  Sigma.Bs      <- Var.Bs(par.comb,L.decoded, inp) # B cov

  Sigma.Eps     <- diag(sd.Eps(par.comb,L.decoded)^2)
  
  Sum.Sigma     <-  Sigma.Bs +Sigma.Eps 

  Sum.Sigma.inv <- solve(Sum.Sigma)
  
  
  log.det.Sum.Sigma <- determinant(Sum.Sigma,logarithm=TRUE)

  if ( log.det.Sum.Sigma$sign < 0 ) 
  {warning("determinant Sigma.Eps+Sigma.B < 0") 
   loglikeli=-Inf}
  else{
  loglikeli <- - 0.5 * length(L.encoded) * log(2*pi) -
    0.5 * as.numeric(log.det.Sum.Sigma$modulus) -
    0.5 * t(y.obs.trans-y.calc.trans) %*% (Sum.Sigma.inv ) %*% 
    (y.obs.trans-y.calc.trans) +
    sum(log(abs(boxcox.deriv)))}
 
  return(loglikeli)
}


# ---------------------------------------------------------------------------
# Covariance matrix for the rain- and self-dependent bias term
# ---------------------------------------------------------------------------

sysanal.Var.Bs <- function(psi,L, inp, i=NA,j=NA, var_prev=NA, compl.var=NA)
{
	L.decoded <- L
	L.encoded <- L
	if ( is.vector(L) )
	{
		L.decoded <- sysanal.decode(L)
	}
	else
	{
		L.encoded <- rownames(L)
	}
	n <- length(L.encoded)
	
	i.local <- i
	if ( is.na(i.local[1]) ) i.local <- 1:n
	j.local <- j
	if ( is.na(j.local[1]) ) j.local <- 1:n
	Var.s <- matrix(0,nrow=length(i.local),ncol=length(j.local))
	vars <- unique(L.decoded$var)
	
	as.var <- array(dim=c(length(i.local)/length(vars),1))
	Dt=L.decoded$val[2]-L.decoded$val[1]
	
	beta <- 1/psi["corrlen"]^1 
	if ( is.na(beta) ) beta <- 0
	exp.dec  <- exp(-2*beta*(Dt))
	
	for ( var in vars ) 
	{
		ind <- which(L.decoded$var==var)
		ind.i <- match(ind,i.local)
		ind.i <- ind.i[!is.na(ind.i)]
		ind.j <- match(ind,j.local)
		ind.j <- ind.j[!is.na(ind.j)]
		
		name.sd.B <- paste("sd.B",var,sep="_")
		var.B <- psi[name.sd.B]^2
		if ( is.na(var.B) ) var.B <- 0
		
		name.ks.B <- paste("ks",var,sep="_")
		ks.B      <- psi[name.ks.B] ; if (is.na(ks.B)) ks.B <- 0
		
		if ( length(ind.i)>0 & length(ind.j)>0 )
		{
			
			if ( is.na(j) ) # Calibration
			{
				if (is.na(var_prev))
				{prim.var =var.B
				 }else{
				 prim.var=var_prev
				 }
				
				as.var[1]  = var.B +(prim.var-var.B-(ks.B*inp[1])^2)*exp.dec+(ks.B*inp[1])^2
				
				for(i in 2:length(ind.i))
				{as.var[i] = var.B +(as.var[i-1]-var.B-(ks.B*inp[i])^2)*exp.dec+(ks.B*inp[i])^2}
				
				SBs       = array(dim=c(length(as.var),length(as.var)))
				lower.ind = lower.tri(SBs,T)
				autoreg   = exp(-beta*abs(outer(L.decoded$val,L.decoded$val,"-")))
				SBs       = as.vector(as.var)*autoreg
				SBs[lower.ind] = 0
				SBs       = SBs+t(SBs)
				diag(SBs) = as.var
				
				Var.s[ind.i,ind.j] <- SBs #Covariance structure of the slow bias
			}
			else 
			{ #SBL1,2
				dist                <- rep(1,length(ind.i)) %*% t(L.decoded$val[j.local[ind.j]]) -
					L.decoded$val[i.local[ind.i]] %*% t(rep(1,length(ind.j)))
				
				Var.s[ind.i,ind.j]  <- (compl.var[ind.i])*exp(-beta*dist)  
				
			} 
			
			rownames(Var.s) <- L.encoded[i.local] #i.local and j.local are the same thing in L1
			colnames(Var.s) <- L.encoded[j.local]  
		}
		else # SBL2
		{Var.s[1,1]  <- (var.B+(var_prev-var.B-(ks.B*inp[i])^2)*exp.dec+ (ks.B*inp[j])^2)}
	}
	return(Var.s)
}



# ---------------------------------------------------------------------------
# Compute Smoothing Distributions
# ---------------------------------------------------------------------------


sysanal.predict.inp.bias.L1 <- function(par,model,L1,y.obs,
                                      Var.Bs, inp,sd.Eps,
                                      L2,
                                        y.calc=NA, par.tr,  ...)
{
  # decode likelihood definitions:
  
  L1.decoded <- L1
  L1.encoded <- L1
  if ( is.vector(L1) )
  {
    L1.decoded <- sysanal.decode(L1)
  }
  else
  {
    L1.encoded <- rownames(L1)
  }
  n1 <- length(L1.encoded)
  
  L.decoded <- L1.decoded
  L.encoded <- L1.encoded
  n2 <- 0
  n <- n1
  
  # calculate results of deterministic model:
  
  if ( is.na(y.calc[1]) )
  {   
    if ( !is.na(L2[1]) ) y.calc    <- model(par,c(L1,L2),...) else y.calc    <- model(par,c(L1),...)
  }
  else
  {
    if ( length(y.calc) != n )
    {
      cat("*** y.calc is not of correct length:",length(y.calc),
          "instead of",n,"\n")
      return(NA)
    }
  }
  
  y.calc.L1 <- y.calc[1:n1] 
  y.calc.L2 <- NA
  
  # transform results and observations:
  
  if(is.na(par.tr["alpha"])) 
  {
  	y.calc.L1.trans  <- sysanal.boxcox(y.calc.L1,par.tr["l1"],par.tr["l2"])
  	if ( !is.na(L2[1]) )	y.calc.L2        <- sysanal.boxcox(y.calc[(n1+1):length(y.calc)],par.tr["l1"],par.tr["l2"])
  	y.obs.trans      <- sysanal.boxcox(y.obs,par.tr["l1"],par.tr["l2"])
  }
  else{
  	y.calc.L1.trans  <- sysanal.logsinh(y.calc.L1,par.tr["alpha"],par.tr["beta"])
  	if ( !is.na(L2[1]) )	y.calc.L2        <- sysanal.logsinh(y.calc[(n1+1):length(y.calc)],par.tr["alpha"],par.tr["beta"])
  	y.obs.trans      <- sysanal.logsinh(y.obs,par.tr["alpha"],par.tr["beta"])
  }
  
  # calculate predictions for layout 1:   
  
  inp = inp[(par["Del.Max"]+1-round(par["Delta"])):(length(inp)-round(par["Delta"]))]
  
  Sigma.Bs.L1       <- Var.Bs(par,L1.decoded,inp) #slow B cov
  Sigma.Eps.L1      <- diag(sd.Eps(par,L1.decoded)^2)
  
  Sum.Sigma.L1       <- Sigma.Bs.L1 + Sigma.Eps.L1
  Sum.Sigma.L1.inv   <- solve(Sum.Sigma.L1) 

  Bs.mean.L1 <- as.numeric(Sigma.Bs.L1 %*% Sum.Sigma.L1.inv %*% (y.obs.trans-y.calc.L1.trans))
  Bs.var.L1  <- diag(Sigma.Bs.L1 %*% Sum.Sigma.L1.inv %*% (Sigma.Eps.L1)) #+Sigma.Bf.L1
  
  
  names(Bs.mean.L1) <- L1.encoded 
  names(Bs.var.L1)  <- L1.encoded

  
  
  Y.mean.L1 <- y.calc.L1.trans + (Sigma.Bs.L1)%*% Sum.Sigma.L1.inv %*% ( y.obs.trans - y.calc.L1.trans ) 
  Y.var.L1  <- diag((Sigma.Bs.L1) %*% Sum.Sigma.L1.inv %*% (Sigma.Eps.L1)) + diag(Sigma.Eps.L1) 
  
  names(Y.mean.L1) <- L1.encoded
  names(Y.var.L1)  <- L1.encoded
  

 
  return(list(y.calc.L1        = y.calc.L1.trans,
              Bs.mean.L1        = Bs.mean.L1,
              Bs.var.L1         = Bs.var.L1,
              Y.mean.L1         = Y.mean.L1,
              Y.var.L1          = Y.var.L1,
              y.calc.L2         = y.calc.L2,
  						fsigma_b2.nL1     = Sigma.Bs.L1[nrow(Sigma.Bs.L1),ncol(Sigma.Bs.L1)]

              ))
}



# ---------------------------------------------------------------------------
# Compute Hindcasting Distributions (vectorial version for L2)
# ---------------------------------------------------------------------------


sysanal.predict.bias.OU <- function(parsamp.L1,model,L1,y.obs, #CMP
                                 predict.bias.cond,
                                 ppt,
                                 probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),
                                 L2 =NA,y.calc=NA, 
																		par.tr,par.fix=NULL, inp,
                                 ...)
{
  # decode likelihood definitions:
  
  L1.decoded <- L1
  L1.encoded <- L1
  if ( !is.na(L1) )
  {
    L1.decoded <- sysanal.decode(L1)
  }
  else
  {
    L1.encoded <- rownames(L1)
  }
  n1 <- length(L1.encoded)
  
  L.decoded <- L1.decoded
  L.encoded <- L1.encoded
  n2 <- 0
  n <- n1

  L2.available <- TRUE
  if ( is.vector(L2) ) { if ( is.na(L2[1]) ) L2.available <- FALSE }   

    L2.decoded <- L2
    L2.encoded <- L2
  if ( !is.na(L2[1]) )
  {
    if ( is.vector(L2) )
    {
      L2.decoded <- sysanal.decode(L2)
    }
    else
    {
      L2.encoded <- rownames(L2)
    }
    n2 <- length(L2.encoded)
    
    if ( !is.na(L1[1]) ) L.decoded <- rbind(L1.decoded,L2.decoded)
    else  L.decoded <- L2.decoded
    if ( !is.na(L1[1]) ) L.encoded <- c(L1.encoded,L2.encoded)
    else  L.encoded <- L2.encoded
    if ( !is.na(L1[1]) ) n <- n1 + n2
    else  n <- n2
  }   else   {     n2 <- 0
        L.decoded <- L1.decoded
        L.encoded <- L1.encoded
   }
    n <- n1 + n2
    
  # initialize result arrays:
  
  y.L1.samp <- matrix(nrow=nrow(parsamp.L1),ncol=n1)
  colnames(y.L1.samp) <- L1.encoded
  Bs.L1.samp <- y.L1.samp
  Y.L1.samp <- y.L1.samp
  neg.var.Bs.L1 <- rep(0,n1)
  names(neg.var.Bs.L1) <- L1.encoded 
  neg.var.Y.L1 <- neg.var.Bs.L1

  
  Bf.L1.samp <- y.L1.samp
  
  y.L2.samp    <- NA
  Bs.L2.samp    <- NA
  Y.L2.samp    <- NA
  neg.var.Bs.L2 <- NA 
  neg.var.Y.L2 <- NA
  
  Bf.L2.samp    <- NA
  neg.var.Bf.L2 <- NA

  if ( !is.na(L2[1]) )
    {y.L2.samp <- matrix(nrow=nrow(parsamp.L1),ncol=n2)
    colnames(y.L2.samp) <- L2.encoded
     neg.var.Bs.L2 <- rep(0,n2)
     names(neg.var.Bs.L2) <- L2.encoded}
    Bs.L2.samp <- y.L2.samp
    Y.L2.samp <- y.L2.samp
    bias.slow.var <- y.L2.samp 
    neg.var.Y.L2 <- neg.var.Bs.L2
 
  # calculate samples:
  
  par.old <- rep(NA,ncol(parsamp.L1))
  num.eval <- 0
  for ( j in 1:nrow(parsamp.L1) )
  {
    par <- c(parsamp.L1[j,],par.fix)
    
    if ( !is.na(L1[1]) ) 
    {
    if ( j==1 | sum((par-par.old)^2) != 0 ) # do not reevaluate if parameter
    {                                       # values stayed the same
      par.old <- par
      num.eval <- num.eval + 1      
      res <- predict.bias.cond(par     = par,
                               model   = model,
                               L1      = L1, 
                               y.obs   = y.obs,
                               L2      = L2, 
                               y.calc  = y.calc,
      												 par.tr  = par.tr,
                               inp     = inp,
                               ...)
    }
    y.L1.samp[j,] <- res$y.calc.L1; # remember: it's transformed
    if ( !is.na(L2[1]) ) y.L2.samp[j,] <- res$y.calc.L2; # remember: it's transformed    
    
    } else { # i.e. if there is just L2 # *********
      
      y.calc.L2    <- model(par,L2,...)
      names(y.calc.L2)= rownames(L2)
      
      # transform results and observations FOR MULTIOUTPUT WITH DIFFERENT TRANSF OF EACH:
      out.vars  <- unique(as.character(L.decoded$var))
      
      y.calc.L2.trans <- rep(NA, length(y.calc.L2))
      
      for ( var in out.vars ) 
      {
        ind.var.L2 = which((sysanal.decode(L2)$var)==var)
        
        if(!is.na(par.tr["l1"])) #if(is.na(par.tr["alpha"])) 
        {
          y.calc.L2.trans [ind.var.L2]  <- sysanal.boxcox(y.calc.L2[ind.var.L2],par.tr["l1"],par.tr["l2"])
          
        } else {
          
          if (!is.na(par.tr["alpha"])) 
          {
            y.calc.L2.trans [ind.var.L2]  <- sysanal.logsinh(y.calc.L2[ind.var.L2],par.tr["alpha"],par.tr["beta"])
            
          } 
          else      {
            print("problem with transf param in pred")
          }
        }
        
      }
      y.L2.samp[j,]  <- y.calc.L2.trans 
      
    } # *********

    # Draw realization/OU paths of the errors
    # -----------------------------------
if ( !is.na(L1[1]) )
{
    for ( i in 1:n1 )
    {
      var <- res$Y.var.L1[i]
      if ( var < 0 ) 
      { 
        cat("* Warning: negative variance of Y:",var,"at",L1.encoded[i],"\n") 
        var <- 0; neg.var.Y.L1[i] <- neg.var.Y.L1[i]+1
      }
      Y.L1.samp[j,i] <- rnorm(1,res$Y.mean.L1[i],sqrt(var)) # Draw realizations for Y
    }
    if ( ! is.na(res$Bs.mean.L1[1]) )
    {
      for ( i in 1:n1 )
      {   
        var <- res$Bs.var.L1[i]
        if ( var < 0 ) 
        { 
          cat("* Warning: negative variance of Bs:",var,"at",L1.encoded[i],"\n") 
          var <- 0; neg.var.Bs.L1[i] <- neg.var.Bs.L1[i]+1 
        }
        Bs.L1.samp[j,i] <- rnorm(1,res$Bs.mean.L1[i],sqrt(var))  # Draw realizations for Bs 
      }
      for ( i in 1:n1 )
      {   
        var <- res$Bf.var.L1[i]
      }
    }
    } 
  }
  neg.var.Bs.L1 <- neg.var.Bs.L1/n1
  neg.var.Y.L1 <- neg.var.Y.L1/n1

  
  # Predictions in L2 
  # --------------------------------------------------

  if ( !is.na(L2[1]) )
  
  {Dt        <- L.decoded$val[2]-L.decoded$val[1]
  var       <- unique(L.decoded$var)
  
  for ( j in 1:nrow(parsamp.L1) ) # for all (selected) MCMC realizations (of several parameter sets) # Add: for every variable
  {                               # draw realizations of the observed system output and the model
    par       <- c(parsamp.L1[j,],par.fix)
    
    beta      <- 1/par["corrlen"]^1; if (is.na(par["corrlen"])) beta <- 0 # !!difference with original!!     if ( is.na(beta) ) beta <- 0
    exp.2dec  <- exp(-2*beta*(Dt))
    exp.dec   <- exp(-beta*(Dt))
    
    name.ks.B <- paste("ks",var,sep="_")
    ks.B      <- par[name.ks.B] ; if (is.na(ks.B)) ks.B <- 0
   
    name.sd.B <- paste("sd.B",var,sep="_")
    sigma_b2  <- par[name.sd.B]^2 ; if ( is.na(sigma_b2) ) sigma_b2  <- 0
    
    name.sd.Eps <- paste("sd.Eps",var,sep="_")
    sd.Eps      <- par[name.sd.Eps] ;  if ( is.na(sd.Eps) ) stop(paste("*** parameter", name.sd.Eps, "not found in sd.Eps.L"))
    
if ( !is.na(L1[1]) )  bs        = Bs.L1.samp[j,n1] 
else                  bs        = rnorm(1,0,par[name.sd.B])
    
    ppt_t = ppt[(par["Del.Max"]+1-round(par["Delta"])):(length(ppt)-round(par["Delta"]))]
    
    for ( i in 1:n2 ) # for all points in time L2
    { 
  	
    	slow_jump_var   = (sigma_b2+ (ppt_t[n1+i]*ks.B)^2)*(1-exp.2dec) # unconditional on previous variance
    
      Bs.L2.samp[j,i] = rnorm(1,bs*exp.dec,sqrt(slow_jump_var)) 

      E.L2.samp       = rnorm(1,0,sd.Eps)
    
      Y.L2.samp[j,i]  = y.L2.samp[j,i] + Bs.L2.samp[j,i] + E.L2.samp # realizations/points in the path
      
      bias.slow.var[j,i] = slow_jump_var
      
      bs             = Bs.L2.samp[j,i] # newB=act_val (dep on prev val=bs, jump_var=slow_jump_var)
      
      }
    }
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # derive quantiles:
  
  y.L1.quant <- matrix(nrow=length(probs),ncol=n1)
  colnames(y.L1.quant) <- L1.encoded
  rownames(y.L1.quant) <- probs
  Btot.L1.quant <- y.L1.quant
  yplusBtot.L1.quant <- y.L1.quant
  Y.L1.quant  <- y.L1.quant
 
if (n1 >0)
{
  
for ( i in 1:n1 )  # for all time points in inference period (we have distributions given by different paths)
  {
    y.L1.quant[,i]         <- quantile(y.L1.samp[,i],probs=probs,na.rm=TRUE)
    Btot.L1.quant[,i]      <- quantile(Bs.L1.samp[,i],probs=probs,na.rm=TRUE) 
    yplusBtot.L1.quant[,i] <- quantile(y.L1.samp[,i]+Bs.L1.samp[,i],probs=probs,na.rm=TRUE) # modif!!
    Y.L1.quant[,i]         <- quantile(Y.L1.samp[,i],probs=probs,na.rm=TRUE)
  }
}  
  y.L2.quant  <- NA
  Btot.L2.quant <- NA
  yplusBtot.L2.quant <- NA
  Y.L2.quant  <- NA
  if ( n2 > 0 )
  {
    y.L2.quant <- matrix(nrow=length(probs),ncol=n2)
    colnames(y.L2.quant) <- L2
    rownames(y.L2.quant) <- probs
    Btot.L2.quant <- y.L2.quant
    yplusBtot.L2.quant <- y.L2.quant
    Y.L2.quant <- y.L2.quant

    for ( i in 1:n2 )
    {
      y.L2.quant[,i]         <- quantile(y.L2.samp[,i],probs=probs,na.rm=TRUE)
      Btot.L2.quant[,i]      <- quantile(Bs.L2.samp[,i],probs=probs,na.rm=TRUE)
      yplusBtot.L2.quant[,i] <- quantile(y.L2.samp[,i]+Bs.L2.samp[,i],probs=probs,na.rm=TRUE)
      Y.L2.quant[,i]         <- quantile(Y.L2.samp[,i],probs=probs,na.rm=TRUE)
    }
  }
  
  return(list(y.L1.samp       = y.L1.samp,
              Y.L1.samp       = Y.L1.samp,
              Bs.L1.samp      = Bs.L1.samp,
              Bs.L2.samp      = Bs.L2.samp,
              Y.L2.samp       = Y.L2.samp,
              y.L2.samp       = y.L2.samp,
              yplusB.L1.quant = yplusBtot.L1.quant,
              yplusB.L2.quant = yplusBtot.L2.quant,
              Y.L1.quant      = Y.L1.quant,
              y.L1.quant      = y.L1.quant,
              Y.L2.quant      = Y.L2.quant,
              y.L2.quant      = y.L2.quant,
              Btot.L1.quant   = Btot.L1.quant,
              Btot.L2.quant   = Btot.L2.quant, 
              bias.slow.var.L2= bias.slow.var,
              neg.var.Bs.L1   = neg.var.Bs.L1,
              neg.var.Y.L1    = neg.var.Y.L1,
              neg.var.Bs.L2   = neg.var.Bs.L2,
              neg.var.Y.L2    = neg.var.Y.L2))
}


################################################################################

# sysanal.logsinh
# ==============

# purpose:
# output transformation


sysanal.logsinh <- function(data,alpha=0,beta=100)
{

	return(beta*log(sinh((alpha+data)*(1/beta))))
	
}

################################################################################

# sysanal.logsinh.deriv
# ====================

# purpose:
# calculate derivative of output transformation


sysanal.logsinh.deriv <- function(data,alpha=0,beta=100)
{
# 	return( 1/tanh(a+b*data)  )
	return( (tanh((alpha+data)*(1/beta)))^(-1) )
}

################################################################################

# sysanal.logsinh.inv
# ==================

# purpose:
# inverse output transformation


sysanal.logsinh.inv <- function(data,alpha=0,beta=100)
{
# 	return ( (asinh(exp(b*data))-a)/b )
	return ( (asinh(exp(data/beta))-alpha/beta)*beta)
	
}

# ---------------------------------------------------------------------------
# likelihood of the rainfall multipliers
# ---------------------------------------------------------------------------

sysanal.loglikeli.mult <- function(par,model,L,y.obs,Var.Bs, inp, sd.Eps,
                                       par.tr,par.fix=NULL, Inp=Inp,...)
{
  if (any (par<0))
  {return(-Inf)}
  
  # decode layout definition: 
  
  L.decoded <- L
  L.encoded <- L
  if ( is.vector(L) )
  {
    L.decoded <- sysanal.decode(L)
  }
  else
  {
    L.encoded <- rownames(L)
  }
  
  # calculate results of deterministic model:

  Inp=Inp*par["beta"]
  y.calc        <- model(par,L.encoded, Inp, ...) 

  # evaluate likelihood function:
  par.comb      <- c(par,par.fix)


  inp = inp[(par.comb["Del.Max"]+1-round(par.comb["Delta"])):(length(inp)-round(par.comb["Delta"]))]
  
  # transform results and observations:
  
  if(is.na(par.tr["alpha"])) 
  {
    y.calc.trans  <- sysanal.boxcox(y.calc,par.tr["l1"],par.tr["l2"])
    y.obs.trans   <- sysanal.boxcox(y.obs,par.tr["l1"],par.tr["l2"])
    boxcox.deriv  <- sysanal.boxcox.deriv(y.obs,par.tr["l1"],par.tr["l2"])
  }  else{
    y.calc.trans  <- sysanal.logsinh(y.calc,par.tr["alpha"],par.tr["beta"])
    y.obs.trans   <- sysanal.logsinh(y.obs,par.tr["alpha"],par.tr["beta"])
    boxcox.deriv  <- sysanal.logsinh.deriv(y.obs,par.tr["alpha"],par.tr["beta"])
  }
  

  Sigma.Eps     <- diag(sd.Eps(par.comb,L.decoded)^2)
  
  Sum.Sigma     <- Sigma.Eps
  
  Sum.Sigma.inv <- solve(Sum.Sigma)
  
  
  log.det.Sum.Sigma <- determinant(Sum.Sigma,logarithm=TRUE)
  
  if ( log.det.Sum.Sigma$sign < 0 ) 
  {warning("determinant Sigma.Eps+Sigma.B < 0") 
   loglikeli=-Inf}
  else{
    loglikeli <- - 0.5 * length(L.encoded) * log(2*pi) -
      0.5 * as.numeric(log.det.Sum.Sigma$modulus) -
      0.5 * t(y.obs.trans-y.calc.trans) %*% (Sum.Sigma.inv ) %*% 
      (y.obs.trans-y.calc.trans) +
      sum(log(abs(boxcox.deriv)))}

  
  return(loglikeli)
}

# ---------------------------------------------------------------------------
# Hindcasting Distributions
# ---------------------------------------------------------------------------


sysanal.predict.mult <- function(parsamp.L1,model,L1,y.obs, 
                                    predict.bias.cond,
                                    ppt,
                                    probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),
                                    L2 =NA,y.calc=NA, 
                                    par.tr,par.fix=NULL, inp, Inp,
                                    ...)
{
  # decode likelihood definitions:
  
  L1.decoded <- L1
  L1.encoded <- L1
  if ( !is.na(L1) )
  {
    L1.decoded <- sysanal.decode(L1)
  }
  else
  {
    L1.encoded <- rownames(L1)
  }
  n1 <- length(L1.encoded)
  
  L.decoded <- L1.decoded
  L.encoded <- L1.encoded
  n2 <- 0
  n <- n1
  
  L2.available <- TRUE
  if ( is.vector(L2) ) { if ( is.na(L2[1]) ) L2.available <- FALSE }   
  
  L2.decoded <- L2
  L2.encoded <- L2
  if ( !is.na(L2[1]) )
  {
    if ( is.vector(L2) )
    {
      L2.decoded <- sysanal.decode(L2)
    }
    else
    {
      L2.encoded <- rownames(L2)
    }
    n2 <- length(L2.encoded)
    
    if ( !is.na(L1[1]) ) L.decoded <- rbind(L1.decoded,L2.decoded)
    else  L.decoded <- L2.decoded
    if ( !is.na(L1[1]) ) L.encoded <- c(L1.encoded,L2.encoded)
    else  L.encoded <- L2.encoded
    if ( !is.na(L1[1]) ) n <- n1 + n2
    else  n <- n2
  }   else   {     n2 <- 0
                   L.decoded <- L1.decoded
                   L.encoded <- L1.encoded
  }
  n <- n1 + n2
  
  # initialize result arrays:
  
  y.L1.samp <- matrix(nrow=nrow(parsamp.L1),ncol=n1)
  colnames(y.L1.samp) <- L1.encoded
  Bs.L1.samp <- y.L1.samp
  Y.L1.samp <- y.L1.samp
  neg.var.Bs.L1 <- rep(0,n1)
  names(neg.var.Bs.L1) <- L1.encoded 
  neg.var.Y.L1 <- neg.var.Bs.L1
  #   neg.var.Bf.L1 <- rep(0,n1)
  #   names(neg.var.Bf.L1) <- L1.encoded
  
  Bf.L1.samp <- y.L1.samp
  
  y.L2.samp    <- NA
  Bs.L2.samp    <- NA
  Y.L2.samp    <- NA
  neg.var.Bs.L2 <- NA 
  neg.var.Y.L2 <- NA
  
  Bf.L2.samp    <- NA
  neg.var.Bf.L2 <- NA
  #   
  if ( !is.na(L2[1]) )
  {y.L2.samp <- matrix(nrow=nrow(parsamp.L1),ncol=n2)
   colnames(y.L2.samp) <- L2.encoded
   neg.var.Bs.L2 <- rep(0,n2)
   names(neg.var.Bs.L2) <- L2.encoded}
  Bs.L2.samp <- y.L2.samp
  Y.L2.samp <- y.L2.samp
  bias.slow.var <- y.L2.samp 
  neg.var.Y.L2 <- neg.var.Bs.L2

  
  # calculate samples:
  
  par.old <- rep(NA,ncol(parsamp.L1))
  num.eval <- 0
  for ( j in 1:nrow(parsamp.L1) )
  {
    par <- c(parsamp.L1[j,],par.fix)
    
    if ( !is.na(L1[1]) ) 
    {
      if ( j==1 | sum((par-par.old)^2) != 0 ) # do not reevaluate if parameter
      {                                       # values stayed the same
        par.old <- par
        num.eval <- num.eval + 1      
        res <- predict.bias.cond(par     = par,
                                 model   = model,
                                 L1      = L1, 
                                 y.obs   = y.obs,
                                 L2      = L2, 
                                 y.calc  = y.calc,
                                 par.tr  = par.tr,
                                 inp     = inp,
                                 Inp     = Inp,
                                 ...)
      }
      y.L1.samp[j,] <- res$y.calc.L1; # remember: it's transformed
      if ( !is.na(L2[1]) ) y.L2.samp[j,] <- res$y.calc.L2; # remember: it's transformed    
      
    } else { # i.e. if there is just L2 # *********
      Unc.Inp=Inp*par["beta"] # implemented for both L1 and L2 but separated
      y.calc.L2    <- model(par,L2,Unc.Inp, ...)
      names(y.calc.L2)= rownames(L2)
      
      # transform results and observations FOR MULTIOUTPUT WITH DIFFERENT TRANSF OF EACH:
      out.vars  <- unique(as.character(L.decoded$var))
      
      y.calc.L2.trans <- rep(NA, length(y.calc.L2))
      
      for ( var in out.vars ) 
      {
        ind.var.L2 = which((sysanal.decode(L2)$var)==var)
        
        if(!is.na(par.tr["l1"])) 
        {
          y.calc.L2.trans [ind.var.L2]  <- sysanal.boxcox(y.calc.L2[ind.var.L2],par.tr["l1"],par.tr["l2"])
          
        } else {
          
          if (!is.na(par.tr["alpha"])) 
          {
            y.calc.L2.trans [ind.var.L2]  <- sysanal.logsinh(y.calc.L2[ind.var.L2],par.tr["alpha"],par.tr["beta"])
            
          } 
          else      {
            print("problem with transf param in pred")
          }
        }
        
      }
      y.L2.samp[j,]  <- y.calc.L2.trans 
      
    } # *********
    
    # Draw realization/OU paths of the errors
    # -----------------------------------
    if ( !is.na(L1[1]) )
    {
      for ( i in 1:n1 )
      {
        var <- res$Y.var.L1[i]
        if ( var < 0 ) 
        { 
          cat("* Warning: negative variance of Y:",var,"at",L1.encoded[i],"\n") 
          var <- 0; neg.var.Y.L1[i] <- neg.var.Y.L1[i]+1
        }
        Y.L1.samp[j,i] <- rnorm(1,res$Y.mean.L1[i],sqrt(var)) # Draw realizations for Y
      }
      if ( ! is.na(res$Bs.mean.L1[1]) )
      {
        for ( i in 1:n1 )
        {   
          var <- res$Bs.var.L1[i]
          if ( var < 0 ) 
          { 
            cat("* Warning: negative variance of Bs:",var,"at",L1.encoded[i],"\n") 
            var <- 0; neg.var.Bs.L1[i] <- neg.var.Bs.L1[i]+1 
          }
          Bs.L1.samp[j,i] <- rnorm(1,res$Bs.mean.L1[i],sqrt(var))  # Draw realizations for Bs 
        }
        for ( i in 1:n1 )
        {   
          var <- res$Bf.var.L1[i]
        }
      }
    } 
  }
  neg.var.Bs.L1 <- neg.var.Bs.L1/n1
  neg.var.Y.L1 <- neg.var.Y.L1/n1
  
  
  # predictions in L2 
  # --------------------------------------------------
  
  if ( !is.na(L2[1]) )
    
  {Dt        <- L.decoded$val[2]-L.decoded$val[1]
   var       <- unique(L.decoded$var)
   
   for ( j in 1:nrow(parsamp.L1) ) # for all (selected) MCMC realizations (of several parameter sets) # Add: for every variable
   {                               # draw realizations of the observed system output and the model
     par       <- c(parsamp.L1[j,],par.fix)
     
     beta      <- 1/par["corrlen"]^1; if (is.na(par["corrlen"])) beta <- 0 # !!difference with original!!     if ( is.na(beta) ) beta <- 0
     exp.2dec  <- exp(-2*beta*(Dt))
     exp.dec   <- exp(-beta*(Dt))
     
     name.ks.B <- paste("ks",var,sep="_")
     ks.B      <- par[name.ks.B] ; if (is.na(ks.B)) ks.B <- 0
    
     name.sd.B <- paste("sd.B",var,sep="_")
     sigma_b2  <- par[name.sd.B]^2 ; if ( is.na(sigma_b2) ) sigma_b2  <- 0
     
     
     name.sd.Eps <- paste("sd.Eps",var,sep="_")
     sd.Eps      <- par[name.sd.Eps] ;  if ( is.na(sd.Eps) ) stop(paste("*** parameter", name.sd.Eps, "not found in sd.Eps.L"))
     
     if ( !is.na(L1[1]) )  bs        = Bs.L1.samp[j,n1] 
     else                  bs        = rnorm(1,0,par[name.sd.B])
     
     ppt_t = ppt[(par["Del.Max"]+1-round(par["Delta"])):(length(ppt)-round(par["Delta"]))]
     
     for ( i in 1:n2 ) # for all points in time L2
     { 
   
       slow_jump_var   = (sigma_b2+ (ppt_t[n1+i]*ks.B)^2)*(1-exp.2dec) # unconditional on previous variance
       
       Bs.L2.samp[j,i] = rnorm(1,bs*exp.dec,sqrt(slow_jump_var)) 
       
       E.L2.samp       = rnorm(1,0,sd.Eps)
       
       Y.L2.samp[j,i]  = y.L2.samp[j,i] + Bs.L2.samp[j,i] + E.L2.samp # realizations/points in the path
       
       
       bias.slow.var[j,i] = slow_jump_var
       
       bs             = Bs.L2.samp[j,i] # newB=act_val (dep on prev val=bs, jump_var=slow_jump_var)
       
     }
   }
  }
  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # derive quantiles:
  
  y.L1.quant <- matrix(nrow=length(probs),ncol=n1)
  colnames(y.L1.quant) <- L1.encoded
  rownames(y.L1.quant) <- probs
  Btot.L1.quant <- y.L1.quant
  yplusBtot.L1.quant <- y.L1.quant
  Y.L1.quant  <- y.L1.quant
  
  if (n1 >0)
  {
    
    for ( i in 1:n1 )  # for all time points in inference period (we have distributions given by different paths)
    {
      y.L1.quant[,i]         <- quantile(y.L1.samp[,i],probs=probs,na.rm=TRUE)
      Btot.L1.quant[,i]      <- quantile(Bs.L1.samp[,i],probs=probs,na.rm=TRUE) 
      yplusBtot.L1.quant[,i] <- quantile(y.L1.samp[,i]+Bs.L1.samp[,i],probs=probs,na.rm=TRUE) # modif!!
      Y.L1.quant[,i]         <- quantile(Y.L1.samp[,i],probs=probs,na.rm=TRUE)
    }
  }  
  y.L2.quant  <- NA
  Btot.L2.quant <- NA
  yplusBtot.L2.quant <- NA
  Y.L2.quant  <- NA
  if ( n2 > 0 )
  {
    y.L2.quant <- matrix(nrow=length(probs),ncol=n2)
    colnames(y.L2.quant) <- L2
    rownames(y.L2.quant) <- probs
    Btot.L2.quant <- y.L2.quant
    yplusBtot.L2.quant <- y.L2.quant
    Y.L2.quant <- y.L2.quant
    
    for ( i in 1:n2 )
    {
      y.L2.quant[,i]         <- quantile(y.L2.samp[,i],probs=probs,na.rm=TRUE)
      Btot.L2.quant[,i]      <- quantile(Bs.L2.samp[,i],probs=probs,na.rm=TRUE)
      yplusBtot.L2.quant[,i] <- quantile(y.L2.samp[,i]+Bs.L2.samp[,i],probs=probs,na.rm=TRUE)
      Y.L2.quant[,i]         <- quantile(Y.L2.samp[,i],probs=probs,na.rm=TRUE)
    }
  }
  
  return(list(y.L1.samp       = y.L1.samp,
              Y.L1.samp       = Y.L1.samp,
              Bs.L1.samp      = Bs.L1.samp,
              Bs.L2.samp      = Bs.L2.samp,
              Y.L2.samp       = Y.L2.samp,
              y.L2.samp       = y.L2.samp,
              yplusB.L1.quant = yplusBtot.L1.quant,
              yplusB.L2.quant = yplusBtot.L2.quant,
              Y.L1.quant      = Y.L1.quant,
              y.L1.quant      = y.L1.quant,
              Y.L2.quant      = Y.L2.quant,
              y.L2.quant      = y.L2.quant,
              Btot.L1.quant   = Btot.L1.quant,
              Btot.L2.quant   = Btot.L2.quant, 
              bias.slow.var.L2= bias.slow.var,
              neg.var.Bs.L1   = neg.var.Bs.L1,
              neg.var.Y.L1    = neg.var.Y.L1,
              neg.var.Bs.L2   = neg.var.Bs.L2,
              neg.var.Y.L2    = neg.var.Y.L2))
}


################################################################################

# Function implementing a simple standard posterior probability density (up to
# a multiplicative normalizing constant). It is based on a likelihood 
# definition provided by the user and a multivariate normal or
# lognormal distribution at the prior.
# We here have three (independent) time series that we want to simultaneously calibrate

sysanal.logposterior.3 <- function(par,model,L.1,y.1,L.2,y.2,L.3,y.3,
                                 prior.dist="lognormal",prior.mean=1,prior.sd=1,
                                 prior.cor=NA,prior.def=NA,
                                 loglikeli=sysanal.loglikeli,
                                 inp.1,Inp.1,
                                 inp.2,Inp.2,
                                 inp.3,Inp.3, 
                                 ...)
{

  logprior  <- sysanal.calcpdf_mv(z=par,dist=prior.dist,mean=prior.mean,
                                  sd=prior.sd,cor=prior.cor,distdef=prior.def)

  if ( is.na(logprior) ) return(NA)
  par["beta"]=par["beta1"]
  loglikeli1 <- loglikeli(par=par,model=model,L=L.1,y=y.1,inp=inp.1,Inp=Inp.1,...)
  par["beta"]=par["beta2"]
  loglikeli2 <- loglikeli(par=par,model=model,L=L.2,y=y.2, inp=inp.2, Inp=Inp.2,...)
  par["beta"]=par["beta3"]
  loglikeli3 <- loglikeli(par=par,model=model,L=L.3,y=y.3, inp=inp.3, Inp=Inp.3,...)
  
  loglikeli <- loglikeli1+loglikeli2+loglikeli3

  return(logprior+loglikeli)
}


################################################################################
# ---------------------------------------------------------------------------
# For conditional predictions with multipliers and iid out err (i.e. calculate smoothing distributions)
# ---------------------------------------------------------------------------


sysanal.predict.mult.L1 <- function(par,model,L1,y.obs,
                                        Var.Bs, inp,sd.Eps,
                                        L2,
                                        y.calc=NA, par.tr,Inp,  ...)
{
  # decode likelihood definitions:
  
  L1.decoded <- L1
  L1.encoded <- L1
  if ( is.vector(L1) )
  {
    L1.decoded <- sysanal.decode(L1)
  }
  else
  {
    L1.encoded <- rownames(L1)
  }
  n1 <- length(L1.encoded)
  
  L.decoded <- L1.decoded
  L.encoded <- L1.encoded
  n2 <- 0
  n <- n1
  
  # calculate results of deterministic model (we assume for now only L1):
  Unc.Inp=Inp*par["beta"]
  
  if ( is.na(y.calc[1]) )
  {   
    if ( !is.na(L2[1]) ) y.calc    <- model(par,c(L1,L2),Unc.Inp, ...) else y.calc    <- model(par,c(L1),Unc.Inp, ...)
  }
  else
  {
    if ( length(y.calc) != n )
    {
      cat("*** y.calc is not of correct length:",length(y.calc),
          "instead of",n,"\n")
      return(NA)
    }
  }
  
  y.calc.L1 <- y.calc[1:n1] 
  y.calc.L2 <- NA
  
  # transform results and observations:
  
  if(is.na(par.tr["alpha"])) 
  {
    y.calc.L1.trans  <- sysanal.boxcox(y.calc.L1,par.tr["l1"],par.tr["l2"])
    if ( !is.na(L2[1]) )	y.calc.L2        <- sysanal.boxcox(y.calc[(n1+1):length(y.calc)],par.tr["l1"],par.tr["l2"])
    y.obs.trans      <- sysanal.boxcox(y.obs,par.tr["l1"],par.tr["l2"])
  }
  else{
    y.calc.L1.trans  <- sysanal.logsinh(y.calc.L1,par.tr["alpha"],par.tr["beta"])
    if ( !is.na(L2[1]) )	y.calc.L2        <- sysanal.logsinh(y.calc[(n1+1):length(y.calc)],par.tr["alpha"],par.tr["beta"])
    y.obs.trans      <- sysanal.logsinh(y.obs,par.tr["alpha"],par.tr["beta"])
  }
  
  # calculate predictions for layout 1:   
  
  inp = inp[(par["Del.Max"]+1-round(par["Delta"])):(length(inp)-round(par["Delta"]))]
  
  Sigma.Bs.L1       <- Var.Bs(par,L1.decoded,inp) #slow B cov
  Sigma.Eps.L1      <- diag(sd.Eps(par,L1.decoded)^2)
  
  Sum.Sigma.L1       <- Sigma.Bs.L1 + Sigma.Eps.L1
  Sum.Sigma.L1.inv   <- solve(Sum.Sigma.L1) 

  Bs.mean.L1 <- as.numeric(Sigma.Bs.L1 %*% Sum.Sigma.L1.inv %*% (y.obs.trans-y.calc.L1.trans))
  Bs.var.L1  <- diag(Sigma.Bs.L1 %*% Sum.Sigma.L1.inv %*% (Sigma.Eps.L1)) #+Sigma.Bf.L1
  
  names(Bs.mean.L1) <- L1.encoded 
  names(Bs.var.L1)  <- L1.encoded
  
  Y.mean.L1 <- y.calc.L1.trans + (Sigma.Bs.L1)%*% Sum.Sigma.L1.inv %*% ( y.obs.trans - y.calc.L1.trans ) 
  Y.var.L1  <- diag((Sigma.Bs.L1) %*% Sum.Sigma.L1.inv %*% (Sigma.Eps.L1)) + diag(Sigma.Eps.L1) 
  
  names(Y.mean.L1) <- L1.encoded
  names(Y.var.L1)  <- L1.encoded
  
  
  
  return(list(y.calc.L1        = y.calc.L1.trans,
              Bs.mean.L1        = Bs.mean.L1,
              Bs.var.L1         = Bs.var.L1,
              Y.mean.L1         = Y.mean.L1,
              Y.var.L1          = Y.var.L1,
              y.calc.L2         = y.calc.L2,
              fsigma_b2.nL1     = Sigma.Bs.L1[nrow(Sigma.Bs.L1),ncol(Sigma.Bs.L1)]

  ))
}


################################################################################
# function transforming a normal autocorrelated process into the rainfall process
################################################################################

tr.norm2ppt <- function(x,par)
{
  y  <- ifelse(x<par["b1"], 0, par["a1"]*(x-par["b1"])^par["alpha1"] )
  b2 <- par["x0"] - ( par["alpha1"]*par["a1"]/(par["alpha2"]*par["a2"])*(par["x0"]-par["b1"])^(par["alpha1"]-1) )^(1/(par["alpha2"]-1))
  c2 <- par["a1"]*(par["x0"]-par["b1"])^par["alpha1"] - par["a2"]*(par["x0"]-b2)^par["alpha2"]
  y  <- ifelse(x<par["x0"],y,par["a2"]*(x-b2)^par["alpha2"]+c2)
  return(y)
}




################################################################################
# Draw a realization of the input observation from a conditioned mvtnd
################################################################################

norm.inp.obs.gen <- function(ti.gr.Ime,ti.gr.1mi, inp.meas, ppt2norm.tra, par, inp.elsew)
{  
  # -----------------------------------------------------------------------
  # This function generate a sample path of the transformed rainfall at the pluviometric station (Eq. 17)
  # For time points where the rain was > 0 it simply gives back the transformed value
  # For time points where the rain was = 0 it draws from a correlated truncated normal distribution
  #
  # Arguments:
  # ti.gr        time grid of the period (e.g. defining at which hour every meas was taken)
  # inp.meas     measured rain over the period
  # ppt2norm.tra function transforming (positive) rain into a (truncated) normal distribution. Already parametrized
  # par          parameters defining: i) the upper truncation of the normal distribution (its a quantile of the norm distr) [thresh],
  #              ii) variance of the difference of two transf rain time series with a distance equal to this one [var.xi], iii) correlation length (e.g. in hr)  {not used}
  # inp.elsew    generated stochastic process representing transformed precipitation over the catchment
  
  #
  # Return a list with:
  # xio          sample/trace of the transf rain at the measurement site.
  #
  #                                first version:         March 2015 - Dario Del Giudice
  # ---------------------------------------------------------------------------------------
  
  xio = inp.meas*NA # stoch input process at the observation site
  
  for (i in 1:length(xio))
  {
    if(inp.meas[i]>0)
    {xio[i] = ppt2norm.tra(inp.meas[i])}
  } # we created xio.I
  
  len.I0 <- length(which(is.na(xio))) # length I0 (time interval where 0 rain was observed)
  
  # Build covariance matrix defining how the process at 2 distinct points in time are correlated
  # -----------------
  
  # lower and upper truncation of the the Truncated Multivariate Normal Distribution
  lo.tr = rep(-Inf, len.I0)
  up.tr = rep(par["thresh"], len.I0)
  
  var.xio  = par["var.xi"] # transf rain units square
  
  # make sure that xi has time information attached (it is solved at the minimal resolution)
  names(inp.elsew) = ti.gr.1mi
  
  xi.I0  <- inp.elsew[which(is.na(xio))] # select times with no measurem
  time.I0  <- as.numeric(ti.gr.Ime[which(is.na(xio))])
  time.I   <- as.numeric(ti.gr.Ime[which(!is.na(xio))]) # time points (in hr) where rain was >0

  Sigma.xio.I0.I = diag(var.xio, nrow=length(time.I0), ncol=length(time.I))

  Sigma.xio.I.I = diag(var.xio, nrow=length(time.I), ncol=length(time.I))
  
  inv.Sigma.xio.I.I = diag(1/var.xio, nrow=length(time.I), ncol=length(time.I))
  
  ind.I.xi = match(as.numeric(ti.gr.Ime[which(!is.na(xio))]),as.numeric(ti.gr.1mi)); if(any(is.na(ind.I.xi))) stop("no match between time grids")
  
  xio.min.xi.I = xio[which(!is.na(xio))]-inp.elsew[ind.I.xi] # compare the 2 transf rain during rain (at appropr resol)

  Sigma.xio.I0.I0 = diag(var.xio, nrow=length(time.I0), ncol=length(time.I0)) 

  mu.xio.I0 = xi.I0 + as.numeric(Sigma.xio.I0.I %*% inv.Sigma.xio.I.I %*% xio.min.xi.I)
  
#   sigma.xio.I0 = Sigma.xio.I0.I0 - Sigma.xio.I0.I %*% inv.Sigma.xio.I.I %*% t(Sigma.xio.I0.I)
#   H.xio.I0 = solve(sigma.xio.I0) # precision matrix, the inverse of the covariance matrix
  H.xio.I0 = diag(1/var.xio, nrow=length(time.I0), ncol=length(time.I0))
 
  xio.I0 <- as.numeric(rtmvnorm(n=1, mean=mu.xio.I0, H=H.xio.I0, lower=lo.tr, upper=up.tr, algorithm="gibbs")) # Gibbs sampling with precision matrix
  
  names(xio.I0) <- time.I0
  
  xio.NA = xio
  
  xio[which(is.na(xio))] = xio.I0
  
  return(xio)
  
}


################################################################################
# Sample posterior out.par, inp.par, inp.proc (i.e. run inference with SIP)
################################################################################

mcmc.sip.3e <- function(out.par.ini ,inp.par.fix ,  inp.par.ini ,t.grid.1mi.1 , t.grid.Ime.1 ,
                        t.grid.1mi.2 ,t.grid.Ime.2 ,t.grid.1mi.3 ,t.grid.Ime.3 ,inp.cand.cov ,inp.logpri ,
                        L.1 , y.1 , L.2 , y.2 , L.3 , y.3 , interv.pts ,ppt2norm.tra ,out.cand.cov ,out.loglikh ,
                        out.logpri ,  meas.rain.1 , meas.rain.2 ,meas.rain.3,sampsize,adapt_par=c(20000000000,200,0.5,0.55),
                        traj.Inp.pr.1=NA,traj.Inp.pr.2=NA,traj.Inp.pr.3=NA, thin=1)
{  
  # -----------------------------------------------------------------------
  # This function samples the posterior distribution of the simulator parameters,
  # error model parameters, and stochastic input process.
  # It uses the Random-Walk Metropolis-within-Gibbs, alias the Metropolis algorithm,
  # to sequentially sample the joint log.posterior in three steps.
  #
  # Arguments:
  # out.loglikh  function defining the output error model of the deterministic
  #              model parameters as well as the parameters of the output error model; it uses generated rain
  # out.logpri   function defining the the prior density of the deterministic
  #              model parameters as well as the parameters of the output error model; it uses generated rain
  # inp.logpri   function defining the the prior density of the parameters of the (true) input OU process
  # inp.pro.den  function definining the density of the input stochastic process
  # inp.dens     function definining the mass of the input observations
  # inp.obs.par  list defining the prior distribution of the parameters of the input observation error model
  # 
  # out.par.ini  vector of values of output parameters at which the chain is to be started.
  # inp.cand.cov covariance matrix of the candidate density (or jump distr.) for input OU process
  # out.cand.cov covariance matrix of the candidate density (or jump distr.) for the output lkh function
  # sampsize     MCMC iterations (length of the chain)
  # adapt.Metr   algorithm to sample the parameters
  # t.grid.cal   vector of temporal points where input data are available
  # inp.trans    function to transform X into rainfall
  # meas.rain    observed input with its true data resolution 
  # traj.Inp.pr  we can provide a good starting point for the stochastic input process
  # ...          additional parameters passed to "mcmc.sip.3e" or its functions
  #
  # Return a list with:
  # post.out.par sample/trace of the target distribution of out parameters as a matrix with sample points in its rows.
  # post.inp.par sample/trace of the target distribution of inp parameters as a matrix with sample points in its rows.
  # out.log.dens matrix with the log-densities of the functions evaluated at the accepted parameters/processes
  # post.inp     sample of the realizations of the inferred input process (each matrix row is a path)
  # inp.log.dens density associated with the input sample
  # hyd.mod.run  number of executions of the deterministic model
  # out.cand.cov jump covariance matrix for output sampling
  # xio          realization of the input process at the observation site
  #
  #                                first version:         January 2015 - Dario Del Giudice
  # ---------------------------------------------------------------------------------------
  
  # initialization
  library(MASS); library(mvtnorm)   # ... for Multivariate Normal Distribution
  library("MHadaptive") # ... for making the cov matrices a symmetric with all positive eigenvalues
  out.nacc  <- 0 ; inp.pro.nacc.1 <- 0; inp.pro.nacc.2 <- 0; inp.pro.nacc.3 <- 0; inp.par.nacc <- 0 
  out.npro  <- 0 ; inp.pro.npro <- 0; inp.npro  <- 0
  reducsize <- floor(sampsize/thin)+1

  post.out.par <- array(dim=c(reducsize,length(out.par.ini))) ;  colnames(post.out.par)  <- names(out.par.ini)
  post.out.par.full <- array(dim=c(sampsize,length(out.par.ini))) ;  colnames(post.out.par.full)  <- names(out.par.ini)
  post.inp.par <- array(dim=c(reducsize,length(inp.par.ini))) ;  colnames(post.inp.par)  <- names(inp.par.ini)
  post.inp.1  <- array(dim=c(reducsize,length(t.grid.1mi.1)));  colnames(post.inp.1) <- t.grid.1mi.1
  post.inp.2  <- array(dim=c(reducsize,length(t.grid.1mi.2)));  colnames(post.inp.2) <- t.grid.1mi.2
  post.inp.3  <- array(dim=c(reducsize,length(t.grid.1mi.3)));  colnames(post.inp.3) <- t.grid.1mi.3
  
  print(paste("Memory need for input samples: ", format(object.size(cbind(post.inp.1,post.inp.2, post.inp.3)),units="Mb")))
  log.dens  <- array(dim=c( reducsize,2))
  inp.log.dens.1  <- array(dim=c( reducsize,1))
  inp.log.dens.2  <- array(dim=c( reducsize,1))
  inp.log.dens.3  <- array(dim=c( reducsize,1))
 
  colnames(inp.log.dens.1) = c("DensInp.1")
  colnames(inp.log.dens.2) = c("DensInp.2")
  colnames(inp.log.dens.3) = c("DensInp.3")
  colnames(log.dens) = c("PriorOutPar","LkhOutPar") 
  out.par   <- out.par.ini
  inp.par   <- inp.par.ini
  inp.par.cand <- inp.par.ini
  out.cand.cov<-makePositiveDefinite(out.cand.cov)
  acc.opt     <- .3

 # define starting paths
if(is.na(traj.Inp.pr.1))  traj.Inp.pr.1 <-rou(tau=inp.par.fix["corrl_inp"],y1=NA,yn=NA,t=as.numeric(t.grid.1mi.1))$y# draw a first realization of the OU process unconditionned 
names(traj.Inp.pr.1) = t.grid.1mi.1
if(is.na(traj.Inp.pr.2))    traj.Inp.pr.2 <-rou(tau=inp.par.fix["corrl_inp"],y1=NA,yn=NA,t=as.numeric(t.grid.1mi.2))$y# draw a first realization of the OU process unconditionned 
  names(traj.Inp.pr.2) = t.grid.1mi.2
if(is.na(traj.Inp.pr.3))    traj.Inp.pr.3 <-rou(tau=inp.par.fix["corrl_inp"],y1=NA,yn=NA,t=as.numeric(t.grid.1mi.3))$y# draw a first realization of the OU process unconditionned 
  names(traj.Inp.pr.3) = t.grid.1mi.3
 

  t.grid.1mi.1  <-as.numeric(t.grid.1mi.1)
  t.grid.Ime.1  <-as.numeric(t.grid.Ime.1)
  t.grid.1mi.2  <-as.numeric(t.grid.1mi.2)
  t.grid.Ime.2  <-as.numeric(t.grid.Ime.2)
  t.grid.1mi.3  <-as.numeric(t.grid.1mi.3)
  t.grid.Ime.3  <-as.numeric(t.grid.Ime.3)
 
  out.lkh.call<- 0
   
  xio.1 <- norm.inp.obs.gen(ti.gr.Ime=t.grid.Ime.1,ti.gr.1mi=t.grid.1mi.1, inp.meas=meas.rain.1, 
                               ppt2norm.tra = ppt2norm.tra, 
                               par = c(inp.par,inp.par.fix), 
                               inp.elsew = traj.Inp.pr.1)
  
  xio.2 <- norm.inp.obs.gen(ti.gr.Ime=t.grid.Ime.2,ti.gr.1mi=t.grid.1mi.2, inp.meas=meas.rain.2, 
                            ppt2norm.tra = ppt2norm.tra, 
                            par = c(inp.par,inp.par.fix), 
                            inp.elsew = traj.Inp.pr.2)
  
  xio.3 <- norm.inp.obs.gen(ti.gr.Ime=t.grid.Ime.3,ti.gr.1mi=t.grid.1mi.3, inp.meas=meas.rain.3, 
                            ppt2norm.tra = ppt2norm.tra, 
                            par = c(inp.par,inp.par.fix), 
                            inp.elsew = traj.Inp.pr.3)
  

  # select the appropriate time points to compare xio and xi
  ind.comp.xi.1 = match(t.grid.Ime.1,t.grid.1mi.1); if(any(is.na(ind.comp.xi.1))) stop("no match between time grids")
  ind.comp.xi.2 = match(t.grid.Ime.2,t.grid.1mi.2); if(any(is.na(ind.comp.xi.2))) stop("no match between time grids")
  ind.comp.xi.3 = match(t.grid.Ime.3,t.grid.1mi.3); if(any(is.na(ind.comp.xi.3))) stop("no match between time grids")
  
  
  Sigma.xio.1 <- diag(inp.par["var.xi"], nrow=length(t.grid.Ime.1), ncol=length(t.grid.Ime.1))
  Sigma.xio.2 <- diag(inp.par["var.xi"], nrow=length(t.grid.Ime.2), ncol=length(t.grid.Ime.2))
  Sigma.xio.3 <- diag(inp.par["var.xi"], nrow=length(t.grid.Ime.3), ncol=length(t.grid.Ime.3))

  inp.dens.1  <-  dmvnorm(xio.1, mean = traj.Inp.pr.1[ind.comp.xi.1], sigma = Sigma.xio.1, log = T)
  inp.dens.2  <-  dmvnorm(xio.2, mean = traj.Inp.pr.2[ind.comp.xi.2], sigma = Sigma.xio.2, log = T)
  inp.dens.3  <-  dmvnorm(xio.3, mean = traj.Inp.pr.3[ind.comp.xi.3], sigma = Sigma.xio.3, log = T)
  
  
  tr.inp.err.1 <- sum(abs(xio.1-traj.Inp.pr.1[ind.comp.xi.1]))
  tr.inp.err.2 <- sum(abs(xio.2-traj.Inp.pr.2[ind.comp.xi.2]))
  tr.inp.err.3 <- sum(abs(xio.3-traj.Inp.pr.3[ind.comp.xi.3]))
  
z=1

  for(i in 1:sampsize)
  {
    z=z+1
    out.npro     <- out.npro+1
    inp.npro     <- inp.npro+1
    
    # Step I: theta, sample posterior ouput parameters
    # -----------------------------------------

# I.0: (re-)calculate posterior with old parameter values but new input time series: 
out.lkh.1  <- out.loglikh(out.par ,L=L.1, y.obs=y.1,sim.inp=inp.trans(traj.Inp.pr.1)*60/1000,t.grid.cal=t.grid.1mi.1)
out.lkh.2  <- out.loglikh(out.par ,L=L.2, y.obs=y.2,sim.inp=inp.trans(traj.Inp.pr.2)*60/1000,t.grid.cal=t.grid.1mi.2)
out.lkh.3  <- out.loglikh(out.par ,L=L.3, y.obs=y.3,sim.inp=inp.trans(traj.Inp.pr.3)*60/1000,t.grid.cal=t.grid.1mi.3)
out.lkh    <- out.lkh.1+out.lkh.2+out.lkh.3
out.post   <- out.lkh + out.logpri(out.par)

    
    # I.i: Generate a candidate for the next sample output posterior parameter 
    out.par.cand <- mvrnorm(1,mu=out.par,Sigma=out.cand.cov)  # Draw candidate point from the proposal/jump density
    
    # I.ii: Decide whether to accept or reject the candidate output posterior parameter 
    out.lkh.1.cand  <- out.loglikh(out.par.cand ,L=L.1, y.obs=y.1,sim.inp=inp.trans(traj.Inp.pr.1)*60/1000,t.grid.cal=t.grid.1mi.1)
    out.lkh.2.cand  <- out.loglikh(out.par.cand ,L=L.2, y.obs=y.2,sim.inp=inp.trans(traj.Inp.pr.2)*60/1000,t.grid.cal=t.grid.1mi.2)
    out.lkh.3.cand  <- out.loglikh(out.par.cand ,L=L.3, y.obs=y.3,sim.inp=inp.trans(traj.Inp.pr.3)*60/1000,t.grid.cal=t.grid.1mi.3)
    
    out.lkh.cand <- out.lkh.1.cand+out.lkh.2.cand+out.lkh.3.cand
    
    out.post.cand <- out.lkh.cand + out.logpri(out.par.cand)    # evaluate the unnormalized posterior at proposed parameter set
    
    out.ratio    <- out.post.cand-out.post  # log of ratio of densities becomes a difference 
    
    if(!is.finite(out.ratio))   out.ratio<- -Inf      # never jump outside the range of the parameter space
    
    if( log(runif(1,0,1)) <= out.ratio)          # Make the jump according to acceptance probability
    {
      out.par   = out.par.cand
      out.post  = out.post.cand
      out.lkh.1 = out.lkh.1.cand
      out.lkh.2 = out.lkh.2.cand
      out.lkh.3 = out.lkh.3.cand
      out.lkh   = out.lkh.1+out.lkh.2+out.lkh.3
      out.nacc  <- out.nacc+1
    }  # else the chain remains there without jumping
    
    post.out.par[(floor(i/thin)+1),]    <- out.par     # store this information
    post.out.par.full[i,]    <- out.par
    out.lkh.call <- out.lkh.call+1
    
    # I.iii: Adapt the jump cov mat for more efficient generation of candidate output posterior parameter 
    
    if(i > adapt_par[1] && i %% adapt_par[2] == 0 && i < (adapt_par[4]*sampsize) )          # adapt the proposal covariance structure
    { #  adapt_par[1]: after which iteration to begin adaptation. adapt_par[2]: frequency with which updating occurs. 
      #  adapt_par[3]: proportion of the previous states to include when updating. adapt_par[4]: when to stop adapting.  
      out.acc.rr <- out.nacc/out.npro
      len <-floor(i*adapt_par[3]):i # select elements to include in the updating
      x <- post.out.par.full[len,] # previously sampled parameters
      N <- length(len)
      out.c <- out.acc.rr/acc.opt
      p_sigma <- (N-1) * var(x)/N  # compute/update the covariance matrix
      p_sigma <-out.c*makePositiveDefinite(p_sigma)   # deal with rounding problems that can de-symmetrize
      out.c     <- 1
      out.nacc  <- 0
      out.npro  <- 0
      if(!(0 %in% p_sigma) )   out.cand.cov<-p_sigma # to avoid jump cov degeneration
    } # end adaptation
    
    # Step II: psi, sample posterior input (observation) parameters
    # -----------------------------------------
    
    # II.i: Generate a candidate for the next sample input (observation) posterior parameter 
    inp.par.cand <- mvrnorm(1,mu=inp.par,Sigma=inp.cand.cov)  # Draw candidate point from the proposal density
    
    # II.ii: Decide whether to accept or reject the candidate input (observation) posterior parameter 
    
    boni.inp.par   <-inp.dens.1+inp.dens.2+inp.dens.3+ inp.logpri(inp.par) # use same rain 
    
    Sigma.xio.1.cand <- diag(inp.par.cand["var.xi"], nrow=length(t.grid.Ime.1), ncol=length(t.grid.Ime.1))
    Sigma.xio.2.cand <- diag(inp.par.cand["var.xi"], nrow=length(t.grid.Ime.2), ncol=length(t.grid.Ime.2))
    Sigma.xio.3.cand <- diag(inp.par.cand["var.xi"], nrow=length(t.grid.Ime.3), ncol=length(t.grid.Ime.3))
    
    inp.obs.d.1.cand.pa <-dmvnorm(xio.1, mean = traj.Inp.pr.1[ind.comp.xi.1], sigma = Sigma.xio.1.cand, log = T)
    inp.obs.d.2.cand.pa <-dmvnorm(xio.2, mean = traj.Inp.pr.2[ind.comp.xi.2], sigma = Sigma.xio.2.cand, log = T)
    inp.obs.d.3.cand.pa <-dmvnorm(xio.3, mean = traj.Inp.pr.3[ind.comp.xi.3], sigma = Sigma.xio.3.cand, log = T)
    
    # fxio|xi(psi cand)*fpsi(psi cand), i.e. observation error model * prior of the parameters ; # fxi(psi cand), i.e. prior rainfall (potential) model [not needed: these param don't modify the OU proc]
    
    boni.inp.cand.par   <- inp.obs.d.1.cand.pa +inp.obs.d.2.cand.pa+inp.obs.d.3.cand.pa + inp.logpri(inp.par.cand)
    
    inp.obs.ratio    <- boni.inp.cand.par-boni.inp.par
    
    if(!is.finite(inp.obs.ratio))   inp.obs.ratio<- -Inf      # never jump outside the range of the parameter space
    if( log(runif(1,0,1)) <= inp.obs.ratio)          # Make the jump according to acceptance probability
    { 
      inp.par       = inp.par.cand
      inp.par.nacc  = inp.par.nacc+1 
      Sigma.xio.1   = Sigma.xio.1.cand
      Sigma.xio.2   = Sigma.xio.2.cand
      Sigma.xio.3   = Sigma.xio.3.cand
    }  # else the chain remains there without jumping
    
    post.inp.par[(floor(i/thin)+1),]    <- inp.par     # store this information
    
    # Step III: xio, sample posterior input process (transf) at the obs site 
    # -----------------------------------------
    
    # III.i: Generate a sample for xio (potential rain at the obs site);apparently no rejection (all proposals are accepted)
 
    xio.1 <- norm.inp.obs.gen(ti.gr.Ime=t.grid.Ime.1,ti.gr.1mi=t.grid.1mi.1, inp.meas=meas.rain.1, 
                              ppt2norm.tra = ppt2norm.tra, 
                              par = c(inp.par,inp.par.fix), 
                              inp.elsew = traj.Inp.pr.1)
    
    xio.2 <- norm.inp.obs.gen(ti.gr.Ime=t.grid.Ime.2,ti.gr.1mi=t.grid.1mi.2, inp.meas=meas.rain.2, 
                              ppt2norm.tra = ppt2norm.tra, 
                              par = c(inp.par,inp.par.fix), 
                              inp.elsew = traj.Inp.pr.2)
    
    xio.3 <- norm.inp.obs.gen(ti.gr.Ime=t.grid.Ime.3,ti.gr.1mi=t.grid.1mi.3, inp.meas=meas.rain.3, 
                              ppt2norm.tra = ppt2norm.tra, 
                              par = c(inp.par,inp.par.fix), 
                              inp.elsew = traj.Inp.pr.3)

    # since xio changes we have to recompute this density 
    inp.dens.1  <-  dmvnorm(xio.1, mean = traj.Inp.pr.1[ind.comp.xi.1], sigma = Sigma.xio.1, log = T)
    inp.dens.2  <-  dmvnorm(xio.2, mean = traj.Inp.pr.2[ind.comp.xi.2], sigma = Sigma.xio.2, log = T)
    inp.dens.3  <-  dmvnorm(xio.3, mean = traj.Inp.pr.3[ind.comp.xi.3], sigma = Sigma.xio.3, log = T)

    
    # Step IV.1: xi, sample posterior input process (transf) at the point of interest (catchment areal average) 
    # -----------------------------------------
    
    # subintervals partitioning of the calibration domain 

leng.subI= min(round(length(t.grid.1mi.1)/2), (i-1)/(sampsize/2-1)*(interv.pts[2]-interv.pts[1])+interv.pts[1]) # avg number of points we want to consider at once (rememb: we use fine input resolution)
if (i>sampsize/2) leng.subI= interv.pts[2]

    # we assume to always partition the interval in at least 2 parts
    nr.subI   = round(length(t.grid.1mi.1)/leng.subI)
    
    break.I   = rep(NA,nr.subI-1)
    for (int in 1:(nr.subI-1)) # define break points of the intervals
    {
      break.I[int] = round(int*leng.subI + rnorm(1,0,.39*1*leng.subI)) # length of the intervals is aleatory with a fix mean
      
      if(int>1 && break.I[int]<=break.I[int-1])     break.I[int]=break.I[int-1]+5 # get a minimal interval always ensured
      
      if(break.I[int]>=length(t.grid.1mi.1)) break.I[int]=length(t.grid.1mi.1)-5
      
      if(break.I[int]<=1) break.I[int]=5
    }
    
    
    if(!all(break.I == cummax(break.I))|break.I[(nr.subI-1)]>=length(t.grid.1mi.1)) 
    {
      print("long jump!") 
      break.I   = rep(NA,nr.subI-1)
      for (int in 1:(nr.subI-1)) # define break points of the intervals
      {
        break.I[int] = round(int*leng.subI) # length of the intervals is fix 
        
        if(int>1 && break.I[int]<=break.I[int-1])     break.I[int]=break.I[int-1]+5 # get a minimal interval always ensured
        
        if(break.I[int]>=length(t.grid.1mi.1)) break.I[int]=length(t.grid.1mi.1)-5
        
        if(break.I[int]<=1) break.I[int]=5
      }
    }
    break.I = c(1,break.I,length(t.grid.1mi.1))  # include the start and ending points
    
    # loop over the subintervals to generate nr.subI partially-new trajectories conditioned on the previous path at the extremities
    # -------------------------------------------
    
    out.lkh.suI.1  = out.lkh.1  # initialize the output lkh using the previously generated and accepted rainfall but updated parameters
    
    for (l in 1:nr.subI)
    {
      cand.traj.Inp.1     = traj.Inp.pr.1
      
      if (l==1) # only condition on the end value
      {
        cand.traj.Inp.1[break.I[l]:break.I[l+1]] <- rou(tau=inp.par.fix["corrl_inp"],y1=NA,yn=traj.Inp.pr.1[break.I[l+1]],
                                                      t=as.numeric(t.grid.1mi.1[break.I[l]:break.I[l+1]]))$y
      } else {
        if (l==nr.subI) # only condition on the start value
        {
          cand.traj.Inp.1[break.I[l]:break.I[l+1]] <- rou(tau=inp.par.fix["corrl_inp"],y1=traj.Inp.pr.1[break.I[l]],yn=NA,
                                                        t=as.numeric(t.grid.1mi.1[break.I[l]:break.I[l+1]]))$y
          
        } else { # condition on both start and end values
          cand.traj.Inp.1[break.I[l]:break.I[l+1]] <- rou(tau=inp.par.fix["corrl_inp"],y1=traj.Inp.pr.1[break.I[l]],
                                                        yn=traj.Inp.pr.1[break.I[l+1]],t=as.numeric(t.grid.1mi.1[break.I[l]:break.I[l+1]]))$y        
        }
      }  
      # implement acceptance/rejection 
      
      real.Pluv.cand    <- inp.trans(cand.traj.Inp.1)*60/1000 # partially updated, partially candidate, partially old X       #      plot.out = F;        if(rnorm(1, mean = 0, sd = 1)>1.65)         plot.out = T # show some proposed rains
      
      # fyo|xi(xi cand), i.e. output error model
      out.lkh.cand.suI.1  <- out.loglikh(out.par,L=L.1 , y.obs=y.1 ,sim.inp=real.Pluv.cand ,t.grid.cal=t.grid.1mi.1 , displ.res=F) # slow because we solve an ode (can be slighly optimized)
      out.lkh.call        <- out.lkh.call+1
      # fxio|xi(xi cand), i.e. inp observation error model [alternatively: calculate density manually]
      inp.obs.cand.suI.1  <- dmvnorm(xio.1, mean = cand.traj.Inp.1[ind.comp.xi.1], sigma = Sigma.xio.1, log = T)
      
      inp.proc.ratio    <- out.lkh.cand.suI.1  + inp.obs.cand.suI.1  - (inp.dens.1 + out.lkh.suI.1 ) 
      inp.pro.npro      <- inp.pro.npro+1
      tr.inp.err.cand.1 <- sum(abs(xio.1-cand.traj.Inp.1[ind.comp.xi.1]))
      
      if( log(runif(1,0,1)) <= inp.proc.ratio )         
      {
      #################### ===============================
      # Plot proposed (blue) vs previous
      #################### ===============================
      par(mfrow=c(2,1), mar=c(0,1.5,0.1,0), oma= c(2,0,0.0,0), mgp = c(2.5, .5, 0))
      
      # i: rain #
      plot(t.grid.1mi.1,cand.traj.Inp.1, col="blue", xaxt = "n", type="l", ylim=c(-3.5,4))
      abline(v=t.grid.1mi.1[break.I], lty="dotted", col="darkorchid", lwd=.1); abline(h=inp.par.fix["thresh"], lty="dashed", col="wheat3", lwd=.1)
      par(new=T)
      plot(t.grid.1mi.1,traj.Inp.pr.1-0.00, col=gray(.61), xaxt = "n", type="l", ylim=c(-3.5,4))
      points(t.grid.1mi.1[ind.comp.xi.1],xio.1, cex=.7, pch="|")
      legend("topright", legend=paste("inp dens=", round(inp.obs.cand.suI.1, digits = 3)), bty="n", inset=c(0,.0), text.col="blue")
      legend("topright", legend=paste("inp dens=", round(inp.dens.1, digits = 3)), bty="n", inset=c(0,.08), text.col=gray(.61))
      legend("bottom", legend=paste("it=",i), bty="n")
      legend("topleft", legend=paste("SumAbsRes=",round(tr.inp.err.cand.1, digits = 3)), bty="n", text.col="blue")
      legend("topleft", legend=paste("SumAbsRes=",round(tr.inp.err.1, digits = 3)), bty="n", inset=c(0,+.08), text.col=gray(.61) )
      
      # ii: discharge #
      
      sim.inp = inp.trans(cand.traj.Inp.1)
      mod.tim.ind <- match(sysanal.decode(L.1)$val,as.numeric(t.grid.1mi.1))
      mod.inp     <- sim.inp[mod.tim.ind]*60/1000 
      y.mod       <- model.function(par=out.par,L=L.1, Inp=mod.inp, dt=dt.mod)
      NS.cand     <- 1-(sum((y.1-y.mod)^2, na.rm = T))/sum((y.1-mean(y.1, na.rm = T))^2, na.rm = T)
      plot(t.grid.Ome.1,y.mod, col="blue", type="l", ylim=c(0,90), xlim=as.numeric(c(t.grid.1mi.1[1],t.grid.1mi.1[length(t.grid.1mi.1)])))
      abline(v=t.grid.1mi.1[break.I], lty="dotted", col="darkorchid", lwd=.1)
      par(new=T)
      sim.inp = inp.trans(traj.Inp.pr.1)
      mod.inp     <- sim.inp[mod.tim.ind]*60/1000 
      y.mod       <- model.function(par=out.par,L=L.1, Inp=mod.inp, dt=dt.mod)
      NS.prev     <- 1-(sum((y.1-y.mod)^2, na.rm = T))/sum((y.1-mean(y.1, na.rm = T))^2, na.rm = T)
      plot(t.grid.Ome.1,y.mod-.0, col=gray(.61), type="l", ylim=c(0,90), xlim=as.numeric(c(t.grid.1mi.1[1],t.grid.1mi.1[length(t.grid.1mi.1)])))
      points(t.grid.Ome.1,y.1, cex=.5)
      legend("topright", legend=paste("out lkh=", round(out.lkh.cand.suI.1, digits = 3)), bty="n", inset=c(0,.0), text.col="blue")
      legend("topright", legend=paste("out lkh=", round(out.lkh.suI.1, digits = 3)), bty="n", inset=c(0,.08), text.col=gray(.61))
      legend("topleft", legend=paste("NS=",round(NS.cand, digits = 3)), bty="n", text.col="blue")
      legend("topleft", legend=paste("NS=",round(NS.prev, digits = 3)), bty="n", inset=c(0,+.08), text.col=gray(.61) )

      
        traj.Inp.pr.1[break.I[l]:break.I[l+1]]= cand.traj.Inp.1[break.I[l]:break.I[l+1]] # as l increases the input process becomes piece-wise more and more updated (some pieces will remain the same)
        inp.dens.1   = inp.obs.cand.suI.1 
        out.lkh.suI.1= out.lkh.cand.suI.1  
        inp.pro.nacc.1   <- inp.pro.nacc.1 +1
        tr.inp.err.1 = tr.inp.err.cand.1 
        legend("top", legend=paste("ACCEPTED"), bty="n", text.col="blue")
        par(new=F)
      
      }	
    } # at the end of the subcycle we have a (gradually & partially) updated input trajectory which we can store
    
    
    # Step IV.2: xi, sample posterior input process (transf) at the point of interest (center of catchment with no obs) 
    # -----------------------------------------
    
    # subintervals partitioning of the calibration domain 
leng.subI= min(round(length(t.grid.1mi.2)/2), (i-1)/(sampsize/2-1)*(interv.pts[2]-interv.pts[1])+interv.pts[1]) # avg number of points we want to consider at once (rememb: we use fine input resolution)
if (i>sampsize/2) leng.subI= interv.pts[2]


    # we assume to always partition the interval in at least 2 parts
    nr.subI   = round(length(t.grid.1mi.2)/leng.subI)
    
    break.I   = rep(NA,nr.subI-1)
    for (int in 1:(nr.subI-1)) # define break points of the intervals
    {
      break.I[int] = round(int*leng.subI + rnorm(1,0,.39*1*leng.subI)) # length of the intervals is NOT aleatory with a fix mean
      
      if(int>1 && break.I[int]<=break.I[int-1])     break.I[int]=break.I[int-1]+5 # get a minimal interval always ensured
      
      if(break.I[int]>=length(t.grid.1mi.2)) break.I[int]=length(t.grid.1mi.2)-5
      
      if(break.I[int]<=1) break.I[int]=5
    }
    
    
    if(!all(break.I == cummax(break.I))|break.I[(nr.subI-1)]>=length(t.grid.1mi.2)) 
    {
      print("long jump!") 
      break.I   = rep(NA,nr.subI-1)
      for (int in 1:(nr.subI-1)) # define break points of the intervals
      {
        break.I[int] = round(int*leng.subI) # length of the intervals is fix 
        
        if(int>1 && break.I[int]<=break.I[int-1])     break.I[int]=break.I[int-1]+5 # get a minimal interval always ensured
        
        if(break.I[int]>=length(t.grid.1mi.2)) break.I[int]=length(t.grid.1mi.2)-5
        
        if(break.I[int]<=1) break.I[int]=5
      }
    }
    break.I = c(1,break.I,length(t.grid.1mi.2))  # include the start and ending points
    
    # loop over the subintervals to generate nr.subI partially new trajectories conditionned on the previous path at the extremities
    # -------------------------------------------
    
    out.lkh.suI.2  = out.lkh.2  # initialize the output lkh using the previously generated and accepted rainfall but updated parameters
    
    for (l in 1:nr.subI)
    {
      cand.traj.Inp.2     = traj.Inp.pr.2
      
      if (l==1) # only condition on the end value
      {
        cand.traj.Inp.2[break.I[l]:break.I[l+1]] <- rou(tau=inp.par.fix["corrl_inp"],y1=NA,yn=traj.Inp.pr.2[break.I[l+1]],
                                                        t=as.numeric(t.grid.1mi.2[break.I[l]:break.I[l+1]]))$y
      } else {
        if (l==nr.subI) # only condition on the start value
        {
          cand.traj.Inp.2[break.I[l]:break.I[l+1]] <- rou(tau=inp.par.fix["corrl_inp"],y1=traj.Inp.pr.2[break.I[l]],yn=NA,
                                                          t=as.numeric(t.grid.1mi.2[break.I[l]:break.I[l+1]]))$y
          
        } else { # condition on both start and end values
          cand.traj.Inp.2[break.I[l]:break.I[l+1]] <- rou(tau=inp.par.fix["corrl_inp"],y1=traj.Inp.pr.2[break.I[l]],
                                                          yn=traj.Inp.pr.2[break.I[l+1]],t=as.numeric(t.grid.1mi.2[break.I[l]:break.I[l+1]]))$y        
        }
      }  
      # implement acceptance/rejection 
      
      real.Pluv.cand    <- inp.trans(cand.traj.Inp.2)*60/1000 # partially updated, partially candidate, partially old X       #      plot.out = F;        if(rnorm(1, mean = 0, sd = 1)>1.65)         plot.out = T # show some proposed rains
      
      # fyo|xi(xi cand), i.e. output error model
      out.lkh.cand.suI.2  <- out.loglikh(out.par,L=L.2 , y.obs=y.2 ,sim.inp=real.Pluv.cand ,t.grid.cal=t.grid.1mi.2 , displ.res=F) # slow because we solve an ode (can be slighly optimized)
      out.lkh.call        <- out.lkh.call+1
      # fxio|xi(xi cand), i.e. inp observation error model [alternatively: calculate density manually]
      inp.obs.cand.suI.2  <- dmvnorm(xio.2, mean = cand.traj.Inp.2[ind.comp.xi.2], sigma = Sigma.xio.2, log = T)
      
      inp.proc.ratio    <- out.lkh.cand.suI.2  + inp.obs.cand.suI.2  - (inp.dens.2 + out.lkh.suI.2 ) 
      inp.pro.npro      <- inp.pro.npro+1
      tr.inp.err.cand.2 <- sum(abs(xio.2-cand.traj.Inp.2[ind.comp.xi.2]))
      
      if( log(runif(1,0,1)) <= inp.proc.ratio )         
      {
        #################### ===============================
        # Plot proposed (blue) vs previous
        #################### ===============================
        par(mfrow=c(2,1), mar=c(0,1.5,0.1,0), oma= c(2,0,0.0,0), mgp = c(2.5, .5, 0))
        
        # i: rain #
        plot(t.grid.1mi.2,cand.traj.Inp.2, col="blue", xaxt = "n", type="l", ylim=c(-3.5,4))
        abline(v=t.grid.1mi.2[break.I], lty="dotted", col="darkorchid", lwd=.1); abline(h=inp.par.fix["thresh"], lty="dashed", col="wheat3", lwd=.1)
        par(new=T)
        plot(t.grid.1mi.2,traj.Inp.pr.2-0.00, col=gray(.61), xaxt = "n", type="l", ylim=c(-3.5,4))
        points(t.grid.1mi.2[ind.comp.xi.2],xio.2, cex=.7, pch="|")
        legend("topright", legend=paste("inp dens=", round(inp.obs.cand.suI.2, digits = 3)), bty="n", inset=c(0,.0), text.col="blue")
        legend("topright", legend=paste("inp dens=", round(inp.dens.2, digits = 3)), bty="n", inset=c(0,.08), text.col=gray(.61))
        legend("bottom", legend=paste("it=",i), bty="n")
        legend("topleft", legend=paste("SumAbsRes=",round(tr.inp.err.cand.2, digits = 3)), bty="n", text.col="blue")
        legend("topleft", legend=paste("SumAbsRes=",round(tr.inp.err.2, digits = 3)), bty="n", inset=c(0,+.08), text.col=gray(.61) )
        
        # ii: discharge #
        
        sim.inp = inp.trans(cand.traj.Inp.2)
        mod.tim.ind <- match(sysanal.decode(L.2)$val,as.numeric(t.grid.1mi.2))
        mod.inp     <- sim.inp[mod.tim.ind]*60/1000 
        y.mod       <- model.function(par=out.par,L=L.2, Inp=mod.inp, dt=dt.mod)
        NS.cand     <- 1-(sum((y.2-y.mod)^2, na.rm = T))/sum((y.2-mean(y.2, na.rm = T))^2, na.rm = T)
        plot(t.grid.Ome.2,y.mod, col="blue", type="l", ylim=c(0,90), xlim=as.numeric(c(t.grid.1mi.2[1],t.grid.1mi.2[length(t.grid.1mi.2)])))
        abline(v=t.grid.1mi.2[break.I], lty="dotted", col="darkorchid", lwd=.1)
        par(new=T)
        sim.inp = inp.trans(traj.Inp.pr.2)
        mod.inp     <- sim.inp[mod.tim.ind]*60/1000 
        y.mod       <- model.function(par=out.par,L=L.2, Inp=mod.inp, dt=dt.mod)
        NS.prev     <- 1-(sum((y.2-y.mod)^2, na.rm = T))/sum((y.2-mean(y.2, na.rm = T))^2, na.rm = T)
        plot(t.grid.Ome.2,y.mod-.0, col=gray(.61), type="l", ylim=c(0,90), xlim=as.numeric(c(t.grid.1mi.2[1],t.grid.1mi.2[length(t.grid.1mi.2)])))
        points(t.grid.Ome.2,y.2, cex=.5)
        legend("topright", legend=paste("out lkh=", round(out.lkh.cand.suI.2, digits = 3)), bty="n", inset=c(0,.0), text.col="blue")
        legend("topright", legend=paste("out lkh=", round(out.lkh.suI.2, digits = 3)), bty="n", inset=c(0,.08), text.col=gray(.61))
        legend("topleft", legend=paste("NS=",round(NS.cand, digits = 3)), bty="n", text.col="blue")
        legend("topleft", legend=paste("NS=",round(NS.prev, digits = 3)), bty="n", inset=c(0,+.08), text.col=gray(.61) )
        
        
        traj.Inp.pr.2[break.I[l]:break.I[l+1]]= cand.traj.Inp.2[break.I[l]:break.I[l+1]] # as l increases the input process becomes piece-wise more and more updated (some pieces will remain the same)
        inp.dens.2   = inp.obs.cand.suI.2 
        out.lkh.suI.2= out.lkh.cand.suI.2  
        inp.pro.nacc.2   <- inp.pro.nacc.2 +1
        tr.inp.err.2 = tr.inp.err.cand.2 
        legend("top", legend=paste("ACCEPTED"), bty="n", text.col="blue")
        par(new=F)
      }  
    } # at the end of the subcycle we have a (gradually & partially) updated input trajectory which we can store
    
  
    
    # Step IV.3: xi, sample posterior input process (transf) at the point of interest (center of catchment with no obs) 
    # -----------------------------------------
    
    # subintervals partitioning of the calibration domain 
leng.subI= min(round(length(t.grid.1mi.3)/2), (i-1)/(sampsize/2-1)*(interv.pts[2]-interv.pts[1])+interv.pts[1]) # avg number of points we want to consider at once (rememb: we use fine input resolution)
if (i>sampsize/2) leng.subI= interv.pts[2]
 
    # we assume to always partition the interval in at least 2 parts
    nr.subI   = round(length(t.grid.1mi.3)/leng.subI)
    
    break.I   = rep(NA,nr.subI-1)
    for (int in 1:(nr.subI-1)) # define break points of the intervals
    {
      break.I[int] = round(int*leng.subI + rnorm(1,0,.39*1*leng.subI)) # length of the intervals is NOT aleatory with a fix mean
      
      if(int>1 && break.I[int]<=break.I[int-1])     break.I[int]=break.I[int-1]+5 # get a minimal interval always ensured
      
      if(break.I[int]>=length(t.grid.1mi.3)) break.I[int]=length(t.grid.1mi.3)-5
      
      if(break.I[int]<=1) break.I[int]=5
    }
    
    
    if(!all(break.I == cummax(break.I))|break.I[(nr.subI-1)]>=length(t.grid.1mi.3)) 
    {
      print("long jump!") 
      break.I   = rep(NA,nr.subI-1)
      for (int in 1:(nr.subI-1)) # define break points of the intervals
      {
        break.I[int] = round(int*leng.subI) # length of the intervals is fix 
        
        if(int>1 && break.I[int]<=break.I[int-1])     break.I[int]=break.I[int-1]+5 # get a minimal interval always ensured
        
        if(break.I[int]>=length(t.grid.1mi.3)) break.I[int]=length(t.grid.1mi.3)-5
        
        if(break.I[int]<=1) break.I[int]=5
      }
    }
    break.I = c(1,break.I,length(t.grid.1mi.3))  # include the start and ending points
    
    # loop over the subintervals to generate nr.subI partially new trajectories conditionned on the previous path at the extremities
    # -------------------------------------------
    
    out.lkh.suI.3  = out.lkh.3  # initialize the output lkh using the previously generated and accepted rainfall but updated parameters
    
    for (l in 1:nr.subI)
    {
      cand.traj.Inp.3     = traj.Inp.pr.3
      
      if (l==1) # only condition on the end value
      {
        cand.traj.Inp.3[break.I[l]:break.I[l+1]] <- rou(tau=inp.par.fix["corrl_inp"],y1=NA,yn=traj.Inp.pr.3[break.I[l+1]],
                                                        t=as.numeric(t.grid.1mi.3[break.I[l]:break.I[l+1]]))$y
      } else {
        if (l==nr.subI) # only condition on the start value
        {
          cand.traj.Inp.3[break.I[l]:break.I[l+1]] <- rou(tau=inp.par.fix["corrl_inp"],y1=traj.Inp.pr.3[break.I[l]],yn=NA,
                                                          t=as.numeric(t.grid.1mi.3[break.I[l]:break.I[l+1]]))$y
          
        } else { # condition on both start and end values
          cand.traj.Inp.3[break.I[l]:break.I[l+1]] <- rou(tau=inp.par.fix["corrl_inp"],y1=traj.Inp.pr.3[break.I[l]],
                                                          yn=traj.Inp.pr.3[break.I[l+1]],t=as.numeric(t.grid.1mi.3[break.I[l]:break.I[l+1]]))$y        
        }
      }  
      # implement acceptance/rejection 
      
      real.Pluv.cand    <- inp.trans(cand.traj.Inp.3)*60/1000 # partially updated, partially candidate, partially old X       #      plot.out = F;        if(rnorm(1, mean = 0, sd = 1)>1.65)         plot.out = T # show some proposed rains
      
      # fyo|xi(xi cand), i.e. output error model
      out.lkh.cand.suI.3  <- out.loglikh(out.par,L=L.3 , y.obs=y.3 ,sim.inp=real.Pluv.cand ,t.grid.cal=t.grid.1mi.3 , displ.res=F) # slow because we solve an ode (can be slighly optimized)
      out.lkh.call        <- out.lkh.call+1
      # fxio|xi(xi cand), i.e. inp observation error model [alternatively: calculate density manually]
      inp.obs.cand.suI.3  <- dmvnorm(xio.3, mean = cand.traj.Inp.3[ind.comp.xi.3], sigma = Sigma.xio.3, log = T)
      
      inp.proc.ratio    <- out.lkh.cand.suI.3  + inp.obs.cand.suI.3  - (inp.dens.3 + out.lkh.suI.3 ) 
      inp.pro.npro      <- inp.pro.npro+1
      tr.inp.err.cand.3 <- sum(abs(xio.3-cand.traj.Inp.3[ind.comp.xi.3]))
      
      if( log(runif(1,0,1)) <= inp.proc.ratio )         
      {
        #################### ===============================
        # Plot proposed (blue) vs previous
        #################### ===============================
        par(mfrow=c(2,1), mar=c(0,1.5,0.1,0), oma= c(2,0,0.0,0), mgp = c(2.5, .5, 0))
        
        # i: rain #
        plot(t.grid.1mi.3,cand.traj.Inp.3, col="blue", xaxt = "n", type="l", ylim=c(-3.5,4))
        abline(v=t.grid.1mi.3[break.I], lty="dotted", col="darkorchid", lwd=.1); abline(h=inp.par.fix["thresh"], lty="dashed", col="wheat3", lwd=.1)
        par(new=T)
        plot(t.grid.1mi.3,traj.Inp.pr.3-0.00, col=gray(.61), xaxt = "n", type="l", ylim=c(-3.5,4))
        points(t.grid.1mi.3[ind.comp.xi.3],xio.3, cex=.7, pch="|")
        legend("topright", legend=paste("inp dens=", round(inp.obs.cand.suI.3, digits = 3)), bty="n", inset=c(0,.0), text.col="blue")
        legend("topright", legend=paste("inp dens=", round(inp.dens.3, digits = 3)), bty="n", inset=c(0,.08), text.col=gray(.61))
        legend("bottom", legend=paste("it=",i), bty="n")
        legend("topleft", legend=paste("SumAbsRes=",round(tr.inp.err.cand.3, digits = 3)), bty="n", text.col="blue")
        legend("topleft", legend=paste("SumAbsRes=",round(tr.inp.err.3, digits = 3)), bty="n", inset=c(0,+.08), text.col=gray(.61) )
        
        # ii: discharge #
        
        sim.inp = inp.trans(cand.traj.Inp.3)
        mod.tim.ind <- match(sysanal.decode(L.3)$val,as.numeric(t.grid.1mi.3))
        mod.inp     <- sim.inp[mod.tim.ind]*60/1000 
        y.mod       <- model.function(par=out.par,L=L.3, Inp=mod.inp, dt=dt.mod)
        NS.cand     <- 1-(sum((y.3-y.mod)^2, na.rm = T))/sum((y.3-mean(y.3, na.rm = T))^2, na.rm = T)
        plot(t.grid.Ome.3,y.mod, col="blue", type="l", ylim=c(0,90), xlim=as.numeric(c(t.grid.1mi.3[1],t.grid.1mi.3[length(t.grid.1mi.3)])))
        abline(v=t.grid.1mi.3[break.I], lty="dotted", col="darkorchid", lwd=.1)
        par(new=T)
        sim.inp = inp.trans(traj.Inp.pr.3)
        mod.inp     <- sim.inp[mod.tim.ind]*60/1000 
        y.mod       <- model.function(par=out.par,L=L.3, Inp=mod.inp, dt=dt.mod)
        NS.prev     <- 1-(sum((y.3-y.mod)^2, na.rm = T))/sum((y.3-mean(y.3, na.rm = T))^2, na.rm = T)
        plot(t.grid.Ome.3,y.mod-.0, col=gray(.61), type="l", ylim=c(0,90), xlim=as.numeric(c(t.grid.1mi.3[1],t.grid.1mi.3[length(t.grid.1mi.3)])))
        points(t.grid.Ome.3,y.3, cex=.5)
        legend("topright", legend=paste("out lkh=", round(out.lkh.cand.suI.3, digits = 3)), bty="n", inset=c(0,.0), text.col="blue")
        legend("topright", legend=paste("out lkh=", round(out.lkh.suI.3, digits = 3)), bty="n", inset=c(0,.08), text.col=gray(.61))
        legend("topleft", legend=paste("NS=",round(NS.cand, digits = 3)), bty="n", text.col="blue")
        legend("topleft", legend=paste("NS=",round(NS.prev, digits = 3)), bty="n", inset=c(0,+.08), text.col=gray(.61) )
        
        
        traj.Inp.pr.3[break.I[l]:break.I[l+1]]= cand.traj.Inp.3[break.I[l]:break.I[l+1]] # as l increases the input process becomes piece-wise more and more updated (some pieces will remain the same)
        inp.dens.3   = inp.obs.cand.suI.3 
        out.lkh.suI.3= out.lkh.cand.suI.3  
        inp.pro.nacc.3   <- inp.pro.nacc.3 +1
        tr.inp.err.3 = tr.inp.err.cand.3 
        legend("top", legend=paste("ACCEPTED"), bty="n", text.col="blue")
        par(new=F)
      }  
    } # at the end of the subcycle we have a (gradually & partially) updated input trajectory which we can store
    
    
    
 
    post.inp.1[(floor(i/thin)+1),] <- traj.Inp.pr.1
    post.inp.2[(floor(i/thin)+1),] <- traj.Inp.pr.2
    post.inp.3[(floor(i/thin)+1),] <- traj.Inp.pr.3
    log.dens[(floor(i/thin)+1),] <- c(out.post-out.lkh,out.lkh) #  aggregated
    inp.log.dens.1[(floor(i/thin)+1),] <- inp.dens.1
    inp.log.dens.2[(floor(i/thin)+1),] <- inp.dens.2
    inp.log.dens.3[(floor(i/thin)+1),] <- inp.dens.3
    if(any(i == round(sampsize*c(seq(.2,.8,.2)))))   print(paste('sampling ',i/sampsize*100,'% complete',sep='')) # info of where we are
  } # end iterations
  
  print(paste("inp pro acc rate =", sum(c(inp.pro.nacc.1,inp.pro.nacc.2,inp.pro.nacc.3))/inp.pro.npro))

  
  
  res <- list(
    post.out.par= post.out.par,
    post.inp.par= post.inp.par,
    post.inp.1  = post.inp.1,
    post.inp.2  = post.inp.2,
    post.inp.3  = post.inp.3,
    out.log.dens = log.dens,
    inp.log.dens.1 = inp.log.dens.1, # fxio|xi
    inp.log.dens.2 = inp.log.dens.2, # fxio|xi
    inp.log.dens.3 = inp.log.dens.3, # fxio|xi
    hyd.mod.run = out.lkh.call,
    out.cand.cov= out.cand.cov,
    xio.1 = xio.1,
    xio.2 = xio.2,
    xio.3 = xio.3
  )
  
  return(res)}



# ---------------------------------------------------------------------------
# Predictions of output/input with SIP in the calibration period
# ---------------------------------------------------------------------------


sysanal.predict.SIP.L1 <- function(parsamp2prop,model,L1,x.L1.samp, par.tr,inp.trans,t.grid.1mi,
                                   probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),dt.mod)
{

  # preallocate results
  
  Yo.L1.samp <- matrix(nrow=nrow(parsamp2prop),ncol=length(L1))
  colnames(Yo.L1.samp) <- L1
  rownames(Yo.L1.samp) <- rownames(parsamp2prop)
  y.L1.samp <- Yo.L1.samp 
  Yo.L1.quant <- matrix(nrow=length(probs),ncol=length(L1))
  colnames(Yo.L1.quant) <- L1
  rownames(Yo.L1.quant) <- probs
  y.L1.quant <- Yo.L1.quant
  x.L1.quant <- matrix(nrow=length(probs),ncol=length(t.grid.1mi))
  rownames(x.L1.quant) <- probs
  
  # loop over sample size
  
  for (i in 1:nrow(parsamp2prop))
  
  {
  # calculate results of deterministic model:
    
    
    sim.inp = inp.trans(x.L1.samp[i,])
    mod.tim.ind <- match(sysanal.decode(L1)$val,as.numeric(t.grid.1mi))
    mod.inp     <- sim.inp[mod.tim.ind]*60/1000 # right units for the model

    y.calc    <- model(par=parsamp2prop[i,],L=L1, Inp=mod.inp, dt=dt.mod)
  

  # transform the results: 
  
  if(!is.na(par.tr["l1"]))  
  {
    y.calc.L1.trans  <- sysanal.boxcox(y.calc,par.tr["l1"],par.tr["l2"])
  } else {
    if (!is.na(par.tr["alpha"])) 
    {
      y.calc.L1.trans  <- sysanal.logsinh(y.calc,par.tr["alpha"],par.tr["beta"])   
    } 
    else      {
      print("problem with transf param in pred")
    }
  }
  
  # add white noise: 
  
  Yo.L1.trans  <- rnorm(n=length(L1), mean = y.calc.L1.trans, sd = parsamp2prop[i,"sd.Eps_Q"])
  
  # transf back:
  
  # backtransformation if needed
  if(!is.na(par.tr["alpha"])) 
  {
    ypE.real     =sysanal.logsinh.inv(Yo.L1.trans,par.tr["alpha"],par.tr["beta"])

  }  else  {
    if(par.tr["l1"]<1)
    {
      ypE.real     =sysanal.boxcox.inv(Yo.L1.trans,par.tr["l1"],par.tr["l2"])

    } else  { print(paste("check",par.tr))
              ypE.real     = Yo.L1.trans
    }  
  }
  
  y.L1.samp[i,]  = y.calc
  Yo.L1.samp[i,] = ypE.real

  } # end of the loop over the sample length

  # compute quantiles 
  
for ( j in 1:length(L1) )  # for all time points in inference period (we have distributions given by different paths)
{
  y.L1.quant[,j]  = quantile(y.L1.samp[,j],probs=probs,na.rm=TRUE) 
  Yo.L1.quant[,j] = quantile(Yo.L1.samp[,j],probs=probs,na.rm=TRUE) 
}
for ( j in 1:length(t.grid.1mi) )  # for all time points in inference period (we have distributions given by different paths)
{
  x.L1.quant[,j] = quantile(x.L1.samp[,j],probs=probs,na.rm=TRUE)
}  
  
  return(list(  Yo.L1.samp = Yo.L1.samp,
                y.L1.samp  = y.L1.samp, 
                Yo.L1.quant = Yo.L1.quant,
                y.L1.quant = y.L1.quant,
                x.L1.quant = x.L1.quant
  ))
}


# ---------------------------------------------------------------------------
# Predictions of output/input with SIP in the validation period
# ---------------------------------------------------------------------------


sysanal.hindc.SIP.L2 <- function(parsamp2prop,model,L2,meas.rain, par.tr,inp.trans,ppt2norm.tra,t.grid.1mi,t.grid.Ime,inp.par.fix,
                                   inp.iter = 100*nrow(parsamp2prop), interv.pts=c(100,50),
                                   probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975), dt.mod,  ...)
{
  
  # infer xi in L2 
  # --------------------------
  
  
  library(MASS); library(mvtnorm)   # ... for Multivariate Normal Distribution
  inp.pro.nacc <- 0;  inp.pro.npro <- 0
  post.inp  <- array(dim=c( inp.iter,length(t.grid.1mi)));  colnames(post.inp) <- t.grid.1mi
  inp.log.dens  <- array(dim=c( inp.iter,1))
  inp.SSE       <- array(dim=c( inp.iter,1))
  colnames(inp.SSE) = c("inp.SSE")
  colnames(inp.log.dens) = c("DensInp")
  
  traj.Inp.pr <-rou(tau=inp.par.fix["corrl_inp"],y1=NA,yn=NA,t=as.numeric(t.grid.1mi))$y# draw a first realization of the OU process unconditionned 
  names(traj.Inp.pr) = t.grid.1mi
  t.grid.1mi  <-as.numeric(t.grid.1mi)
  t.grid.Ime  <-as.numeric(t.grid.Ime)
  
  inp.par  = parsamp2prop[1,]
  
  xio <- norm.inp.obs.gen(ti.gr.Ime=t.grid.Ime,ti.gr.1mi=t.grid.1mi, inp.meas=meas.rain, 
                          ppt2norm.tra = ppt2norm.tra, 
                          par = c(inp.par,inp.par.fix), 
                          inp.elsew = traj.Inp.pr)
  
  time.dist  <- abs(outer(t.grid.Ime,t.grid.Ime,"-")) #in hr
  rownames(time.dist) <- t.grid.Ime
  colnames(time.dist) <- t.grid.Ime
  
  # select the appropriate time points to compare xio and xi
  ind.comp.xi = match(t.grid.Ime,t.grid.1mi); if(any(is.na(ind.comp.xi))) stop("no match between time grids")
  
  Sigma.xio <- diag(inp.par["var.xi"], nrow=length(t.grid.Ime), ncol=length(t.grid.Ime))
  inv.Sigma.xio <- solve(Sigma.xio)
  log.det.Sum.Sigma <- determinant(Sigma.xio,logarithm=TRUE)
  inp.dens       <-   - 0.5 * length(t.grid.Ime) * log(2*pi) -
    0.5 * as.numeric(log.det.Sum.Sigma$modulus) -
    0.5 * t(xio-traj.Inp.pr[ind.comp.xi]) %*% (inv.Sigma.xio) %*% 
    (xio-traj.Inp.pr[ind.comp.xi]) 
  
  tr.inp.err <- sum(abs(xio-traj.Inp.pr[ind.comp.xi]))
  
  
  # start iterations to assimilate the input  
  
  for(i in 1:inp.iter)
  {
    if(rnorm(1, mean = 0, sd = 1)>1.95) print(i) # display current iteration just sometimes
    
    
    # Step II: psi, sample posterior input (observation) parameters
    # -----------------------------------------
    
    inp.par   <- parsamp2prop[runif(1, min = 1, max =nrow(parsamp2prop)),]
    Sigma.xio <- diag(inp.par["var.xi"], nrow=length(t.grid.Ime), ncol=length(t.grid.Ime))
    
    # Step III: xio, sample posterior input process (transf) at the obs site 
    # -----------------------------------------
    
    # III.i: Generate a sample for xio (potent rain at the obs site)
    
    xio <- norm.inp.obs.gen(ti.gr.Ime=t.grid.Ime,ti.gr.1mi=t.grid.1mi, inp.meas=meas.rain, 
                            ppt2norm.tra = ppt2norm.tra, 
                            par = c(inp.par,inp.par.fix), 
                            inp.elsew = traj.Inp.pr)
    
    # since xio changes we have to recompute this density !!
    inp.dens   <-dmvnorm(xio, mean = traj.Inp.pr[ind.comp.xi], sigma = Sigma.xio, log = T)
    
    # Step IV: xi, sample posterior input process (transf) at the point of interest (center of catchment with no obs) 
    # -----------------------------------------
    
    # subintervals partitioning of the calibration domain 
   leng.subI= min(round(length(t.grid.1mi)/2), (i-1)/(sampsize/2-1)*(interv.pts[2]-interv.pts[1])+interv.pts[1]) # avg number of points we want to consider at once (rememb: we use fine input resolution)
    if (i>sampsize/2) leng.subI= interv.pts[2]
    
    # we assume to always partition the interval in at least 2 parts
    nr.subI   = round(length(t.grid.1mi)/leng.subI)
    
    break.I   = rep(NA,nr.subI-1)
    for (int in 1:(nr.subI-1)) # define break points of the intervals
    {
      break.I[int] = round(int*leng.subI + rnorm(1,0,.39*leng.subI)) # length of the intervals is aleatory with a fix mean
      
      if(int>1 && break.I[int]<=break.I[int-1])     break.I[int]=break.I[int-1]+5 # get a minimal interval always ensured
      
      if(break.I[int]>=length(t.grid.1mi)) break.I[int]=length(t.grid.1mi)-5
      
      if(break.I[int]<=1) break.I[int]=5
    }
    
    if(!all(break.I == cummax(break.I))|break.I[(nr.subI-1)]>=length(t.grid.1mi)) 
    {
      print("long jump!"); print(break.I)
      break.I   = rep(NA,nr.subI-1)
      for (int in 1:(nr.subI-1)) # define break points of the intervals
      {
        break.I[int] = round(int*leng.subI) # length of the intervals is aleatory with a fix mean
        
        if(int>1 && break.I[int]<=break.I[int-1])     break.I[int]=break.I[int-1]+5 # get a minimal interval always ensured
        
        if(break.I[int]>=length(t.grid.1mi)) break.I[int]=length(t.grid.1mi)-5
        
        if(break.I[int]<=1) break.I[int]=5
      }
    }
    
    break.I = c(1,break.I,length(t.grid.1mi))  # include the start and ending points
    
    # loop over the subintervals to generate nr.subI partially new trajectories conditionned on the previous path at the extremities
    # -------------------------------------------
    
    for (l in 1:nr.subI)
    {
      cand.traj.Inp     = traj.Inp.pr
      if (l==1) # only condition on the end value
      {
        cand.traj.Inp[break.I[l]:break.I[l+1]] <- rou(tau=inp.par.fix["corrl_inp"],y1=NA,yn=traj.Inp.pr[break.I[l+1]],
                                                      t=as.numeric(t.grid.1mi[break.I[l]:break.I[l+1]]))$y
      } else {
        if (l==nr.subI) # only condition on the start value
        {
          cand.traj.Inp[break.I[l]:break.I[l+1]] <- rou(tau=inp.par.fix["corrl_inp"],y1=traj.Inp.pr[break.I[l]],yn=NA,
                                                        t=as.numeric(t.grid.1mi[break.I[l]:break.I[l+1]]))$y
          
        } else { # condition on both start and end values
          cand.traj.Inp[break.I[l]:break.I[l+1]] <- rou(tau=inp.par.fix["corrl_inp"],y1=traj.Inp.pr[break.I[l]],
                                                        yn=traj.Inp.pr[break.I[l+1]],t=as.numeric(t.grid.1mi[break.I[l]:break.I[l+1]]))$y        
        }
      }  
      # implement acceptance/rejection 
      
      # fxio|xi(xi cand), i.e. observation error model [alternatively: calculate density manually]
      
      inp.obs.cand.suI   <-dmvnorm(xio, mean = cand.traj.Inp[ind.comp.xi], sigma = Sigma.xio, log = T)
      
      inp.proc.ratio    <- inp.obs.cand.suI  - (inp.dens ) 
      inp.pro.npro      <- inp.pro.npro+1
      tr.inp.err.cand   <- sum(abs(xio-cand.traj.Inp[ind.comp.xi]))
      
      
      if( log(runif(1,0,1)) <= inp.proc.ratio )         
      {
        #################### ===============================
        # Plot proposed (blue) vs previous
        #################### ===============================      
        # i: rain #
        par(mfrow=c(1,1), mar=c(0,1.5,0.1,0), oma= c(2,0,0.0,0), mgp = c(2.5, .5, 0))
        
        plot(t.grid.1mi,cand.traj.Inp, col="blue", xaxt = "n", type="l", ylim=c(-3.5,4))
        abline(v=t.grid.1mi[break.I], lty="dotted", col="darkorchid", lwd=.1); abline(h=inp.par.fix["thresh"], lty="dashed", col="wheat3", lwd=.1)
        par(new=T)
        plot(t.grid.1mi,traj.Inp.pr-0.00, col=gray(.61), xaxt = "n", type="l", ylim=c(-3.5,4))
        points(t.grid.1mi[ind.comp.xi],xio, cex=.7, pch="|")
        legend("topright", legend=paste("inp dens=", round(inp.obs.cand.suI, digits = 3)), bty="n", inset=c(0,.0), text.col="blue")
        legend("topright", legend=paste("inp dens=", round(inp.dens, digits = 3)), bty="n", inset=c(0,.08), text.col=gray(.61))
        legend("bottom", legend=paste("it=",i), bty="n")
        legend("topleft", legend=paste("SumAbsRes=",round(tr.inp.err.cand, digits = 3)), bty="n", text.col="blue")
        legend("topleft", legend=paste("SumAbsRes=",round(tr.inp.err, digits = 3)), bty="n", inset=c(0,+.08), text.col=gray(.61) )
        
        
        traj.Inp.pr[break.I[l]:break.I[l+1]]= cand.traj.Inp[break.I[l]:break.I[l+1]] # as l increases the input process becomes piece-wise more and more updated (some pieces will remain the same)
        inp.dens       = inp.obs.cand.suI 
        inp.pro.nacc   = inp.pro.nacc +1
        tr.inp.err     = tr.inp.err.cand 
        
      }  
    } # at the end of the subcycle we have a (gradually & partially) updated input trajectory which we can store
    
    
    post.inp[i,] <- traj.Inp.pr
    inp.log.dens[i,] <- inp.dens
    inp.SSE[i,] <- tr.inp.err 
    
  } # end iterations to sample the posterior input
  
  print(paste("inp pro acc rate =", sum(c(inp.pro.nacc))/inp.pro.npro))
  
  
  x.L2.samp <- post.inp[(max(1,(inp.iter+1-nrow(parsamp2prop))):inp.iter),]
  
  # ---------------------------------------------------------------
  # preallocate results
  
  Yo.L2.samp <- matrix(nrow=nrow(parsamp2prop),ncol=length(L2))
  colnames(Yo.L2.samp) <- L2
  rownames(Yo.L2.samp) <- rownames(parsamp2prop)
  y.L2.samp <- Yo.L2.samp 
  Yo.L2.quant <- matrix(nrow=length(probs),ncol=length(L2))
  colnames(Yo.L2.quant) <- L2
  rownames(Yo.L2.quant) <- probs
  y.L2.quant <- Yo.L2.quant
  x.L2.quant <- matrix(nrow=length(probs),ncol=length(t.grid.1mi))
  rownames(x.L2.quant) <- probs
  
  # loop over sample size
  
  for (i in 1:min(nrow(parsamp2prop),nrow(x.L2.samp)))
    
  {
    # calculate results of deterministic model:
    
    sim.inp = inp.trans(x.L2.samp[i,])
    mod.tim.ind <- match(sysanal.decode(L2)$val,as.numeric(t.grid.1mi))
    mod.inp     <- sim.inp[mod.tim.ind]*60/1000 # right units for the model
    
    y.calc    <- model(par=parsamp2prop[i,],L=L2, Inp=mod.inp, dt=dt.mod)
    
    
    # transform the results: 
    
    if(!is.na(par.tr["l1"]))  
    {
      y.calc.L2.trans  <- sysanal.boxcox(y.calc,par.tr["l1"],par.tr["l2"])
    } else {
      if (!is.na(par.tr["alpha"])) 
      {
        y.calc.L2.trans  <- sysanal.logsinh(y.calc,par.tr["alpha"],par.tr["beta"])   
      } 
      else      {
        print("problem with transf param in pred")
      }
    }
    
    # add white noise: 
    
    Yo.L2.trans  <- rnorm(n=length(L2), mean = y.calc.L2.trans, sd = parsamp2prop[i,"sd.Eps_Q"])
    
    # transf back:
    
    # backtransformation if needed
    if(!is.na(par.tr["alpha"])) 
    {
      ypE.real     =sysanal.logsinh.inv(Yo.L2.trans,par.tr["alpha"],par.tr["beta"])
      
    }  else  {
      if(par.tr["l1"]<1)
      {
        ypE.real     =sysanal.boxcox.inv(Yo.L2.trans,par.tr["l1"],par.tr["l2"])
        
      } else  { print(paste("check",par.tr))
                ypE.real     = Yo.L2.trans
      }  
    }
    
    y.L2.samp[i,]  = y.calc
    Yo.L2.samp[i,] = ypE.real
    
  } # end of the loop over the sample length
  
  
  # compute quantiles 
  
  for ( j in 1:length(L2) )  # for all time points in inference period (we have distributions given by different paths)
  {
    y.L2.quant[,j]  = quantile(y.L2.samp[,j],probs=probs,na.rm=TRUE) 
    Yo.L2.quant[,j] = quantile(Yo.L2.samp[,j],probs=probs,na.rm=TRUE) 
  }
  for ( j in 1:length(t.grid.1mi) )  # for all time points in inference period (we have distributions given by different paths)
  {
    x.L2.quant[,j] = quantile(x.L2.samp[,j],probs=probs,na.rm=TRUE)
  }  
  
  return(list(  Yo.L2.samp = Yo.L2.samp, # (real output space)
                y.L2.samp  = y.L2.samp,
                x.L2.quant = x.L2.quant,
                Yo.L2.quant = Yo.L2.quant,
                y.L2.quant = y.L2.quant,
                x.L2.samp = x.L2.samp
  ))
}

# ===============================================================================
# Draw from Ornstein-Uhlenbeck process
# ===============================================================================

# First version: Peter Reichert, Jan. 25, 2014
# Last revision: Peter Reichert, June 24, 2014

# This function draws a realization of an Ornstein-Uhlenbeck process.

# Arguments: 
#   mean:     Asymptotic mean of the process.
#   sd:       Asymptotic standard deviation of the process.
#   tau:      Correlation time of the process.
#   y1:       Optional starting value of the process (if NA, the value is drawn 
#             from the asymptotic distribution).
#   yn:       Optional end value of the process (if not NA, and n = length(t) > 1
#             the process is conditioned to this value at the last time point)
#   t:        Time points at which the process should be sampled.
#             (Note: the value at t[1] will be the starting value, y1, if 
#             provided, the value at t[n] will be the end value, yn, if provided
#             and n is the number of time points in t.)
#   log:      Indicator whether the log of the variable should be an
#             Ornstein-Uhlenbeck process (log=TRUE) rather than the
#             variable itself.
#             (Note: mean and sd are interpreted in original units
#             also for log=TRUE.)

rou <- function(mean=0,sd=1,tau=0.1,y1=NA,yn=NA,t=0:1000/1000,log=FALSE)
{
  # consistency checks:
  
  n <- length(t)
  if ( n > 1 )
  {
    if ( min(diff(t)) <= 0 )
    {
      warning("t must be strictly increasing")
      return(NA)
    }
  }
  if ( sd < 0 )
  {
    warning("sd must be non-negative")
    return(NA)
  }
  if ( tau <= 0 )
  {
    warning("tau must be positive")
    return(NA)
  }
  
  # calculate exp of log if log=TRUE (note: mean and sd are in original units):
  
  if ( log ) 
  {
    if ( mean <= 0 )
    {
      warning("if log=TRUE, mean must be positive")
      return(NA)
    }
    meanlog <- log(mean/(sqrt(1+sd*sd/(mean*mean))))
    sdlog   <- sqrt(log(1+sd*sd/(mean*mean)))
    res <- rou(mean = meanlog,
               sd   = sdlog,
               tau  = tau,
               y1   = log(y1),
               yn   = log(yn),
               t    = t,
               log  = FALSE)
    res$y <- exp(res$y)
    return(res)
  }
  else
  {  
    # draw and return sample:
    
    y <- rep(NA,n)
    if ( !is.na(y1) ) 
    {
      y[1] <- y1
    }
    else
    {
      if ( !is.na(yn) )
      {
        y[1] <- rnorm(n    = 1,
                      mean = mean + (yn-mean)*exp(-(t[n]-t[1])/tau),
                      sd   = sd * sqrt(1-exp(-2*(t[n]-t[1])/tau)))
      }
      else
      {
        y[1] <- rnorm(n=1,mean=mean,sd=sd)
      }
    }
    
    if ( n > 1 )
    {
      for ( i in 2:n )
      {
        if ( is.na(yn) ) # if we want to condition only on the last first
        {
          y[i] <- rnorm(n    = 1,
                        mean = mean + (y[i-1]-mean)*exp(-(t[i]-t[i-1])/tau),
                        sd   = sd * sqrt(1-exp(-2*(t[i]-t[i-1])/tau)))
        }
        else # if we want to condition on the last value
        {
          y[i] <- rnorm(n    = 1,
                        mean = mean + (y[i-1]-mean)*
                          exp(-(t[i]-t[i-1])/tau)*(1-exp(-2*(t[n]-t[i])/tau))/
                          (1-exp(-2*(t[n]-t[i-1])/tau))+
                          (yn-mean)*exp(-(t[n]-t[i])/tau)*(1-exp(-2*(t[i]-t[i-1])/tau))/ # addiert
                          (1-exp(-2*(t[n]-t[i-1])/tau)), # addiert
                        sd   = sd * sqrt((1-exp(-2*(t[i]-t[i-1])/tau))*
                                           (1-exp(-2*(t[n]-t[i])/tau))/
                                           (1-exp(-2*(t[n]-t[i-1])/tau))))
          #          browser()
        }
      }
    }
    return(list(x=t,y=y))
  }
}


# ===============================================================================
# Calculate probability density of an Ornstein-Uhlenbeck process
# ===============================================================================

# First version: Peter Reichert, Nov. 05, 2014
# Last revision: Peter Reichert, Nov. 05, 2014

# This function returns the height of the probability distribution 
# of an Ornstein-Uhlenbeck process at each point.

# Arguments: 
#   t:        time points of a time series.
#   y:        time series corresponding to the points in time t.
#   mean:     Asymptotic mean of the process.
#   sd:       Asymptotic standard deviation of the process.
#   tau:      Correlation time of the process.
#   log:      If true, the log of the density is returned)


dou <- function(t,y,mean=0,sd=1,tau=0.1,y1=NA,yn=NA,log=TRUE)
{
  n <- length(t)
  if ( length(y) != n ) return(NA)
  
  log.dou <- dnorm(x=y[1],mean=mean,sd=sd,log=TRUE)
  for ( i in 1:(n-1) )
  {
    log.dou <- log.dou + dnorm(x    = y[i+1],
                               mean = mean+(y[i]-mean)*exp(-(t[i+1]-t[i])/tau),
                               sd   = sd*sqrt(1-exp(-2*(t[i+1]-t[i])/tau)),
                               log  = TRUE)
  }
  
  if ( log ) return(log.dou)
  return(exp(log.dou))
}

# #################################################################
# Adaptive MCMC sampling based on Haario et al. (2001)
# #################################################################

Metro_Hastings_0815 <-
  function(li_func,pars,prop_sigma=NULL,iterations=50000,burn_in=1,adapt_par=c(100,100,0.5,0.45),quiet=FALSE,...)
  {  library("MHadaptive")
    
    if (!is.finite(li_func(pars, ...))) 

      stop("Seed parameter values <pars> are not in the defined parameter space.  Try new starting values for <pars>.")
    
    if(!is.null(dim(prop_sigma)))
    {
      if( ( dim(prop_sigma)[1] != length(pars) ||  dim(prop_sigma)[2] != length(pars) ) && !is.null(prop_sigma) )
        stop("prop_sigma not of dimension length(pars) x length(pars)")
    }
    
    if(is.null(prop_sigma)) #if no proposal matrix given, estimate the Fisher information to use as the diagonal (start in the right variance scale)
    {
      if(length(pars)!=1)
      {
        fit<-optim(pars,li_func,control=list("fnscale"=-1),method = "Nelder-Mead", hessian=TRUE,...)
        fisher_info<-solve(-fit$hessian)
        prop_sigma<-sqrt(diag(fisher_info))
        prop_sigma<-diag(prop_sigma)
      }else{
        prop_sigma<-1+pars/2
      }
    }
    
    prop_sigma<-makePositiveDefinite(prop_sigma)
    mu<-pars
    pi_X<-li_func(pars,...)			        # initial likelihood evaluation
    k_X<-pars
    trace          <-array(dim=c(iterations,length(pars)))
    # -*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*-  
    colnames(trace)    = names(pars)
    log.post           = array(dim=c(iterations,1))
    colnames(log.post) = "log.post"
    acc.opt           = .3
    # -*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*-  
    deviance<-array(dim=iterations) 
    announce<-floor(seq(iterations/10,iterations,length.out=10))
    nacc  <- 0
    npro  <- 0

    for(i in 1:iterations)
    {
      
      npro  <- npro+1
      k_Y<-mvrnorm(1,mu=k_X,Sigma=prop_sigma)	# Draw proposal point
      
      pi_Y<-li_func(k_Y,...)		        # evaluate likelihood at proposal point
      
      a_X_Y = (pi_Y)-(pi_X)		        # Compare relative likelihoods
      if(!is.finite(a_X_Y))           # -*--*--*--*--*--*--*--*--*-
        a_X_Y<--Inf                          # never jump outside the range of the parameter space
      if( log(runif(1,0,1)) <= a_X_Y)	        # Make the jump according to MH probability
      {
        k_X = k_Y
        pi_X = pi_Y
        nacc  <- nacc+1
      }	
      
      trace[i,]<-k_X	
      log.post[i,] <-(pi_X)   # -*--*--*--*--*--*--*--*--*-
      deviance[i]<-(-2*pi_X)                      # Store the deviance for calculating DIC
      if(i > adapt_par[1] && i %% adapt_par[2] == 0 && i < (adapt_par[4]*iterations) )	        # adapt the proposal covariance structure
      {   
        acc.rr  <- nacc/npro
        len<-floor(i*adapt_par[3]):i
        x<-trace[len,]
        N<-length(len)
                    
        c = acc.rr/acc.opt

        p_sigma <- (N-1) * var(x)/N
        p_sigma <-c*makePositiveDefinite(p_sigma)   # To deal with rounding problems that can de-symmetrize
        c     <- 1
        nacc  <- 0
        npro  <- 0
        if(!(0 %in% p_sigma) ) 
          prop_sigma<-p_sigma
      }
      if(!quiet && i %in% announce)
        print(paste('updating: ',i/iterations*100,'%',sep=''))
    }
    trace<-trace[burn_in:iterations,]
    DIC<-NULL
    ## Calculate the DIC
    D_bar<-mean(deviance[burn_in:iterations])
    if(length(pars)>1)
    {
      theta_bar<-sapply(1:length(pars),function(x){mean( trace[,x] )})
    }else
      theta_bar<-mean( trace )
    
    D_hat<-li_func(theta_bar,...)
    pD<-D_bar-D_hat
    DIC<-D_hat + 2*pD
    
    accept_rate<-length(unique(trace[,1]))/(iterations-burn_in) 

    val<-list("trace"=trace,"prop_sigma"=prop_sigma,"DIC"=DIC,'acceptance_rate'=accept_rate,
              "unnor.log.post"=log.post)
    class(val)<-"MHposterior"
    return(val)
  }


# ###############################################
# Lumped conceptual hydrological model
# ###############################################

lin_1box_ww <- function(par,L, Inp,dt)
{
  # model:
  # ------
  #
  #   flow through a 1 cascade reservoir model:
  #
  #   dS1                1
  #   ----  =  A * P  -  -- * S1 + bf
  #    dt                K
  #
  #             1                2*pi*t            4*pi*t            2*pi*t            4*pi*t
  #   flow  =  -- * S1 + si1*sin(------) + si2*sin(------) + co1*cos(------) + co2*cos(------)
  #             K                  24                24                24                24
  #
  # state variables:
  # ----------------
  #
  #   S1      cascade state 1
  #
  # parameters:
  # -----------
  #
  #   A       catchment area [m2]
  #   K       reservoir residence time [hr]
  #   si1,si2,co1,co2     parameters of trigonometric function to describe dry weather variation
  #   P       Inp [m/hr]
  #   bf      constant groundwater inflow in the sewer [m3/hr]
  
  # decode layout:
  # --------------
  time=sysanal.decode(L)$val
  Dt = time[2]-time[1]
  
  co1=-par["nco1"]; # if (is.na(co1)) print("par missing")
  co2=par["co2"];   # if (is.na(co2)) print("par missing")
  si1=-par["nsi1"]; # if (is.na(si1)) print("par missing")
  si2=-par["nsi2"]; # if (is.na(si2)) print("par missing")
  
  WWF = si1*sin(2*pi*time/24)+si2*sin(4*pi*time/24)+co1*cos(2*pi*time/24)+co2*cos(4*pi*time/24)
  
  # STEP I: start simulations inustus_temp times earlier
  # burn-in phase
  #---------------
  inustus_temp = 20 # time steps
  #   dt i.e. dt*60 min (t.s. of rain and calculation)
  #   dt.o=4/60
  start=max(dt, time[1]-inustus_temp*(dt))
  end  =time[length(time)]
  time =seq(start,end, by=dt)
  if ((length(time)-inustus_temp)<length(L))
  {
    end = end + dt*(length(L)-(length(time)-inustus_temp))
    time =seq(start,end, by=dt)
  }

  Inp.elon  =c(rep(0,inustus_temp),Inp) # m/h
  
  # calculate results:
  # ------------------
  S1_ini <- 4; # close to dry weather analitical sol
  
  K=par["k"];       # if (is.na(K)) print("par missing")

  Inp.fl  = par["A"]*Inp.elon+par["bf"] # m3/h
  
  S0 <- S1_ini
  S  <- rep(NA,length(time))
  SWF<- S
  
  # iii
  for (t in 1:length(time))
  {
    
    S[t]  <- Inp.fl[t]*K + (S0-Inp.fl[t]*K)*exp(-1/K*Dt) # analitical sol linear reservoir
    
    SWF[t]<- Inp.fl[t] - (S[t]-S0)/Dt  # m3/min
    
    S0    <- S[t] 
    
  }
  
  flow=  SWF[(inustus_temp+1):length(SWF)]/3.6 + WWF #l/s
  res <- as.numeric(flow)
  res[which(res<0)]  <- 1e-7 # eliminate neg val
  # save the output at the same dt as flow data
  names(res) <- L

  return(res)
}
