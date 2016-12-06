fast.smacof <- function(w1,delta, ndim = 2,verbose = FALSE, init=NULL,
                        modulus = 1, itmax = 1000, eps = 1e-6)  
{
    require(shapes)
    require(protoclass)
    #require(smacof)
    ## w1 is the value such that W is then of the form:
    ##       ee^T-I on diagonal blocks
    ##       w1 I off diagonal
    ## delta ... matrix or list of dissmilarity matrices to be jointly embedded
    ##           vertices are assumed matched

    m<-length(delta)
    nvert<-nrow(delta[[1]])
  # ndim ... number of dimensions
  # init ... matrix with starting values of dimension nvert*m \times ndim
  # type ... either "ratio", "interval", "ordinal", "spline" (replaces metric)
  # ties ... ties for pava (primary, secondary, tertiary)
  # modulus ... modulus for nonmetric update
  # itmax ... maximum number of iterations
  # eps ... change in loss function
  # spline.degree ... degree of the spline in case a spline transformation is chosen
  # spline.intKnots ... number of interior knots for the spline in case a spline transformation is chosen
  
  ## --- sanity checks
  # make the delta's into diss objects, 
  # can parallelize here
    diss<-delta
    for(i in 1:m){
        if ((is.matrix(diss[[i]])) || (is.data.frame(diss[[i]]))) {
            diss[[i]] <- strucprep(diss[[i]])  #if data are provided as dissimilarity matrix
            attr(diss[[i]], "Labels") <- c((nvert*(i-1)+1):(nvert*i))
        }
    }
    
    p <- ndim                                     
    n <- nvert*m
    NN<-n*(n-1)/2
    if (p > (n - 1)) stop("Maximum number of dimensions is nvert*m-1!")
    
    nn <- nvert*(nvert-1)/2
    m <- length(diss)
    
    wgthsd<-as.dist(matrix(1,nvert,nvert)-diag(nvert))
    wgthso<-w1*diag(nvert)
    
  # can parallelize here
    dhat<-list()
    for(i in 1:m){
        dhat[[i]] <- normDissN(diss[[i]], wgthsd, 1)*sqrt(NN/(m*nn))  #normalize each dissimilarity to nvert(nvert-1)/2
    }
    if (is.null(init)){   
  # to initialize, instead of embedding everything at once using torgerson
  # which would be inefficient speed-wise, let's embed each separately
  # and procrustes fit them to each other and call that the initialization?
  # this can be done more elegantly, but i think this is quick and dirty
        x<-list()
        xmean <- torgerson(sqrt(Reduce('+',diss)), p=p)
        for(i in 1:m){
            x[[i]] <- torgerson(sqrt(diss[[i]]), p=p) # intialize each x
            xx     <- procOPA(xmean,x[[i]],scale=FALSE,reflect=TRUE)
            x[[i]] <- xx$Bhat
        } 
    } else{
        x <- list()
        for(i in 1:m){x[[i]]<-as.matrix(init)[((i-1)*nvert+1):(nvert*i),]}
    }
  
    itel <- 1  #iteration number
  # make the distance matrix list (can be parallelized)
    d    <- list()
    for (i in 1:m){
        d[[i]]<-list()
        for (j in i:m){
            d[[i]][[j]]<-dist2(x[[i]],x[[j]])
            if (i==j){d[[i]][[j]]<-as.dist(d[[i]][[j]])}
        }
    }
    numer   <-0
    denom   <-0
    for(i in 1:m){numer <-numer+sum(d[[i]][[i]]*dhat[[i]])}
    for(i in 1:m){
        for(j in i:m){
            if(i==j){denom<-denom+sum(wgthsd*d[[i]][[i]]^2)}
            else{denom<-denom+sum(w1*diag(d[[i]][[j]]^2))}
        }
    }
    
    lb   <- numer/denom      #normalization tr(X'VX); 
    for(i in 1:m){
        x[[i]]    <- lb*x[[i]]      #modify x with lb-factor
        for(j in i:m){
            d[[i]][[j]]    <- lb*d[[i]][[j]]  #modify d with lb-factor
        }
    }
    c1   <- (nvert+w1)/(nvert*(nvert+m*w1))
    c2   <- w1/(nvert*(nvert+m*w1))
    
  # calculate stress (can be parallelized)
    sold<-0
    for(i in 1:m){
        for(j in i:m){
            if(i==j){sold<-sold+sum(wgthsd*(dhat[[i]]-d[[i]][[j]])^2)}
            else{sold<-sold+sum(w1*diag(d[[i]][[j]]^2))}
        }
    }
    sold <- sold/(NN)         #stress (to be minimized in repeat loop)
    
  #--------------- begin majorization --------------------
    repeat {    #majorization loop
    # make B*x's
        bx <- list()
        for(i in 1:m){
            b <- bmat(dhat[[i]],wgthsd,d[[i]][[i]])
            bx[[i]] <- b%*%x[[i]]
        }
    # make updated embedding
        y <- list()
        ysum<-c2*Reduce("+",bx)
        for(i in 1:m){
            y[[i]]<-ysum+(c1-c2)*bx[[i]]
        }
    # make distance matrices for y
        e    <- list()
        for (i in 1:m){
            e[[i]]<-list()
            for (j in i:m){
                if (i==j){e[[i]][[j]]<-dist(y[[i]])}
                if(i != j){e[[i]][[j]]<-dist2(y[[i]],y[[j]])}}
        }

    
    # calculate new stress (can be parrallelized)
        snon<-0
        for(i in 1:m){
            for(j in i:m){
                if(i==j){snon<-snon+sum(wgthsd*(dhat[[i]]-e[[i]][[j]])^2)}
                else{snon<-snon+sum(w1*diag(e[[i]][[j]])^2)}
            }
        }
        snon <- snon/NN
        
    #print out intermediate stress
        if (verbose) cat("Iteration: ",formatC(itel,width=3, format="d"),
                         " Stress (normalized):", formatC(c(snon),digits=8,width=12,format="f"),
                         " Difference: ", formatC(sold-snon,digits=8,width=12,format="f"),"\n")
    
        if (((sold-snon)<eps) || (itel == itmax)) break()
        x <- y                           #update configurations
        d <- e                           #update configuration distances
        sold <- snon                     #update stress
        itel <- itel+1	                 #increase iterations
    }
  #------------------ end majorization --------------- 
  
    stress <- sqrt(snon)               #stress normalization
  
    confdiss<-e
    for(i in 1:m){
        for(j in i:m){
            if(i==j){confdiss<-normDissN(e[[i]][[j]],wgthsd,1)}
            else{confdiss<-normDissN(e[[i]][[j]],wgthso,1)}
        }
    }
    
    if (itel == itmax) warning("Iteration limit reached! Increase itmax argument!") 
  
  #return configurations, configuration distances, normalized observed distances 
    result <- list(delta = diss, obsdiss = dhat, confdiss = confdiss, conf=y, stress=stress, niter = itel)
    result 
}




strucprep <-function(x)
{
  distvec <- as.vector(x[lower.tri(x)])
  n <- dim(x)[1]
  dissim <- structure(distvec, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE)
}

transform <-function (Target, x, w = rep(1,length(x$x)), normq = 0){
  # min ||Result-Target|| s.t. transformation 
  #    x   object of type "optScal" (S3 class)
  #    Target: unconstrained vector of target values
  #    w: vector nonnegative weights
  #
  #    x$trans=none     no transformation 
  #          linear   linear transformation
  #          nominal  nominal transformation
  #          ordinalp ordinal primary approach to ties (untie ties)
  #          ordinals secondary approach to ties (keep ties tied)
  #          ordinalt tertiary approach to ties (keep ties tied)
  #          spline   I-spline transformation
  #          mspline  monotone I-spline transformation
  #          power    use dis, do not change anything
  #    x$missing: 
  #       none     missing values (NA) are deleted, that is, their 
  #                weights w are considered to be 0
  #       single   missing values (NA) are considered a single category that does not need to 
  #                follow the restrictions of trans
  #       multiple each missing value (NA) receives its own category and does follow the 
  #                restrictions of trans
  
  #    normq >0       sum of squares equal to normq
  n <- length(x$x)
  b <- NULL
  iord3 <- x$iord
  Result <- rep(0,n)
  if (x$missing == "none") {     # Only operate on nonmissings
    ind_act <- x$iord_nonmiss
    nties_act <- x$nties_nonmis
  } else if(x$missing %in% c("single","multiple")){ # Last tieblock contains missings
    ind_act <- x$iord
    nties_act <- length(x$ties)
  }
  
  y  <- rep(0,nties_act)   
  w2 <- rep(0,nties_act)
  Temp <-                       # Make y as the weighted mean (in order of x$iord)
    .C("weightedMean",
       as.numeric(y), 
       as.numeric(w2), 
       as.numeric(Target), 
       as.numeric(w), 
       as.integer(x$iord), 
       as.integer(x$ties), 
       as.integer(n), 
       as.integer(nties_act))
  y <- Temp[[1]]
  w2 <- Temp[[2]]
  
  if (n > x$n_nonmis & (x$missing %in% c("single","multiple"))) {                 # Fill the estimates 
    Result[x$iord_mis] <- y[(x$nties_nonmis+1):length(x$ties)]
  }
  if (x$trans == "none") {
    # no transformation
    Result[x$iord_nonmis] <- x$x[x$iord_nonmis]
  } else if (x$trans %in% c("linear","interval","mlinear","minterval","mspline","spline")){ 
    # linear transformation subject to pos constraints
    ind_ties_nonmis <- 1:x$nties_nonmis
    w3   <- w2[ind_ties_nonmis]
    y3   <- y[ind_ties_nonmis]
    #A    <- crossprod(x$base,(matrix(w3,x$nties_nonmis,2)*x$base))
    #b    <- crossprod(x$base,w3*y3)
    #z <- solve(A,b)  # linear least-squares solution, no positivity constraints
    #z <- nnls(A,b)   # linear least-squares solution with positivity constraints
    
    ncoef <- ncol(x$base)
    A <- matrix(w3^.5,x$nties_nonmis,ncoef)*x$base
    f <- w3^.5*y3
    if (x$trans %in% c("spline","interval","linear")) {
      f <- crossprod(A,f)
      A <- crossprod(A)
      b <- solve(A,f)         # nonmonotone spline, interval, or linear
    } else {
      nnls.res <- nnls(A,f)   # monotone spline, interval, or linear
      b <- nnls.res$x
    }    
    Result[x$iord_nonmis] <- rep(x$base %*% b,x$ties[ind_ties_nonmis])
  } else if (x$trans %in% c("nominals","nominal")) {
    # nominal transformation secondary approach to ties (keep ties tied)
    Result[x$iord_nonmis] <- rep(y,x$ties[1:nties_act])
  } else if (x$trans == "ordinalp") {
    # ordinal transformation primary approach to ties (untie ties)
    iord3 <- order(x$x,Target)
    if (n>x$n_nonmi){
      iord3_nonmis <- iord3[-((x$n_nonmis+1):n)]
    } else {
      iord3_nonmis <- iord3
    }
    Temp  <- .C("wmonreg", as.numeric(Target[iord3_nonmis]), 
                as.numeric(w[iord3_nonmis]), 
                as.integer(x$n_nonmis))
    Result[iord3_nonmis] <- Temp[[1]]
  } else if (x$trans %in% c("ordinals","ordinalt","ordinal")) {
    # ordinal transformation secondary approach to ties (keep ties tied)
    #--------
    Temp <- .C("wmonreg", as.numeric(y), 
               as.numeric(w2), 
               as.integer(x$nties_nonmis)) 
    ycon   <- Temp[[1]]
    ind_ties_nonmis <- 1:x$nties_nonmis
    if (x$trans %in% c("ordinals","ordinal")) {    
      Result[x$iord_nonmis] <- rep(ycon[ind_ties_nonmis],x$ties[1:x$nties_nonmis])      
    } else { # trans == "ordinalt", tertiary approach to ties
      Result[x$iord_nonmis] <- Target[x$iord_nonmis] + 
        rep(ycon[ind_ties_nonmis] - y[ind_ties_nonmis], x$ties[1:x$nties_nonmis])      
    }
  } 
  
  if (normq>0){             # Normalize to length normq
    Result <- Result*(normq/sum(w*Result^2))^(.5)
  }
  return(list(res = Result, b = b, iord.prim = iord3))
}
transPrep <- function (x, trans = "ordinals", spline.intKnots = 4, spline.degree = 2,
                       missing = "none"){
  # Prepares for transformation 
  # and gives an initial form of the data
  #    trans:
  #       none     no transformation 
  #       linear   linear transformation
  #       interval   linear transformation
  #       nominal  nominal transformation
  #       ordinalp ordinal primary approach to ties (untie ties)
  #       ordinals secondary approach to ties (keep ties tied)
  #       ordinalt tertiary approach to ties (keep ties tied)
  #       spline   I-spline transformation
  #       mspline  monotone I-spline transformation
  #    missing: 
  #       none     missing values (NA) are deleted, that is, their 
  #                weights w are considered to be 0
  #       single   missing values (NA) are considered a single category that does not need to 
  #                follow the restrictions of trans
  #       multiple each missing value (NA) receives its own category and does follow the 
  #                restrictions of trans
  #            
  #
  n  <- length(x)
  knotSeq <- NULL
  base <- NULL
  #xSorted <- sort(as.vector(x), index.return = TRUE, na.last = TRUE)
  iord <- order(as.vector(x),  na.last = TRUE)  # Missing values are ordered last
  y    <- as.vector(x)[iord]
  y[is.na(y)] <- Inf       # Replace NA by Inf for the computation of tie blocks
  # Find tie blocks
  indTieBlock <-c(1,(2:n)[!y[-n]==y[-1]])
  ties <- c(indTieBlock[-1],n+1)-indTieBlock
  # Determine number of nonmissing data
  n_nonmis <- ifelse(is.infinite(y[n]), n - ties[length(ties)], n)
  iord_nonmis <- iord[1:n_nonmis]             # Order permutation for nonmissings
  if (n_nonmis < n){                          # If there are missings
    iord_mis <- iord[(n_nonmis+1):n]
    nties_nonmis <- length(ties) - 1
    if (missing=="multiple"){                 # add tieblocks of 1 for each missing value
      ties <- c(ties[-length(ties)],rep(1,n-n_nonmis))
    }
  } else {
    nties_nonmis <- length(ties)
    iord_mis <- NULL
  }
  x_unique <- x[iord[cumsum(ties[1:nties_nonmis])]]
  
  # Set xInit initial 
  if (trans %in% c("none","linear","interval")) {
    base <- cbind(rep(1,nties_nonmis),x_unique-x_unique[1])
    xInit <- rep(0,n)
    xInit[iord_nonmis] <- rep(x_unique,ties[1:nties_nonmis])
  } else if (trans %in% c("ordinalp","ordinals","ordinalt","ordinal","nominal","nominals","nominalp")){
    i <- 1:nties_nonmis
    xInit <- rep(0,n)
    xInit[iord_nonmis] <- rep(i,ties[1:nties_nonmis])
  } else if (trans %in% c("spline","mspline")){    
    if (ties[1]!=n){
      res <- splineSetUp(x_unique, spline.intKnots + 2, spline.degree)
      #base <- matrix(0, n, ncol(res$base))
      #base[iord[1:n_nonmis],] <- res$base
      base    <- res$base
      knotSeq <- res$knotSeq
    } else {
      base <- matrix(1, n_nonmis, 1)
    }
    xInit <- rowSums(base)
  }
  return(list(x = x, 
              x_unique = x_unique,
              n = n,
              n_nonmis = n_nonmis,
              trans = trans, 
              spline.allKnots = spline.intKnots + 2, 
              spline.degree = spline.degree, 
              spline.knotSeq = knotSeq,
              xInit = xInit, 
              iord = iord, 
              ties = ties, 
              nties_nonmis = nties_nonmis,
              base = base, 
              missing = missing,
              iord_nonmis = iord_nonmis,
              iord_mis = iord_mis,
              #factor = fact,
              class = "optScal"))
}
normDiss <-function(diss,wghts)
{
  return(diss/sqrt(sum(wghts*diss^2)))
}
normDissN <-function(diss,wghts,m)
{
  N <- length(diss)*m
  dissnorm <- diss/sqrt(sum(wghts*diss^2))*sqrt(N)
  return(dissnorm)
}
summary.smacofSP <-function(object, ...)
  {
    cat("\n")
    cat("Configurations:\n")
    print(round(object$conf,4))
    #cat("\nConfiguration dissimilarities: \n")
    #print(round(object$confdiss,4))
  }
summary.smacofR <-function(object, ...)
  {
    cat("\n")
    cat("Subjects configurations (rows):\n")
    print(round(object$conf.row,4))
    cat("\n")
    cat("Objects configurations (columns):\n")
    print(round(object$conf.col,4))
    
    cat("\n\n")
    cat("Stress per point rows:\n")
    spp.perc.row <- object$spp.row/sum(object$spp.row)*100
    sppmat.row <- cbind(sort(object$spp.row), sort(spp.perc.row))
    colnames(sppmat.row) <- c("SPP","SPP(%)")
    print(round(sppmat.row, 4))
    cat("\n")
    
    cat("Stress per point columns:\n")
    spp.perc.col <- object$spp.col/sum(object$spp.col)*100
    sppmat.col <- cbind(sort(object$spp.col), sort(spp.perc.col))
    colnames(sppmat.row) <- c("SPP","SPP(%)")
    print(round(sppmat.row, 4))
    cat("\n") 
  }
summary.smacofID<-function(object, ...)
  {
    cat("\n")
    cat("Group Stimulus Space (Joint Configurations):\n")
    print(round(object$gspace, 4))
    
    cat("\n\n")
    cat("Stress per point:\n")
    
    spp.perc <- object$spp/sum(object$spp)*100
    sppmat <- cbind(sort(object$spp), sort(spp.perc))
    colnames(sppmat) <- c("SPP","SPP(%)")
    print(round(sppmat, 4))
    cat("\n")
    
  }
summary.smacofB<-function(object, ...)
  {
    cat("\n")
    cat("Configurations:\n")
    print(round(object$conf,4))
    
    cat("\n\n")
    cat("Stress per point:\n")
    
    spp.perc <- object$spp/sum(object$spp)*100
    sppmat <- cbind(sort(object$spp), sort(spp.perc))
    colnames(sppmat) <- c("SPP","SPP(%)")
    print(round(sppmat, 4))
    cat("\n")
  }
initWeights <-function(diss) {
    if (!is.list(diss)) {
      n<-attr(diss,"Size")
      return(as.dist(matrix(1,n,n)))
    }
    n <- attr(diss[[1]],"Size"); m<-length(diss)
    return(repList(as.dist(matrix(1,n,n)),m))
  }
torgerson<-function(diss, p=p)
  {
    #diss ... dissimilarity matrix
    #p ... number of dimensions
    
    #------------------- begin subroutines -------------------
    #Torgerson's double centering
    doubleCenter <- function(x) {
      n <- dim(x)[1]
      m <- dim(x)[2]
      s <- sum(x)/(n*m)
      xr <- rowSums(x)/m
      xc <- colSums(x)/n
      return((x-outer(xr,xc,"+"))+s)
    }
    #-------------------- end subroutines --------------------
    z <- eigen(-doubleCenter(as.matrix(diss)^2)/2,symmetric=TRUE)
    v <- pmax(z$values,0)
    if (p == 1) normdiag <- cbind(sqrt(v[1])) else normdiag <- diag(sqrt(v[1:p]))
    return(z$vectors[,1:p]%*%normdiag)
  }
`bmat` <-
  function(diss, wgths, d, eps = 1E-12)
  {
    z <- ifelse(d < eps, 1, 0)
    b <- as.matrix((wgths*diss*(1-z))/(d+z))
    r <- rowSums(b) 
    return(diag(r)-b)
  }


