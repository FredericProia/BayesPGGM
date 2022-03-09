library(GeneralizedHyperbolic)
library(mvtnorm)
library(control)
library(ks)
library(expm)

# Gibbs samplers for the Bayesian PGGMs described in the paper
# "A Bayesian approach for partial Gaussian graphical models with sparsity" (E. Okome Obiang, P. Jézéquel, F. Proïa), 2022.
GibbsPGGM = function(Y, X, ListHyp, type = "s", NGrp = c(), shr = "a", std = FALSE, initl = 1, Or = NA, N = 10000, Nd = 5000){
  
  # Arguments :
  # Y : matrix of (n x q) outcomes
  # X : matrix of (n x p) predictors
  # ListHyp : list of hyperparameters (a, b) or (a1, b1, a2, b2) for the prior spike probability
  # type : "s" (Sparse, default), "gs" (Group-Sparse), "sgs" (Sparse-Group-Sparse), "ns" (No Sparsity)
  # NGrp : vector of group sizes (empty if no group structure, default)
  # shr : "g" (global shrinkage) or "a" (adaptative shrinkage, default)
  # std : boolean, standardization of X and Y before running (false, default)
  # initl : initial value of the second parameter in the Gamma prior for the shrinkage parameter (1, default)
  # Or : list with Oy and Lam (and Nu for (sgs)) to run the oracle mode, else NA (default) 
  # N : number of iterations of the sampler (10000, default)
  # Nd : number of iterations of the sampler used for estimation (5000, default), others are burn-ins
  
  # Will return :
  # EstD : estimation of Delta based on posterior median
  # EstOy : estimation of Omega_y based on posterior mean
  # ProbInc : for (gs), probability of inclusion of each group (NA for (s) and (sgs))
  
  p = ncol(X) # Number of predictors
  q = ncol(Y) # Number of outcomes
  n = nrow(X) # Number of observations
  
  # No standardization for oracle mode
  std = ifelse(!is.na(Or[1]), FALSE, std)
  
  # Standardization
  if (std){
    normx = sqrt(colSums(X^2)/n)
    X = sweep(X, 2L, normx, "/", check.margin = FALSE)
    normy = sqrt(colSums(Y^2)/n)
    Y = sweep(Y, 2L, normy, "/", check.margin = FALSE)
  }
  
  # Each predictor is given a group for (gs) and (sgs), according to the sizes
  if ((type == "gs") | (type == "sgs")){
    m = length(NGrp) # Number of groups
    IndGrp = cumsum(c(1, NGrp))
    GrpIndiv = rep(0, p)
    for (g in 1:m){
      GrpIndiv[IndGrp[g]:(IndGrp[g+1]-1)] = g
    }
  }
  
  # No sparsity is equivalent to sparsity with pi=0
  ns = FALSE
  if (type == "ns"){
    type = "s"
    ns = TRUE
  }
  
  # Sparse estimation (s) 
  if (type == "s"){
    
    # Univariate outcomes (s)
    if (q == 1){
      
      # Hyperparameters
      l = rep(initl, p)
      u = 2
      v = 1/2
      a = ListHyp$a
      b = ListHyp$b
      
      # Preliminary computations
      XX = t(X)%*%X
      YY = norm(Y, "2")^2
      YX = t(X)%*%Y

      # Initial values (arbitrary or fixed in oracle mode)
      if (is.na(Or)){
        Oy = 1
        Lam = rep(0, p)
        for (i in 1:p){
          Lam[i] = rgamma(1, shape=1, rate=l[i])
        }
      } else{
        Oy = Or$Oy
        Lam = Or$Lam
      }
      Pi = 0 # If no sparsity
      if (!ns){
        Pi = rbeta(1, a, b)
      }
      D = rep(0, p)
      SLam = Lam
      
      # The N simulations of Delta and Omega_y will be kept
      MEstD = matrix(0, N, p)
      MEstOy = rep(0, N)
      
      # Gibbs sampler
      for (k in 1:N){
        
        # Simulation of Delta and Lambda (except in oracle mode)
        N0 = 0
        for (i in 1:p){
          H = Oy*YX[i] + t(D)%*%XX[, i] - D[i]*XX[i, i]
          S = Lam[i]/(1 + Lam[i]*XX[i, i])
          pr = Pi/(Pi + (1-Pi)*(1+Lam[i]*XX[i, i])^(-1/2)*exp(H^2*S/(2*Oy)))
          if (is.nan(pr)){
            pr = 0
          }
          
          # U(0,1) in [O, pr] : exactly zero
          # Else : non-zero, slab distribution
          if (runif(1) < pr){
            N0 = N0+1
            D[i] = 0
            if (is.na(Or)){
              Lam[i] = rexp(1, l[i])
            }
          } else{
            D[i] = rnorm(1, -H*S, sqrt(S*Oy))
            if (is.na(Or)){
              Lam[i] = rgig(1, chi=D[i]^2/Oy, psi=2*l[i], lambda=1/2)
            }
          }
        }
        
        # Sum of Lam will be useful to update the shrinkage hyperparameter
        SLam = SLam+Lam
        
        # Simulation of Omega_y (except in oracle mode)
        if (is.na(Or)){
          if (N0 < p){
            Oy = rgig(1, chi=t(D)%*%(XX + diag(1/Lam))%*%D, psi=YY+1/v, lambda=(n-p+N0+u)/2)
          } else{
            Oy = rgamma(1, (n+u)/2, rate=(YY+1/v)/2)
          }
        }
          
        # Simulation of Pi (except if no sparsity)
        if (!ns){
          Pi = rbeta(1, N0+a, p-N0+b)
        }
        
        # EM step to update the shrinkage hyperparameter
        if (shr == "a"){ # With adaptative shrinkage
          l = 1/(SLam/k)
        } else if (shr == "g"){ # With global shrinkage
          l = rep(p/sum(SLam/k), p)
        }
    
        # Registration of the k-th generation of the sampler
        MEstD[k, ] = D
        MEstOy[k] = Oy
      }
      
      # Only the Nd last iterations are used to estimate Delta and Omega_y
      MEstD = MEstD[(N-Nd+1):N, ]
      MEstOy = MEstOy[(N-Nd+1):N]
      
      # Posterior median for Delta and posterior mean for Omega_y
      PostMedD = apply(MEstD, 2, median)
      PostMeanOy = mean(MEstOy)
      
      # Denormalization if necessary
      if (std){
        PostMeanOy = PostMeanOy/normy^2
        PostMedD = PostMedD/(normx*normy)
      }
    }
    
    # Multivariate outcomes (s)
    else if (q > 1){
      
      # Hyperparameters
      l = rep(initl, p)
      u = q
      V = diag(q)/q
      a = ListHyp$a
      b = ListHyp$b
      
      # Preliminary computations
      XX = t(X)%*%X
      YY = t(Y)%*%Y
      YX = t(Y)%*%X
      SPhi = sqrtm(YY + solve(V))
      
      # Initial values (arbitrary or fixed in oracle mode)
      if (is.na(Or)){
        Oy = diag(q)
        Lam = rep(0, p)
        for (i in 1:p){
          Lam[i] = rgamma(1, shape=(q+1)/2, rate=l[i])
        }
      } else{
        Oy = Or$Oy
        Lam = Or$Lam
      }
      iOy = solve(Oy)
      Pi = 0 # If no sparsity
      if (!ns){
        Pi = rbeta(1, a, b)
      }
      D = matrix(0, nrow=q, ncol=p)
      SLam = Lam
  
      # The (vectorized) N simulations of Delta and Omega_y will be kept
      MEstD = matrix(0, N, q*p)
      MEstOy = matrix(0, N, q^2)
      
      # Gibbs sampler
      for (k in 1:N){
        
        # Simulation of Delta and Lambda (except in oracle mode)
        N0 = 0
        for (i in 1:p){
          H = Oy%*%YX[, i] + D%*%XX[, i] - D[, i]*XX[i, i]
          S = Lam[i]/(1 + Lam[i]*XX[i, i])
          pr = Pi/(Pi + (1-Pi)*(1+Lam[i]*XX[i, i])^(-q/2)*exp(S*t(H)%*%iOy%*%H/2))
          if (is.nan(pr)){
            pr = 0
          }

          # U(0,1) in [O, pr] : exactly zero
          # Else : non-zero, slab distribution
          if (runif(1) < pr){
            N0 = N0+1
            D[, i] = 0
            if (is.na(Or)){
              Lam[i] = rgamma(1, (q+1)/2, rate=l[i])
            }
          } else{
            D[, i] = rmvnorm(1, -S*H, sigma = S*Oy)
            if (is.na(Or)){
              Lam[i] = rgig(1, chi=t(D[, i])%*%iOy%*%D[, i], psi=2*l[i], lambda=1/2)
            }
          }
        }
        
        # Sum of Lam will be useful to update the shrinkage hyperparameter
        SLam = SLam+Lam
        
        # Simulation of Omega_y (except in oracle mode)
        # Warning : the mode of the posterior MGIG distribution is used (waiting for a better alternative)
        # It is computed by solving a Riccati equation
        if (is.na(Or)){
          A = (n-p+N0+u-q-1)/2*diag(q)
          Psi = D%*%(XX + diag(1/Lam))%*%t(D)
          Oy = care(A, SPhi, Psi)$X
          iOy = solve(Oy)
        }
        
        # Simulation of Pi (except if no sparsity)
        if (!ns){
          Pi = rbeta(1, N0+a, p-N0+b)
        }
        
        # EM step to update the shrinkage hyperparameter
        if (shr == "a"){ # With adaptative shrinkage
          l = ((q+1)/2)/(SLam/k)
        } else if (shr == "g"){ # With global shrinkage
          l = rep((p*(q+1)/2)/sum(SLam/k), p)
        }
        
        # Registration of the k-th generation of the sampler
        MEstD[k, ] = vec(D)
        MEstOy[k, ] = vec(Oy)
      }
      
      # Only the Nd last iterations are used to estimate Delta and Omega_y
      MEstD = MEstD[(N-Nd+1):N, ]
      MEstOy = MEstOy[(N-Nd+1):N, ]
      
      # Posterior median for Delta and posterior mean for Omega_y
      PostMedD = invvec(apply(MEstD, 2, median), p, q)
      PostMeanOy = invvec(apply(MEstOy, 2, mean), q, q)
      
      # Denormalization if necessary
      if (std){
        PostMeanOy = solve(diag(normy))%*%PostMeanOy%*%solve(diag(normy))
        PostMedD = solve(diag(normy))%*%PostMedD%*%solve(diag(normx))
      }
    }
    
    # Probability of inclusion only for (gs)
    ProbInc = NA
  }
  
  # Group-sparse estimation (gs) 
  else if (type == "gs"){
    
    # Univariate outcomes (gs)
    if (q == 1){
      
      # Hyperparameters
      l = rep(initl, m)
      u = 2
      v = 1/2
      a = ListHyp$a
      b = ListHyp$b
      
      # Preliminary computations
      XX = t(X)%*%X
      YY = norm(Y, "2")^2
      YX = t(X)%*%Y
      
      # Initial values (arbitrary or fixed in oracle mode)
      if (is.na(Or)){
        Oy = 1
        Lam = rep(0, m)
        for (g in 1:m){
          Lam[g] = rgamma(1, shape=(NGrp[g]+1)/2, rate=l[g])
        }
      } else{
        Oy = Or$Oy
        Lam = Or$Lam
      }
      Pi = rbeta(1, a, b)
      D = rep(0, p)
      SLam = Lam

      # The N simulations of Delta and Omega_y will be kept
      MEstD = matrix(0, N, p)
      MEstOy = rep(0, N)
      
      # Probability of inclusion of each group
      ProbInc = rep(0, m)
      
      # Gibbs sampler
      for (k in 1:N){
        
        # Simulation of Delta and Lambda (except in oracle mode)
        N0 = 0
        G0 = 0
        DLam = c()
        for (g in 1:m){
          Indg = IndGrp[g]:(IndGrp[g+1]-1)
          H = Oy*YX[Indg] + XX[Indg, ]%*%D - XX[Indg, Indg]%*%D[Indg]
          S = Lam[g]*solve(diag(NGrp[g]) + Lam[g]*XX[Indg, Indg])
          pr = Pi/(Pi + (1-Pi)*det(diag(NGrp[g]) + Lam[g]*XX[Indg, Indg])^(-1/2)*exp(t(H)%*%S%*%H/(2*Oy)))
          if (is.nan(pr)){
            pr = 0
          }
          
          # Update of the probability of inclusion of this group
          if (k > N-Nd){
            ProbInc[g] = ProbInc[g]+pr
          }

          # U(0,1) in [O, pr] : exactly zero
          # Else : non-zero, slab distribution
          if (runif(1) < pr){
            N0 = N0+NGrp[g]
            G0 = G0+1
            D[Indg] = 0
            if (is.na(Or)){
              Lam[g] = rgamma(1, shape=(NGrp[g]+1)/2, rate=l[g])
            }
          } else{
            D[IndGrp[g]:(IndGrp[g+1]-1)] = rmvnorm(1, -S%*%H, sigma=Oy*S)
            if (is.na(Or)){
              Lam[g] = rgig(1, chi=norm(D[IndGrp[g]:(IndGrp[g+1]-1)], "2")^2/Oy, psi=2*l[g], lambda=1/2)
            }
          }
          DLam = c(DLam, rep(Lam[g], NGrp[g]))
        }
        
        # Sum of Lam will be useful to update the shrinkage hyperparameter
        SLam = SLam+Lam
        
        # Simulation of Omega_y (except in oracle mode)
        if (is.na(Or)){
          if (N0 < p){
            Oy = rgig(1, chi=t(D)%*%(XX + diag(1/DLam))%*%D, psi=YY+1/v, lambda=(n-p+N0+u)/2)
          } else{
            Oy = rgamma(1, (n+u)/2, rate=(YY+1/v)/2)
          }
        }
        
        # Simulation of Pi
        Pi = rbeta(1, G0+a, m-G0+b)
        
        # EM step to update the shrinkage hyperparameter
        if (shr == "a"){ # With adaptative shrinkage
          l = ((NGrp+1)/2)/(SLam/k)
        } else if (shr == "g"){ # With global shrinkage
          l = rep(((p+m)/2)/sum(SLam/k), m)
        }

        # Registration of the k-th generation of the sampler
        MEstD[k, ] = D
        MEstOy[k] = Oy
      }
      
      # Only the Nd last iterations are used to estimate Delta and Omega_y
      MEstD = MEstD[(N-Nd+1):N, ]
      MEstOy = MEstOy[(N-Nd+1):N]
      
      # Posterior median for Delta and posterior mean for Omega_y
      PostMedD = apply(MEstD, 2, median)
      PostMeanOy = mean(MEstOy)
      
      # Probability of inclusion of each group
      ProbInc = ProbInc/Nd
      
      # Denormalization if necessary
      if (std){
        PostMeanOy = PostMeanOy/normy^2
        PostMedD = PostMedD/(normx*normy)
      }
    }
    
    # Multivariate outcomes (gs)
    else if (q > 1){

      # Hyperparameters
      l = rep(initl, m)
      u = q
      V = diag(q)/q
      a = ListHyp$a
      b = ListHyp$b
      
      # Preliminary computations
      XX = t(X)%*%X
      YY = t(Y)%*%Y
      YX = t(Y)%*%X
      SPhi = sqrtm(YY + solve(V))
      
      # Initial values (arbitrary or fixed in oracle mode)
      if (is.na(Or)){
        Oy = diag(q)
        Lam = rep(0, m)
        for (g in 1:m){
          Lam[g] = rgamma(1, shape=(q*NGrp[g]+1)/2, rate=l[g])
        }
      } else{
        Oy = Or$Oy
        Lam = Or$Lam
      }
      iOy = solve(Oy)
      Pi = rbeta(1, a, b)
      D = matrix(0, nrow=q, ncol=p)
      SLam = Lam

      # The (vectorized) N simulations of Delta and Omega_y will be kept
      MEstD = matrix(0, N, q*p)
      MEstOy = matrix(0, N, q^2)
      
      # Probability of inclusion of each group
      ProbInc = rep(0, m)
      
      # Gibbs sampler
      for (k in 1:N){
        
        # Simulation of Delta and Lambda (except in oracle mode)
        N0 = 0
        G0 = 0
        DLam = c()
        for (g in 1:m){
          Indg = IndGrp[g]:(IndGrp[g+1]-1)
          H = Oy%*%YX[, Indg] + D%*%XX[, Indg] - D[, Indg]%*%XX[Indg, Indg]
          S = Lam[g]*solve(diag(NGrp[g]) + Lam[g]*XX[Indg, Indg])
          pr = Pi/(Pi + (1-Pi)*det(diag(NGrp[g]) + Lam[g]*XX[Indg, Indg])^(-q/2)*exp(sum(diag(t(H)%*%iOy%*%H%*%S/2))))
          if (is.nan(pr)){
            pr = 0
          }
          
          # Update of the probability of inclusion of this group
          if (k > N-Nd){
            ProbInc[g] = ProbInc[g]+pr
          }
          
          # U(0,1) in [O, pr] : exactly zero
          # Else : non-zero, slab distribution
          if (runif(1) < pr){
            N0 = N0+NGrp[g]
            G0 = G0+1
            D[, Indg] = 0
            if (is.na(Or)){
              Lam[g] = rgamma(1, shape=(q*NGrp[g]+1)/2, rate=l[g])
            }
          } else{
            VD = rmvnorm(1, vec(-H%*%S), sigma=S%x%Oy)
            D[, Indg] = invvec(VD, NGrp[g], q)
            if (is.na(Or)){
              Lam[g] = rgig(1, chi=sum(diag(t(D[, Indg])%*%iOy%*%D[, Indg])), psi=2*l[g], lambda=1/2)
            }
          }
          DLam = c(DLam, rep(Lam[g], NGrp[g]))
        }
        
        # Sum of Lam will be useful to update the shrinkage hyperparameter
        SLam = SLam+Lam
        
        # Simulation of Omega_y (except in oracle mode)
        # Warning : the mode of the posterior MGIG distribution is used (waiting for a better alternative)
        # It is computed by solving a Riccati equation
        if (is.na(Or)){
          A = (n-p+N0+u-q-1)/2*diag(q)
          Psi = D%*%(XX + diag(1/DLam))%*%t(D)
          Oy = care(A, SPhi, Psi)$X
          iOy = solve(Oy)
        }
        
        # Simulation of Pi
        Pi = rbeta(1, G0+a, m-G0+b)
        
        # EM step to update the shrinkage hyperparameter
        if (shr == "a"){ # With adaptative shrinkage
          l = ((q*NGrp+1)/2)/(SLam/k)
        } else if (shr == "g"){ # With global shrinkage
          l = rep(((q*p+m)/2)/sum(SLam/k), m)
        }
        
        # Registration of the k-th generation of the sampler
        MEstD[k, ] = vec(D)
        MEstOy[k, ] = vec(Oy)
      }
      
      # Only the Nd last iterations are used to estimate Delta and Omega_y
      MEstD = MEstD[(N-Nd+1):N, ]
      MEstOy = MEstOy[(N-Nd+1):N, ]
      
      # Posterior median for Delta and posterior mean for Omega_y
      PostMedD = invvec(apply(MEstD, 2, median), p, q)
      PostMeanOy = invvec(apply(MEstOy, 2, mean), q, q)
      
      # Probability of inclusion of each group
      ProbInc = ProbInc/Nd
      
      # Denormalization if necessary
      if (std){
        PostMeanOy = solve(diag(normy))%*%PostMeanOy%*%solve(diag(normy))
        PostMedD = solve(diag(normy))%*%PostMedD%*%solve(diag(normx))
      }
    }
    
  }
  
  # Sparse-group-sparse estimation (sgs) 
  else if (type == "sgs"){
    
    # Univariate outcomes (sgs)
    if (q == 1){
      
      # Hyperparameters
      gam = rep(initl, m)
      l = rep(initl, p)
      u = 2
      v = 1/2
      a1 = ListHyp$a1
      a2 = ListHyp$a2
      b1 = ListHyp$b1
      b2 = ListHyp$b2
      
      # Preliminary computations
      XX = t(X)%*%X
      YY = norm(Y, "2")^2
      YX = t(X)%*%Y
      
      # Initial values (arbitrary or fixed in oracle mode)
      if (is.na(Or)){
        Oy = 1
        if (shr == "a"){
          Lam = rep(1, m) # Lam=1 in adaptative shrinkage
          Nu = rep(0, p)
          for (i in 1:p){
            Nu[i] = rgamma(1, shape=1, rate=l[i])
          }
        } else if (shr == "g"){
          Nu = rep(1, p) # Nu=1 in global shrinkage
          Lam = rep(0, m)
          for (g in 1:m){
            Lam[g] = rgamma(1, shape=(NGrp[g]+1)/2, rate=gam[g])
          }
        }
      } else{
        Oy = Or$Oy
        if (shr == "a"){
          Lam = rep(1, m)
          Nu = Or$Nu
        } else if (shr == "g"){
          Lam = Or$Lam
          Nu = rep(1, p)
        }
      }
      Pi1 = rbeta(1, a1, b1)
      Pi2 = rbeta(1, a2, b2)
      D = rep(0, p)
      Dln = rep(p, 0)
      for (g in 1:m){
        Dln[IndGrp[g]:(IndGrp[g+1]-1)] = Lam[g]*Nu[IndGrp[g]:(IndGrp[g+1]-1)]
      }
      SLam = Lam
      SNu = Nu
      
      # The N simulations of Delta and Omega_y will be kept
      MEstD = matrix(0, N, p)
      MEstOy = rep(0, N)
      
      # Gibbs sampler
      for (k in 1:N){
        
        # Simulation of Delta, Lambda and Nu (except in oracle mode)
        N0 = 0
        G0 = 0
        J0 = 0
        ag = 1
        Igi = rep(0, p)
        for (i in 1:p){
          g = GrpIndiv[i]
          H = Oy*YX[i] + t(D)%*%XX[, i] - D[i]*XX[i, i]
          S = Nu[i]*Lam[g]/(1 + Nu[i]*Lam[g]*XX[i, i])
          r = (1-Pi1)*Pi2*Igi[i] + Pi1*(1-Igi[i])
          pr = r/(r + (1-Pi1)*(1-Pi2)*(1+Nu[i]*Lam[g]*XX[i, i])^(-1/2)*exp(H^2*S/(2*Oy)))
          if (is.nan(pr)){
            pr = 0
          }
          
          # U(0,1) in [O, pr] : exactly zero
          # Else : non-zero, slab distribution
          # Nu updated only for adaptative shrinkage (except in oracle mode)
          if (runif(1) < pr){
            N0 = N0+1
            D[i] = 0
            if ( (shr == "a") & (is.na(Or)) ){
              Nu[i] = rexp(1, l[i])
            }
          } else{
            D[i] = rnorm(1, -H*S, sqrt(S*Oy))
            if ( (shr == "a") & (is.na(Or)) ){
              Nu[i] = rgig(1, chi=D[i]^2/(Oy*Lam[g]), psi=2*l[i], lambda=1/2)
            }
          }
        
          # Treatments for the group just explored
          # Lam updated only for global shrinkage (except in oracle mode)
          if (i == p){
            g = g+1
          }
          if (g != ag){
            Indg = IndGrp[g-1]:(IndGrp[g]-1)
            for (j in Indg){
              Indgj = Indg[Indg != j]
              Igi[j] = (sum(D[Indgj] != 0) > 0)
            }
            if (sum(D[Indg] != 0) > 0){
              N0g = sum(D[Indg] == 0)
              J0 = J0+N0g
              Nug = Nu[Indg]
              if ( (shr == "g") & (is.na(Or)) ){
                Lam[g-1] = rgig(1, chi=t(D[Indg])%*%diag(1/Nug)%*%D[Indg]/Oy, psi=2*gam[g-1], lambda=(NGrp[g-1]+N0g)/2)
              }
            } else{
              G0 = G0+1
              if ( (shr == "g") & (is.na(Or)) ){
                Lam[g-1] = rgamma(1, shape=NGrp[g-1], rate=gam[g-1])
              }
            }
            Dln[Indg] = Lam[g-1]*Nu[Indg]
            ag = g
          }
        }
        
        # Sums of Lam and Nu will be useful to update the shrinkage hyperparameters
        SLam = SLam+Lam
        SNu = SNu+Nu
        
        # Simulation of Omega_y (except in oracle mode)
        if (is.na(Or)){
          if (N0 < p){
            Oy = rgig(1, chi=t(D)%*%(XX + diag(1/Dln))%*%D, psi=YY+1/v, lambda=(n-p+N0+u)/2)
          } else{
            Oy = rgamma(1, (n+u)/2, rate=(YY+1/v)/2)
          }
        }
        
        # Simulation of Pi1 and Pi2
        Pi1 = rbeta(1, G0+a1, m-G0+b1)
        Pi2 = rbeta(1, J0+a2, p-N0+b2)
        
        # EM step to update the shrinkage hyperparameters
        if (shr == "a"){ # With adaptative shrinkage
          l = 1/(SNu/k)
        } else if (shr == "g"){ # With global shrinkage
          gam = ((NGrp+1)/2)/(SLam/k)
        }
        
        # Registration of the k-th generation of the sampler
        MEstD[k, ] = D
        MEstOy[k] = Oy
      }
      
      # Only the Nd last iterations are used to estimate Delta and Omega_y
      MEstD = MEstD[(N-Nd+1):N, ]
      MEstOy = MEstOy[(N-Nd+1):N]
      
      # Posterior median for Delta and posterior mean for Omega_y
      PostMedD = apply(MEstD, 2, median)
      PostMeanOy = mean(MEstOy)
      
      # Denormalization if necessary
      if (std){
        PostMeanOy = PostMeanOy/normy^2
        PostMedD = PostMedD/(normx*normy)
      }
    }
    
    # Multivariate outcomes (sgs)
    else if (q > 1){
      
      # Hyperparameters
      gam = rep(initl, m)
      l = rep(initl, p)
      u = q
      V = diag(q)/q
      a1 = ListHyp$a1
      a2 = ListHyp$a2
      b1 = ListHyp$b1
      b2 = ListHyp$b2
      
      # Preliminary computations
      XX = t(X)%*%X
      YY = t(Y)%*%Y
      YX = t(Y)%*%X
      SPhi = sqrtm(YY + solve(V))
      
      # Initial values (arbitrary or fixed in oracle mode)
      if (is.na(Or)){
        Oy = diag(q)
        if (shr == "a"){
          Lam = rep(1, m) # Lam=1 in adaptative shrinkage
          Nu = rep(0, p)
          for (i in 1:p){
            Nu[i] = rgamma(1, shape=(q+1)/2, rate=l[i])
          }
        } else if (shr == "g"){
          Nu = rep(1, p) # Nu=1 in global shrinkage
          Lam = rep(0, m)
          for (g in 1:m){
            Lam[g] = rgamma(1, shape=(q*NGrp[g]+1)/2, rate=gam[g])
          }
        }
      } else{
        Oy = Or$Oy
        if (shr == "a"){
          Lam = rep(1, m)
          Nu = Or$Nu
        } else if (shr == "g"){
          Lam = Or$Lam
          Nu = rep(1, p)
        }
      }
      iOy = solve(Oy)
      Pi1 = rbeta(1, a1, b1)
      Pi2 = rbeta(1, a2, b2)
      D = matrix(0, nrow=q, ncol=p)
      Dln = rep(p, 0)
      for (g in 1:m){
        Dln[IndGrp[g]:(IndGrp[g+1]-1)] = Lam[g]*Nu[IndGrp[g]:(IndGrp[g+1]-1)]
      }
      SLam = Lam
      SNu = Nu
      
      # The (vectorized) N simulations of Delta and Omega_y will be kept
      MEstD = matrix(0, N, q*p)
      MEstOy = matrix(0, N, q^2)
      
      # Gibbs sampler
      for (k in 1:N){
        
        # Simulation of Delta, Lambda and Nu (except in oracle mode)
        N0 = 0
        G0 = 0
        J0 = 0
        ag = 1
        Igi = rep(0, p)
        for (i in 1:p){
          g = GrpIndiv[i]
          H = Oy%*%YX[, i] + D%*%XX[, i] - D[, i]*XX[i, i]
          S = Nu[i]*Lam[g]/(1 + Nu[i]*Lam[g]*XX[i, i])
          r = (1-Pi1)*Pi2*Igi[i] + Pi1*(1-Igi[i])
          pr = r/(r + (1-Pi1)*(1-Pi2)*(1+Nu[i]*Lam[g]*XX[i, i])^(-q/2)*exp(S*t(H)%*%iOy%*%H/2))
          if (is.nan(pr)){
            pr = 0
          }
          
          # U(0,1) in [O, pr] : exactly zero
          # Else : non-zero, slab distribution
          # Nu updated only for adaptative shrinkage (except in oracle mode)
          if (runif(1) < pr){
            N0 = N0+1
            D[, i] = 0
            if ( (shr == "a") & (is.na(Or)) ){
              Nu[i] = rgamma(1, shape=(q+1)/2, rate=l[i])
            }
          } else{
            D[, i] = rmvnorm(1, -S*H, sigma = S*Oy)
            if ( (shr == "a") & (is.na(Or)) ){
              Nu[i] = rgig(1, chi=t(D[, i])%*%iOy%*%D[, i], psi=2*l[i], lambda=1/2)
            }
          }
          
          # Treatments for the group just explored
          # Lam updated only for global shrinkage (except in oracle mode)
          if (i == p){
            g = g+1
          }
          if (g != ag ){
            Indg = IndGrp[g-1]:(IndGrp[g]-1)
            for (j in Indg){
              Indgj = Indg[Indg != j]
              Igi[j] = (sum(D[, Indgj] != 0) > 0)
            }
            if (sum(D[, Indg] != 0) > 0){
              N0g = sum(colSums(abs(D[, Indg])) == 0)
              J0 = J0+N0g
              Nug = Nu[Indg]
              if ( (shr == "g") & (is.na(Or)) ){
                Lam[g-1] = rgig(1, chi=sum(diag(diag(1/Nug)%*%t(D[, Indg])%*%iOy%*%D[, Indg])), psi=2*gam[g-1], lambda=(q*N0g+1)/2)
              }
            } else{
              G0 = G0+1
              if ( (shr == "g") & (is.na(Or)) ){
                Lam[g-1] = rgamma(1, shape=(q*NGrp[g-1]+1)/2, rate=gam[g-1])
              }
            }
            Dln[Indg] = Lam[g-1]*Nu[Indg]
            ag = g
          }
        }
        
        # Sums of Lam and Nu will be useful to update the shrinkage hyperparameters
        SLam = SLam+Lam
        SNu = SNu+Nu
        
        # Simulation of Omega_y (except in oracle mode)
        # Warning : the mode of the posterior MGIG distribution is used (waiting for a better alternative)
        # It is computed by solving a Riccati equation
        if (is.na(Or)){
          A = (n-p+N0+u-q-1)/2*diag(q)
          Psi = D%*%(XX + diag(1/Dln))%*%t(D)
          Oy = care(A, SPhi, Psi)$X
          iOy = solve(Oy)
        }

        # Simulation of Pi1 and Pi2
        Pi1 = rbeta(1, G0+a1, m-G0+b1)
        Pi2 = rbeta(1, J0+a2, p-N0+b2)
        
        # EM step to update the shrinkage hyperparameters
        if (shr == "a"){ # With adaptative shrinkage
          l = ((q+1)/2)/(SNu/k)
        } else if (shr == "g"){ # With global shrinkage
          gam = ((q*NGrp+1)/2)/(SLam/k)
        }
        
        # Registration of the k-th generation of the sampler
        MEstD[k, ] = vec(D)
        MEstOy[k, ] = vec(Oy)
      }
      
      # Only the Nd last iterations are used to estimate Delta and Omega_y
      MEstD = MEstD[(N-Nd+1):N, ]
      MEstOy = MEstOy[(N-Nd+1):N, ]
      
      # Posterior median for Delta and posterior mean for Omega_y
      PostMedD = invvec(apply(MEstD, 2, median), p, q)
      PostMeanOy = invvec(apply(MEstOy, 2, mean), q, q)
      
      # Denormalization if necessary
      if (std){
        PostMeanOy = solve(diag(normy))%*%PostMeanOy%*%solve(diag(normy))
        PostMedD = solve(diag(normy))%*%PostMedD%*%solve(diag(normx))
      }
    }
    
    # Probability of inclusion only for (gs)
    ProbInc = NA
  }
  
  # Return :
  # EstD : estimation of Delta
  # EstOy : estimation of Omega_y
  # ProbInc : for (gs), probability of inclusion of each group
  return(list(EstD = PostMedD, EstOy = PostMeanOy, ProbInc = ProbInc))
}