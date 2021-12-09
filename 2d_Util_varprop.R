mse <- function(x){
  # mean sqaured error
  mean(x^2)
}

mae <- function(x){
  # mean absolute error
  mean(abs(x))
}

rank_calculation <- function(x, cut_off=1e-1){
  # calculate the rank of the final estimated W
  W <- x$par
  sum(svd(W)$d>cut_off)
}

truncate_A <- function(A, d){
  # approximate matrix A using first d svd
  tmp <- svd(A)
  A_approx <- tmp$u[,1:d]%*%diag(tmp$d[1:d])%*%t(tmp$v[,1:d])
  return(A_approx)
}

randomLambda <- function(R, spectral_norm=0.9){
  # function to simulate a R matrix for the R-rank transition operator A
  Lambda <- matrix(rnorm(R^2), R, R)
  Lambda_svd <- svd(Lambda)
  Lambda <- Lambda_svd$u%*%diag(Lambda_svd$d/max(Lambda_svd$d)*spectral_norm)%*%t(Lambda_svd$v)
  return(Lambda)
}

randomLambda_Aue <- function(R, sigma_l, spectral_norm){
  # function to simulate a R matrix for the R-rank transition operator A based on Aue (2014)
  sd <- sigma_l%*%t(sigma_l)
  Lambda <- matrix(rnorm(R^2, c(sd)), R, R)
  Lambda_svd <- svd(Lambda)
  Lambda <- Lambda_svd$u%*%diag(Lambda_svd$d/max(Lambda_svd$d)*spectral_norm)%*%t(Lambda_svd$v)
  return(Lambda)
}

cosine_fun_2d <- function(s,k=1){
  # cosine function for RKHS basis
  #point = #matrix(s,nrow = )
   #print(s)
  if(k==1){
    return(cos(0*(sum(s))))
  }
  else{
    return(sqrt(2)*cos(2*(k-1)*pi*sum(s)))
  }
}

A_transition <- function(x1,x2,Lambda){
  # function to calculate the true transition operator at A(x1,x2)
  R <- dim(Lambda)[1]
  u_x1 <- u_x2 <- c()
  for(r in 1:R){
    u_x1 <- c(u_x1, cosine_fun_2d(x1,k=r))# need test
    u_x2 <- c(u_x2, cosine_fun_2d(x2,k=r))
  }
    #print(length(u_x1))
    #print(dim(Lambda))
  value <- c(t(u_x1)%*%Lambda%*%t(t(u_x2)))
  return(value)
}

A_transition_vec <- Vectorize(A_transition, vectorize.args=c('x1','x2'))

far_Xt <- function(s, x_coef, R){
  # function to recover the true FAR at time s
  Xt_value <- rep(0, nrow(s))
  for(k in 1:R){
      u.q = c()
      for(row in 1:nrow(s)){
          tmp.u <- cosine_fun_2d(s[row,], k=k)
          u.q = c(u.q,tmp.u)
      } 
    tmp <- u.q*c(x_coef[k])
    Xt_value <- Xt_value + tmp
  }
  return(Xt_value)
}

far_Xt.fast <- function(s, R){
    Umat = c()
  for(k in 1:R){
      tmp = sapply(1:nrow(s),function(row) cosine_fun_2d(s[row,],k = k))
      Umat = cbind(Umat,tmp)             
  }
  return(Umat)
}


kernel <- function(x1,x2,sigma.kernel){
   #sigma = sqrt(0.001)
    #sigma = sqrt(0.1)#0.1
  # function for reproducing kernel
  # return value k(x1,x2) = 1+k1(x1)k1(x2)+k2(x1)*k2(x2)-k4(abs(x1-x2)) [see Gu(2013) in Chapter 2.3.3]
  #val <- 1+(x1-0.5)*(x2-0.5) #1+k1(x1)k1(x2)
  #val <- val+((x1-0.5)^2-1/12)/2*((x2-0.5)^2-1/12)/2 #k2(x1)*k2(x2)
  #wk <- abs(x1-x2)
  #val <- val-((wk-0.5)^4-(wk-0.5)^2/2+7/240)/24
  val = exp(-norm(x2-x1,type = "2")^2/(2*sigma.kernel^2))
  return(val)
}
                   
kernel.poly = function(x1,x2){
  # function for reproducing kernel
  # return value k(x1,x2) = 1+k1(x1)k1(x2)+k2(x1)*k2(x2)-k4(abs(x1-x2)) [see Gu(2013) in Chapter 2.3.3]
  val <- 1+(x1-0.5)*(x2-0.5) #1+k1(x1)k1(x2)
  val <- val+((x1-0.5)^2-1/12)/2*((x2-0.5)^2-1/12)/2 #k2(x1)*k2(x2)
  wk <- abs(x1-x2)
  val <- val-((wk-0.5)^4-(wk-0.5)^2/2+7/240)/24
  return(val)
}                   

penloss_func <- function(W, X, Z, K_sqr, lambda_pen){
  # the penalized loss function f(W)+lambda*nuclearnorm(W)
  n <- dim(X)[1]
  fW <- sum((X-K_sqr%*%W%*%Z/n)^2) # f(W)
  penW <- lambda_pen*sum(abs(svd(W)$d))
  return(fW+penW)
}

grad_f <- function(W, K, gradf_part1, Z2){
  # gradient of the likelihood function f(W)
  return(gradf_part1+K%*%W%*%Z2)
}

Q_func <- function(W, W_prev, X, Z, Z2, K, K_sqr, gradf_part1, lambda_pen, stepL){
  # the approximate penalized loss function in JiYe(2009)
  n <- dim(X)[1]
  sum((X-K_sqr%*%W_prev%*%Z/n)^2)+tr(t(W-W_prev)%*%grad_f(W_prev,K,gradf_part1,Z2)) + stepL/2*sum((W-W_prev)^2) +
    lambda_pen*sum(abs(svd(W)$d))
}

argmin_Qfunc <- function(W_prev, K, lambda_pen, stepL){
  # function to obtain the minimizer of Q func in JiYe(2009)
  C <- W_prev-1/stepL*grad_f(W_prev,K,gradf_part1,Z2)
  thres <- lambda_pen/stepL
  tmp <- svd(C)
  d <- pmax(0, tmp$d-thres)
  d <- diag(d[d>0])
  dimension_d <- dim(d)[1]
  return(tmp$u[,1:dimension_d]%*%d%*%t(tmp$v[,1:dimension_d]))
}

tr <- function(W){
  # trace function
  sum(diag(W))
}

trace_norm_optim <- function(W_init, X, Z, Z2, K, K_sqr, gradf_part1, lambda_pen, stepL0, gamma=1.05, niter=1000, tol=10^(-8)){
  # optimization routine in JiYe(2009) Algorithm 2 and has exactly the same result as trace_norm_optim
  # This is a clear representation of JiYe(2009).
  stepL_prev <- stepL0
  W_prev <- W_init # intital guess of W
  Y <- W_init # Z matrix in JiYe(2009)
  alpha_prev <- 1
  value_F_prev <- Inf
  stop_counter <- 0
  for(iter in 1:niter){
    # print(iter)
    stepL <- stepL_prev
    passed <- F
    while (!passed){
      W <- argmin_Qfunc(Y,K,lambda_pen,stepL)
      value_F <- penloss_func(W,X,Z,K_sqr,lambda_pen)
      value_Q <- Q_func(W,Y,X,Z,Z2,K,K_sqr,gradf_part1,lambda_pen,stepL)
      passed <- value_F<=value_Q # can simplify here
      if(!passed){
        stepL <- stepL*gamma
      }
    }
    stepL_prev <- stepL
    W_current <- W
    alpha_current <- (1+sqrt(1+4*alpha_prev^2))/2
    Y <- W_current+((alpha_prev-1)/alpha_current)*(W_current-W_prev)
    W_prev <- W_current
    alpha_prev <- alpha_current
    if(value_F/value_F_prev>1-tol){
      stop_counter <- stop_counter + 1
    }
    value_F_prev <- value_F
    if(stop_counter>10){
      break
    }
  }
  # summarize optimization result
  result <- list()
  result$par <- W_current
  result$iter <- iter
  result$value <- value_F
  return(result)
}

trace_norm_optim_fast <- function(W_init, X, Z, Z2, K, K_sqr, gradf_part1, lambda_pen, stepL0, gamma=1.01, niter=1e5, tol=10^(-6)){
  # optimization routine in JiYe(2009) Algorithm 2 (with algebraic simplication) [it has the same result as trace_norm_optim]
  stepL_prev <- stepL0
  W_prev <- W_init # intital guess of W
  Y <- W_init # Z matrix in JiYe(2009)
  alpha_prev <- 1
  value_F_prev <- Inf
  stop_counter <- 0 # counter for if function value decreases less than tolerance relatively
  n <- dim(X)[1]
  n_tracks <- 50 # track the optimization process
  stepL_sequence <- value_F_sequence <- stop_counter_sequence <- rank_W_sequence <- c() # record the function values for each iteration
  for(iter in 1:niter){
    stepL <- stepL_prev
    if(iter%%(niter/n_tracks)==1){
      time1 <- proc.time()
    }
    # gradient at Y
    passed <- F
    grad_f_atY <- gradf_part1+K%*%Y%*%Z2 # grad_f(Y,K,gradf_part1,Z2)
    f_atY <- sum((X-K_sqr%*%Y%*%Z/n)^2)
    while (!passed){
      # W <- argmin_Qfunc(Y,K,lambda_pen,stepL)
      C <- Y-1/stepL*grad_f_atY
      thres <- lambda_pen/stepL
      tmp <- svd(C)
      d <- pmax(0, tmp$d-thres)
      positive_d <- d[d>0]
      dimension_d <- length(positive_d)
      if(dimension_d==0){
        W <- diag(rep(0, dim(Y)[1]))
      }else{
        # d <- diag(positive_d, nrow=dimension_d)
        # W <- tmp$u[,1:dimension_d]%*%d%*%t(tmp$v[,1:dimension_d])
        W <- (tmp$u[,1:dimension_d]*rep(positive_d, rep(n, dimension_d)))%*%t(tmp$v[,1:dimension_d])
      }
      
      # value_F <- penloss_func(W,X,Z,K_sqr,lambda_pen)
      fW <- sum((X-K_sqr%*%W%*%Z/n)^2) # f(W) computationally consuming
      penW <- lambda_pen*sum(positive_d) # penW <- lambda_pen*sum(abs(svd(W)$d)) nuclear norm of W
      value_F <- fW+penW
      
      # value_Q <- Q_func(W,Y,X,Z,Z2,K,K_sqr,gradf_part1,lambda_pen,stepL)
      diff_WY <- W-Y
      value_Q <- f_atY + sum(diag(t(diff_WY)%*%grad_f_atY)) + stepL/2*sum((diff_WY)^2) + penW
      
      passed <- value_F<=value_Q
      if(!passed){
        stepL <- stepL*gamma
      }
    }
    stepL_prev <- stepL
    W_current <- W
    alpha_current <- (1+sqrt(1+4*alpha_prev^2))/2
    Y <- W_current+((alpha_prev-1)/alpha_current)*(W_current-W_prev)
    W_prev <- W_current
    alpha_prev <- alpha_current
    # if(value_F-value_F_prev > -tol){ # absolute error 
    if(value_F > value_F_prev*(1-tol)){ # relative error
      stop_counter <- stop_counter + 1
    }
    value_F_prev <- value_F
    
    # trace the optimization process
    if(iter%%(niter/n_tracks)==0){
      value_F_sequence <- c(value_F_sequence, value_F_prev)
      stepL_sequence <- c(stepL_sequence, stepL)
      stop_counter_sequence <- c(stop_counter_sequence, stop_counter)
      rank_W_sequence <- c(rank_W_sequence, dimension_d)
      #print(iter)
      #print(value_F)
      #print(stepL)
      #print(stop_counter)
      #print(dimension_d)
      time2 <- proc.time()
      print(as.numeric(time2-time1)[1])
    }
    if(stop_counter>10){
      break
    }
  }
  # summarize optimization result
  result <- list()
  result$par <- W_current
  result$iter <- iter
  result$value <- value_F
  result$value_F_sequence <- value_F_sequence
  result$stepL_sequence <- stepL_sequence
  result$stop_counter_sequence <- stop_counter_sequence
  result$rank_W_sequence <- rank_W_sequence
  return(result)
}

cv_rkhs <- function(cv_folds=5, X_all, Xl_all, K, K_sqr, lambda_pen, stepL0, gamma=1.05, niter=1000, tol=10^(-6), A_true, A_fro_est, A_var_est,n,init_zero){
  # function for estimating forecating error via cross validation for a given lambda_pen  
  # X_all: observation matrix
  # Xl_all: lagged observation matrix
  total_train <- dim(X_all)[2]
  cv_indices <- createFolds(1:total_train, k=cv_folds)
  
  cv_err <- c()
  for(cv_fold_index in 1:cv_folds){
    cv_index_test <- sort(cv_indices[[cv_fold_index]])
    #print(cv_index_test)
    cv_index_train <- sort(unlist(cv_indices)[!unlist(cv_indices)%in%cv_index_test])
    #print(cv_index_train)
    # CV Training
    X <- X_all[,cv_index_train]
    Xl <- Xl_all[,cv_index_train]
    Z <- K_sqr%*%Xl
    gradf_part1 <- -2/n*K_sqr%*%X%*%t(Z)
    Z2 <- 2/n^2*Z%*%t(Z)
    
    # Optimization of the penalized loss function (we try three different initial points)
    #W_init1 <- solve(K_sqr)%*%A_true%*%solve(K_sqr) # initial using true A
    #opt_result1 <- trace_norm_optim_fast(W_init=W_init1, X=X, Z=Z, Z2=Z2, K=K, K_sqr=K_sqr, gradf_part1=gradf_part1,
                                         #lambda_pen=lambda_pen, stepL0=stepL0, gamma=gamma_rkhs, niter=niter_rkhs, tol=tolerance_rkhs)
    
    #W_init2 <- solve(K_sqr)%*%A_var_est%*%solve(K_sqr) # initial using Aue(2014) estimation
    #opt_result2 <- trace_norm_optim_fast(W_init=W_init2, X=X, Z=Z, Z2=Z2, K=K, K_sqr=K_sqr, gradf_part1=gradf_part1,
                                         #lambda_pen=lambda_pen, stepL0=stepL0, gamma=gamma_rkhs, niter=niter_rkhs, tol=tolerance_rkhs)
    
    # Summarize optimization result
    #opt_result <- list(opt_result1, opt_result2) # combine all optimization results
    #fvalues <- c(opt_result1$value, opt_result2$value)
    #niters <- c(opt_result1$iter, opt_result2$iter)
      
      if(init_zero == TRUE){
      W_init <- solve(K_sqr)%*%matrix(0,n,n)%*%solve(K_sqr) # initial using Aue(2014) estimation
    }else if(init_zero == "Fro"){
      W_init <- solve(K_sqr)%*%A_fro_est%*%solve(K_sqr) # initial using Fro estimation
    }
    else if (init_zero == "aue"){
        W_init <- solve(K_sqr)%*%A_var_est%*%solve(K_sqr) # initial using Fro estimationv 
    }
    opt_result <- trace_norm_optim_fast(W_init=W_init, X=X, Z=Z, Z2=Z2, K=K, K_sqr=K_sqr, gradf_part1=gradf_part1,
                                        lambda_pen=lambda_pen, stepL0=stepL0, gamma=gamma_rkhs, niter=niter_rkhs,
                                        tol=tolerance_rkhs) 
    
    W_est <- opt_result$par#opt_result[[which.min(fvalues)]]$par
    D_est <- solve(K_sqr)%*%W_est%*%solve(K_sqr)
    A_est <- K%*%D_est%*%K
    
    # one step ahead prediction
    X_true <- X_all[,cv_index_test]
    X_forecast_rkhs <- A_est%*%Xl_all[,cv_index_test]/n ##A_est is a n by n matrix, nothing to do with T
    cv_err <- c(cv_err, mean(apply(X_true-X_forecast_rkhs, 2, mse)))
  }
  return(cv_err)
}

cv_rkhs_Rcpp <- function(cv_folds=5, X_all, Xl_all, K, K_sqr, lambda_pen, stepL0, gamma=1.05, niter=1000, tol=10^(-6), A_true, A_var_est,init_zero=F){
  # function for estimating forecating error via cross validation for a given lambda_pen
  # X_all: observation matrix
  # Xl_all: lagged observation matrix
  n <- dim(X_all)[1]
  total_train <- dim(X_all)[2]
  cv_indices <- createFolds(1:total_train, k=cv_folds)
  
  cv_err <- c()
  for(cv_fold_index in 1:cv_folds){
    print(paste('CV fold', cv_fold_index))
    cv_index_test <- sort(cv_indices[[cv_fold_index]])
    cv_index_train <- sort(unlist(cv_indices)[!unlist(cv_indices)%in%cv_index_test])
    # CV Training
    X <- X_all[,cv_index_train]
    Xl <- Xl_all[,cv_index_train]
    Z <- K_sqr%*%Xl
    gradf_part1 <- -2/n*K_sqr%*%X%*%t(Z)
    Z2 <- 2/n^2*Z%*%t(Z)
    
    # Optimization of the penalized loss function (we try three different initial points)
    if(init_zero){
      W_init <- solve(K_sqr)%*%matrix(0,n,n)%*%solve(K_sqr) # initial using Aue(2014) estimation
    }else{
      W_init <- solve(K_sqr)%*%A_true%*%solve(K_sqr) # initial using Aue(2014) estimation
    }
    opt_result <- trace_norm_optim_Rcpp(W_init=W_init, X=X, Z=Z, Z2=Z2, K=K, K_sqr=K_sqr, gradf_part1=gradf_part1,
                                        lambda_pen=lambda_pen, stepL0=stepL0, gamma=gamma_rkhs, niter=niter_rkhs,
                                        tol=tolerance_rkhs, n=n)
    
    # Summarize optimization result
    W_est <- opt_result$par
    D_est <- solve(K_sqr)%*%W_est%*%solve(K_sqr)
    A_est <- K%*%D_est%*%K
    
    # one step ahead prediction
    X_true <- X_all[,cv_index_test]
    X_forecast_rkhs <- A_est%*%Xl_all[,cv_index_test]/n
    cv_err <- c(cv_err, mean(apply(X_true-X_forecast_rkhs, 2, mse)))
  }
  return(cv_err)
}
                   
                   
far2.generate = function(totalT_train,n,R,Lambda,sample_points,noise_range){
  burnin <- 100 # burnin stage for simulation of FAR(1)
  #noise_range <- 0.1 # functional noise, Z is U(-noise_range, noise_range)
  sampling_noise <- 0
  totalT_train <- totalT_train # number of curves
  totalT_test <- round(totalT_train*0.2)
  totalT <- totalT_train + totalT_test
  #n <- current_simulation_setting[2] # number of sampling points per curve
  #R <- current_simulation_setting[3] # rank of transition operator A
  #spectral_norm <- current_simulation_setting[4] # the spectral norm of the diagonal transition operator A
  init_zero <- T # whether use 0 matrix as initial value of optimization
  #Lambda <- diag(rep(spectral_norm, R))
  # Simulate sampling points and transition operator A
    #sample_points <- seq(0,1,length.out=n) # s1,s2,...,sn
    # sample_points <- sort(runif(n,0,1))
    
    # Plot contour of true transition operator A at (s1,s2,...,sn) by (s1,s2,...,sn)
    A_true <- sapply(1:n,function(j) sapply(1:n,function(i) A_transition(sample_points[i,],sample_points[j,],Lambda)))
                                            #outer(sample_points, sample_points, A_transition_vec, Lambda=Lambda)
  far_coef <- c() # xt coefficient
    for(t in 1:(totalT+burnin)){
      z <- runif(R, min=-noise_range, max=noise_range)
      if(t==1){
        far_coef <- rbind(far_coef, z)
      }else{
        tmp <- Lambda%*%far_coef[(t-1),]+z
          #print(tmp)
        far_coef <- rbind(far_coef, c(tmp))
      }
    }
    far_coef <- far_coef[-c(1:burnin),]
    rownames(far_coef) <- NULL
    far_coef_onestep <- t(Lambda%*%t(far_coef)) # best one step ahead estimation of far_coef
    #print(dim(far_coef))
    # Sample FAR time series at (s1,s2,...,sn)
    #far_sample <- c()
    #for(t in 1:totalT){
      #tmp <- far_Xt(sample_points, far_coef[t,], R=R)  
      #tmp <- tmp+rnorm(n, sd=sampling_noise) # add iid sampling noise to each observation point ##？
      #far_sample <- rbind(far_sample, tmp)
    #}
    U = far_Xt.fast(s = sample_points,R = R)
    #print(U) 
    #print(far_coef)                                        
    far_sample = far_coef%*%t(U)                                        
    #print(length(tmp)) 
    #print(dim(far_sample))                                         
    X_true <- far_sample[(totalT_train+1):totalT,] # later used for prediction
    mylist = list("X.test" = X_true, "Xtotal" = far_sample,"Atrue" = A_true)
    return(mylist)
    }
                                            
#ridge
ridge.vectorize.d.2 = function(K,X,Xl,U1,U2,D1,D2,S0,sqrtK,lambda){
    n = nrow(K)
    #sqrtK = sqrtm(K)
    #S0 = sqrtK%*%Xl/n
    #decomp.S0 = svd(S0)
    #U1 = decomp.S0$u
    #D1 = decomp.S0$d
    #decomp.K = svd(K)
    #U2 = decomp.K$u
    #D2 = decomp.K$d
    W = matrix(0,nrow = n,ncol = n)
    for(i in 1:n){
        for(j in 1:n){
            W = W + U2[,j]%*%t(U2[,j])%*%sqrtK%*%X%*%t(Xl)%*%sqrtK%*%U1[,i]%*%t(U1[,i])/(lambda + D1[i]*D2[j])
        }
    }
    #W.seq = sapply(1:n,function(j) sapply(1:n,function(i) W.func(U2 = U2,sqrtK = sqrtK,X = X,Xl = Xl,
                                                         #U1 = U1,lambda = lambda,D1 = D1,D2 = D2,i = i,j = j)))
    #nc = ncol(W.seq)
    #nr = nrow(W.seq)                                      
    #for (time in 0:(nr/nc -1)){
        #tmp = W.seq[(nc*time+1):(nc*(time+1)),]
        #W = W + tmp
    #}                                      
    W = W/n
    return(W)
}
ridge.GCV = function(K,X,Xl,lambda,lambda.kernel){
    n = nrow(K)
    T = ncol(X) +1
    numerator.mat = matrix(0,nrow = n,ncol = T-1)
    trace.A = 0
    K = K + lambda.kernel*diag(n)
    sqrtK = sqrtm(K)
      K_sqr = sqrtK
      n = nrow(K)
    sqrtK = sqrtm(K)
    S0 = sqrtK%*%Xl/n
    decomp.S0 = svd(S0%*%t(S0))
    U1 = decomp.S0$u
    D1 = decomp.S0$d
    decomp.K = svd(K)
    U2 = decomp.K$u
    D2 = decomp.K$d
    for(i in 1:n){
        for(j in 1:n){
            numerator.mat = numerator.mat + sqrtK%*%U2[,j]%*%t(U2[,j])%*%sqrtK%*%X%*%t(S0)%*%U1[,i]%*%t(U1[,i])%*%S0/(lambda + D1[i]*D2[j])
            trace.A = trace.A + tr(U1[,i]%*%t(U1[,i])%*%S0%*%t(S0))*tr(U2[,j]%*%t(U2[,j])%*%K)/(lambda + D1[i]*D2[j])
        }
    }
    #num.seq = sapply(1:n,function(j) sapply(1:n,function(i) num.func(sqrtK = sqrtK,U2 = U2,S0 = S0,U1 = U1,lambda = lambda,
                                                                    #D1 = D1,D2 = D2,i = i,j = j)))
    #trA.seq = sapply(1:n,function(j) sapply(1:n,function(i) trace.A.func(U1 = U1,S0 = S0,U2 = U2,K = K,lambda = lambda,
                                                                     #D1 = D1,D2 = D2,i = i,j = j)))
    #nc = ncol(num.seq)
    #print(dim(num.seq))                                        
    #nr = nrow(num.seq)                                      
    #for (time in 0:(nr/(T-1) -1)){
        #tmp = num.seq[((T-1)*time+1):((T-1)*(time+1)),]
        #numerator.mat = numerator.mat + t(tmp)
    #}                                                                                
    numerator = norm(X - numerator.mat, type = "F")^2/(n*T)
    denominator = (1 - trace.A/(n*T))^2
    GCV.value = numerator/denominator
    return(GCV.value)
}
                                            
# fPCA
fd.2d = function(coef, basisobj1, basisobj2){
    if (is.null(coef) && is.null(basisobj1) && is.null(basisobj2)){
        basisobj1 <- basisfd()
        basisobj2 <- basisfd()
    } 
    if (is.null(coef)) 
        coef <- rep(0, basisobj[["nbasis"]])
    {
        if (!is.numeric(coef)) 
            stop("'coef' is not numeric.")
        else if (is.vector(coef)) {
            coef <- as.matrix(coef)
            if (identical(type, "constant")) 
                coef <- t(coef)
            coefd <- dim(coef)
            ndim <- length(coefd)
        }
        else if (is.matrix(coef)) {
            coefd <- dim(coef)
            ndim <- length(coefd)
        }
        else if (is.array(coef)) {
            coefd <- dim(coef)
            ndim <- length(coefd)
        }
        else stop("Type of 'coef' is not correct")
    }
    if (ndim > 3) 
        stop("'coef' not of dimension 1, 2 or 3")
    if (ndim > 1) 
        nrep <- coefd[2]
    else nrep <- 1
    if (ndim > 2) 
        nvar <- coefd[3]
    else nvar <- 1
    fdobj <- list(coefs = coef, basis1 = basisobj1, basis2 = basisobj2)
    oldClass(fdobj) <- "fd"
    fdobj
}

smooth.basis.2d = function(argvals=sample_points, y=t(far_sample[1:totalT_train,]), fdParob1=basis1, fdParob2=basis2, lambda = 0){
    basisobj1 <- fdParob1
    basisobj2 <- fdParob2
    nbasis1 <- basisobj1$nbasis
    nbasis2 <- basisobj2$nbasis
    Phi1mat = eval.basis(argvals[,1],basisobj1)
    Phi2mat = eval.basis(argvals[,2],basisobj2)
    #print(Phi1mat)
    #print(Phi2mat)
    Phimat = (Phi2mat %x% matrix(1,nrow = 1, ncol = ncol(Phi1mat)))*(matrix(1,nrow = 1, ncol = ncol(Phi2mat)) %x% Phi1mat)
    Smat = t(Phimat) %*% Phimat
    #print(eigen(Smat)$values)
    if (lambda > 0) {
        Jmat1 <- inprod(basisobj1, basisobj1)
        Jmat2 = inprod(basisobj2, basisobj2)
        DJ1mat <- eval.penalty(basisobj1, 1)
        D2J1mat = eval.penalty(basisobj1, 2)
        DJ2mat <- eval.penalty(basisobj2, 1)
        D2J2mat = eval.penalty(basisobj2, 2)
        Kmat = Jmat2 %x% D2J1mat + 2*DJ2mat %x% DJ1mat + D2J2mat %x% Jmat1
        Smat = Smat + lambda*Kmat
    }
    Smat = (Smat + t(Smat))/2
    #print(lambda)
    #print(Smat)
    #print(dim(Smat))
    Rmat = chol(Smat)
    Rmat.inv = solve(Rmat)
    coef = solve(Rmat, t(Rmat.inv)%*%t(Phimat)%*%y)
    coef <- as.matrix(coef)
    nbasis = nbasis1*nbasis2
    #fdnames1
    #fdnames2
    #print(is.null(basisobj1))
    fdobj <- fd.2d(coef[1:nbasis, ], basisobj1, basisobj2)
    smoothlist <- list(fd = fdobj, lambda = lambda)
}

eval.fd.2d = function(sample_points, x_eigenfunc){
    basis1 = x_eigenfunc$basis1
    basis2 = x_eigenfunc$basis2
    coef = x_eigenfunc$coef
    Phi1mat = eval.basis(sample_points[,1],basis1)
    Phi2mat = eval.basis(sample_points[,2],basis2)
    Phimat = (Phi2mat %x% matrix(1,nrow = 1, ncol = ncol(Phi1mat)))*(matrix(1,nrow = 1, ncol = ncol(Phi2mat)) %x% Phi1mat)
    X = Phimat %*% coef
    return(X)
}

fPCA2D = function (fdobj, nharm = 2, centerfns = FALSE,lambda) 
{
    if (!(is.fd(fdobj))) 
        stop("First argument is neither a functional data or a functional parameter object.")
    #meanfd <- mean.fd(fdobj)
    if (centerfns) {
        fdobj <- center.fd(fdobj)
    }
    coef <- fdobj$coefs
    coefd <- dim(coef)
    #print(coefd)
    ndim <- length(coefd)
    nrep <- coefd[2]
    #coefnames <- dimnames(coef)
    if (nrep < 2) 
        stop("PCA not possible without replications.")
    basisobj1 <- fdobj$basis1
    basisobj2 <- fdobj$basis2
    nbasis1 <- basisobj1$nbasis
    nbasis2 <- basisobj2$nbasis
    #type <- basisobj$type
    #harmbasis = basisobj
    #harmbasis <- harmfdPar$fd$basis# harmbasis1,2
    #nhbasis <- harmbasis$nbasis
    #Lfdobj <- harmfdPar$Lfd
    #lambda <- harmfdPar$lambda
    #if (ndim == 3) {
        #nvar <- coefd[3]
        #ctemp <- matrix(0, nvar * nbasis, nrep)
        #for (j in 1:nvar) {
            #index <- 1:nbasis + (j - 1) * nbasis
            #ctemp[index, ] <- coef[, , j]
        #}
    #}
    #else {
        nvar <- 1
        ctemp <- coef
    #}
    Jmat1 <- inprod(basisobj1, basisobj1)
    Jmat2 = inprod(basisobj2, basisobj2)
    Jmat = Jmat1 %x% Jmat2
    Lmat = Jmat
    if (lambda > 0) {
        DJ1mat <- eval.penalty(basisobj1, 1)
        D2J1mat = eval.penalty(basisobj1, 2)
        DJ2mat <- eval.penalty(basisobj2, 1)
        D2J2mat = eval.penalty(basisobj2, 2)
        Kmat = Jmat2 %x% D2J1mat + 2*DJ2mat %x% DJ1mat + D2J2mat %x% Jmat1
        Lmat <- Lmat + lambda * Kmat
    }
    Lmat <- (Lmat + t(Lmat))/2
    Mmat <- chol(Lmat)
    Mmatinv <- solve(Mmat)
    Wmat <- crossprod(t(ctemp))/nrep
    #Jmat = inprod(harmbasis, basisobj)
    MIJW = crossprod(Mmatinv, Jmat)
    if (nvar == 1) {
        Cmat = MIJW %*% Wmat %*% t(MIJW)
    }
    #else {
        #Cmat = matrix(0, nvar * nhbasis, nvar * nhbasis)
        #for (i in 1:nvar) {
            #indexi <- 1:nbasis + (i - 1) * nbasis
            #for (j in 1:nvar) {
                #indexj <- 1:nbasis + (j - 1) * nbasis
                #Cmat[indexi, indexj] <- MIJW %*% Wmat[indexi, 
                  #indexj] %*% t(MIJW)
            #}
        #}
    #}
    Cmat <- (Cmat + t(Cmat))/2
    result <- eigen(Cmat)
    eigvalc <- result$values
    eigvecc <- as.matrix(result$vectors)#[, 1:nharm])
    sumvecc <- apply(eigvecc, 2, sum)
    #eigvecc[, eigvalc < 0] <- -eigvecc[, eigvalc < 0]
    #print("eigvalc")
    #print(eigvalc)
    #print(which(eigvalc<0))
    #print(dim(eigvecc))
    eigvecc[, which(eigvalc<0)] <- -eigvecc[, which(eigvalc<0)]
    eigvecc = eigvecc[, 1:nharm]
    #eigvecc[, sumvecc < 0] <- -eigvecc[, sumvecc < 0]
    eigvalc <- abs(eigvalc)
    #varprop <- eigvalc[1:nharm]/sum(eigvalc)
    varprop <- eigvalc/sum(eigvalc)#eigvalc[1:nharm]/sum(eigvalc)
    if (nvar == 1) {
        harmcoef <- Mmatinv %*% eigvecc
    }
    #else {
        #harmcoef <- array(0, c(nbasis, nharm, nvar))
        #for (j in 1:nvar) {
            #index <- 1:nbasis + (j - 1) * nbasis
            #temp <- eigvecc[index, ]
            #harmcoef[, , j] <- Mmatinv %*% temp
        #}
    #}
    #harmnames <- rep("", nharm)
    #for (i in 1:nharm) harmnames[i] <- paste("PC", i, sep = "")
    #if (length(coefd) == 2) 
        #harmnames <- list(coefnames[[1]], harmnames, "values")
    #if (length(coefd) == 3) 
        #harmnames <- list(coefnames[[1]], harmnames, coefnames[[3]])
    harmfd <- fd.2d(harmcoef,basisobj1, basisobj2)
    harmscr <- t(ctemp)%*%Jmat%*%harmcoef#t(t(harmcoef)%*%Jmat%*%ctemp)
    #if (nvar == 1) {
        #harmscr <- inprod(fdobj, harmfd)
    #}
    #else {
        #harmscr <- array(0, c(nrep, nharm, nvar))
        #coefarray <- fdobj$coefs
        #harmcoefarray <- harmfd$coefs
        #for (j in 1:nvar) {
            #fdobjj <- fd(as.matrix(coefarray[, , j]), basisobj)
            #harmfdj <- fd(as.matrix(harmcoefarray[, , j]), basisobj)
            #harmscr[, , j] <- inprod(fdobjj, harmfdj)
        #}
    #}
    #print("varprop")
    #print(varprop)
    pcafd <- list(harmfd, eigvalc, harmscr, varprop)
    class(pcafd) <- "pca.fd"
    names(pcafd) <- c("harmonics", "values", "scores", "varprop")
    return(pcafd)
}


fFPE = function(npca_component.max,fpca.obj,totalT_train){
    fFPE.seq <- c() # selection criteria for number of pca components in Aue (2014)
    for(npca_component in 1:npca_component.max){
      x_scores_train <- fpca.obj$scores[,1:npca_component]
        #print(fpca.obj$scores)
      # Vector AR
      if(npca_component>1){
        var_x <- VAR(x_scores_train,p=1,type='none')
        A_var <- c()
        for(npca_index in 1:npca_component){
            #print(sum(var_x$varresult[[npca_index]]$coefficients))
          A_var <- rbind(A_var, var_x$varresult[[npca_index]]$coefficients)  
        }
          #print(tr(A_var))
        residual <- t(x_scores_train[-1,])-A_var%*%t(x_scores_train[-totalT_train,]) 
        SigmaZ <- residual%*%t(residual)/(totalT_train-1)
          #print("sigmaZ")
      }else{
        ar_x <- ar(x_scores_train, order.max=1, aic=F, demean=F)
          #print("tran.A")
          #print(ar_x$var.pred)
        A_var <- ar_x$ar
        SigmaZ <- matrix(ar_x$var.pred, 1, 1)
      }
        #print("d")
        #print(npca_component)
        #print(tr(SigmaZ))
        #print(sum(fpca.obj$values[-c(1:npca_component)]))
      fFPE_tmp <- (totalT_train+npca_component)/(totalT_train-npca_component)*tr(SigmaZ) + sum(fpca.obj$values[-c(1:npca_component)])
      fFPE.seq <- c(fFPE.seq, fFPE_tmp)
    }
    return(fFPE.seq)
}
                                            
## aue 2d
aue = function(sample_points,x_curve_train,lambda,nbasis,totalT_train){                                                   
npca_component.max <- nbasis   # the maximum number of npca_component   
fpca_max <- fPCA2D(fdobj = x_curve_train$fd, nharm=npca_component.max, centerfns=FALSE,lambda = lambda) 
varprop =  cumsum(fpca_max$varprop)
    
#fFPE.seq = fFPE(npca_component.max = npca_component.max,fpca.obj = fpca_max,totalT_train = totalT_train)   
    #print(fFPE.seq)
    # final selected fPCA model based on fFPE
    npca_component_selected <- which(varprop >0.999)[1]#which.min(fFPE.seq)
    #print(npca_component_selected)
    fpca <- fPCA2D(fdobj = x_curve_train$fd, nharm=npca_component_selected, centerfns=F,lambda = lambda)
    #print(fpca$varprop)
    #fpcatest <- fPCA2D(fdobj = x_curve_train$fd, nharm=2, centerfns=F,lambda = x_curve_train$lambda)
    #print(fpcatest$varprop)
#pca.fd(x_curve_train$fd, nharm=npca_component_selected, harmfdPar=fdPar(x_curve_train$fd), centerfns=F)
    #print("score")
    #print(fpca_max$scores)
    x_scores_train <- fpca_max$scores[,1:npca_component_selected]
    #print(x_scores_train)
    # selected Vector AR
    if(npca_component_selected>1){
      var_x <- VAR(x_scores_train,p=1,type='none')
      A_var <- c()
      for(npca_index in 1:npca_component_selected){
        A_var <- rbind(A_var, var_x$varresult[[npca_index]]$coefficients)
      }
      residual <- t(x_scores_train[-1,])-A_var%*%t(x_scores_train[-totalT_train,])
    }else{
      ar_x <- ar(x_scores_train, order.max=1, aic=F, demean=F)
      A_var <- ar_x$ar
    }
    
    # estimating transition operator A
    x_eigenfunc <- fpca$harmonics
    
    fpcalist = list("eigenfunc" = x_eigenfunc, "varmat" = A_var)
    return(fpcalist)
}

cv_fpca.2d <- function(cv_folds = 5,lambda_pen=lambda_pen,lambda_smooth,sample_points,far_sample,basis1,basis2,npoints_approx){
  # function for estimating forecating error via cross validation for a given lambda_pen  
  # X_all: observation matrix
  # Xl_all: lagged observation matrix
    
    nbasis1 = basis1$nbasis
    nbasis2 = basis2$nbasis
    nbasis = nbasis1*nbasis2
    
  total_train <- dim(far_sample)[1]
  cv_indices <- createFolds(1:total_train, k=cv_folds)
  
  cv_err <- c()
  for(cv_fold_index in 1:cv_folds){
    cv_index_test <- sort(cv_indices[[cv_fold_index]])
    cv_index_train <- sort(unlist(cv_indices)[!unlist(cv_indices)%in%cv_index_test])
      
    x_curve_train <- smooth.basis.2d(argvals=sample_points, y=t(far_sample[cv_index_train,]), fdParob1=basis1, fdParob2=basis2, lambda = lambda_smooth)
    x_curve_test <- smooth.basis.2d(argvals=sample_points, y=t(far_sample[cv_index_test,]), fdParob1=basis1, fdParob2=basis2, lambda = lambda_smooth)  
    aue.output = aue(sample_points = sample_points,x_curve_train = x_curve_train,lambda = lambda_pen,nbasis = nbasis,totalT_train = length(cv_index_train))  
    x_eigenfunc = aue.output$eigenfunc
    A_var = aue.output$varmat
      
    # one step ahead prediction
    #npoints_approx <- 1e3
    approx_points <- cbind(rep(1:npoints_approx,each = npoints_approx),rep(1:npoints_approx,npoints_approx))/npoints_approx
      #seq(0,1,length.out=npoints_approx)
      #print( t(eval.fd.2d(approx_points, x_eigenfunc)))
    x_scores_test <- t(eval.fd.2d(approx_points, x_eigenfunc))%*%eval.fd.2d(approx_points, x_curve_test$fd)/nrow(approx_points)
    if(is.null(dim(A_var))){
      x_score_forecast <- A_var*x_scores_test
    }else{
      x_score_forecast <- A_var%*%x_scores_test
    }
    X_forecast_aue <- t(eval.fd.2d(sample_points, x_eigenfunc)%*%x_score_forecast)  
      
    X_true <- far_sample[cv_index_test,]#X_all[,cv_index_test]
      
      ########### need test
      
    cv_err <- c(cv_err, mean(apply(X_true-X_forecast_aue, 2, mse)))
      
  }
  return(mean(cv_err))
}

far.tr = function(X,lambda,X.test){
  #alpha n^2 * T
  n = ncol(X)
  Time = nrow(X) - 1
  X.response = X[2:(Time+1),]
  X.variable = X[1:Time,]
  vec.X.variable = as.vector(matrix(X.variable,nrow = Time*n, ncol = 1))
  #K = outer(vec.X.variable,vec.X.variable,kernel.poly) + 10^(-1)*diag(Time*n)
  #print(K)
  #print(dim(K))  
  #output = sapply(1:n, function(j) sapply(1:n, function(k) inner.sam(X.response[,j],k,K,alpha[seq(j,nrow(alpha),by = n),],eps,lambda,niter)))
  output = sapply(1:n, function(j) inner.far.tr(X.response[,j],X.variable,lambda))
  #print("output")
  #print(output)
  #print(dim(output))
  A_est = n*t(output)
  X.plus.seq <- X.test%*%t(A_est)/n
  #X.plus = apply(X.plus.seq, 2, sum)
  rt = list("X.pred" = X.plus.seq,"A.est" = A_est)
  return(rt)
}

inner.far.tr = function(X.j, X.variable,lambda){
    #print(dim(alpha.seq))
  Time = length(X.j)#nrow(alpha.seq)
  n = ncol(X.variable)
  #print(dim(X.bcd))
  #print(X.j)
  group = rep(1:1, each = n)
  fit = gglasso(x = X.variable, y = X.j, group = group, loss = "ls", intercept = FALSE, nlambda = 100)
  fit.cv <- cv.gglasso(x = X.variable, y = X.j, group = group, nfolds = 10)
  # Best lambda
  best_lambda_fit.cv <- fit.cv$lambda.1se
  # Final coefficients of variables
  beta = coef(object = fit, s = best_lambda_fit.cv)
  coef = beta[-1]
  #fit = grpreg(X.variable, X.j, group, lambda = lambda,penalty="grLasso")
  #plot(fit)
  #fit$beta
  #coef = fit$beta#coef(fit, lambda=lambda)
  #print(coef)
  #print(alpha.seq)
  return(coef)
}
                  
                  
basis.region = function(nregion,R){
    u.q = c()
    for(k in 1:R){
        if (k == 1){
            tmp.u = 1/4
        }
        else{
            lambda = 2*(k-1)*pi
            S = 2*cos((0.5 + 0.5*mod(nregion, 3))*lambda) - cos((1 + 0.5*mod(nregion, 3))*lambda) - cos((0.5*mod(nregion, 3))*lambda)
            tmp.u = sqrt(2)*S/(lambda^2)
        }
      u.q = c(u.q,tmp.u)    
  }
    return(u.q)
}

Lambda.low.rank = function(Nregion,R){
    Lambda = matrix(0,nrow = R, ncol = R)
    for (nregion in 1:Nregion){
        #print(nregion)
        u.tilde = basis.region(nregion,R)
        #print(u.tilde %o% u.tilde)
        Lambda = Lambda + u.tilde %o% u.tilde
    }
    return(Lambda)
}
data.generation.low.rank = function(totalT_train,n,R,sample_points,noise_range,Nregion){
    burnin <- 100 # burnin stage for simulation of FAR(1)
  #noise_range <- 0.1 # functional noise, Z is U(-noise_range, noise_range)
  sampling_noise <- 0
  totalT_train <- totalT_train # number of curves
  totalT_test <- round(totalT_train*0.2)
  totalT <- totalT_train + totalT_test
  init_zero <- T # whether use 0 matrix as initial value of optimization
    # Plot contour of true transition operator A at (s1,s2,...,sn) by (s1,s2,...,sn)
    Lambda = c()
    for (nregion in 1:Nregion){
        u.tilde = basis.region(nregion,R)
        Lambda = Lambda + u.tilde %o% u.tilde
    }
    A_true <- sapply(1:n,function(j) sapply(1:n,function(i) A_transition(sample_points[i,],sample_points[j,],Lambda)))
                                            #outer(sample_points, sample_points, A_transition_vec, Lambda=Lambda)
  far_coef <- c() # xt coefficient
    for(t in 1:(totalT+burnin)){
      z <- runif(R, min=-noise_range, max=noise_range)
      if(t==1){
        far_coef <- rbind(far_coef, z)
      }else{
        tmp <- Lambda%*%far_coef[(t-1),]+z
          #print(tmp)
        far_coef <- rbind(far_coef, c(tmp))
      }
    }                                        
    far_coef <- far_coef[-c(1:burnin),]
    rownames(far_coef) <- NULL
    far_coef_onestep <- t(Lambda%*%t(far_coef)) # best one step ahead estimation of far_coef
    U = far_Xt.fast(s = sample_points,R = R)
    noise.coef <- matrix(runif(R*totalT, min=-noise_range, max=noise_range),nrow = totalT, ncol = R)                                       
    far_noise = noise.coef %*% t(U)                                        
    #far_sample = far_coef%*%t(U)                                                                                
    X_true <- far_sample[(totalT_train+1):totalT,] # later used for prediction
    mylist = list("X.test" = X_true, "Xtotal" = far_sample,"Atrue" = A_true, "Lambda" = Lambda)
    return(mylist)
}                  
                  
basis.region.mul = function(nregion,R,Nregion,ver,hor){
    u.q = c()
    for(k in 1:R){
        if (k == 1){
            tmp.u = 1/Nregion
        }
        else{
            lambda = 2*(k-1)*pi
            #S = 2*cos((0.5 + 0.5*mod(nregion, 3))*lambda) - cos((1 + 0.5*mod(nregion, 3))*lambda) - cos((0.5*mod(nregion, 3))*lambda)
            S = (sin(lambda*hor[2])-sin(lambda*hor[1]))*(sin(lambda*ver[2])-sin(lambda*ver[1])) - (cos(lambda*hor[2])-cos(lambda*hor[1]))*(cos(lambda*ver[2])-cos(lambda*ver[1]))
            tmp.u = sqrt(2)*S/(lambda^2)
        }
      u.q = c(u.q,tmp.u)    
  }
    return(u.q)
}

Lambda.low.rank.mul = function(Nregion,R){
    Lambda = matrix(0,nrow = R, ncol = R)
    grid.on.side = sqrt(Nregion)
    length.of.grid = 1/grid.on.side
    for (nregion in 1:Nregion){
        ver = length.of.grid*c(ceiling(nregion/grid.on.side)-1,ceiling(nregion/grid.on.side))
        hor = length.of.grid*c(mod(nregion,grid.on.side)-1,mod(nregion,grid.on.side))
        #print(nregion)
        u.tilde = basis.region.mul(nregion,R,Nregion,ver,hor)
        #print(u.tilde %o% u.tilde)
        Lambda = Lambda + u.tilde %o% u.tilde
    }
    return(Lambda)
}
data.generation.low.rank.mul = function(totalT_train,n,R,sample_points,noise_range,Nregion){
    burnin <- 100 # burnin stage for simulation of FAR(1)
  #noise_range <- 0.1 # functional noise, Z is U(-noise_range, noise_range)
  sampling_noise <- 0
  totalT_train <- totalT_train # number of curves
  totalT_test <- round(totalT_train*0.2)
  totalT <- totalT_train + totalT_test
  init_zero <- T # whether use 0 matrix as initial value of optimization
    # Plot contour of true transition operator A at (s1,s2,...,sn) by (s1,s2,...,sn)
    Lambda = c()
    grid.on.side = sqrt(Nregion)
    length.of.grid = 1/grid.on.side
    for (nregion in 1:Nregion){
        ver = length.of.grid*c(ceiling(nregion/grid.on.side)-1,ceiling(nregion/grid.on.side))
        hor = length.of.grid*c(mod(nregion,grid.on.side)-1,mod(nregion,grid.on.side))
        u.tilde = basis.region.mul(nregion,R,Nregion,ver,hor)
        Lambda = Lambda + u.tilde %o% u.tilde
    }
    A_true <- sapply(1:n,function(j) sapply(1:n,function(i) A_transition(sample_points[i,],sample_points[j,],Lambda)))
                                            #outer(sample_points, sample_points, A_transition_vec, Lambda=Lambda)
  far_coef <- c() # xt coefficient
    for(t in 1:(totalT+burnin)){
      z <- runif(R, min=-noise_range, max=noise_range)
      if(t==1){
        far_coef <- rbind(far_coef, z)
      }else{
        tmp <- Lambda%*%far_coef[(t-1),]+z
          #print(tmp)
        far_coef <- rbind(far_coef, c(tmp))
      }
    }                                        
    far_coef <- far_coef[-c(1:burnin),]
    rownames(far_coef) <- NULL
    far_coef_onestep <- t(Lambda%*%t(far_coef)) # best one step ahead estimation of far_coef
    U = far_Xt.fast(s = sample_points,R = R)
    noise.coef <- matrix(runif(R*totalT, min=-noise_range, max=noise_range),nrow = totalT, ncol = R)                                       
    far_noise = noise.coef %*% t(U)                                        
    #far_sample = far_coef%*%t(U)                                                                                
    X_true <- far_sample[(totalT_train+1):totalT,] # later used for prediction
    mylist = list("X.test" = X_true, "Xtotal" = far_sample,"Atrue" = A_true, "Lambda" = Lambda)
    return(mylist)
}                  
                  
