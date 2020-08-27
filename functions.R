
#---------------------------------------------------------------------------------------------#
# functions definitions
#---------------------------------------------------------------------------------------------#
library(DTRlearn2)
library(prodlim)


## function to simulate from RBM by Gibbs sampling
# each categorical Yj takes value from 1: maxcat[j]
gib_rbm_general2 = function(nsim, ny, nz, maxcat, a, b, w) {
  
  ysim = matrix(0, ncol=ny, nrow=nsim)
  zsim = matrix(0, ncol=nz, nrow=nsim)
  z_p = matrix(0, ncol=nz, nrow=nsim)
  
  # set random initial values of y 
  ysim[1,] = rbinom(ny, 1, 0.5)
  
  for(m in 1:(nsim-1)) {
    z_p[m,] = exp(b + t(ysim[m,]) %*% w) / (1 + exp(b + t(ysim[m,]) %*% w)) # note based on ysim from the same iteration, hence this is the poterior mean
    zsim[m,] = rbinom(nz, 1, z_p[m,])
    y_p = matrix(0 , nrow=ny, ncol = max(maxcat)+1)
    y_p[, 1] = 1 # exp(0) for Y=0
    wz = w %*% as.matrix(zsim[m,])
    ##
    a_iter = 1
    for(j in 1:ny) {
      for(p in 1:maxcat[j]) {
        y_p[j,p+1] = exp( a[a_iter] + wz[j] * p )  # don't foget to *p (the value of Y)
        a_iter = a_iter + 1
      }
    }
    y_p = y_p / apply(y_p, 1, sum) # row divided by scalar 
    for(j in 1:ny) {
      ysim[(m+1),j] = which.max( rmultinom(1, 1, y_p[j,]) ) - 1
    }
  }
  z_p[nsim,] = exp(b + ysim[nsim,] %*% w) / (1 + exp(b + ysim[nsim,] %*% w) )
  zsim[nsim,] = rbinom(nz, 1, z_p[nsim])
  
  list(ysim=ysim, zsim=zsim, z_p = z_p)
}

## fucntion to subsample n from the start iteration, and subsample every m    
subsample = function(dat, start, n, m) { # subsample every m iterations 
  ncol = dim(dat)[2]
  output = matrix(numeric(), ncol=ncol, nrow=n)
  index = start
  for(i in 1:n) {
    output[i, ] = dat[index, ]
    index = index + m
  }
  output 
}

## function to calculate the log-likelihood of the RBM model
cal_loglik2 = function(a, b, w, ys, zs, n_ys, n_zs, obs_y0, maxcat) {  #add maxcat
  nobs = dim(obs_y0)[1]
  a_ys = rep(0, n_ys)
  n_a = length(a)      
  ny = dim(w)[1]   
  sum_indY = rep(0, n_a)   
  a_iter = 1
  for(j in 1:ny) {
    for(p in 1:maxcat[j]) {
      a_ys = a_ys + a[a_iter] * (ys[,j] == p)
      sum_indY[a_iter] = sum(obs_y0[,j]==p)
      a_iter = a_iter + 1
    }
  }
  b_zs = b %*% t(zs)
  wyz = as.matrix(ys) %*% w %*% t(zs)
  C_terms = exp( matrix(rep(a_ys, n_zs), ncol=n_zs) + matrix(rep(b_zs, each=n_ys), ncol=n_zs) + wyz )
  C = sum(C_terms)
  E0 = exp( matrix(rep(b_zs, nobs), nrow=nobs, byrow=T) + obs_y0 %*% w %*% t(zs)  )  # add b_zs to each line of it
  loglik = - nobs * log(C)  +  sum(a*sum_indY) +  sum(log(apply(E0, 1, sum)))
  loglik 
}


## function to simulate the full graphical model  
sim_graph2 = function(seed, n_tot, nsim, start, m,  ny, nz, nx, max_cat, a, b, w, theta0, theta1, theta2, theta3, theta4, theta5) {
  set.seed(seed)
  y0sim = matrix(0, ncol=ny, nrow=n_tot)
  p_z0 = matrix(0, ncol=nz, nrow=n_tot)
  z0sim = matrix(0, ncol=nz, nrow=n_tot)
  y1sim = matrix(0, ncol=ny, nrow=n_tot)
  z1sim = matrix(0, ncol=nz, nrow=n_tot)
  
  # simulate Y0 and Z0 first 
  sim_y0z0 = gib_rbm_general2(nsim, ny, nz, max_cat, a, b, w)
  y0sim = subsample(sim_y0z0$ysim, start = start, n=n_tot, m=m)
  z0sim = subsample(sim_y0z0$zsim, start = start, n=n_tot, m=m)
  z0_p = subsample(sim_y0z0$z_p, start=start, n=n_tot, m=m)
  
  # simulate X 
  X = matrix(rnorm(n_tot*nx, 0, sd=1), nrow=n_tot)
  # simulate treatment trt
  trt = 2 * rbinom(n_tot, 1, 0.5) - 1 
  
  # simulate Z1 based on Z0, X, trt
  theta_linear = theta0 + theta1 %*% t(z0sim) + theta2 %*% t(X) + matrix(theta3,nrow=nz) %*% trt + theta4 %*% t(trt * z0sim) + theta5 %*% t(trt * X)
  p_z1 = t( exp(theta_linear) / (1 + exp(theta_linear)) )
  z1sim = matrix(rbinom(n_tot*nz, 1, p_z1), ncol=nz) 
  
  # simualte Y1 based on Z1, based on RBM 
  y1sim = predict_y3 (max_cat, a, w, z1sim)   
  
  y0sim_sum = rowSums(y0sim)
  y1sim_sum = rowSums(y1sim)
  z0sim_sum = rowSums(z0sim)
  z1sim_sum = rowSums(z1sim)
  
  # counterfactural trt and simulate Z1 under the counterfactural trt 
  trt_c = (-1) * trt 
  theta_linear_c = theta0 + theta1 %*% t(z0sim) + theta2 %*% t(X) + matrix(theta3,nrow=nz) %*% trt_c + theta4 %*% t(trt_c * z0sim) + theta5 %*% t(trt_c * X)
  p_z1_c = t( exp(theta_linear_c) / (1 + exp(theta_linear_c)) )
  z1sim_c = matrix(rbinom(n_tot*nz, 1, p_z1_c), ncol=nz)
  
  # the optimal treatment for each person based on Z1sum 
  trt_better = rowSums(z1sim) < rowSums(z1sim_c)
  trt_c_better = rowSums(z1sim) > rowSums(z1sim_c)
  z1sum_eq = rowSums(z1sim) == rowSums(z1sim_c)
  
  opt_trt = 2 * rbinom(n_tot, 1, 0.5) - 1
  opt_trt[trt_better] = trt[trt_better]
  opt_trt[trt_c_better] = trt_c[trt_c_better]
  opt_trt2 = opt_trt
  opt_trt2[z1sum_eq] = 0 
  
  opt_z1sum = rowSums(z1sim) * trt_better + rowSums(z1sim_c) * trt_c_better + rowSums(z1sim) * z1sum_eq
  
  # Z1sum under 1 and -1 
  z1sum_both = matrix(NA, nrow=n_tot, ncol=2)
  z1sum_both[,1] = rowSums(z1sim) * (trt==1) + rowSums(z1sim_c) * (trt_c==1)
  z1sum_both[,2] = rowSums(z1sim) * (trt==-1) + rowSums(z1sim_c) * (trt_c==-1)
  
  
  list(nobs=n_tot, X=X, trt=trt, y0_obs=y0sim, y1_obs=y1sim, z0=z0sim, z0_p=z0_p, z1=z1sim, y0sum=y0sim_sum, y1sum=y1sim_sum, z0sum=z0sim_sum, z1sum=z1sim_sum,  a=a, b=b, w=w, theta0=theta0, theta1=theta1, theta2=theta2, theta3=theta3, theta4=theta4, theta5=theta5, opt_trt=opt_trt, opt_trt2=opt_trt2, opt_z1sum=opt_z1sum, z1sum_both=z1sum_both, z1_p=p_z1, seed=seed)  
}   

## function to simulate Ys based on conditional distribution given Zs
predict_y3 = function(maxcat, a, w, z) { 
  if(is.matrix(z)) nobs = dim(z)[1]    # z: each line is 1 obs
  if(is.vector(z)) nobs = length(z)
  ny = dim(w)[1]
  y_pred = matrix(numeric(), nrow=nobs, ncol=ny)
  
  wz =  as.matrix(z) %*% t(w)
  for(i in 1:nobs) {
    y_p = matrix(0, nrow=ny, ncol = max(maxcat)+1)
    y_p[,1] = 1
    a_iter = 1
    for(j in 1:ny) {
      for(p in 1:maxcat[j]) {
        y_p[j,p+1] = exp( a[a_iter] + wz[i,j] * p )  # don't foget to *p 
        a_iter = a_iter + 1
      }
      y_pred[i,j] = which.max( rmultinom(1, 1, y_p[j,]) ) - 1  # differs here, with randomness
    }
  }
  y_pred
}

## function to run one iteration of RBM 
do_rbm_train2 = function(a, b, w, momentum_a=0, momentum_b=0, momentum_w=0, ys, zs, n_ys, n_zs, obs_y0, maxcat, momentum, lrate, fix_w=NULL, fix_wsign=NULL)  {
  nobs = dim(obs_y0)[1]
  n_a = sum(maxcat)
  
  a_ys = rep(0, n_ys)
  a_iter = 1
  for(j in 1:ny) {
    for(p in 1:maxcat[j]) {
      a_ys = a_ys + a[a_iter] * (ys[,j] == p)
      a_iter = a_iter + 1
    }
  }
  b_zs = b %*% t(zs)
  wyz = as.matrix(ys) %*% w %*% t(zs)  # yTWz (ys on the column)
  
  C_terms = exp( matrix(rep(a_ys, n_zs), ncol=n_zs) + matrix(rep(b_zs, each=n_ys), ncol=n_zs) + wyz )
  C = sum(C_terms)
  
  #----------- Calculate the derivatives (in the vector form)
  # derivative wrt to alpha  
  sum_Cy = rep(0, n_a)      # for first term in gradient a
  sum_indY = rep(0, n_a)    # for second term in log-likelihood 
  
  C_terms_sumz = apply(C_terms, 1, sum)  # C_terms to sum over z 
  a_iter = 1
  for(j in 1:ny) {
    for(p in 1:maxcat[j]) {
      sum_Cy[a_iter] = sum( C_terms_sumz * (ys[,j] == p) ) # for first term in gradient a
      sum_indY[a_iter] = sum(obs_y0[,j]==p)
      a_iter = a_iter + 1
    }
  }
  d_a = - nobs/C * sum_Cy + sum_indY
  
  # calculate exp(beta^T z0 + y0i^T w z0) for all z0, all subjects i. Note here use zs for z0 as all possible zs.
  E0 = exp( matrix(rep(b_zs, nobs), nrow=nobs, byrow=T) + obs_y0 %*% w %*% t(zs)  )  # add b_zs to each line of it
  
  # derivative wrt to beta 
  d_b = - nobs/C * apply(apply(C_terms, 2, sum)*zs, 2, sum) + apply( E0 %*% as.matrix(zs) / apply(E0, 1, sum), 2, sum )   # first term: sum over y first 
  
  # derivative wrt to w 
  E0_y0z_E0 = list() 
  for(i in 1:nobs) {
    E0_y0z_E0[[i]] = obs_y0[i,] %*% t(E0[i,]) %*% as.matrix(zs) / apply(E0, 1, sum)[i]
  }
  d_w = - nobs/C  *  t(ys) %*% C_terms %*% as.matrix(zs) + #do sum at the same time
    Reduce('+', E0_y0z_E0)  
  
  #-------- calucalte the log likelihood 
  loglik = - nobs * log(C)  +  sum(a * sum_indY) +  sum(log(apply(E0, 1, sum)))
  
  #-------- gradient ascent with momentum
  momentum_a = momentum * momentum_a + d_a
  momentum_b = momentum * momentum_b + d_b
  momentum_w = momentum * momentum_w + d_w
  
  a = a + momentum_a * lrate
  b = b + momentum_b * lrate
  w = w + momentum_w * lrate
  
  if(!is.null(fix_w)) {
    for(k in 1:nz) {
      w_temp = w[fix_w[k], k]
      if(w_temp * fix_wsign[k] <0) {
        print(c("flip", k, w_temp))
        w[,k] = -1 * w[,k]
        b[k] = -1 * b[k]
      }
    }
  }
  
  list(loglik=loglik, a=a, b=b, w=w, momentum_a=momentum_a, momentum_b=momentum_b, momentum_w=momentum_w, lrate=lrate)
}

## function to run iterations of RBMs under one lrate
rbm_train2 = function(ny, nz, n_cat, a0, b0, w0,  obs_y0, lrate, lrate_scale = 1, momentum, niter, batchsize, outlast=F, fix_w=NULL, fix_wsign=NULL) {
  nobs = dim(obs_y0)[1]
  numbatches = nobs / batchsize
  maxcat = n_cat - 1
  n_a = sum(maxcat)
  
  a = a0
  b = b0
  w = w0
  
  n_ys = prod(n_cat)   
  n_zs = 2^nz          
  
  ys_matrix_cat = matrix(NA, nrow=max(n_cat), ncol=length(n_cat))
  for(i in 1:ny) {
    ys_matrix_cat[1:n_cat[i],i] = 0:(n_cat[i]-1)
  }
  ys = na.omit( expand.grid(data.frame(ys_matrix_cat)) )
  
  zs_matrix_cat = matrix(rep(0:1, nz), ncol=nz)
  zs = na.omit( expand.grid(data.frame(zs_matrix_cat)) )
  
  loglike = rep(0, niter)
  loglike[1] = cal_loglik2(a, b, w, ys, zs, n_ys, n_zs, obs_y0, maxcat) # can be slow
  loglik_opt = loglike[1]
  print(loglike[1])
  
  # set initial values
  momentum_a = rep(0, n_a)
  momentum_b = rep(0, nz)
  momentum_w = matrix(0, ncol = nz, nrow = ny)
  a_opt = a
  b_opt = b
  w_opt = w
  
  ### iterate the process of gradient asent 
  for(m in 2:niter) {  # each iteration on the whole sample 
    
    randperm <- sample(1:nobs, nobs)  
    if (numbatches >= 1) {
      for (l in 1:numbatches) {  # gradient ascent on every batch
        batch_obs_y0 <- obs_y0[randperm[((l - 1) * batchsize + 1):(l * batchsize)], ]
        batch_n = dim(obs_y0)[1]
        rbm <- do_rbm_train2 (a, b, w, momentum_a, momentum_b, momentum_w, ys, zs, n_ys, n_zs, batch_obs_y0, maxcat, momentum, lrate/batch_n, fix_w=fix_w, fix_wsign = fix_wsign)  
        momentum_a = rbm$momentum_a; momentum_b = rbm$momentum_b; momentum_w = rbm$momentum_w;
        a = rbm$a; b = rbm$b; w = rbm$w;
      }
    }
    if (numbatches > as.integer(numbatches)) {  
      batch_obs_y0 <- obs_y0[randperm[(as.integer(numbatches) * batchsize): nobs], ]
      rbm <- do_rbm_train2 (a, b, w, momentum_a, momentum_b, momentum_w, ys, zs, n_ys, n_zs, batch_obs_y0, maxcat, momentum, lrate/batch_n, fix_w=fix_w, fix_wsign = fix_wsign)  
      momentum_a = rbm$momentum_a; momentum_b = rbm$momentum_b; momentum_w = rbm$momentum_w;
      a = rbm$a; b = rbm$b; w = rbm$w;
    }
    
    loglik = cal_loglik2(a, b, w, ys, zs, n_ys, n_zs, obs_y0, maxcat)
    loglike[m] = loglik 
    print(c(m, loglik))
    
    if(loglik > loglik_opt)  { 
      loglik_opt = loglik
      a_opt = a
      b_opt = b
      w_opt = w
    }
    lrate = lrate * lrate_scale
  }
  list(loglik_opt=loglik_opt, a=a_opt, b=b_opt, w=w_opt)  # oglike=loglike,  lrate=lrate
}


#------------------- functions to run RBM with MCMC
## fit RBM with contrastive divergence
rbm.train_cat3 <- function(ny,nz, n_cat, a0, b0, w0, obs_y0, 
                           lrate, lrate_scale = 1, momentum=0.5, niter, batchsize=100, cd=1, 
                           cal_likelihood=F, outlast=F, momentum_a=NULL, momentum_b=NULL, momentum_w=NULL) {  
  if (!is.matrix(obs_y0))   stop("obs_y0 must be a matrix!")
  input_dim <- ncol(obs_y0)
  maxcat = n_cat - 1
  numepochs = niter
  
  n_a = length(a0)
  
  best = list()  # to store the best rbm achieved 
  
  if (cal_likelihood == T) { 
    n_ys = prod(n_cat)  # number of possible y values
    n_zs = 2^nz         # number of possible z values 
    
    # create the full sets of y and z (each row is one possible value) 
    ys_matrix_cat = matrix(NA, nrow=max(n_cat), ncol=length(n_cat))
    for(i in 1:ny) {
      ys_matrix_cat[1:n_cat[i],i] = 0:(n_cat[i]-1)
    }
    ys = na.omit( expand.grid(data.frame(ys_matrix_cat)) )  # vector memory exhausted (limit reached?) here with ny=17
    print("ys generated: dimension")
    print(dim(ys))
    
    zs_matrix_cat = matrix(rep(0:1, nz), ncol=nz)
    zs = na.omit( expand.grid(data.frame(zs_matrix_cat)) )
    print("zs generated: dimension")
    print(dim(zs))
    
    loglik = numeric(numepochs)
    loglik_opt = cal_loglik2 (a0, b0, w0, ys, zs, n_ys, n_zs, obs_y0, maxcat) # change to "cal_loglik2". was t(w0)
    print(loglik_opt)
  }

  
  rbm <- list(
    size = c(input_dim, nz),
    W = t(w0),
    vW = matrix(rep(0,nz*input_dim), c(nz,input_dim)),  # record the momentum from last iteration
    B = a0,
    vB = rep(0,n_a),  
    C = b0,
    vC = rep(0,nz),
    learningrate = lrate,
    learningrate_scale = lrate_scale,
    momentum = momentum,
    cd=cd,
    e = numeric(0)  # MSE of reconstruction error between sample and model simulated in each minibatch
  )
  if(!is.null(momentum_a)) rbm$vB = momentum_a    # add
  if(!is.null(momentum_b))  rbm$vC = momentum_b
  if(!is.null(momentum_w)) rbm$vW = momentum_w
  
  m <- nrow(obs_y0);
  numbatches <- m / batchsize;
  s <- 0
  for(i in 1:numepochs){
    # print(i)
    randperm <- sample(1:m,m)
    if(numbatches >= 1){
      for(l in 1 : numbatches){
        s <- s + 1
        batch_x <- obs_y0[randperm[((l-1)*batchsize+1):(l*batchsize)], ] 
        rbm <- do.rbm.train_mod2 (rbm,batch_x,s, ny, n_cat)
      }
    }
    #last fraction of sample
    if(numbatches > as.integer(numbatches)){
      batch_x <- obs_y0[randperm[(as.integer(numbatches)*batchsize):m], ]      
      s <- s + 1
      rbm <- do.rbm.train_mod2 (rbm,batch_x,s, ny, n_cat)
    }
    rbm$learningrate <- rbm$learningrate * rbm$learningrate_scale;
    
    # add 
    if(cal_likelihood==T) {
      loglik[i] = cal_loglik2 (rbm$B, rbm$C, t(rbm$W), ys, zs, n_ys, n_zs, obs_y0=obs_y0, maxcat)
      print(loglik[i])
      if (outlast==F & loglik[i] > loglik_opt) { 
        loglik_opt = loglik[i]
        best = rbm
      }
    }
  } # end of empochs
  if(length(best)==0) best = rbm  # for either no better likelhood found or no likelihood computed 
  if(cal_likelihood==T) {
    best$loglik_opt = loglik_opt
    best$loglik = loglik
  }
  # make output the same format as the exact 
  best_output <- list(
    a = best$B,
    b = best$C,
    w = t(best$W),
    cd = best$cd,
    recon_error = best$e,  # add followings
    momentum_a = best$vB,  
    momentum_b = best$vC, 
    momentum_w = best$vW
  )
  best_output
}


## function to run each gradient ascent iteration
do.rbm.train_mod2 <- function(rbm, batch_x, s, ny, n_cat){
  m <- nrow(batch_x)
  v1 <- batch_x       # start from the observed data 
  h1 <- binary.state(rbm.up_mod(rbm, v1))  # simulate based on the current params and observed data
  
  vn <- v1
  hn <- h1   
  for(c in 1:rbm$cd){  # iterate CD
    vn <- rbm.down_mod2 (rbm, hn, ny, n_cat)
    hn <- binary.state( rbm.up_mod (rbm, vn)) 
  }
  
  dW <- (t(h1) %*% v1 - t(hn)  %*% vn) / m  
  dW <- rbm$learningrate * dW
  rbm$vW <- rbm$vW * rbm$momentum + dW  # add previous momentum
  dW <- rbm$vW  
  rbm$W <- rbm$W + dW
  
  dB = numeric(length(rbm$B))
  a_iter = 1
  max_cat = n_cat - 1
  for(j in 1:ny) {
    for(p in 1:max_cat[j]){
      dB[a_iter] = mean(v1[,j]==p)  - mean(vn[,j]==p)
      a_iter = a_iter + 1 
    }
  }
  dB <- rbm$learningrate * dB
  rbm$vB <- rbm$vB * rbm$momentum + dB
  dB <- rbm$vB
  rbm$B <- rbm$B + dB
  
  dC <- colMeans(h1 - hn) 
  dC <- rbm$learningrate * dC
  rbm$vC <- rbm$vC * rbm$momentum + dC
  dC <- rbm$vC
  rbm$C <- rbm$C + dC
  
  rbm$e[s] <- sum((v1 - vn)^2)/m
  rbm
}


## Infer hidden units state by visible units
rbm.up_mod <- function(rbm, v){
  sum <- t( t(v %*% t(rbm$W)) + rbm$C )
  h <- 1/(1+exp(-sum))
  h
}


## Simulate visible vector by hidden units states
rbm.down_mod2 <- function(rbm, h, ny, n_cat){   
  maxcat = n_cat - 1  # note!! not + 1
  wz = h %*% rbm$W  
  v = t( apply(wz, 1, function(x) sim_ycat2 (rbm$B, ny, maxcat, x)) )
}

## function to calculate for each observation the probability of each category for each Y, and simulate Ys based on that
sim_ycat2 = function(a, ny, maxcat, wz_i) {
  y_p = matrix(0, nrow=ny, ncol = max(maxcat)+1)
  y_p[,1] = 1
  a_iter = 1
  for(j in 1:ny) {
    for(p in 1:maxcat[j]) {
      y_p[j, p+1] = exp( a[a_iter] + wz_i[j] * p )
      a_iter = a_iter + 1
    }
  }
  y_p = y_p / apply(y_p, 1, sum)
  v_sim = numeric(ny)
  # simulate based on the probabilities
  for(j in 1:ny) {
    v_sim[j] = which.max( rmultinom(1, 1, y_p[j,]) ) - 1
  }
  v_sim
}


## simulate hidden state based on probability
binary.state <- function(h){
  p <- matrix( runif(length(h),min=0,max=1) ,nrow=nrow(h),ncol=ncol(h))
  h[h>p] <- 1
  h[h<=p] <- 0
  h
}



#--------- function to run a series of RBMs with specified learning rates (lrates)
run_rbms = function(lrates=c(2, 1, 0.5), niters=c(rep(20,3)), ny, nz, n_cat, a0, b0, w0,  
                    obs_y0, lrate_scale = 1, momentum=0.5, batchsize=100, 
                    MCMC=F, cd=NULL, fix_w=NULL, fix_wsign=NULL, cal_likelihood=F, outlast=F) {
  
  steps = length(lrates)
  print(c("lrate", lrates[1]))
  
  if(MCMC == F)  {
    model = rbm_train2 (ny=ny, nz=nz,  n_cat=n_cat, a0=a0, b0=b0, w0=w0,  obs_y0=obs_y0, lrates[1], lrate_scale = lrate_scale, momentum = momentum, niter=niters[1], batchsize=batchsize, outlast=outlast, fix_w=fix_w, fix_wsign=fix_wsign)
  }
  else  {
    model = rbm.train_cat3 (ny=ny, nz=nz,  n_cat=n_cat, a0=a0, b0=b0, w0=w0,  obs_y0=obs_y0, lrates[1], lrate_scale = lrate_scale, momentum = momentum, niter=niters[1], batchsize=batchsize, cd=cd, cal_likelihood=cal_likelihood, outlast=outlast)
  }
  
  
  for(u in 2:steps) {
    print(c("step", u, "lrate", lrates[u]))
    if(MCMC == F) { 
      model = rbm_train2 (ny=ny, nz=nz,  n_cat=n_cat, a0=model$a, b0=model$b, w0=model$w,  obs_y0=obs_y0, lrates[u], lrate_scale = lrate_scale, momentum = momentum, niter=niters[u], batchsize=batchsize, outlast=outlast, fix_w=fix_w, fix_wsign=fix_wsign)
    }
    else {
      model = rbm.train_cat3 (ny=ny, nz=nz,  n_cat=n_cat, a0=model$a, b0=model$b, w0=model$w,  obs_y0=obs_y0, lrates[u], lrate_scale = lrate_scale, momentum = momentum, niter=niters[u], batchsize=batchsize, cd=cd, cal_likelihood=cal_likelihood, outlast=outlast)
    }
  }
  model
}



## function to calculate gtilde(Y1)
# g_multiplier: is a vector of length nz s.t. g_multiplier^T z1 is the objective quantity to optimize
cal_g_tilda_Y1 = function(a, w, ny, nz, n_cat, obs_y1, g_multiplier) {   
  
  max_cat = n_cat - 1
  n_ys = prod(n_cat)
  n_zs = 2^nz
  
  zs_matrix_cat = matrix(rep(0:1, nz), ncol=nz)
  zs = na.omit( expand.grid(data.frame(zs_matrix_cat)) )
  
  # since condition on Z1, Y1 are independent. Hence, the joint conditional probability is the product of the individual ones. For each possible Z1, can construct a conditional probability matrix for Y1 (same structure as ys_matrix_cat)
  exp_a_pwzs = list()  # list length n_zs, one for each possible zs
  w_zs = w %*% t(zs)
  for(t in 1:n_zs) {
    exp_a_pwzs[[t]] = matrix(NA, nrow=max(n_cat), ncol=length(n_cat))
    exp_a_pwzs[[t]][1,] = 1
    a_iter = 1
    for(j in 1:ny) {
      for(p in 1:max_cat[j]) {
        exp_a_pwzs[[t]][p+1,j] = exp(a[a_iter] + p * w_zs[j,t])
        a_iter = a_iter + 1
      }
    }
  }
  prob_y_zs = matrix(numeric(), nrow=n_ys, ncol=n_zs)  
  
  for(t in 1:n_zs) {
    prob_y_zs_t_ind = t( t(exp_a_pwzs[[t]]) / colSums(exp_a_pwzs[[t]], na.rm=T) )  # temp var for condprob of Y given zs=zs[t,]
    prob_y_zs_t = na.omit( expand.grid(data.frame(prob_y_zs_t_ind)) ) # temp var for cond prob of each Y (col=ny) for each ys(row=n_ys), for zs[t,]
    prob_y_zs[,t] = apply(prob_y_zs_t, 1, prod)
    # sum(prob_y_zs[,t])  # 1
  }
  ginv_prob_y_zs = ginv(prob_y_zs)  
  remove(prob_y_zs)
  
  G = matrix( colSums(t(zs) * g_multiplier), nrow=1)
  G_sigma_ginv = G %*% ginv_prob_y_zs
  
  # if(outputZ1s==F) {  # output the g_multiplier applied to each Z1
  #   G = matrix( colSums(t(zs) * g_multiplier), nrow=1)  
  #   G_sigma_ginv = G %*% ginv_prob_y_zs
  # }
  # if(outputZ1s==T | rescaleZs==T)  G_sigma_ginv = t(zs) %*% ginv_prob_y_zs  # output each individual Z1k
  # remove(ginv_prob_y_zs)
  
  prod = c(1, cumprod(n_cat)[-ny])
  y1_index = colSums(t(obs_y1) * prod) + 1
  
  g_tilda_Y1 = G_sigma_ginv [, y1_index]
  # if(rescaleZs==T)  {g_tilda_Y1 = apply(2 * g_tilda_Y1 - 1, 2, function(x) sum(x * g_multiplier)) }
  
  remove(G_sigma_ginv)
  g_tilda_Y1
}


## function to predict Z (prob of being 1) based on Y 
predict_z = function(b, w, y) {  # y each line as 1 obs
  z_p = exp(b + y %*% w) / (1 + exp(b + y %*% w) )
  z_p
}





