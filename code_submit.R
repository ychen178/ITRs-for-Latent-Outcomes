
#---------------------------------------------------------------------------------------------#
# This is the organized code for the simulations in the paper 
# "Learning Individualized Treatment Rules for Multiple-Domain Latent Outcomes"
# Yuan Chen 
# June 2019
#---------------------------------------------------------------------------------------------#
library(DTRlearn2)
library(prodlim)  # for row.match() used in gtilde(Y1) calculation




#----------- function definitions 
## function to simulate from RBM by Gibbs sampling
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

## function to simulate the full graphical model  
sim_graph2 = function(seed, n_tot, nsim, start, m,  ny, nz, nx, max_cat, a, b, w, theta0, theta1, theta2, theta3, theta4, theta5, rbmfun, pred_yfun) {
  set.seed(seed)
  y0sim = matrix(0, ncol=ny, nrow=n_tot)
  p_z0 = matrix(0, ncol=nz, nrow=n_tot)
  z0sim = matrix(0, ncol=nz, nrow=n_tot)
  y1sim = matrix(0, ncol=ny, nrow=n_tot)
  z1sim = matrix(0, ncol=nz, nrow=n_tot)
  
  # simulate Y0 and Z0 first 
  sim_y0z0 = rbmfun(nsim, ny, nz, max_cat, a, b, w)
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
  y1sim = pred_yfun (max_cat, a, w, z1sim)   
  
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
  # z1sum_eq_n = sum(z1sum_eq) 
  
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

## function to run one iteration of rbm 
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
  loglike[1] = cal_loglik2(a, b, w, ys, zs, n_ys, n_zs, obs_y0, maxcat) # very slow
  loglik_opt = loglike[1]
  print(loglike[1])
  
  # set initial values
  momentum_a = rep(0, n_a)  # changed here
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
  list(loglike=loglike, loglik_opt=loglik_opt, a=a_opt, b=b_opt, w=w_opt, lrate=lrate)
}

## function to solve iterations of RBMs under series of lrate
run_rbms = function(lrates=c(2, rep(1,2)), niters=c(20, rep(20,20)), ny, nz, n_cat, a0, b0, w0,  obs_y0, lrate_scale = 1, momentum=0.5, batchsize=100, outlast=F, fix_w=NULL, fix_wsign=NULL) {
  
  steps = length(lrates)
  print(c("lrate", lrates[1]))
  
  model = rbm_train2 (ny=ny, nz=nz,  n_cat=n_cat, a0=a0, b0=b0, w0=w0,  obs_y0=obs_y0, lrates[1], lrate_scale = lrate_scale, momentum = momentum, niter=niters[1], batchsize=batchsize, outlast=outlast, fix_w=fix_w, fix_wsign=fix_wsign)
  
  for(u in 2:steps) {
    print(c("step", u, "lrate", lrates[u]))
    model = rbm_train2 (ny=ny, nz=nz,  n_cat=n_cat, a0=model$a, b0=model$b, w0=model$w,  obs_y0=obs_y0, lrates[u], lrate_scale = lrate_scale, momentum = momentum, niter=niters[u], batchsize=batchsize, outlast=outlast, fix_w=fix_w, fix_wsign=fix_wsign)
  }
  model
}

## function to calculate the log-likelihood  
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

## function to calculate gtilde(Y1)
cal_g_tilda_Y1 = function(a, w,  ny, nz, n_cat, obs_y1, g_multiplier=rep(1,nz), outputZ1s=F, rescaleZs=F) {   
  
  max_cat = n_cat - 1
  n_ys = prod(n_cat)
  n_zs = 2^nz
  
  zs_matrix_cat = matrix(rep(0:1, nz), ncol=nz)
  zs = na.omit( expand.grid(data.frame(zs_matrix_cat)) )
  
  # since condition on Z1, Y1 are independent. Hence, the joint conditional probability is the product of the invidual ones. For each possible Z1, can construct a conditional probability matrix for Y1 (same structure as ys_matrix_cat)
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
    # print(t)
    prob_y_zs_t_ind = t( t(exp_a_pwzs[[t]]) / colSums(exp_a_pwzs[[t]], na.rm=T) )  # temp var for condprob of Y given zs=zs[t,]
    prob_y_zs_t = na.omit( expand.grid(data.frame(prob_y_zs_t_ind)) ) # temp var for cond prob of each Y (col=ny) for each ys(row=n_ys), for zs[t,]
    prob_y_zs[,t] = apply(prob_y_zs_t, 1, prod)
    # sum(prob_y_zs[,t])  # 1
  }
  ginv_prob_y_zs = ginv(prob_y_zs)  
  remove(prob_y_zs)
  
  if(outputZ1s==F) {  # output the Z1sum
    G = matrix( colSums(t(zs) * g_multiplier), nrow=1)  
    G_sigma_ginv = G %*% ginv_prob_y_zs
  }
  if(outputZ1s==T | rescaleZs==T)  G_sigma_ginv = t(zs) %*% ginv_prob_y_zs  # output each individual Z1k
  remove(ginv_prob_y_zs)
  
  prod = c(1, cumprod(n_cat)[-ny])
  y1_index = colSums(t(obs_y1) * prod) + 1
  
  g_tilda_Y1 = G_sigma_ginv [, y1_index]
  if(rescaleZs==T)  {g_tilda_Y1 = apply(2 * g_tilda_Y1 - 1, 2, function(x) sum(x * g_multiplier)) }
  
  remove(G_sigma_ginv)
  g_tilda_Y1
}

## function to create pairwise interactions  
create_int = function(dat) {
  nobs = dim(dat)[1]
  ncol = dim(dat)[2]
  dat_int = matrix(numeric(), nrow=nobs, ncol=ncol*(ncol-1)/2+ncol)  # change ncol here
  col = 1
  for(k in 1:ncol) {
    for(t in k:ncol) {
      dat_int[,col] = dat[,k] * dat[,t]
      col = col + 1
    }
  }
  dat_int
}

## function to predict Z (prob of being 1) based on Y 
predict_z = function(b, w, y) {  # y each line as 1 obs
  z_p = exp(b + y %*% w) / (1 + exp(b + y %*% w) )
  z_p
}
#------------------------------


# set up the simulation settings  
ny = 9
nz = 3
a_true = c(-1, -0.5, -1, -1,-2, -0.5,-1, -1,-1, 0.5,-0.5,-0.5, -0.5,-1,-2, -0.5,-1,-2)
b_true = c(-1, 1, -0.5)
w_true = matrix(c(1,0,-1,  1.2,0.5,0, 0.8,1,-0.5, 0.3,0,-0.5, 0.1,0.3,0, 0.2,0,-1.2, 0.2,-0.1,-1, 0.5,-0.2,0, 0.3,0,-1.5), ncol=3, nrow = 9, byrow=T)
fix_w = c(2,3,1)
fix_wsign = c(1,1,-1)
theta0_true = c(-2, -2.5, -1)
theta1_true = matrix(c(1,0.5,0,  0.5,1.5,-1,  0,-1,1), byrow=T, nrow=nz) # on Z0 
theta2_true = matrix(c(-0.5,1,  -1,-0.5, -0.5,0.5), nrow=nz)
theta3_true = c(1, 0.5, 0.5)
theta4_true = matrix(c(-1,0.5,0,  -0.5,-1,0,  -0.2,0,0.5), byrow=T, nrow=nz) # on Z0, no Z0sum
theta5_true = matrix(c(0.5,-0.5,  -1,0,  0,-0.5), byrow=T, nrow=nz)  # on X, 

# set up the parameters  
n_cat = c(2, 2, 2, 3, 3, 3, 4, 4, 4)
maxcat = n_cat - 1
max_cat = n_cat - 1
n_a = sum(n_cat) - ny    # number of parameters in a: no param for yj=0
n_ys = prod(n_cat)  # number of possible y values
n_zs = 2^nz         # number of possible z values
ys = expand.grid(0:(n_cat[1]-1), 0:(n_cat[2]-1), 0:(n_cat[3]-1), 0:(n_cat[4]-1), 0:(n_cat[5]-1), 0:(n_cat[6]-1), 0:(n_cat[7]-1), 0:(n_cat[8]-1), 0:(n_cat[9]-1))  
zs = expand.grid(0:1, 0:1, 0:1)
zs0zs1_matrix_cat = matrix(rep(0:1, nz*2), ncol=nz*2)
zs0zs1 = na.omit( expand.grid(data.frame(zs0zs1_matrix_cat)) )

# simulate TEST SET with sample size 10000
n_test = 10000
dat_t = sim_graph2 (seed=1, n_tot=n_test, nsim=270000, start=10000, m=25,  ny=9, nz=nz, nx=2, max_cat=maxcat, a=a_true, b=b_true, w=w_true, theta0=theta0_true, theta1=theta1_true, theta2=theta2_true, theta3=theta3_true, theta4=theta4_true, theta5=theta5_true, rbmfun=gib_rbm_general2, pred_yfun=predict_y3)

H_y_t = cbind(rep(1,dat_t$nobs), dat_t$X, dat_t$y0_obs)   
y0_int_t = create_int(dat_t$y0_obs)
H_y_int_t = cbind(rep(1,dat_t$nobs), dat_t$X, dat_t$y0_obs, y0_int_t)  


vz = matrix(NA, nrow=1000, ncol=13)
correct = matrix(NA, nrow=1000, ncol=13)  

loss = "logit"
n_tot = 400    # 800, 1200
batchsize = 80 # 160, 240

# loss = "hinge"
# c = 2^c(-12:4)
 

for(q in 1: 100) {

  dat = sim_graph2 (seed=q, n_tot=n_tot, nsim=31000, start=10000, m=25,  ny=9, nz=nz, nx=2, max_cat=max_cat, a=a_true, b=b_true, w=w_true, theta0=theta0_true, theta1=theta1_true, theta2=theta2_true, theta3=theta3_true, theta4=theta4_true, theta5=theta5_true, rbmfun=gib_rbm_general2, pred_yfun=predict_y3)  
  
  nobs=n_tot
  
  a = runif(n_a, -0.1, 0.1)
  b = runif(nz, min = -0.1, max = 0.1)
  w = matrix(runif(ny*nz, min = -0.1, max = 0.1), nrow = ny)
  for(k in 1:nz) {
    w[fix_w[k], k] = fix_wsign[k]  # set those starting values to be 1/-1
  }

  rbm = run_rbms(lrates=c(rep(2), rep(1,2)), niters=c(20, rep(20,20)), ny=ny, nz=nz,  n_cat=n_cat, a0=a, b0=b, w0=w,  obs_y0=dat$y0_obs, lrate_scale = 1, momentum=0.5, batchsize=batchsize, outlast=F)  
  
  gtilda_Y1 = cal_g_tilda_Y1 (rbm$a, rbm$w,  ny=9, nz=nz, n_cat, dat$y1_obs) 
  
  ##-----  O-learning with ysum as outcome and y0 items in H
  H_y = cbind(rep(1,nobs), dat$X, dat$y0_obs)   
  y0_int = create_int(dat$y0_obs)
  H_y_int = cbind(rep(1,dat$nobs), dat$X, dat$y0_obs, y0_int)  
  
  owl_y = owl(H_y, dat$trt, RR=-dat$y1sum, n=nobs, K=1, pi=rep(0.5,nobs), augment=F, c=c, loss=loss)
  
  predz_y = predict(owl_y, H=H_y_t, AA=dat_t$trt, RR=dat_t$z1sum, K=1, pi=rep(0.5, 10000))
  vz[q,1] = predz_y$valuefun
  correct[q, 1] = 1 - sum(predz_y$treatment[[1]] * dat_t$opt_trt2==-1) / 10000
  
  owl_y_int = owl(H_y_int, dat$trt, RR=-dat$y1sum, n=nobs, K=1, pi=rep(0.5,nobs), augment=F, c=c, loss=loss)
  predy_y_int = predict(owl_y_int, H=H_y_int_t, AA=dat_t$trt, RR=dat_t$y1sum, K=1, pi=rep(0.5, 10000))
  predz_y_int = predict(owl_y_int, H=H_y_int_t, AA=dat_t$trt, RR=dat_t$z1sum, K=1, pi=rep(0.5, 10000))
  vz[q,11] = predz_y_int$valuefun
  correct[q, 11] = 1 - sum(predz_y_int$treatment[[1]] * dat_t$opt_trt2==-1) / 10000
  
  ##-----  O-learning with gtilda_Y1 as outcome and predicted E[z0|Y0] in H
  pred_z0 = predict_z(rbm$b, rbm$w, dat$y0_obs)       
  H_z = cbind(rep(1,nobs), dat$X, pred_z0)      
  predz0_int = create_int(pred_z0)
  H_z_int = cbind(rep(1,nobs), dat$X, pred_z0, predz0_int)       
  
  pred_z0_t = predict_z(rbm$b, rbm$w, dat_t$y0_obs)
  H_z_t =  cbind(rep(1,n_test), dat_t$X, pred_z0_t)
  predz0_int_t = create_int(pred_z0_t)
  H_z_int_t =  cbind(rep(1,n_test), dat_t$X, pred_z0_t, predz0_int_t)
  
  owl_z = owl(H_z, dat$trt, RR = -gtilda_Y1, n=nobs, K=1, pi=rep(0.5,nobs), augment=F, c=c, loss=loss)
  owl_z_int = owl(H_z_int, dat$trt, RR = -gtilda_Y1, n=nobs, K=1, pi=rep(0.5,nobs), augment=F, c=c, loss=loss)
  
  predz_z = predict(owl_z, H=H_z_t, AA=dat_t$trt, RR=dat_t$z1sum, K=1, pi=rep(0.5, 10000))
  vz[q,2] = predz_z$valuefun
  correct[q,2] = 1 - sum(predz_z$treatment[[1]] * dat_t$opt_trt2==-1) / 10000
  
  predz_z_int = predict(owl_z_int, H=H_z_int_t, AA=dat_t$trt, RR=dat_t$z1sum, K=1, pi=rep(0.5, 10000))
  vz[q,12] = predz_z_int$valuefun
  correct[q, 12] = 1 - sum(predz_z_int$treatment[[1]] * dat_t$opt_trt2==-1) / 10000
  
  ##---------  O-learning with true gtilda(Y1) as outcome and true E[z0|Y0] in H
  gtilda_Y1_true = cal_g_tilda_Y1 (a_true, w_true,  ny=9, nz=nz, n_cat, dat$y1_obs)
  
  H_z_true2 = cbind(rep(1,nobs), dat$X, dat$z0_p)
  H_z_t_true2 =  cbind(rep(1,n_test), dat_t$X, dat_t$z0_p)  
  z0_int2 = create_int(dat$z0_p)
  z0_int2_t = create_int(dat_t$z0_p)
  H_z_true2_int = cbind(rep(1,nobs), dat$X, dat$z0_p, z0_int2, rowSums(z0_int2))
  H_z_t_true2_int =  cbind(rep(1,n_test), dat_t$X, dat_t$z0_p, z0_int2_t, rowSums(z0_int2_t))
  
  owl_z_true3 = owl(H_z_true2_int, dat$trt, RR = -gtilda_Y1_true, n=nobs, K=1, pi=rep(0.5,nobs), augment=F, c=c, loss=loss)  # note use negative dat$z1sum
  predz_z_true3 = predict(owl_z_true3, H=H_z_t_true2_int, AA=dat_t$trt, RR=dat_t$z1sum, K=1, pi=rep(0.5, 10000))
  vz[q, 3] = predz_z_true3$valuefun
  correct[q, 3] = 1 - sum(predz_z_true3$treatment[[1]] * dat_t$opt_trt2==-1) / 10000
  
  owl_z_true3_noint = owl(H_z_true2, dat$trt, RR = -gtilda_Y1_true, n=nobs, K=1, pi=rep(0.5,nobs), augment=F, c=c, loss=loss)  # note use negative dat$z1sum
  predz_z_true3_noint = predict(owl_z_true3_noint, H=H_z_t_true2, AA=dat_t$trt, RR=dat_t$z1sum, K=1, pi=rep(0.5, 10000))
  vz[q, 13] = predz_z_true3_noint$valuefun
  correct[q, 13] = 1 - sum(predz_z_true3_noint$treatment[[1]] * dat_t$opt_trt2==-1) / 10000

  print(c("q", vz[q,]))
}


 

