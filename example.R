
# example code to run the algorithm
source("functions.R")

# set up the simulation parameters  
ny = 9    # number of observed items
nz = 2    # number of latent Z's
nx = 2    # number of X

n_cat = c(2, 2, 2, 3, 3, 3, 4, 4, 4)  # number of categories for each Y
maxcat = n_cat - 1

# a_true = c(-1, -0.5, -1, -1,-2, -0.5,-1, -1,-1, 0.5,-0.5,-0.5, -0.5,-1,-2, -0.5,-1,-2)  # number of parameters in a: no param for yj=0
# b_true = c(-1, 1, -0.5)
# w_true = matrix(c(1,0,-1,  1.2,0.5,0, 0.8,1,-0.5, 0.3,0,-0.5, 0.1,0.3,0, 0.2,0,-1.2, 0.2,-0.1,-1, 0.5,-0.2,0, 0.3,0,-1.5), ncol=3, nrow = 9, byrow=T)
fix_w = c(1,2)         # to control the direction of one entry of W for each latent Z
fix_wsign = c(3,3)    # signs of those entries

# theta0_true = c(-2, -2.5, -1)
# theta1_true = matrix(c(1,0.5,0,  0.5,1.5,-1,  0,-1,1), byrow=T, nrow=nz) # on Z0 
# theta2_true = matrix(c(-0.5,1,  -1,-0.5, -0.5,0.5), nrow=nz)
# theta3_true = c(1, 0.5, 0.5)
# theta4_true = matrix(c(-1,0.5,0,  -0.5,-1,0,  -0.2,0,0.5), byrow=T, nrow=nz) # on Z0, no Z0sum
# theta5_true = matrix(c(0.5,-0.5,  -1,0,  0,-0.5), byrow=T, nrow=nz)  # on X, 

theta0_true = c(-2, -2.5)
theta1_true = matrix(c(1,0.5,  0.5,1.5), byrow=T, nrow=nz) # on Z0
theta2_true = matrix(c(-0.5,1,  -1,-0.5), nrow=nz) # on X
theta3_true = c(1, 0.5)
theta4_true = matrix(c(-1,0.5, -0.5,-1), byrow=T, nrow=nz) # on Z0 A, no Z0sum
theta5_true = matrix(c(0.5,-0.5,  -1,0), byrow=T, nrow=nz)  # on X A,


# try new RBM setting
a_true = c(-1, -0.5, -1, -1,-2, -0.5,-1, -1,-1, 0.5,-0.5,-0.5, -0.5,-1,-2, -0.5,-1,-2)  # number of parameters in a: no param for yj=0
b_true = c(-2, 0)
w_true = matrix(c(1.5,-0.5,  0,1.5, -0.5,0.5, 1.5,0.5, -1,0.5, 2,-1, 0.5,-0.5, 0.5,0, -1,0.5), ncol=2, nrow = 9, byrow=T)

fix_w = c(4,2)         # to control the direction of one entry of W for each latent Z
fix_wsign = c(3,3)    # signs of those entries


# create training data with sample size of 400, with RBMs simulated from MCMC wtih burn-in of 10000, take one obs every 25 iterations
nobs = 400
dat = sim_graph2 (seed=1, n_tot=nobs, nsim=20000, start=10000, m=25,  ny=ny, nz=nz, nx=nx, max_cat=maxcat, 
                  a=a_true, b=b_true, w=w_true, 
                  theta0=theta0_true, theta1=theta1_true, theta2=theta2_true, 
                  theta3=theta3_true, theta4=theta4_true, theta5=theta5_true)
  table(dat$opt_trt)
  table(dat$opt_trt2)
  
  apply(dat$z0==1, 2, sum)


# initialization parameters for the RBM model
a = runif(length(a_true), -0.01, 0.01)
b = runif(nz, min = -0.01, max = 0.01)
w = matrix(runif(ny*nz, min = -0.01, max = 0.01), nrow = ny)
for(k in 1:nz) {
  w[fix_w[k], k] = fix_wsign[k]  # set those starting values to be 1/-1
}

# fit the RBM
rbm = run_rbms(lrates=c(rep(1, 2), rep(0.5,4)), niters=rep(50, 6), ny=ny, nz=nz,  n_cat=n_cat, 
               a0=a, b0=b, w0=w,  obs_y0=dat$y0_obs, 
               lrate_scale = 1, momentum=0.5, batchsize=50, fix_w=fix_w, fix_wsign=fix_wsign,
               MCMC=F, cd=1, cal_likelihood=F, outlast=F)
    rbm$a
    a_true
    rbm$w
    w_true
    apply(dat$z0, 2, sum)

# construct surrogate predictors based on fitted RBM
pred_z0 = predict_z(rbm$b, rbm$w, dat$y0_obs)  
    apply((pred_z0 > 0.5) != dat$z0, 2, sum)
H_z = cbind(dat$X, pred_z0)
    str(H_z)
    str(dat$X)

# construct the surrogate outcome based on the fitted RBM
gtilda_Y1 = cal_g_tilda_Y1 (rbm$a, rbm$w,  ny=9, nz=nz, n_cat, dat$y1_obs) 


# run AOL with surrogate outcome and predictor (minimize the surrogate outcome), using the "DTRlearn2" package
owl_z = owl(H_z, dat$trt, RR = -gtilda_Y1, n=nobs, K=1, pi=rep(0.5,nobs), augment=F, loss="logit")




