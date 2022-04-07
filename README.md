# Learning Optimal Individualized Treatment Rules (ITRs) for Latent Outcomes

For many mental disorders, latent mental status inferred from multiple-domain psychological or clinical symptoms may perform as a better characterization of the underlying disorder status than a simple summary score of the symptoms, and they may also serve as more reliable and representative features to differentiate treatment responses.

We provide a new paradigm for learning optimal individualized treatment rules (ITRs) by modeling patients' latent mental status.  We first learn the multi-domain latent states at baseline from the observed symptoms under a restricted Boltzmann machine (RBM) model. We then optimize a value function defined by the latent states after treatment by exploiting a transformation of the observed symptoms based on the RBM without modeling the relationship between the latent mental states before and after treatment.


### Reference
**Yuan Chen**, Donglin Zeng, Yuanjia Wang (2021). Learning Individualized Treatment Rules for Multiple-Domain Latent Outcomes. *Journal of the American Statistical Association*, 116(533), 269-282.


### Code 
In ```functions.R```, we provide the code to simulate data from a RBM model and code to simulate training / test data with observed individual item scores before and after treatment, as well code for implementing the RBM model. Two algorithms are implemented for fitting the RBM model: 1) to maximize the exact likeilhood using gradient ascent, 2) using MCMC approximation for the gradients with the contrastive divergence algorithm. 

An example is provided in ```example.R``` to illustrate the use of the functions, including simulating data, fitting RBM model, constructing surrogate outcome and predictors and learnining the optimal ITRs for the latent outcomes.

