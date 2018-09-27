# FastLZeroSpikeInference: A package for estimating spike times from calcium imaging data using an L0 penalty 

This package implements an algorithm for deconvolving calcium imaging data
for a single neuron in order to estimate the times at which the neuron
spikes.

This algorithm solves the optimization problems

### AR(1) model

<img src="un-constr.png" alt="alt text" width="600" height="80">

for the global optimum, where y_t is the observed fluorescence at the tth timepoint.

### Constrained AR(1) model

<img src="constr.png" alt="alt text" width="600" height="140">

for the global optimum, where y_t is the observed fluorescence at the tth timepoint.

We introduce the constant EPS > 0, typically on the order of 10^-10, to avoid 
arbitrarily small calcium concentrations that would result in numerical  
instabilities. In practice, this means that the estimated calcium concentration 
decays according to the AR(1) model for values greater than EPS and is equal to EPS thereafter.

When estimating the spikes, it is not necessary to explicitly compute the 
calcium concentration. Therefore, if only the spike times are required, the user
can avoid this computation cost by setting the compute_fitted_values boolean to false. 
By default, the calcium concentration is not estimated. 

Given the set of estimated spikes produced from the estimate_spike, the calcium concentration
can be estimated with the estimate_calcium function.

For additional information see: 
 
Jewell, Hocking, Fearnhead, and Witten (2018). [Fast Nonconvex Deconvolution of Calcium Imaging Data](https://arxiv.org/abs/1802.07380)

Jewell and Witten (2017). [Exact Spike Train Inference Via L0 Optimization](https://arxiv.org/abs/1703.08644)

R examples 
```r
sim <- simulate_ar1(n = 500, gam = 0.95, poisMean = 0.009, sd = 0.05, seed = 1)
plot(sim)
 
## Fits for a single tuning parameter

# AR(1) model
fit <- estimate_spikes(dat = sim$fl, gam = 0.95, lambda = 1)
print(fit)

# compute fitted values from prev. fit
fit <- estimate_calcium(fit)
plot(fit)

# or
fit <- estimate_spikes(dat = sim$fl, gam = 0.95, lambda = 1, estimate_calcium = T)
plot(fit)
 
# Constrained AR(1) model
fit <- estimate_spikes(dat = sim$fl, gam = 0.95, lambda = 1, constraint = T, estimate_calcium = T)
print(fit)
plot(fit)
```

Install 
-----

If ``devtools`` is installed type 

```r
devtools::install_github("jewellsean/FastLZeroSpikeInference")
```

Usage
----

Once installed type 
```{r}
library(FastLZeroSpikeInference)
```


Alpha Python Instructions
---

Below are instructions for an alpha-release of a simple c-types python wrapper for this C++ code. Only Unix type systems are currently supported. 

NB: ubuntu requires `sudo apt install clang g++` 

Within terminal, clone this repo and run the make script: 

```
git clone "https://github.com/jewellsean/FastLZeroSpikeInference.git"
cd FastLZeroSpikeInference/python
./make.sh
```

An example using this code can be viewed [here](https://github.com/jewellsean/FastLZeroSpikeInference/blob/master/examples/python/simple_example.py).

References
-----
Jewell, Hocking, Fearnhead, and Witten (2018). [Fast Nonconvex Deconvolution of Calcium Imaging Data](https://arxiv.org/abs/1802.07380)

Jewell and Witten (2017). [Exact Spike Train Inference Via L0 Optimization](https://arxiv.org/abs/1703.08644)

Hocking, T. D., Rigaill, G., Fearnhead, P., & Bourque, G. (2017). [A log-linear time algorithm for constrained changepoint detection](https://arxiv.org/abs/1703.03352)
