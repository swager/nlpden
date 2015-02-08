# nlpden
Density estimation by quadratic programming, as proposed by Wager (2014).

To install this package in R, run the following commands:

```R
install.packages("devtools")
library(devtools) 
install_github("swager/nlpden")
```

Example usage:

```R
library(nlpden)

# Make some data...
mu = c(runif(N/2, min = -3, max = 3), rep(0, N/2))
X = mu + rnorm(N)
  
#  Run nlpden
f.hat = nlpden(X)
plot(f.hat, type = "l")
```

#### References
Stefan Wager. <b>A Geometric Approach to Density Estimation with Additive Noise</b>. <i>Statistica Sinica</i>, 24(2), 2014. 
