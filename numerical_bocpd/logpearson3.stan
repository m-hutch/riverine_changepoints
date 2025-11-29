
functions {
  
 real logpearson3_lpdf(real y, real shape, real scale, real location) {
    if (y < location) {
      return negative_infinity(); // invalid region
    }
    real x = y - location; // shifted
    return gamma_lpdf(x | shape, 1/scale);
  }

}

data {
  int<lower=1> N;
  array[N] real<lower=0> y;
  real<lower=0> shape_prior;
  real<lower=0> scale_prior;
  real location_prior;
}

parameters {
  real<lower=0> shape;
  real<lower=0> scale;
  real location;
}
model {
  // Priors
  shape ~ lognormal(shape_prior, 2);
  scale ~ lognormal(scale_prior, 2);
  location ~ normal(location_prior, 2);

  // Likelihood
  for (n in 1:N) {
    target += logpearson3_lpdf(y[n] | shape, scale, location);
  }
}
