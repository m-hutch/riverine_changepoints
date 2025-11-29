#include <Rcpp.h>
using namespace Rcpp;

// adaptive hazard function: can be constant or depend on run length
inline double hazard(int run_length, double lambda, bool adaptive) {
  // constant
  if (!adaptive) return 1.0 / lambda;
  // adaptive hazard: increases slightly with run length
  return std::min(1.0, (1.0 / lambda) * (1.0 + 0.01 * run_length));
}

// [[Rcpp::export]]
List bocpd_pearsonIII(NumericVector data,
                      double prior_shape,
                      double prior_scale,
                      double location,
                      double hazard_lambda,
                      bool adaptive_hazard = false) {
  
  int T = data.size();
  NumericMatrix run_length_probs(T + 1, T + 1);
  run_length_probs(0, 0) = 1.0;
  
  // Track shape and scale for each run length
  std::vector<double> shape(T + 1, prior_shape);
  std::vector<double> scale(T + 1, prior_scale);
  
  NumericVector changepoint_probs(T);
  
  for (int t = 0; t < T; ++t) {
    NumericVector new_rl_probs(T + 1, 0.0);
    
    for (int rl = 0; rl <= t; ++rl) {
      double shifted_x = data[t] - location;
      if (shifted_x <= 0) shifted_x = 1e-8; // avoid invalid values
      
      double pred_shape = shape[rl];
      double pred_scale = scale[rl];
      
      // Pearson Type III uses Gamma PDF on shifted data
      double pred_prob = R::dgamma(shifted_x, pred_shape, pred_scale, false);
      
      // Growth: no changepoint
      new_rl_probs[rl + 1] += run_length_probs(rl, t) * pred_prob * (1.0 - hazard(rl, hazard_lambda, adaptive_hazard));
    }
    
    // Changepoint probability
    double cp_prob = 0.0;
    for (int rl = 0; rl <= t; ++rl) {
      double shifted_x = data[t] - location;
      if (shifted_x <= 0) shifted_x = 1e-8;
      
      double pred_shape = shape[rl];
      double pred_scale = scale[rl];
      double pred_prob = R::dgamma(shifted_x, pred_shape, pred_scale, false);
      
      cp_prob += run_length_probs(rl, t) * pred_prob * hazard(rl, hazard_lambda, adaptive_hazard);
    }
    new_rl_probs[0] = cp_prob;
    
    // Normalize
    double total = std::accumulate(new_rl_probs.begin(), new_rl_probs.end(), 0.0);
    for (int rl = 0; rl <= t + 1; ++rl) {
      run_length_probs(rl, t + 1) = new_rl_probs[rl] / total;
    }
    
    // Update sufficient statistics for Pearson III
    for (int rl = 0; rl <= t + 1; ++rl) {
      if (rl == 0) {
        shape[rl] = prior_shape + 1;
        scale[rl] = prior_scale; // scale remains constant
      } else {
        shape[rl] = shape[rl - 1] + 1;
        // FUTURE WORK: update scale adaptively
      }
    }
    
    changepoint_probs[t] = run_length_probs(0, t + 1);
  }
  
  return List::create(
    Named("run_length_probs") = run_length_probs,
    Named("changepoint_probs") = changepoint_probs
  );
}