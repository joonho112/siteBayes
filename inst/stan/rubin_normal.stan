functions {
  //----------------------------------------------------------------
  // 1) Helper function for scalar-valued priors (tau, sigma_tau, etc.)
  //----------------------------------------------------------------
  real realprior_lpdf(real x, int family, vector pars) {
    /*
      family-coded logic:
        0 => uniform
        1 => normal
        2 => cauchy
        3 => student_t
        5 => lognormal
        6 => half_student_t
        7 => half_cauchy
        9 => (placeholder for multinormal)
        10 => (placeholder for lkj)
    */
    if (family == 1) {
      // e.g. pars = [loc, scale, 0] => Normal(loc, scale)
      return normal_lpdf(x | pars[1], pars[2]);

    } else if (family == 2) {
      // e.g. pars = [location, scale, 0] => Cauchy(loc, scale)
      return cauchy_lpdf(x | pars[1], pars[2]);

    } else if (family == 3) {
      // e.g. pars = [df, loc, scale] => Student_t(df, loc, scale)
      return student_t_lpdf(x | pars[1], pars[2], pars[3]);

    } else if (family == 5) {
      // e.g. pars = [mu, sigma, 0] => Lognormal(mu, sigma)
      return lognormal_lpdf(x | pars[1], pars[2]);

    } else if (family == 6) {
      // half_student_t => T[0,∞). Typically declared real<lower=0>.
      // Student_t with location=0, scale=pars[2], df=pars[1].
      return student_t_lpdf(x | pars[1], 0, pars[2]);

    } else if (family == 7) {
      // half_cauchy => T[0,∞).
      // cauchy with loc=pars[1], scale=pars[2], restricted to x>0 in parameters block
      return cauchy_lpdf(x | pars[1], pars[2]);

    } else if (family == 0) {
      // uniform => e.g. pars=[lower, upper, 0]
      return uniform_lpdf(x | pars[1], pars[2]);

    } else {
      // fallback or unrecognized => treat as wide uniform
      return uniform_lpdf(x | -1.0e6, 1.0e6);
    }
  }

  //----------------------------------------------------------------
  // 2) Helper function for vector-valued priors (beta, if any)
  //----------------------------------------------------------------
  real vecprior_lpdf(vector v, int family, vector pars) {
    /*
      Typically 1 => normal, 3 => student_t, 0 => uniform, etc.
      If the user sets others, you can expand the logic accordingly.
    */
    if (family == 1) {
      // e.g. pars = [loc, scale, 0]
      return normal_lpdf(v | pars[1], pars[2]);

    } else if (family == 3) {
      // e.g. pars = [df, loc, scale]
      return student_t_lpdf(v | pars[1], pars[2], pars[3]);

    } else if (family == 0) {
      // vector uniform => typically not standard in Stan
      // fallback to a wide normal
      return normal_lpdf(v | 0, 1.0e4);

    } else {
      // fallback
      return normal_lpdf(v | 0, 1.0e4);
    }
  }
}

data {
  // 1) Basic meta-analysis data
  int<lower=0> J;                // Number of sites
  vector[J]     tau_hat_j;       // Observed site-level effects
  vector<lower=0>[J] se_j;       // Standard errors for each site

  // 2) Covariates (for fixed effects)
  int<lower=0> Nc;               // Number of covariates
  matrix[J, Nc] X;               // Covariate matrix

  // 3) Priors, encoded
  int prior_hypermean_fam;       // (for the overall effect 'tau')
  int prior_hypersd_fam;         // (for the between-site SD 'sigma_tau')
  int prior_beta_fam;
  vector[3] prior_hypermean_val;
  vector[3] prior_hypersd_val;
  vector[3] prior_beta_val;
}

transformed data {
  // No special transformations here
}

parameters {
  real tau;                       // overall average treatment effect
  real<lower=0> sigma_tau;        // between-site SD
  vector[J] eta;                  // non-centered random effects
  vector[Nc] beta;                // fixed-effect coefficients
}

transformed parameters {
  // The "true" site-level parameters
  vector[J] tau_j;

  // The fixed-effects contribution
  vector[J] fe_j;

  // 1) partial-pooling site effects
  for (j in 1:J) {
    tau_j[j] = tau + sigma_tau * eta[j];
  }

  // 2) if we have covariates, fe_j = X * beta, else 0
  if (Nc > 0) {
    fe_j = X * beta;
  } else {
    fe_j = rep_vector(0.0, J);
  }
}

model {
  // 1) Priors on beta
  target += vecprior_lpdf(beta | prior_beta_fam, prior_beta_val);

  // 2) Priors on tau (overall effect) and sigma_tau (between-site SD)
  target += realprior_lpdf(tau       | prior_hypermean_fam, prior_hypermean_val);
  target += realprior_lpdf(sigma_tau | prior_hypersd_fam,   prior_hypersd_val);

  // 3) Non-centered random effects
  eta ~ normal(0, 1);

  // 4) Likelihood: measure tau_hat_j around (tau_j + fe_j) with se_j
  for (j in 1:J) {
    tau_hat_j[j] ~ normal(tau_j[j] + fe_j[j], se_j[j]);
  }
}

generated quantities {
  // 1) per-site log-likelihood
  vector[J] log_lik;
  // 2) posterior predictions
  vector[J] posterior_pred;
  // 3) shrinkage factor for each site
  vector[J] shrinkage;
  // 4) geometric mean-based measure
  real log_se_squared_mean = 0;
  real I_level;

  // fill them in
  for(j in 1:J) {
    // log-likelihood contribution
    log_lik[j] = normal_lpdf(tau_hat_j[j] | tau_j[j] + fe_j[j], se_j[j]);

    // posterior predictive
    posterior_pred[j] = normal_rng(tau_j[j] + fe_j[j], se_j[j]);

    // shrinkage factor: ratio of sigma_tau^2 to (sigma_tau^2 + se_j[j]^2)
    shrinkage[j] = sigma_tau^2 / (sigma_tau^2 + se_j[j]^2);

    // accumulate sum for geometric mean of se_j[j]^2
    log_se_squared_mean += log(se_j[j]^2);
  }
  // finalize geometric mean-based measure
  log_se_squared_mean /= J;
  I_level = sigma_tau^2 / (sigma_tau^2 + exp(log_se_squared_mean));
}
