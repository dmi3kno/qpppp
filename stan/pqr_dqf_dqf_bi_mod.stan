// begin file stan/gexp_pdf_pdf_mod.stan
functions{
  // standard quantile function for the error term
real fsld_s_sf(real u, real eta, real k, real dlt){
    real Qslogis = (1-dlt)*log(u)-dlt*log1m(u);
    real Qslogis05 = (1-dlt)*log(0.5)-dlt*log1m(0.5);
    real Qunif = k*u;
    real Qunif05 = k*0.5;
    return eta*(Qslogis+Qunif-Qslogis05-Qunif05);
  }
  // quantile function for prior
  real fsld_s_qf(real u, real a, real eta, real k, real dlt){
    real Qslogis = (1-dlt)*log(u)-dlt*log1m(u);
    real Qunif = k*u;
    return a + eta*(Qslogis+Qunif);
  }
  // quantile function for prior
  real fld_s_qf(real u, real a, real eta, real k){
    // this is centered with location 0
    real Qlogis = logit(u);
    real Qunif = k*u;
    return a + eta*(Qlogis+Qunif);
  }
// this is regression quantile funciton
 real rqf_s_qf(real u, real alph, real bt, real eta, real k, real dlt, real x){
    if(is_nan(u)) reject("p can not be nan!");
    if(u>=1) return positive_infinity();
    if(u<=0) return negative_infinity();
   real S =  fsld_s_sf(u, eta, k, dlt); // location is zero
   return alph+bt*sqrt(x)+sqrt(x)*S;
 }

// this is quantile density for the error term DQF FSLD
real fsld_s_ldqf_lpdf(real u, real eta, real k, real dlt){
  return -log(eta*((1-dlt)/u+dlt/(1-u)+k));
}

// This function is the algebra system with a signature:
// vector algebra_system (vector y, vector theta, data vector x)
real rootfun(real u0, vector params, data real y_r){
  return y_r - rqf_s_qf(u0, params[1], params[2], params[3], params[4], params[5], params[6]);
}

real iqf_brent01(vector params, data real y_r, data real rel_tol, data real max_steps, int verbose){
    real u0=rel_tol;
    real u1=1-rel_tol;
    int steps_taken=0;
    real f0=rootfun(u0, params, y_r);
    real f1=rootfun(u1, params, y_r);
    real tmp=0;
    // assuming rootfun is non-increasing
    // if target function at both ends is positive, the true root is above
    // if target function at both ends is negative, the true root is below
    if(f1>0 && f0>0) return(1-rel_tol);
    if(f1<0 && f0<0) return(rel_tol);

   if(is_inf(f0) || is_inf(f1)) reject("Tolerance is too small!");

    // swap for non-decreasing function
    if (abs(f0)<abs(f1)){
      tmp=u0;
      u0=u1;
      u1=tmp;
      tmp=f0;
      f0=f1;
      f1=tmp;
    }
    real u2=u0;
    real f2=f0;
    int mflag=1;
    real n=0;
    real d=0;

    while(steps_taken<max_steps && abs(u1-u0)>rel_tol){
      f0=rootfun(u0, params, y_r);
      f1=rootfun(u1, params, y_r);
      f2=rootfun(u2, params, y_r);

      if(f0!=f2 && f1 !=f2){
        real l0=(u0*f1*f2)/((f0-f1)*(f0-f2));
        real l1=(u1*f0*f2)/((f1-f0)*(f1-f2));
        real l2=(u2*f1*f0)/((f2-f0)*(f2-f1));
        n = l0+l1+l2;
      } else {
        n = u1-(f1*(u1-u0)/(f1-f0));
      }

      if((n<(3*u0+u1)/4 || n>u1) ||
      (mflag==1 && abs(n-u1)>=abs(u1-u2)/2.0) ||
      (mflag==0 && abs(n-u1)>=abs(u2-d)/2.0) ||
      (mflag==1 && abs(u1-u2)<rel_tol ) ||
      (mflag==0 && abs(u2-d)<rel_tol)
      ){
        n=(u0+u1)/2.0;
        mflag=1;
      } else {
        mflag=0;
      }
      real fn=rootfun(n, params, y_r);
      d=u2;
      u2=u1;
      if(f0*fn<0){
        u1=n;
      } else {
        u0=n;
      }
      if(abs(f0)<abs(f1)){
        tmp=u0;
        u0=u1;
        u1=tmp;
      }
      steps_taken+=1;
    }
   if(verbose) print("Steps taken ", steps_taken);
    return u1;

  }

}// end of functions block

data {
  int<lower=0> N; // number of data points
  vector[N] y; // data
  vector[N] x; // covariate
  real intrc_a;
  real intrc_eta;
  real intrc_k; // hyperparameter alpha
  real slope_a;
  real slope_eta;
  real slope_k;
  real slope_dlt;
  real rel_tol;
  real f_tol;
  int max_steps;
  int verbose;
}
parameters {
  //real alph; // strictly positive
  //real bt; // strictly positive
  real<lower=0, upper=1> v; // for indirect intercept
  real<lower=0, upper=1> w; // for indirect slope
  real<lower=0> eta; // strictly positive scale
  real<lower=0> k; // strictly positive flatness
  real<lower=0, upper=1> dlt; // strictly positive
}
transformed parameters{
  real alph = fld_s_qf(v, intrc_a, intrc_eta, intrc_k);
  real bt = fsld_s_qf(w, slope_a, slope_eta, slope_k, slope_dlt);
}

model {
  vector[N] u;
  target += exponential_lpdf(eta | 0.1);
  target += exponential_lpdf(k | 1);
  target += beta_lpdf(dlt | 2,1);
  for (i in 1:N){
   u[i] = iqf_brent01(to_vector({alph, bt, eta, k, dlt, x[i]}), y[i], rel_tol, max_steps, verbose);
   target += fsld_s_ldqf_lpdf(u[i]|eta,k, dlt)-log(sqrt(x[i]));
  }
}
// end file stan/gexp_pdf_pdf_mod.stan
