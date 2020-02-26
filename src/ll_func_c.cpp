//
//  logl.h
//  cal
//
//  Created by David Hodgson on 03/02/2020.
//  Copyright © 2020 David Hodgson. All rights reserved.
//

#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
#include <boost/random.hpp>
#include <boost/math/distributions.hpp>
#include "ascent/Ascent.h"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins("cpp14")]]

using namespace Rcpp;
using namespace std;
using namespace Eigen;
using namespace asc;

typedef vector< double > num_vec;       //General-purpose numerical vector

std::random_device dev;
std::mt19937 engine(dev());
typedef boost::mt19937 PRNG_s;
PRNG_s rng(engine()); //Generate non-static random numbers (pick different numbers from prior distribution each run)

double logistic_dist(double a, double b, char r, double x = 0)
{
  if (r =='r')
  {
    boost::random::uniform_real_distribution<> u(0,1);
    double t = u(rng);
    return a + b*log(t/(1-t));
    
  }
  else if ( r =='m')
  {boost::math::logistic_distribution<> u(a,b); return mean(u);}
  else
  {boost::math::logistic_distribution<> u(a,b); return pdf(u,x);}
}

double poisson_cdf(double l, double a, double x)
{
  if( l == 0.0 || a == 0.0)
  {
    boost::math::poisson_distribution<> p(0.000001); return cdf(p,x);
  }
  else
  {
    boost::math::poisson_distribution<> p(l*a); return cdf(p,x);
  }
}

long double stirl(double n)
{
  if (n == 0)
    return 0;
  else {
    double x = n + 1;
    return (x - 0.5)*log(x) - x + 0.5*log(2*PI) + 1.0/(12*x) - 1.0/(360.0*pow(x,3)); // https://en.wikipedia.org/wiki/Stirling%27s_approximation#Speed_of_convergence_and_error_estimates
  }
}

namespace model
{

num_vec prop_init_ex(double l, double a1, double a2)
{
  num_vec prop(2);
  prop[0] = abs(poisson_cdf(l,a2,0)-poisson_cdf(l,a1,0))/((a2-a1)*l);
  prop[1] = 1 - prop[0];
  return prop;
}

inline double trans_par(double x, double a, double b)
{
  if (x<-100)
    return a;
  else if (x>100)
    return b;
  else
    return a + (b-a)*exp(x)/(exp(x) + 1);
}

//Initial conditions
num_vec init_cond(List dem, NumericVector age_lim, num_vec popsize, NumericVector pars)
{
  int A =  age_lim.size();
    
  num_vec comp_size_init;
  double a1, a2;
  double I1 = trans_par(pars["I1"], 0, 1);
  double I2 = trans_par(pars["I2"], 0, 1);
  double I3 = trans_par(pars["I3"], 0.25, 1);
  double pI1, pI2;
  
  for (int a = 0; a < A ; a++){
    if (a < 1)
    {
      a1 = age_lim(a); a2 =  age_lim(a+1);
    }
    
    else
    {
      a1 =  age_lim(a); a2 = 90;
    }

    num_vec prop_init = prop_init_ex(I3, a1, a2);
    pI1 = prop_init[0]; pI2 = prop_init[1];
    
    comp_size_init.push_back(popsize[a]*pI1*(1-I1*1.0)*(1-I2));  //Sus
    comp_size_init.push_back(popsize[a]*pI1*I1*1.0);    //Inf
    comp_size_init.push_back(popsize[a]*pI1*(1-I1*1.0)*I2);     //Rec
    
    comp_size_init.push_back(popsize[a]*pI2*(1-I1*1.0)*(1-I2));  //Sus
    comp_size_init.push_back(popsize[a]*pI2*I1*1.0);    //Inf
    comp_size_init.push_back(popsize[a]*pI2*(1-I1*1.0)*I2);     //Rec
    
    comp_size_init.push_back(0);     //Rec
  }
  return comp_size_init;
}

struct ODE_dynamics
{
  double t_start, t_burn, t_end, dt;
  bool get_incidence;
  
  ODE_dynamics(const double t_start_i, const double t_burn_i, const double t_end_i, const double dt_i): t_start(t_start_i), t_burn(t_burn_i), t_end(t_end_i) , dt(dt_i) {}
};

class ODE_desc
{
  
public:
  NumericVector pars;
  DataFrame data;
  
  double d1 = trans_par(pars["d1"], 0, 1);
  double ga0 = 1.0/(trans_par(pars["ga0"], 2, 20));
  double ga1 = 1;
  double alpha = trans_par(pars["alpha"], 0, 1);
  double om = 1.0/(trans_par(pars["om"], 60, 365));
  
  double a = trans_par(pars["a"], 0, 1);
  double b = trans_par(pars["b"], 0, 1);
  double psi = trans_par(pars["psi"], 0, 0.1*365);
  double phi = trans_par(pars["phi"], 100, 250);
  double th = 0;
  int A;
  double mu_t; //data here
  num_vec eta;

  NumericMatrix poly;
  NumericVector age_lim;
  num_vec popsize;
  List dem;
        
  ODE_desc(NumericMatrix poly_t, List dem_t, NumericVector age_lim_t, num_vec popsize_t, NumericVector pars_t): poly(poly_t), dem(dem_t), age_lim(age_lim_t), popsize(popsize_t), pars(pars_t) {

      A = age_lim.size();
      mu_t = dem["mu"];
      int tot = dem["pop"];
      
      for (int i = 1; i < A; i++)
      {
          eta.push_back( 1.0/(365.0*(age_lim[i] - age_lim[i-1])) );
      }
      eta.push_back(mu_t/(tot - (mu*365)*age_lim[A-1]) );
  };
  
  //Variables in the ODE solver
  double S0, I0, I_t, R0, S1, I1, R1, S0p, I0p, R0p, S1p, I1p, R1p, N, mu, beta, d;

  void operator() (  num_vec &x , num_vec &dxdt , const double  t )
  {
    
    int t_d = (int)t%365;
    
    beta = a + (b-a)*exp(-(t_d-phi)*(t_d-phi)/(2*psi*psi));
    
    for (int a = 0; a < A ; a++)
    {
      mu = a==0 ? mu_t : 0;
      
      I_t = 0.0;
      double I_t_temp = 0;
      for (int b = 0; b < A ; b++)
      {
        I_t_temp += (x[7*b+1]+x[7*b+4]*alpha)/popsize[b]*poly(b,a);
      }
      //I_t = min(I_t_temp, 1.0);
      I_t = I_t_temp;
      
      // SIR for age group and strain.
      S0 = x[7*a+0]; I0 = x[7*a+1]; R0 = x[7*a+2]; N = popsize[a];
      S1 = x[7*a+3]; I1 = x[7*a+4]; R1 = x[7*a+5];
      if (a == 0)
      {
        S0p = 0; I0p = 0; R0p = 0;
        S1p = 0; I1p = 0; R1p = 0;
      }
      else
      {
        S0p = x[7*(a-1)+0]; I0p = x[7*(a-1)+1]; R0p = x[7*(a-1)+2];
        S1p = x[7*(a-1)+3]; I1p = x[7*(a-1)+4]; R1p = x[7*(a-1)+5];
      }
      
      // ODEs transmission with aging
      dxdt[7*a+0] = mu - beta*I_t*S0        - S0*eta[a+1] + S0p*eta[a];
      dxdt[7*a+1] = beta*I_t*S0 - ga0*I0      - I0*eta[a+1] + I0p*eta[a];
      dxdt[7*a+2] = ga0*I0         - om*R0       - R0*eta[a+1] + R0p*eta[a];
      dxdt[7*a+3] = - d1*beta*I_t*S1  + om*R0  + om*R1     - S1*eta[a+1] + S1p*eta[a];
      dxdt[7*a+4] = d1*beta*I_t*S1 - ga0*ga1*I1      - I1*eta[a+1] + I1p*eta[a];
      dxdt[7*a+5] = ga0*ga1*I1        - om*R1        - R1*eta[a+1] + R1p*eta[a];
      
      // Number of new infections in age group a due to strain s
      dxdt[7*a+6] = beta*I_t*(S0 + d1*S1) + N*th*0.01;
    }
  }
};
};

using namespace model;


double llikelihood_func(VectorXd inc, NumericVector pars, int t_w, DataFrame y){
  
  int A = 2;
  VectorXd ep_t(A);
  ep_t << trans_par(pars["ep1"],0,0.01), trans_par(pars["ep2"],0,0.01);
  
  double ll = 0;
  for (int a = 0; a < A; a++)
  {
    NumericVector d = y[a];
    double N = inc(a);
    double x = d[t_w];
    if (x > N)
    {
      return log(0);
    }
    
    ll += stirl(N) - stirl(x) - stirl(N-x) + x*log(ep_t(a)) + (N-x)*log(1-ep_t(a));
  }
  return ll;
}

inline double check_stability(num_vec x, double t)
{
  int A = 2;
  long double X_w = 0;
  for (int j = 0; j < A*7; j++)
  {
    if (x[j] < 0)
    {
      return log(0);
    }
    else
      X_w += x[j];
  }
  
  if (std::isinf(X_w) || std::isnan(X_w))
    return log(0);
  else
    return 1.0;
}

// [[Rcpp::export]]
double eval_ll_cpp(DataFrame data, NumericMatrix poly, List dem, NumericVector age_lim, VectorXd pars_val)
{
  NumericVector pars = NumericVector::create(
    _["ga0"] = pars_val(0),
    _["om"] = pars_val(1),
    _["a"] = pars_val(2),
    _["b"] = pars_val(3),
    _["phi"] = pars_val(4),
    _["psi"] = pars_val(5),
    _["I1"] = pars_val(6),
    _["I2"] = pars_val(7),
    _["I3"] = pars_val(8),
    _["d1"] = pars_val(9),
    _["ep1"] = pars_val(10),
    _["ep2"] = pars_val(11),
    _["alpha"] = pars_val(12)
  );
  
  int A =  age_lim.size();
  num_vec popsize(A, 0);
  double mu = dem["mu"];
  double tot = dem["pop"];
  
  for (int i = 0; i < A-1; i++)
  {
    popsize[i] = (mu*365)*(age_lim[i+1]-age_lim[i]);
  }
  popsize[A-1] = tot - (mu*365)*age_lim[A-1];
  
  ODE_dynamics ode_par(0, 52*7, 52*7*8, 1.0);
  
  EulerT<state_t> integrator;
  SystemT<state_t, system_t> System;
  asc::Recorder recorder;
  
  VectorXd inc(A);
  VectorXd inc_tot(A);
  inc_tot = VectorXd::Zero(A);
  
  
        
  ODE_desc ODE_desc_inst(poly, dem, age_lim, popsize, pars);
  num_vec x0 = model::init_cond(dem, age_lim, popsize, pars);

  double t = ode_par.t_start;
  int t_d = 0;
  int t_w = 0;
  double ll = 0;
  double val;
  

  while (t < ode_par.t_end)
  {
    // Burn in section if required
    integrator(ODE_desc_inst, x0, t, 1.0);

    val = check_stability(x0, t);
    if (val < 0.0)
    {
      return log(0);
    }
    if (t > ode_par.t_burn)
    {
      if (t_d == 0)
      {
        for (int a = 0; a < A; a++)
          x0[6 + 7*a] = 0.0; //Incidence at t_d = 0;
      }
      if (t_d%7 == 0 && t_d > 0)
      {
        for (int a = 0; a < A; a++)
        {
          inc(a) = x0[6 + 7*a]; //Incidence at t_d = 7;
          x0[6 + 7*a] = 0.0;
        }
        ll += llikelihood_func(inc, pars, t_w, data);

        if (std::isinf(ll))
        {
          return log(0);
        }
        t_w ++;
      }
      t_d++;
    }
  }
  return ll; // RETURN THE LOG LIKELIHOOD
}
