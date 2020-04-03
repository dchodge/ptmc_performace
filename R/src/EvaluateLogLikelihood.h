//
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



class EvaluateLogLikelihood
{
public:
    // Need to be defined in constructor
    int A;
    num_vec populationPerAgeGroup, eta, modelIncidencePerTime;
    double dailyBirthRate, totPopulation;
    double t_start, t_burn, t_end, dt;

    // Defined later but of size A
    NumericVector ageStratification, ep_t;
    double currentODETime;
    bool loglikelihoodError;
    
    EvaluateLogLikelihood(double dailyBirthRate_t, double totPopulation_t, NumericVector ageStratification_t): dailyBirthRate(dailyBirthRate_t), totPopulation(totPopulation_t), ageStratification(ageStratification_t){
       
        A = ageStratification.size();
        eta.push_back(0);
        for (int i = 0; i < A-1; i++){
            populationPerAgeGroup.push_back(dailyBirthRate*365*(ageStratification[i+1]-ageStratification[i]));
            eta.push_back( 1.0/(365.0*(ageStratification[i+1] - ageStratification[i])) );
            modelIncidencePerTime.push_back(0);
        }            
        modelIncidencePerTime.push_back(0);
        populationPerAgeGroup.push_back(totPopulation - (dailyBirthRate*365)*ageStratification[A-1]);
        eta.push_back(dailyBirthRate/(totPopulation - (dailyBirthRate*365)*ageStratification[A-1]) );

        t_start = 0;
        t_burn = 52*7;
        t_end = 52*7*8;
        dt = 1;
        currentODETime = 0;
        loglikelihoodError = false;
    }
    
    double evaluateLogLikelihoodCpp(VectorXd currentParamValues);
    
    // define in R after constructure
    NumericMatrix contactMatrix;
    NumericMatrix observedData;
    
    // Related to parameter values
    NumericVector parameterValuesTransformed;
    
    void transformParameterValuesforODE(VectorXd currentParamValues){

        this->parameterValuesTransformed = NumericVector::create(
                                                                 _["ga0"] = 1.0/logisticTransform(currentParamValues(0), 2, 20),
                                                                 _["om"] = 1.0/logisticTransform(currentParamValues(1), 60, 365),
                                                                 _["a"] = logisticTransform(currentParamValues(2), 0, 1.0),
                                                                 _["b"] = logisticTransform(currentParamValues(3), 0, 1),
                                                                 _["phi"] = logisticTransform(currentParamValues(4), 100, 250),
                                                                 _["psi"] = logisticTransform(currentParamValues(5), 0, 0.1*365),
                                                                 _["I1"] = logisticTransform(currentParamValues(6), 0, 1),
                                                                 _["I2"] = logisticTransform(currentParamValues(7), 0, 1),
                                                                 _["I3"] = logisticTransform(currentParamValues(8), 0.25, 1),
                                                                 _["d1"] = logisticTransform(currentParamValues(9), 0, 1),
                                                                 _["ep1"] = logisticTransform(currentParamValues(10), 0,0.01),
                                                                 _["ep2"] = logisticTransform(currentParamValues(11), 0,0.01),
                                                                 _["alpha"] = logisticTransform(currentParamValues(12), 0, 1)
                                                                 );

      this->ep_t = NumericVector::create(this->parameterValuesTransformed["ep1"], this->parameterValuesTransformed["ep2"]);
    }
    
    inline double logisticTransform(double x, double a, double b){
        if (x<-100){return a;}
        else if (x>100){return b;}
        else {return a + (b-a)*exp(x)/(exp(x) + 1);}
    }
    
    num_vec generateInitialStates(){
        num_vec initialStates;
        double a1, a2;
        double I1 = this->parameterValuesTransformed["I1"];
        double I2 = this->parameterValuesTransformed["I2"];
        double I3 = this->parameterValuesTransformed["I3"];
        
        for (int a = 0; a < this->A ; a++){
            if (a < A-1){
                a1 = this->ageStratification(a); a2 =  this->ageStratification(a+1);
            }
            else{
                a1 =  this->ageStratification(a); a2 = 90;
            }
            
            num_vec propEachExposureGroup = initialProportionExposure(I3, a1, a2);
            
            initialStates.push_back(this->populationPerAgeGroup[a]*propEachExposureGroup[0]*(1-I1*1.0)*(1-I2));
            initialStates.push_back(this->populationPerAgeGroup[a]*propEachExposureGroup[0]*I1*1.0);
            initialStates.push_back(this->populationPerAgeGroup[a]*propEachExposureGroup[0]*(1-I1*1.0)*I2);
            
            initialStates.push_back(this->populationPerAgeGroup[a]*propEachExposureGroup[1]*(1-I1*1.0)*(1-I2));
            initialStates.push_back(this->populationPerAgeGroup[a]*propEachExposureGroup[1]*I1*1.0);
            initialStates.push_back(this->populationPerAgeGroup[a]*propEachExposureGroup[1]*(1-I1*1.0)*I2);
            initialStates.push_back(0);
        }
        return initialStates;
    }
    
    num_vec initialProportionExposure(double l, double a1, double a2){
        num_vec prop(this->A);
        prop[0] = abs(poisson_cdf(l,a2,0)-poisson_cdf(l,a1,0))/((a2-a1)*l);
        prop[1] = 1 - prop[0];
        return prop;
    }
    
    double poisson_cdf(double l, double a, double x){
      if( l == 0.0 || a == 0.0){
        boost::math::poisson_distribution<> p(0.000001); return cdf(p,x);
      }
      else{
        boost::math::poisson_distribution<> p(l*a); return cdf(p,x);
      }
    }
    
    double evaluateTimeStepLogLikelihood(int weekNo){
        double ll = 0;
        double estimatedLogBinomialCoeff;
        for (int a = 0; a < this->A; a++){
           // NumericVector d = this->observedData[weekNo, a];
            double dataNewInfections = this->observedData(weekNo, a);
            if (dataNewInfections > this->modelIncidencePerTime[a]){
                return log(0);
            }
            
            estimatedLogBinomialCoeff = stirlingApproximation(this->modelIncidencePerTime[a]) - stirlingApproximation(dataNewInfections) - stirlingApproximation(this->modelIncidencePerTime[a]-dataNewInfections);
            ll += estimatedLogBinomialCoeff + dataNewInfections*log(this->ep_t(a)) + (this->modelIncidencePerTime[a]-dataNewInfections)*log(1-this->ep_t(a));
        }
        if (std::isinf(ll)){
          this->loglikelihoodError=true;
        }
        return ll;
    }
    
    long double stirlingApproximation(double n){
        if (n == 0)
            return 0;
        else {
            double x = n + 1;
            return (x - 0.5)*log(x) - x + 0.5*log(2*PI) + 1.0/(12*x) - 1.0/(360.0*pow(x,3)); // https://en.wikipedia.org/wiki/Stirling%27s_approximation#Speed_of_convergence_and_error_estimates
        }
    }
    
    inline void checkStability(num_vec x)
    {
        long double X_w = 0;
        // check state values are never negative
        for (int j = 0; j < this->A*7; j++){
            if (x[j] < 0) {
//Rcout << "Negative number at position: " << j << ". Value: " << x[j] << ". At time: " << this->currentODETime << endl;
              this->loglikelihoodError=true; 
              return; }
            else {X_w += x[j]; }
        }
        // Check state values are non-infinte numeric values
        if (std::isinf(X_w) || std::isnan(X_w)) {
            this->loglikelihoodError=true;
            return;
        }
    }
};

class ODE_desc
{
  EvaluateLogLikelihood* finELL;
  
public:
  double S0, I0, I_t, R0, S1, I1, R1, S0p, I0p, R0p, S1p, I1p, R1p, N, mu, beta, d;
  double d1, ga0, ga1, alpha, om, a, b, psi, phi;
  double dailyBirthRate;
  NumericMatrix contactMatrix;
  int A;
  num_vec populationPerAgeGroup, eta;
  
  ODE_desc(EvaluateLogLikelihood* finELL_t): finELL(finELL_t){
    
    ga0 = finELL->parameterValuesTransformed["ga0"];
    ga1 = 1;
    om = finELL->parameterValuesTransformed["om"];
    a = finELL->parameterValuesTransformed["a"];
    b = finELL->parameterValuesTransformed["b"];
    psi = finELL->parameterValuesTransformed["psi"];
    phi = finELL->parameterValuesTransformed["phi"];
    d1 = finELL->parameterValuesTransformed["d1"];
    alpha = finELL->parameterValuesTransformed["alpha"];
    A = finELL->A;
    eta = finELL->eta;
    populationPerAgeGroup = finELL->populationPerAgeGroup;
    contactMatrix = finELL->contactMatrix;
    dailyBirthRate = finELL->dailyBirthRate;
   };
  
  
  void operator() (  num_vec &x , num_vec &dxdt , const double  t )
  {
    int t_d = (int)t%365;
    beta = a + (b-a)*exp(-(t_d-phi)*(t_d-phi)/(2*psi*psi));
    
    for (int a = 0; a < A ; a++)
    {
      mu = a==0 ? dailyBirthRate : 0;
      
      I_t = 0.0;
      for (int b = 0; b < A ; b++){
        I_t += (x[7*b+1]+x[7*b+4]*alpha)/populationPerAgeGroup[b]*contactMatrix(a,b);
      }
      //I_t = min(I_t, 1.0);
      
      S0 = x[7*a+0]; I0 = x[7*a+1]; R0 = x[7*a+2]; N = populationPerAgeGroup[a];
      S1 = x[7*a+3]; I1 = x[7*a+4]; R1 = x[7*a+5];
      
      if (a == 0){
        S0p = 0; I0p = 0; R0p = 0; S1p = 0; I1p = 0; R1p = 0;
      }
      else{
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
      dxdt[7*a+6] = beta*I_t*(S0 + d1*S1);
    }
  }
};

double EvaluateLogLikelihood::evaluateLogLikelihoodCpp(VectorXd currentParamValues)
{
  transformParameterValuesforODE(currentParamValues);
  this->loglikelihoodError = false;
  
  EulerT<state_t> integrator;
  SystemT<state_t, system_t> System;
  asc::Recorder recorder;

  num_vec x0 = generateInitialStates();
  ODE_desc ODE_desc_inst(this);
  this->currentODETime = this->t_start;
  int dayNoAfterBurn = 0;
  int weekNo = 0;
  double valueLogLikelihood = 0;

  while (this->currentODETime < this->t_end){
    integrator(ODE_desc_inst, x0, this->currentODETime, this->dt);
    checkStability(x0);
    if (this->loglikelihoodError){
      return log(0);
    }
    
    if (this->currentODETime > this->t_burn)
    {
      if (dayNoAfterBurn == 0){
        for (int a = 0; a < this->A; a++)
          x0[6 + 7*a] = 0.0; //Incidence at t_d = 0;
      }
      if (dayNoAfterBurn%7 == 0 && dayNoAfterBurn > 0)
      {
        for (int a = 0; a < this->A; a++){
          this->modelIncidencePerTime[a] = x0[6 + 7*a]; //Incidence at t_d = 7;
          x0[6 + 7*a] = 0.0;
        }
        valueLogLikelihood += evaluateTimeStepLogLikelihood(weekNo);
        if (this->loglikelihoodError){
          return log(0);
        }
        weekNo ++;
      }
      dayNoAfterBurn++;
    }
  }
  return valueLogLikelihood; // RETURN THE LOG LIKELIHOOD
}
