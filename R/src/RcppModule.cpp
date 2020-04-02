//
//
//  Created by David Hodgson on 03/02/2020.
//  Copyright © 2020 David Hodgson. All rights reserved.
//

#include <Rcpp.h>
#include "EvaluateLogLikelihood.h"

RCPP_MODULE(EvaluateLogLikelihoodModule) {
    using namespace Rcpp;
    
    class_<EvaluateLogLikelihood>( "EvaluateLogLikelihood" )
    .constructor<double, double, NumericVector>()
    .field( "contactMatrix", &EvaluateLogLikelihood::contactMatrix )
    .field( "observedData", &EvaluateLogLikelihood::observedData )
    //.method( "evaluateLogLikelihoodCpp", &EvaluateLogLikelihood::evaluateLogLikelihoodCpp )
    ;
    /*.constructor<double,double, NumericVector>()
    .field( "contactMatrix", &EvaluateLogLikelihood::contactMatrix )
    .field( "observedData", &EvaluateLogLikelihood::observedData )
    .method( "evaluateLogLikelihood", &EvaluateLogLikelihood::evaluateLogLikelihoodCpp )
    ;*/
}
