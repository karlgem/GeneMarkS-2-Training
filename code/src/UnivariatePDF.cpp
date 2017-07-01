//
//  UnivariatePDF.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/5/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "UnivariatePDF.hpp"

#include <algorithm>
#include <stdexcept>
#include <stdlib.h>
#include <sstream>
#include <math.h>

using namespace gmsuite;
using namespace std;

// default constructor
UnivariatePDF::UnivariatePDF() {
    
}


UnivariatePDF::UnivariatePDF(const vector<double> &weights, bool lnspace, double pcounts, bool norm) {
    
    construct(weights, lnspace, pcounts, norm);
}

// construct distribution from weights
void UnivariatePDF::construct(const vector<double> &weights, bool lnspace, double pcounts, bool norm) {
    this->probabilities = weights;
    
    // check if conversion necessary
    if (lnspace)
        convertFromLinearToLogspace();
    
    // add pseudocounts
    if (pcounts > 0) {
        for (size_t n = 0; n < probabilities.size(); n++)
            probabilities[n] += pcounts;
    }
    
    if (norm)
        normalize();
    
    // compute cumulative distribution (for easy sampling)
    computeCumulativeDistribution();
}

/**
 * Sample a value from this distribution
 */
size_t UnivariatePDF::sample() const {
    
    if (cumulative.size() == 0)
        throw std::domain_error("Cannot sample from an empty distribution");
    
    // generate a uniform random number between 0 and max(cumulative distribution), and select corresponding sample
    double u = rand() / (double) RAND_MAX;
    double b = cumulative[probabilities.size()-1];
    
    u = u * b;
    
    size_t s = cumulative.size();
    
    // find location i where u < cumulative[i]
    for (size_t i = 0; i < cumulative.size(); i++) {
        if (u <= cumulative[i]) {
            s = i;
            break;
        }
    }
    
    return s;
}



// convert values in <probabilities> from log-space to linear space
void UnivariatePDF::convertFromLinearToLogspace() {
    
    if (probabilities.size() == 0)
        return;
    
    // apply exponential
    for (size_t i = 0; i < probabilities.size(); i++)
        probabilities[i] = exp(probabilities[i]);
}

void UnivariatePDF::computeCumulativeDistribution() {
    
    if (probabilities.size() == 0)
        return;
    
    // get minimum value
    double minVal = *std::min_element(probabilities.begin(), probabilities.end());
    
    // make sure minimum value is not negative
    if (minVal < 0)
        throw std::domain_error("Cannot normalize weights with negative values.");
    
    cumulative.resize(probabilities.size(), 0);
    cumulative[0] = probabilities[0];
    
    for (size_t i = 1; i < probabilities.size(); i++) {
        cumulative[i] = probabilities[i] + cumulative[i-1];
    }
}






void UnivariatePDF::normalize() {
    
    double total = 0;
    for (size_t i = 0; i < probabilities.size(); i++) {
        total += probabilities[i];
    }
    
    if (total > 0) {
        for (size_t i = 0; i < probabilities.size(); i++)
            probabilities[i] /= total;
    }
    
}



const double& UnivariatePDF::operator[] (size_t pos) const {
    return probabilities.at(pos);
}




size_t UnivariatePDF::size() const {
    return probabilities.size();
}


string UnivariatePDF::toString() const {
    stringstream ssm; ssm << fixed;
    
    for (size_t n = 0; n < probabilities.size(); n++) {
        ssm << n << "\t" << probabilities[n] << endl;
    }
    
    return ssm.str();
}





UnivariatePDF::localization_metric_t UnivariatePDF::localization(size_t windowLength) const {
    
    if (windowLength > this->probabilities.size())
        throw invalid_argument("Window Length cannot be larger than probability vector size.");
    
    localization_metric_t result;
    result.windowBegin = 0;
    result.windowLength = windowLength;
    result.windowTotal = 0;
    
    // make sure window length is less than the probability vector size
    
    // loop over distribution, accumulating results and computing the max.
    
    // for every location
    for (size_t n = 0; n < this->probabilities.size()-windowLength+1; n++) {
        
        // sum over all elements in window
        double accum = 0;
        for (size_t i = 0; i < windowLength; i++) {
            accum += this->probabilities[i+n];
        }
        
        // check if this location has max
        if (accum >= result.windowTotal) {
            result.windowTotal = accum;
            result.windowBegin = n;
        }
    }
    
    return result;
}
















