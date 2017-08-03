//
//  UnivariatePDF.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/5/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef UnivariatePDF_hpp
#define UnivariatePDF_hpp

#include <stdio.h>
#include <vector>
#include <string>

using std::string;
using std::vector;

namespace gmsuite {
    
    /**
     * @class UnivariatePDF
     * @brief This class implements a linear distribution
     *
     * Operations such as sample and sample-in-log-space allow simple usage
     */
    class UnivariatePDF {
        
    public:
        
        /**
         * Default constructor
         */
        UnivariatePDF();
        
        
        /**
         * Constructor: create a distribution from set of weights. These can be simple counts
         * to actual weights, from which a normalized probably distribution will
         * be created.
         *
         * @param weights the set of weights
         * @param lnspace if set, then the weights are converted from natural logspace back
         * into linear space, before the normalization is created.
         * @param normalize whether the weights should be normalized to sum to one
         */
        UnivariatePDF(const vector<double> &weights, bool lnspace=false, double pcounts=0, bool normalize=true);
        
        /**
         * Construct distribution from set of weights. These can be simple counts
         * to actual weights, from which a normalized probably distribution will
         * be created.
         *
         * @param weights the set of weights
         * @param lnspace if set, then the weights are converted from natural logspace back
         * into linear space, before the normalization is created.
         * @param normalize whether the weights should be normalized to sum to one
         */
        void construct(const vector<double> &weights, bool lnspace=false, double pcounts=0, bool normalize=true);
        
        
        /**
         * Sample a value from this distribution
         */
        size_t sample() const;
        
        
        const double& operator[] (size_t pos) const;
        
        size_t size() const;
        
        string toString() const;
        
        
        /**
         * Compute a localization metric of the distribution.
         */
        typedef struct {
            size_t windowBegin;         // the starting position of the window
            size_t windowLength;        // the length of the window (must be > 0)
            double windowTotal;         // the total amount of probability in the window
        } localization_metric_t;
        
        localization_metric_t localization(size_t windowLength = 0) const;
        
        
    private:
        
        vector<double> probabilities;           /**< the distribution's probabilities */
        vector<double> cumulative;              /**< cumulative distribution (for sampling) */
        
        /**
         * Convert values from log-space to linear space
         */
        void convertFromLinearToLogspace();
        
        /**
         * Compute the cumulative distribution for easy sampling
         */
        void computeCumulativeDistribution();
        
        /**
         * Normalize probability distribution
         */
        void normalize();
        
    };
}


#endif /* UnivariatePDF_hpp */
