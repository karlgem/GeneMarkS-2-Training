//
//  LabelsParser.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 10/19/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef LabelsParser_hpp
#define LabelsParser_hpp

#include <stdio.h>
#include <vector>

#include "Label.hpp"

using std::vector;

namespace gmsuite {
    
    /**
     * @class LabelsParser
     * @brief parse and filter labels
     */
    class LabelsParser {
        
    public:
        
        typedef enum {FGIO, NFGIO, AMBIG} operon_status_t;
        
        static void partitionBasedOnOperonStatus(const vector<Label*> &labels, size_t fgioThresh, size_t nfgioThresh,
                                                 vector<operon_status_t> &status);
        
    };
}

#endif /* LabelsParser_hpp */
