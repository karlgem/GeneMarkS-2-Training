//
//  OldGMS2ModelFile.hpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 11/7/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#ifndef OldGMS2ModelFile_hpp
#define OldGMS2ModelFile_hpp

#include <stdio.h>
#include <string>
#include <vector>

using std::string;
using std::vector;
using std::pair;

namespace gmsuite {
    
    class OldGMS2ModelFile {
        
    public:
        
        /**
         * Create an model file with a filename
         */
        OldGMS2ModelFile(std::string filename);
        
        vector<pair<string, double> > getNoncoding() const;
        
        size_t getNDEC() const;
        
        
    private:
        
        std::string filename;
        
    };
    
}

#endif /* OldGMS2ModelFile_hpp */
