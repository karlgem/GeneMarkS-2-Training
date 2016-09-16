//
//  ModuleUtilities.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/15/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "ModuleUtilities.hpp"

#include "LabelFile.hpp"
#include "SequenceFile.hpp"
#include "SequenceParser.hpp"

using namespace std;
using namespace gmsuite;


// constructor
ModuleUtilities::ModuleUtilities(const OptionsUtilities& opt) : options(opt) {
    
}

// Run module according to the provided options.
void ModuleUtilities::run() {
    
    if (options.utility == EXTRACT_UPSTR) {
        runExtractUpstream();
    }
    else            // unrecognized utility to run
        throw invalid_argument("Unknown utility function " + options.utility);
    
}



void ModuleUtilities::runExtractUpstream() {
    
    // read sequence file
    SequenceFile sequenceFile (options.extractUpstreamUtility.fn_sequence, SequenceFile::READ);
    Sequence strSequence = sequenceFile.read();
    
    // read label file
    LabelFile labelFile (options.extractUpstreamUtility.fn_label, LabelFile::READ);
    vector<Label*> labels;
    labelFile.read(labels);
    
    // create numeric sequence
    AlphabetDNA alph;
    CharNumConverter cnc (&alph);
    
    NumSequence numSequence (strSequence, cnc);
    
    // extract upstream regions from numeric sequence
    vector<NumSequence> upstreams;
    SequenceParser::extractUpstreamSequences(numSequence, labels, cnc, options.extractUpstreamUtility.length, upstreams, options.extractUpstreamUtility.allowOverlaps, options.extractUpstreamUtility.minimumGeneLength);
    
    // convert upstream numeric sequences back to string
    vector<Sequence> strUpstreams (upstreams.size());
    
    for (size_t n = 0; n < upstreams.size(); n++) {
        strUpstreams[n] = Sequence(cnc.convert(upstreams[n].begin(), upstreams[n].end()));
    }
    
    
    // write sequences to file
    SequenceFile outputFile(options.extractUpstreamUtility.fn_output, SequenceFile::WRITE);
    outputFile.write(strUpstreams);
    
    
}
