//
//  OptionsUtilities.cpp
//  GeneMark Suite
//
//  Created by Karl Gemayel on 9/15/16.
//  Copyright © 2016 Karl Gemayel. All rights reserved.
//

#include "OptionsUtilities.hpp"

#include <vector>
#include <fstream>
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/any.hpp>

using namespace std;
using namespace gmsuite;
using namespace boost::xpressive;
namespace po = boost::program_options;

OptionsUtilities::OptionsUtilities(string mode) : Options(mode) {
    
}

#define STR_EXTRACT_UPSTR "extract_upstream"
#define STR_LABELS_SIMILARITY_CHECK "labels-sim"
#define STR_START_MODEL_INFO "start-model-info"
#define STR_MATCH_SEQ_TO_UPSTREAM "match-seq-to-upstream"
#define STR_MATCH_SEQ_TO_NONCODING "match-seq-to-noncoding"
#define STR_EMIT_NONCODING "emit-noncoding"
#define STR_COUNT_NUM_ORF "count-num-orfs"
#define STR_EXTRACT_SC_PER_OPERON_STATUS "extract-sc-per-operon-status"
#define STR_EXTRACT_SC_PER_MOTIF_STATUS "extract-sc-per-motif-status"
#define STR_COMPUTE_GC "compute-gc"
#define STR_SEPARATE_FGIO_AND_IG "separate-fgio-and-ig"
#define STR_EXTRACT_START_CONTEXT "extract-start-context"
#define STR_DNA_TO_AA "dna-to-aa"

namespace gmsuite {
    // convert string to utility_t
    std::istream& operator>>(std::istream& in, OptionsUtilities::utility_t& unit) {
        std::string token;
        in >> token;
        
        if (token == STR_EXTRACT_UPSTR)                 unit = OptionsUtilities::EXTRACT_UPSTR;
        else if (token == STR_START_MODEL_INFO)         unit = OptionsUtilities::START_MODEL_INFO;
        else if (token == STR_MATCH_SEQ_TO_UPSTREAM)    unit = OptionsUtilities::MATCH_SEQ_TO_UPSTREAM;
        else if (token == STR_MATCH_SEQ_TO_NONCODING)   unit = OptionsUtilities::MATCH_SEQ_TO_NONCODING;
        else if (token == STR_LABELS_SIMILARITY_CHECK)  unit = OptionsUtilities::LABELS_SIMILARITY_CHECK;
        else if (token == STR_EMIT_NONCODING)           unit = OptionsUtilities::EMIT_NON_CODING;
        else if (token == STR_COUNT_NUM_ORF)           unit = OptionsUtilities::COUNT_NUM_ORF;
        else if (token == STR_EXTRACT_SC_PER_OPERON_STATUS) unit = OptionsUtilities::EXTRACT_SC_PER_OPERON_STATUS;
        else if (token == STR_EXTRACT_SC_PER_MOTIF_STATUS) unit = OptionsUtilities::EXTRACT_SC_PER_MOTIF_STATUS;
        else if (token == STR_COMPUTE_GC)               unit = OptionsUtilities::COMPUTE_GC;
        else if (token == STR_SEPARATE_FGIO_AND_IG)     unit = OptionsUtilities::SEPARATE_FGIO_AND_IG;
        else if (token == STR_EXTRACT_START_CONTEXT)    unit = OptionsUtilities::EXTRACT_START_CONTEXT;
        else if (token == STR_DNA_TO_AA)                unit = OptionsUtilities::DNA_TO_AA;
        else
            throw boost::program_options::invalid_option_value(token);
        
        return in;
    }
    
//    // validate utility option
//    void validate(boost::any& v,
//                  const vector<string>& values,
//                  OptionsUtilities::utility_t * target_type, int) {
//        
//        static sregex expr = sregex::compile("^\\s*(\\S+)\\s*");
//        
//        po::validators::check_first_occurrence(v);  // check no previous assignment to v
//        const string& s = po::validators::get_single_string(values);  // get value
//        
//        smatch match;
//        if (regex_search(s,match,expr))
//    }
    
}

// parse CMD options
bool OptionsUtilities::parse(int argc, const char *argv[]) {
    
    
    try {
        vector<string> config_fnames;                           // holds names of all configuration files (specified by user)
        
        // Declare a group of options that will be allowed
        // only on the command line (CML)
        po::options_description generic("General Options");
        generic.add_options()
        ("version", "Print version string")
        ("help,h", "Display help message")
        ("config", po::value(&config_fnames), "Config file where options may be specified (can be specified more than once)")
        ;
        
        // Declare a group of options that will be allowed
        // on both CML and in the config files
        po::options_description config("Configuration");
        config.add_options()
        ("verbose,v", po::value<int>(&verbose)->default_value(0), "Verbose level")
        ;
        
        // Create set of hidden arguments (which can correspond to positional arguments). This is used
        // to add positional arguments, while not putting their description in the "options" section.
        // Hidden options are allowed in both CML and config files
        po::options_description hidden;
        hidden.add_options()
        ("mode", po::value<string>(&mode)->required(), "Program Mode")
        ("utility", po::value<utility_t>(&utility)->required(), "Utility function")
        ("subargs", po::value<std::vector<std::string> >(), "Arguments for command")
        ;
        
        // Congregate options into further groups
        po::options_description cmdline_options;
        cmdline_options.add(generic).add(config).add(hidden);
        
        po::options_description config_file_options;
        config_file_options.add(config).add(hidden);
        
        po::options_description visible("Allowed options");
        visible.add(generic).add(config);
        
        // Specify positional arguments
        po::positional_options_description pos;
        pos.add("mode",1);
        pos.add("utility",1);
        pos.add("subargs",-1);
        
        
        
        // create storage component for storing names and values of arguments
        po::variables_map vm;
        
        // store command-line options
        po::parsed_options parsed = po::command_line_parser(argc, argv).          // pass in input
            options(cmdline_options).                     // specify options list
            positional(pos).                              // specify which are positional
            allow_unregistered().                         // allow unregistered args since we haven't yet set-up the utility-based ones
            run();
        
        po::store(parsed,vm);
        
        if (!vm.count("utility")){
                cout << make_usage_string(basename(argv[0]), cmdline_options, pos) << endl;
                return false;
        }
        
        // try parsing arguments.
        po::notify(vm);
        
        // Utility: Extract upstream
        if (utility == EXTRACT_UPSTR) {
            
            po::options_description utilDesc(string(STR_EXTRACT_UPSTR) + " options");
            utilDesc.add_options()
                ("sequence,s", po::value<string>(&extractUpstreamUtility.fn_sequence)->required(), "Sequence filename")
                ("label,l", po::value<string>(&extractUpstreamUtility.fn_label)->required(), "Label filename")
                ("output,o", po::value<string>(&extractUpstreamUtility.fn_output)->required(), "Output filename")
                ("length", po::value<size_t>(&extractUpstreamUtility.length)->required(), "Upstream length")
                ("allow-overlap-with-cds", "If set, then upstream (non-coding) regions are allowed to overlap with coding regions. If not set, these sequences are ignored.")
                ("min-gene-length", po::value<size_t>(&extractUpstreamUtility.minimumGeneLength)->default_value(0), "Minimum gene length")
            ;
            
            cmdline_options.add(utilDesc);
            
            // Collect all the unrecognized options from the first pass. This will include the
            // (positional) mode and command name, so we need to erase them
            vector<string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
            opts.erase(opts.begin());       // erase mode
            opts.erase(opts.begin());       // erase command name
            
            // Parse again...
            po::store(po::command_line_parser(opts).options(utilDesc).run(), vm);
            
            // get remaining parameters whose values were not assigned in add_options() above
            extractUpstreamUtility.allowOverlaps = vm.count("allow-overlap-with-cds") > 0;
        }
        // Similarity Check
        else if (utility == LABELS_SIMILARITY_CHECK) {
            po::options_description utilDesc(string(STR_LABELS_SIMILARITY_CHECK) + " options");
            addProcessOptions_LabelsSimilarityCheck(labelsSimilarityCheck, utilDesc);
            
            cmdline_options.add(utilDesc);
            
            // Collect all the unrecognized options from the first pass. This will include the
            // (positional) mode and command name, so we need to erase them
            vector<string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
            opts.erase(opts.begin());       // erase mode
            opts.erase(opts.begin());       // erase command name
            
            // Parse again...
            po::store(po::command_line_parser(opts).options(utilDesc).run(), vm);
            
        }
        // Emit Noncoding
        else if (utility == EMIT_NON_CODING) {
            po::options_description utilDesc(string(STR_EMIT_NONCODING) + " options");
            addProcessOptions_EmitNonCoding(emitNonCoding, utilDesc);
            
            cmdline_options.add(utilDesc);
            
            // Collect all the unrecognized options from the first pass. This will include the
            // (positional) mode and command name, so we need to erase them
            vector<string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
            opts.erase(opts.begin());       // erase mode
            opts.erase(opts.begin());       // erase command name
            
            // Parse again...
            po::store(po::command_line_parser(opts).options(utilDesc).run(), vm);
            
        }
        // Emit Noncoding
        else if (utility == COUNT_NUM_ORF) {
            po::options_description utilDesc(string(STR_COUNT_NUM_ORF) + " options");
            addProcessOptions_CountNumORF(countNumORF, utilDesc);
            
            cmdline_options.add(utilDesc);
            
            // Collect all the unrecognized options from the first pass. This will include the
            // (positional) mode and command name, so we need to erase them
            vector<string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
            opts.erase(opts.begin());       // erase mode
            opts.erase(opts.begin());       // erase command name
            
            // Parse again...
            po::store(po::command_line_parser(opts).options(utilDesc).run(), vm);
            
        }
        // Utility: Start-Model Info
        else if (utility == START_MODEL_INFO) {
            
            po::options_description utilDesc(string(STR_START_MODEL_INFO) + " options");
            utilDesc.add_options()
                ("sequence,s", po::value<string>(&startModelInfoUtility.fn_sequence)->required(), "Sequence filename")
                ("label,l", po::value<string>(&startModelInfoUtility.fn_label)->required(), "Label filename")
                ("num-sim-noncoding", po::value<size_t>(&startModelInfoUtility.numOfSimNonCoding)->default_value(1000), "Number of simulated non-coding sequences.")
            ;
            
            // gms2 training options
            po::options_description gms2training("GMS2 Training");
            OptionsGMS2Training::addProcessOptions(startModelInfoUtility.optionsGMS2Training, gms2training);
            
            utilDesc.add(gms2training);
            
            cmdline_options.add(utilDesc);
            
            // Collect all the unrecognized options from the first pass. This will include the
            // (positional) mode and command name, so we need to erase them
            vector<string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
            opts.erase(opts.begin());       // erase mode
            opts.erase(opts.begin());       // erase command name
            
            // Parse again...
            po::store(po::command_line_parser(opts).options(utilDesc).run(), vm);
            
            // get remaining parameters whose values were not assigned in add_options() above
            startModelInfoUtility.allowOverlaps = vm.count("allow-overlap-with-cds") > 0;
            
            
            
        }
        // Utility: Match sequence to upstream
        else if (utility == MATCH_SEQ_TO_UPSTREAM) {
            po::options_description utilDesc(string(STR_MATCH_SEQ_TO_UPSTREAM) + " options");
            utilDesc.add_options()
                ("match-to", po::value<string>(&matchSeqWithUpstream.matchTo)->required(), "Sequence to be matched to")
            ;
            
            // add extract upstream process
            po::options_description extractUpstreamOptions(string(STR_EXTRACT_UPSTR) + " options");
            addProcessOptions_ExtractUpstream(matchSeqWithUpstream, extractUpstreamOptions);
            
            utilDesc.add(extractUpstreamOptions);
            
            cmdline_options.add(utilDesc);
            
            // Collect all the unrecognized options from the first pass. This will include the
            // (positional) mode and command name, so we need to erase them
            vector<string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
            opts.erase(opts.begin());       // erase mode
            opts.erase(opts.begin());       // erase command name
            
            // Parse again...
            po::store(po::command_line_parser(opts).options(utilDesc).run(), vm);
            
            // get remaining parameters whose values were not assigned in add_options() above
            extractUpstreamUtility.allowOverlaps = vm.count("allow-overlap-with-cds") > 0;
            
        }
        // Utility: Match sequence to noncoding
        else if (utility == MATCH_SEQ_TO_NONCODING) {
            po::options_description utilDesc(string(STR_MATCH_SEQ_TO_NONCODING) + " options");
            utilDesc.add_options()
                ("match-to", po::value<string>(&matchSeqWithNoncoding.matchTo)->required(), "Sequence to be matched to")
            ;
            
            // add start-model-info
            po::options_description startModelInfoOptions(string(STR_START_MODEL_INFO) + " options");
            addProcessOptions_StartModelInfo(matchSeqWithNoncoding, startModelInfoOptions);
            
            utilDesc.add(startModelInfoOptions);
            
            cmdline_options.add(utilDesc);
            
            // Collect all the unrecognized options from the first pass. This will include the
            // (positional) mode and command name, so we need to erase them
            vector<string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
            opts.erase(opts.begin());       // erase mode
            opts.erase(opts.begin());       // erase command name
            
            // Parse again...
            po::store(po::command_line_parser(opts).options(utilDesc).run(), vm);
            
            // get remaining parameters whose values were not assigned in add_options() above
            matchSeqWithNoncoding.allowOverlaps = vm.count("allow-overlap-with-cds") > 0;
        }
        else if (utility == EXTRACT_SC_PER_OPERON_STATUS) {
            po::options_description utilDesc (string(STR_EXTRACT_SC_PER_OPERON_STATUS) + " options");
            addProcessOptions_ExtractStartContextPerOperonStatus(extractStartContextPerOperonStatus, utilDesc);
            
            cmdline_options.add(utilDesc);
            
            // Collect all the unrecognized options from the first pass. This will include the
            // (positional) mode and command name, so we need to erase them
            vector<string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
            opts.erase(opts.begin());       // erase mode
            opts.erase(opts.begin());       // erase command name
            
            // Parse again...
            po::store(po::command_line_parser(opts).options(utilDesc).run(), vm);
            
            // get remaining parameters whose values were not assigned in add_options() above
            extractStartContextPerOperonStatus.allowOverlaps = vm.count("allow-overlap-with-cds") > 0;
            
            
            
        }
        else if (utility == EXTRACT_SC_PER_MOTIF_STATUS) {
            po::options_description utilDesc (string(STR_EXTRACT_SC_PER_MOTIF_STATUS) + " options");
            addProcessOptions_ExtractStartContextPerMotifStatus(extractStartContextPerMotifStatus, utilDesc);
            
            cmdline_options.add(utilDesc);
            
            // Collect all the unrecognized options from the first pass. This will include the
            // (positional) mode and command name, so we need to erase them
            vector<string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
            opts.erase(opts.begin());       // erase mode
            opts.erase(opts.begin());       // erase command name
            
            // Parse again...
            po::store(po::command_line_parser(opts).options(utilDesc).run(), vm);
            
            // get remaining parameters whose values were not assigned in add_options() above
            extractStartContextPerMotifStatus.allowOverlaps = vm.count("allow-overlap-with-cds") > 0;
            
            
            
        }
        else if (utility == COMPUTE_GC) {
            po::options_description utilDesc (string(STR_COMPUTE_GC) + " options");
            addProcessOptions_ComputeGC(computeGC, utilDesc);
            
            cmdline_options.add(utilDesc);
            
            // Collect all the unrecognized options from the first pass. This will include the
            // (positional) mode and command name, so we need to erase them
            vector<string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
            opts.erase(opts.begin());       // erase mode
            opts.erase(opts.begin());       // erase command name
            
            // Parse again...
            po::store(po::command_line_parser(opts).options(utilDesc).run(), vm);
        }
        else if (utility == SEPARATE_FGIO_AND_IG) {
            po::options_description utilDesc (string(STR_SEPARATE_FGIO_AND_IG) + " options");
            addProcessOptions_SeparateFGIOAndIG(separateFGIOAndIG, utilDesc);
            
            cmdline_options.add(utilDesc);
            
            // Collect all the unrecognized options from the first pass. This will include the
            // (positional) mode and command name, so we need to erase them
            vector<string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
            opts.erase(opts.begin());       // erase mode
            opts.erase(opts.begin());       // erase command name
            
            // Parse again...
            po::store(po::command_line_parser(opts).options(utilDesc).run(), vm);
        }
        // Extract start-contxt
        else if (utility == EXTRACT_START_CONTEXT) {
            po::options_description utilDesc (string(STR_EXTRACT_START_CONTEXT) + " options");
            addProcessOptions_ExtractStartContext(extractStartContext, utilDesc);
            
            cmdline_options.add(utilDesc);
            
            // Collect all the unrecognized options from the first pass. This will include the
            // (positional) mode and command name, so we need to erase them
            vector<string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
            opts.erase(opts.begin());       // erase mode
            opts.erase(opts.begin());       // erase command name
            
            // Parse again...
            po::store(po::command_line_parser(opts).options(utilDesc).run(), vm);
        }
        // Convert DNA to AA
        else if (utility == DNA_TO_AA) {
            po::options_description utilDesc (string(STR_DNA_TO_AA) + " options");
            addProcessOptions_DNAToAA(dnaToAA, utilDesc);
            
            cmdline_options.add(utilDesc);
            
            // Collect all the unrecognized options from the first pass. This will include the
            // (positional) mode and command name, so we need to erase them
            vector<string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
            opts.erase(opts.begin());       // erase mode
            opts.erase(opts.begin());       // erase command name
            
            // Parse again...
            po::store(po::command_line_parser(opts).options(utilDesc).run(), vm);
        }
//        else                                                                        // unrecognized utility
//            throw po::invalid_option_value(utility);
        
        // if help specified, print usage message and quit
        if (vm.count("help")) {
            cout << make_usage_string(basename(argv[0]), cmdline_options, pos) << endl;
            return false;
        }
        
        // update all values and make sure required are provided
        po::notify(vm);
    }
    catch (exception &ex) {
        cerr << "Error: " << ex.what() << endl;
        return false;
    }
    
    return true;
    
    
}


void OptionsUtilities::addProcessOptions_ExtractUpstream(ExtractUpstreamUtility &options, po::options_description &processOptions) {
    
    processOptions.add_options()
        ("sequence,s", po::value<string>(&options.fn_sequence)->required(), "Sequence filename")
        ("label,l", po::value<string>(&options.fn_label)->required(), "Label filename")
        ("output,o", po::value<string>(&options.fn_output)->required(), "Output filename")
        ("length", po::value<size_t>(&options.length)->required(), "Upstream length")
        ("allow-overlap-with-cds", "If set, then upstream (non-coding) regions are allowed to overlap with coding regions. If not set, these sequences are ignored.")
        ("min-gene-length", po::value<size_t>(&options.minimumGeneLength)->default_value(0), "Minimum gene length")
    ;
}

void OptionsUtilities::addProcessOptions_ExtractStartContext(ExtractStartContext &options, po::options_description &processOptions) {
    
    processOptions.add_options()
    ("sequence,s", po::value<string>(&options.fn_sequence)->required(), "Sequence filename")
    ("label,l", po::value<string>(&options.fn_label)->required(), "Label filename")
    ("output,o", po::value<string>(&options.fn_output)->required(), "Output filename")
    ("left", po::value<int>(&options.leftRelativeToStart)->required(), "Left index relative to start")
    ("right", po::value<int>(&options.rightRelativeToStart)->required(), "Right index relative to start")
    ("allow-overlap-with-cds", po::bool_switch(&options.allowOverlaps), "If set, contexts can overlap CDS")
    ("min-gene-length", po::value<size_t>(&options.minimumGeneLength), "Minimum gene length associated with start-context")
    ("output-fasta", po::bool_switch(&options.outputFastaDefs)->default_value(false), "If set, output is in FASTA format")
    ;
}



void OptionsUtilities::addProcessOptions_StartModelInfo(StartModelInfoUtility &options, po::options_description &processOptions) {
    
    processOptions.add_options()
        ("sequence,s", po::value<string>(&options.fn_sequence)->required(), "Sequence filename")
        ("label,l", po::value<string>(&options.fn_label)->required(), "Label filename")
        ("num-sim-noncoding", po::value<size_t>(&options.numOfSimNonCoding)->default_value(1000), "Number of simulated non-coding sequences.")
    ;
    
    // gms2 training options
    po::options_description gms2training("GMS2 Training");
    OptionsGMS2Training::addProcessOptions(options.optionsGMS2Training, gms2training);
    
    processOptions.add(gms2training);
}



void OptionsUtilities::addProcessOptions_LabelsSimilarityCheck(LabelsSimilarityCheck &options, po::options_description &processOptions) {
    
    processOptions.add_options()
        (",A", po::value<string>(&options.fn_labelsA)->required(), "Set A labels file")
        (",B", po::value<string>(&options.fn_labelsB)->required(), "Set B labels file")
    ;
    
}


void OptionsUtilities::addProcessOptions_EmitNonCoding(EmitNonCoding &options, po::options_description &processOptions) {
    
    processOptions.add_options()
        ("mod,m", po::value<string>(&options.fn_mod)->required(), "Model file containing non-coding model")
        ("out,o", po::value<string>(&options.fn_out)->required(), "Name of output file containing sequence")
        ("length,l", po::value<NumSequence::size_type>(&options.length)->required(), "Length of generated non-coding sequence")
        ("order", po::value<unsigned>(&options.order)->default_value(numeric_limits<int>::infinity()))
    ;
    
}


void OptionsUtilities::addProcessOptions_CountNumORF(CountNumORF &options, po::options_description &processOptions) {
    processOptions.add_options()
        ("seq,s", po::value<string>(&options.fn_sequence)->required(), "Sequence file")
        ("mod,m", po::value<string>(&options.fn_mod)->required(), "Model file")
        ("print-seq", po::bool_switch(&options.printSeq)->default_value(false), "Print orf sequences and positions")
    ;
    
}


void OptionsUtilities::addProcessOptions_ExtractStartContextPerOperonStatus(ExtractStartContextPerOperonStatus &options, po::options_description &processOptions) {
    processOptions.add_options()
    ("sequence,s", po::value<string>(&options.fn_sequence)->required(), "Sequence filename")
    ("label,l", po::value<string>(&options.fn_label)->required(), "Label filename")
    ("allow-overlap-with-cds", "If set, then upstream (non-coding) regions are allowed to overlap with coding regions. If not set, these sequences are ignored.")
    ("min-gene-length", po::value<size_t>(&options.minimumGeneLength)->default_value(0), "Minimum gene length")
    ("dist-thresh-fgio", po::value<size_t>(&options.distThreshFGIO)->default_value(22), "Minimum distance between genes to be classified as FGIO")
    ("dist-thres-ig", po::value<size_t>(&options.distThreshIG)->default_value(22), "Maximum distance between genes to be classified as IG")
    ;
}


void OptionsUtilities::addProcessOptions_ExtractStartContextPerMotifStatus(ExtractStartContextPerMotifStatus &options, po::options_description &processOptions) {
    processOptions.add_options()
    ("sequence,s", po::value<string>(&options.fn_sequence)->required(), "Sequence filename")
    ("label,l", po::value<string>(&options.fn_label)->required(), "Label filename")
    ("allow-overlap-with-cds", "If set, then upstream (non-coding) regions are allowed to overlap with coding regions. If not set, these sequences are ignored.")
    ("min-gene-length", po::value<size_t>(&options.minimumGeneLength)->default_value(0), "Minimum gene length")
    ("dist-thresh-fgio", po::value<size_t>(&options.distThreshFGIO)->default_value(22), "Minimum distance between genes to be classified as FGIO")
    ("dist-thres-ig", po::value<size_t>(&options.distThreshIG)->default_value(22), "Maximum distance between genes to be classified as IG")
    ("match-to", po::value<string>(&options.matchTo)->default_value("TAAGGAGGTGA"), "Match to 16S tail")
    ("allow-ag-subsitution", po::bool_switch(&options.allowAGSubstitution)->default_value(true), "Allow AG substitution in 16S match.")
    ("match-thresh", po::value<size_t>(&options.matchThresh)->default_value(4), "Match threshold")
    ;
}


void OptionsUtilities::addProcessOptions_ComputeGC(ComputeGC &options, po::options_description &processOptions) {
    processOptions.add_options()
    ("sequence,s", po::value<string>(&options.fn_sequence)->required(), "Sequence filename")
    ("label,l", po::value<string>(&options.fn_label)->default_value(""), "Label filename (for genes-based GC)")
    ("per-gene", po::bool_switch(&options.perGene)->default_value(false), "If set, per-gene GC is printed")
    ;
}

void OptionsUtilities::addProcessOptions_SeparateFGIOAndIG(SeparateFGIOAndIG &options, po::options_description &processOptions) {
    processOptions.add_options()
    ("label,l", po::value<string>(&options.fn_label)->required(), "Label filename")
    ("fgio", po::value<string>(&options.fnout_fgio)->default_value(""), "Filename for FGIO labels")
    ("ig", po::value<string>(&options.fnout_ig)->default_value(""), "Filename for IG labels")
    ("dist-thresh-fgio", po::value<size_t>(&options.distThreshFGIO)->default_value(25), "Minimum distance between genes to be classified as FGIO")
    ("dist-thres-ig", po::value<size_t>(&options.distThreshIG)->default_value(22), "Maximum distance between genes to be classified as IG")
    ;
}

void OptionsUtilities::addProcessOptions_DNAToAA(DNAToAA &options, po::options_description &processOptions) {
    processOptions.add_options()
    ("sequence,s", po::value<string>(&options.fnseqs)->required(), "Sequence filename")
    ("output-fasta", po::bool_switch(&options.outputFastaDefs)->default_value(false), "If set, output is in FASTA format")
    ;
}


