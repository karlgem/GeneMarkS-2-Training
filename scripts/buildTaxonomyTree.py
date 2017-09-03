# Author: Karl Gemayel <karl@gatech.edu>

import argparse
from collections import defaultdict
import matplotlib
matplotlib.use('Agg');
import matplotlib.pyplot as plt
import numpy as np
import pdb
import re
import subprocess
import PIL
import os.path
import os

##### Command line arguments
parser = argparse.ArgumentParser(description='Find number of leaderless in taxonomy.');
parser.add_argument('list', type=str, help='Genome list');
parser.add_argument('taxonomy', type=str, help="Taxonomy information");
parser.add_argument('datapath', type=str, help="Path to data folder");
parser.add_argument('names', type=str, help="File containing taxid names");
parser.add_argument('--gtype', type=str, choices=['group-' + x for x in ['a','b','c','d','e', 'a2']], help="Genome type");
parser.add_argument('--run-dir', type=str, default="gms2", help="GMS2 directory");
parser.add_argument('--type-from-class-file', action='store_true', default=False, dest='typeFromClassFile');
parser.add_argument('--class-file', type=str, default="class", help="name of class file");
parser.add_argument('--verbose', action='store_true', default=False, dest='verbose');

args = parser.parse_args();

# groupTagPair = {
#     'group-a' : 'archaea-promoter',
#     'group-b' : 'bacteria-promoter',
#     'group-c' : 'class-c',
#     'group-d' : 'pure-rbs',
#     'group-e' : 'upstream-signature',
#     'group-a2': 'archaea-promoter-2'
# }
groupTagPair = {
    'group-a' : 'group-a',
    'group-b' : 'group-b',
    'group-c' : 'group-c',
    'group-d' : 'group-d',
    'group-e' : 'group-e',
    'group-a2': 'group-a'
}


##### Read input file
def readGenomeList(fname):

    data = list();

    f = open(fname, "r");

    species = "";
    for line in f:
        line = line.strip().split();

        d = {"gcfid" : line[0], "gcode" : line[1], "taxid" : line[2] };
        data += [d];

    f.close();

    return data;


##### Read taxonomy information
def readTaxonomyFile(fname):

    data = defaultdict(dict);

    f = open(fname, "r");

    for line in f:
        line = line.strip().split('|');

        d = {"parent" : line[1].strip(), "type" : line[2].strip() };

        data[line[0].strip()] = d;

    f.close();

    return data;



def GetFullAncestry(taxid, taxonomy):

    ancestry = [taxid];

    currentTaxid = taxid;

    while currentTaxid != "1":
        current = taxonomy[currentTaxid];
        parentTaxid = current["parent"];

        ancestry.insert(0, parentTaxid);

        currentTaxid = parentTaxid;

    return ancestry;


def BuildRecordTree(taxonomy):

    def newNode(taxid, parent, name):
        aNode = dict(); #defaultdict(dict);
        aNode["taxid"] = taxid;
        aNode["name"] = name;
        aNode["parent"] = parent;
        aNode["children"] = defaultdict(dict);
        aNode["total"] = 0;
        aNode["total-of-type"] = 0;


        if args.gtype == "group-a" or args.gtype == "group-a2":
            aNode["sum-of-percent-of-leaderless-in-fgio-in-prediction"]   = 0;
            aNode["sum-of-percent-of-leaderless-in-all-genes-in-prediction"] = 0;

            aNode["sum-of-percent-of-fgio-in-training"] = 0;
            aNode["sum-of-percent-of-fgio-in-prediction"] = 0

        elif args.gtype == "group-b":
            aNode["sum-of-percent-of-leaderless-in-fgio-in-prediction"] = 0;
            aNode["sum-of-percent-of-leaderless-in-fgio-in-training"]   = 0;

            aNode["sum-of-percent-of-fgio-in-training"]   = 0;
            aNode["sum-of-percent-of-fgio-in-prediction"] = 0;

            aNode["sum-of-percent-of-leaderless-in-all-genes-in-training"]   = 0;
            aNode["sum-of-percent-of-leaderless-in-all-genes-in-prediction"] = 0;



        return aNode;

    root = newNode("1", "1", taxonomy["1"]["name"]);

    # for every specie
    counter = 0;
    totalTaxIds = len(taxonomy);

    for currentTaxid in taxonomy:

        counter += 1;
        if args.verbose:
            out = str(counter) + "/" + str(totalTaxIds);
            if counter == totalTaxIds:
                print out
            else:
                print out + "\r",

        # get its ancestry
        ancestry = GetFullAncestry(currentTaxid, taxonomy);

        # fill in the tree
        tmpNode = root;
        for tmpTaxid  in ancestry:
            tmpName = taxonomy[tmpTaxid]["name"];
            # skip root
            if tmpTaxid == "1":
                continue;

            # go to child if it exists
            if tmpTaxid in tmpNode["children"]:
                tmpNode = tmpNode["children"][tmpTaxid];
            # otherwise, create new node
            else:
                tmpNode["children"][tmpTaxid] = newNode(tmpTaxid, tmpNode["taxid"], tmpName);
                tmpNode = tmpNode["children"][tmpTaxid];

    return root;


def GetGenomeInfo(genomeList, dataPath):

    current = 0;
    totalGenomes = len(genomeList);

    for genome in genomeList:
        current += 1;

        if args.verbose:
            out = str(current) + "/" + str(totalGenomes);
            if current == totalGenomes:
                print out
            else:
                print out + "\r",

        # get path to mod file
        pathToGenome = os.path.join(dataPath, genome["gcfid"]);
        pathToGenomeRun = os.path.join(pathToGenome, "runs", args.run_dir);
        pathToMod = os.path.join(pathToGenomeRun, "GMS2.mod");

        cmd = "ls " + str(pathToGenomeRun) + "/itr_*.lst | sort -t _ -k 2 -n | tail -n 2 | head -n 1";
        # print cmd
        # proc = subprocess.Popen("pwd", stdout=subprocess.PIPE)
        # trainingLST = proc.stdout.read()
        trainingLST = os.popen(cmd).read().strip();
        # print trainingLST

        genomeType = "";
        numLeaderlessGenes = "0";
        numFGIO = "0";
        consensusPromoter = "";
        consensusRBS = "";

        # check if Mod file exists
        if os.path.isfile(pathToMod):

            file = open(pathToMod, "r");

            # for line in file:
            while True:
                line = file.readline()
                if not line:
                    break;

                # Genome type
                if re.search("GENOME_TYPE", line):
                    genomeType = line.strip().split()[1];

                # Number of leaderless genes
                if re.search("NUM_LEADERLESS", line):
                    numLeaderlessGenes = line.strip().split()[1];

                # Number of first genes in operon
                if re.search("NUM_FGIO", line):
                    numFGIO = line.strip().split()[1];

                # collect promoter matrix
                if (args.gtype == "group-a" and genomeType == groupTagPair["group-a"]) or (args.gtype == "group-b" and genomeType == groupTagPair["group-b"]) or (args.gtype == "group-a2" and genomeType == groupTagPair["group-a2"]):
                    if  re.search("\$PROMOTER_MAT",  line):
                        lineA = file.readline(); lineA = lineA.strip().split();
                        lineC = file.readline(); lineC = lineC.strip().split();
                        lineG = file.readline(); lineG = lineG.strip().split();
                        lineT = file.readline(); lineT = lineT.strip().split();

                        for i in range(1, len(lineA)):
                            currMax = float(lineA[i]);
                            letter = "A";
                            if currMax <= float(lineC[i]):
                                currMax = float(lineC[i]);
                                letter = "C";
                            if currMax <= float(lineG[i]):
                                currMax = float(lineG[i]);
                                letter = "G";
                            if currMax <= float(lineT[i]):
                                currMax = float(lineT[i]);
                                letter = "T";

                            consensusPromoter += letter;

                # get RBS consensus
                if re.search("\$RBS_MAT", line):
                    lineA = file.readline(); lineA = lineA.strip().split();
                    lineC = file.readline(); lineC = lineC.strip().split();
                    lineG = file.readline(); lineG = lineG.strip().split();
                    lineT = file.readline(); lineT = lineT.strip().split();

                    for i in range(1, len(lineA)):
                        currMax = float(lineA[i]);
                        letter = "A";
                        if currMax <= float(lineC[i]):
                            currMax = float(lineC[i]);
                            letter = "C";
                        if currMax <= float(lineG[i]):
                            currMax = float(lineG[i]);
                            letter = "G";
                        if currMax <= float(lineT[i]):
                            currMax = float(lineT[i]);
                            letter = "T";

                        consensusRBS += letter;

                    #if consensus[0] == "T" and consensus[1] == "A" and consensus[5] == "T":
                    #    consensus = "TANNNT";


            file.close();

        if genomeType != "":
            if genomeType == "archaea-promoter-2":
                genomeType = "archaea-promoter";
            genome["genome-type"] = genomeType

        if args.typeFromClassFile:
            pathToClassFile = os.path.join(dataPath, genome["gcfid"], "runs", args.run_dir, args.class_file);
            file = open(pathToClassFile, "r");
            answer = file.readline().strip();
            file.close();

            if answer == "yes":
                genome["genome-type"] = "bacteria-promoter";
            elif answer == "no":
                genome["genome-type"] = "other";

        if numLeaderlessGenes != "":
            genome["total-leaderless-genes-in-fgio-in-training"] = numLeaderlessGenes;
            genome["total-leaderless-genes-in-all-genes-in-training"] = numLeaderlessGenes;

        if numFGIO != "":
            genome["fgio-in-training"] = numFGIO;

        if consensusPromoter != "":
            genome["consensus-promoter"] = consensusPromoter;

        if args.gtype == genomeType and consensusRBS != "":
            genome["consensus-rbs"] = consensusRBS


        # count total genes in training
        def getTotalGenes(fname):

            file = open(fname, "r");

            totalGenes = 0;
            for line in file:

                currGene = {};      # empty dictionary

                m = re.search("\s*\d+\s+([+-])\s+(\d+)\s+(\d+)\s+(\d+)\s+(atypical|native)", line)
                if m:
                    currGene['length'] = int(m.group(4));
                    currGene['model']  = m.group(5);

                    # ignore short genes
                    if currGene['length'] > 300:
                        totalGenes += 1
                    # allow short native genes for group A
                    elif (args.gtype == "group-a" or args.gtype == "group-a2") and currGene['model'] == "native":
                        totalGenes += 1

            file.close();

            return totalGenes;

        totalGenes = getTotalGenes(trainingLST);
        genome["total-genes-in-training"] = totalGenes;


def GetTotalNumberOfGenes(genomeList, dataPath):

    def readGeneFile(fname):

        file = open(fname, "r");

        # read file into list
        genesList = list();
        for line in file:

            currGene = {};      # empty dictionary

            m = re.search("\s*\d+\s+([+-])\s+(\d+)\s+(\d+)\s+(\d+)\s+(atypical|native)", line)
            if m:
                currGene['strand'] = m.group(1);
                currGene['left']   = int(m.group(2));
                currGene['right']  = int(m.group(3));
                currGene['length'] = int(m.group(4));
                currGene['model']  = m.group(5);

                # get start model type
                m2 = re.search("atypical|native\s+[ACGT]+\s+\d+\s+(\d+)", line)
                if m2:
                    currGene['start-model-type'] = m2.group(1)

                # ignore short genes
                if currGene['length'] > 300:
                    genesList += [currGene];
                # allow short native genes for group A
                elif (args.gtype == "group-a" or args.gtype == "group-a2") and currGene['model'] == "native":
                    genesList += [currGene];

        file.close();

        return genesList;

    current = 0;
    totalGenomes = len(genomeList);
    for genome in genomeList:

        current += 1;

        if args.verbose:
            out = str(current) + "/" + str(totalGenomes);
            if current == totalGenomes:
                print out
            else:
                print out + "\r",

        # get path to LST file
        pathToLST = os.path.join(dataPath, genome["gcfid"], "runs", args.run_dir, "gms2.lst");

        genomeType = "";

        totalGenes = 0;
        totalPredictedLeaderlessInFGIO = 0;
        totalPredictedLeaderlessInAllGenes = 0;
        FGIOInPrediction = 0;

        prevRight = 0;
        prevLeft = 0;
        prevStrand = "";

        # check if file exists
        if os.path.isfile(pathToLST):

            # read genes
            genesList = readGeneFile(pathToLST);

            # go through list and gather total number of genes, genes predicted
            # as leaderless, and FGIO

            totalGenes = len(genesList);

            for n in range(len(genesList)):

                currGene = genesList[n];

                ##### Get FGIO

                isFGIO = False;

                # positive strand
                if currGene['strand'] == "+":
                    if n == 0:
                        isFGIO = True
                    else:
                        prevGene = genesList[n-1];

                        # if previous gene on negative strand
                        if prevGene['strand'] == '-':
                            isFGIO = True;      # gene is FGIO
                        else:
                            # if on same strand but too far away
                            if currGene['left'] - prevGene['right'] > 25:
                                isFGIO = True;
                elif currGene['strand'] == "-":
                    if n == totalGenes-1:
                        isFGIO = True;
                    else:
                        prevGene = genesList[n+1];

                        # if upstream gene on negative strand
                        if prevGene['strand'] == '+':
                            isFGIO = True;      # gene is FGIO
                        else:
                            # if on same strand but too far away
                            if prevGene['left'] - currGene['right'] > 25:
                                isFGIO = True;
                else:
                    print "Warning: invalid strand '" + currGene['strand'] + "'"

                if isFGIO == True:
                    FGIOInPrediction += 1;

                    ##### Get leaderless
                    if 'start-model-type' in currGene and currGene['start-model-type'] == '2':
                        totalPredictedLeaderlessInFGIO += 1;

                if 'start-model-type' in currGene and currGene['start-model-type'] == '2':
                    totalPredictedLeaderlessInAllGenes += 1;

        genome["total-genes-in-prediction"] = totalGenes;
        genome['fgio-in-prediction'] = FGIOInPrediction;
        # print genome['gcfid'], FGIOInPrediction, totalGenes, FGIOInPrediction / float(totalGenes)
        genome['total-leaderless-genes-in-fgio-in-prediction'] = totalPredictedLeaderlessInFGIO;
        genome['total-leaderless-genes-in-all-genes-in-prediction'] = totalPredictedLeaderlessInAllGenes;



def GetTaxidNames(taxonomy, fnNames):

    # open file containing taxid to names mapping
    f = open(fnNames, "r")

    for line in f:
        line = line.strip().split('|');
        taxid = line[0].strip();
        name = line[1].strip();
        nameClass = line[3].strip();

        if taxid in taxonomy:
            # if taxid name already assigned, skip it
            if "name" not in taxonomy[taxid]:
                taxonomy[taxid]["name"] = name;
            # if taxid already assigned but new name has class "scientific name", then assign it
            elif nameClass == "scientific name":
                taxonomy[taxid]["name"] = name;

    f.close();




def BFS(root):

    queue = [root];

    while queue:
        vertex = queue.pop(0);
        print vertex["taxid"];

        for neighborTaxid in vertex["children"]:
            queue.append(vertex["children"][neighborTaxid]);




def DFS(root, path=None):

    if path is None:
        path = list();

    path.append(root["name"]);

    for neighborTaxid in root["children"]:
        DFS(root["children"][neighborTaxid], path)

    # print leaf
    if len(root["children"]) >= 0:
        percentLeaderless = 0;
        if root["total"] > 0:
            percentLeaderless = 100 * float(root["num leaderless"]) / float(root["total"]);

            unrollPath = "";
            for x in path:
                unrollPath += "\"" + x + "\"";
            print unrollPath + "\"" + str(root["num leaderless"]) + "\"" + str(root["total"]) + "\"" + str(percentLeaderless);


    del path[-1];


# if args.gtype == "group-a":
#             aNode["sum-of-percent-of-leaderless-in-fgio-in-prediction"]   = 0;

#             aNode["sum-of-percent-of-fgio-in-prediction"] = 0

#         elif args.gtype == "group-b":
#             aNode["sum-of-percent-of-leaderless-in-fgio-in-prediction"] = 0;
#             aNode["sum-of-percent-of-leaderless-in-fgio-in-training"]   = 0;

#             aNode["sum-of-percent-of-fgio-in-training"]   = 0;
#             aNode["sum-of-percent-of-fgio-in-prediction"] = 0;

#             aNode["sum-of-percent-of-leaderless-in-all-genes-training"]   = 0;
#             aNode["sum-of-percent-of-leaderless-in-all-genes-prediction"] = 0;

def UpdateRecordTree(genomeList, taxonomyTree, taxonomy):

    def updateNodeWithGenomeInfo(node, genome):

        node['total'] += 1;

        if genome['genome-type'] == groupTagPair[args.gtype]:

            node['total-of-type'] += 1;

            # common options for group-a AND group-b
            if genome['genome-type'] == groupTagPair["group-a"] or genome['genome-type'] == groupTagPair["group-b"] or genome['genome-type'] == groupTagPair["group-a2"] :
                # Training: fgio / all
                if "fgio-in-training" in genome and genome["total-genes-in-prediction"] > 0:
                    currNode["sum-of-percent-of-fgio-in-training"] += 100 * float(genome["fgio-in-training"]) / float(genome["total-genes-in-training"]);

                # Prediction: fgio / all
                if "fgio-in-prediction" in genome and float(genome["fgio-in-prediction"]) > 0:
                    currNode["sum-of-percent-of-fgio-in-prediction"] += 100*float(genome["fgio-in-prediction"]) / float(genome["total-genes-in-prediction"])

                # Prediction: leaderless / FGIO
                if genome["fgio-in-prediction"] > 0:
                    currNode["sum-of-percent-of-leaderless-in-fgio-in-prediction"] += 100 * float(genome["total-leaderless-genes-in-fgio-in-prediction"]) / float(genome["fgio-in-prediction"])

                # Prediction: leaderless / all
                if genome["total-genes-in-prediction"] > 0:
                    currNode["sum-of-percent-of-leaderless-in-all-genes-in-prediction"] += 100 * float(genome["total-leaderless-genes-in-all-genes-in-prediction"]) / float(genome["total-genes-in-prediction"])


            # group-b options only
            if genome['genome-type'] == groupTagPair["group-b"]:

                # Training: leaderless / FGIO
                if genome["total-genes-in-training"] > 0:
                    currNode["sum-of-percent-of-leaderless-in-fgio-in-training"] += 100 * float(genome["total-leaderless-genes-in-fgio-in-training"]) / float(genome["fgio-in-training"])

                # Training: leaderless / all
                if genome["total-genes-in-prediction"] > 0:
                    currNode["sum-of-percent-of-leaderless-in-all-genes-in-training"] += 100 * float(genome["total-leaderless-genes-in-all-genes-in-training"]) / float(genome["total-genes-in-training"])



    current = 0;
    totalGenomes = len(genomeList);

    # pdb.set_trace()

    for genome in genomeList:

        current += 1;

        if args.verbose:
            out = str(current) + "/" + str(totalGenomes);
            if current == totalGenomes:
                print out
            else:
                print out + "\r",

        # pdb.set_trace()

        genomeTaxid = genome["taxid"];
        ancestry = GetFullAncestry(genomeTaxid, taxonomy);

        if "genome-type" not in genome:
            continue;

        isBacteriaLeaderless = (genome["genome-type"] == args.gtype);
        isClass= (genome["genome-type"] == args.gtype);

        # go through tree using ancestry
        currNode = taxonomyTree
        for currTaxid in ancestry:

            # if root, don't loop further
            if currTaxid == "1":
                updateNodeWithGenomeInfo(currNode, genome);
                continue;

            if currTaxid in currNode["children"]:
                currNode = currNode["children"][currTaxid];

                # if final child, get consensus sequences
                if currTaxid == ancestry[-1]:
                    if "consensus" in genome:
                        currNode["consensus"] = genome["consensus"];
                    if "consensus-rbs" in genome:
                        currNode["consensus-rbs"] = genome["consensus-rbs"];

                updateNodeWithGenomeInfo(currNode, genome);


def PrintTree(root, depth=None):

    def addNewLevel(depth, datapoint):
        """ Draw a new level in the tree with the given depth and datapoint"""
        oneLevel = "    |";
        depthLevel = oneLevel * depth;

        startPositionOfNumbers = 120;

        # add name and spaces to result
        result = "";
        if depth > 0:
            result += depthLevel + "__ ";
        result += datapoint["name"] + (" " * (startPositionOfNumbers - len(depthLevel) - len(datapoint["name"]) - len("__ ")))

        # Percent of type
        percentOfType = 0;
        if datapoint["total"] > 0:
            percentOfType = 100 * datapoint["total-of-type"] / float(datapoint["total"]);
            percentOfType = round(percentOfType, 1)


        # Add total, total-of-type, and percent-of-type
        result += str(datapoint["total"]) + "\t" + str(datapoint["total-of-type"]) + "\t" + str(percentOfType);

        # add info for classes A and B
        if args.gtype == 'group-a' or args.gtype == 'group-b' or args.gtype == "group-a2":

            # sum of percentage of fgio
            percentOfFGIOInTraining   = 0;
            percentOfFGIOInPrediction = 0;
            if datapoint['total-of-type'] > 0:
                percentOfFGIOInTraining   = round(float(datapoint["sum-of-percent-of-fgio-in-training"])   / float(datapoint["total-of-type"]) , 1);
                percentOfFGIOInPrediction = round(float(datapoint["sum-of-percent-of-fgio-in-prediction"]) / float(datapoint["total-of-type"]) , 1);

            # result += "\t" + str(percentOfFGIOInTraining);
            result += "\t" + str(percentOfFGIOInPrediction);

            # percentage of leaderless in fgio
            if args.gtype == "group-b":
                percentOfLeaderlessInFGIOInTraining = 0;
                if float(datapoint['total-of-type']) > 0:
                    percentOfLeaderlessInFGIOInTraining = round(float(datapoint["sum-of-percent-of-leaderless-in-fgio-in-training"]) / float(datapoint['total-of-type']), 1);

                # result += "\t" + str(percentOfLeaderlessInFGIOInTraining);

            percentOfLeaderlessInFGIOInPrediction = 0;
            if float(datapoint["total-of-type"]) > 0:
                percentOfLeaderlessInFGIOInPrediction = round(float(datapoint["sum-of-percent-of-leaderless-in-fgio-in-prediction"]) / float(datapoint["total-of-type"]), 1);

            result += "\t" + str(percentOfLeaderlessInFGIOInPrediction);


            # percentage of leaderless in all genes
            if args.gtype == "group-b":
                percentOfLeaderlessInAllGenesInTraining = 0;
                if float(datapoint['total-of-type']) > 0:
                    percentOfLeaderlessInAllGenesInTraining = round(float(datapoint["sum-of-percent-of-leaderless-in-all-genes-in-training"]) / float(datapoint['total-of-type']), 1);

                # result += "\t" + str(percentOfLeaderlessInAllGenesInTraining);

            percentOfLeaderlessInAllGenesInPrediction = 0;
            if float(datapoint["total-of-type"]) > 0:
                percentOfLeaderlessInAllGenesInPrediction = round(float(datapoint["sum-of-percent-of-leaderless-in-all-genes-in-prediction"]) / float(datapoint["total-of-type"]), 1);

            result += "\t" + str(percentOfLeaderlessInAllGenesInPrediction);


        # result += "\t" + str(round(average,1));

        if "consensus-rbs" in datapoint:
            result += "\t" + datapoint["consensus-rbs"];

        if "consensus-promoter" in datapoint:
            result += "\t" + datapoint["consensus-promoter"];

        return result;

    def stopExpansion(datapoint, depth=None):
#        if datapoint["total"] <=5:
#            return True;
        if datapoint["total-of-type"] == 0:
            return True;

        if depth != None:
            if depth >=2000:
                return True;

        return False

    if depth is None:
        depth = 0;

    if root["total"] == 0:
        return;

    print addNewLevel(depth, root);

    #pdb.set_trace();
    #for neighborTaxid in root["children"]:
    #for neighborTaxid in sorted(root["children"].iteritems(), key=lambda (k,v) : (v,k)):
    for neighborTaxid in sorted(root["children"].iteritems(), key=lambda (k,v) : v["total-of-type"], reverse=True):
        if not stopExpansion(root["children"][neighborTaxid[0]], depth):
            PrintTree(root["children"][neighborTaxid[0]], depth+1)




if args.verbose:
    print "Reading genome list... "

genomeList = readGenomeList(args.list);

# get genome types
if args.verbose:
    print "Reading genome info... "
GetGenomeInfo(genomeList, args.datapath);

if args.verbose:
    print "Reading gene list... "
GetTotalNumberOfGenes(genomeList, args.datapath);

# read taxonomy information
if args.verbose:
    print "Reading taxonomy file... "
taxonomy = readTaxonomyFile(args.taxonomy);

# add name information to taxonomy
if args.verbose:
    print "Reading taxid names... "
GetTaxidNames(taxonomy, args.names);

# build taxonomy tree
if args.verbose:
    print "Building empty record tree... "
root = BuildRecordTree(taxonomy);


# update taxonomy tree with leaderless information
if args.verbose:
    print "Updating record tree... "
UpdateRecordTree(genomeList, root, taxonomy);

#DFS(root);

if args.verbose:
    print "Printing record tree... "
PrintTree(root);



