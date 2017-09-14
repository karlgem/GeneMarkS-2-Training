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
from TaxonomyTree import TaxonomyTree

##### Command line arguments
parser = argparse.ArgumentParser(description='Find number of leaderless in taxonomy.');
parser.add_argument('list', type=str, help='Genome list');
parser.add_argument('taxonomy', type=str, help="Taxonomy information");
parser.add_argument('datapath', type=str, help="Path to data folder");
parser.add_argument('names', type=str, help="File containing taxid names");
parser.add_argument('--group', type=str, choices=['A', 'B', 'C', 'D', 'E', 'C2'], help="Genome group");
parser.add_argument('--run-dir', type=str, default="gms2", help="GMS2 directory");
parser.add_argument('--verbose', action='store_true', default=False, dest='verbose');

args = parser.parse_args();

groupTagPair = {
    'A' : 'group-a',
    'B' : 'group-b',
    'C' : 'group-c',
    'D' : 'group-d',
    'E' : 'group-e',
    'C2': 'group-c2'
}

groupQuery = groupTagPair[args.group];


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
                # if (args.group == "A" and genomeType == "group-a") or (args.group == "B" and genomeType == "group-b") or (groupQuery == "group-c" and genomeType == "group-c2"):
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


            file.close();

        if genomeType != "":
            genome["genome-type"] = genomeType

        if numLeaderlessGenes != "":
            genome["total-leaderless-genes-in-fgio-in-training"] = numLeaderlessGenes;
            genome["total-leaderless-genes-in-all-genes-in-training"] = numLeaderlessGenes;

        if numFGIO != "":
            genome["fgio-in-training"] = numFGIO;

        if groupQuery == genomeType and consensusPromoter != "":
            genome["consensus-promoter"] = consensusPromoter;

        if groupQuery == genomeType and consensusRBS != "":
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
                    elif groupQuery == "group-a" and currGene['model'] == "native":
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
                elif groupQuery == "group-a" and currGene['model'] == "native":
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

taxtree = TaxonomyTree(taxonomy, args.group)    


# update taxonomy tree with leaderless information
if args.verbose:
    print "Updating record tree... "
for genome in genomeList:
    taxtree.UpdateTree(genome);

if args.verbose:
    print "Printing record tree... "
taxtree.PrintTree();



