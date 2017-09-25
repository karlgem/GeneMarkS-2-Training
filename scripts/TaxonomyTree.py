from collections import defaultdict
import pdb

class TaxonomyTree:

    taxonomyInfo = None;
    root = None;
    group = None;

    def __init__(self, taxonomyInfo, group):
        """ Initialize the taxonomy tree

            @param taxonomyInfo a dict of dict where the key
            is a taxid, and the value is anoter dictionary whose keys are
            the parent taxid ('parent'), name of genus species ('name')
         """

        self.group = group;
        self.taxonomyInfo = taxonomyInfo;
        self.root = self.BuildTree(taxonomyInfo);
        
    def BuildTree(self, taxonomy):

        root = self.NewNode("1", "1", taxonomy["1"]["name"]);

        # for every specie
        counter = 0;
        totalTaxIds = len(taxonomy);

        for currentTaxid in taxonomy:

            # get its ancestry
            ancestry = self.GetFullAncestry(currentTaxid, taxonomy);

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
                    tmpNode["children"][tmpTaxid] = self.NewNode(tmpTaxid, tmpNode["taxid"], tmpName);
                    tmpNode = tmpNode["children"][tmpTaxid];

        return root;

    # create a new node with default values and taxid, parent, and name
    def NewNode(self, taxid, parent, name):
        aNode = dict(); #defaultdict(dict);
        aNode["taxid"] = taxid;
        aNode["name"] = name;
        aNode["parent"] = parent;
        aNode["children"] = defaultdict(dict);
        aNode["total"] = 0;
        aNode["total-of-type"] = 0;


        if self.group == "A":
            aNode["sum-of-percent-of-leaderless-in-fgio-in-prediction"]   = 0;
            aNode["sum-of-percent-of-leaderless-in-all-genes-in-prediction"] = 0;

            aNode["sum-of-percent-of-fgio-in-training"] = 0;
            aNode["sum-of-percent-of-fgio-in-prediction"] = 0

        elif self.group == "B":
            aNode["sum-of-percent-of-leaderless-in-fgio-in-prediction"] = 0;
            aNode["sum-of-percent-of-leaderless-in-fgio-in-training"]   = 0;

            aNode["sum-of-percent-of-fgio-in-training"]   = 0;
            aNode["sum-of-percent-of-fgio-in-prediction"] = 0;

            aNode["sum-of-percent-of-leaderless-in-all-genes-in-training"]   = 0;
            aNode["sum-of-percent-of-leaderless-in-all-genes-in-prediction"] = 0;

        return aNode;

    ##### Get ancestry of taxid's: root to current
    def GetFullAncestry(self, taxid, taxonomy):

        ancestry = [taxid];

        currentTaxid = taxid;

        while currentTaxid != "1":
            current = taxonomy[currentTaxid];
            parentTaxid = current["parent"];

            ancestry.insert(0, parentTaxid);

            currentTaxid = parentTaxid;

        return ancestry;

    def UpdateTree(self, genome):

        genomeTaxid = genome["taxid"];
        ancestry = self.GetFullAncestry(genomeTaxid, self.taxonomyInfo);

        # go through tree using ancestry
        currNode = self.root
        for currTaxid in ancestry:

            # if root, don't loop further
            if currTaxid == "1":
                self.UpdateNode(currNode, genome);
                continue;

            if currTaxid in currNode["children"]:
                currNode = currNode["children"][currTaxid];

                # if final child, get consensus sequences
                if currTaxid == ancestry[-1]:
                    if "consensus" in genome:
                        currNode["consensus"] = genome["consensus"];
                    if "consensus-rbs" in genome:
                        currNode["consensus-rbs"] = genome["consensus-rbs"];

                self.UpdateNode(currNode, genome);

    def UpdateNode(self, node, genome):

        node['total'] += 1;

        genomeGroup = genome["genome-type"];

        # make sure genome group matches group represented by the tree
        if genomeGroup == self.group:

            node['total-of-type'] += 1;

            # common options for group-a AND group-b
            if genomeGroup == "A" or genomeGroup == "B":
                
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
            if genomeGroup == "B":

                # Training: leaderless / FGIO
                if genome["total-genes-in-training"] > 0:
                    currNode["sum-of-percent-of-leaderless-in-fgio-in-training"] += 100 * float(genome["total-leaderless-genes-in-fgio-in-training"]) / float(genome["fgio-in-training"])

                # Training: leaderless / all
                if genome["total-genes-in-prediction"] > 0:
                    currNode["sum-of-percent-of-leaderless-in-all-genes-in-training"] += 100 * float(genome["total-leaderless-genes-in-all-genes-in-training"]) / float(genome["total-genes-in-training"])

    def PrintTree(self):

        self.PrintTreeRecursive(self.root);

    def PrintTreeRecursive(self, currentNode, depth=None):
        if depth == None:
            depth = 0;

        # pdb.set_trace()

        if currentNode == None:
            return;

        if currentNode["total"] == 0:
            return;

        # add current level
        print self.NewLevel(depth, currentNode);

        # loop over children sorted by total-of-type
        for childTaxID in sorted(currentNode["children"].iteritems(), key=lambda (k,v) : v["total-of-type"], reverse=True):
            self.PrintTreeRecursive(currentNode["children"][childTaxID[0]], depth+1)



    def NewLevel(self, depth, datapoint):
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
        if self.group == 'A' or self.group == 'B':

            # sum of percentage of fgio
            percentOfFGIOInTraining   = 0;
            percentOfFGIOInPrediction = 0;
            if datapoint['total-of-type'] > 0:
                percentOfFGIOInTraining   = round(float(datapoint["sum-of-percent-of-fgio-in-training"])   / float(datapoint["total-of-type"]) , 1);
                percentOfFGIOInPrediction = round(float(datapoint["sum-of-percent-of-fgio-in-prediction"]) / float(datapoint["total-of-type"]) , 1);

            # result += "\t" + str(percentOfFGIOInTraining);
            result += "\t" + str(percentOfFGIOInPrediction);

            # percentage of leaderless in fgio
            if self.group == "B":
                percentOfLeaderlessInFGIOInTraining = 0;
                if float(datapoint['total-of-type']) > 0:
                    percentOfLeaderlessInFGIOInTraining = round(float(datapoint["sum-of-percent-of-leaderless-in-fgio-in-training"]) / float(datapoint['total-of-type']), 1);

                # result += "\t" + str(percentOfLeaderlessInFGIOInTraining);

            percentOfLeaderlessInFGIOInPrediction = 0;
            if float(datapoint["total-of-type"]) > 0:
                percentOfLeaderlessInFGIOInPrediction = round(float(datapoint["sum-of-percent-of-leaderless-in-fgio-in-prediction"]) / float(datapoint["total-of-type"]), 1);

            result += "\t" + str(percentOfLeaderlessInFGIOInPrediction);


            # percentage of leaderless in all genes
            if self.group == "B":
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

























