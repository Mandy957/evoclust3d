####################################################################################################
#                                                                                                  #
#           PROJECT Protein Adaptation                                                             #
#           CLASS ScopeAlgorithmTreeSet                                                            #
#           PROGRAMMER Jeremy Adams                                                                #
#           STARTED 25-07-13                                                                       #
#           LASTMOD 06-11-14                                                                       #
#                                                                                                  #
#           DESCRIPTION class representation of all branch segments to be used for statistical     #
#                       analysis. Main stats and P-value generating module                         #
#                                                                                                  #
####################################################################################################

from scopealgorithmtrios import ScopeAlgorithmTrios
from pdbatom import PDBAtom
from fastmltree import FastMLTree
from fasequence import FASequence
from numpy import array

import os
import re
import math
import sets
import random
import numpy
from scipy.stats import percentileofscore
from staticmethods import *


class ScopeAlgorithmTriosSet:
    
    "CONSTRUCTOR"
    def __init__(self , DataDIR , reportPATH , scoringPATH , outPATH):
        
        """
        Class attributes:
        DataDIR (String): Directory where main program is held
        Directory (String): Directory to read in report files and output the final PValue file
        ProteinFamilyName (String): Protein family descriptor to use in random distribution file generation
        
        ScopeXMLFile (String): Path to report (mutation mapping) file
        ModdedTreeFile (String): Path to newick syntax tree file with branch names according to first module
        ScoringMatrixXMLFile (String): Path to XML format file of all scoring keys (to be used in random distributions)
        DistributionPath (String): Path to the random distribution directory
        
        HydroPATH (String): Path to hydropathyindex file
        MassPATH (String): Path to sidechainmass file
        
        Hydro_D (Dict): Key is one letter AA code, Value is its hydropathy index
        Mass_D (Dict): Key is one letter AA code, Value is its side chain mass value
        
        Tree (FastMLTree obj): Tree object with the renamed branches
        NodeSequenceKey_L (List): List containing the names of all nodes
        NodeToSequence_D (Dict): Key is the node name, Value is the ancestral or extant sequence at that node
        
        BranchToAlgorithm_D (Dict): Key is the Branch key name, Value is a ScopeAlgorithm instance for that branch segment
        PDBContents_D (Dict): Key is PDB ID, value is dictionaries of all atom and residue information for that PDB file
        PDBXMLContents_D (Dict): Key is PDB ID, value is dictionaries of all atom and residue information for that PDBXML file
        
        ScoringMatrixCoverageKeys_D (Dict): Key is PDB ID, value is list of chain and position keys that correspond to successfully aligned regions
        AccsToMutationCount_D (Dict): Key is PDB ID, value is a list of integers for all branch segments with those number of mutations
        AccsToDistanceCount_D (Dict): Key is PDB ID, value is a list of integers for all branch segments with those number of mutations that can be joined by pairwise distances
        
        MassChanges_L (List): all mass change calculations that have happened anywhere on the tree
        HydroChanges_L (List): all hydropathy index change calculations that have happened anywhere on the tree
        
        RandomDistributions_D (Dict): Dictionary structure pointing to various number arrays based on the PDB ID used and the number of items drawn before averaging
        BranchToPValues_D (Dict): Dictionary structure pointing to the four P-Values for an ancestral, derived, PDB alignment triad
        """
        
        self.DataDIR = DataDIR
        
        #gets input/output directory and protein family name
        
        
        #gets all path information for the relevant input files
        self.ScopeXMLFile = reportPATH
        self.ScoringMatrixXMLFile = scoringPATH
        self.outPATH = outPATH
        
        #paths to more input files
        #self.HydroPATH = self.DataDIR+"misc/hydropathyindex"
        #self.MassPATH = self.DataDIR+"misc/sidechainmass"
        
        #makes a dictionary out of hydropathy index and mass input files
        #self.Hydro_D = {line.split()[0] : float(line.replace("\n","").split()[1]) for line in open(self.HydroPATH,"r").readlines()}
        #self.Mass_D = {line.split()[0] : float(line.replace("\n","").split()[1]) for line in open(self.MassPATH,"r").readlines()}
        
        #create the FastMLTree object and set branch lengths
        
        #parse out sequence information
        NodeToSequence_LD = self.getNodeToSequence_LD()
        self.NodeSequenceKey_L = NodeToSequence_LD[0]
        self.NodeToSequence_D = NodeToSequence_LD[1]
        
        #parse out report XML file and create ScopeAlgorithm instances
        self.BranchToAlgorithm_D = self.getBranchToAlgorithm_D()
        
        #set PDB content dictionary and PDBXML content dictionary
        
        PDB_L = []
        #for each branch key, check for new PDB ID keys
        for BranchKey in self.BranchToAlgorithm_D.keys():
            AccKeysSearch = re.compile("<PDBs>(.+?)</PDBs>").search(self.BranchToAlgorithm_D[BranchKey].alignmentSet)
            if AccKeysSearch:
                AccKeys_L = AccKeysSearch.group(1).split(";")
                #for each PDB ID key found
                for AccKey in AccKeys_L:        
                    #only executes if the PDB ID key has not already been added to the dictionary
                    if AccKey in set(PDB_L):
                        pass
                    else:
                        PDB_L.append(AccKey)
        
        PDBAndPDBXMLContents_Dicts = getAllPDBFileDicts(PDB_L)
        self.PDBContents_D = PDBAndPDBXMLContents_Dicts[0]
        self.PDBXMLContents_D = PDBAndPDBXMLContents_Dicts[1]
        [self.setPDBAndPDBXMLContentDictionaries(self.BranchToAlgorithm_D[Key]) for Key in self.BranchToAlgorithm_D.keys()]
        
        #parses out ScoringMatrixCoverage file
        self.ScoringMatrixCoverageKeys_D = self.getScoringMatrixCoverageKeys()
        self.ScoringMatrixPDBXMLMatchedKeys_D = self.getScoringMatrixPDBXMLMatchedKeys_D()
        
        #gets the indices for PDB IDs to be used in SAS and distance random distribution generation
        self.AccsToMutationCount_D = self.getNCoveredMutations_D()
        
        #self.AccsToDistanceCount_D = self.getNDistances_D()
        
        #get list of mass and hydropathy index change values for use in random distributions
        #BranchSegmentMutations_L = self.getAllBranchSegmentMutations()
        #self.MassChanges_L = BranchSegmentMutations_L[0]
        #self.HydroChanges_L = BranchSegmentMutations_L[1]
        
        self.RandomDistributions_D = self.getRandomDistributions_D() #create all random distributions
        #print self.RandomDistributions_D
        self.BranchToPValues_D = self.getAllBranchSegmentPValues() #get all PValues
        
        self.output() #output to PValue file
    
    "gets list and dictionary of node keys to sequences"
    def getNodeToSequence_LD(self):
        ret = {}
        retKey_L = []
        
        #parses all Seq headers in the xml report file, and makes a dictionary and list entry for each one
        AllSequences = re.findall("<Seq>.+?</Seq>" , open(self.ScopeXMLFile , "r").read())
        for Seq in AllSequences:
            SeqKey = re.compile("<H>(.+?)</H>").search(Seq).group(1)
            SeqObj = FASequence(SeqKey , re.compile("<S>(.+?)</S>").search(Seq).group(1))
            ret[SeqKey] = SeqObj
            retKey_L.append(SeqKey)
        return [retKey_L,ret]
    
    "gets dictionary of branch keys to ScopeAlgorithm objects representing those branches"
    def getBranchToAlgorithm_D(self):
        ret = {}
        key = "NA"
        val = open(self.ScopeXMLFile , "r").read()
        
        ret[key] = ScopeAlgorithmTrios(val)
        
        return ret
    
    
    

    
    "sets the ScopeAlgorithm instance's PDBContents_D and PDBXMLContents_D"
    def setPDBAndPDBXMLContentDictionaries(self,SA):
        SA.PDBContents_D = self.PDBContents_D
        SA.PDBXMLContents_D = self.PDBXMLContents_D
    
    "gets dictionary of covered residues in each PDB accession"
    def getScoringMatrixCoverageKeys(self):
        ret = {}
        #finds all PDBIDs and adds their key as a separate entry in the dictionary
        for CoverageKeySet in re.findall("<Coverage>.+?</Coverage>" , re.compile("(<Coverages>.+?</Coverages>)",re.DOTALL).search(open(self.ScoringMatrixXMLFile,"r").read()).group(1)):
            ret[re.compile("<ID>(.+?)</ID>",).search(CoverageKeySet).group(1)] = re.compile("<Keys>(.+?)</Keys>").search(CoverageKeySet).group(1).split(",") #get all covered keys
        return ret
    
    "gets dictionary of residues that are covered by the alignment and also covered by the parse PDB structure file"
    def getScoringMatrixPDBXMLMatchedKeys_D(self):
        Ret = {}
        
        for AccKey in self.PDBXMLContents_D.keys():
            
            Ret[AccKey] = []
            for PosKey in self.ScoringMatrixCoverageKeys_D[AccKey]:
                #checks if the coverage key is also in the PDBXML residue dictionary
                if PosKey in sets.Set(self.PDBXMLContents_D[AccKey]["XMLResidue_D"].keys()):
                    Ret[AccKey].append(PosKey)
        return Ret
    
    "gets dictionary of PDB IDs and the list of integers used to draw SAS random distributions"
    def getNCoveredMutations_D(self):
        
        AccsToMutationCount_D = {}
        for BranchKey in self.BranchToAlgorithm_D.keys():
            SA = self.BranchToAlgorithm_D[BranchKey]
            
            #for each ScopeAlgorithm instance, if mutations are present
            if SA.mutationsPresent:
                AccsToMutationCountForSingleSA_D = {}
                AccsToKey_D = {}
                
                #for each mutated site
                for MutationXMLKey in SA.mutationsXMLkey_L:
                    
                    #checks if that mutation has coverage, then gets the accession and position of that mutation
                    if SA.mutationScore_D[MutationXMLKey]["Coverage"]:
                        AllAccPos = SA.getAccessionPosition_L(SA.mutationsXML_D[MutationXMLKey])
                        for AccPos in AllAccPos:
                            Acc = AccPos[0]
                            Pos = AccPos[1]
                            if Acc in AccsToMutationCount_D.keys():
                                pass
                            else:
                                AccsToMutationCount_D[Acc] = []
                            
                            #adds the PDB ID to the count dictionary for a single ScopeAlgorithm
                            if Acc in AccsToMutationCountForSingleSA_D.keys():
                                pass
                            else:
                                AccsToMutationCountForSingleSA_D[Acc] = 0
                                AccsToKey_D[Acc] = []
                        
                            #adds one count to that PDB ID
                            AccsToMutationCountForSingleSA_D[Acc] += 1
                            AccsToKey_D[Acc].append(Pos)
                
                #sets the ScopeAlgorithm AccsToMutationCount to the SingleSA_D
                SA.AccsToMutationCount = AccsToMutationCountForSingleSA_D
                SA.AccsToKey_D = AccsToKey_D
                
                #adds the index to the overall dictionary if it is not already in there
                for Key in AccsToMutationCountForSingleSA_D.keys():
                    if AccsToMutationCountForSingleSA_D[Key] in sets.Set(AccsToMutationCount_D[Key]):
                        pass
                    else:
                        #only adds the index if it is greater than 2 mutations
                        if AccsToMutationCountForSingleSA_D[Key] >= 2:
                            AccsToMutationCount_D[Key].append(AccsToMutationCountForSingleSA_D[Key])
                            
        return AccsToMutationCount_D
    
    "gets dictionary of PDB IDs and the list of integers used to draw pairwise distance random distributions"
    """
    def getNDistances_D(self):
        AccsToDistanceCount_D = {}
        for BranchKey in self.BranchToAlgorithm_D.keys():
            SA = self.BranchToAlgorithm_D[BranchKey]
            
            #checks each ScopeAlgorithm instance if it has mutations present
            if SA.mutationsPresent:
                AccsToDistanceCountForSingleSA_D = {}
                AccsToDistanceKeyList_D = {}
                
                #checks each mutation for PDB coverage and adds it to the appropriate dictioanry (according to PDB ID)
                for MutationXMLKey in SA.mutationsXMLkey_L:
                    if SA.mutationScore_D[MutationXMLKey]["Coverage"]:
                        DistanceMarkers_L = SA.getAllDistanceMarkers(SA.mutationsXML_D[MutationXMLKey])
                        for DistanceMarker in DistanceMarkers_L:
                            DistancePositionAndAccession_L = SA.getDistancePositionKeysAndAccession_L(DistanceMarker)
                            if DistancePositionAndAccession_L[1] in AccsToDistanceKeyList_D.keys():
                                pass
                            else:
                                AccsToDistanceKeyList_D[DistancePositionAndAccession_L[1]] = []
                            
                            #fills graph list of highest set of positions that can be grouped together under the structural alignment bounds
                            PDBXMLMarkers = SA.getDistancePositionPDBXMLMarkers(DistancePositionAndAccession_L[0])
                            
                            if PDBXMLMarkers.split(";")[0] == PDBXMLMarkers.split(";")[1]:
                                pass
                            else:
                                if PDBXMLMarkers in sets.Set(AccsToDistanceKeyList_D[DistancePositionAndAccession_L[1]]):
                                    pass
                                else:
                                    if ";".join([PDBXMLMarkers.split(";")[1] , PDBXMLMarkers.split(";")[0]]) in sets.Set(AccsToDistanceKeyList_D[DistancePositionAndAccession_L[1]]):
                                        pass
                                    else:
                                        AccsToDistanceKeyList_D[DistancePositionAndAccession_L[1]].append(PDBXMLMarkers)
                
                #finds the most efficient (highest possible number) of graph connections for a PDB structure within the alignment bounds
                AccsToHighestCouplingKeyList_D = {}
                for AccKey in AccsToDistanceKeyList_D.keys():
                    PositionKeyToPartnersList_D = {}
                    
                    for PDBXMLMarker in AccsToDistanceKeyList_D[AccKey]:
                        KeyA = PDBXMLMarker.split(";")[0]
                        KeyB = PDBXMLMarker.split(";")[1]
                        
                        if KeyA in PositionKeyToPartnersList_D.keys():
                            pass
                        else:
                            PositionKeyToPartnersList_D[KeyA] = []
                        if KeyB in PositionKeyToPartnersList_D.keys():
                            pass
                        else:
                            PositionKeyToPartnersList_D[KeyB] = []
                        
                        PositionKeyToPartnersList_D[KeyA].append(KeyB)
                        PositionKeyToPartnersList_D[KeyB].append(KeyA)
                    
                    HighestPairing = 0
                    for Key in PositionKeyToPartnersList_D.keys():
                        if len(PositionKeyToPartnersList_D[Key]) > HighestPairing:
                            HighestPairing = len(PositionKeyToPartnersList_D[Key])
                    
                    HighestPairing_D = {Key:PositionKeyToPartnersList_D[Key] for Key in PositionKeyToPartnersList_D.keys() if len(PositionKeyToPartnersList_D[Key]) == HighestPairing}
                    if len(HighestPairing_D) > 2:
                        TempKey = HighestPairing_D.keys()[0]
                        HighestPairing_L = HighestPairing_D[TempKey]+[TempKey]
                        
                        if AccKey in AccsToDistanceCount_D.keys():
                            pass
                        else:
                            AccsToDistanceCount_D[AccKey] = []
                        if len(HighestPairing_L) in sets.Set(AccsToDistanceCount_D[AccKey]):
                            pass
                        else:
                            AccsToDistanceCount_D[AccKey].append(len(HighestPairing_L))
                        AccsToDistanceCountForSingleSA_D[AccKey] = len(HighestPairing_L)
                        AccsToHighestCouplingKeyList_D[AccKey] = HighestPairing_L
                    
                self.BranchToAlgorithm_D[BranchKey].AccsToDistanceCount = AccsToDistanceCountForSingleSA_D
                self.BranchToAlgorithm_D[BranchKey].AccsToDistanceKey_D = AccsToHighestCouplingKeyList_D
                
        return AccsToDistanceCount_D
    """
    
    "gets hydropathy index and mass change lists for all mutations that have occurred in the tree"
    """
    def getAllBranchSegmentMutations(self):
        MassRet = []
        HydroRet = []
        
        #gets all mutations from all branches and all states
        M_L = [[self.BranchToAlgorithm_D[B].getMutationType(self.BranchToAlgorithm_D[B].mutationsXML_D[MutationXMLKey]) for MutationXMLKey in self.BranchToAlgorithm_D[B].mutationsXMLkey_L if self.BranchToAlgorithm_D[B].mutationScore_D[MutationXMLKey]["Coverage"]] for B in self.BranchToAlgorithm_D.keys() if self.BranchToAlgorithm_D[B].mutationsPresent]
        Mut_L = []
        for M in M_L:
            for Mut in M:
                Mut_L.append(Mut)
        
        #makes lists of hydropathy index changes and mass changes from the list of state changes
        HydroRet = [self.getHydroDif(Mut[0],Mut[1]) for Mut in Mut_L]
        MassRet = [self.getMassDif(Mut[0],Mut[1]) for Mut in Mut_L]

        return [MassRet,HydroRet]
    """
        
    ####################################################################################################
    
    "gets squared difference in side chain mass between ancestral and derived sequence states"
    def getMassDif(self,StateA,StateB):
        return math.pow(self.Mass_D[StateA] - self.Mass_D[StateB] , 2)
    "gets squared difference in hydropathy index between ancestral and derived sequence states"
    def getHydroDif(self,StateA,StateB):
        return math.pow(self.Hydro_D[StateA] - self.Hydro_D[StateB] , 2)
    
    "gets list of all possible combinations of distances in a list of mutated sites from the same PDB structure"
    def getCombinatorialListOfPairwiseDistances(self,AccKey,AccsA_L):
        
        AccsB_L = AccsA_L[1:]
        count = 0
        Distances_L = []
        SA = self.BranchToAlgorithm_D[self.BranchToAlgorithm_D.keys()[0]]
        
        #combinations of pairwise distances
        for AccA in AccsA_L:
            for AccB in AccsB_L:
                #gets both points
                APDBXMLLine = SA.getPDBXMLLine(AccKey,AccA)
                BPDBXMLLine = SA.getPDBXMLLine(AccKey,AccB)
                
                if APDBXMLLine and BPDBXMLLine:    
                    APoint = SA.getAlphaCarbonPoint(APDBXMLLine)
                    BPoint = SA.getAlphaCarbonPoint(BPDBXMLLine)
                
                    #as long as they are not null, calculate the magnitude of distance between them and add it to a list
                    if APoint and BPoint:
                        Distances_L.append(SA.getDistanceMagnitude(APoint , BPoint))
            AccsB_L = AccsB_L[1:]
        
        return Distances_L
    
    "general method for writing/retrieving a random distribution array"
    def getAnyAverageRandomDist(self,AveragedNumbers_L):
        FinalNumbers_L = [AveragedNumber for AveragedNumber in AveragedNumbers_L if math.isnan(AveragedNumber) == False]
        return array(FinalNumbers_L)
    "gets the average of numbers in a list"
    def getAveragedData(self,Numbers_L):
        return numpy.mean([Number for Number in Numbers_L if Number != None])
    "gets a random sample of integers to be used as random indices to draw numbers for the random distributions"
    def getRandomSampleOfIntegers(self,MaxLength,Index):
        return random.sample(range(MaxLength),Index)
    
    "get SAS random distribution for a single PDB ID and index"
    def getRelativeSASRandDistForIndex(self,AccKey,Index):
        Ret = None
        #print AccKey
        #print Index
        #gets 10000 SAS averages of randomly selected indices on the protein structure within the alignment bounds
        try:
        
            Ret = self.getAnyAverageRandomDist(\
                    [self.getAveragedData(\
                        [self.BranchToAlgorithm_D[self.BranchToAlgorithm_D.keys()[0]].getRelativeGlobalSAS(self.PDBXMLContents_D[AccKey]["XMLResidue_D"][self.ScoringMatrixPDBXMLMatchedKeys_D[AccKey][xIndex]]) \
                            for xIndex in self.getRandomSampleOfIntegers(len(self.ScoringMatrixPDBXMLMatchedKeys_D[AccKey]),Index)]) for i in range(0,10000)])
        except Exception as e:
            print e
            
        return Ret
    
    "get distance random distribution for a single PDB ID and index"
    def getRelativeDistanceRandDistForIndex(self,AccKey,Index):
        Ret = None
        
        #gets 10000 distance averages of randomly selected indices on the protein structure within the alignment bounds
        try:
            Ret = self.getAnyAverageRandomDist(\
                    [self.getAveragedData(\
                        self.getCombinatorialListOfPairwiseDistances(AccKey,[self.ScoringMatrixPDBXMLMatchedKeys_D[AccKey][xIndex]\
                            for xIndex in self.getRandomSampleOfIntegers(len(self.ScoringMatrixPDBXMLMatchedKeys_D[AccKey]),Index)])) for i in range(0,10000)])
        except Exception as e:
            print e
        return Ret
    
    "get hydropathy index change random distribution for a single PDB ID and index"
    def getRelativeHydroRandDistForIndex(self,AccKey,Index):
        
        #gets 10000 hydropathy index change averages of randomly selected indices on the protein structure within the alignment bounds
        try:
            Ret = self.getAnyAverageRandomDist(
            [self.getAveragedData(
                [self.HydroChanges_L[xIndex] for xIndex in self.getRandomSampleOfIntegers(len(self.HydroChanges_L),Index)]) for a in range(0,10000)])
        except Exception as e:
            print e
        return Ret
    
    "get mass change random distribution for a single PDB ID and index"
    def getRelativeMassRandDistForIndex(self,AccKey,Index):
        Ret = None
        
        #gets 10000 mass change averages of randomly selected indices on the protein structure within the alignment bounds
        try:
            
            Ret = self.getAnyAverageRandomDist(
            [self.getAveragedData(
                [self.MassChanges_L[xIndex] for xIndex in self.getRandomSampleOfIntegers(len(self.MassChanges_L),Index)]) for a in range(0,10000)])
        except Exception as e:
            print e
        return Ret
    
    "gets SAS random distributions for all indices for one PDB ID"
    def getRelativeSASRandDistAllIndices(self,AccKey,Indices):
        return {str(Index) : self.getRelativeSASRandDistForIndex(AccKey,Index) for Index in Indices}
    "gets distance random distributions for all indices for one PDB ID"
    def getDistanceRandDistAllIndices(self,AccKey,Indices):
        return {str(Index) : self.getRelativeDistanceRandDistForIndex(AccKey,Index) for Index in Indices}
    "gets hydropathy index change random distributions for all indices for one PDB ID"
    def getRelativeHydroRandDistAllIndices(self,AccKey,Indices):
        return {str(Index) : self.getRelativeHydroRandDistForIndex(AccKey,Index) for Index in Indices}
    "gets mass change random distributions for all indices for one PDB ID"
    def getRelativeMassRandDistAllIndices(self,AccKey,Indices):
        return {str(Index) : self.getRelativeMassRandDistForIndex(AccKey,Index) for Index in Indices}
    
    "gets all random distributions for each criterion and each PDB ID key and each index (ie. highest method for random distributions)"
    def getRandomDistributions_D(self):
        return {"SAS":{AccKey : self.getRelativeSASRandDistAllIndices(AccKey,self.AccsToMutationCount_D[AccKey]) for AccKey in self.AccsToMutationCount_D.keys()},\
                "Dist": {AccKey : self.getDistanceRandDistAllIndices(AccKey,self.AccsToMutationCount_D[AccKey]) for AccKey in self.AccsToMutationCount_D.keys()}}
                #"Hydro": {AccKey : self.getRelativeHydroRandDistAllIndices(AccKey,self.AccsToMutationCount_D[AccKey]) for AccKey in self.AccsToMutationCount_D.keys()},\
                #"Mass" : {AccKey : self.getRelativeMassRandDistAllIndices(AccKey,self.AccsToMutationCount_D[AccKey]) for AccKey in self.AccsToMutationCount_D.keys()}}
    
    ####################################################################################################
    
    "general method for retrieving a P-value"
    def getGeneralPValue(self,Criterion,Acc,Index,Point):
        Ret = None
        
        #does not execute if there are only 0.0's in the Random distribution (faulty distribution)
        if len(set(list(self.RandomDistributions_D[Criterion][Acc][str(Index)]))) == 1:
            pass
        #gets CDF function (percentile) using the observed point as the input and distribution as the background
        else:
            Ret = percentileofscore(self.RandomDistributions_D[Criterion][Acc][str(Index)],Point)
        
        return Ret
    
    "get SAS p-value for a particular ancestral, derived, triple alignment triad"
    def getSASPValue(self,SA,AccKey):
        PVal = None
        Avg = None
        Ret = None
        
        if str(SA.AccsToMutationCount[AccKey]) in self.RandomDistributions_D["SAS"][AccKey].keys():
            try:
                RSAS_L = []
                
                for AccPos in SA.getAllPositionKeysAccordingToAccession(AccKey):
                    pdbxmlLine = SA.getPDBXMLLine(AccPos[0],AccPos[1])
                    if pdbxmlLine:
                        SASToAdd = SA.getRelativeGlobalSAS(pdbxmlLine)
                        RSAS_L.append(SASToAdd)
                
                Avg = self.getAveragedData(RSAS_L)
                
                if numpy.isnan(Avg):
                    pass
                else:
                    PVal = self.getGeneralPValue("SAS",AccKey,SA.AccsToMutationCount[AccKey],Avg)  / 100.0
                    Ret = [Avg,PVal]
                
                
            except Exception as e:
                pass
                
        return Ret
    
    "get distance p-value for a particular ancestral, derived, triple alignment triad"
    def getDistPValue(self,SA,AccKey):
        PVal = None
        Avg = None
        Ret = None
            
        if AccKey in SA.AccsToMutationCount.keys():
            if str(SA.AccsToMutationCount[AccKey]) in self.RandomDistributions_D["Dist"][AccKey].keys():
                try:
                    Pairwise_L = self.getCombinatorialListOfPairwiseDistances(AccKey,SA.AccsToKey_D[AccKey])
                    Avg = self.getAveragedData(Pairwise_L)
                    if numpy.isnan(Avg):
                        pass
                    else:
                        PVal = self.getGeneralPValue("Dist",AccKey,SA.AccsToMutationCount[AccKey],Avg) / 100.0
                        Ret = [Avg,PVal]
                    #Ret = self.getGeneralPValue("Dist",AccKey,SA.AccsToMutationCount[AccKey],self.getAveragedData(self.getCombinatorialListOfPairwiseDistances(AccKey,SA.AccsToKey_D[AccKey]))) / 100.0
                    
                    
                except Exception as e:
                    print e
        return Ret
    
    "gets PValues for all four criteria for one ancestral,derived, PDB alignment triad"
    def getAllPValuesForAccession(self,SA,AccKey):
        return {"SAS":self.getSASPValue(SA,AccKey),\
                "Dist":self.getDistPValue(SA,AccKey)}
                
    "gets all PValues for one ancestral, derived alignment pair"
    def getAllPValuesForBranchSegment(self,SA):
        return {Acc:self.getAllPValuesForAccession(SA,Acc) for Acc in self.AccsToMutationCount_D.keys() if Acc in SA.AccsToMutationCount.keys()}
    "gets all PValues for all ancestral, derived alignment pairs"
    def getAllBranchSegmentPValues(self):
        return {BranchKey : self.getAllPValuesForBranchSegment(self.BranchToAlgorithm_D[BranchKey]) for BranchKey in self.BranchToAlgorithm_D.keys() if self.BranchToAlgorithm_D[BranchKey].mutationsPresent}
    
    
    ####################################################################################################
    
    "writes PValue information to the appropriate file"
    def output(self):
        AllOutput_L = ["Branch               PDB   #M   Msas  Psas  Mdis  Pdis"] #header line
        
        #for each ancestral, derived, PDB alignment triad
        for BranchKey in self.BranchToPValues_D.keys():
            
            for AccKey in self.BranchToPValues_D[BranchKey].keys():
                
                #if there is a PValue to this triad, then it will format an appropriate output string
                OutputString = None
                if self.BranchToPValues_D[BranchKey][AccKey]["SAS"] and self.BranchToPValues_D[BranchKey][AccKey]["Dist"]:
                    if self.BranchToPValues_D[BranchKey][AccKey]["Dist"][0] != 0.0:
                        
                        OutputString = "%s%s %s%s %s%s %s%s %s%s %s%s %s" %   (BranchKey, " "*(20-len(BranchKey)),\
                                                                               AccKey," "*(5-len(AccKey)),\
                                                                               str(len(self.BranchToAlgorithm_D[BranchKey].getAllMutationXMLAccordingToAccession(AccKey))), " "*(4-len(str(len(self.BranchToAlgorithm_D[BranchKey].getAllMutationXMLAccordingToAccession(AccKey))))),\
                                                                               str(round(self.BranchToPValues_D[BranchKey][AccKey]["SAS"][0] , 3)), " "*(5-len(str(round(self.BranchToPValues_D[BranchKey][AccKey]["SAS"][0] , 3)))),\
                                                                               str(round(self.BranchToPValues_D[BranchKey][AccKey]["SAS"][1] , 3)), " "*(5-len(str(round(self.BranchToPValues_D[BranchKey][AccKey]["SAS"][1] , 3)))),\
                                                                               str(round(self.BranchToPValues_D[BranchKey][AccKey]["Dist"][0] , 3)), " "*(5-len(str(round(self.BranchToPValues_D[BranchKey][AccKey]["Dist"][0] , 3)))),\
                                                                               str(round(self.BranchToPValues_D[BranchKey][AccKey]["Dist"][1] , 3)))
                        AllOutput_L.append(OutputString)
                                                                   
        #writes PValues to output file and displays end prompt.
         
        with open(self.outPATH,"w") as w:
            w.write("\n".join(AllOutput_L))
            
        print "Done.\nAll P-Values written to %s" % (self.outPATH)
    