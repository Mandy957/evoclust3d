####################################################################################################
#                                                                                                  #
#           PROJECT Protein Adaptation                                                             #
#           CLASS MercatorParser                                                                   #
#           PROGRAMMER Jeremy Adams                                                                #
#           STARTED 16-12-14                                                                       #
#           LASTMOD 16-12-14                                                                       #
#                                                                                                  #
#           DESCRIPTION creates scope input sets from the orthologous genes found by mercator      #
#                                                                                                  #
####################################################################################################

import sys
import os
import re
import random
from cogent import LoadTree
from staticmethods import *

class MercatorParser:
    
    """
    Class attributes:
    DIR (String): Path to the Genome general directory
    GFF (String): Path to GFF annotation files
    Fasta (String): Path to the fasta file directory containing the genomes
    ProjectName (String): String reflecting the specific desired merctor project (and directory in the Mercator Projects folder)
    ProjectDIR (String): Path to the specific mercator project of interest
    """
    
    "CONSTRUCTOR"
    def __init__(self, DataDIR , ProjectName):
        
        #sets all directory variables
        self.DataDIR = DataDIR
        self.DIR = self.DataDIR+"Genomes/"
        self.GFF = self.DIR+"GFF/"
        self.Fasta = self.DIR+"ModdedComplete/"
        
        self.ProjectName = ProjectName
        if self.ProjectName.endswith("/"):
            self.ProjectName = self.ProjectName[:-1]
        
        self.ProjectDIR = self.DIR+"MercatorProjects/"+self.ProjectName+"/"
        self.MapPATH = self.ProjectDIR+"map"
        
        #sets species names and parses phylogenetic tree
        self.SpeciesNames_L = open("%sgenomes" % (self.ProjectDIR) , "r").read().replace("\n","").split()
        self.TreePath = self.ProjectDIR+"SpeciesTree.nwk"
        self.OrigTree = LoadTree(self.TreePath)
        self.SpeciesToProteins_D = {spec : {} for spec in self.SpeciesNames_L}
        
        #parses the GFF file for each species
        for SpeciesName in self.SpeciesNames_L:
            CDSIDToCDSCoords_D = self.parseGFF("%s%s.gff" % (self.GFF,SpeciesName))
            Fasta = readFasta("%s%s.fa" % (self.Fasta,SpeciesName))
            FastaKey_L = Fasta[0]
            Fasta_D = Fasta[1]
            
            self.SpeciesToProteins_D[SpeciesName] = self.JoinCDSAndTranslate(CDSIDToCDSCoords_D,Fasta_D)
        
        #parse the map and then create input sets for each line in the map
        self.OutTreeDIR = self.ProjectDIR+"Trees/"
        self.OutAlignmentDIR = self.ProjectDIR+"Alignments/"
        if os.path.exists(self.OutTreeDIR):
            pass
        else:
            os.system("mkdir %s" % (self.OutTreeDIR))
        
        if os.path.exists(self.OutAlignmentDIR):
            pass
        else:
            os.system("mkdir %s" % (self.OutAlignmentDIR))
        
        parsedMap = self.parseMap()
        self.parsedMapKey_L = parsedMap[0]
        self.parsedMap_D = parsedMap[1]
        self.createInputSets()
        
    "parses a GFF file for a single species, gets location and sequences of proteins"
    def parseGFF(self, PATH):
        
        allCDSLines_L = []
        for line in open(PATH,"r").readlines():
            if line[0] != "#":
                if line.split()[2] == "CDS":
                    allCDSLines_L.append(line.replace("\n",""))
        
        CDSIDToCDSCoords_D = {}
        
        previouscdsID = re.compile("ID=(cds.+?);").search(allCDSLines_L[0]).group(1)
        CDSIDToCDSCoords_D[previouscdsID] = []
        
        for line in allCDSLines_L:
            
            ls = line.split()
            
            cdsID = re.compile("ID=(cds.+?);").search(line).group(1)
                
            if cdsID == previouscdsID:
                pass
            else:
                CDSIDToCDSCoords_D[cdsID] = []
                previouscdsID = cdsID
            
            CDSIDToCDSCoords_D[cdsID].append([ls[0],ls[3],ls[4],ls[6],ls[7]])
        
        s = 0
        
        for key in CDSIDToCDSCoords_D.keys():
            s += len(CDSIDToCDSCoords_D[key])
        
        return CDSIDToCDSCoords_D
    
    "joins a list of CDS sequences and translates them into a single protein sequence"
    def JoinCDSAndTranslate(self, CDSIDToCDSCoords_D , Fasta_D):
        Ret = {}
        goodcount = 0
        badcount = 0
        for Key in CDSIDToCDSCoords_D.keys():
            
            Sign = ""
            Contig = CDSIDToCDSCoords_D[Key][0][0]
            DNASeq = []
            NoN = True
            
            for item in CDSIDToCDSCoords_D[Key]:
                
                miniSeq = Fasta_D[item[0]][int(item[1])-1:int(item[2])]
                if re.compile("N").search(miniSeq):
                    NoN = False
                
                if NoN:
                    if item[3] == "+":
                        DNASeq.append(miniSeq)
                        Sign = "+"
                    else:
                        DNASeq.append(reverseComplement(miniSeq))
                        Sign = "-"
            
            
            if NoN:
                DNASeq = "".join(DNASeq)
                ProteinSeq = translate(convertTomRNA(DNASeq))
            
                if ProteinSeq[0] == "M" and ProteinSeq[-1] == "*":
                    if ProteinSeq.count("*") == 1:
                        ProteinSeq = ProteinSeq[:-1]
                        
                        SortedStartCoords_L = sorted([int(x[1]) for x in CDSIDToCDSCoords_D[Key]] , cmp=numeric_compare)
                        SortedEndCoords_L = sorted([int(x[2]) for x in CDSIDToCDSCoords_D[Key]] , cmp=numeric_compare)
                        
                        start = str(SortedStartCoords_L[0])
                        end = str(SortedEndCoords_L[-1])
                        
                        if Contig in Ret.keys():
                            pass
                        else:
                            Ret[Contig] = {"+":{},"-":{}}
                        Ret[Contig][Sign]["%s-%s" % (start,end)] = ProteinSeq
        return Ret
    
    "parses all lines of the ortholog alignment map from mercator"
    def parseMap(self):
        Key_L = []
        D = {}
        
        mapLines_L = [line.replace("\n","") for line in open(self.MapPATH,"r").readlines()]
        
        for line in mapLines_L:
            
            ls = line.split()
            Key_L.append(ls[0])
            D[ls[0]] = {"NAs" : [] , "Coords" : {Key : [] for Key in self.SpeciesNames_L}}
            
            statsForOneSpecies = []
            i = 0
            specCount = 0
            
            for elem in ls[1:]:
                if i == 0:
                    if elem == "NA":
                        D[ls[0]]["NAs"].append(self.SpeciesNames_L[specCount])
                
                if i < 4:
                    statsForOneSpecies.append(elem)
                    i += 1
                else:
                    D[ls[0]]["Coords"][self.SpeciesNames_L[specCount]] = statsForOneSpecies
                    specCount += 1
                    i = 0
                    statsForOneSpecies = []
                    statsForOneSpecies.append(elem)
                    i += 1
                
            
            D[ls[0]]["Coords"][self.SpeciesNames_L[specCount]] = statsForOneSpecies
            
        return [Key_L, D]
    
    "creates input sets for the scope algorithm for all lines of the map"   
    def createInputSets(self):
        
        for mapKey in self.parsedMapKey_L[107:]:
            SpeciesToProteinsForThisMapLine_D = {Key : {} for Key in self.SpeciesNames_L}
            
            for speciesKey in self.SpeciesNames_L:
                if speciesKey in set(self.parsedMap_D[mapKey]["NAs"]):
                    pass
                else:
                    
                    mapSpeciesSegment = self.parsedMap_D[mapKey]["Coords"][speciesKey]
                    contigKey = mapSpeciesSegment[0]
                    signKey = mapSpeciesSegment[3]
                    startBoundarystr = mapSpeciesSegment[1]
                    finishBoundarystr = mapSpeciesSegment[2]
                    
                    if contigKey == "NA" or signKey == "NA" or startBoundarystr == "NA" or finishBoundarystr == "NA":
                        pass
                    else:
                        
                        startBoundary = int(startBoundarystr)
                        finishBoundary = int(finishBoundarystr)
                        
                        if contigKey in self.SpeciesToProteins_D[speciesKey].keys():
                            for StartFinishKey in self.SpeciesToProteins_D[speciesKey][contigKey][signKey].keys():
                                
                                s = StartFinishKey.split("-")
                                StartState = int(s[0])
                                FinishState = int(s[1])
                                
                                if StartState > startBoundary and FinishState < finishBoundary:
                                    distanceFromStart = str(StartState - startBoundary)
                                    length = str(FinishState - StartState + 1)
                                    distanceFromEnd = str(finishBoundary - FinishState)
                                    
                                    if signKey == "-":
                                        distanceFromEnd = str(StartState - startBoundary)
                                        distanceFromStart = str(finishBoundary - FinishState)
                                    
                                    SpeciesToProteinsForThisMapLine_D[speciesKey]["%s-%s-%s" % (distanceFromStart,length,distanceFromEnd)] = self.SpeciesToProteins_D[speciesKey][contigKey][signKey][StartFinishKey]
            
            self.createInputSetsForOneMapLine(SpeciesToProteinsForThisMapLine_D , mapKey)
    
    "extends the above method by creating all the input sets belonging in one line of the map"
    def createInputSetsForOneMapLine(self, SpeciesToProteinsMapLine_D , mapKey):
        
        subCount = 1
        
        Used_D = {SpeciesKey : {PosKey : False for PosKey in SpeciesToProteinsMapLine_D[SpeciesKey].keys()} for SpeciesKey in SpeciesToProteinsMapLine_D.keys()}
        
        while(self.checkUsed_D(Used_D)):
        
            LowestNonZero = self.getLowestNonZero(Used_D)
            LowestNonZeroSpecies = LowestNonZero[0]
            
            anchorPos = ""
            B= True
            
            for PosKey in SpeciesToProteinsMapLine_D[LowestNonZeroSpecies].keys():
                if B:
                    if Used_D[LowestNonZeroSpecies][PosKey] == False:
                        anchorPos = PosKey
                        B = False
            
            anchorPosSplit = [int(Pos) for Pos in anchorPos.split("-")]
            
            Ortholog_D = {LowestNonZeroSpecies : anchorPos}
            
            for SpeciesKey in SpeciesToProteinsMapLine_D.keys():
                if SpeciesKey == LowestNonZeroSpecies:
                    pass
                else:
                    
                    PositionsRemaining_L = self.checkIfPositionsRemaining(Used_D[SpeciesKey])
                    RemainingCount = PositionsRemaining_L[0]
                    Positions_L = PositionsRemaining_L[1]
                    
                    if RemainingCount != 0:
                        BestTestPosition = ""
                        BestTestPositionScore = float("inf")
                        
                        for TestPosition in Positions_L:
                            TestPosSplit = [int(Pos) for Pos in TestPosition.split("-")]
                            
                            TestScoreA = abs(TestPosSplit[0] - anchorPosSplit[0])
                            TestScoreB = abs(TestPosSplit[2] - anchorPosSplit[2])
                            TestScore = TestScoreA + TestScoreB
                            
                            if TestScore < BestTestPositionScore:
                                BestTestPositionScore = TestScore
                                BestTestPosition = TestPosition
                            
                        Ortholog_D[SpeciesKey] = BestTestPosition
            
            
            for SPKey in Ortholog_D.keys():
                Used_D[SPKey][Ortholog_D[SPKey]] = True
            
            clustalInput_L = []
            
            for SPKey in Ortholog_D.keys():
                clustalInput_L.append(">"+SPKey)
                clustalInput_L.append(SpeciesToProteinsMapLine_D[SPKey][Ortholog_D[SPKey]])
                
            clustalInput = "\n".join(clustalInput_L)
            clustalOutput = executeClustalW(clustalInput)
            alignedSequences_D = parseClustalW(clustalOutput)
            
            self.writeOneOrthologSet(alignedSequences_D , mapKey , subCount)
            subCount += 1
            
    "gets the species with the least amount of proteins in the set (but with more than 0)"
    def getLowestNonZero(self,Used_D):
        Lowest = float("inf")
        LowestKey = ""
        for Key in Used_D.keys():
            Count = 0
            for PosKey in Used_D[Key].keys():
                if Used_D[Key][PosKey] == False:
                    Count += 1
            
            if Count != 0:
                if Count < Lowest:
                    Lowest = Count
                    LowestKey = Key
        
        return [LowestKey,Lowest]
    
    "checks the Used_D to see if the species has any proteins left to be aligned with others"
    def checkIfPositionsRemaining(self, Used_DSubset):
        Count = 0
        PosKeys = []
        for key in Used_DSubset.keys():
            if Used_DSubset[key] == False:
                Count += 1
                PosKeys.append(key)
        
        return [Count,PosKeys]
                        
    "checks the dictionary to see if there are any more ortholog sets for the given map line"
    def checkUsed_D(self, Used_D):
        
        numFullyUsed = 0
        
        for Key in Used_D.keys():
            L = [Used_D[Key][pos] for pos in Used_D[Key].keys()]
            
            if len(set(L)) == 1:
                if list(set(L))[0] == True:
                    numFullyUsed += 1
            
            elif len(set(L)) == 0:
                numFullyUsed += 1
        
        numLeft = len(Used_D.keys()) - numFullyUsed
        
        Ret = True
        if numLeft < 3:
            Ret = False
        
        return Ret
        
    def calculatePercentIdentity(self,seq_D):
        seqs_L = seq_D.values()
        
        m = 0
        
        for i in range(0,len(seqs_L[0])):
            states_L = [seq[i] for seq in seqs_L]
            if len(set(states_L)) == 1:
                if states_L[0] != "-":
                    m += 1
        
        return float(float(m) / float(len(seqs_L[0])))
    
    "gets a set of orthologous proteins, then writes the appropriate tree and alignment fasta files for scope"  
    def writeOneOrthologSet(self, alignment_D,mapLine,subCount):
        
        outputAlignment = []
        for speciesKey in self.SpeciesNames_L:
            if speciesKey in alignment_D.keys():
                outputAlignment.append(">%s\n%s" % (speciesKey,alignment_D[speciesKey]))
        
        outputAlignment = "\n".join(outputAlignment)
        
        LongToShort_D = {}
        ShortToLong_D = {}
        
        for speciesKey in self.SpeciesNames_L:
            if len(speciesKey) > 10:
                LongToShort_D[speciesKey] = speciesKey[0:10]
                ShortToLong_D[speciesKey[0:10]] = speciesKey
        
        copyTree = self.OrigTree
        
        for speciesKey in self.SpeciesNames_L:
            if speciesKey in alignment_D.keys():
                pass
            else:
                copyTree.remove(speciesKey)
                copyTree.prune()
        
        newTreeAsNewick = copyTree.getNewick()
        
        for key in LongToShort_D.keys():
            newTreeAsNewick = newTreeAsNewick.replace(key, LongToShort_D[key])
        
        
        PhyMLResults = usePhyMLForBranchLengths(outputAlignment,newTreeAsNewick)
        Success = PhyMLResults[0]
        returnNewick = PhyMLResults[1]
        
        if Success:
            
            for key in ShortToLong_D.keys():
                returnNewick = returnNewick.replace(key,ShortToLong_D[key]).replace("\n","")
                
            open("%s%s-%s.fa" % (self.OutAlignmentDIR,mapLine,str(subCount)) , "w").write(outputAlignment)
            open("%s%s-%s.nwk" %(self.OutTreeDIR,mapLine,str(subCount)) , "w").write(returnNewick)
        
        print "%s-%s" % (mapLine,str(subCount))
        
        
        