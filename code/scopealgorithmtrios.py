####################################################################################################
#                                                                                                  #
#           PROJECT Protein Adaptation                                                             #
#           CLASS ScopeAlgorithm                                                                   #
#           PROGRAMMER Jeremy Adams                                                                #
#           STARTED 25-07-13                                                                       #
#           LASTMOD 06-11-14                                                                       #
#                                                                                                  #
#           DESCRIPTION class representation of a single branch segment and associated observed    #
#                       mutation properties used for statistical analysis                          #
#                                                                                                  #
####################################################################################################

import os
import math
import re
import sys
import random
from numpy import array

class ScopeAlgorithmTrios:                    
    
    "dictionary of accessible surface area when the amino acid exists independent of a protein"
    def getIdealizedSASA_D(self):
        return {'A': 102.77, 'R': 279.20, 'N': 174.75, 'D': 172.69,\
                'C': 143.36, 'Q': 220.26, 'E': 214.78, 'G':  47.92,\
                'H': 221.38, 'I': 241.44, 'L': 205.91, 'K': 246.36,\
                'M': 213.64, 'F': 224.78, 'P': 170.15, 'S': 126.89,\
                'T': 182.33, 'W': 243.47, 'Y': 261.73, 'V': 204.59,\
                'X': 1.0}    
    
    "CONSTRUCTOR"
    def __init__(self, alignmentSet):
        
        """
        Class attributes:
        alignmentSet (String): Contains XML format information of what the mutations in this branch segment aligned to
        BranchName (String): Key corresponding to the branch name of the branch segment as it appears in the modified Newick file
        pdbsPresent (Bool): False if the reference sequence did not align to anything, True if it did
        pdbIdsHaveFiles (Bool): False if there are no files in the PDB directory for the hit structures, True if there are
        mutationsPresent (Bool): False if there are no mutations, True if there are
        pdbAccessionKeysSearch (Search Obj): Contains regex search information on what PDB ID keys are contained
        nMutations (Int): number of mutations for this particular branch length
        mutationsXML_L (List): list of all the mutation information lines for the branch segment
        pdbAccessionKeys_L (List): list of the PDB IDs that the branch segment aligned to
        
        mutationsXML_D (Dict): Key is the mutation position, value is the full mutation XML line
        mutationsXMLkey_L (List): List of all mutation position keys
        mutationType_D (Dict): Key is the mutation position, value is the ancestral and derived states for that residue position
        mutationQuality_D (Dict): Key is the mutation position, value is whether the reconstruction probabilities are satisfactory or unsatisfactory
        mutationScore_D (Dict): Key is mutation position, value is whether or not there is coverage for this residue
        """
        
        #starting conditions
        self.alignmentSet = alignmentSet
        
        self.pdbsPresent = False
        self.pdbIdsHaveFiles = True
        self.mutationsPresent = False
        self.pdbAccessionKeysSearch = re.compile("<PDBs>(.+?)</PDBs>").search(alignmentSet)
        
        self.nMutations = 0
        
        #finds all mutation lines in the XML file
        self.mutationsXML_L = re.findall("<M>.+?</R>" , alignmentSet)
        
        #if the branch segment has valid PDb IDs, then it makes a list out of these IDs
        if self.pdbAccessionKeysSearch:
            self.pdbAccessionKeys_L = self.pdbAccessionKeysSearch.group(1).split(";")
        
        #checks to see if there are mutations in this branch segment
        if len(self.mutationsXML_L) == 0:
            pass
        else:
            self.mutationsPresent = True
            
            #gets all the details from the mutation XML and places it in appropriate dictionaries/lists
            mutationsXML_H = self.get_mutationsXML_H(self.mutationsXML_L)
            self.mutationsXML_D = mutationsXML_H['mutationXML_D']
            self.mutationsXMLkey_L = mutationsXML_H['mutationXMLkey_L']
            self.mutationType_D = mutationsXML_H['mutationType_D']
            self.mutationQuality_D = mutationsXML_H['mutationQuality_D']
            self.mutationScore_D = {}
            
            self.nMutations = len(self.mutationsXMLkey_L)
            
            #checks to see if each mutation has a full mutation XML line and therefore is covered or not
            for mutationsXMLkey in self.mutationsXMLkey_L:
                self.mutationScore_D[mutationsXMLkey] = {}
                self.mutationScore_D[mutationsXMLkey]['Coverage'] = True
                
                if re.compile("<R>NOCOVERAGE</R>").search(self.mutationsXML_D[mutationsXMLkey]):
                    self.mutationScore_D[mutationsXMLkey]['Coverage'] = False

    "gets all dictionaries pertaining to the mutations"                 
    def get_mutationsXML_H(self, mutationsXML_L):
        mutationXML_D = {}
        mutationXMLkey_L = []
        mutationType_D = {}
        mutationQuality_D = {}
        
        #for each mutation, get the original line, key, type, and quality score
        for mutation in mutationsXML_L:
            key = re.compile("<M>(.+?)</M>").search(mutation).group(1).split("|")[0]
            substitution = re.compile("<M>(.+?)</M>").search(mutation).group(1).split("|")[1]
            quality = re.compile("<Q>(.+?)</Q>").search(mutation).group(1)
            
            mutationXML_D[key] = mutation
            mutationXMLkey_L.append(key)
            mutationType_D[key] = substitution
            mutationQuality_D[key] = quality
        
        return {'mutationXML_D' : mutationXML_D,\
                'mutationXMLkey_L' : mutationXMLkey_L,\
                'mutationType_D' : mutationType_D,\
                'mutationQuality_D' : mutationQuality_D}
    
    "gets mutation type of a mutationXML line"
    def getMutationType(self,mutationXML):
        
        return re.compile("<M>(.+?)</M>").search(mutationXML).group(1).split("|")[1]
    
    "gets PDB ID, chain, and position of a mutationXML line"
    def getAccessionPosition_L(self,mutationXML):
        Ret = None
        base = re.compile("<R>(.+?)</R>").search(mutationXML).group(1)
        
        if base == "NOCOVERAGE":
            pass
        else:
            allpdbaccs = base.split(",")
            Ret = [[pdb.split("|")[0].lower() , pdb.split("|")[1]] for pdb in allpdbaccs]
            
        return Ret
    
    "gets the appropriate PDBXML line using accession and position information"
    def getPDBXMLLine(self,accession,position):
        Ret = None
        if position in self.PDBXMLContents_D[accession]['XMLResidue_D'].keys():
            Ret = self.PDBXMLContents_D[accession]['XMLResidue_D'][position]
        else:
            pass
            #print accession+" "+position + " not found"
        
        return Ret
    
    "gets the amino acid of the PDB residue described in the PDBXML line"
    def getPDBXMLResidueType(self,pdbXMLLine ):
        return re.compile("<t>(.+?)</t>").search(pdbXMLLine).group(1)
    
    "gets both global and local SAS values from a PDB XML line"
    def getPDBXMLSAS(self, pdbXMLLine):
        return re.compile("<s>(.+?)</s>").search(pdbXMLLine).group(1).split(";")
    "get global SAS from a PDB XML line"
    def getGlobalSAS(self, pdbXMLLine):
        return float(self.getPDBXMLSAS(pdbXMLLine)[0])
    "get local SAS from a PDB XML line"
    def getLocalSAS(self, pdbXMLLine):
        return float(self.getPDBXMLSAS(pdbXMLLine)[1])
    
    "gets absolute global SAS for a pdbXML residue and divides it by the ideal to get the relative global sas"
    def getRelativeGlobalSAS(self,pdbXMLLine):
        return self.getGlobalSAS(pdbXMLLine) / self.getIdealizedSASA_D()[self.getPDBXMLResidueType(pdbXMLLine)]
    "gets absolute local SAS for a pdbXML residue and divides it by the ideal to get the relative local sas"
    def getRelativeLocalSAS(self,pdbXMLLine):
        return self.getLocalSAS(pdbXMLLine) / self.getIdealizedSASA_D()[self.getPDBXMLResidueType(pdbXMLLine)]
    
    
    "gets all keys for pairwise distances between residues"
    def getAllDistanceMarkers(self,mutationXMLLine):
        return [re.compile("<v>(.+?)</v>").search(SubsetString).group(1) for SubsetString in re.findall("<v>.+?</v>",mutationXMLLine)]
    "gets a split list of 2 indices where one index is the Accession, and the other index is the position"
    def getDistancePositionKeysAndAccession_L(self,DistanceMarker):
        return DistanceMarker.split(";")
    "gets distance position PDBXML Markers"
    def getDistancePositionPDBXMLMarkers(self,DistancePositionKeys):
        return ";".join([DistancePositionKeys.split("->")[0].split(":")[1] , DistancePositionKeys.split("->")[1].split(":")[1]])
    
    "gets a list of all x,y,z coordinate sets for a residue in a PDBXML file"
    def getPDBXMLPositions(self,pdbXMLLine):
        return re.compile("<p>(.+?)</p>").search(pdbXMLLine).group(1).split(";")
    "makes a 3-D floating point coordinate out of the string format numbers for a coordinate point"
    def makePoint(self,pdbXMLPointSplit):
        Ret = None
        if pdbXMLPointSplit[0] == "NA" or pdbXMLPointSplit[1] == "NA" or pdbXMLPointSplit[2] == "NA":
            pass
        else:
            Ret = [float(pdbXMLPointSplit[0]) , float(pdbXMLPointSplit[1]) , float(pdbXMLPointSplit[2])]
        return Ret
    
    "makes a 3-D point for the centroid coordinates"
    def getCentroidPoint(self):
        return self.makePoint(self.getPDBXMLPositions(pdbXMLLine)[0].split(","))
    "makes a 3-D point for the alpha carbon coordinates"
    def getAlphaCarbonPoint(self , pdbXMLLine):
        return self.makePoint(self.getPDBXMLPositions(pdbXMLLine)[1].split(","))
    "makes a 3-D point for the amine nitrogen coordinates"
    def getAminoNitrogenPoint(self,pdbXMLLine):
        return self.makePoint(self.getPDBXMLPositions(pdbXMLLine)[2].split(","))
    "makes a 3-D point for the carbonyl carbon coordinates"
    def getCarbonylCarbonPoint(self):
        return self.makePoint(self.getPDBXMLPositions(pdbXMLLine)[3].split(","))
    "makes a 3-D point for the side chain start coordinates"
    def getSideChainStartPoint(self):
        return self.makePoint(self.getPDBXMLPositions(pdbXMLLine)[4].split(","))
    "makes a 3-D point for the side chain end coordinates"
    def getSideChainEndPoint(self):
        return self.makePoint(self.getPDBXMLPositions(pdbXMLLine)[5].split(","))
    
    "computes a 1-D vector between two points in 1-D space"
    def computeVector(self , pointa , pointb):
        return pointb - pointa
    
    "computes the magnitude of distance between two 3-D points"
    def getDistanceMagnitude(self , pointa , pointb):
        return math.sqrt(math.pow(self.computeVector(pointa[0] , pointb[0]) , 2) +\
                         math.pow(self.computeVector(pointa[1] , pointb[1]) , 2) +\
                         math.pow(self.computeVector(pointa[2] , pointb[2]) , 2))
    
    "checks if the mutated site aligns with a real PDBXML residue, or a filler X residue"
    def checkIfMutationAssociatesToRealResidue(self , mutationXML):
        ret = False
        
        #gets accession and position of a PDBXML line
        accession_position = self.getAccessionPosition_L(mutationXML)
        accession = accession_position[0]
        position = accession_position[1]
        
        #if the position is in the PDBXML file
        if position in self.pdbAccessionsToPDBXMLFileContents_D[accession]['XMLResidue_D'].keys():
            pdbXMLLine = self.pdbAccessionsToPDBXMLFileContents_D[accession]['XMLResidue_D'][position]
            pdbXMLResidueType = self.getPDBXMLResidueType(pdbXMLLine)
            #if the residue type is "X", it is not a real residue
            if pdbXMLResidueType == "X":
                pass
            else:
                ret = True
        
        return ret
    
    "gets all accession-position lists for mutations provided that mutation associates to the PDB ID provided"    
    def getAllPositionKeysAccordingToAccession(self,AccKey):
        ret = []
        for mutationsXMLKey in self.mutationsXMLkey_L:
            if self.getAccessionPosition_L(self.mutationsXML_D[mutationsXMLKey]):
                for AccPos in self.getAccessionPosition_L(self.mutationsXML_D[mutationsXMLKey]):
                    if AccPos[0].lower() == AccKey:
                        ret.append(AccPos)
        return ret
        
        
        #return [self.getAccessionPosition_L(self.mutationsXML_D[mutationsXMLKey])[1] for mutationsXMLKey in self.mutationsXMLkey_L if self.getAccessionPosition_L(self.mutationsXML_D[mutationsXMLKey])[0] == AccKey]
    "gets all mutation XML lines provided those mutations associate to the PDB ID provided"
    def getAllMutationXMLAccordingToAccession(self,AccKey):
        ret = []
        for mutationsXMLKey in self.mutationsXMLkey_L:
            if self.getAccessionPosition_L(self.mutationsXML_D[mutationsXMLKey]):
                for AccPos in self.getAccessionPosition_L(self.mutationsXML_D[mutationsXMLKey]):
                    if AccPos[0].lower() == AccKey:
                        ret.append(self.mutationsXML_D[mutationsXMLKey])
        return ret
        #return [self.mutationsXML_D[mutationsXMLKey] for mutationsXMLKey in self.mutationsXMLkey_L if self.getAccessionPosition_L(self.mutationsXML_D[mutationsXMLKey])[0] == AccKey]
    "gets all mutation XML keys provided those keys associate to the PDB ID provided"
    def getAllMutationXMLKeysAccordingToAccession(self,AccKey):
        return [mutationsXMLKey for mutationsXMLKey in self.mutationsXMLkey_L if self.getAccessionPosition_L(self.mutationsXML_D[mutationsXMLKey])[0] == AccKey]
            
    "modifies the temperature factor of a PDB file so that when opened in pymol, it will be colored according to the mutations"
    def createPDBColoredFile(self, accession , outputPath):
        check_D = {}
        
        #checks if mutations are present
        if self.mutationsPresent:
            #gets the chain and position keys in the PDBXML file that correspond to mutated sites
            for mutationKey in self.mutationsXMLkey_L:
                PDBRes = re.compile("<R>(.+?)</R>").search(self.mutationsXML_D[mutationKey]).group(1)
                PDBAcc = PDBRes.split("|")[0]
                if PDBAcc == accession:
                    chain = re.compile("([A-Z]+)").search(PDBRes.split("|")[1]).group(1)
                    resNum = re.compile("([0-9]+)").search(PDBRes.split("|")[1]).group(1)
                    check_D[":%s.%s" % (resNum , chain)] = mutationKey
                    
        
        end_L = []
        
        #for each atom
        for atomKey in self.PDBContents_D[accession]['AtomKey_L']:
            resKey = re.compile("^(:.+?\..+?)-").search(atomKey).group(1)
            
            #if the atom key is in the dictionary of mutated sites, the temp factor is set to 1.0 (will appear red) and added to the list of lines
            #otherwise the temp factor will be set to 0.0 (will appear blue)
            
            if resKey in check_D.keys():
                appropriateScore_D = self.mutationScore_D[check_D[resKey]]

                newNum = 1.0
                if newNum < 0.0:
                    newNum = 0.0
                
                self.PDBContents_D[accession]['Atoms_D'][atomKey].Info['new_temp_factor'] = newNum
                self.PDBContents_D[accession]['Atoms_D'][atomKey].create_new_line()
                end_L.append(self.PDBContents_D[accession]['Atoms_D'][atomKey].Info['new_line'])
                
            else:
                self.PDBContents_D[accession]['Atoms_D'][atomKey].Info['new_temp_factor'] = 0.0
                self.PDBContents_D[accession]['Atoms_D'][atomKey].create_new_line()
                end_L.append(self.PDBContents_D[accession]['Atoms_D'][atomKey].Info['new_line'])
        
            
        #writes to the output path specified
        with open(outputPath , "w") as w:
            w.write('\n'.join(end_L))
    
        

