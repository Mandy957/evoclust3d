####################################################################################################
#                                                                                                  #
#           PROJECT Protein Adaptation                                                             #
#           CLASS WholeTreeOrthologousSubgroup                                                     #
#           PROGRAMMER Jeremy Adams                                                                #
#           STARTED 08-01-14                                                                       #
#           LASTMOD 07-08-14                                                                       #
#                                                                                                  #
#           DESCRIPTION Class that represents a single protein family with associated phylogenetic #
#                       tree and protein sequences. Ancestral sequences are reconstructed and      #
#                       properly placed here                                                       #
#                                                                                                  #
####################################################################################################

import re
import os
import sets

from fasequence import FASequence
from wholetreetriplealignmentset import WholeTreeTripleAlignmentSet
from staticmethods import *

class WholeTreeOrthologousSubgroup:
    
    """
    Class attributes:
    FastMLTree (FastMLTree class object): Tree object that parses newick tree and prepares it for mutation analysis in evolutionary stages
    LeafKey_L (List): List of terminal node keys in the FastMLTree
    FastaPATH (String): Path to user-defined sequence file in FASTA format
    FastaContent (String): String representing all of the content in the FASTA file
    FastaKey_L (List): List of strings representing the headers of sequences
    Fasta_D (Dict): Key is the Fasta sequence header, value is the sequence itself
    
    RepresentativeSequence (FASequence class object): Sequence within all reconstructed sequences that best represents the protein family: used to BLAST to PDB
    
    FastMLResults_D (Dict): Results of FastML ancestral reconstruction
    FastMLReconstructedFASTA (String): Fasta format string of all sequences (including ancestral) resulting from FastML
    FastMLProb_L (List): List of lines of reconstruction probabilites from FastML
    
    FastMLNodeNameToReconstructedSeq_D (Dict): Key is the FastML node name, value is string representing sequence at that point in evolutionary history
    TreeNodeToReconstructedSeq_D (Dict): Key is the Original (Cogent) node name, value is string representing sequence at that point in evolutionary history
    Probability_D (Dict): Key is Original node name, value is list of probabilities for each sequence index (residue) in the sequence
    
    WholeTreeTripleAlignmentSet (WholeTreeTripleAlignmentSet class object): performs lower-scale analysis, including the mapping of residues in the sequences to PDB structure residues
    ExitString (String): string to be used in final output 
    """
    
    "CONSTRUCTOR"
    def __init__(self, FastMLTree , FastaPATH , FastaKey_L , Fasta_D):
        self.FastMLTree = FastMLTree
        self.LeafKey_L = self.FastMLTree.LeafKey_L
        self.FastaPATH = FastaPATH
        self.FastaContent = open(FastaPATH,"r").read()
        self.FastaKey_L = FastaKey_L
        self.Fasta_D = Fasta_D
        
        #builds the consensus/representative sequence to be used for BLASTing
        self.RepresentativeSequence = self.buildConsensusSequence()
        
        #executes FastML to get reconstructed sequences and probabilities
        self.FastMLResults_D = executeFastML(self.FastaContent , self.FastMLTree.FastMLInputTreeString , False)
        self.FastMLReconstructedFASTA = self.FastMLResults_D['OutputFASTA']
        self.FastMLProb_L = self.FastMLResults_D['OutputProbLines_L']
        
        if self.FastMLReconstructedFASTA == "NOPATH":
            self.FastMLReconstructedFASTA = self.getFakeFASTA()
        
        #parses out the reconstructed sequences according to FastML node names and then places these sequences according to Original (Cogent) node names
        FH = getInputTempFile(self.FastMLReconstructedFASTA)
        ReadFastMLOutputFasta = readFasta(FH.name)
        self.FastMLNodeNameToReconstructedSeq_D = ReadFastMLOutputFasta[1]
        self.TreeNodeToReconstructedSeq_D = self.prepareNodeToItsSequenceDict()
        self.ReconstructedFASTAString = self.prepareReconstructedFASTAString()
        
        #get probability dict
        self.Probability = self.prepareNodeToItsProbabilityDict()
        self.Probability_D = self.Probability[0]
        self.ReconstructedProbabilityString = self.Probability[1]
        
        #performs lower-scale analysis
        self.WholeTreeTripleAlignmentSet = WholeTreeTripleAlignmentSet(self.RepresentativeSequence , self.TreeNodeToReconstructedSeq_D , self.FastMLTree.BranchKey_L , self.Probability_D)
        
        #gets output xml string and matrix graphics string
        self.ExitString = self.getExitString()
        self.MatrixGraphicsString = self.WholeTreeTripleAlignmentSet.MatrixGraphicsXML
    
    "retrieves the most representative sequence out of all extant sequences in the family"
    def buildConsensusSequence(self):
        potential_states = ['-','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X'] #residue possibilities
        
        sequences = []
        
        for key in sorted(self.FastaKey_L):
            sequences.append([key , self.Fasta_D[key]])
        
        consensus_sequence = []
        
        #for each state (residue) of all the sequences, the most prominent state is added into the consensus sequence
        for state in range(0,len(sequences[0][1])):
            states = []
            state_counts = {key : 0 for key in potential_states}
            for sequence in sequences:
                states.append(sequence[1][state])
            
            for state in states:
                state_counts[state] += 1
            
            consensus_state = ""
            consensus_count = 0
            for state_key in state_counts.keys():
                if state_counts[state_key] > consensus_count:
                    consensus_state = state_key
                    consensus_count = state_counts[state_key]
            
            consensus_sequence.append(consensus_state)
            
        #"fictional" sequence composing the most prominent residue at each state
        consensus_sequence = ''.join(consensus_sequence)
        
        
        identity_count_hash = {sequence[0] : 0 for sequence in sequences}
        #finds the real sequence that matches most closely to the fictional consensus sequence
        for sequence in sequences:
            key = sequence[0]
            seq = sequence[1]
            
            for i in range(0,len(consensus_sequence)):
        
                if seq[i] == consensus_sequence[i]:
                    identity_count_hash[key] += 1
        
        representative_key = ""
        representative_count = 0
        for key in identity_count_hash.keys():
            if identity_count_hash[key] > representative_count:
                representative_key = key
                representative_count = identity_count_hash[key]
        
        return FASequence(representative_key , self.Fasta_D[representative_key])
    
    def getFakeFASTA(self):
        ret = []
        
        for NodeKey in self.FastMLTree.FastMLToOriginalMatchedNodes_D.keys():
            ret.append(">%s\n%s" % (NodeKey,self.Fasta_D[self.FastaKey_L[0]]))
        return "\n".join(ret)
        
    
    "moves reconstructed sequences under the FastML header to reconstructed sequences under the appropriate Cogent header"
    def prepareNodeToItsSequenceDict(self):
        return {self.FastMLTree.FastMLToOriginalMatchedNodes_D[FastMLNodeKey] : self.FastMLNodeNameToReconstructedSeq_D[FastMLNodeKey] \
                for FastMLNodeKey in self.FastMLTree.FastMLToOriginalMatchedNodes_D.keys()}
    
    "prepares a FASTA format file containing sequence node headers and the corresponding sequence"
    def prepareReconstructedFASTAString(self):
        L = []
        for NodeKey in self.FastMLTree.NodeKey_L:
            L.append(">"+NodeKey)
            L.append(re.sub("(.{60})" , r"\1\n" , self.TreeNodeToReconstructedSeq_D[NodeKey]))
        return "\n".join(L)
        
        
    
    "get dictionary of the residue type as the key and the column in which it appears for the probability listings of FastML"
    def getResidueToProbLineIndexDict(self):
        return {"A":2,"C":3,"D":4,"E":5,"F":6,"G":7,"H":8,"I":9,"K":10,"L":11,"M":12,"N":13,"P":14,"Q":15,"R":16,"S":17,"T":18,"V":19,"W":20,"Y":21}
    
    "gets dictionary of the reliability/confidence of each reconstructed state for each state in an ancestral sequence"
    def prepareNodeToItsProbabilityDict(self):
        ResidueToProbLineIndex_D = self.getResidueToProbLineIndexDict() #gets residue to index dict
        ret = {}
        retLines_L = [self.FastMLProb_L[0].replace("\n","")]
        #for each line from the file
        for line in self.FastMLProb_L[1:]:
            lineSplit = line.replace("\n","").split(",")
            OriginalNodeName = self.FastMLTree.FastMLToOriginalMatchedNodes_D[lineSplit[0]]
            
            if OriginalNodeName in ret.keys():
                pass
            else:
                ret[OriginalNodeName] = []
            
            #gets the residue at the proper index, then finds the reconstruction probability of that residue and adds it to the list for that node
            toAppend = ""
            
            sequenceIndex = int(lineSplit[1]) - 1
            residue = self.TreeNodeToReconstructedSeq_D[OriginalNodeName][sequenceIndex]
            
            if residue == "-":
                pass
            else:
                toAppend = lineSplit[ResidueToProbLineIndex_D[residue]]
            
            ret[OriginalNodeName].append(toAppend)
            retLines_L.append(",".join([self.FastMLTree.FastMLToOriginalMatchedNodes_D[lineSplit[0]] , ",".join(lineSplit[1:])]))
        
        #leaf sequences are considered to have 100% confidence
        for LeafKey in self.LeafKey_L:
            ret[LeafKey] = [1.0 for i in range(0,len(self.TreeNodeToReconstructedSeq_D[LeafKey]))]
        
        return [ret,"\n".join(retLines_L)]
    
    "gets output string for final XML output"        
    def getExitString(self):
        ret = ["\t\t<OSG>\n" , "\t\t\t<Sequences>\n"] #initial lines
        ret = ret + ["\t\t\t\t<Seq><H>%s</H><S>%s</S></Seq>\n" % (key , self.TreeNodeToReconstructedSeq_D[key]) for key in self.FastMLTree.NodeKey_L] #all reconstructed sequences
        ret = ret + ["\t\t\t</Sequences>\n" , "\t\t\t<Alignment_set>\n" , self.WholeTreeTripleAlignmentSet.ExitString , "\t\t\t</Alignment_set>\n" , "\t\t</OSG>\n"] #closing lines with the bulk of the lower-scale analysis
        return "".join(ret)
            