####################################################################################################
#                                                                                                  #
#           PROJECT Protein Adaptation                                                             #
#           CLASS WholeTreeTripleAlignmentSet                                                      #
#           PROGRAMMER Jeremy Adams                                                                #
#           STARTED 08-01-14                                                                       #
#           LASTMOD 06-11-14                                                                       #
#                                                                                                  #
#           DESCRIPTION Class that aligns the representative sequence to PDB structures, and       #
#                       instantiates mutation cataloguing and analysis for each branch segment     #
#                                                                                                  #
####################################################################################################

import re

from fasequence import FASequence
from adalignment import ADAlignment
from referencetopdb2dscoringmatrix import ReferenceToPDB2DScoringMatrix
from triplealignment import TripleAlignment

from staticmethods import *

class WholeTreeTripleAlignmentSet:
    
    """
    Class attributes:
    RepresentativeSeqWithGaps (FASequence Obj): Sequence containing information on the most representative extant sequence of the entire set
    TreeNodeToReconstructedSeq_D (Dict): Key is the tree node name, value is the sequence at that node
    BranchKey_L (List): list of all branch segments denoted by ancestral and derived branch keys
    Probability_D (Dict): dictionary of reconstruction probabilities for all nodes
    
    RepresentativeSeqWithoutGaps (FASequence Obj): Same sequence as the with gap sequence although with gaps removed
    RepresentativeSeqScoringMatrix (Reference2PDBScoringMatrix Obj): 2-D grid object of the RepresentativeSeqWithGaps object
    PDBFile_L (list): list of PDB ID accessions used for this protein family
    TripleAlignments_D (Dict): Key is the branch key, value is a TripleAlignment object representing that branch
    
    ExitString (String): Report XML for all branches
    MatrixGraphicsXML (String): Scoring keys for PDB alignment boundaries
    """
    
    "CONSTRUCTOR"
    def __init__(self , RepresentativeSeqWithGaps , TreeNodeToReconstructedSeq_D , BranchKey_L , Probability_D):
        
        #variables from WholeTreeOrthologousSubgroup
        self.RepresentativeSeqWithGaps = RepresentativeSeqWithGaps
        self.TreeNodeToReconstructedSeq_D = TreeNodeToReconstructedSeq_D
        self.BranchKey_L = BranchKey_L
        self.Probability_D = Probability_D
        
        #gets gaps from the with gap seq and prepares the sequence without gaps
        GapIndicesAndWithoutGapSeq_D = self.getListOfGapIndices()
        self.RepresentativeSeqWithoutGaps = GapIndicesAndWithoutGapSeq_D["WithoutGapSeq"]
        self.GapIndices = GapIndicesAndWithoutGapSeq_D["GapIndices"]
        
        #BLASTs the sequence without gaps to the PDB
        self.RepresentativeSeqWithoutGaps.completeMasterList()
        
        #adds gaps to the with gap PDBMasterList and gets the scoring matrix for it
        self.RepresentativeSeqWithGaps.PDBMasterList = self.gapCorrectPDBMasterList(self.RepresentativeSeqWithoutGaps.PDBMasterList , self.GapIndices , len(self.RepresentativeSeqWithGaps.info['sequence']))
        
        self.RepresentativeSeqScoringMatrix = ReferenceToPDB2DScoringMatrix(self.RepresentativeSeqWithGaps.referenceSequencePrint())
            
        self.PDBFile_L = self.getPDBFile_L()
        
        #performs TripleAlignments for each branch in the overall tree
        self.TripleAlignments_D = {}
        for branchKey in self.BranchKey_L:
            
            branchSplit = branchKey.split(">>")
            
            if branchSplit[0] in self.TreeNodeToReconstructedSeq_D.keys() and branchSplit[1] in self.TreeNodeToReconstructedSeq_D.keys():
                adAlignment = ADAlignment(self.TreeNodeToReconstructedSeq_D[branchSplit[0]] , self.TreeNodeToReconstructedSeq_D[branchSplit[1]] , Probability_D[branchSplit[0]] , Probability_D[branchSplit[1]])
                self.TripleAlignments_D[branchKey] = TripleAlignment(adAlignment.mutations,\
                                                     self.RepresentativeSeqScoringMatrix,\
                                                     self.PDBFile_L,\
                                                     "\t\t\t\t\t\t\t")
                
        #prepares report XML and scoring matrix XML
        self.ExitString = self.getExitString()
        self.MatrixGraphicsXML = self.RepresentativeSeqScoringMatrix.scoringMatrixXMLPrint()
    
    "gets list of gaps and gap lengths in the sequence to prepare for reinsertion"
    def getListOfGapIndices(self):
        
        seq_with_gaps = self.RepresentativeSeqWithGaps.info['sequence']
        missing_gaps_in_the_without_gap_seq = []
        
        without_gap_seq = []
        
        a_gap_has_been_found = False
        the_gap_has_ended = False
        gap_start = 0
        
        #for each state in the sequence
        for i in range(0,len(seq_with_gaps)):
            #if in gap mode
            if a_gap_has_been_found:
                if seq_with_gaps[i] != "-": #ends the gap when the "-" is not found
                    gap_end = i
                    
                    missing_gaps_in_the_without_gap_seq.append([gap_start , gap_end , len(without_gap_seq)])
                    without_gap_seq.append(seq_with_gaps[i])
                    a_gap_has_been_found = False
            
                elif i == len(seq_with_gaps) - 1:
                    missing_gaps_in_the_without_gap_seq.append([gap_start , len(seq_with_gaps) , len(without_gap_seq)])
                
            #if not in gap mode and looking for a gap
            else:
                if seq_with_gaps[i] == "-": #starts the gap when the "-" is found
                    gap_start = i
                    a_gap_has_been_found = True
        
                else:
                    without_gap_seq.append(seq_with_gaps[i])
                   
        return {'WithoutGapSeq' : FASequence('seq_without_gaps' , ''.join(without_gap_seq)),\
                'GapIndices' : missing_gaps_in_the_without_gap_seq}
    
    "reinsert gaps and adjust alignment coordinates appropriately in the PDB sequence"
    def gapCorrectPDBMasterList(self , pdb_master_list , gap_indices , query_length):
        
        corrected_list = []
        
        #for each target sequence in the master list
        for pdb_hit in pdb_master_list:
            new_qseq = []
            new_sseq = []
            old_seq_index_counter = 0
            
            #gets old information from BLASTing with the without gap seq
            old_qseq = get_qseq(pdb_hit)
            old_sseq = get_sseq(pdb_hit)
            
            old_qstart = int(get_qstart(pdb_hit))
            old_seq_blast_numbering_index_counter = old_qstart
            old_qend = int(get_qend(pdb_hit))
            
            #for each gap found in the with gap seq
            #adds gaps to the sseq (PDB) and qseq (representative sequence without gaps) wherever there were gap indices
            for gap in gap_indices:
                gap_position_in_the_without_gap_seq = gap[2] + 1
                
                if gap_position_in_the_without_gap_seq >= old_qstart and gap_position_in_the_without_gap_seq <= old_qend:
                    blast_alignment_seq_before_gaps = old_qseq[old_seq_index_counter:old_seq_index_counter+(gap_position_in_the_without_gap_seq-old_seq_blast_numbering_index_counter)]
                    gaps_in_the_blast_alignment = blast_alignment_seq_before_gaps.count("-")
                    
                    if gaps_in_the_blast_alignment == 0:
                        pass
                    else:
                        blast_alignment_seq_is_not_satisfied = True
                        
                        while blast_alignment_seq_is_not_satisfied:
                            blast_alignment_candidate = old_qseq[old_seq_index_counter:old_seq_index_counter+(gap_position_in_the_without_gap_seq-old_seq_blast_numbering_index_counter)+gaps_in_the_blast_alignment]
                            new_amount_of_gaps_in_the_blast_alignment = blast_alignment_candidate.count("-")
                            
                            if new_amount_of_gaps_in_the_blast_alignment == gaps_in_the_blast_alignment:
                                blast_alignment_seq_is_not_satisfied = False
                                blast_alignment_seq_before_gaps = blast_alignment_candidate
                            else:
                                gaps_in_the_blast_alignment = new_amount_of_gaps_in_the_blast_alignment
                        
                        
                    new_qseq.append(blast_alignment_seq_before_gaps)
                    new_sseq.append(old_sseq[old_seq_index_counter:old_seq_index_counter+(gap_position_in_the_without_gap_seq-old_seq_blast_numbering_index_counter)+gaps_in_the_blast_alignment])
                    
                    new_qseq.append("-" * (gap[1] - gap[0]))
                    new_sseq.append("-" * (gap[1] - gap[0]))
                    
                    new_qseq.append(old_qseq[old_seq_index_counter+(gap_position_in_the_without_gap_seq-old_seq_blast_numbering_index_counter)])
                    new_sseq.append(old_sseq[old_seq_index_counter+(gap_position_in_the_without_gap_seq-old_seq_blast_numbering_index_counter)])
                    
                    old_seq_index_counter = old_seq_index_counter+(gap_position_in_the_without_gap_seq-old_seq_blast_numbering_index_counter)+gaps_in_the_blast_alignment+1
                        
                    old_seq_blast_numbering_index_counter = gap_position_in_the_without_gap_seq+1
            
            new_qseq.append(old_qseq[old_seq_index_counter:])
            new_sseq.append(old_sseq[old_seq_index_counter:])
            
            #places the new sequence lists into a single string
            new_qseq = ''.join(new_qseq)
            new_sseq = ''.join(new_sseq)
            
            
            #finds different start and end coordinates for the with gap seq alignments by adjusting the original without gap seq start and end coordinates
            new_qstart_not_found = True
            new_qend_not_found = True
            
            new_qstart = 0
            new_qend = 0
            new_sstart = 0
            new_send = 0
            
            if len(gap_indices) > 0 :
            
                for i in range(0,len(gap_indices)):
                    if new_qstart_not_found:
                        if (gap_indices[i][2] + 1) >= old_qstart:
                            for a in range(0,i):
                                new_qstart += (gap_indices[a][1] - gap_indices[a][0])
                            
                            new_qstart += old_qstart
                            
                            new_qstart_not_found = False
                            
                    else:
                        if new_qend_not_found:
                            if (gap_indices[i][2] + 1) >= old_qend:
                                for a in range(0,i):
                                    new_qend += (gap_indices[a][1] - gap_indices[a][0])
                                
                                new_qend += old_qend
                                new_qend_not_found = False
            else:
                
                new_qstart = old_qstart
                new_qend = old_qend
                        
            
            #creates a similar line as in the without gap seq, but with the correctly adjusted gaps
            list_item = []
            list_item.append(get_qcovs(pdb_hit))
            list_item.append(get_evalue(pdb_hit))
            list_item.append(get_bitscore(pdb_hit))
            list_item.append(get_sseqid(pdb_hit))
            list_item.append(str(query_length))
            list_item.append(str(new_qstart))
            list_item.append(str(new_qend))
            list_item.append(new_qseq)
            list_item.append(get_slen(pdb_hit))
            list_item.append(get_sstart(pdb_hit))
            list_item.append(get_send(pdb_hit))
            list_item.append(new_sseq)
            list_item.append(get_qcovhsp(pdb_hit))
            
            corrected_list.append(','.join(list_item))
                
        return corrected_list
    
    "gets a list of PDB ID accessions from the overall sequence alignment"
    def getPDBFile_L(self):
        
        ret = []
        
        #reads the processed XML of the scoring matrix and finds all unique PDB IDs
        if self.RepresentativeSeqScoringMatrix.ProcessedXML:
            for pdb in self.RepresentativeSeqScoringMatrix.ProcessedXML['pdb_files']:
                
                accession = pdb['pdb_id'].split("|")[0].lower()
                this_file_has_not_yet_been_read = True
                
                for key in ret:
                    if key == accession:
                        this_file_has_not_yet_been_read = False
                
                if this_file_has_not_yet_been_read:
                    ret.append(accession)
        
        return ret
    
    "gets the main report XML string using all of the TripleAlignment reports"
    def getExitString(self):
        
        #adds the Branch tag to each TripleAlignment report string
        Ret = ["\t\t\t\t<Branches>\n"]
            
        for branchKey in self.BranchKey_L:
            Ret = Ret + ["\t\t\t\t\t<Branch>\n" , "\t\t\t\t\t\t<Branch_name>%s</Branch_name>\n" % (branchKey) , "\t\t\t\t\t\t<Muts>\n" ,\
                         self.TripleAlignments_D[branchKey].exit['exit_string'] , "\t\t\t\t\t\t</Muts>\n" , "\t\t\t\t\t</Branch>\n"]
                 
        Ret.append("\t\t\t\t</Branches>\n")
        return "".join(Ret)
          