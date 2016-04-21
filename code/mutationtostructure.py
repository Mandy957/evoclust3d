from staticmethods import *
from fasequence import FASequence
from referencetopdb2dscoringmatrix import ReferenceToPDB2DScoringMatrix
from triplealignment import TripleAlignment

class MutationToStructure:
    
    def __init__(self,FastaPATH, InterestHeader, IngroupHeader, OutgroupHeader):
        
        Fasta = readFasta(FastaPATH)
        self.interestSeq = Fasta[1][InterestHeader]
        self.ingroupSeq = Fasta[1][IngroupHeader]
        self.outgroupSeq = Fasta[1][OutgroupHeader]
        
        interestSeqObj = FASequence(InterestHeader , self.interestSeq)
        ingroupSeqObj = FASequence(IngroupHeader , self.ingroupSeq)
        outgroupSeqObj = FASequence(OutgroupHeader , self.outgroupSeq)
        
        self.RepresentativeSeqWithGaps = interestSeqObj
        
        GapIndicesAndWithoutGapSeq_D = self.getListOfGapIndices()
        self.RepresentativeSeqWithoutGaps = GapIndicesAndWithoutGapSeq_D["WithoutGapSeq"]
        self.GapIndices = GapIndicesAndWithoutGapSeq_D["GapIndices"]
        
        self.RepresentativeSeqWithoutGaps.completeMasterList()
        
        self.RepresentativeSeqWithGaps.PDBMasterList = self.gapCorrectPDBMasterList(self.RepresentativeSeqWithoutGaps.PDBMasterList , self.GapIndices , len(self.RepresentativeSeqWithGaps.info['sequence']))
        
        self.RepresentativeSeqScoringMatrix = ReferenceToPDB2DScoringMatrix(self.RepresentativeSeqWithGaps.referenceSequencePrint())
            
        self.PDBFile_L = self.getPDBFile_L()
        
        self.MutationString = self.getMutationString()
        
        self.triplealn = TripleAlignment(self.MutationString,\
                                         self.RepresentativeSeqScoringMatrix,\
                                         self.PDBFile_L,\
                                         "")
        
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
    
    def getMutationString(self):
        
        mutationString = []
        
        for i in range(0,len(self.interestSeq)):
            interestS = self.interestSeq[i]
            ingroupS = self.ingroupSeq[i]
            outgroupS = self.outgroupSeq[i]
            
            if interestS == "-" or ingroupS == "-" or outgroupS == "-":
                pass
            else:
                if ingroupS == outgroupS:
                    if interestS != ingroupS:
                        mutationString.append("%s|%s|%s|%s|%s" % (str(i+1) , outgroupS , interestS , "1.0" , "1.0"))
                        
        
        ret = ""
        #builds a single string out of the whole list of mutations
        if len(mutationString) == 0:
            ret = None
        else:
            ret = ';'.join(mutationString)
                    
        return ret  
        
        
        
          
        