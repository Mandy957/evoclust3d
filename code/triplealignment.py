####################################################################################################
#                                                                                                  #
#           PROJECT Protein Adaptation                                                             #
#           CLASS TripleAlignment                                                                  #
#           PROGRAMMER Jeremy Adams                                                                #
#           STARTED 25-07-13                                                                       #
#           LASTMOD 06-11-14                                                                       #
#                                                                                                  #
#           DESCRIPTION class representation of an ancestral, derived, alignment pair with all PDB #
#                       matchings. Used to align mutations to the PDB structure and get mutation   #
#                       to structure matches.                                                      #
#                                                                                                  #
####################################################################################################

import re

class TripleAlignment:
    
    "CONSTRUCTOR"
    def __init__(self , mutation_string , scoring_matrix , pdbfile_L , totaltabs):
        
        """
        Class attributes:
        exit (Dict): contains information on whether the alignment process was succesfully conducted, and the output alignment string
        """
        
        self.exit = self.initialize(mutation_string , scoring_matrix , pdbfile_L , totaltabs)
        
    "main execution of the class"
    def initialize(self , mutation_string , scoring_matrix , pdbfile_L , totaltabs):
        
        exit_status = False
        
        #output string for this branch segment
        exit_string = []
        exit_string.append("%s<PDBs>%s</PDBs>" % (totaltabs , ';'.join(pdbfile_L)))
        
        #if there are no mutations from the AD alignment, then nothing is done here
        if mutation_string == None:
            exit_string.append(totaltabs+"<Mut>No Mutations</Mut>\n")
            exit_string = '\n'.join(exit_string)
        else:
            
            #gets reconstruction quality from the AD alignment
            qualityAdjustedMutations_L = self.completeQualityInformation(mutation_string)
            
            #does not execute if there are no good mutations
            if qualityAdjustedMutations_L == None:
                pass
            else:
                
                #maps all mutations to the aligned PDB structures
                allMutationsWithPDBMatchedResidues = self.mapAllMutationsToPDBResidues(qualityAdjustedMutations_L , scoring_matrix , pdbfile_L)
                
                #does not execute further if there are no mutations 
                if allMutationsWithPDBMatchedResidues == None:
                    
                    exit_string.append(totaltabs+"<Mut>SAS Information Unclear</Mut>\n")
                    exit_string = '\n'.join(exit_string)
                    
                else:
                    
                    exit_status = True
                    
                    allMutationsWithPDBMatchedResidues = self.positionCorrectMutations(allMutationsWithPDBMatchedResidues , scoring_matrix , pdbfile_L)
                    
                    
                    for mutation in allMutationsWithPDBMatchedResidues:
                        exit_string.append(totaltabs+mutation)
                    
                    exit_string = '\n'.join(exit_string)+"\n"
                """    
                    #completes the distance information for final completion of mutation XML data
                    finalMutations = self.completeDistanceInformation(allMutationsWithPDBMatchedResidues , scoring_matrix , pdbfile_L)
                    
                    
                    for mutation in finalMutations:
                        exit_string.append(totaltabs+mutation)
                    
                    exit_string = '\n'.join(exit_string)+"\n"
                """
        
        #sets the status and final string as the class exit variable
        return {'exit_status' : exit_status,\
                'exit_string' : exit_string}
    
    "get quality information for the group of mutations"
    def completeQualityInformation(self , mutation_string):
        
        qualityAdjustedMutations_L = []
        b = "OK"
        
        if mutation_string == None:
            b = None
        else:
            
            #puts quality as U (unsat) if probabilities are lower than 0.75, or S (sat) if higher
            allMutations_L = mutation_string.split(";")
            for mutation in allMutations_L:
                mutationSplit = mutation.split("|")
                quality = "U"
                
                if float(mutationSplit[3]) > 0.75 and float(mutationSplit[4]) > 0.75:
                    quality = "S"
                
                qualityAdjustedMutations_L.append("<M>%s|%s%s</M><Q>%s</Q>" % (mutationSplit[0],mutationSplit[1],mutationSplit[2],quality))
            
        ret = None
        
        if b == "OK":
            ret = qualityAdjustedMutations_L
        return ret
    
    "gets a dictionary of important coordinate mappings between the AD pair and the PDB sequence"
    def getStatsFromCorrectRefSeqXML(self, correct_pdb_accession_chain , refseq_xml):
        
        all_refseq_xmls = re.findall("<PDB_alignment>.+?</PDB_alignment>" , refseq_xml , re.DOTALL)
        
        correct_refseq_xml_not_yet_found = True
        correct_refseq_xml = ""
        
        #gets the correct refseq XML
        for ref_seq_xml in all_refseq_xmls:
            if correct_refseq_xml_not_yet_found:
                refseq_xml_pdb_accession_chain = re.compile("<PDB_id>(.+?)</PDB_id>").search(ref_seq_xml).group(1)
                
                
                if correct_pdb_accession_chain == refseq_xml_pdb_accession_chain:
                    correct_refseq_xml_not_yet_found = False
                    correct_refseq_xml = ref_seq_xml
        
        #gets the correct positions and sequence string for the reference (representative) sequence
        refseq_start = int(re.compile("<Alignment_start_query>(.+?)</Alignment_start_query>").search(correct_refseq_xml).group(1))
        refseq_end = int(re.compile("<Alignment_end_query>(.+?)</Alignment_end_query>").search(correct_refseq_xml).group(1))
        refseq = re.compile("<Aligned_query_sequence>(.+?)</Aligned_query_sequence>").search(correct_refseq_xml).group(1)
        #gets the correct positions and sequence string for the PDB target sequence
        pdbseq_start = int(re.compile("<Alignment_start_subject>(.+?)</Alignment_start_subject>").search(correct_refseq_xml).group(1))
        pdbseq_end = int(re.compile("<Alignment_end_subject>(.+?)</Alignment_end_subject>").search(correct_refseq_xml).group(1))
        pdbseq = re.compile("<Aligned_subject_sequence>(.+?)</Aligned_subject_sequence>").search(correct_refseq_xml).group(1)
        
        return {'refseq_start' : refseq_start,\
                'refseq_end' : refseq_end,\
                'refseq' : refseq,\
                'pdbseq_start' : pdbseq_start,\
                'pdbseq_end' : pdbseq_end,\
                'pdbseq' : pdbseq}
    
    "gets the correct PDB residue for a position on the representative sequence"
    def getCorrectPositionAndResidueInPDBSeq(self, xmlseq_hash , position):
        refseq_start = xmlseq_hash['refseq_start']
        refseq_end = xmlseq_hash['refseq_end']
        refseq = xmlseq_hash['refseq']
        pdbseq_start = xmlseq_hash['pdbseq_start']
        pdbseq_end = xmlseq_hash['pdbseq_end']
        pdbseq = xmlseq_hash['pdbseq']
        
        pdbstep = 0
                    
        pdbresidue = ""
        pdb_correct_position = ""
        
        #continues along sequence until the starting position of the alignment to the PDB seq is reached
        for i in range(0,len(refseq)):
            if (i + refseq_start) == position:
                if pdbseq[i] == "-":
                    
                    #if the mutated site aligns to a gap in the PDB sequence, then it will find the closest non-gap residue
                    forward_limit_has_not_been_found = True
                    forward_count = 0
                    while forward_limit_has_not_been_found:
                        if pdbseq[i+forward_count] == "-":
                            forward_count += 1
                        else:
                            forward_limit_has_not_been_found = False
                    
                    reverse_limit_has_not_been_found = True
                    reverse_count = 0
                    while reverse_limit_has_not_been_found:
                        if pdbseq[i-reverse_count] == "-":
                            reverse_count += 1
                        else:
                            reverse_limit_has_not_been_found = False
                    
                    if reverse_count > forward_count:
                        pdbresidue = pdbseq[i+forward_count]
                        pdb_correct_position = pdbseq_start + pdbstep + 1
                    else:
                        pdbresidue = pdbseq[i-reverse_count]
                        pdb_correct_position = pdbseq_start + pdbstep - 1
                    
                    
                #the correct position on the PDB sequence when it does not align to a gap in the PDB sequence
                else:
                    pdbresidue = pdbseq[i]
                    pdb_correct_position = pdbseq_start + pdbstep
            
            #increments for non gap positions
            if pdbseq[i] == "-":
                pass
            else:
                pdbstep += 1
        
        Success = True
        if pdb_correct_position == "":
            Success = False
        
        return {'success' : Success,\
                'pdb_residue' : pdbresidue,\
                'pdb_correct_position' : pdb_correct_position}
    
    "maps all mutated sites to residues and gets PDB information"
    def mapAllMutationsToPDBResidues(self, qualityAdjustedMutations_L , scoring_matrix , pdb_file_hash):
        
        allMutationsWithPDBMatchedResidues_L = []
        
        #for all mutations
        for mutation in qualityAdjustedMutations_L:
            
            #gets the position of the mutation as an integer
            mutationStats = re.compile("<M>(.+?)</M>").search(mutation).group(1)
            mutationStatsSplit = mutationStats.split("|")
            
            position = int(mutationStatsSplit[0])
            position_index = (int(position)-1)
            from_state = mutationStatsSplit[1][0]
            to_state = mutationStatsSplit[1][1]
            
            quality = re.compile("<Q>(.+?)</Q>").search(mutation).group(1)
            
            
            correct_pdb_accession_chain = scoring_matrix.Matrix[position_index][position_index]['pdb']
            other_pdb_accession_chains = []
            
            for cell in scoring_matrix.Matrix[position_index]:
                if cell['pdb'] != "none":
                    if cell['pdb'] != correct_pdb_accession_chain:
                        if cell['pdb'] in set(other_pdb_accession_chains):
                            pass
                        else:
                            other_pdb_accession_chains.append(cell['pdb'])
            
            #if mutation does not associate to a PDB
            if correct_pdb_accession_chain == 'none':
                
                allMutationsWithPDBMatchedResidues_L.append("<M>%s|%s%s</M><Q>%s</Q><R>NOCOVERAGE</R>" % (str(position) , from_state , to_state , quality))
                
            #if mutation does associate to a PDB
            else:
                Res_L = []
                #gets position on the PDB sequence that matches the position on the AD sequence
                xmlseqhash = self.getStatsFromCorrectRefSeqXML(correct_pdb_accession_chain , scoring_matrix.ReferenceSequenceXML)
                pdbposition_hash = self.getCorrectPositionAndResidueInPDBSeq(xmlseqhash , position)
                
                if pdbposition_hash['success']:
                
                    #solvent_accessibility information and general position information
                    accession = correct_pdb_accession_chain.split("|")[0].lower()
                    chain = correct_pdb_accession_chain.split("|")[1]
                    Res_L.append(correct_pdb_accession_chain+str(pdbposition_hash['pdb_correct_position']))
                
                for pdb_accession_chain in other_pdb_accession_chains:
                    xmlseqhash = self.getStatsFromCorrectRefSeqXML(pdb_accession_chain , scoring_matrix.ReferenceSequenceXML)
                    pdbposition_hash = self.getCorrectPositionAndResidueInPDBSeq(xmlseqhash , position)
                    
                    if pdbposition_hash['success']:
                    
                        #solvent_accessibility information and general position information
                        accession = pdb_accession_chain.split("|")[0].lower()
                        chain = pdb_accession_chain.split("|")[1]
                        Res_L.append(pdb_accession_chain+str(pdbposition_hash['pdb_correct_position']))
                
                if len(Res_L) > 0:
                    
        
                    allMutationsWithPDBMatchedResidues_L.append("<M>%s|%s%s</M><Q>%s</Q><R>%s</R>" % (str(position) , from_state , to_state, quality,\
                                                                                                    ",".join(Res_L)))
                else:
                    allMutationsWithPDBMatchedResidues_L.append("<M>%s|%s%s</M><Q>%s</Q><R>NOCOVERAGE</R>" % (str(position) , from_state, to_state , quality))
        
        return allMutationsWithPDBMatchedResidues_L
    
    "method to readjust aligned residues to the PDB structure when there are overlapping residues"
    def positionCorrectMutations(self, allMutationsWithPDBMatchedResidues , scoring_matrix , pdbfile_L):
        MutKeys_L = []
        PDBAccToMutKeysToPositions_D = {}
        
        #looks at all mutations and excludes the mutations without coverage by the PDB
        for Mut in allMutationsWithPDBMatchedResidues:
            if re.compile("NOCOVERAGE").search(Mut):
                pass
            else:
                #get the mutation position and the position it currently aligns to in the PDB
                MutKey = re.compile("<M>(.+?)</M>").search(Mut).group(1).split("|")[0]
                MutKeys_L.append(MutKey)
                AllRes_L = re.compile("<R>(.+?)</R>").search(Mut).group(1).split(",")
                
                for Res in AllRes_L:
                    Accession = Res.split("|")[0]
                    Position = Res.split("|")[1]
                
                    
                    #adds to the dictionary where the key is the mutation position on the primary sequence, and the value is the position on the PDB structure
                    if Accession in PDBAccToMutKeysToPositions_D.keys():
                        pass
                    else:
                        PDBAccToMutKeysToPositions_D[Accession] = {}
                    PDBAccToMutKeysToPositions_D[Accession][MutKey] = Position
        
        #does the analysis for each aligned PDB structure
        for AccKey in PDBAccToMutKeysToPositions_D.keys():
            RelevantMutKeys_L = []
            for MutKey in MutKeys_L:
                if MutKey in PDBAccToMutKeysToPositions_D[AccKey].keys():
                    RelevantMutKeys_L.append(MutKey)
            
            #print AccKey
            #print len(set(RelevantMutKeys_L))
            #print len(set([PDBAccToMutKeysToPositions_D[AccKey][MutKey] for MutKey in RelevantMutKeys_L]))
            #print "*" * 50
            
            
            #if the set of mutations is equal to the set of aligned PDB residues, then there is no overlap at all
            if len(set(RelevantMutKeys_L)) == len(set([PDBAccToMutKeysToPositions_D[AccKey][MutKey] for MutKey in RelevantMutKeys_L])):
                pass
            else: #if there is overlap
                NeedToAdjust = True
                print AccKey
                print NeedToAdjust
                PositionTest_L = []
                NeedToAdjustPositions = []
                
                #gets the positions where there is overlap and adds them to a list
                for RMutKey in RelevantMutKeys_L:
                    if PDBAccToMutKeysToPositions_D[AccKey][RMutKey] in set(PositionTest_L):
                        if PDBAccToMutKeysToPositions_D[AccKey][RMutKey] in set(NeedToAdjustPositions):
                            pass
                        else:
                            NeedToAdjustPositions.append(PDBAccToMutKeysToPositions_D[AccKey][RMutKey])
                    else:
                        PositionTest_L.append(PDBAccToMutKeysToPositions_D[AccKey][RMutKey])
                
                #for each position that needs to be adjusted, the following method takes place
                for AdjustPosition in NeedToAdjustPositions:
                    ThereIsOverlap = True
                    NewAdjustPosition = AdjustPosition
                    
                    #finds the number of residues that are overlapping in the region of interest
                    nNeedsReadjusts = len([PDBAccToMutKeysToPositions_D[AccKey][RMutKey] for RMutKey in RelevantMutKeys_L if PDBAccToMutKeysToPositions_D[AccKey][RMutKey] == NewAdjustPosition])
                    #if the readjust has been correct by a previous readjust, then this will correct for that
                    if nNeedsReadjusts <= 1:
                        ThereIsOverlap = False
                    
                    #variables to check the status of the loop
                    AdjustRange = []
                    Encapsulated_L = []
                    ReadjustDone = False
                    
                    #continuously checks if there are still residues in the overlap range
                    while ThereIsOverlap:
                        print "stuck in the new while loop"
                        #the adjust range is the number of residues that need to be adjusted spread out on either end from the central adjust position
                        AdjustN = int(re.compile("([0-9]+)").search(NewAdjustPosition).group(1))
                        LowerLimit = AdjustN - (nNeedsReadjusts / 2)
                        UpperLimit = AdjustN + (nNeedsReadjusts /2)
                        
                        Encapsulated_L = []
                        
                        #checks how many residues are encapsulated by the adjust range
                        for RMutKey in RelevantMutKeys_L:
                            TestPos = int(re.compile("([0-9]+)").search(PDBAccToMutKeysToPositions_D[AccKey][RMutKey]).group(1))
                            if TestPos >= LowerLimit and TestPos <= UpperLimit:
                                Encapsulated_L.append(RMutKey)
                        #if the number of residues in the adjust range is the same as the original number that needed readjusting, then this readjust is finished
                        if len(Encapsulated_L) == nNeedsReadjusts:
                            ThereIsOverlap = False
                            AdjustRange = [LowerLimit , UpperLimit]
                            ReadjustDone = True
                        #if more residues are encapsulated by the adjust range, then we need to reset the numbers and iterate over the while loop again
                        else:
                            NewLowerLimit = PDBAccToMutKeysToPositions_D[AccKey][Encapsulated_L[0]]
                            NewUpperLimit = PDBAccToMutKeysToPositions_D[AccKey][Encapsulated_L[-1]]
                            
                            Chain = None
                            NewL = None
                            NewU = None
                            
                            NormalChainSearch = re.compile("([a-zA-Z]+)").search(NewLowerLimit)
                            
                            if NormalChainSearch:
                            
                                Chain = re.compile("([a-zA-Z]+)").search(NewLowerLimit).group(1)
                                NewL = int(re.compile("([0-9]+)").search(NewLowerLimit).group(1))
                                NewU = int(re.compile("([0-9]+)").search(NewUpperLimit).group(1))
                                
                            else:
                                Chain = NewLowerLimit[0]
                                NewL = int(re.compile("([0-9]+)").search(NewLowerLimit[1:]).group(1))
                                NewU = int(re.compile("([0-9]+)").search(NewUpperLimit[1:]).group(1))
                            
                            NewP = (NewL + NewU) / 2
                            
                            NewAdjustPosition = Chain + str(NewP)
                            nNeedsReadjusts = len(Encapsulated_L)
                    
                    #if the readjust is completed, the original dict will be changed to reflect the readjusted positions
                    if ReadjustDone:
                        point = AdjustRange[0]
                        for EncapsulatedMutKey in Encapsulated_L:
                            
                            Chain = ""
                            NewL = None
                            NewU = None
                            NormalChainSearch = re.compile("([a-zA-Z]+)").search(PDBAccToMutKeysToPositions_D[AccKey][EncapsulatedMutKey])
                            
                            oldChain = None
                            if NormalChainSearch:
                                oldChain = re.compile("([a-zA-Z]+)").search(PDBAccToMutKeysToPositions_D[AccKey][EncapsulatedMutKey]).group(1)
                            else:
                                oldChain = PDBAccToMutKeysToPositions_D[AccKey][EncapsulatedMutKey][0]
                            PDBAccToMutKeysToPositions_D[AccKey][EncapsulatedMutKey] = oldChain + str(point)
                        
                            point +=1
        
        Ret_L = []
        elemCount = 0
        
        
        #print PDBAccToMutKeysToPositions_D
        
        #adds the list of mutation lines in the same order that it came in, but includes the readjusted positions
        for DesiredAccKey in PDBAccToMutKeysToPositions_D.keys():
            
            for mutLine in allMutationsWithPDBMatchedResidues:
                if re.compile("NOCOVERAGE").search(mutLine):
                    Ret_L.append(mutLine)
                else:
                    NewMutInfo = ""
                    CorrectAccKey = ""
                    B = False
                    
                    if MutKeys_L[elemCount] in PDBAccToMutKeysToPositions_D[DesiredAccKey].keys():
                        NewMutInfo = PDBAccToMutKeysToPositions_D[DesiredAccKey][MutKeys_L[elemCount]]
                        B = True
                
                    #print DesiredAccKey
                    #print NewMutInfo
                    #print "*" * 10
                    
                    NewMutLine = mutLine
                    
                    if B:
                        deletePattern = DesiredAccKey+ re.escape("|") + "[A-Za-z0-9]+"
                        replacePattern = DesiredAccKey+ "|" + PDBAccToMutKeysToPositions_D[DesiredAccKey][MutKeys_L[elemCount]]
                
                        NewMutLine = re.sub(deletePattern , replacePattern , mutLine)
                    
                    Ret_L.append(NewMutLine)
                                
                    elemCount += 1
            
            allMutationsWithPDBMatchedResidues = Ret_L
            elemCount = 0
            Ret_L = []
        
        #remove duplicates of the same accession position on the same line
        
        Ret_L = []
        for mutLine in allMutationsWithPDBMatchedResidues:
            newLine = mutLine
            
            if re.compile("NOCOVERAGE").search(newLine):
                pass
            else:
                allRes = re.compile("<R>(.+?)</R>").search(newLine).group(1).split(",")
                newAllRes = set(allRes)
                
                #print newLine
                #print len(allRes)
                #print len(newAllRes)
                
                if len(allRes) == len(newAllRes):
                    pass
                else:
                    newString = ",".join(list(newAllRes))
                    
                    newLine = re.sub("<R>.+?</R>","<R>%s</R>" % (newString) , newLine)
                    
            Ret_L.append(newLine)
            
            """
            for DesiredAccKey in PDBAccToMutKeysToPositions_D.keys():
                #print DesiredAccKey
                checkPattern = "("+DesiredAccKey+ re.escape("|") + "[A-Za-z0-9]+)"
                theCheckSearch = re.compile(checkPattern).search(newLine)
                if theCheckSearch:
                    checkString = theCheckSearch.group(1)
                    
                    checknumber = newLine.count(checkString)
                    
                    if checknumber > 1:
                        newLine = re.sub(checkString, "" , newLine , checknumber-1)
            """
            
                        
            
            
            
        
        
        return Ret_L
    
    "gets relevant inter-residue PDB keys for distance information"
    def completeDistanceInformation(self , allMutationsWithPDBMatchedResidues , scoring_matrix , pdb_file_hash):
        
        
        allMutationsWithDistanceInfo_L = []
        
        #starts with a list of mutations that have completed SAS and position information
        for mutation_a in allMutationsWithPDBMatchedResidues:
            
            #skips residues without coverage
            if re.compile("NOCOVERAGE").search(mutation_a):
                allMutationsWithDistanceInfo_L.append("<Mut>%s<D>NOCOVERAGE</D></Mut>" % (mutation_a))
            else:
                
                #prepares distance information to be placed in line
                mutationWithDistanceInfo_S = []
                mutationWithDistanceInfo_S.append(mutation_a)
                mutationWithDistanceInfo_S.append("<D>")
                
                #gets position of mutation a
                mutation_a_position = int(re.compile("<M>(.+?)</M>").search(mutation_a).group(1).split("|")[0])
                mutation_a_position_index = mutation_a_position - 1
                
                #for each pairwise combination with mutation a
                for mutation_b in allMutationsWithPDBMatchedResidues:
                    
                    #gets position of mutation b
                    mutation_b_position = int(re.compile("<M>(.+?)</M>").search(mutation_b).group(1).split("|")[0])
                    mutation_b_position_index = mutation_b_position - 1
                    
                    correct_pdb_accession_chain = scoring_matrix.Matrix[mutation_a_position_index][mutation_b_position_index]['pdb']
                    
                    #skips if there is no accession number at this point
                    if correct_pdb_accession_chain == 'none':
                        pass
                        
                    else:
                        
                        accession = correct_pdb_accession_chain.split("|")[0].lower()
                        chain = correct_pdb_accession_chain.split("|")[1]
                        
                        #gets the corresponding PDB positions for both mutations a and b
                        xmlseqhash = self.getStatsFromCorrectRefSeqXML(correct_pdb_accession_chain , scoring_matrix.ReferenceSequenceXML)
                        pdbposition_hash_a = self.getCorrectPositionAndResidueInPDBSeq(xmlseqhash , mutation_a_position)
                        pdbposition_hash_b = self.getCorrectPositionAndResidueInPDBSeq(xmlseqhash , mutation_b_position)
                        
                        pdb_position_a = chain + str(pdbposition_hash_a['pdb_correct_position'])
                        pdb_position_b = chain + str(pdbposition_hash_b['pdb_correct_position'])
                        #places the position of mutation a on the PDB structure and the corresponding position of mutation B on the structure
                        mutationWithDistanceInfo_S.append("<v>%s:%s->%s:%s;%s</v>" % (str(mutation_a_position) , pdb_position_a,\
                                                                                          str(mutation_b_position) , pdb_position_b,\
                                                                                          accession))
                mutationWithDistanceInfo_S.append("</D>")
                allMutationsWithDistanceInfo_L.append("<Mut>%s</Mut>" % ''.join(mutationWithDistanceInfo_S))
                    
        return allMutationsWithDistanceInfo_L
        
