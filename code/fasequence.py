####################################################################################################
#                                                                                                  #
#           PROJECT Protein Adaptation                                                             #
#           CLASS FASequence                                                                       #
#           PROGRAMMER Jeremy Adams                                                                #
#           STARTED 02-10-13                                                                       #
#           LASTMOD 12-08-14                                                                       #
#                                                                                                  #
#           DESCRIPTION Class representing a FASTA sequence that can be used as input for a BLASTP #
#                       search                                                                     #
#                                                                                                  #
####################################################################################################

import re
from staticmethods import *
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Application import ApplicationError

class FASequence:
    
    "CONSTRUCTOR"
    def __init__(self,header,sequence):
        
        """
        Class attributes:
        Info (Dict): dictionary containing header and sequence information
        PDBMasterList (List): List containing all blast hits that make a unique contribution to the coverage of the sequence by the PDB
        """
        
        self.info = {}
        self.info['header'] = header
        self.info['sequence'] = sequence
        self.info['query_coverage_threshold'] = 0.95
        self.PDBMasterList = []
        
    
    "returns a 1-dimensional list of dictionaries of equal length to the sequence proper. Each dictionary contains a very high evalue and a blank pdb accession number to be overwritten"
    def getSeqAsCellDictList(self):
        return [getSingleCell() for letter in self.info['sequence']]
    "appends a BLAST PDB record to the master list if it passes the criteria to be included"
    def addToMasterList(self,table_row):
        self.PDBMasterList.append(table_row)
    "returns a float representing the percentage of the query sequence that has been convered by one or more subject sequences"
    def getQueryCoverage(self,seq_array):
        qlen = 0.0
        cover = 0.0
        for position in seq_array:
            if position['pdb'] != 'none':
                cover += 1.0
            
            qlen +=1.0
        
        return float(cover / qlen)
    
    "checks if a new incoming evalue is better (lower) than the evalue currently sitting at that spot, if it is, then the new evalue and pdb number will take its spot"
    def fillOneSpaceInDictList(self,seq_array_index,qstart,qend,evalue,pdb,count,boolean):
        
        position = {}
        change_occurred = boolean
        
        if float(evalue) < 0.000001:
            if seq_array_index['evalue'] > float(evalue): #if the current cell's evalue is greater than the challenging evalue
                if count >= int(qstart) and count <= int(qend): #if the alignment hist covers over this part in the sequence then it is sufficient to change the cell value
                    position['evalue'] = float(evalue) #new position with new evalue
                    position['pdb'] = pdb #new position with new pdb accession
                    change_occurred = True #a change did occur
                else: #no change occurred, so it returns the exact same thing
                    position['evalue'] = seq_array_index['evalue']
                    position['pdb'] = seq_array_index['pdb']
            else: #no change occurred
                position['evalue'] = seq_array_index['evalue']
                position['pdb'] = seq_array_index['pdb']
        else: #no change occurred
            position['evalue'] = seq_array_index['evalue']
            position['pdb'] = seq_array_index['pdb']
                    
        return (change_occurred , position)
    
    "checks each dict list index for the length of the entire list using a single PDB hit"
    def fillUpDictList(self,seq_array,qstart,qend,evalue,sseqid,qlen):
                
        ret_list = []

        count = 1
        #starts off the first index and check if the new evalue will replace the old, at first a change has not occurred
        change_occurred , first_iteration = self.fillOneSpaceInDictList(seq_array[count-1],qstart,qend,evalue,sseqid,count,False)
        
        ret_list.append(first_iteration)
        
        count += 1
        #checks over each index in the entire dict list to change values, keeping "True" if even one cell has changed
        for seq_index in seq_array[1:]:
            change_occurred , iteration = self.fillOneSpaceInDictList(seq_array[count-1],qstart,qend,evalue,sseqid,count,change_occurred)
            ret_list.append(iteration)
            count += 1
            
        return (change_occurred , ret_list)
    
    "performs a BLASTP operation to the protein data bank and returns all hits and information contained therein"
    def blastPToPDB(self , dbPATH):
        
        retBool = True
        OFH = getOutputTempFile()
        try:
            #wraps the input sequence into an input temporary file, then executes the blastP
            IFH = getInputTempFile(self.info['sequence'])
            
            
            cline = NcbiblastpCommandline(query=IFH.name, db=dbPATH,\
                                          out=OFH.name,\
                                          outfmt='"10 qcovs evalue bitscore sseqid qlen qstart qend qseq slen sstart send sseq qcovhsp"')
            
            cline()
            IFH.close()
        except ApplicationError:
            retBool = False
            IFH.close()
            OFH.close()
        
        #returns the table and whether or not the blastP ran successfully
        ret = {'Success' : retBool}
        if retBool:
            ret['Table'] = OFH.read()
            OFH.close()
        return ret
    
    "adds a slash to a directory path if it is needed"
    def addSlashIfNeeded(self,PATH):
        Ret = PATH
        if PATH.endswith("/"):
            pass
        else:
            Ret = PATH+"/"
        return Ret
    
    
    "the main hierarchy process of getting pdb ids and adding all significant pdb records to the master list"
    def completeMasterList(self):
        
        #executes blast and gets a table of all the records
        blast_list = (self.blastPToPDB(self.addSlashIfNeeded(sys.argv[2])+'misc/pdbaa')['Table'].split('\n'))[:-1]
        #uses the first blast record to fill up the sequence dict list
        if len(blast_list) > 0:
            self.addToMasterList(blast_list[0])
        
            useless_bool , seq_array = self.fillUpDictList(self.getSeqAsCellDictList(),\
                                                              get_qstart(self.PDBMasterList[0]),\
                                                              get_qend(self.PDBMasterList[0]),\
                                                              get_evalue(self.PDBMasterList[0]),\
                                                              get_pdb_accession_number(get_sseqid(self.PDBMasterList[0])),\
                                                              get_qlen(self.PDBMasterList[0]))
            
            
            #if the sequence dict list is not sufficiently covered then it will complete the master list for all blast hit records
            if self.getQueryCoverage(seq_array) < self.info['query_coverage_threshold']:
                
                if len(blast_list) > 1:
                    
                    for blast_row in blast_list[1:]:
                        
                        change_occurred , seq_array = self.fillUpDictList(seq_array,\
                                                                                            get_qstart(blast_row),\
                                                                                            get_qend(blast_row),\
                                                                                            get_evalue(blast_row),\
                                                                                            get_pdb_accession_number(get_sseqid(blast_row)),\
                                                                                            get_qlen(blast_row))
                        if change_occurred:
                            self.addToMasterList(blast_row)
    
    "returns a single string representing fully formatted xml output of all blast pdb hits for this particular sequence"
    def referenceSequencePrint(self):
        
        qq = []
        
        qq.append("<Reference_sequence>\n")
            
        qq.append("\t<Header>"+self.info['header']+"</Header>\n")
        qq.append("\t<Sequence>"+self.info['sequence']+"</Sequence>\n")
            
        qq.append("\t<Sequence_length>"+str(len(self.info["sequence"]))+"</Sequence_length>\n")
        
        qq.append("\t<PDB_associations>\n")
        
        #relevant stats for each hit in the PDB master list
        for row in self.PDBMasterList:
            qq.append("\t\t<PDB_alignment>\n")
            
            qq.append("\t\t\t<PDB_id>"+get_pdb_accession_number(get_sseqid(row))+"|"+get_pdb_chain(get_sseqid(row))+"</PDB_id>\n")
            qq.append("\t\t\t<E_value>"+get_evalue(row)+"</E_value>\n")
            qq.append("\t\t\t<Alignment_start_query>"+get_qstart(row)+"</Alignment_start_query>\n")
            qq.append("\t\t\t<Alignment_end_query>"+get_qend(row)+"</Alignment_end_query>\n")
            qq.append("\t\t\t<Aligned_query_sequence>"+get_qseq(row)+"</Aligned_query_sequence>\n")
            qq.append("\t\t\t<Alignment_start_subject>"+get_sstart(row)+"</Alignment_start_subject>\n")
            qq.append("\t\t\t<Alignment_end_subject>"+get_send(row)+"</Alignment_end_subject>\n")
            qq.append("\t\t\t<Aligned_subject_sequence>"+get_sseq(row)+"</Aligned_subject_sequence>\n")
            qq.append("\t\t</PDB_alignment>\n")
            
        qq.append("\t</PDB_associations>\n")
    
        qq.append("</Reference_sequence>")
        
        return ''.join(qq)
    