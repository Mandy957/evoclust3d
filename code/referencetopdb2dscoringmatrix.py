####################################################################################################
#                                                                                                  #
#           PROJECT Protein Adaptation                                                             #
#           CLASS ReferenceToPDB2DScoringMatrix                                                    #
#           PROGRAMMER Jeremy Adams                                                                #
#           STARTED 25-07-13                                                                       #
#           LASTMOD 06-11-14                                                                       #
#                                                                                                  #
#           DESCRIPTION class representation of the two-dimensional grid created by the sequence's #
#                       relation to itself and the PDB structure mappings contained therein        #
#                                                                                                  #
####################################################################################################

import re
import tempfile
from staticmethods import *

class ReferenceToPDB2DScoringMatrix:
    
    "CONSTRUCTOR"
    def __init__(self,reference_sequence_xml):
        
        """
        Class attributes:
        ReferenceSequenceXML (String): String containing information of what PDB structures the reference sequence BLASTed to and where
        SuccessfullyAligned (Bool): False if the sequence did not blast to any PDB entries, True if it did
        ProcessedXML (Dict): Contains key length and PDB hit information about the reference sequence xml
        PDBCoverageKeys (Dict): Key is the PDB ID, Value is a list of integers and chain IDs stating where on the structure the sequence validly BLASTed to
        Matrix (List of Lists): 2-D matrix where each individual square contains the best (lowest e-value) PDB structure containing those two coordinate points
        """
        
        #takes input from argument
        self.ReferenceSequenceXML = reference_sequence_xml
        self.SuccessfullyAligned = True
        
        #if reference_sequence_xml is blank, the sequence did not BLAST to anything and hence it was not successfully aligned
        if self.ReferenceSequenceXML == "":
            self.SuccessfullyAligned = False
            self.ProcessedXML = self.process_reference_sequence_xml(self.ReferenceSequenceXML)
            self.Matrix = []
            self.PDBCoverageKeys = {}
        #if there were BLAST hits, then the sequence was successfully aligned
        else:
            
            self.ProcessedXML = self.processReferenceSequenceXML(self.ReferenceSequenceXML)
            self.Matrix = self.computeMatrix()
            self.PDBCoverageKeys = self.getPDBCoverageKeys()
    
    "creates the 2-D matrix containing PDB entries in each cell"
    def computeMatrix(self):
        
        #gets blank 2-D square grid with length the same as the sequence, cells are initialized to an e-value of infinity
        matrix = []
        matrix = self.getAllRows()
        
        #for each PDB that was hit
        for pdb in self.ProcessedXML['pdb_files']:
            #challenger stats are those of the current PDB structure
            challenger_evalue = float(pdb['e_value'])
            challenger_pdb = pdb['pdb_id']
            challenger_start = int(pdb['start'])
            challenger_end = int(pdb['end'])
            
            #check each cell in the matrix
            row_count = 1
            for row in matrix:
                cell_count = 1
                
                #checks if the PDB structure alignment falls within the bounds of that cell in the grid
                if row_count >= challenger_start and row_count <= challenger_end:
                    for cell in row:
                        if cell_count >= challenger_start and cell_count <= challenger_end:
                            
                            #checks if the challenger e-value is better than what is already there
                            if challenger_evalue < matrix[row_count-1][cell_count-1]['evalue']:
                                #if the challenger evalue is lower than what was there, than it is replace by the challenger evalue and PDB ID
                                matrix[row_count-1][cell_count-1]['evalue'] = challenger_evalue
                                matrix[row_count-1][cell_count-1]['pdb'] = challenger_pdb
                                
                    
                        cell_count += 1
                    
                row_count += 1
              
        return matrix
    
    "gets a single row of single blank cells the length of the sequence"
    def getSingleRow(self):
        return [getSingleCell() for letter in range(0,int(self.ProcessedXML['length']))]
    "gets all rows for the length of the sequence"
    def getAllRows(self):
        return [self.getSingleRow() for letter in range(0,int(self.ProcessedXML['length']))]
    "processes part of the reference sequence xml for a single PDB ID, giving its start and end coordinates, and associated evalue"
    def processSinglePDBAlignment(self,xml):
        return {'pdb_id' : re.compile("<PDB_id>(.+?)</PDB_id>").search(xml).group(1),\
                'e_value' : re.compile("<E_value>(.+?)</E_value>").search(xml).group(1),\
                'start' : re.compile("<Alignment_start_query>(.+?)</Alignment_start_query>").search(xml).group(1),\
                'end' : re.compile("<Alignment_end_query>(.+?)</Alignment_end_query>").search(xml).group(1)}
    
    "processes each PDB ID into a dictionary for all PDB IDs that were hit"
    def processAllPDBAlignments(self,xml):
        xml_list = re.findall("<PDB_alignment>.+?</PDB_alignment>" , xml , re.DOTALL)
        return [self.processSinglePDBAlignment(xml_subset) for xml_subset in xml_list]
        
    "gets the length of the sequence and all of the PDB information from the reference sequence xml"
    def processReferenceSequenceXML(self,xml):
        
        length = re.compile("<Sequence_length>(.+?)</Sequence_length>").search(xml).group(1)
        return {'length' : length,\
                'pdb_files' : self.processAllPDBAlignments(xml)}
    
    "gets a dictionary of keys to later be used in the random distribution processing portion of the program"
    def getPDBCoverageKeys(self):
        PDBCoverageKeys_D = {}
        
        #one key in the dictionary for each PDB ID
        for PDBAlignment in re.findall("<PDB_alignment>.+?</PDB_alignment>",self.ReferenceSequenceXML,re.DOTALL):
            #gets the accession and key
            AccChain = re.compile("<PDB_id>(.+?)</PDB_id>").search(PDBAlignment).group(1).split("|")
            
            if AccChain[0].lower() in PDBCoverageKeys_D.keys():
                pass
            else:
                PDBCoverageKeys_D[AccChain[0].lower()] = []
            
            #gets the start and end alignment coordinates for the PDB hit
            PDBStart = int(re.compile("<Alignment_start_subject>(.+?)</Alignment_start_subject>").search(PDBAlignment).group(1))
            PDBEnd = int(re.compile("<Alignment_end_subject>(.+?)</Alignment_end_subject>").search(PDBAlignment).group(1))
            
            #adds one key for each index between the start and end coordinates inclusive
            for i in range(PDBStart,PDBEnd+1):
                PDBCoverageKeys_D[AccChain[0].lower()].append(AccChain[1]+str(i))
        
        return PDBCoverageKeys_D
    
    "prints out the coverage keys to a file"
    def scoringMatrixXMLPrint(self):
        
        qq = ["<ScoringMatrix>\n"]
        
        qq.append("\t"+re.compile("(<PDB_associations>.+?</PDB_associations>)" , re.DOTALL).search(self.ReferenceSequenceXML).group(1))
        qq.append("\n")
        
        """
        qq.append("\t<CellHits>\n")
        for r in range(0,len(self.Matrix)):
            for c in range(0,len(self.Matrix[r])):
                
                if self.Matrix[r][c]["pdb"] != "none":
                    qq.append("<H>%s,%s,%s</H>\n" % (str(r+1) , str(c+1) , self.Matrix[r][c]["pdb"]))
            
        qq.append("\t</CellHits>\n")
        """
        
        #start xml tag
        qq.append("\t<Coverages>\n")
        
        #for each PDB ID that was hit, add the PDB ID key
        for PDBKey in self.PDBCoverageKeys.keys():
            t = []
            qq.append("\t\t<Coverage><ID>%s</ID><Keys>%s</Keys></Coverage>\n" % (PDBKey , ",".join([CoverageKey for CoverageKey in self.PDBCoverageKeys[PDBKey]])))
            
        #end xml tag
        qq.append("\t</Coverages>\n")
        
        qq.append("</ScoringMatrix>")
        return ''.join(qq)
    