####################################################################################################
#                                                                                                  #
#           PROJECT Protein Adaptation                                                             #
#           CLASS WholeTreeGeneScope                                                               #
#           PROGRAMMER Jeremy Adams                                                                #
#           STARTED 03-12-13                                                                       #
#           LASTMOD 07-08-14                                                                       #
#                                                                                                  #
#           DESCRIPTION Top level class for the performance of the protein adaptation detection    #
#                       method/program                                                             #
#                                                                                                  #
####################################################################################################

from staticmethods import *
from fastmltree import FastMLTree
from wholetreeorthologoussubgroup import WholeTreeOrthologousSubgroup

import re
import sets
import os


class WholeTreeGeneScope:
    
    """
    Class attributes:
    FastaPATH (String): Path to user-defined sequence file in FASTA format
    UserTreePATH (String): Path to user-defined tree file in newick format
    ProjectName (String): user-defined project name to determine output file names
    
    ExitStatus (Bool): True if the program can execute fully, False if there is an error
    ExitString (String): The final XML format output string, or error string if there is an error
    OSGString (String): String returned by the orthologous subgroup that contains most mutation analysis information
    
    FastMLTree (FastMLTree class object): Tree object that parses newick tree and prepares it for mutation analysis in evolutionary stages
    FastaKey_L (List): List of strings representing the headers of sequences
    Fasta_D (Dict): Key is the Fasta sequence header, value is the sequence itself
    
    OSG (WholeTreeOrthologousSubgroup class object): Object that performs the bulk of the protein adaptation detection analysis
    """
    
    "CONSTRUCTOR"
    def __init__(self , FastaPATH , UserTreePATH , ProjectName):
        self.FastaPATH = FastaPATH
        self.UserTreePATH = UserTreePATH
        self.ProjectName = ProjectName
        
        #declaration of output variables
        self.ExitStatus = False
        self.ExitString = ""
        self.OSGString = ""
        
        #parses the tree file according to FastMLTree methods
        self.FastMLTree = FastMLTree(self.UserTreePATH , True)
        
        if self.FastMLTree.Parsed:
            #opens the Fasta file and gets the sequences as a dictionary
            ReadFasta = readFasta(self.FastaPATH)
            self.FastaKey_L = ReadFasta[0]
            self.Fasta_D = ReadFasta[1]
            
            #validates the sequences with the tree
            ValidSeqsWithTree = self.ValidateSeqsWithTree()
            
            #if sequence headers match with node headers
            if ValidSeqsWithTree[0]:
                
                #instantiates the WholeTreeOrthologousSubgroup class object
                self.OSG = WholeTreeOrthologousSubgroup(self.FastMLTree , self.FastaPATH , self.FastaKey_L , self.Fasta_D)
                
                #output string from the WholeTreeOrthologousSubgroup
                self.ExitStatus = True
                #print "A"
                #if the analysis worked
                if self.ExitStatus:
                    #print "B"
                    #prepares the final XML output string
                    exit_string = []
                    exit_string.append("<Group>\n")
                    exit_string.append("\t<Group_id>NA</Group_id>\n\t<Number_OSGs>1</Number_OSGs>\n")
                    exit_string.append("\t<OSGs>\n")
                    exit_string.append(self.OSG.ExitString)
                    exit_string.append("\t</OSGs>\n")
                    exit_string.append("</Group>\n")
                    self.ExitString = ''.join(exit_string)
                    
                    os.system("mkdir -p %s" % self.ProjectName)
                    
                    #writes the protein adaptation XML file
                    with open("%s/Report.xml" % (self.ProjectName) , "w") as w:
                        w.write(self.ExitString)
                    
                    #writes a new modified newick tree file with internal node names according to the cogent convention
                    with open("%s/ModdedTree.nwk" % (self.ProjectName) , "w") as w:
                        w.write(self.FastMLTree.CogentTree.getNewick(with_distances=True).replace("'",""))
                    
                    #writes an XML file containing information pertaining to the ReferenceToPDB2DScoringMatrix object
                    with open("%s/ScoringMatrix.xml" % (self.ProjectName) , "w") as w:
                        w.write(self.OSG.MatrixGraphicsString)
                        
                    #writes FASTA file of reconstructed sequences
                    with open("%s/AncestralSeqs.fa" % (self.ProjectName) , "w") as w:
                        w.write(self.OSG.ReconstructedFASTAString)
                    
                    #writes text file of reconstruction probabilities
                    with open("%s/AncestralProb.txt" % (self.ProjectName) , "w") as w:
                        w.write(self.OSG.ReconstructedProbabilityString)
                    
                    #Final output message
                    print "Done.\n" + \
                          "Main output XML file written to %s/Report.xml\n" % (self.ProjectName) + \
                          "Modified tree Newick file written to %s/ModdedTree.nwk\n" % (self.ProjectName) + \
                          "Scoring Matrix XML file written to %s/ScoringMatrix.xml\n"  % (self.ProjectName) + \
                          "FASTA format file of reconstructed sequences written to %s/AncestralSeqs.fa\n" % (self.ProjectName) + \
                          "Text file of reconstructed sequence probabilities written to %s/AncestralProb.txt" % (self.ProjectName)
                
            #if the tree terminal nodes were not matched with the Fasta sequence headers 
            else:
                self.ExitString = self.getTreeMatchedToFastaErrorMessage()
                print self.ExitString
        
        #if the FastMLTree object was not properly parsed
        else:
            self.ExitString = self.getTreeErrorMessage()
            print self.ExitString
        
    "Method to validate if the input data is compatible to perform the full protein adaptation analysis"
    def ValidateSeqsWithTree(self):
        Ret = False
        
        InFaButNotInTree = []
        InTreeButNotInFa = []
        
        #print "-".join(sorted(self.FastaKey_L))
        #print "*" * 50
        #print "-".join(sorted(self.FastMLTree.LeafKey_L))
        
        if "-".join(sorted(self.FastaKey_L)) == "-".join(sorted(self.FastMLTree.LeafKey_L)):
            Ret = True
        else:
            
            for FaKey in self.FastaKey_L:
                if FaKey in set(self.FastMLTree.LeafKey_L):
                    pass
                else:
                    InFaButNotInTree.append(FaKey)
            
            for TreeKey in self.FastMLTree.LeafKey_L:
                if TreeKey in set(self.FastaKey_L):
                    pass
                else:
                    InTreeButNotInFa.append(TreeKey)
        
        #print "\n".join(L)
            
        return [Ret , InFaButNotInTree , InTreeButNotInFa]
    
    "Returns an error message explaining something went wrong with parsing the phylogenetic tree"
    def getTreeErrorMessage(self):
        return "ERROR: Problem parsing the phylogenetic tree.\n" +\
               "Please make sure that the tree file is in valid Newick syntax, and that branch lengths are included."
    
    "Returns an error message explaining that the terminal node names in the tree are not matched to the headers in the Fasta file"
    def getTreeMatchedToFastaErrorMessage(self):
        return "ERROR: Terminal node names in the tree file do not match sequence headers in the FASTA file.\n" +\
               "Please make sure that there is an exact correspondence between terminal nodes in the tree file and sequences in the FASTA file " +\
               "and that the names in each match exactly."
    

    
    