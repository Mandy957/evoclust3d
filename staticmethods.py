####################################################################################################
#                                                                                                  #
#           PROJECT Protein Adaptation                                                             #
#           FILE staticmethods.py                                                                  #
#           PROGRAMMER Jeremy Adams                                                                #
#           STARTED 13-09-13                                                                       #
#           LASTMOD 12-08-14                                                                       #
#                                                                                                  #
#           DESCRIPTION File containing non-class specific methods, including methods that are     #
#                       used in a wide array of general instances                                  #
#                                                                                                  #
####################################################################################################

import re
import tempfile
import time
import os
import string
import random
import sets
import scipy.stats
import pexpect
import sys
import numpy
from numpy import array
from pdbatom import PDBAtom
import re
import math
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO


##################################################RANDOM METHODS####################################

"adds a slash to a directory path if it is needed"
def addSlashIfNeeded(PATH):
    Ret = PATH
    if PATH.endswith("/"):
        pass
    else:
        Ret = PATH+"/"
    return Ret

"used to generate a random string of characters to create temporary directories"
def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

##################################################FASTA METHODS#####################################

"takes in a path to a Fasta file and parses out all of the information"
def readFasta(PATH):
    F = open(PATH,"r")
    AllLines_L = F.readlines()
    F.close()
    
    ret_D = {}  #dictionary where key is the header, value is the sequence
    retKey_L = []   #list with all of the dictionary keys in the proper order
    
    AllLinesInOneSequence = [AllLines_L[0].replace("\n","")]
    
    for line in AllLines_L[1:]:
        lineR = line.replace("\n","")
        
        if lineR.startswith(">"):   #previous sequence is done, new sequence is beginning
            ret_D[AllLinesInOneSequence[0][1:]] = "".join(AllLinesInOneSequence[1:])
            retKey_L.append(AllLinesInOneSequence[0][1:])
            
            AllLinesInOneSequence = [lineR]
            
        else:   #this line is part of the same sequence as the previous line
            AllLinesInOneSequence.append(lineR)
    
    ret_D[AllLinesInOneSequence[0][1:]] = "".join(AllLinesInOneSequence[1:])
    retKey_L.append(AllLinesInOneSequence[0][1:])
    
    return [retKey_L , ret_D]

##################################################BLAST RESULTS METHODS#############################

"returns regex pattern to retrieve all information from BLAST PDB result row"
def master_list_regex_pattern():
    return re.compile("(.+?),(.+?),(.+?),(.+?),(.+?),(.+?),(.+?),(.+?),(.+?),(.+?),(.+?),(.+?),(.+?)")
"gets query coverage of the blast hit"
def get_qcovs(table_row):
    return master_list_regex_pattern().search(table_row).group(1)
"gets blast hit expect value"
def get_evalue(table_row):
    return master_list_regex_pattern().search(table_row).group(2)
"gets blast hit bitscore"""
def get_bitscore(table_row):
    return master_list_regex_pattern().search(table_row).group(3)
"gets blast hit subject sequence identifier"""
def get_sseqid(table_row):
    return master_list_regex_pattern().search(table_row).group(4)
"gets blast hit query sequence length"""
def get_qlen(table_row):
    return master_list_regex_pattern().search(table_row).group(5)
"gets blast hit position of query sequence where alignment begins"""
def get_qstart(table_row):
    return master_list_regex_pattern().search(table_row).group(6)
"gets blast hit position of query sequence where alignment ends"
def get_qend(table_row):
    return master_list_regex_pattern().search(table_row).group(7)
"gets blast hit sequence of query that is aligned to the subject"
def get_qseq(table_row):
    return master_list_regex_pattern().search(table_row).group(8)
"gets blast hit subject sequence length"
def get_slen(table_row):
    return master_list_regex_pattern().search(table_row).group(9)
"gets blast hit position of subject sequence where alignment begins"
def get_sstart(table_row):
    return master_list_regex_pattern().search(table_row).group(10)
"gets blast hit position of subject sequence where alignemnt ends"
def get_send(table_row):
    return master_list_regex_pattern().search(table_row).group(11)
"gets blast hit sequence of subject that is aligned to the query"""
def get_sseq(table_row):
    return master_list_regex_pattern().search(table_row).group(12)
"gets blast hit query coverage high score predict"
def get_qcovhsp(table_row):
    return master_list_regex_pattern().search(table_row).group(13)
"gets blast hit pdb accession number"""
def get_pdb_accession_number(sseqid):
    
    p = re.escape("|")
    
    pattern = re.compile('.+?'+p+".+?"+p+".+?"+p+"(.+?)"+p+".+?")
     
    return pattern.search(sseqid).group(1)
"gets blast hit pdb chain"
def get_pdb_chain(sseqid):
    p = re.escape("|")
    
    pattern = re.compile('.+?'+p+".+?"+p+".+?"+p+".+?"+p+"(.+?)")
    
    return pattern.search(sseqid).group(1)
    
##################################################SCORING MATRIX METHODS############################

"gets a single small dict representing one cell in the 2-D scoring matrix"
def getSingleCell():
        return {'evalue'  :   100000,\
                'pdb'  :   'none'}

##################################################TEMPORARY FILE METHODS############################

"gets a temporary file with the input content written to it"
def getInputTempFile(Content):
    handle = tempfile.NamedTemporaryFile(mode='r+',suffix='',prefix='', delete=True)
    handle.write(Content)
    handle.seek(0)
    return handle
"gets a temporary file with no content written to it so that an external program can write to this file"
def getOutputTempFile():
    handle = tempfile.NamedTemporaryFile(mode='r+',suffix='',prefix='', delete=True)
    return handle

##################################################FASTML METHODS####################################

"executes the FastML program using the specified alignment file and tree file specified"
def executeFastML(AlignmentContents , TreeString , NodesOrSequences):
    
    randDIR = "/tmp/%s/" % (id_generator()) #generates a random directory that FastML will write all files into
    #temporary files for FASTA and tree contents
    AlignmentFH = getInputTempFile(AlignmentContents)
    TreeFH = getInputTempFile(TreeString)
    #executes FastML
    os.system("perl %sFastML/www/fastml/FastML_Wrapper.pl --MSA_File %s --seqType aa --outDir %s --Tree %s --jointReconstruction no > /dev/null" %  (addSlashIfNeeded(sys.argv[2]),\
                                                                                                                                                     AlignmentFH.name,\
                                                                                                                                                     randDIR,\
                                                                                                                                                     TreeFH.name))
    
    ret = {}
    
    if NodesOrSequences: #gets the Newick format file with nodes named according to FastML convention
        ret = open("%stree.newick.txt" % (randDIR) , "r").read()
    else:   #gets the contents of the reconstructed sequence FASTA file and the probability information
        
        OutputProbLines_L = open("%sprob.marginal.csv" % (randDIR) , "r").readlines()
        ret["OutputProbLines_L"] = OutputProbLines_L 
        
        if os.path.exists("%sseq.marginal_IndelAndChars.txt" % (randDIR)):    
            OutputFASTA = open("%sseq.marginal_IndelAndChars.txt" % (randDIR) , "r").read()
            ret["OutputFASTA"] = OutputFASTA
        else:
            ret["OutputFASTA"] = "NOPATH"
               
    os.system("rm -rf %s" % (randDIR)) #removes all the files that were just generated
    
    return ret

##################################################TREE METHODS######################################

"Renames internal nodes for later parsing by pyCogent"
def fixUpFileForCogent(tree_file):
    
    #reads in tree string
    TreeString = open(tree_file , "r").read().replace("\n","")
    
    Letters_L = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]
    LetterACount = 0
    LetterBCount = 0
    i = 0
    NotThroughTheString = True
    
    #while loop moves one space along the tree string until it reaches the end
    while NotThroughTheString:
        #close bracket symbolizes end of clade and hence is the spot to name a node
        if TreeString[i] == ")":
            
            #new node name is created
            newInternalNodeName = Letters_L[LetterACount]+Letters_L[LetterBCount]
            lengthToColon = len(re.compile("^(.+?)[:;]").search(TreeString[i:]).group(1)) - 1
            
            #replaces the original node name (if any) with the replace name
            TreeString = TreeString[:i+1]+ newInternalNodeName + TreeString[i+lengthToColon+1:]
            
            #changes what the new internal node name will be
            if LetterBCount == 25:
                LetterBCount = 0
                LetterACount += 1
            else:
                LetterBCount += 1
                
        #check to end while loop
        if i == len(TreeString) - 1:
            NotThroughTheString = False
        i += 1
    
    TreeString = TreeString.replace("-","")
    
    #replaces actual distances with 0.0 distance for convenience
    if re.compile("\).+?:[0-9.]+?;$").search(TreeString):
        pass
    else:
        TreeString = TreeString[:-1] + ":0.0" + TreeString[-1]
    
    #returns temp file containing this modified string
    modifiedFile = getInputTempFile(TreeString)
    
    return [modifiedFile,TreeString]

"Fully parses cogent tree object to give all branch segments and nodes"
def completeNodesLeavesBranches(CogentTreeObj):
    
    NodeKey_L = CogentTreeObj.getNodeNames() #all node names
    LeafKey_L = [str(tip).split(":")[0].replace("'","") for tip in CogentTreeObj.iterTips()] #only the tips
    UpperKey_L = []
    
    #all the nodes that are not tips are upper keys
    for NodeKey in NodeKey_L:
        if NodeKey in sets.Set(LeafKey_L):
            pass
        else:
            UpperKey_L.append(NodeKey)
    
    TopKey = ""
    Nodes_D = {Key : {'immediate' : [] , 'terminal' : []} for Key in NodeKey_L} #each node is a sub dict containing lists of immediate derived and termini under that node
    
    for NodeKey in NodeKey_L:
        Node = CogentTreeObj.getNodeMatchingName(NodeKey)
        #if the node has no ancestors, it is the root or topkey
        if len(Node.ancestors()) == 0:
            TopKey = NodeKey
        else:
            #each node that is not the root has 1 ancestor, the ancestors immediate derived is the current node
            ImmediateAncestorName = Node.ancestors()[0].Name
            Nodes_D[ImmediateAncestorName]['immediate'].append(NodeKey)
    
    #for each node, get the terminial nodes underneath it
    for NodeKey in NodeKey_L:
        Node = CogentTreeObj.getNodeMatchingName(NodeKey)
        
        for terminalChild in Node.iterTips():
            Nodes_D[NodeKey]['terminal'].append(terminalChild.Name)
        
    BranchKey_L = []
    #all branches are all collections of immediate ancestor to immediate derived relationships
    for AncestorNodeKey in Nodes_D.keys():
        for ChildNode in Nodes_D[AncestorNodeKey]['immediate']:
            BranchKey_L.append(AncestorNodeKey+">>"+ChildNode)
    #returns entire dict
    return {"NodeKey_L" : NodeKey_L,\
            "LeafKey_L" : LeafKey_L,\
            "UpperKey_L" : UpperKey_L,\
            "TopKey" : TopKey,\
            "BranchKey_L" : BranchKey_L,\
            "Nodes_D" : Nodes_D}

##################################################PDB Dictionary METHODS############################

"gets atom key list and atom specification dictionary for a list of PDB file atom lines"
def get_pdbatom_H(atom_array):
    atom_D = {}
    atomkey_L = []
    
    for atom in atom_array:
        #sets the key
        s = ':%s.%s-%s' % (atom.Info['residue_num'],\
                           atom.Info['chain'],\
                           atom.Info['name'])
        
        atom_D[s] = atom #adds the atom into the dictionary
        atomkey_L.append(s) #adds the key to the list
    
    return {'atom_D' : atom_D,\
            'atomkey_L' : atomkey_L}

"gets residue key list and residue specification dictionary for a list of PDB file residues"
def get_pdbresidue_H(atom_D , atomkey_L):
    residue_D = {}
    residuekey_L = []
    
    residue_check_pattern = re.compile("^(:.+?\..+?)-")
    
    atoms_per_residue = []
    test_residue_check = residue_check_pattern.search(atomkey_L[0]).group(1)
    
    #checks each atom
    for atomkey in atomkey_L:
        current_residue_check = residue_check_pattern.search(atomkey).group(1)
        
        #if the atom belongs to the same residue as the last atom, they go in the same list
        if current_residue_check == test_residue_check:
            atoms_per_residue.append(atom_D[atomkey])
        
        #otherwise, they go into separate lists and the former list gets added to the dictionary
        else:
            residuekey_L.append(test_residue_check)
            residue_D[test_residue_check] = atoms_per_residue
            
            atoms_per_residue = []
            atoms_per_residue.append(atom_D[atomkey])
            test_residue_check = current_residue_check
    
    residuekey_L.append(test_residue_check)
    residue_D[test_residue_check] = atoms_per_residue
    
    return {'residue_D' : residue_D,\
            'residuekey_L' : residuekey_L}

"gets residue key list and residue specification dictionary for all lines in a PDBXML file"
def get_pdbXML_residue_H(residue_L):
    XMLresidue_D = {}
    XMLresiduekey_L = []
    
    #for each residue, the key is the chain+position combination
    for residue in residue_L:
        key = re.compile("<n>(.+?)</n>").search(residue).group(1)
        XMLresidue_D[key] = residue
        XMLresiduekey_L.append(key)
    
    return {'XMLresidue_D' : XMLresidue_D,\
            'XMLresiduekey_L' : XMLresiduekey_L}


"gets PDBContents dictionary and PDBXMLContents dictionary"
def getAllPDBFileDicts(PDB_L):
    AccToPDBContents_D = {}
    AccToPDBXMLContents_D = {}
    
    
    for AccKey in PDB_L:
                 
        #sets dictionary spaces for the new PDB ID key
        AccToPDBContents_D[AccKey] = {}
        AccToPDBXMLContents_D[AccKey] = {}
        
        #checks if the paths to PDB file and PDBXML file exist
        pdbIDHasFile = True
        
        if os.path.exists("%s%s/pdb%s.ent" % (addSlashIfNeeded(sys.argv[2])+"pdb/" , AccKey[1:3] , AccKey)):
            pass
        else:
            pdbIDHasFile = False
        if os.path.exists("%s%s/pdb%s.xml" % (addSlashIfNeeded(sys.argv[2])+"pdbxml/" , AccKey[1:3] , AccKey)):
            pass
        else:
            pdbIDHasFile = False
    
        if pdbIDHasFile:
            pdbF = open("%s%s/pdb%s.ent" % (addSlashIfNeeded(sys.argv[2])+"pdb/" , AccKey[1:3] , AccKey) , "r")
            pdbContent = pdbF.readlines()
            pdbF.close()
            #gets a list of atom lines for each line in the PDB file
            AccToPDBContents_D[AccKey]['Atoms_L'] = [PDBAtom(line) for line in pdbContent if line.startswith("ATOM")]
            
            #gets atom key list and atom dictionary for each atom in the list
            pdbatoms_H = get_pdbatom_H(AccToPDBContents_D[AccKey]['Atoms_L'])
            AccToPDBContents_D[AccKey]['AtomKey_L'] = pdbatoms_H['atomkey_L']
            AccToPDBContents_D[AccKey]['Atoms_D'] = pdbatoms_H['atom_D']
            
            #gets residue key list and residue dictionary for each residue in the PDB file
            pdbresidues_H = get_pdbresidue_H(AccToPDBContents_D[AccKey]['Atoms_D'],\
                                                 AccToPDBContents_D[AccKey]['AtomKey_L'])
            AccToPDBContents_D[AccKey]['ResidueKey_L'] = pdbresidues_H['residue_D']
            AccToPDBContents_D[AccKey]['Residue_D'] = pdbresidues_H['residuekey_L']
            
            #gets a list of residue lines for each line in the PDBXML file
            pdbxmlF = open("%s%s/pdb%s.xml" % (addSlashIfNeeded(sys.argv[2])+"pdbxml/", AccKey[1:3], AccKey) , "r")
            pdbxmlContent = pdbxmlF.read()
            pdbxmlF.close()
            
            AccToPDBXMLContents_D[AccKey]['Residues_L'] = [line for line in re.compile("<Residues>\n(.+?)\n\t</Residues>" , re.DOTALL).search(pdbxmlContent).group(1).split("\n")]
            
            #gets residue key list and residue dictionary for each residue in the PDBXML file
            pdbXMLresidues_H = get_pdbXML_residue_H(AccToPDBXMLContents_D[AccKey]['Residues_L'])
            
            AccToPDBXMLContents_D[AccKey]['XMLResidueKey_L'] = pdbXMLresidues_H['XMLresiduekey_L']
            AccToPDBXMLContents_D[AccKey]['XMLResidue_D'] = pdbXMLresidues_H['XMLresidue_D']
                        
    return [AccToPDBContents_D , AccToPDBXMLContents_D]

##################################################SIGNIFICANCE METHODS##############################

"method to combine a list of p-values into a single p-value"
def getCombinedPValue(PValue_L):
    
    V = []
    
    for PVal in PValue_L:
        if PVal == 0.0:
            V.append(float("0.001"))
        else:
            V.append(float(PVal))
    
    S = sum([math.log(i) for i in V])
    ChiSquare = -2.0 * S
    L = 2*len(PValue_L)
    P = 1.0 - scipy.stats.chi2.cdf(ChiSquare,L)
    
    return P

##################################################MERCATOR METHODS##################################

"reverse complements a string of DNA"
def reverseComplement(DNASeq):
    D = {"A":"T","T":"A","C":"G","G":"C"}
    return "".join([D[i] for i in DNASeq[::-1]])

"converts DNA (with T) to an mRNA (with U)"
def convertTomRNA(DNASeq):
    return DNASeq.replace("T","U")

"uses the genetic code to translate an mRNA sequence to an amino acid sequence"
def translate(mRNASeq):
    D = {"UUU":"F","UCU":"S","UAU":"Y","UGU":"C",\
         "UUC":"F","UCC":"S","UAC":"Y","UGC":"C",\
         "UUA":"L","UCA":"S","UAA":"*","UGA":"*",\
         "UUG":"L","UCG":"S","UAG":"*","UGG":"W",\
         "CUU":"L","CCU":"P","CAU":"H","CGU":"R",\
         "CUC":"L","CCC":"P","CAC":"H","CGC":"R",\
         "CUA":"L","CCA":"P","CAA":"Q","CGA":"R",\
         "CUG":"L","CCG":"P","CAG":"Q","CGG":"R",\
         "AUU":"I","ACU":"T","AAU":"N","AGU":"S",\
         "AUC":"I","ACC":"T","AAC":"N","AGC":"S",\
         "AUA":"I","ACA":"T","AAA":"K","AGA":"R",\
         "AUG":"M","ACG":"T","AAG":"K","AGG":"R",\
         "GUU":"V","GCU":"A","GAU":"D","GGU":"G",\
         "GUC":"V","GCC":"A","GAC":"D","GGC":"G",\
         "GUA":"V","GCA":"A","GAA":"E","GGA":"G",\
         "GUG":"V","GCG":"A","GAG":"E","GGG":"G"}
    
    L = int(len(mRNASeq)/3) * 3
    
    return "".join([D[mRNASeq[i:i+3]] for i in range(0,L,3)])

"sort by numerical value"
def numeric_compare(x,y):
    return x-y

"executes clustalw to perform a multiple sequence alignment"
def executeClustalW(inputString):
    inFH = getInputTempFile(inputString)
    outFH = getOutputTempFile()
    
    clustalw_cline = ClustalwCommandline(infile=inFH.name , outfile=outFH.name)
    clustalw_cline()
    
    output = outFH.read()
    return output

"parses clustalw results"
def parseClustalW(clustalOutput):
    
    D = {}
    
    for line in clustalOutput.split("\n"):
        if len(line) == 0:
            pass
        elif line.startswith("CLUSTAL"):
            pass
        elif re.compile("[\.:\*]").search(line):
            pass
        elif re.compile("^\s+$").search(line):
            pass
        else:
            
            ls = line.split()
            if ls[0] in D.keys():
                pass
            else:
                D[ls[0]] = []
            
            D[ls[0]].append(ls[1])
    
    R = {Key : "".join(D[Key]) for Key in D.keys()}
    return R

"executes phyml in Mercator to get the branch lengths between orthologous protein sequences"
def usePhyMLForBranchLengths(alignmentFasta,newicktree):
    #converts the alignment file in FASTA to PHYLIP format
    alnFastaInFH = getInputTempFile(alignmentFasta)
    alnPhylipOutFH = getOutputTempFile()
    
    input_handle = open(alnFastaInFH.name,"rU")
    output_handle = open(alnPhylipOutFH.name,"w")
    
    alignments = AlignIO.parse(input_handle, "fasta")
    AlignIO.write(alignments, output_handle , "phylip")
    
    input_handle.close()
    output_handle.close()
    treeInFH = getInputTempFile(newicktree)
    
    #removes the output files to be created if they already exist
    if os.path.exists(alnPhylipOutFH.name+"_phyml_stats.txt"):
        os.system("rm %s_phyml_stats.txt" % (alnPhylipOutFH.name))
    if os.path.exists(alnPhylipOutFH.name+"_phyml_tree.txt"):
        os.system("rm %s" % (alnPhylipOutFH.name+"_phyml_tree.txt"))
    
    #spawns a process and executes all the correct input keys when prompted
    child = pexpect.spawn("phyml")
    child.expect(". Enter the sequence file name >")
    child.sendline(alnPhylipOutFH.name)
    child.sendline("D")
    child.sendline("+")
    child.sendline("+")
    child.sendline("O")
    child.sendline("U")
    child.sendline("Y")
    child.expect(". Enter the name of the input tree file >")
    child.sendline(treeInFH.name)
    
    #checks for how long the calculation has been going for, if it has stalled, then the process will be halted
    argc = 0
    startTime = time.time()
    B = True
    while child.isalive() and B:
        argc += 1
        if argc % 1000 == 0:
            newTime = time.time() - startTime
            if newTime > 60.0:
                B = False
    
    Ret = None
    if B:
        Ret = open(alnPhylipOutFH.name+"_phyml_tree.txt" , "r").read()
    #returns the branch lengthed tree if the process was successful
    return [B , Ret]

def usePhyML(alignmentFasta):
    #converts the alignment file in FASTA to PHYLIP format
    alnFastaInFH = getInputTempFile(alignmentFasta)
    alnPhylipOutFH = getOutputTempFile()
    
    input_handle = open(alnFastaInFH.name,"rU")
    output_handle = open(alnPhylipOutFH.name,"w")
    
    alignments = AlignIO.parse(input_handle, "fasta")
    AlignIO.write(alignments, output_handle , "phylip")
    
    input_handle.close()
    output_handle.close()
    
    #print alnPhylipOutFH.read()
    
    #removes the output files to be created if they already exist
    if os.path.exists(alnPhylipOutFH.name+"_phyml_stats.txt"):
        os.system("rm %s_phyml_stats.txt" % (alnPhylipOutFH.name))
    if os.path.exists(alnPhylipOutFH.name+"_phyml_tree.txt"):
        os.system("rm %s" % (alnPhylipOutFH.name+"_phyml_tree.txt"))
    
    #print alnPhylipOutFH.read()
    
    #spawns a process and executes all the correct input keys when prompted
    child = pexpect.spawn("phyml")
    child.expect(". Enter the sequence file name >")
    child.sendline(alnPhylipOutFH.name)
    child.sendline("D")
    child.sendline("Y")
    
    #checks for how long the calculation has been going for, if it has stalled, then the process will be halted
    argc = 0
    startTime = time.time()
    B = True
    while child.isalive() and B:
        argc += 1
        if argc % 1000 == 0:
            newTime = time.time() - startTime
            if newTime > 60.0:
                B = False
    
    Ret = None
    if B:
        Ret = open(alnPhylipOutFH.name+"_phyml_tree.txt" , "r").read()
    #returns the branch lengthed tree if the process was successful
    return [B , Ret]
    
    
def entropy(prob_L):
    s = 0.0
    
    for prob in prob_L:
        if prob == 0.0:
            pass
        else:
            s -= prob * (math.log(prob))
    return s

###NCBI TAXONOMY####

def getNodeNameToNodeNumberAndReverse_D(PATH):
    allLines_L = [l.replace("\n","").split("\t|\t") for l in open(PATH,"r").readlines()]
    
    NameToNumber_D = {}
    NumberToName_D = {}
    
    for ls in allLines_L:
        NameToNumber_D[ls[1]] = ls[0]
        NumberToName_D[ls[0]] = ls[1]
    
    return {"NameToNumber" : NameToNumber_D,\
            "NumberToName" : NumberToName_D}

def getChildToParent_D(PATH):
    allLines_L = [l.replace("\n","").split("\t|\t") for l in open(PATH,"r").readlines()]
    
    ChildToParent_D = {}
    
    for ls in allLines_L:
        ChildToParent_D[ls[0]] = ls[1]
    
    return ChildToParent_D

def getListOfNodesToOrigin(nodeOfInterest , ChildToParent_D):
    Nodes_L = []
    
    nodeCheck = nodeOfInterest
    B = True
    
    while B:
        parentNode = ChildToParent_D[nodeCheck]
        
        if parentNode == "1":
            B = False
        
        nodeCheck = parentNode
        Nodes_L.append(parentNode)
    
    return Nodes_L
    
def parseInparanoidTable(PATH , SpecA , SpecB):
    initial_D = {}
    F = open(PATH,"r")
    allLines_L = F.readlines()
    F.close()
    for line in allLines_L:
        ls = line.replace("\n","").split()
        
        if ls[0] in initial_D.keys():
            pass
        else:
            initial_D[ls[0]] = []
        
        initial_D[ls[0]].append(ls)
    
    Ret = []
    for key in initial_D.keys():
        SpecALines = []
        SpecBLines = []
        
        for ls in initial_D[key]:
            if ls[2] == SpecA:
                SpecALines.append(ls[4])
            elif ls[2] == SpecB:
                SpecBLines.append(ls[4])
        
        Ret.append({SpecA : SpecALines , SpecB : SpecBLines})
    
    return Ret
    
def makePoint(pdbXMLPointSplit):
        Ret = None
        if pdbXMLPointSplit[0] == "NA" or pdbXMLPointSplit[1] == "NA" or pdbXMLPointSplit[2] == "NA":
            pass
        else:
            Ret = [float(pdbXMLPointSplit[0]) , float(pdbXMLPointSplit[1]) , float(pdbXMLPointSplit[2])]
        return Ret

def computeVector(pointa , pointb):
    return pointb - pointa

def compute3dVector(pointa,pointb):
    return [computeVector(pointa[0],pointb[0]) , computeVector(pointa[1],pointb[1]) , computeVector(pointa[2],pointb[2])]
  
"computes the magnitude of distance between two 3-D points"
def getDistanceMagnitude(pointa , pointb):
    return math.sqrt(math.pow(computeVector(pointa[0] , pointb[0]) , 2) +\
                     math.pow(computeVector(pointa[1] , pointb[1]) , 2) +\
                     math.pow(computeVector(pointa[2] , pointb[2]) , 2))

def getDistanceMagnitudeForVector(veca):
    return math.sqrt(math.pow(veca[0],2) + math.pow(veca[1],2) + math.pow(veca[2],2))
    
def getCombinatorialListOfPairwiseDistances(PointsA_L):
        
    PointsB_L = PointsA_L[1:]
    count = 0
    Distances_L = []
    
    
    #combinations of pairwise distances
    for PointsA in PointsA_L:
        for PointsB in PointsB_L:
            #gets both points
            
            #as long as they are not null, calculate the magnitude of distance between them and add it to a list
            if PointsA and PointsB:
                Distances_L.append(getDistanceMagnitude(PointsA , PointsB))
        PointsB_L = PointsB_L[1:]
        
    return Distances_L

def getCombinatorialListOfCoplanarities(VecA_L):
    
    VecB_L = VecA_L[1:]
    count = 0
    Coplanar_L = []
    
    
    #combinations of pairwise distances
    for VecA in VecA_L:
        for VecB in VecB_L:
            #gets both points
            
            #as long as they are not null, calculate the magnitude of distance between them and add it to a list
            if VecA and VecB:
                vecAMag = getDistanceMagnitudeForVector(VecA)
                vecBMag = getDistanceMagnitudeForVector(VecB)
                Dot = numpy.dot(vecAMag,vecBMag)
                costheta = Dot /(vecAMag * vecBMag)
                arccos = numpy.arccos(costheta)
                deg = numpy.rad2deg(arccos)
                
                if numpy.isnan(deg):
                    pass
                else:
                
                    Coplanar_L.append(deg)
        VecB_L = VecB_L[1:]
        
    return Coplanar_L





def getSASD():
    return  {'A': 102.77, 'R': 279.20, 'N': 174.75, 'D': 172.69,\
            'C': 143.36, 'Q': 220.26, 'E': 214.78, 'G':  47.92,\
            'H': 221.38, 'I': 241.44, 'L': 205.91, 'K': 246.36,\
            'M': 213.64, 'F': 224.78, 'P': 170.15, 'S': 126.89,\
            'T': 182.33, 'W': 243.47, 'Y': 261.73, 'V': 204.59,\
            'X': 1.0}

def getBlosum62_D():
    return {"A":[["AR",0.0],["AN",0.0542132467845],["AD",0.102536362544],["AC",0.158672917418],["AQ",0.185358817166],["AE",0.226229114076],["AG",0.289578074288],["AH",0.391513403053],["AI",0.415674960933],["AL",0.471691309052],["AK",0.546940738069],["AM",0.617261690107],["AF",0.643827383099],["AP",0.68193292463],["AS",0.725087149898],["AT",0.811515807188],["AW",0.873302079577],["AY",0.886043995673],["AV",0.923668710182]],\
            "R":[["RA",0.0],["RN",0.0557961153037],["RD",0.102313497464],["RC",0.158604478535],["RQ",0.184337498454],["RE",0.249907212669],["RG",0.314363478906],["RH",0.380675491773],["RI",0.411975751577],["RL",0.456142521341],["RK",0.534331312631],["RM",0.640356303353],["RF",0.670171965854],["RP",0.716565631572],["RS",0.762835580849],["RT",0.830508474576],["RW",0.876902140294],["RY",0.893603859953],["RV",0.935172584436]],\
            "N":[["NA",0.0],["NR",0.0467496220491],["ND",0.0904756367019],["NC",0.172345621584],["NQ",0.205140132574],["NE",0.259681358297],["NG",0.319688335853],["NH",0.427956739156],["NI",0.473543435283],["NL",0.530875683219],["NK",0.59658099779],["NM",0.670543086405],["NF",0.70252354925],["NP",0.745319223165],["NS",0.786021630422],["NT",0.854285382021],["NW",0.911152459588],["NY",0.920339574369],["NV",0.952203744621]],\
            "D":[["DA",0.0],["DR",0.0605628323175],["DN",0.119569446246],["DC",0.210867591752],["DQ",0.23395149786],["DE",0.276228764103],["DG",0.381662559979],["DH",0.456361042666],["DI",0.489300998573],["DL",0.531318895085],["DK",0.60822202049],["DM",0.681364284788],["DF",0.698223317339],["DP",0.73103358838],["DS",0.782648164959],["DT",0.850992089223],["DW",0.902995720399],["DY",0.913240824796],["DV",0.944365192582]],\
            "C":[["CA",0.0],["CR",0.0356283100626],["CN",0.0690097897609],["CD",0.114267372813],["CQ",0.142834216017],["CE",0.174129353234],["CG",0.259027443428],["CH",0.330284063553],["CI",0.383726528647],["CL",0.430749478414],["CK",0.526079281014],["CM",0.599582731504],["CF",0.628470550473],["CP",0.670357887979],["CS",0.710479858771],["CT",0.770983790724],["CW",0.843684801797],["CY",0.891670678864],["CV",0.934520943669]],\
            "Q":[["QA",0.0],["QR",0.0385356454721],["QN",0.098605916355],["QD",0.15176243908],["QC",0.188711322679],["QE",0.210812648759],["QG",0.317352374476],["QH",0.385356454721],["QI",0.414484869092],["QL",0.456307378443],["QK",0.52680494163],["QM",0.600476028562],["QF",0.635611470022],["QP",0.672447013487],["QS",0.725376855945],["QT",0.801428085685],["QW",0.853337866939],["QY",0.897993879633],["QV",0.956250708376]],\
            "E":[["EA",0.0],["ER",0.0634633911368],["EN",0.126204238921],["ED",0.188342967245],["EC",0.286247591522],["EQ",0.349951830443],["EG",0.463150289017],["EH",0.520472061657],["EI",0.548410404624],["EL",0.584898843931],["EK",0.593208092486],["EM",0.682080924855],["EF",0.700264932563],["EP",0.731575144509],["ES",0.787090558767],["ET",0.852360308285],["EW",0.90655105973],["EY",0.919677263969],["EV",0.952673410405]],\
            "G":[["GA",0.0],["GR",0.118468846046],["GN",0.193350097793],["GD",0.323414361554],["GC",0.403883766415],["GQ",0.465912265996],["GE",0.549734562727],["GH",0.6162335848],["GI",0.644314054205],["GL",0.681056160939],["GK",0.682592903046],["GM",0.731628946633],["GF",0.756496227997],["GP",0.791422184968],["GS",0.843392008941],["GT",0.897177982677],["GW",0.945794914781],["GY",0.959905001397],["GV",0.99804414641]],\
            "H":[["HA",0.0],["HR",0.0262539184953],["HN",0.0592998955068],["HD",0.110501567398],["HC",0.14367816092],["HQ",0.187173458725],["HE",0.220741901776],["HG",0.251044932079],["HI",0.277298850575],["HL",0.323275862069],["HK",0.406870428422],["HM",0.4644723093],["HF",0.559952978056],["HP",0.619252873563],["HS",0.691353187043],["HT",0.763714733542],["HW",0.823537095089],["HY",0.879179728318],["HV",0.956504702194]],\
            "I":[["IA",0.0],["IR",0.0542870456664],["IN",0.0958760484623],["ID",0.153308480895],["IC",0.191053122088],["IQ",0.22518639329],["IE",0.26817334576],["IG",0.303471575023],["IH",0.334109972041],["IL",0.375116495806],["IK",0.515843429637],["IM",0.567218080149],["IF",0.599021435228],["IP",0.657269338304],["IS",0.70246971109],["IT",0.756290773532],["IW",0.808946877912],["IY",0.81931500466],["IV",0.863350419385]],\
            "L":[["LA",0.0],["LR",0.0784559468605],["LN",0.157663867653],["LD",0.228474746209],["LC",0.302794836446],["LQ",0.377240255671],["LE",0.455194886577],["LG",0.46384258679],["LH",0.465221205665],["LI",0.545431758366],["LK",0.696829176589],["LM",0.721268329365],["LF",0.76450683043],["LP",0.837949617747],["LS",0.874044366462],["LT",0.87868153904],["LW",0.935956886828],["LY",0.951121694448],["LV",0.998621381125]],\
            "K":[["KA",0.0],["KR",0.0683650812201],["KN",0.168517003623],["KD",0.242842117565],["KC",0.308753067664],["KQ",0.362276498773],["KE",0.438237700129],["KG",0.524482879514],["KH",0.565501928246],["KI",0.617038681781],["KL",0.668575435316],["KM",0.691363795723],["KF",0.718709828211],["KP",0.76089751081],["KS",0.815122122239],["KT",0.882084842819],["KW",0.939932219236],["KY",0.9543064158],["KV",0.98936543181]],\
            "M":[["MA",0.0],["MR",0.0235406902429],["MN",0.0492117596932],["MD",0.0785044737963],["MC",0.0923519386451],["MQ",0.111525351513],["ME",0.144546229229],["MG",0.160630592245],["MH",0.179590967192],["MI",0.257456327226],["ML",0.286536003409],["MK",0.323285044738],["MF",0.348210481466],["MP",0.411695781849],["MS",0.442479761398],["MT",0.494674051981],["MW",0.54825308905],["MY",0.873242437154],["MV",0.924051981253]],\
            "F":[["FA",0.0],["FR",0.0391503025812],["FN",0.0854637520069],["FD",0.13091268371],["FC",0.162158824256],["FQ",0.194392985056],["FE",0.234531307892],["FG",0.26664196616],["FH",0.297517599111],["FI",0.353587748549],["FL",0.41533901445],["FK",0.487711498086],["FM",0.532295912066],["FP",0.60590342102],["FS",0.636902556502],["FT",0.696307274299],["FW",0.747807830061],["FY",0.804989502285],["FV",0.916265283438]],\
            "P":[["PA",0.0],["PR",0.050792303339],["PN",0.103706847765],["PD",0.153225806452],["PC",0.209535936616],["PQ",0.244906621392],["PE",0.310979060555],["PG",0.376202603282],["PH",0.42883418223],["PI",0.506932654216],["PL",0.561827956989],["PK",0.602574985852],["PM",0.668222976797],["PF",0.709111488398],["PS",0.744623655914],["PT",0.808290888512],["PW",0.880305602716],["PY",0.906338426712],["PV",0.942416525184]],\
            "S":[["SA",0.0],["SR",0.0828627405785],["SN",0.145902961853],["SD",0.213553071338],["SC",0.274288348508],["SQ",0.317736544889],["SE",0.395067419615],["SG",0.457531404863],["SH",0.501901578887],["SI",0.565748530598],["SL",0.618992739426],["SK",0.623256886021],["SM",0.689293534632],["SF",0.745764665207],["SP",0.801198570935],["ST",0.8530598133],["SW",0.931082171257],["SY",0.942261150167],["SV",0.978218278207]],\
            "T":[["TA",0.0],["TR",0.0581908751274],["TN",0.10064530737],["TD",0.156005887015],["TC",0.201403826559],["TQ",0.252688780709],["TE",0.304539793954],["TG",0.355485112646],["TH",0.394882825767],["TI",0.446733839013],["TL",0.497905581343],["TK",0.549643382769],["TM",0.60568323333],["TF",0.662628778444],["TP",0.709838107098],["TS",0.767462923129],["TW",0.844107324805],["TY",0.883505037926],["TV",0.915317559153]],\
            "W":[["WA",0.0],["WR",0.0129095116307],["WN",0.0293508707831],["WD",0.0389721105834],["WC",0.0485933503836],["WQ",0.08500791621],["WE",0.132992327366],["WG",0.146267202533],["WH",0.158567774936],["WI",0.21044939715],["WL",0.221288515406],["WK",0.23602484472],["WM",0.251004749726],["WF",0.622579466569],["WP",0.67896723907],["WS",0.701376202655],["WT",0.713189623676],["WY",0.755571793935],["WV",0.948970892705]],\
            "Y":[["YA",0.0],["YR",0.0365654205607],["YN",0.0758177570093],["YD",0.107827102804],["YC",0.135864485981],["YQ",0.167056074766],["YE",0.227102803738],["YG",0.259112149533],["YH",0.291004672897],["YI",0.360163551402],["YL",0.404322429907],["YK",0.448598130841],["YM",0.483644859813],["YF",0.539369158879],["YP",0.644626168224],["YS",0.67441588785],["YT",0.710864485981],["YW",0.743691588785],["YV",0.929205607477]],\
            "V":[["VA",0.0],["VR",0.0737676261508],["VN",0.134832770073],["VD",0.182729285631],["VC",0.232839995339],["VQ",0.280503437828],["VE",0.325719613099],["VG",0.371868080643],["VH",0.373616128656],["VI",0.412422794546],["VL",0.549120149167],["VK",0.550518587577],["VM",0.561589558327],["VF",0.644796643748],["VP",0.723808413938],["VS",0.771238783359],["VT",0.793264188323],["VW",0.880200442839],["VY",0.929145787204]]}

def getMutation(Blosum62,From):
    sub_L = Blosum62[From]
    randnum = random.random()
    
    mutType = ""
    keepGoing = True
    for elem in sub_L:
        if keepGoing:
            if randnum >= elem[1]:
                mutType = elem[0]
            else:
                keepGoing = False
    
    
    return mutType
    


