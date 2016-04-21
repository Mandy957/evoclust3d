####################################################################################################
#                                                                                                  #
#           PROJECT Protein Adaptation                                                             #
#           CLASS FastMLTree                                                                       #
#           PROGRAMMER Jeremy Adams                                                                #
#           STARTED 06-08-14                                                                       #
#           LASTMOD 06-08-14                                                                       #
#                                                                                                  #
#           DESCRIPTION Takes in user-defined protein tree, parses all nodes and branches, and     #
#                       matches said nodes to equivalent nodes from FastML Results                 #
#                                                                                                  #
####################################################################################################

from cogent import LoadTree
from staticmethods import *
import re
import sets

class FastMLTree:
    
    """
    Class attributes:
    Parsed (Bool): an indication of whether or not the user-defined tree was successfully parsed, if it was not, then the rest of the analysis is not performed
    TreePath (String): absolute path to tree file
    NeedsToBeCogentModded (Bool): whether or not placeholder names for the internal nodes need to be created
    CogentTree (Object LoadTree): pyCogent Class object containing parsed newick syntax tree
    FastMLInputTreeString (String): representation of tree in newick with internal node names removed
    FastMLOutputTreeString (String): representation of tree in newick with internal nodes named according to FastML naming convention
    FastMLToOriginalMatchedNodes_D (Dict): Key is the node name in the cogent convention, value is the node name in the FastML convention
    
    NodeKey_L (List): List of all node name keys
    LeafKey_L (List): List of all terminal node name keys
    UpperKey_L (List): List of all internal (non-terminal) node name keys
    TopKey (String): root node name key
    BranchKey_L (List): List of all paths (from ancestral to immediate derived) along the tree
    Nodes_D (Dict): Key is the node name, value is a sub-dict containing immediate derived nodes and terminal nodes under the node
    """
    
    "CONSTRUCTOR"
    def __init__(self, TreePath , NeedsToBeCogentModded):
        self.Parsed = True #used to determine if the full analysis can be conducted
        
        try:
            self.TreePath = TreePath
            self.NeedsToBeCogentModded = NeedsToBeCogentModded
            
            self.CogentTree = None
            
            #if the internal nodes need to be renamed, then it is done according to the "FixUpFileForCogent" method
            if self.NeedsToBeCogentModded:
                cogentFixUp = fixUpFileForCogent(self.TreePath)
                self.CogentTreeFile = cogentFixUp[0]
                self.CogentInputTreeString = cogentFixUp[1]
                
                
                self.CogentTree = LoadTree(self.CogentTreeFile.name)
                
            else:
                
                self.CogentTree = LoadTree(self.TreePath)
            
            #prepares an input string for FastML
            self.FastMLInputTreeString = self.FixUpFileForFastML(self.CogentTree)
            
            
            #executes method to fully parse tree, then sets all returned variables as class variables
            CogentNodesLeavesBranches = completeNodesLeavesBranches(self.CogentTree)
            self.NodeKey_L = CogentNodesLeavesBranches['NodeKey_L']
            self.LeafKey_L = CogentNodesLeavesBranches['LeafKey_L']
            self.UpperKey_L = CogentNodesLeavesBranches['UpperKey_L']
            self.TopKey = CogentNodesLeavesBranches['TopKey']
            self.BranchKey_L = CogentNodesLeavesBranches['BranchKey_L']
            self.Nodes_D = CogentNodesLeavesBranches['Nodes_D']
            
            
            
            
            
            #print self.LeafKey_L
            #executes quick run of FastML to get FastML's naming convention of internal nodes
            
            self.FastMLOutputTreeString = executeFastML(self.getTempFASTAFile() , self.FastMLInputTreeString , True)
            
            
            #prepares the FastMLToOriginalMatchedNodes_D
            self.MatchNodes()
            
        except Exception as e:
            
            self.Parsed = False
        
    
    "Removes internal node names so that FastML adds its own naming convention"
    def FixUpFileForFastML(self, CogentTree):
        #gets the tree string for the cogent object
        TreeString = CogentTree.getNewick(with_distances=True).replace("'","")
        
        i = 0
        NotThroughTheString = True
        #while loop moves one space along tree string until it gets to the end
        while NotThroughTheString:
            #when a close bracket is found, it signifies the end of an internal node
            if TreeString[i] == ")":
                if TreeString[i+1] == ";":
                    pass
                else:
                    #tree string replaces the name of the internal node with nothing
                    lengthToColon = len(re.compile("^(.+?)[:;]").search(TreeString[i:]).group(1)) - 1
                    
                    TreeString = TreeString[:i+1]+ TreeString[i+lengthToColon+1:]
            #check to end while loop
            if i == len(TreeString) - 1:
                NotThroughTheString = False
            i += 1
        
        return TreeString
    
    "Prepares simple FastaFile to be given to FastML"   
    def getTempFASTAFile(self):
        retString_L = []
        
        #FastaFile will have the sequence "GREAT" for each terminal sequence
        for LeafKey in self.LeafKey_L:
            retString_L.append(">"+LeafKey)
            retString_L.append("GREAT")
        
        return '\n'.join(retString_L)
    
    "Corrects for instances where FastML anomalously renames terminal nodes"
    def correctForFastMLNameChanges(self):
        
        #gets lists of terminal names in the FastML input and output strings (in the same order)
        FastMLInputNames = [re.compile("^(.+?):").search(TaxString).group(1) for TaxString in re.findall("[A-Za-z0-9_./]+:[.0-9]+",self.FastMLInputTreeString)]
        #print FastMLInputNames
        FastMLOutputNames = [re.compile("^(.+?):").search(TaxString).group(1) for TaxString in re.findall("[A-Za-z0-9_./]+:[.0-9]+",self.FastMLOutputTreeString)]
        FastMLOutputNames = [Name for Name in FastMLOutputNames if re.compile("^N[0-9]+$").search(Name) == None]
        
        #when equivalent node names are not the same, then the output string node name is renamed according to the input string node name
        for i in range(0,len(FastMLInputNames)):
            if FastMLInputNames[i] != FastMLOutputNames[i]:
                self.FastMLOutputTreeString = re.sub("([,\(\)])%s:" % (FastMLOutputNames[i]) , r"\1%s:" % (FastMLInputNames[i]) , self.FastMLOutputTreeString)
    
    "Matches original (cogent) node names with how the nodes are named in FastML" 
    def MatchNodes(self):
        #print "YAY"
        self.correctForFastMLNameChanges() #performs the correction on the output string if necessary
        #print "NAY"
        TerminiStringToNodeName_D = {}
        #a termini string is prepared for each internal node, that is, all termini under the internal node sorted an placed into a single string
        
        for NodeKey in self.UpperKey_L:
            TerminiStringToNodeName_D['-'.join(sorted(self.Nodes_D[NodeKey]['terminal']))] = NodeKey
        
        #prepares a cogent tree object for the fastML output
        FH = getInputTempFile(self.FastMLOutputTreeString)
        
        FastMLCogentTree = LoadTree(FH.name)
        
        
        self.FastMLToOriginalMatchedNodes_D = {}
        
        #for each cogent node in the FastML cogent tree
        for FastMLCogentNodeKey in FastMLCogentTree.getNodeNames():
            
            #a termini string is prepared for the fastML node
            FastMLCogentNode = FastMLCogentTree.getNodeMatchingName(FastMLCogentNodeKey)
            FastMLTermini_L = [tip.Name for tip in FastMLCogentNode.iterTips()]
            
            #if it has more than 0 termini under the node
            if len(FastMLTermini_L) > 0:
                #A fastML termini string is prepared, and this termini string will be the same termini string as the equivalent cogent node
                FastMLTerminiString = '-'.join(sorted(FastMLTermini_L))
                self.FastMLToOriginalMatchedNodes_D[FastMLCogentNodeKey] = TerminiStringToNodeName_D[FastMLTerminiString]
                
            #if it has no termini under it, then the node itself is a terminus and has the same name in FastML and Cogent
            else:
                self.FastMLToOriginalMatchedNodes_D[FastMLCogentNodeKey] = FastMLCogentNodeKey
    
    "Sets branch lengths of each node"
    def setBranchLengths(self):
        
        self.BranchLength_D = {}
        #gets the distance between a node and its immediate ancestor
        for NodeNameKey in self.NodeKey_L:
            HigherNode = self.CogentTree.getNodeMatchingName(NodeNameKey)
            
            for ImmediateNeighbourNodeNameKey in self.Nodes_D[NodeNameKey]['immediate']:
                LowerNode = self.CogentTree.getNodeMatchingName(ImmediateNeighbourNodeNameKey)
                
                self.BranchLength_D[ImmediateNeighbourNodeNameKey] = HigherNode.distance(LowerNode)
            