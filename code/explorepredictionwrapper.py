####################################################################################################
#                                                                                                  #
#           PROJECT Protein Adaptation                                                             #
#           CLASS ExplorePredictionWrapper                                                         #
#           PROGRAMMER Jeremy Adams                                                                #
#           STARTED 09-12-14                                                                       #
#           LASTMOD 09-12-14                                                                       #
#                                                                                                  #
#           DESCRIPTION class to generate potentially multiple set of figures for all ancestral,   #
#                       derived, PDB sequence triads that are significant hits                     #
#                                                                                                  #
####################################################################################################

from exploreprediction import ExplorePrediction
from staticmethods import getCombinedPValue
import os

class ExplorePredictionWrapper:
    
    """
    Class attributes:
    Directory (String): Directory to the analyzed protein family of interest
    PValueFile (String): Path to the P-value file
    SigPValues_L (List): List of all derived clades of interest that had a significant combined P-value
    """
    
    "CONSTRUCTOR"
    def __init__(self, Directory):
        
        #sets the directory and checks if it ends with a "/"
        self.Directory = Directory
        
        if self.Directory.endswith("/"):
            pass
        else:
            self.Directory = self.Directory+"/"
        
        #makes the figures directory if it does not already exist
        if os.path.exists(self.Directory+"Figures"):
            pass
        else:
            os.system("mkdir " +self.Directory+"Figures")
        
        self.PValueFile = self.Directory+"PValues.txt" #sets the path to the p-value file
        self.SigPValues_L = self.getSigPValues() #gets significant clades of interest
        self.explorePredictions() #instantiates an explore prediction class for each clade that was a significant hit
        
    "method to get clades that are significant hits according to the algorithm"
    def getSigPValues(self):
        R = []
        #for each p-value line in the file
        for line in [line.replace("\n","") for line in open(self.PValueFile,"r").readlines()][1:]:
            #executes the combined p-value method and checks if the combined value is less than 0.05
            ls = line.split()
            PValue = float(ls[-1])
            
            
            if PValue <= 0.05:
                R.append([ls[0].split(">>")[1] , ls[1]])
        return R
    
    "method to instantiate an ExplorePrediction class for each significant hit"
    def explorePredictions(self):
        
        for pred in self.SigPValues_L:
            ExplorePrediction(self.Directory , pred[0] , pred[1])
        
        
        
        
            
            
            
            
        
    
    
    