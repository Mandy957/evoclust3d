####################################################################################################
#                                                                                                  #
#           PROJECT Protein Adaptation                                                             #
#           CLASS GenomePrep                                                                       #
#           PROGRAMMER Jeremy Adams                                                                #
#           STARTED 14-11-14                                                                       #
#           LASTMOD 14-11-14                                                                       #
#                                                                                                  #
#           DESCRIPTION class representation of all work necessary to extract orthologous          #
#                       relationships out of a series of whole genomes using Mercator              #
#                                                                                                  #
####################################################################################################

import sys
import os
from cogent import LoadTree
from staticmethods import *

class GenomePrep:
    
    """
    Class attributes:
    ExitString (String): contains a string to provide the user with info as to why the program was not executed successfully
    ExitStatus (Bool): True if executed successfully, false if not
    DIR (String): Path to the Genome general directory
    IncompleteDIR (String): Path to incomplete genome directories
    RawDIR (String): Path to unmodified complete genomes
    ModdedDIR (String): Path to modified complete genomes
    SDB (String): Path to modified complete genomes in SDB format
    GFF (String): Path to GFF annotation files
    MercatorProjects (String): Path to mercator projects directory
    ProjectName (String): String reflecting the specific desired merctor project (and directory in the Mercator Projects folder)
    ProjectDIR (String): Path to the specific mercator project of interest
    SpeciesNamePATH (String): Path to the species names file
    SpeciesTreePATH (String): Path to the species tree file
    """
    
    "CONSTRUCTOR"
    def __init__(self , DataDIR , ProjectName):
        
        #establishes all correct paths to necessary files
        self.ExitString = ""
        self.ExitStatus = False
        
        self.DataDIR = DataDIR
        self.DIR = self.DataDIR+"Genomes/"
        
        self.IncompleteDIR = self.DIR+"Incomplete/"
        self.RawDIR = self.DIR+"RawComplete/"
        self.ModdedDIR = self.DIR+"ModdedComplete/"
        self.SDB = self.DIR+"SDB/"
        
        self.GFF = self.DIR+"GFF/"
        self.MercatorProjects = self.DIR+"MercatorProjects/"
        
        self.ProjectName = ProjectName
        if self.ProjectName.endswith("/"):
            self.ProjectName = self.ProjectName[:-1]
        self.ProjectDIR = self.MercatorProjects+self.ProjectName+"/"
        
        self.SpeciesNamePATH = self.ProjectDIR+"SpeciesNames.txt"
        self.SpeciesTreePATH = self.ProjectDIR+"SpeciesTree.nwk"
        
        #checks if the project directory and the species name and tree files exist
        if os.path.exists(self.ProjectDIR):
            if os.path.exists(self.SpeciesNamePATH) and os.path.exists(self.SpeciesTreePATH):
                
                #reads the species name file
                self.SpeciesName_L = self.readSpeciesFile()
                
                #reads the species tree file and formats it into cogent format, then parses it accordng to the static method
                SPCogentFixUp = fixUpFileForCogent(self.SpeciesTreePATH)
                self.SPCogentTreeFile = SPCogentFixUp[0]
                self.SPCogentInputTreeString = SPCogentFixUp[1]
                
                self.SpeciesCogentTree = LoadTree(self.SPCogentTreeFile.name)
                
                NodesLeaveBranches = completeNodesLeavesBranches(self.SpeciesCogentTree)
                self.SPTopKey = NodesLeaveBranches["TopKey"]
                self.SPNodes_D = NodesLeaveBranches["Nodes_D"]
                self.SPNodeKey_L = NodesLeaveBranches["NodeKey_L"]
                self.SPUpperKey_L = NodesLeaveBranches["UpperKey_L"]
                self.SPLeafKey_L = NodesLeaveBranches["LeafKey_L"]
                self.SPBranchKey_L = NodesLeaveBranches["BranchKey_L"]
                
                #checks if the species file and tree file match
                if set(self.SpeciesName_L) == set(self.SPLeafKey_L):
                    
                    #modifies genomes that need to be modified to be ready for mercator
                    self.NeedsToMod_L = self.getNeedsToMod()
                    self.checkIncompletes_D = self.checkForIncompletes()
                    
                    #checks if program can proceed (ie. all genomes are present even if they are as yet incomplete)
                    if self.checkIncompletes_D["Proceed"]:
                        
                        #preps all incomplete genomes that need to be prepped
                        for Spec in self.NeedsToMod_L:
                            self.performFullPrep(Spec)
                        
                        #checks for the presence of all GFF files
                        self.checkGFFs_D = self.checkForGFFs()
                        
                        #if all the checks are passed, then the appropriate inputs are made and mercator is run
                        if self.checkGFFs_D["Proceed"]:
                            self.makeMercatorProjectInput()
                            self.runMercator()
                            
        #appropriate error message is released for the error that took place               
                        else:
                            self.ExitString = self.exitMissingGFF()
                    else:
                        self.ExitString = self.exitMissingGenomes()
                else:
                    self.ExitString = self.exitUnmatchedFiles()
                
            else:
                self.ExitString = self.exitNoInputFiles()
            
        else:
            self.ExitString = self.exitNoProjectDIR()
                
        print self.ExitString
    
    "reads all species in the species file"
    def readSpeciesFile(self):
        return [line.replace("\n","") for line in open(self.SpeciesNamePATH,"r").readlines()]
    "gets a list of all species that do not have complete genomes in SDB format"
    def getNeedsToMod(self):
        return [spec for spec in self.SpeciesName_L if os.path.exists(self.SDB+spec+".sdb") == False]
    "checks for all genomes in the project of interest that are missing from the appropriate directory"
    def checkForIncompletes(self):
        Missing = [spec for spec in self.NeedsToMod_L if os.path.exists(self.IncompleteDIR+spec) == False]
        B = True
        if len(Missing) > 0:
            B = False
        return {"Proceed" : B , "Missing" : Missing}
    
    "checks the gff directory for the presence of all needed GFF files"
    def checkForGFFs(self):
        Missing = [spec for spec in self.SpeciesName_L if os.path.exists(self.GFF+spec+".gff") == False]
        B = True
        if len(Missing) > 0 :
            B = False
        return {"Proceed":B,"Missing":Missing}
    
    "modifies genomes in separate files into one genome file" 
    def incompleteToRawComplete(self, species):
        
        print "Concatenating input genome files for %s" % (species)
        Out = []
        
        for F in sorted(os.listdir(self.IncompleteDIR+species)):        
            Out.append(open(self.IncompleteDIR+species+"/"+F,"r").read())
            open(self.RawDIR+species+".fa" , "w").write("\n".join(Out))
        
    "modifies the fasta header for the fasta sequences in the completed genome file"
    def rawCompleteToModdedComplete(self, species):
        print "Modifying header names for %s" % (species)
        
        F = readFasta(self.RawDIR+species+".fa")
        
        NewHeaderNames = []
        Output_L = []
        
        if len(F[0]) == len(F[1]):
            B = True
            
            for Key in F[0]:
                if B:
                    NewHeader = Key.split("|")[3]
                    
                    
                    if NewHeader in set(NewHeaderNames):
                        B = False
                    else:
                        NewHeaderNames.append(NewHeader)
                        
                    Output_L.append(">%s\n%s" % (NewHeader , re.sub("(.{70})",r"\1\n",F[1][Key])))
            
            if B:
                open(self.ModdedDIR+species+".fa" , "w").write("\n\n".join(Output_L))
            
            else:
                print "New scaffold name duplicate"
        else:
            print "Error: Non-unique scaffold names"
    
    "converts complete genomes into SDB format that mercator uses"
    def moddedCompleteToSDB(self, species):
        print "Converting %s FASTA file to SDB" % (species)
        
        INF = species + ".fa"
        OUTF = species + ".sdb"
        
        os.chdir(self.ModdedDIR)
        
        os.system("fa2sdb %s < %s" % (OUTF,INF))
        
        os.system("mv %s %s" % (OUTF,self.SDB))
    
    "performs the 3 above methods to fully convert a genome from separate files into a genome in SDB format"
    def performFullPrep(self, species):
        self.incompleteToRawComplete(species)
        self.rawCompleteToModdedComplete(species)
        self.moddedCompleteToSDB(species)
    
    "executes the makeMercatorInput program"   
    def makeMercatorProjectInput(self):
        os.system("makeMercatorInput --genome-dir=%s --gff-dir=%s --out-dir=%s %s" % (self.SDB,self.GFF,self.ProjectDIR," ".join(self.SpeciesName_L)))
    "executes the mercator program"
    def runMercator(self):
        os.chdir(self.ProjectDIR)
        os.system("mercator %s" % (" ".join(self.SpeciesName_L)))
    
    "various error messages"
    def exitNoProjectDIR(self):
        return """Error: The project directory you specified (%s) is not in the MercatorProjects directory.""" % (self.ProjectName)
    def exitNoInputFiles(self):
        return """Error: The files SpeciesNames.txt and SpeciesTree.nwk are not present in the %s project directory""" % (self.ProjectName)
    def exitUnmatchedFiles(self):
        return """Error: Species names in the SpeciesNames.txt file and the SpeciesTree.nwk file are unmatched"""
    def exitMissingGenomes(self):
        return """Error: There are missing genomes for the following species: %s""" % (" ".join(self.checkIncompletes_D["Missing"]))
    def exitMissingGFF(self):
        return """Error: There are missing GFF files for the following species: %s""" % (" ".join(self.checkGFFs_D["Missing"]))
    
    
    