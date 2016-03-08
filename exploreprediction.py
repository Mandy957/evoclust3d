####################################################################################################
#                                                                                                  #
#           PROJECT Protein Adaptation                                                             #
#           CLASS ExplorePrediction                                                                #
#           PROGRAMMER Jeremy Adams                                                                #
#           STARTED 09-12-14                                                                       #
#           LASTMOD 09-12-14                                                                       #
#                                                                                                  #
#           DESCRIPTION class to generate relevant figures for an ancestral, derived, PDB,         #
#                       sequence triad                                                             #
#                                                                                                  #
####################################################################################################

from cogent import LoadTree
from scopealgorithm import ScopeAlgorithm
from fastmltree import FastMLTree
from staticmethods import getOutputTempFile, getAllPDBFileDicts, id_generator
import time
import re
import sets
import pymol
import os
import cairosvg

class ExplorePrediction:
    
    "CONSTRUCTOR"
    def __init__(self, Directory , DerivedoI , PDBoI):
        
        """
        Class attributes:
        Figures_L (List): list of all the figure types that will be created
        FiguresSVG_D (Dict): key is the type of figure, value is the SVG syntax that will draw the figure
        DerivedoInterest (String): Derived node of interest that the figure will be based on
        PDBoInterest (String): PDB structure that the derived node sequence aligned to and achieved a significant hit on
        """
        
        #initial setup of what figures will be created
        self.Figures_L = ["TreeAndStates" , "Alignment" , "Structurecartoon" , "Structuresurface"]
        self.FigureSVG_D = {Key : [] for Key in self.Figures_L}
        
        self.Directory = Directory
        if self.Directory.endswith("/"):
            pass
        else:
            self.Directory = self.Directory+"/"
            
        self.DerivedoInterest = DerivedoI
        self.PDBoInterest = PDBoI
        
        print self.Directory
        print self.DerivedoInterest
        print self.PDBoInterest
        
        #output directory where files will be written
        self.OutputDirectory = "%sFigures/%s-%s/" % (self.Directory,self.DerivedoInterest,self.PDBoInterest)
        
        if os.path.exists(self.OutputDirectory):
            pass
        else:
            os.system("mkdir " +self.OutputDirectory)
        
        #paths to relevant input files
        self.ReportPATH = self.Directory+"Report.xml"
        self.TreePATH = self.Directory+"ModdedTree.nwk"
        self.MatrixPATH = self.Directory+"ScoringMatrix.xml"
        
        #parses the report file for sequences and branch relationships
        self.NodeToSeq_D = {re.compile("<H>(.+?)</H>").search(Seq).group(1) : re.compile("<S>(.+?)</S>").search(Seq).group(1) for Seq in re.findall("<Seq>.+?</Seq>", open(self.ReportPATH,"r").read())}
        self.BranchToAlgorithm_D = {re.compile("<Branch_name>(.+?)</Branch_name>").search(Branch).group(1) : ScopeAlgorithm(Branch) for Branch in re.findall("<Branch>.+?</Branch>",open(self.ReportPATH,"r").read(),re.DOTALL)}
        self.RectCount = 0
        
        #dimensions
        self.TreeFigWIDTH = 750
        self.TreeFigHEIGHT = 500
        self.TreeFigXOffset = 25
        self.TreeFigYOffset = 50
        
        #loads and parses tree, gets evolutionary distances for proper branch lengths
        self.CogentTree = LoadTree(self.TreePATH)
        self.FastMLTree = FastMLTree(self.TreePATH , False)
        self.FastMLTree.setBranchLengths()
        self.LongestDistance = self.getLongestEvoDistance()
        self.EvoDistance_D = {Key : self.getEvoDistance(Key) for Key in self.NodeToSeq_D.keys() if Key != self.FastMLTree.TopKey}
        self.EvoDistance_D[self.FastMLTree.TopKey] = 0.0
        self.ModdedEvoDistance_D = self.modEvoDistance()
        self.TreeCoords_D = self.setTreeCoords()
        
        FurthestPosition = 0.0
        FurthestClade = ""
        
        #gets the furthest evolutionary distance
        for Key in self.FastMLTree.LeafKey_L:
            Val = self.TreeCoords_D[Key][0] + (12*len(Key))
            if Val > FurthestPosition:
                FurthestPosition = Val
                FurthestClade = Key
        
        self.BranchoInterest = ""
        
        for Key in self.FastMLTree.BranchKey_L:
            if Key.split(">>")[1] == self.DerivedoInterest:
                self.BranchoInterest = Key
        
        #gets all relevant information for the states portion of the figure
        self.StateIndices_L = [int(X)-1 for X in self.BranchToAlgorithm_D[self.BranchoInterest].getAllMutationXMLKeysAccordingToAccession(self.PDBoInterest)]
        self.LeafStates_D = {Key : [self.NodeToSeq_D[Key][state] for state in self.StateIndices_L] for Key in self.FastMLTree.LeafKey_L}
        self.StateColour_D = self.getStateToHex()
        
        self.StateInc = 25.0
        
        self.StateFigHEIGHT = 500
        self.StateFigWIDTH = self.StateInc * (len(self.StateIndices_L)) + 50
        self.StateFigXOffset = self.TreeFigXOffset+self.TreeFigWIDTH + (12*len(FurthestClade)) + 25
        self.StateFigYOffset = 50
        #creates the states and tree figure
        self.FigureSVG_D["TreeAndStates"].append(self.getSVGHeader(self.TreeFigHEIGHT+(self.TreeFigYOffset*2) , self.StateFigXOffset+self.StateFigWIDTH+self.TreeFigXOffset))
        self.makeTreeFig()
        self.makeStatesFig()
        self.FigureSVG_D["TreeAndStates"].append("</svg>")
        
        self.TreeAndStatesFOutPATH = self.OutputDirectory+"TreeAndStates.png"
        TreeStateFOut = open(self.TreeAndStatesFOutPATH , "w")
        cairosvg.svg2png(bytestring="\n".join(self.FigureSVG_D["TreeAndStates"]),write_to=TreeStateFOut)
        TreeStateFOut.close()
        
        
        LongestCladeName = ""
        for Key in self.FastMLTree.LeafKey_L:
            if len(Key) > len(LongestCladeName):
                LongestCladeName = Key
        
        #gets all relevant information for the alignment cartoon portion of the figure
        self.MatrixInfo = self.parseScoringMatrix()
        
        self.AlnInc = 11.0
        
        self.AlignmentFigWIDTH = self.AlnInc * len(self.MatrixInfo["Sseq"]) + self.AlnInc + (8*len(LongestCladeName))
        
        self.AlignmentFigHEIGHT = self.AlnInc * (len(self.FastMLTree.LeafKey_L) + 1) + self.AlnInc
        self.AlignmentFigXOffset = self.AlnInc
        self.AlignmentFigYOffset = self.AlnInc
        
        self.FigureSVG_D["Alignment"].append(self.getSVGHeader(self.AlignmentFigHEIGHT,self.AlignmentFigWIDTH))
        self.makeAlignmentFig()
        self.FigureSVG_D["Alignment"].append("</svg>")
        
        self.AlignmentFOutPATH = self.OutputDirectory+"Alignment.png"
        AlignmentFOut = open(self.AlignmentFOutPATH , "w")
        cairosvg.svg2png(bytestring="\n".join(self.FigureSVG_D["Alignment"]),write_to=AlignmentFOut)
        AlignmentFOut.close()
        
        #relevant information for the structure file in PDB format
        self.ColouredStructureFile = self.getColoredStructureFile()
        self.StructureFOutPATH = self.OutputDirectory+"Structure.pdb"
        open(self.StructureFOutPATH,"w").write(self.ColouredStructureFile.read())
        
        self.TotalFigWIDTH = 1000
        self.TotalFigHEIGHT = 600
        
        self.TotalElement_L = [self.getSVGHeader(self.TotalFigHEIGHT,self.TotalFigWIDTH)]
        self.TotalElement_L.append('''\t<image x="0" y="0" width="1000" height="500" xlink:href="file://%s"/>''' % (self.TreeAndStatesFOutPATH))
        self.TotalElement_L.append('''\t<image x="0" y="500" width="1000" height="100" xlink:href="file://%s"/>''' % (self.AlignmentFOutPATH))
        self.TotalElement_L.append("</svg>")
    
    "gets the header for any SVG format file"
    def getSVGHeader(self , FrameHEIGHT , FrameWIDTH):
        return """<?xml version="1.0" standalone="no"?>

<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"
"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">

<svg xmlns:xlink="http://www.w3.org/1999/xlink" xmlns='http://www.w3.org/2000/svg' version='1.1'
    width='%s' height='%s'>
""" % (str(FrameWIDTH) , str(FrameHEIGHT))
    
    "Dictionary where the key is the amino acid character and the value is the background colour"
    def getStateToHex(self):
        return {"A":"80B3E6","C":"E68080","D":"CC4DCC","E":"CC4DCC","F":"80B3E6",\
                "G":"E6994D","H":"1AB3B3","I":"80B3E6","K":"E6331A","L":"80B3E6",\
                "M":"80B3E6","N":"1ACC1A","P":"CCCC00","Q":"1ACC1A","R":"E6331A",\
                "S":"1ACC1A","T":"1ACC1A","V":"80B3E6","W":"80B3E6","Y":"1AB3B3",\
                "-":"FFFFFF","X":"FFFFFF"}
    
    "returns the total evolutionary distance from the origin to the node of interest"
    def getEvoDistance(self,startingToNodeKey):
        distance = 0.0
        rootNodeHasNotBeenReached = True
        ToNodeKey = startingToNodeKey
         
        while rootNodeHasNotBeenReached:
            
            distance += self.FastMLTree.BranchLength_D[ToNodeKey]
            
            branchUpHasNotBeenFound = True
            
            for BranchKey in self.FastMLTree.BranchKey_L:
                if branchUpHasNotBeenFound:
                    if re.compile(">>"+ToNodeKey+"$").search(BranchKey):
                    
                        branchUpHasNotBeenFound = False
                        ToNodeKey = BranchKey.split(">>")[0]
            
            if ToNodeKey == self.FastMLTree.TopKey:
                rootNodeHasNotBeenReached = False
        
        return distance
    
    "gets the node with the longest evolutionary distance from the origin"
    def getLongestEvoDistance(self):
        longestDistance = 0.0
        
        for LeafKey in self.FastMLTree.LeafKey_L:
            
            distance = self.getEvoDistance(LeafKey)
                       
            if distance > longestDistance:
                longestDistance = distance
            
        
        return longestDistance
    
    "modifies evolutionary distance into a different format"
    def modEvoDistance(self):
        Ret = {}
        
        for Key in self.EvoDistance_D.keys():
            if Key == self.FastMLTree.TopKey:
                Ret[Key] = self.EvoDistance_D[Key]
                
            else:
                if self.EvoDistance_D[Key] == 0:
                    Ret[Key] = self.EvoDistance_D[Key]
                else:
                    Ret[Key] = self.EvoDistance_D[Key]
        return Ret
    
    "sets tree node coordinates (horizontal and vertical) for the SVG image"
    def setTreeCoords(self):
        
        Lines_L = self.CogentTree.asciiArt().split("\n")
        MaxVert = 0
        VertCoord_D = {}
        
        for i in range(0,len(Lines_L)):
            
            if re.compile("[a-zA-Z0-9_\.@]+").search(Lines_L[i]):
                Leaves = re.findall("([a-zA-Z0-9_\.@]+)" , Lines_L[i])
                
                for Leaf in Leaves:
                    
                    VertCoord_D[Leaf] = i
                    MaxVert = i
        
        TreeCoords_D = {Key : [(self.ModdedEvoDistance_D[Key] / self.LongestDistance) * self.TreeFigWIDTH + self.TreeFigXOffset , float(float(VertCoord_D[Key]) / float(MaxVert)) * self.TreeFigHEIGHT + self.TreeFigYOffset] for Key in self.NodeToSeq_D.keys()}
        return TreeCoords_D
    
    "adds node names at each node vertex"
    def addNodeNamesAtNodePoints(self):
        for Key in self.FastMLTree.LeafKey_L:
            
            xy = self.TreeCoords_D[Key]
            xStart = str(xy[0])
            yStart = str(xy[1])
            self.FigureSVG_D["TreeAndStates"].append('''\t<text x='%s' y='%s' text-anchor='left' font-size='20' font-family='Courier' style="fill: #000000;"  >%s</text>''' % (xStart,yStart,Key))
    
    "adds the vertical lines of the tree image"
    def addVerticalLines(self):
        
        for branchKey in self.FastMLTree.BranchKey_L:
            fro = branchKey.split(">>")[0]
            to = branchKey.split(">>")[1]
            
            froXY = self.TreeCoords_D[fro]
            toXY = self.TreeCoords_D[to]
            
            if branchKey == self.BranchoInterest:
                self.FigureSVG_D["TreeAndStates"].append('''\t<line class='axis' x1='%s' y1='%s' x2='%s' y2='%s' style="stroke:rgb(255,0,0);stroke-width:1 " />''' % (str(froXY[0]) , str(froXY[1]) , str(froXY[0]) , str(toXY[1])))
            else:
                self.FigureSVG_D["TreeAndStates"].append('''\t<line class='axis' x1='%s' y1='%s' x2='%s' y2='%s' style="stroke:rgb(0,0,0);stroke-width:1 " />''' % (str(froXY[0]) , str(froXY[1]) , str(froXY[0]) , str(toXY[1])))
    
    "adds the horizontal lines of the tree image"
    def addHorizontalLines(self):
        
        for branchKey in self.FastMLTree.BranchKey_L:
            
            fro = branchKey.split(">>")[0]
            to = branchKey.split(">>")[1]
            
            froXY = self.TreeCoords_D[fro]
            toXY = self.TreeCoords_D[to]
            
            if branchKey == self.BranchoInterest:
                self.FigureSVG_D["TreeAndStates"].append('''\t<line class='axis' x1='%s' y1='%s' x2='%s' y2='%s' style="stroke:rgb(255,0,0);stroke-width:1 " />''' % (str(froXY[0]) , str(toXY[1]) , str(toXY[0]) , str(toXY[1])))
            else:
                self.FigureSVG_D["TreeAndStates"].append('''\t<line class='axis' x1='%s' y1='%s' x2='%s' y2='%s' style="stroke:rgb(0,0,0);stroke-width:1 " />''' % (str(froXY[0]) , str(toXY[1]) , str(toXY[0]) , str(toXY[1])))
            
    "does all methods necessary to make the tree image"
    def makeTreeFig(self):
        self.addNodeNamesAtNodePoints()
        self.addVerticalLines()
        self.addHorizontalLines()
              
    "adds the rows for the mutated states in each sequence"
    def addStateRows(self):
        inc =  self.StateInc
        vertInc = float(self.StateFigHEIGHT / float(len(self.LeafStates_D)))
        
        lowestY = float("inf")
        
        for Key in self.TreeCoords_D.keys():
            if self.TreeCoords_D[Key][1] < lowestY:
                lowestY = self.TreeCoords_D[Key][1]
        
        stateY = lowestY - (1.5*vertInc)
        
        stateX = 0.0 + self.StateFigXOffset
        for i in self.StateIndices_L:
            
            self.FigureSVG_D["TreeAndStates"].append('''\t<text x='%s' y='%s' text-anchor='middle' font-size='16' font-family='Courier' transform="rotate(90, %s, %s)" style="fill: #000000;"  >%s</text>''' % (str(stateX),str(stateY),str(stateX),str(stateY),str(i+1)))
            
            stateX += inc
            
        
        for Key in self.LeafStates_D.keys():
            X = 0.0 + self.StateFigXOffset
            
            for State in self.LeafStates_D[Key]:
                Y = self.TreeCoords_D[Key][1]
                
                
                RectX = X - (float(inc/2.0)) 
                RectY = Y - (float(vertInc/2.0)) - 5.0
                
                self.FigureSVG_D["TreeAndStates"].append('''\t<rect class='r%s' x='%s' y='%s' width='%s' height='%s' style="fill:#%s" />''' % (str(self.RectCount),\
                                                                                                                                 str(RectX),str(RectY),\
                                                                                                                                 str(inc),str(vertInc),\
                                                                                                                                 self.StateColour_D[State]))
                self.FigureSVG_D["TreeAndStates"].append('''\t<text x='%s' y='%s' font-size='20' font-family='Courier' text-anchor='middle' style="fill: #000000;"  >%s</text>''' % (str(X),str(Y),State))
                
                X += inc
        
    "executes the method to make the states figure"
    def makeStatesFig(self):
        self.addStateRows()
    
    "parses the scoring matrix for alignment to the PDB sequence information"
    def parseScoringMatrix(self):
        allAlignments_L = re.findall("<PDB_alignment>.+?</PDB_alignment>",open(self.MatrixPATH,"r").read() , re.DOTALL)
        KeyAln = ""
        NotFound = True
        
        for Alignment in allAlignments_L:
            if NotFound:
                PDBID = re.compile("<PDB_id>(.+?)</PDB_id>").search(Alignment).group(1).split("|")[0]
                if self.PDBoInterest.upper() == PDBID:
                    NotFound = False
                    KeyAln = Alignment
                    self.ChainoInterest = re.compile("<PDB_id>(.+?)</PDB_id>").search(Alignment).group(1).split("|")[1].lower()
        
        return {"Qstart" : int(re.compile("<Alignment_start_query>(.+?)</Alignment_start_query>").search(KeyAln).group(1))-1,\
                "Qend" : int(re.compile("<Alignment_end_query>(.+?)</Alignment_end_query>").search(KeyAln).group(1))-1,\
                "Sstart" : int(re.compile("<Alignment_start_subject>(.+?)</Alignment_start_subject>").search(KeyAln).group(1))-1,\
                "Send" : int(re.compile("<Alignment_end_subject>(.+?)</Alignment_end_subject>").search(KeyAln).group(1))-1,\
                "Sseq" : re.compile("<Aligned_subject_sequence>(.+?)</Aligned_subject_sequence>").search(KeyAln).group(1)}
        
    "makes the cartoon of all aligned sequences in the protein family"
    def makeAlignmentFig(self):
        AllSeqs_L = [self.MatrixInfo["Sseq"]] + [self.NodeToSeq_D[Key][self.MatrixInfo["Qstart"] : self.MatrixInfo["Qstart"]+len(self.MatrixInfo["Sseq"])] for Key in self.FastMLTree.LeafKey_L]
        l1 = len(AllSeqs_L[0])
        AllHeaders_L = [self.PDBoInterest] + self.FastMLTree.LeafKey_L
        l2 = 0
        
        for Header in AllHeaders_L:
            if len(Header) > l2:
                l2 = len(Header)
        
        l = l1
        
        xinc = self.AlnInc  
        yinc = self.AlnInc
        
        Y = self.AlignmentFigYOffset
        
        for i in range(0,len(AllSeqs_L)):
            
            X = 0.0 + self.AlignmentFigXOffset
            
            for State in AllSeqs_L[i]:
                
                RectX = X - (float(xinc/2.0)) 
                RectY = Y - (float(yinc/2.0)) - 5.0
                
                self.FigureSVG_D["Alignment"].append('''\t<rect class='r%s' x='%s' y='%s' width='%s' height='%s' style="fill:#%s" />''' % (str(self.RectCount),\
                                                                                                                                 str(RectX),str(RectY),\
                                                                                                                                 str(xinc),str(yinc),\
                                                                                                                                 self.StateColour_D[State]))
                self.FigureSVG_D["Alignment"].append('''\t<text x='%s' y='%s' text-anchor='middle' font-size='10' font-family='Courier' style="fill: #000000;"  >%s</text>''' % (str(X),str(Y),State))
                
                X += xinc
            
            self.FigureSVG_D["Alignment"].append('''\t<text x='%s' y='%s' text-anchor='left' font-size='10' font-family='Courier' style="fill: #000000;"  >%s</text>''' % (str(X+self.AlnInc),str(Y),AllHeaders_L[i]))
            
            Y += yinc    
        
    "gets a PDB format file with the temperature factors coloured to reflect mutated sites"
    def getColoredStructureFile(self):
        NotFound = True
        DesiredBranchKey = ""
        for BranchKey in self.FastMLTree.BranchKey_L:
            if BranchKey.split(">>")[1] == self.DerivedoInterest:
                DesiredBranchKey = BranchKey
                NotFound = False
        
        PDBAndPDBXMLContents = getAllPDBFileDicts([self.PDBoInterest])
        SA = self.BranchToAlgorithm_D[DesiredBranchKey]
        SA.PDBContents_D = PDBAndPDBXMLContents[0]
        SA.PDBXMLContents_D = PDBAndPDBXMLContents[1]
        
        FH = getOutputTempFile()
        SA.createPDBColoredFile(self.PDBoInterest,FH.name)
        
        return FH
    

    
    