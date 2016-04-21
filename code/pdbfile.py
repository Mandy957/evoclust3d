####################################################################################################
#                                                                                                  #
#           PROJECT Protein Adaptation                                                             #
#           CLASS PDBFile                                                                          #
#           PROGRAMMER Jeremy Adams                                                                #
#           STARTED 25-07-13                                                                       #
#           LASTMOD 08-08-14                                                                       #
#                                                                                                  #
#           DESCRIPTION Class that parses an entire PDB structure and provides information about   #
#                       the chains, residue positions, and ASA                                     #
#                                                                                                  #
####################################################################################################

import os
import re
import numpy
import math

from pdbatom import PDBAtom
from pdbresiduetypes import *

from staticmethods import *


class PDBFile:
    
    "CONSTRUCTOR"
    def __init__(self,filename):
        
        """
        Class attributes:
        Accession (String): path to the pdb file
        AllContent (String): string containing everything in the PDB file
        Atoms (List): all PDBAtom objects for all atoms in the PDB file in the order they are listed (alternate position and water molecules excluded)
        ChainGroupedAtoms (List): list within list: each sub list contains all PDBAtom objects for a particular chain
        ChainGroupedResidues (List): list within list: each sub list contains all PDBResidue objects which contain all PDBAtom objects within a particular chain
        """
        
        self.Features_L = ["HELIX","SHEET","SSBOND","LINK"]
        self.Features_D = {key : [] for key in self.Features_L}
        self.FeatureResidue_D = {key : [] for key in self.Features_L}
        
        self.Accession = filename
        self.AllContent = open(self.Accession,'rt').read()
        
        self.completeFeatures_D()
        self.completeFeatureResidue_D()
        
        #get all PDBAtom objects from the file
        self.Atoms = []
        self.Atoms = self.readAllAtoms()
        
        self.runGlobalSAS() #perform SAS analysis on the entire protein structure
        self.ChainGroupedAtoms = self.groupAtomsByChain() #place atoms within their respective chains
        self.runLocalSAS() #perform SAS analysis on each individual chain
        self.ChainGroupedResidues = self.groupAtomsByResidue() #place atoms within their respective residues
        
        
        
    def completeFeatures_D(self):
        for line in self.AllContent.split("\n"):
            ls = line.split()
            if len(ls) > 0:
                if ls[0] in set(self.Features_L):
                    self.Features_D[ls[0]].append(line)
    def completeFeatureResidue_D(self):
        self.completeHelix()
        self.completeSheet()
        self.completeSSBond()
        self.completeLink()
    
    def completeHelix(self):
        for line in self.Features_D["HELIX"]:
            chainStart = line[19].replace(" ","")
            resStart = int(line[21:25].replace(" ",""))
            chainEnd = line[31].replace(" ","")
            resEnd = int(line[33:37].replace(" ",""))
            
            for i in range(resStart,resEnd+1):
                self.FeatureResidue_D["HELIX"].append("%s%s" % (chainStart,str(i)))
            
    def completeSheet(self):
        for line in self.Features_D["SHEET"]:
            chainStart = line[21].replace(" ","")
            resStart = int(line[22:26].replace(" ",""))
            chainEnd = line[32].replace(" ","")
            resEnd = int(line[33:37].replace(" ",""))
            
            for i in range(resStart,resEnd+1):
                self.FeatureResidue_D["SHEET"].append("%s%s" % (chainStart,str(i)))
        
    def completeSSBond(self):
        for line in self.Features_D["SSBOND"]:
            chainStart = line[15]
            resStart = int(line[17:21].replace(" ",""))
            chainEnd = line[29].replace(" ","")
            resEnd = int(line[31:35].replace(" ",""))
            
            self.FeatureResidue_D["SSBOND"].append(chainStart+str(resStart))
            self.FeatureResidue_D["SSBOND"].append(chainEnd+str(resEnd))
            
    def completeLink(self):
        for line in self.Features_D["LINK"]:
            chainStart = line[21]
            resStart = int(line[22:26].replace(" ",""))
            chainEnd = line[51].replace(" ","")
            resEnd = int(line[52:56].replace(" ",""))
            
            self.FeatureResidue_D["LINK"].append(chainStart+str(resStart))
            self.FeatureResidue_D["LINK"].append(chainEnd+str(resEnd))
            
    "Returns an instance of PDBAtom with parameters corresponding to flat file info"
    def readAtom(self,line):
        return PDBAtom(line)
                       
    "Returns array of PDBAtoms corresponding to all atoms in flat file"
    def readAllAtoms(self):
        
        nmr_pattern = re.compile("MODEL(.+?)ENDMDL" , re.DOTALL) #if it is an NMR structure, only take the first model
        
        if nmr_pattern.search(self.AllContent):
            self.AllContent = nmr_pattern.search(self.AllContent).group(1)
        
        #gets all lines, then only atoms, then only non-water molecule atoms, then only atoms without an alternate position
        all_lines = re.split("\n" , self.AllContent)
        all_lines = [line for line in all_lines if line.startswith("ATOM")]
        all_lines = [line for line in all_lines if line[77] != "H"]
        all_lines = [line for line in all_lines if line[16] == " "]
        
        return [self.readAtom(line) for line in all_lines]
    
    "gets the SAS of an atom from the POPS report"
    def getAtomSAS(self,pops_line):
        return pops_line[31:37].replace(" ","")
    
    "gets the SAS of each atom in the entire protein structure"
    def runGlobalSAS(self):
        
        #makes temp input and output files
        all_atom_lines = "\n".join([atom.Info["original_line"] for atom in self.Atoms])
        
        IFH = getInputTempFile(all_atom_lines)
        OFH = getOutputTempFile()
        
        #executes pops
        os.system("pops --pdb " + IFH.name + ' --atomOut --noHeaderOut --noTotalOut --popsOut ' + OFH.name + " > /dev/null")
            
        #reads pops output and ascribes each SAS to the atom it belongs to
        pops_lines = OFH.read().split("\n")[0:-1]
        pops_completion = True
        if len(pops_lines) == 0:
            pops_completion = False
        
        for i in range(0,len(self.Atoms)):
            
            if pops_completion:
                self.Atoms[i].GlobalSAS = self.getAtomSAS(pops_lines[i])
            else:
                self.Atoms[i].GlobalSAS = "0.0"
    
    "gets lists of atom objects grouped by chain identifier"
    def groupAtomsByChain(self):
        
        all_chains = []
        one_chain = []
        check = self.Atoms[0].Info['chain']
        #if atom has the same chain id as the previous, then the chain grows, otherwise the chain is complete and another chain starts
        for atom in self.Atoms:
            
            if atom.Info['chain'] == check:
                one_chain.append(atom)
                
            else:
                all_chains.append(one_chain)
                
                one_chain = []
                one_chain.append(atom)
                
                check = atom.Info['chain']
        
        all_chains.append(one_chain)
        
        return all_chains
    
    "gets the SAS of each atom when the chain exists by itself"
    def runLocalSAS(self):
        
        for chain in self.ChainGroupedAtoms:
            all_atoms = []
            
            #gets input and output temp file for each chain
            for atom in chain:
                all_atoms.append(atom.Info['original_line'])
            
            IFH = getInputTempFile('\n'.join(all_atoms))
            OFH = getOutputTempFile()
            
            #executes pops
            os.system("pops --pdb " + IFH.name + ' --atomOut --noHeaderOut --noTotalOut --popsOut ' + OFH.name + " > /dev/null")
            
            #reads pops and ascribes sas value to the proper atom
            pops_lines = OFH.read().split("\n")[0:-1]
            pops_completion = True
            if len(pops_lines) == 0:
                pops_completion = False
            
            for a in range(0,len(chain)):
                if pops_completion:
                    chain[a].LocalSAS = self.getAtomSAS(pops_lines[a])
                else:
                    chain[a].LocalSAS = '0.0'
    
    "Contructs PDBResidue objects for each grouping of atoms that belong to the same residue"
    def groupAtomsByResidue(self):
        
        all_chains = []
        one_residue = []
        #for each chain
        for chain in self.ChainGroupedAtoms:
            #checks if the next atom has the same residue number as previous, if it does, they are part of the same residue
            all_residues_per_chain = []
            check = int(chain[0].Info['residue_num'])
            
            chain_id = chain[0].Info['chain']
            
            for atom in chain:
                
                challenger = int(atom.Info['residue_num'])
                #if same residue number, it is part of the same residue
                if challenger == check:
                    one_residue.append(atom)
                #if the new residue number is just one more than the last, it is a new residue with no residues in between
                elif challenger == check + 1:
                    
                    
                    all_residues_per_chain.append(self.getSpecifiedResidue(one_residue))
                    
                    one_residue = []
                    one_residue.append(atom)
                    check = challenger
                
                #if the new residue number is more than one than the last, there are fictional residues that need to be added in
                else:
                    
                    all_residues_per_chain.append(self.getSpecifiedResidue(one_residue))
                    one_residue = []
                    one_residue.append(atom)
                    
                    #one fictional residue for every gap number between the current real residue and the last real residue
                    for i in range(check+1,challenger):
                        
                        all_residues_per_chain.append(PDBFictionalResidue([PDBAtom("ATOM      1  N   XXX %s" % (chain_id)+" "*(4-len(str(i)))+str(i)+"      00.000  00.000  00.000  1.00 0.000           N  "),\
                                                                           PDBAtom("ATOM      1  C   XXX %s" % (chain_id)+" "*(4-len(str(i)))+str(i)+"      00.000  00.000  00.000  1.00 0.000           C  "),\
                                                                           PDBAtom("ATOM      1  CA  XXX %s" % (chain_id)+" "*(4-len(str(i)))+str(i)+"      00.000  00.000  00.000  1.00 0.000           C  "),\
                                                                           PDBAtom("ATOM      1  O   XXX %s" % (chain_id)+" "*(4-len(str(i)))+str(i)+"      00.000  00.000  00.000  1.00 0.000           O  ")]))
                    check = challenger
                
            all_residues_per_chain.append(self.getSpecifiedResidue(one_residue))
            one_residue = []
        
            all_chains.append(all_residues_per_chain)
             
        return all_chains
    
    "gets the proper residue object that reflects the residue type of the atoms grouped together"
    def getSpecifiedResidue(self,one_residue):
        
        resi = one_residue[0].Info['residue_name']
        chain = one_residue[0].Info['chain']
        num = one_residue[0].Info['residue_num']
        
        #default is a fictional residue if the observed residue is not one of the standard types
        ret = PDBFictionalResidue([PDBAtom("ATOM      1  N   XXX %s" % (chain)+" "*(4-len(num))+num+"      00.000  00.000  00.000  1.00 0.000           N  "),\
                                   PDBAtom("ATOM      1  C   XXX %s" % (chain)+" "*(4-len(num))+num+"      00.000  00.000  00.000  1.00 0.000           C  "),\
                                   PDBAtom("ATOM      1  CA  XXX %s" % (chain)+" "*(4-len(num))+num+"      00.000  00.000  00.000  1.00 0.000           C  "),\
                                   PDBAtom("ATOM      1  O   XXX %s" % (chain)+" "*(4-len(num))+num+"      00.000  00.000  00.000  1.00 0.000           O  ")])
        
        #checks atom residue name and gets the appropriate object
        if resi == 'ALA':
            ret = PDBAlanine(one_residue)
        elif resi == 'ARG':
            ret = PDBArginine(one_residue)
        elif resi == 'ASN':
            ret = PDBAsparagine(one_residue)
        elif resi == 'ASP':
            ret = PDBAsparticAcid(one_residue)
        elif resi == 'CYS':
            ret = PDBCysteine(one_residue)
        elif resi == 'GLU':
            ret = PDBGlutamicAcid(one_residue)
        elif resi == 'GLN':
            ret = PDBGlutamine(one_residue)
        elif resi == 'GLY':
            ret = PDBGlycine(one_residue)
        elif resi == 'HIS':
            ret = PDBHistidine(one_residue)
        elif resi == 'ILE':
            ret = PDBIsoleucine(one_residue)
        elif resi == 'LEU':
            ret = PDBLeucine(one_residue)
        elif resi == 'LYS':
            ret = PDBLysine(one_residue)
        elif resi == 'MET':
            ret = PDBMethionine(one_residue)
        elif resi == 'PHE':
            ret = PDBPhenylalanine(one_residue)
        elif resi == 'PRO':
            ret = PDBProline(one_residue)
        elif resi == 'SER':
            ret = PDBSerine(one_residue)
        elif resi == 'THR':
            ret = PDBThreonine(one_residue)
        elif resi == 'TRP':
            ret = PDBTryptophan(one_residue)
        elif resi == 'TYR':
            ret = PDBTyrosine(one_residue)
        elif resi == 'VAL':
            ret = PDBValine(one_residue)
        
        return ret
    
    "returns the primary sequence of a single chain"
    def getPrimarySequence(self,chain):
        
        #joins all the one letter symbols together
        sequence = []
        
        for residue in chain:
            sequence.append(residue.Info['one_letter']) 
        return ''.join(sequence)
    
    def vecMagnitude(self,vec):
        return math.sqrt(math.pow(vec[0],2)+\
                         math.pow(vec[1],2)+\
                         math.pow(vec[2],2))
    
    def getAngle(self,a,b):
        ret = "NA"
        costheta = numpy.dot(a,b)
        arccos = numpy.arccos(costheta)
        deg = numpy.rad2deg(arccos)
        
        if numpy.isnan(deg):
            pass
        else:
            ret = deg
        
        return ret
        
   
    def getDihedral(self,a,b,c,d):
        R ="NA"
        a = makePoint(a.split(","))
        b = makePoint(b.split(","))
        c = makePoint(c.split(","))
        d = makePoint(d.split(","))
        u = compute3dVector(a,b)
        v = compute3dVector(b,c)
        w = compute3dVector(c,d)
        uxv = numpy.cross(u,v)
        vxw = numpy.cross(v,w)
        
        uxvmag = self.vecMagnitude(uxv)
        vxwmag = self.vecMagnitude(vxw)
        
        uxvnorm = (uxv[0]/uxvmag , uxv[1]/uxvmag , uxv[2]/uxvmag)
        vxwnorm = (vxw[0]/vxwmag , vxw[1]/vxwmag , vxw[2]/vxwmag)
        
        sign2 = compute3dVector(c,d)
        signdot = numpy.dot(uxvnorm,sign2)
        angle = self.getAngle(uxvnorm,vxwnorm)
        
        if angle == "NA":
            pass
        else:
            
            if signdot < 0.0:
                angle = angle * -1.0
            R = angle
        return R
        
    "prepares the PDBXML document string"
    def xmlWrite(self):
        qq = [] #list of strings
        
        #basic specs
        qq.append("<PDBFile>\n")
        qq.append("\t<Accession>"+self.Accession+"</Accession>\n")
        
        #gets the primary sequence information for each chain
        for chain in self.ChainGroupedResidues:
            qq.append("\t<Chain>\n")
            qq.append("\t\t<Chain_name>"+chain[0].Atoms[0].Info['chain']+"</Chain_name>\n")
            qq.append("\t\t<Primary_sequence>"+\
                      "<Start>"+chain[0].Identity['residue_number']+"</Start>"+\
                      "<seq>"+self.getPrimarySequence(chain)+"</seq>"+\
                      "<End>"+chain[-1].Identity['residue_number']+"</End>"+\
                      "</Primary_sequence>\n")
            qq.append("\t</Chain>\n")
        qq.append("\t<Residues>\n")
        
        #more specific residue information (SAS,BFactor,Position)
        for chain in self.ChainGroupedResidues:
            for i in range(0,len(chain)):
                residue = chain[i]

                #Residue Type and Chain Position
                oneRes = []
                oneRes.append("\t\t<R><t>%s</t><n>%s%s</n>" % (residue.Info['one_letter'],\
                                                              residue.Identity['member_of_chain'],\
                                                              residue.Identity['residue_number']))
                
                #Residue Solvent Accessibility
                if residue.SideChainIsInitialized:
                    oneRes.append("<s>%s;%s</s>" % (residue.AbsoluteSideChainGlobalSAS,\
                                                    residue.AbsoluteSideChainLocalSAS))
                else:
                    oneRes.append("<s>0.0;0.0</s>")
                
                oneRes.append("<p>")
                
                #Residue Average Positions
                if residue.SideChainIsInitialized and residue.BackboneIsInitialized:
                    oneRes.append("%s;" % residue.Positions['centroid'])
                else:
                    oneRes.append("NA,NA,NA;")
                
                if residue.BackboneIsInitialized:
                    oneRes.append("%s;%s;%s;" % (residue.Positions['backbone']['c_alpha'],\
                                                 residue.Positions['backbone']['n_amino'],\
                                                 residue.Positions['backbone']['c_carboxyl']))
                else:
                    oneRes.append("NA,NA,NA;NA,NA,NA;NA,NA,NA;")
                
                if residue.SideChainIsInitialized:
                    oneRes.append("%s;%s" % (residue.Positions['side_chain']['scs'],\
                                             residue.Positions['side_chain']['sce']))
                else:
                    oneRes.append("NA,NA,NA;NA,NA,NA")
                
                oneRes.append("</p>")
                
                #B-Factor
                
                oneRes.append("<B>")
                #BACKBONE ONLY
                if residue.Info['one_letter'] == "X":
                    oneRes.append("NA")
                else:
                    if residue.BackboneIsInitialized:
                        AVGBFactor = residue.getAverageBackboneBFactor()
                        oneRes.append(str(round(AVGBFactor,2)))
                    else:
                        oneRes.append("NA")
                
                #SIDE CHAIN ONLY
                if residue.Info['one_letter'] == "X":
                    oneRes.append(";NA")
                else:
                    if residue.SideChainIsInitialized:
                        AVGBFactor = residue.getAverageSideChainBFactor()
                        oneRes.append(";"+str(round(AVGBFactor,2)))
                    else:
                        oneRes.append(";NA")
                
                #OVERALL
                if residue.Info['one_letter'] == "X":
                    oneRes.append(";NA")
                else:
                    if residue.BackboneIsInitialized and residue.SideChainIsInitialized:
                        AVGBFactor = residue.getAverageTotalBFactor()
                        oneRes.append(";"+str(round(AVGBFactor,2)))
                    else:
                        oneRes.append(";NA")
                    
                oneRes.append("</B>")
                
                #Secondary structure and link feature information
                
                oneRes.append("<F>")
                Features = []
                for F in self.Features_L:
                    if residue.Identity['member_of_chain']+residue.Identity['residue_number'] in set(self.FeatureResidue_D[F]):
                        Features.append(F)
                
                if len(Features) == 0:
                    Features = "None"
                else:
                    Features = "-".join(Features)
                oneRes.append(Features)
                
                
                oneRes.append("</F>")
                
                #dihedral angle
                
                oneRes.append("<A>")
                
                #phi dihedral angle
                phi = "NA"
                if i == 0:
                    pass
                else:
                    previousRes = chain[i-1]
                    if previousRes.Info['one_letter'] != "X" and residue.Info['one_letter'] != "X":
                        if previousRes.BackboneIsInitialized and residue.BackboneIsInitialized:
                        
                            a = previousRes.Positions['backbone']['c_carboxyl']
                            b = residue.Positions['backbone']['n_amino']
                            c = residue.Positions['backbone']['c_alpha']
                            d = residue.Positions['backbone']['c_carboxyl']
                            phi = str(round(self.getDihedral(a,b,c,d),4))
                    
                    
                
                
                #psi dihedral angle
                psi = "NA"
                if i == len(chain)-1:
                    pass
                else:
                    nextRes = chain[i+1]
                    if nextRes.Info['one_letter'] != "X" and residue.Info['one_letter'] != "X":
                        if nextRes.BackboneIsInitialized and residue.BackboneIsInitialized:
                        
                            a = residue.Positions['backbone']['n_amino']
                            b = residue.Positions['backbone']['c_alpha']
                            c = residue.Positions['backbone']['c_carboxyl']
                            d = nextRes.Positions['backbone']['n_amino']
                            psi = str(round(self.getDihedral(a,b,c,d),4))
                            
                oneRes.append(";".join([phi,psi]))
                oneRes.append("</A>")
                
                #Backbone vector
                
                oneRes.append("<V>")
                Vector = "NA,NA,NA"
                
                if residue.BackboneIsInitialized:
                    
                    C = makePoint(residue.Positions['backbone']['c_carboxyl'].split(","))
                    N = makePoint(residue.Positions['backbone']['n_amino'].split(","))
                    Vector = ",".join([str(computeVector(N[0],C[0])) , str(computeVector(N[1],C[1])) , str(computeVector(N[2],C[2]))])
                oneRes.append(Vector)
                oneRes.append("</V>")
                
                oneRes.append("</R>\n")
                
                qq.append(''.join(oneRes))
                
        qq.append("\t</Residues>\n")
        qq.append("</PDBFile>\n")
        return ''.join(qq)
        
    
    