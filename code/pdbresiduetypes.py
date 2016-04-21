####################################################################################################
#                                                                                                  #
#           PROJECT Protein Adaptation                                                             #
#           FILE pdbresiduetypes.py                                                                #
#           PROGRAMMER Jeremy Adams                                                                #
#           STARTED 25-07-13                                                                       #
#           LASTMOD 08-08-14                                                                       #
#                                                                                                  #
#           DESCRIPTION File containing multiple classes outlining specific amino acid types,      #
#                       their component Atoms and structural characteristics as given in the PDB   #
#                                                                                                  #
####################################################################################################

from pdbatom import PDBAtom
import sys
import re
import abc
import math
import numpy

####################################################################################################
#                                                                                                  #
#           CLASS PDBResidue                                                                       #
#                                                                                                  #
#           DESCRIPTION General class encompassing the shared backbone Atoms of all amino acids    #
#                       Contains methods for computing vector points and side chain ASA            #
#                                                                                                  #
####################################################################################################


class PDBResidue:
    
    "CONSTRUCTOR"
    def __init__(self,AtomList):
        
        """
        Class attributes:
        Atoms (List): contains original listing of PDBAtom objects for that particular residue
        Identity (Dict): contains information to chain membership and residue
        SpecificAtoms (Dict): key is the name of a specific atom name/position within the residue, value is PDBAtom object corresponding to that atom
        Positions (Dict): key is the name/title of a 3D position, value is the coordinates of that 3D position represented as a string
        Info(Dict): contains information on amino acid type's one letter and three letter codes
        
        BackboneIsInitialized (Bool): True if all the necessary atoms of the backbone are present in the PDBResidue instance
        SideChainIsInitialized (Bool): True if all necessary atoms of the amino acid type's side chain are present in the PDBResidue instance
        AbsoluteSideChainGlobalSAS (float) : Sum of the SAS of all atoms in the side chain when the entire structure is present
        AbsoluteSideChainLocalSAS (float) : Sum of the SAS of all atoms in the side chain when only the particular chain is present
        """
        
        self.Atoms = []
        self.Atoms = AtomList
        
        self.Identity = {}
        self.Identity['member_of_chain'] = self.Atoms[0].Info['chain']
        self.Identity['residue_number'] = self.Atoms[0].Info['residue_num']
        
        #initializes the backbone of any PDBResidue
        self.SpecificAtoms = {}
        self.SpecificAtoms['backbone'] = {}
        self.SpecificAtoms['backbone']['c_alpha'] = self.getSpecificAtom(self.Atoms,'CA')
        self.SpecificAtoms['backbone']['n_amino'] = self.getSpecificAtom(self.Atoms,'N')
        self.SpecificAtoms['backbone']['c_carboxyl'] = self.getSpecificAtom(self.Atoms,'C')
        self.SpecificAtoms['backbone']['o_carboxyl'] = self.getSpecificAtom(self.Atoms,'O')
        
        self.Positions = {}
        self.Positions['backbone'] = {}
        
        #checks if all backbone atoms are present
        self.BackboneIsInitialized = self.checkIfProperlyInitialized(self.SpecificAtoms['backbone'])
        
        #gets the positions of key points in the backbone if all atoms are present
        if self.BackboneIsInitialized:
            self.Positions['backbone']['c_alpha'] = self.SpecificAtoms['backbone']['c_alpha'].getPointAsString()
            self.Positions['backbone']['n_amino'] = self.SpecificAtoms['backbone']['n_amino'].getPointAsString()
            self.Positions['backbone']['c_carboxyl'] = self.SpecificAtoms['backbone']['c_carboxyl'].getPointAsString()
    
    "gets a desired atom according to the atom's unique name"
    def getSpecificAtom(self, aas , desired_name):
    
        ret_atom = PDBAtom("ATOM      1  XXX XXX X   X      0.0000   0.000   0.000  1.00 0.000           X  ")
        #returns the atom with the specified name
        for atom in aas:
            if atom.Info['name'] == desired_name:
                ret_atom = atom
                
        return ret_atom
    
    "checks if either the backbone or the side chain has all real residues in it"
    def checkIfProperlyInitialized(self,specific_atom_hash):
        
        ret = True
        #if any atom has "Error" for its name, then the atom set is not initialized
        for atom_key in specific_atom_hash.keys():
            if specific_atom_hash[atom_key].Info['name'] == 'XXX':
                ret = False
        
        return ret
    
    "totals the side chain SAS for all atoms in the side chain for SAS values when the whole protein is present"
    def getAbsoluteSideChainGlobalSAS(self):
        sas = 0.0
        
        if self.BackboneIsInitialized and self.SideChainIsInitialized:
            for key in self.SpecificAtoms['side_chain'].keys():
                sas += float(self.SpecificAtoms['side_chain'][key].GlobalSAS)
                
        return sas
    
    "totals the side chain SAS for all atoms in the side chain for SAS values when only that chain is present"
    def getAbsoluteSideChainLocalSAS(self):
        sas = 0.0
        
        if self.BackboneIsInitialized and self.SideChainIsInitialized:
            for key in self.SpecificAtoms['side_chain'].keys():
                sas += float(self.SpecificAtoms['side_chain'][key].LocalSAS)
            
        return sas
    
    "gets the average position of a list of coordinates in a single dimension"
    def compute1Dcentroid(self,number_L):
        theSum = 0.0
        
        for number in number_L:
            theSum += number
        
        return round(float(float(theSum) / float(len(number_L))) , 3)
    
    "gets the average 3D position of all atoms within the residue"
    def compute3Dcentroid(self):
        
        x_L = []
        y_L = []
        z_L = []
        
        #gets lists of 1D coordinates for X,Y,Z dimensions
        for firstKey in self.SpecificAtoms.keys():
            for secondKey in self.SpecificAtoms[firstKey].keys():
                point = self.SpecificAtoms[firstKey][secondKey].getPointAsFloat()
                x_L.append(point[0])
                y_L.append(point[1])
                z_L.append(point[2])
        
        #3D centroid is the collection of 1D centroids
        return "%s,%s,%s" % (str(self.compute1Dcentroid(x_L)),\
                             str(self.compute1Dcentroid(y_L)),\
                             str(self.compute1Dcentroid(z_L)))
    
    "finds the average position of two atoms"
    def computeCentroidOf2Atoms(self, atom_a , atom_b):
        pointa = atom_a.getPointAsFloat()
        pointb = atom_b.getPointAsFloat()
        
        x_L = [pointa[0] , pointb[0]]
        y_L = [pointa[1] , pointb[1]]
        z_L = [pointa[2] , pointb[2]]
        
        return "%s,%s,%s" % (str(self.compute1Dcentroid(x_L)),\
                             str(self.compute1Dcentroid(y_L)),\
                             str(self.compute1Dcentroid(z_L)))
    
    "gets the average BFactor of all atoms in the backbone"
    def getAverageBackboneBFactor(self):
        BB_L = [self.SpecificAtoms['backbone'][Key].Info['old_temp_factor'] for Key in self.SpecificAtoms['backbone'].keys()]
        return numpy.mean(BB_L)
    
    "gets the average BFactor of all atoms in the sidechain"
    def getAverageSideChainBFactor(self):
        SC_L = [self.SpecificAtoms['side_chain'][Key].Info['old_temp_factor'] for Key in self.SpecificAtoms['side_chain'].keys()]
        return numpy.mean(SC_L)
    
    "gets the average BFactor of all atoms"
    def getAverageTotalBFactor(self):
        BB_L = [self.SpecificAtoms['backbone'][Key].Info['old_temp_factor'] for Key in self.SpecificAtoms['backbone'].keys()]
        SC_L = [self.SpecificAtoms['side_chain'][Key].Info['old_temp_factor'] for Key in self.SpecificAtoms['side_chain'].keys()]
        Total_L = BB_L + SC_L
        return numpy.mean(Total_L)
    
####################################################################################################
#                                                                                                  #
#           CLASS PDBAlanine                                                                       #
#                                                                                                  #
#           DESCRIPTION specific class inheriting from PDBResidue with unique side chain atom      #
#                       identities pertaining to amino acid type alanine                           #
#                                                                                                  #
####################################################################################################

class PDBAlanine(PDBResidue):
    
    "CONSTRUCTOR"
    def __init__(self,atom_array):
        PDBResidue.__init__(self,atom_array)
        
        self.Info = {}
        self.Info['one_letter'] = 'A'
        self.Info['three_letter'] = 'ALA'
        
        #initialize the side chain
        self.SpecificAtoms['side_chain'] = {}
        self.SpecificAtoms['side_chain']['c_beta'] = self.getSpecificAtom(self.Atoms,'CB')
        
        #check if side chain is initialized
        self.SideChainIsInitialized = self.checkIfProperlyInitialized(self.SpecificAtoms['side_chain'])
        
        #compute centroid
        if self.SideChainIsInitialized and self.BackboneIsInitialized:
            self.Positions['centroid'] = self.compute3Dcentroid()
        
        if self.SideChainIsInitialized:
            #compute global and local SAS
            self.AbsoluteSideChainGlobalSAS = self.getAbsoluteSideChainGlobalSAS()
            self.AbsoluteSideChainLocalSAS = self.getAbsoluteSideChainLocalSAS()
            #get appropriate side chain coordinates
            self.Positions['side_chain'] = {}
            self.Positions['side_chain']['scs'] = self.SpecificAtoms['backbone']['c_alpha'].getPointAsString()
            self.Positions['side_chain']['sce'] = self.SpecificAtoms['side_chain']['c_beta'].getPointAsString()

####################################################################################################
#                                                                                                  #
#           CLASS PDBArginine                                                                      #
#                                                                                                  #
#           DESCRIPTION specific class inheriting from PDBResidue with unique side chain atom      #
#                       identities pertaining to amino acid type arginine                          #
#                                                                                                  #
####################################################################################################

class PDBArginine(PDBResidue):
    
    "CONSTRUCTOR"
    def __init__(self,atom_array):
        PDBResidue.__init__(self,atom_array)
        
        self.Info = {}
        self.Info['one_letter'] = 'R'
        self.Info['three_letter'] = 'ARG'
        
        #initialize the side chain
        self.SpecificAtoms['side_chain'] = {}
        self.SpecificAtoms['side_chain']['c_beta'] = self.getSpecificAtom(self.Atoms,'CB')
        self.SpecificAtoms['side_chain']['c_gamma'] = self.getSpecificAtom(self.Atoms,'CG')
        self.SpecificAtoms['side_chain']['c_delta'] = self.getSpecificAtom(self.Atoms,'CD')
        self.SpecificAtoms['side_chain']['n_epsilon'] = self.getSpecificAtom(self.Atoms,'NE')
        self.SpecificAtoms['side_chain']['c_zeta'] = self.getSpecificAtom(self.Atoms,'CZ')
        self.SpecificAtoms['side_chain']['n_eta1'] = self.getSpecificAtom(self.Atoms,'NH1')
        self.SpecificAtoms['side_chain']['n_eta2'] = self.getSpecificAtom(self.Atoms,'NH2')
        
        #check if side chain is initialized
        self.SideChainIsInitialized = self.checkIfProperlyInitialized(self.SpecificAtoms['side_chain'])
        
        #compute centroid
        if self.SideChainIsInitialized and self.BackboneIsInitialized:
            self.Positions['centroid'] = self.compute3Dcentroid()
        
        if self.SideChainIsInitialized:
            #compute global and local SAS
            self.AbsoluteSideChainGlobalSAS = self.getAbsoluteSideChainGlobalSAS()
            self.AbsoluteSideChainLocalSAS = self.getAbsoluteSideChainLocalSAS()
            #get appropriate side chain coordinates
            self.Positions['side_chain'] = {}
            self.Positions['side_chain']['scs'] = self.SpecificAtoms['side_chain']['c_delta'].getPointAsString()
            self.Positions['side_chain']['sce'] = self.computeCentroidOf2Atoms(self.SpecificAtoms['side_chain']['n_eta1'],\
                                                                               self.SpecificAtoms['side_chain']['n_eta2'])

####################################################################################################
#                                                                                                  #
#           CLASS PDBAsparagine                                                                    #
#                                                                                                  #
#           DESCRIPTION specific class inheriting from PDBResidue with unique side chain atom      #
#                       identities pertaining to amino acid type asparagine                        #
#                                                                                                  #
####################################################################################################

class PDBAsparagine(PDBResidue):
    
    "CONSTRUCTOR"
    def __init__(self,atom_array):
        PDBResidue.__init__(self,atom_array)
        
        self.Info = {}
        self.Info['one_letter'] = 'N'
        self.Info['three_letter'] = 'ASN'
        
        #initialize the side chain
        self.SpecificAtoms['side_chain'] = {}
        self.SpecificAtoms['side_chain']['c_beta'] = self.getSpecificAtom(self.Atoms,'CB')
        self.SpecificAtoms['side_chain']['c_gamma'] = self.getSpecificAtom(self.Atoms,'CG')
        self.SpecificAtoms['side_chain']['o_delta1'] = self.getSpecificAtom(self.Atoms,'OD1')
        self.SpecificAtoms['side_chain']['n_delta2'] = self.getSpecificAtom(self.Atoms,'ND2')
        
        #check if side chain is initialized
        self.SideChainIsInitialized = self.checkIfProperlyInitialized(self.SpecificAtoms['side_chain'])
        
        #compute centroid
        if self.SideChainIsInitialized and self.BackboneIsInitialized:
            self.Positions['centroid'] = self.compute3Dcentroid()
        
        if self.SideChainIsInitialized:
            #compute global and local SAS
            self.AbsoluteSideChainGlobalSAS = self.getAbsoluteSideChainGlobalSAS()
            self.AbsoluteSideChainLocalSAS = self.getAbsoluteSideChainLocalSAS()
            #get appropriate side chain coordinates
            self.Positions['side_chain'] = {}
            self.Positions['side_chain']['scs'] = self.SpecificAtoms['side_chain']['c_beta'].getPointAsString()
            self.Positions['side_chain']['sce'] = self.computeCentroidOf2Atoms(self.SpecificAtoms['side_chain']['o_delta1'],\
                                                                               self.SpecificAtoms['side_chain']['n_delta2'])

####################################################################################################
#                                                                                                  #
#           CLASS PDBAsparticAcid                                                                  #
#                                                                                                  #
#           DESCRIPTION specific class inheriting from PDBResidue with unique side chain atom      #
#                       identities pertaining to amino acid type aspartic acid                     #
#                                                                                                  #
####################################################################################################

class PDBAsparticAcid(PDBResidue):
    
    "CONSTRUCTOR"
    def __init__(self,atom_array):
        PDBResidue.__init__(self,atom_array)
        
        self.Info = {}
        self.Info['one_letter'] = 'D'
        self.Info['three_letter'] = 'ASP'
        
        #initialize the side chain
        self.SpecificAtoms['side_chain'] = {}
        self.SpecificAtoms['side_chain']['c_beta'] = self.getSpecificAtom(self.Atoms,'CB')
        self.SpecificAtoms['side_chain']['c_gamma'] = self.getSpecificAtom(self.Atoms,'CG')
        self.SpecificAtoms['side_chain']['o_delta1'] = self.getSpecificAtom(self.Atoms,'OD1')
        self.SpecificAtoms['side_chain']['o_delta2'] = self.getSpecificAtom(self.Atoms,'OD2')
        
        #check if side chain is initialized
        self.SideChainIsInitialized = self.checkIfProperlyInitialized(self.SpecificAtoms['side_chain'])
        
        #compute centroid
        if self.SideChainIsInitialized and self.BackboneIsInitialized:
            self.Positions['centroid'] = self.compute3Dcentroid()
        
        if self.SideChainIsInitialized:
            #compute global and local SAS
            self.AbsoluteSideChainGlobalSAS = self.getAbsoluteSideChainGlobalSAS()
            self.AbsoluteSideChainLocalSAS = self.getAbsoluteSideChainLocalSAS()
            #get appropriate side chain coordinates
            self.Positions['side_chain'] = {}
            self.Positions['side_chain']['scs'] = self.SpecificAtoms['side_chain']['c_beta'].getPointAsString()
            self.Positions['side_chain']['sce'] = self.computeCentroidOf2Atoms(self.SpecificAtoms['side_chain']['o_delta1'],\
                                                                               self.SpecificAtoms['side_chain']['o_delta2'])

####################################################################################################
#                                                                                                  #
#           CLASS PDBCysteine                                                                      #
#                                                                                                  #
#           DESCRIPTION specific class inheriting from PDBResidue with unique side chain atom      #
#                       identities pertaining to amino acid type cysteine                          #
#                                                                                                  #
####################################################################################################

class PDBCysteine(PDBResidue):
    
    "CONSTRUCTOR"
    def __init__(self,atom_array):
        PDBResidue.__init__(self,atom_array)
        
        self.Info = {}
        self.Info['one_letter'] = 'C'
        self.Info['three_letter'] = 'CYS'
        
        #initialize the side chain
        self.SpecificAtoms['side_chain'] = {}
        self.SpecificAtoms['side_chain']['c_beta'] = self.getSpecificAtom(self.Atoms,'CB')
        self.SpecificAtoms['side_chain']['s_gamma'] = self.getSpecificAtom(self.Atoms,'SG')
        
        #check if side chain is initialized
        self.SideChainIsInitialized = self.checkIfProperlyInitialized(self.SpecificAtoms['side_chain'])
        
        #compute centroid
        if self.SideChainIsInitialized and self.BackboneIsInitialized:
            self.Positions['centroid'] = self.compute3Dcentroid()
        
        if self.SideChainIsInitialized:
            #compute global and local SAS
            self.AbsoluteSideChainGlobalSAS = self.getAbsoluteSideChainGlobalSAS()
            self.AbsoluteSideChainLocalSAS = self.getAbsoluteSideChainLocalSAS()
            #get appropriate side chain coordinates
            self.Positions['side_chain'] = {}
            self.Positions['side_chain']['scs'] = self.SpecificAtoms['backbone']['c_alpha'].getPointAsString()
            self.Positions['side_chain']['sce'] = self.SpecificAtoms['side_chain']['s_gamma'].getPointAsString()

####################################################################################################
#                                                                                                  #
#           CLASS PDBGlutamicAcid                                                                  #
#                                                                                                  #
#           DESCRIPTION specific class inheriting from PDBResidue with unique side chain atom      #
#                       identities pertaining to amino acid type glutamic acid                     #
#                                                                                                  #
####################################################################################################

class PDBGlutamicAcid(PDBResidue):
    
    "CONSTRUCTOR"
    def __init__(self,atom_array):
        PDBResidue.__init__(self,atom_array)
        
        self.Info = {}
        self.Info['one_letter'] = 'E'
        self.Info['three_letter'] = 'GLU'
        
        #initialize the side chain
        self.SpecificAtoms['side_chain'] = {}
        self.SpecificAtoms['side_chain']['c_beta'] = self.getSpecificAtom(self.Atoms,'CB')
        self.SpecificAtoms['side_chain']['c_gamma'] = self.getSpecificAtom(self.Atoms,'CG')
        self.SpecificAtoms['side_chain']['c_delta'] = self.getSpecificAtom(self.Atoms,'CD')
        self.SpecificAtoms['side_chain']['o_epsilon1'] = self.getSpecificAtom(self.Atoms,'OE1')
        self.SpecificAtoms['side_chain']['o_epsilon2'] = self.getSpecificAtom(self.Atoms,'OE2')
        
        #check if side chain is initialized
        self.SideChainIsInitialized = self.checkIfProperlyInitialized(self.SpecificAtoms['side_chain'])
        
        #compute centroid
        if self.SideChainIsInitialized and self.BackboneIsInitialized:
            self.Positions['centroid'] = self.compute3Dcentroid()
        
        if self.SideChainIsInitialized:
            #compute global and local SAS
            self.AbsoluteSideChainGlobalSAS = self.getAbsoluteSideChainGlobalSAS()
            self.AbsoluteSideChainLocalSAS = self.getAbsoluteSideChainLocalSAS()
            #get appropriate side chain coordinates
            self.Positions['side_chain'] = {}
            self.Positions['side_chain']['scs'] = self.SpecificAtoms['side_chain']['c_gamma'].getPointAsString()
            self.Positions['side_chain']['sce'] = self.computeCentroidOf2Atoms(self.SpecificAtoms['side_chain']['o_epsilon1'],\
                                                                               self.SpecificAtoms['side_chain']['o_epsilon2'])

####################################################################################################
#                                                                                                  #
#           CLASS PDBGlutamine                                                                     #
#                                                                                                  #
#           DESCRIPTION specific class inheriting from PDBResidue with unique side chain atom      #
#                       identities pertaining to amino acid type glutamine                         #
#                                                                                                  #
####################################################################################################

class PDBGlutamine(PDBResidue):
    
    "CONSTRUCTOR"
    def __init__(self,atom_array):
        PDBResidue.__init__(self,atom_array)
        
        self.Info = {}
        self.Info['one_letter'] = 'Q'
        self.Info['three_letter'] = 'GLN'
        
        #initialize the side chain
        self.SpecificAtoms['side_chain'] = {}
        self.SpecificAtoms['side_chain']['c_beta'] = self.getSpecificAtom(self.Atoms,'CB')
        self.SpecificAtoms['side_chain']['c_gamma'] = self.getSpecificAtom(self.Atoms,'CG')
        self.SpecificAtoms['side_chain']['c_delta'] = self.getSpecificAtom(self.Atoms,'CD')
        self.SpecificAtoms['side_chain']['o_epsilon1'] = self.getSpecificAtom(self.Atoms,'OE1')
        self.SpecificAtoms['side_chain']['n_epsilon2'] = self.getSpecificAtom(self.Atoms,'NE2')
        
        #check if side chain is initialized
        self.SideChainIsInitialized = self.checkIfProperlyInitialized(self.SpecificAtoms['side_chain'])
        
        #compute centroid
        if self.SideChainIsInitialized and self.BackboneIsInitialized:
            self.Positions['centroid'] = self.compute3Dcentroid()
        
        if self.SideChainIsInitialized:
            #compute global and local SAS
            self.AbsoluteSideChainGlobalSAS = self.getAbsoluteSideChainGlobalSAS()
            self.AbsoluteSideChainLocalSAS = self.getAbsoluteSideChainLocalSAS()
            #get appropriate side chain coordinates
            self.Positions['side_chain'] = {}
            self.Positions['side_chain']['scs'] = self.SpecificAtoms['side_chain']['c_gamma'].getPointAsString()
            self.Positions['side_chain']['sce'] = self.computeCentroidOf2Atoms(self.SpecificAtoms['side_chain']['o_epsilon1'],\
                                                                               self.SpecificAtoms['side_chain']['n_epsilon2'])

####################################################################################################
#                                                                                                  #
#           CLASS PDBGlycine                                                                       #
#                                                                                                  #
#           DESCRIPTION specific class inheriting from PDBResidue with unique side chain atom      #
#                       identities pertaining to amino acid type glycine                           #
#                                                                                                  #
####################################################################################################

class PDBGlycine(PDBResidue):
    
    "CONSTRUCTOR"
    def __init__(self,atom_array):
        PDBResidue.__init__(self,atom_array)
        
        self.Info = {}
        self.Info['one_letter'] = 'G'
        self.Info['three_letter'] = 'GLY'
        
        #initialize the side chain
        self.SpecificAtoms['side_chain'] = {}
        
        #check if side chain is initialized
        self.SideChainIsInitialized = False
        
        #compute centroid
        if self.BackboneIsInitialized:
            self.SideChainIsInitialized = True
            
            self.Positions['centroid'] = self.compute3Dcentroid()
            
            #compute global and local SAS
            self.AbsoluteSideChainGlobalSAS = float(self.SpecificAtoms['backbone']['c_alpha'].GlobalSAS)
            self.AbsoluteSideChainLocalSAS = float(self.SpecificAtoms['backbone']['c_alpha'].LocalSAS)
            #get appropriate side chain coordinates
            self.Positions['side_chain'] = {}
            self.Positions['side_chain']['scs'] = self.SpecificAtoms['backbone']['c_alpha'].getPointAsString()
            self.Positions['side_chain']['sce'] = self.computeCentroidOf2Atoms(self.SpecificAtoms['backbone']['n_amino'],\
                                                                               self.SpecificAtoms['backbone']['c_carboxyl'])
    
    "overrides the normal method because glycine has no side chain"
    def getAverageSideChainBFactor(self):
        return self.SpecificAtoms['backbone']['c_alpha'].Info['old_temp_factor']

####################################################################################################
#                                                                                                  #
#           CLASS PDBHistidine                                                                     #
#                                                                                                  #
#           DESCRIPTION specific class inheriting from PDBResidue with unique side chain atom      #
#                       identities pertaining to amino acid type histidine                         #
#                                                                                                  #
####################################################################################################

class PDBHistidine(PDBResidue):
    
    "CONSTRUCTOR"
    def __init__(self,atom_array):
        PDBResidue.__init__(self,atom_array)
        
        self.Info = {}
        self.Info['one_letter'] = 'H'
        self.Info['three_letter'] = 'HIS'
        
        #initialize the side chain
        self.SpecificAtoms['side_chain'] = {}
        self.SpecificAtoms['side_chain']['c_beta'] = self.getSpecificAtom(self.Atoms,'CB')
        self.SpecificAtoms['side_chain']['c_gamma'] = self.getSpecificAtom(self.Atoms,'CG')
        self.SpecificAtoms['side_chain']['n_delta1'] = self.getSpecificAtom(self.Atoms,'ND1')
        self.SpecificAtoms['side_chain']['c_delta2'] = self.getSpecificAtom(self.Atoms,'CD2')
        self.SpecificAtoms['side_chain']['c_epsilon1'] = self.getSpecificAtom(self.Atoms,'CE1')
        self.SpecificAtoms['side_chain']['n_epsilon2'] = self.getSpecificAtom(self.Atoms,'NE2')
        
        #check if side chain is initialized
        self.SideChainIsInitialized = self.checkIfProperlyInitialized(self.SpecificAtoms['side_chain'])
        
        #compute centroid
        if self.SideChainIsInitialized and self.BackboneIsInitialized:
            self.Positions['centroid'] = self.compute3Dcentroid()
        
        if self.SideChainIsInitialized:
            #compute global and local SAS
            self.AbsoluteSideChainGlobalSAS = self.getAbsoluteSideChainGlobalSAS()
            self.AbsoluteSideChainLocalSAS = self.getAbsoluteSideChainLocalSAS()
            #get appropriate side chain coordinates
            self.Positions['side_chain'] = {}
            self.Positions['side_chain']['scs'] = self.SpecificAtoms['side_chain']['c_gamma'].getPointAsString()
            self.Positions['side_chain']['sce'] = self.computeCentroidOf2Atoms(self.SpecificAtoms['side_chain']['c_epsilon1'],\
                                                                               self.SpecificAtoms['side_chain']['n_epsilon2'])

####################################################################################################
#                                                                                                  #
#           CLASS PDBIsoleucine                                                                    #
#                                                                                                  #
#           DESCRIPTION specific class inheriting from PDBResidue with unique side chain atom      #
#                       identities pertaining to amino acid type isoleucine                        #
#                                                                                                  #
####################################################################################################

class PDBIsoleucine(PDBResidue):
    
    "CONSTRUCTOR"
    def __init__(self,atom_array):
        PDBResidue.__init__(self,atom_array)
        
        self.Info = {}
        self.Info['one_letter'] = 'I'
        self.Info['three_letter'] = 'ILE'
        
        #initialize the side chain
        self.SpecificAtoms['side_chain'] = {}
        self.SpecificAtoms['side_chain']['c_beta'] = self.getSpecificAtom(self.Atoms,'CB')
        self.SpecificAtoms['side_chain']['c_gamma1'] = self.getSpecificAtom(self.Atoms,'CG1')
        self.SpecificAtoms['side_chain']['c_gamma2'] = self.getSpecificAtom(self.Atoms,'CG2')
        self.SpecificAtoms['side_chain']['c_delta1'] = self.getSpecificAtom(self.Atoms,'CD1')
        
        #check if side chain is initialized
        self.SideChainIsInitialized = self.checkIfProperlyInitialized(self.SpecificAtoms['side_chain'])
        
        #compute centroid
        if self.SideChainIsInitialized and self.BackboneIsInitialized:
            self.Positions['centroid'] = self.compute3Dcentroid()
        
        if self.SideChainIsInitialized:
            #compute global and local SAS
            self.AbsoluteSideChainGlobalSAS = self.getAbsoluteSideChainGlobalSAS()
            self.AbsoluteSideChainLocalSAS = self.getAbsoluteSideChainLocalSAS()
            #get appropriate side chain coordinates
            self.Positions['side_chain'] = {}
            self.Positions['side_chain']['scs'] = self.SpecificAtoms['side_chain']['c_beta'].getPointAsString()
            self.Positions['side_chain']['sce'] = self.SpecificAtoms['side_chain']['c_delta1'].getPointAsString()

####################################################################################################
#                                                                                                  #
#           CLASS PDBLeucine                                                                       #
#                                                                                                  #
#           DESCRIPTION specific class inheriting from PDBResidue with unique side chain atom      #
#                       identities pertaining to amino acid type leucine                           #
#                                                                                                  #
####################################################################################################

class PDBLeucine(PDBResidue):
    
    "CONSTRUCTOR"
    def __init__(self,atom_array):
        PDBResidue.__init__(self,atom_array)
        
        self.Info = {}
        self.Info['one_letter'] = 'L'
        self.Info['three_letter'] = 'LEU'
        
        #initialize the side chain
        self.SpecificAtoms['side_chain'] = {}
        self.SpecificAtoms['side_chain']['c_beta'] = self.getSpecificAtom(self.Atoms,'CB')
        self.SpecificAtoms['side_chain']['c_gamma'] = self.getSpecificAtom(self.Atoms,'CG')
        self.SpecificAtoms['side_chain']['c_delta1'] = self.getSpecificAtom(self.Atoms,'CD1')
        self.SpecificAtoms['side_chain']['c_delta2'] = self.getSpecificAtom(self.Atoms,'CD2')
        
        #check if side chain is initialized
        self.SideChainIsInitialized = self.checkIfProperlyInitialized(self.SpecificAtoms['side_chain'])
        
        #compute centroid
        if self.SideChainIsInitialized and self.BackboneIsInitialized:
            self.Positions['centroid'] = self.compute3Dcentroid()
        
        if self.SideChainIsInitialized:
            #compute global and local SAS
            self.AbsoluteSideChainGlobalSAS = self.getAbsoluteSideChainGlobalSAS()
            self.AbsoluteSideChainLocalSAS = self.getAbsoluteSideChainLocalSAS()
            #get appropriate side chain coordinates
            self.Positions['side_chain'] = {}
            self.Positions['side_chain']['scs'] = self.SpecificAtoms['side_chain']['c_beta'].getPointAsString()
            self.Positions['side_chain']['sce'] = self.computeCentroidOf2Atoms(self.SpecificAtoms['side_chain']['c_delta1'],\
                                                                               self.SpecificAtoms['side_chain']['c_delta2'])

####################################################################################################
#                                                                                                  #
#           CLASS PDBLysine                                                                        #
#                                                                                                  #
#           DESCRIPTION specific class inheriting from PDBResidue with unique side chain atom      #
#                       identities pertaining to amino acid type lysine                            #
#                                                                                                  #
####################################################################################################

class PDBLysine(PDBResidue):
    
    "CONSTRUCTOR"
    def __init__(self,atom_array):
        PDBResidue.__init__(self,atom_array)
        
        self.Info = {}
        self.Info['one_letter'] = 'K'
        self.Info['three_letter'] = 'LYS'
        
        #initialize the side chain
        self.SpecificAtoms['side_chain'] = {}
        self.SpecificAtoms['side_chain']['c_beta'] = self.getSpecificAtom(self.Atoms,'CB')
        self.SpecificAtoms['side_chain']['c_gamma'] = self.getSpecificAtom(self.Atoms,'CG')
        self.SpecificAtoms['side_chain']['c_delta'] = self.getSpecificAtom(self.Atoms,'CD')
        self.SpecificAtoms['side_chain']['c_epsilon'] = self.getSpecificAtom(self.Atoms,'CE')
        self.SpecificAtoms['side_chain']['n_zeta'] = self.getSpecificAtom(self.Atoms,'NZ')
        
        #check if side chain is initialized
        self.SideChainIsInitialized = self.checkIfProperlyInitialized(self.SpecificAtoms['side_chain'])
        
        #compute centroid
        if self.SideChainIsInitialized and self.BackboneIsInitialized:
            self.Positions['centroid'] = self.compute3Dcentroid()
        
        if self.SideChainIsInitialized:
            #compute global and local SAS
            self.AbsoluteSideChainGlobalSAS = self.getAbsoluteSideChainGlobalSAS()
            self.AbsoluteSideChainLocalSAS = self.getAbsoluteSideChainLocalSAS()
            #get appropriate side chain coordinates
            self.Positions['side_chain'] = {}
            self.Positions['side_chain']['scs'] = self.SpecificAtoms['side_chain']['c_gamma'].getPointAsString()
            self.Positions['side_chain']['sce'] = self.SpecificAtoms['side_chain']['n_zeta'].getPointAsString()

####################################################################################################
#                                                                                                  #
#           CLASS PDBMethionine                                                                    #
#                                                                                                  #
#           DESCRIPTION specific class inheriting from PDBResidue with unique side chain atom      #
#                       identities pertaining to amino acid type methionine                        #
#                                                                                                  #
####################################################################################################

class PDBMethionine(PDBResidue):
    
    "CONSTRUCTOR"
    def __init__(self,atom_array):
        PDBResidue.__init__(self,atom_array)
        
        self.Info = {}
        self.Info['one_letter'] = 'M'
        self.Info['three_letter'] = 'MET'
        
        #initialize the side chain
        self.SpecificAtoms['side_chain'] = {}
        self.SpecificAtoms['side_chain']['c_beta'] = self.getSpecificAtom(self.Atoms,'CB')
        self.SpecificAtoms['side_chain']['c_gamma'] = self.getSpecificAtom(self.Atoms,'CG')
        self.SpecificAtoms['side_chain']['s_delta'] = self.getSpecificAtom(self.Atoms,'SD')
        self.SpecificAtoms['side_chain']['c_epsilon'] = self.getSpecificAtom(self.Atoms,'CE')
        
        #check if side chain is initialized
        self.SideChainIsInitialized = self.checkIfProperlyInitialized(self.SpecificAtoms['side_chain'])
        
        #compute centroid
        if self.SideChainIsInitialized and self.BackboneIsInitialized:
            self.Positions['centroid'] = self.compute3Dcentroid()
        
        if self.SideChainIsInitialized:
            #compute global and local SAS
            self.AbsoluteSideChainGlobalSAS = self.getAbsoluteSideChainGlobalSAS()
            self.AbsoluteSideChainLocalSAS = self.getAbsoluteSideChainLocalSAS()
            #get appropriate side chain coordinates
            self.Positions['side_chain'] = {}
            self.Positions['side_chain']['scs'] = self.SpecificAtoms['backbone']['c_alpha'].getPointAsString()
            self.Positions['side_chain']['sce'] = self.SpecificAtoms['side_chain']['s_delta'].getPointAsString()

####################################################################################################
#                                                                                                  #
#           CLASS PDBPhenylalanine                                                                 #
#                                                                                                  #
#           DESCRIPTION specific class inheriting from PDBResidue with unique side chain atom      #
#                       identities pertaining to amino acid type phenylalanine                     #
#                                                                                                  #
####################################################################################################

class PDBPhenylalanine(PDBResidue):
    
    "CONSTRUCTOR"
    def __init__(self,atom_array):
        PDBResidue.__init__(self,atom_array)
        
        self.Info = {}
        self.Info['one_letter'] = 'F'
        self.Info['three_letter'] = 'PHE'
        
        #initialize the side chain
        self.SpecificAtoms['side_chain'] = {}
        self.SpecificAtoms['side_chain']['c_beta'] = self.getSpecificAtom(self.Atoms,'CB')
        self.SpecificAtoms['side_chain']['c_gamma'] = self.getSpecificAtom(self.Atoms,'CG')
        self.SpecificAtoms['side_chain']['c_delta1'] = self.getSpecificAtom(self.Atoms,'CD1')
        self.SpecificAtoms['side_chain']['c_delta2'] = self.getSpecificAtom(self.Atoms,'CD2')
        self.SpecificAtoms['side_chain']['c_epsilon1'] = self.getSpecificAtom(self.Atoms,'CE1')
        self.SpecificAtoms['side_chain']['c_epsilon2'] = self.getSpecificAtom(self.Atoms,'CE2')
        self.SpecificAtoms['side_chain']['c_zeta'] = self.getSpecificAtom(self.Atoms,'CZ')
        
        #check if side chain is initialized
        self.SideChainIsInitialized = self.checkIfProperlyInitialized(self.SpecificAtoms['side_chain'])
        
        #compute centroid
        if self.SideChainIsInitialized and self.BackboneIsInitialized:
            self.Positions['centroid'] = self.compute3Dcentroid()
        
        if self.SideChainIsInitialized:
            #compute global and local SAS
            self.AbsoluteSideChainGlobalSAS = self.getAbsoluteSideChainGlobalSAS()
            self.AbsoluteSideChainLocalSAS = self.getAbsoluteSideChainLocalSAS()
            #get appropriate side chain coordinates
            self.Positions['side_chain'] = {}
            self.Positions['side_chain']['scs'] = self.SpecificAtoms['side_chain']['c_gamma'].getPointAsString()
            self.Positions['side_chain']['sce'] = self.SpecificAtoms['side_chain']['c_zeta'].getPointAsString()

####################################################################################################
#                                                                                                  #
#           CLASS PDBProline                                                                       #
#                                                                                                  #
#           DESCRIPTION specific class inheriting from PDBResidue with unique side chain atom      #
#                       identities pertaining to amino acid type proline                           #
#                                                                                                  #
####################################################################################################

class PDBProline(PDBResidue):
    
    "CONSTRUCTOR"
    def __init__(self,atom_array):
        PDBResidue.__init__(self,atom_array)
        
        self.Info = {}
        self.Info['one_letter'] = 'P'
        self.Info['three_letter'] = 'PRO'
        
        #initialize the side chain
        self.SpecificAtoms['side_chain'] = {}
        self.SpecificAtoms['side_chain']['c_beta'] = self.getSpecificAtom(self.Atoms,'CB')
        self.SpecificAtoms['side_chain']['c_gamma'] = self.getSpecificAtom(self.Atoms,'CG')
        self.SpecificAtoms['side_chain']['c_delta'] = self.getSpecificAtom(self.Atoms,'CD')
        
        #check if side chain is initialized
        self.SideChainIsInitialized = self.checkIfProperlyInitialized(self.SpecificAtoms['side_chain'])
        
        #compute centroid
        if self.SideChainIsInitialized and self.BackboneIsInitialized:
            self.Positions['centroid'] = self.compute3Dcentroid()
        
        if self.SideChainIsInitialized:
            #compute global and local SAS
            self.AbsoluteSideChainGlobalSAS = self.getAbsoluteSideChainGlobalSAS()
            self.AbsoluteSideChainLocalSAS = self.getAbsoluteSideChainLocalSAS()
            #get appropriate side chain coordinates
            self.Positions['side_chain'] = {}
            self.Positions['side_chain']['scs'] = self.SpecificAtoms['backbone']['c_alpha'].getPointAsString()
            self.Positions['side_chain']['sce'] = self.computeCentroidOf2Atoms(self.SpecificAtoms['side_chain']['c_gamma'],\
                                                                               self.SpecificAtoms['side_chain']['c_delta'])

####################################################################################################
#                                                                                                  #
#           CLASS PDBSerine                                                                        #
#                                                                                                  #
#           DESCRIPTION specific class inheriting from PDBResidue with unique side chain atom      #
#                       identities pertaining to amino acid type serine                            #
#                                                                                                  #
####################################################################################################

class PDBSerine(PDBResidue):
    
    "CONSTRUCTOR"
    def __init__(self,atom_array):
        PDBResidue.__init__(self,atom_array)
        
        self.Info = {}
        self.Info['one_letter'] = 'S'
        self.Info['three_letter'] = 'SER'
        
        #initialize the side chain
        self.SpecificAtoms['side_chain'] = {}
        self.SpecificAtoms['side_chain']['c_beta'] = self.getSpecificAtom(self.Atoms,'CB')
        self.SpecificAtoms['side_chain']['o_gamma'] = self.getSpecificAtom(self.Atoms,'OG')
        
        #check if side chain is initialized
        self.SideChainIsInitialized = self.checkIfProperlyInitialized(self.SpecificAtoms['side_chain'])
        
        #compute centroid
        if self.SideChainIsInitialized and self.BackboneIsInitialized:
            self.Positions['centroid'] = self.compute3Dcentroid()
        
        if self.SideChainIsInitialized:
            #compute global and local SAS
            self.AbsoluteSideChainGlobalSAS = self.getAbsoluteSideChainGlobalSAS()
            self.AbsoluteSideChainLocalSAS = self.getAbsoluteSideChainLocalSAS()
            #get appropriate side chain coordinates
            self.Positions['side_chain'] = {}
            self.Positions['side_chain']['scs'] = self.SpecificAtoms['backbone']['c_alpha'].getPointAsString()
            self.Positions['side_chain']['sce'] = self.SpecificAtoms['side_chain']['o_gamma'].getPointAsString()

####################################################################################################
#                                                                                                  #
#           CLASS PDBThreonine                                                                     #
#                                                                                                  #
#           DESCRIPTION specific class inheriting from PDBResidue with unique side chain atom      #
#                       identities pertaining to amino acid type threonine                         #
#                                                                                                  #
####################################################################################################

class PDBThreonine(PDBResidue):
    
    "CONSTRUCTOR"
    def __init__(self,atom_array):
        PDBResidue.__init__(self,atom_array)
        
        self.Info = {}
        self.Info['one_letter'] = 'T'
        self.Info['three_letter'] = 'THR'
        
        #initialize the side chain
        self.SpecificAtoms['side_chain'] = {}
        self.SpecificAtoms['side_chain']['c_beta'] = self.getSpecificAtom(self.Atoms,'CB')
        self.SpecificAtoms['side_chain']['o_gamma1'] = self.getSpecificAtom(self.Atoms,'OG1')
        self.SpecificAtoms['side_chain']['c_gamma2'] = self.getSpecificAtom(self.Atoms,'CG2')
        
        #check if side chain is initialized
        self.SideChainIsInitialized = self.checkIfProperlyInitialized(self.SpecificAtoms['side_chain'])
        
        #compute centroid
        if self.SideChainIsInitialized and self.BackboneIsInitialized:
            self.Positions['centroid'] = self.compute3Dcentroid()
        
        if self.SideChainIsInitialized:
            #compute global and local SAS
            self.AbsoluteSideChainGlobalSAS = self.getAbsoluteSideChainGlobalSAS()
            self.AbsoluteSideChainLocalSAS = self.getAbsoluteSideChainLocalSAS()
            #get appropriate side chain coordinates
            self.Positions['side_chain'] = {}
            self.Positions['side_chain']['scs'] = self.SpecificAtoms['backbone']['c_alpha'].getPointAsString()
            self.Positions['side_chain']['sce'] = self.computeCentroidOf2Atoms(self.SpecificAtoms['side_chain']['c_gamma2'],\
                                                                               self.SpecificAtoms['side_chain']['o_gamma1'])

####################################################################################################
#                                                                                                  #
#           CLASS PDBTryptophan                                                                    #
#                                                                                                  #
#           DESCRIPTION specific class inheriting from PDBResidue with unique side chain atom      #
#                       identities pertaining to amino acid type tryptophan                        #
#                                                                                                  #
####################################################################################################

class PDBTryptophan(PDBResidue):
    
    "CONSTRUCTOR"
    def __init__(self,atom_array):
        PDBResidue.__init__(self,atom_array)
        
        self.Info = {}
        self.Info['one_letter'] = 'W'
        self.Info['three_letter'] = 'TRP'
        
        #initialize the side chain
        self.SpecificAtoms['side_chain'] = {}
        self.SpecificAtoms['side_chain']['c_beta'] = self.getSpecificAtom(self.Atoms,'CB')
        self.SpecificAtoms['side_chain']['c_gamma'] = self.getSpecificAtom(self.Atoms,'CG')
        self.SpecificAtoms['side_chain']['c_delta1'] = self.getSpecificAtom(self.Atoms,'CD1')
        self.SpecificAtoms['side_chain']['c_delta2'] = self.getSpecificAtom(self.Atoms,'CD2')
        self.SpecificAtoms['side_chain']['n_epsilon1'] = self.getSpecificAtom(self.Atoms,'NE1')
        self.SpecificAtoms['side_chain']['c_epsilon2'] = self.getSpecificAtom(self.Atoms,'CE2')
        self.SpecificAtoms['side_chain']['c_epsilon3'] = self.getSpecificAtom(self.Atoms,'CE3')
        self.SpecificAtoms['side_chain']['c_zeta2'] = self.getSpecificAtom(self.Atoms,'CZ2')
        self.SpecificAtoms['side_chain']['c_zeta3'] = self.getSpecificAtom(self.Atoms,'CZ3')
        self.SpecificAtoms['side_chain']['c_eta2'] = self.getSpecificAtom(self.Atoms,'CH2')
        
        #check if side chain is initialized
        self.SideChainIsInitialized = self.checkIfProperlyInitialized(self.SpecificAtoms['side_chain'])
        
        #compute centroid
        if self.SideChainIsInitialized and self.BackboneIsInitialized:
            self.Positions['centroid'] = self.compute3Dcentroid()
        
        if self.SideChainIsInitialized:
            #compute global and local SAS
            self.AbsoluteSideChainGlobalSAS = self.getAbsoluteSideChainGlobalSAS()
            self.AbsoluteSideChainLocalSAS = self.getAbsoluteSideChainLocalSAS()
            #get appropriate side chain coordinates
            self.Positions['side_chain'] = {}
            self.Positions['side_chain']['scs'] = self.SpecificAtoms['side_chain']['c_delta1'].getPointAsString()
            self.Positions['side_chain']['sce'] = self.computeCentroidOf2Atoms(self.SpecificAtoms['side_chain']['c_zeta3'],\
                                                                               self.SpecificAtoms['side_chain']['c_eta2'])

####################################################################################################
#                                                                                                  #
#           CLASS PDBTyrosine                                                                      #
#                                                                                                  #
#           DESCRIPTION specific class inheriting from PDBResidue with unique side chain atom      #
#                       identities pertaining to amino acid type tyrosine                          #
#                                                                                                  #
####################################################################################################

class PDBTyrosine(PDBResidue):
    
    "CONSTRUCTOR"
    def __init__(self,atom_array):
        PDBResidue.__init__(self,atom_array)
        
        self.Info = {}
        self.Info['one_letter'] = 'Y'
        self.Info['three_letter'] = 'TYR'
        
        #initialize the side chain
        self.SpecificAtoms['side_chain'] = {}
        self.SpecificAtoms['side_chain']['c_beta'] = self.getSpecificAtom(self.Atoms,'CB')
        self.SpecificAtoms['side_chain']['c_gamma'] = self.getSpecificAtom(self.Atoms,'CG')
        self.SpecificAtoms['side_chain']['c_delta1'] = self.getSpecificAtom(self.Atoms,'CD1')
        self.SpecificAtoms['side_chain']['c_delta2'] = self.getSpecificAtom(self.Atoms,'CD2')
        self.SpecificAtoms['side_chain']['c_epsilon1'] = self.getSpecificAtom(self.Atoms,'CE1')
        self.SpecificAtoms['side_chain']['c_epsilon2'] = self.getSpecificAtom(self.Atoms,'CE2')
        self.SpecificAtoms['side_chain']['c_zeta'] = self.getSpecificAtom(self.Atoms,'CZ')
        self.SpecificAtoms['side_chain']['o_eta'] = self.getSpecificAtom(self.Atoms,'OH')
        
        #check if side chain is initialized
        self.SideChainIsInitialized = self.checkIfProperlyInitialized(self.SpecificAtoms['side_chain'])
        
        #compute centroid
        if self.SideChainIsInitialized and self.BackboneIsInitialized:
            self.Positions['centroid'] = self.compute3Dcentroid()
        
        if self.SideChainIsInitialized:
            #compute global and local SAS
            self.AbsoluteSideChainGlobalSAS = self.getAbsoluteSideChainGlobalSAS()
            self.AbsoluteSideChainLocalSAS = self.getAbsoluteSideChainLocalSAS()
            #get appropriate side chain coordinates
            self.Positions['side_chain'] = {}
            self.Positions['side_chain']['scs'] = self.SpecificAtoms['side_chain']['c_gamma'].getPointAsString()
            self.Positions['side_chain']['sce'] = self.SpecificAtoms['side_chain']['c_zeta'].getPointAsString()

####################################################################################################
#                                                                                                  #
#           CLASS PDBValine                                                                        #
#                                                                                                  #
#           DESCRIPTION specific class inheriting from PDBResidue with unique side chain atom      #
#                       identities pertaining to amino acid type valine                            #
#                                                                                                  #
####################################################################################################

class PDBValine(PDBResidue):
    
    "CONSTRUCTOR"
    def __init__(self,atom_array):
        PDBResidue.__init__(self,atom_array)
        
        self.Info = {}
        self.Info['one_letter'] = 'V'
        self.Info['three_letter'] = 'VAL'
        
        #initialize the side chain
        self.SpecificAtoms['side_chain'] = {}
        self.SpecificAtoms['side_chain']['c_beta'] = self.getSpecificAtom(self.Atoms,'CB')
        self.SpecificAtoms['side_chain']['c_gamma1'] = self.getSpecificAtom(self.Atoms,'CG1')
        self.SpecificAtoms['side_chain']['c_gamma2'] = self.getSpecificAtom(self.Atoms,'CG2')
        
        #check if side chain is initialized
        self.SideChainIsInitialized = self.checkIfProperlyInitialized(self.SpecificAtoms['side_chain'])
        
        #compute centroid
        if self.SideChainIsInitialized and self.BackboneIsInitialized:
            self.Positions['centroid'] = self.compute3Dcentroid()
        
        if self.SideChainIsInitialized:
            #compute global and local SAS
            self.AbsoluteSideChainGlobalSAS = self.getAbsoluteSideChainGlobalSAS()
            self.AbsoluteSideChainLocalSAS = self.getAbsoluteSideChainLocalSAS()
            #get appropriate side chain coordinates
            self.Positions['side_chain'] = {}
            self.Positions['side_chain']['scs'] = self.SpecificAtoms['backbone']['c_alpha'].getPointAsString()
            self.Positions['side_chain']['sce'] = self.computeCentroidOf2Atoms(self.SpecificAtoms['side_chain']['c_gamma1'],\
                                                                               self.SpecificAtoms['side_chain']['c_gamma2'])

####################################################################################################
#                                                                                                  #
#           CLASS PDBFictionalResidue                                                              #
#                                                                                                  #
#           DESCRIPTION specific class inheriting from PDBResidue where a residue is unresolved in #
#                       a PDB structure. Serves largely as a placeholder                           #
#                                                                                                  #
####################################################################################################

class PDBFictionalResidue(PDBResidue):
    
    "CONSTRUCTOR"
    def __init__(self,atom_array):
        PDBResidue.__init__(self,atom_array)
        
        #all fields are given null or false values because the structure is not resolved at this point
        
        self.SpecificAtoms['side_chain'] = {}
        
        self.SideChainIsInitialized = True
        self.Info = {}
        self.Info['one_letter'] = 'X'
        self.Info['three_letter'] = 'XXX'
        
        self.AbsoluteSideChainGlobalSAS = "0.0"
        self.AbsoluteSideChainLocalSAS = "0.0"
        
        self.Positions['centroid'] = "0,0,0"
        self.Positions['backbone']['c_alpha'] = "0,0,0"
        self.Positions['backbone']['n_amino'] = "0,0,0"
        self.Positions['backbone']['c_carboxyl'] = "0,0,0"
        
        self.Positions['side_chain'] = {}
        self.Positions['side_chain']['scs'] = "0,0,0"
        self.Positions['side_chain']['sce'] = "0,0,0"
