####################################################################################################
#                                                                                                  #
#           PROJECT Protein Adaptation                                                             #
#           CLASS PDBAtom                                                                          #
#           PROGRAMMER Jeremy Adams                                                                #
#           STARTED 25-07-13                                                                       #
#           LASTMOD 08-08-14                                                                       #
#                                                                                                  #
#           DESCRIPTION Low level class representation the specification of a single atom from     #
#                       a PDB file                                                                 #
#                                                                                                  #
####################################################################################################

class PDBAtom:
    
    "CONSTRUCTOR"
    def __init__(self,OrigLine):
        
        """
        Class attributes:
        Info (Dict): contains specifications about the element, including position, element, host chain, etc.
        """
        
        self.Info = {}
        self.Info['name'] = OrigLine[13:16].replace(" ","")
        self.Info['residue_name'] = OrigLine[17:20].replace(" ","")
        self.Info['chain'] = OrigLine[21].replace(" ","")
        self.Info['residue_num'] = OrigLine[22:26].replace(" ","")
        self.Info['x'] = float(OrigLine[30:38].replace(" ",""))
        self.Info['y'] = float(OrigLine[38:46].replace(" ",""))
        self.Info['z'] = float(OrigLine[46:54].replace(" ",""))
        self.Info['occupancy'] = OrigLine[56:60].replace(" ","")
        self.Info['old_temp_factor'] = float(OrigLine[61:66].replace(" " , ""))
        self.Info['new_temp_factor'] = 0.00
        self.Info['original_line'] = OrigLine
    
    "produces an artificial line with the temperature factor altered to highlight mutations"
    def create_new_line(self):
        self.Info['new_temp_string'] = str(round(self.Info['new_temp_factor'] , 2))
        self.Info['new_line'] = self.Info['original_line'][0:61] + " " * (5 - len(self.Info['new_temp_string'])) + self.Info['new_temp_string'] + self.Info['original_line'][66:-1]
    
    "gets the 3-D coordinates of the atom in string format"
    def getPointAsString(self):
        return "%s,%s,%s" % (str(self.Info['x']) , (self.Info['y']) , (self.Info['z']))
    
    "gets the 3-D coordinates of the atom in float format"
    def getPointAsFloat(self):
        return [self.Info['x'],self.Info['y'],self.Info['z']]
    
    