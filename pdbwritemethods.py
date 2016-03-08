####################################################################################################
#                                                                                                  #
#           PROJECT Protein Adaptation                                                             #
#           FILE PDBWriteMethods                                                                   #
#           PROGRAMMER Jeremy Adams                                                                #
#           STARTED 07-01-14                                                                       #
#           LASTMOD 10-01-15                                                                       #
#                                                                                                  #
#           DESCRIPTION methods containing ways to write different combinations of PDBXML files    #
#                                                                                                  #
####################################################################################################

import re
import os
from pdbfile import PDBFile

"writes a single instance of PDBFile with a given directory and PDB accession"
def writeSinglePDBXMLFile(pdbDIR , pdbXMLDIR , d , acc):
    
    ret = "Success"
    try:
        pdb = PDBFile("%s%s/pdb%s.ent" % (pdbDIR , d , acc))
    
        with open("%s%s/pdb%s.xml" % (pdbXMLDIR , d , acc) , "w") as xml:
            xml.write(pdb.xmlWrite())
            xml.flush()
            os.fsync(xml)
        
    except Exception as e:
        ret = "Failure"
        
    return ret
    