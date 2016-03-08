####################################################################################################
#                                                                                                  #
#           PROJECT Protein Adaptation                                                             #
#           Script UpdatePDBXML                                                                    #
#           PROGRAMMER Jeremy Adams                                                                #
#           STARTED 27-11-14                                                                       #
#           LASTMOD 27-11-14                                                                       #
#                                                                                                  #
#           DESCRIPTION script to convert new PDB coordinate files to their corresponding          #
#                       XML format and place the new XML file in the appropriate directory         #
#                                                                                                  #
####################################################################################################

import os
import sys

"adds a slash to a directory path if it is needed"
def addSlashIfNeeded(PATH):
    Ret = PATH
    if PATH.endswith("/"):
        pass
    else:
        Ret = PATH+"/"
    return Ret

USAGE = """python updatepdbxml.py --data-dir [datadir]

Where:
datadir: the absolute path to the directory containing this script
"""

#performs all necessary checks, prints the USAGE string if there is incorrect input
if len(sys.argv) != 3:
    print USAGE
else:
    if sys.argv[1] != "--data-dir":
        print USAGE
    else:
        
        if os.path.exists(addSlashIfNeeded(sys.argv[2])+"pdb/") and os.path.exists(addSlashIfNeeded(sys.argv[2])+"pdbxml/"):
            sys.path.append(addSlashIfNeeded(sys.argv[2])+"src")
            
            from pdbwritemethods import *
            
            PDBDir = addSlashIfNeeded(sys.argv[2])+"pdb/"
            PDBXMLDir = addSlashIfNeeded(sys.argv[2])+"pdbxml/"
            
            #gets all the subdirectories for the PDB and PDBXML directories
            PDBSubDIR_L = sorted(os.listdir(PDBDir))
            PDBXMLSubDIR_L = sorted(os.listdir(PDBXMLDir))
            
            #makes subdirectory in the PDBXML directory if it exists in the PDB directory
            for PDBSubDIR in PDBSubDIR_L:
                if PDBSubDIR in set(PDBXMLSubDIR_L):
                    pass
                else:
                    print "making directory " + PDBSubDIR
                    os.system("mkdir " +PDBXMLDir+PDBSubDIR)
            
            FilesToWrite = []
            
            #gets a list of all files to be converted to PDBXML format
            for PDBSubDIR in PDBSubDIR_L:
                for PDBFile in sorted(os.listdir(PDBDir+PDBSubDIR)):
                    PDBFileToXML = PDBFile.replace(".ent",".xml")
                    
                    if os.path.exists(PDBXMLDir+PDBSubDIR+"/"+PDBFileToXML):
                        pass
                    else:
                        FilesToWrite.append(PDBFileToXML.replace("pdb","").replace(".xml",""))
            
            #executes the method to write the files
            for File in FilesToWrite:
                
                R = writeSinglePDBXMLFile(PDBDir , PDBXMLDir , File[1:3],File)
                
                if R == "Success":
                    print File + " was successfully converted to PDBXML format"
                else:
                    print File + " was not converted to PDBXML format"
        
        else:
            print USAGE
        