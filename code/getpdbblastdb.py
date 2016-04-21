####################################################################################################
#                                                                                                  #
#           PROJECT Protein Adaptation                                                             #
#           Script GetPDBBlastDB                                                                   #
#           PROGRAMMER Jeremy Adams                                                                #
#           STARTED 10-01-15                                                                       #
#           LASTMOD 10-01-15                                                                       #
#                                                                                                  #
#           DESCRIPTION script to download the PDB sequences in Fasta format and make a blast      #
#                       database out of it                                                         #
#                                                                                                  #
####################################################################################################

import urllib2
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

ftp = "ftp://ftp.ncbi.nlm.nih.gov/blast/db/pdbaa.tar.gz"

USAGE = """python getpdbblastdb.py --data-dir [datadir]

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
        if os.path.exists(addSlashIfNeeded(sys.argv[2])+"pdb/"):
            
            os.system("wget %s" % (ftp)) #fetches the file from the FTP
            os.system("gunzip pdbaa.tar.gz")
            os.system("tar -xvf pdbaa.tar")
            os.system("mv pdbaa.* %smisc" % addSlashIfNeeded(sys.argv[2]))
            
        else:
            print USAGE