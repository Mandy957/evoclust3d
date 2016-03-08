####################################################################################################
#                                                                                                  #
#           PROJECT Protein Adaptation                                                             #
#           Script UpdatePDB                                                                       #
#           PROGRAMMER Jeremy Adams                                                                #
#           STARTED 27-11-14                                                                       #
#           LASTMOD 27-11-14                                                                       #
#                                                                                                  #
#           DESCRIPTION script to download new PDB coordinate files to the appropriate directory   #
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


ftp = "ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/"


USAGE = """python updatepdb.py --data-dir [datadir]

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
            pdbDIR = addSlashIfNeeded(sys.argv[2]+"pdb/")


            allDir_L = [line.replace("\n","").split()[-1] for line in urllib2.urlopen(ftp).readlines()] #gets all directories from the FTP site
            unmadeDir_L = [Dir for Dir in allDir_L if os.path.exists(pdbDIR+Dir) == False] #gets list of directories that have not been made locally
            
            #makes local directories for the previously unmade directories
            for unmadeDir in unmadeDir_L:
                print "making directory "+unmadeDir
                os.system("mkdir %s%s" % (pdbDIR,unmadeDir))
            
            #for each directory in the FTP
            for Dir in allDir_L:
                print "Looking in directory "+Dir
                allFile_L = [line.replace("\n","").split()[-1] for line in urllib2.urlopen(ftp+Dir).readlines()] #gets all the gzipped pdb files
                
                #for each file within the directory
                for File in allFile_L:
                    fileWOGZ = File.replace(".gz","")
                    
                    #checks if the file already exists locally
                    if os.path.exists(pdbDIR+Dir+"/"+fileWOGZ):
                        pass
                    
                    #if the file does not exist locally
                    else:
                        print "fetching file" + File
                        os.system("wget %s" % (ftp+Dir+"/"+File)) #fetches the file from the FTP
                        os.system("gunzip %s" % (File)) #unzips the file
                        os.system("mv %s %s" % (fileWOGZ , pdbDIR+Dir)) #moves the file to the correct local subdirectory
        
        else:
            print USAGE