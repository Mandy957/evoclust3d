####################################################################################################
#                                                                                                  #
#           PROJECT Protein Adaptation                                                             #
#           SCRIPT Adapt                                                                           #
#           PROGRAMMER Jeremy Adams                                                                #
#           STARTED 10-01-15                                                                       #
#           LASTMOD 10-01-15                                                                       #
#                                                                                                  #
#           DESCRIPTION checks command line input and then executes all stages of the algorithm    #
#                                                                                                  #
####################################################################################################

import sys
import os

def addSlashIfNeeded(PATH):
    Ret = PATH
    if PATH.endswith("/"):
        pass
    else:
        Ret = PATH+"/"
    return Ret

USAGE = """python adaptation3d.py --data-dir [datadir] --aln [aln_file] --tree [tree_file] --out-dir [output_dir]

Where:
datadir: the absolute path to the directory containing this script
aln_file: the path to the input alignment file
tree_file: the path to the input tree file
output_dir: the path to the output directory to be written
"""

#performs all necessary checks, prints the USAGE string if there is incorrect input
if len(sys.argv) != 9:
    print USAGE
else:
    
    if sys.argv[1] == "--data-dir" and sys.argv[3] == "--aln" and sys.argv[5] == "--tree" and sys.argv[7] == "--out-dir":
        if os.path.exists(sys.argv[4]) and os.path.exists(sys.argv[6]):
            if os.path.exists(addSlashIfNeeded(sys.argv[2])+"src/wholetreegenescope.py"):
                
                #adds the source directory to the path
                sys.path.append(addSlashIfNeeded(sys.argv[2])+"src")
                
                #from wholetreegenescope import WholeTreeGeneScope
                from scopealgorithmtreeset import ScopeAlgorithmTreeSet
                #from explorepredictionwrapper import ExplorePredictionWrapper
                
                datadir = addSlashIfNeeded(sys.argv[2])
                
                #instantiates instances of all classes needed for the scope algorithm
                
                #if os.path.exists(addSlashIfNeeded(sys.argv[8])+"Report.xml"):
                #    pass
                #else:
                #WholeTreeGeneScope(sys.argv[4] , sys.argv[6] , sys.argv[8])
                
                
                if os.path.exists(addSlashIfNeeded(sys.argv[8])+"PValues.txt"):
                    pass
                else:
                    ScopeAlgorithmTreeSet(sys.argv[8] , datadir)
                
                """
                if os.path.exists(addSlashIfNeeded(sys.argv[8])+"Figures"):
                    pass
                else:
                    ExplorePredictionWrapper(sys.argv[8])
                """
            
            else:
                print USAGE
        else:
            print USAGE
    else:
        print USAGE
