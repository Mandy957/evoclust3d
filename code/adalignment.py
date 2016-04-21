####################################################################################################
#                                                                                                  #
#           PROJECT Protein Adaptation                                                             #
#           CLASS ADAlignment                                                                      #
#           PROGRAMMER Jeremy Adams                                                                #
#           STARTED 25-07-13                                                                       #
#           LASTMOD 12-08-14                                                                       #
#                                                                                                  #
#           DESCRIPTION class representation of the aligned ancestral and derived sequences. Used  #
#                       to obtain mutation positions between the sequences                         #
#                                                                                                  #
####################################################################################################

class ADAlignment:
    
    "CONSTRUCTOR"
    def __init__(self , FromSeq_S , ToSeq_S , FromProb_L , ToProb_L):
        
        """
        Class attributes:
        FromSeq_S (String): string containing the sequence of the ancestral state
        ToSeq_S (String): string containing the sequence of the derived state
        FromProb_L (List): list of probabilities of reconstructed states for the ancestral sequence, one probability per state in the sequence
        ToProb_L (List): list of probabilities of reconstructed states for the derived sequence
        """
        
        self.FromSeq_S = FromSeq_S
        self.ToSeq_S = ToSeq_S
        self.FromProb_L = FromProb_L
        self.FromProb_L = FromProb_L
        self.ToProb_L = ToProb_L
        
        #finds the mutations between the ancestral and derived
        self.mutations = self.findMutations()
    
    "builds a list wherein each element represents where a mutation has occurred"
    def findMutations(self):
        
        mutation_string = []
        #for each state in the aligned sequence pair
        for i in range(0,len(self.FromSeq_S)):
            
            #do nothing if there is a gap or a non-standard residue
            if self.FromSeq_S[i] == "X" or self.FromSeq_S[i] == "-" or self.FromSeq_S[i] == "*" or self.ToSeq_S[i] == "X" or self.ToSeq_S[i] == "-" or self.ToSeq_S[i] == "*":
                pass
            else:
                #if there is an amino acid mismatch, an element containing the position, amino acid types, and probabilities are added to the list
                if self.FromSeq_S[i] != self.ToSeq_S[i]:
                    mutation_string.append("%s|%s|%s|%s|%s" % (str(i+1) , self.FromSeq_S[i] , self.ToSeq_S[i] , self.FromProb_L[i] , self.ToProb_L[i]))
        
        
        ret = ""
        #builds a single string out of the whole list of mutations
        if len(mutation_string) == 0:
            ret = None
        else:
            ret = ';'.join(mutation_string)
                    
        return ret