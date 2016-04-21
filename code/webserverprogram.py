import os
import re
import sys
import random
from staticmethods import *
from mutationtostructure import MutationToStructure
from scopealgorithmtriosset import ScopeAlgorithmTriosSet
from scipy.stats import percentileofscore

class WebServerProgram:
    
    def __init__(self, interest , ingroup , outgroup , cusOrOma):
        self.WSDIR = "/home/j7adams/Adaptation3D/"
        self.interest = interest
        self.ingroup = ingroup
        self.outgroup = outgroup
        self.cusOrOma = cusOrOma
        self.ProjDIR = ""
        if self.cusOrOma == "oma":
            self.ProjDIR = self.WSDIR+"projects/oma/"+"-".join([self.interest,self.ingroup,self.outgroup])
        elif self.cusOrOma == "custom":
            self.ProjDIR = self.WSDIR+"projects/custom/"+"-".join([self.interest,self.ingroup,self.outgroup])
        self.SAS_D = getSASD()
        
        self.makeDirIfNeeded(self.ProjDIR)
        OG = None
        
        if self.cusOrOma == "oma":
            OG = self.getPerfectOrthologyGroupsForOma()
        elif self.cusOrOma == "custom":
            OG = self.getPerfectOrthologyGroupsForCustom()
        
        
        self.GroupKey_L = OG[0]
        self.Groups_D = OG[1]
        NeededProtein_D = OG[2]
        
        self.SpeciesToFasta_D = {sp : self.readSpeciesFasta(sp , NeededProtein_D[sp]) for sp in [self.interest,self.ingroup,self.outgroup]}
        
        self.makeDirIfNeeded(self.ProjDIR+"/Alignments/")
        self.getAlignments()
        self.makeDirIfNeeded(self.ProjDIR+"/Reports/")
        self.makeDirIfNeeded(self.ProjDIR+"/Tables/")
        self.makeDirIfNeeded(self.ProjDIR+"/ScoringMatrices/")
        self.makeDirIfNeeded(self.ProjDIR+"/PValues/")
        self.runMutationToStructure()
        
        #self.makeDirIfNeeded(self.ProjDIR+"/Summary/")
        #self.getDatasetSummary()
        #self.getRandomSummary()
        
        
    def makeDirIfNeeded(self,DIR):
        if os.path.exists(DIR):
            pass
        else:
            os.mkdir(DIR)
    
    def getPerfectOrthologyGroupsForOma(self):
        L = [self.interest , self.ingroup , self.outgroup]
        S = set(L)
        
        Ret_D = {sp : [] for sp in L}
        Groups_D = {}
        GroupKey_L = []
        with open(self.WSDIR+"oma/raw/oma-groups.txt" , "r") as GroupF:
            for line in GroupF:
                D = {sp : [] for sp in L}
                ls = line.split()
                Group = ls[0]
                
                for seq in ls[1:]:
                    species = seq[:5]
                    if species in S:
                        D[species].append(seq)
                
                if len(D[self.interest]) == 1 and len(D[self.ingroup]) == 1 and len(D[self.outgroup]) == 1:
                    GroupKey_L.append(Group)
                    Groups_D[Group] = ([D[self.interest][0] , D[self.ingroup][0] , D[self.outgroup][0]])
                    
                    for sp in L:
                        Ret_D[sp].append(D[sp][0])
        
        return [GroupKey_L , Groups_D , Ret_D]
    
    def getPerfectOrthologyGroupsForCustom(self):
        Paths_L = [self.WSDIR+"custom/sqltables/sqltable.%s-%s"% (self.interest,self.ingroup) , self.WSDIR+"custom/sqltables/sqltable.%s-%s"% (self.interest,self.outgroup)]
        D = {}
        
        for Path in Paths_L:
            Specs = re.compile("sqltable\.(.+?)$").search(Path).group(1).split("-")
            
            Res = parseInparanoidTable(Path , Specs[0] , Specs[1])
            D["%s-%s" % (Specs[0] , Specs[1])] = Res
        
        taken_D = {i : False for i in range(0,len(D["%s-%s" % (self.interest,self.outgroup)]))}
        
        for b in range(0,len(D["%s-%s" % (self.interest,self.ingroup)])):
            NcoriicepsFirstHeaders = set(D["%s-%s" % (self.interest,self.ingroup)][b][self.interest])
            
            MatchNotFound = True
            for a in range(0,len(D["%s-%s" % (self.interest,self.outgroup)])):
                if MatchNotFound:
                    secondline = D["%s-%s" % (self.interest,self.outgroup)][a]
                
                    NcoriicepsSecondHeaders = secondline[self.interest]
                    
                    ThisWorks = False
                
                    for secondHeader in NcoriicepsSecondHeaders:
                        if secondHeader in NcoriicepsFirstHeaders:
                            ThisWorks = True
                    
                    if ThisWorks:
                        if taken_D[a] == False:
                            MatchNotFound = False
                            taken_D[a] = True
                            D["%s-%s" % (self.interest,self.ingroup)][b][self.outgroup] = D["%s-%s" % (self.interest,self.outgroup)][a][self.outgroup]
        
        GroupKey_L = []
        Groups_D = {}
        Ret_D = {sp : [] for sp in [self.interest,self.ingroup,self.outgroup]}
        i = 1
        L = [self.interest,self.ingroup,self.outgroup]
        for line in D["%s-%s" % (self.interest,self.ingroup)]:
            if self.interest in line.keys() and self.ingroup in line.keys() and self.outgroup in line.keys():
                if len(line[self.interest]) == 1 and len(line[self.ingroup]) == 1 and len(line[self.outgroup]) == 1:
                    GroupKey_L.append(str(i))
                    Groups_D[str(i)] = [line[self.interest][0] , line[self.ingroup][0] , line[self.outgroup][0]]
                    Ret_D[self.interest].append(line[self.interest][0])
                    Ret_D[self.ingroup].append(line[self.ingroup][0])
                    Ret_D[self.outgroup].append(line[self.outgroup][0])
                    
                    i += 1
        
        return [GroupKey_L , Groups_D , Ret_D]    
        
    def readSpeciesFasta(self,species , NeededProtein_L):
        PATH = self.WSDIR+self.cusOrOma+"/proteins/"+species
        S = set(NeededProtein_L)
        
        Ret = {}
        AllLinesInOneSequence = []
        
        with open(PATH , "r") as f:
            AllLinesInOneSequence.append(f.readline().replace("\n",""))
            
            for line in f:
                
                if line.startswith(">"):
                    
                    HeaderLine = AllLinesInOneSequence[0].replace(">","")
                    if self.cusOrOma == "custom":
                        HeaderLine = HeaderLine.split()[0]
                    
                    if HeaderLine in S:
                        Ret[HeaderLine] = "".join(AllLinesInOneSequence[1:])
                        
                    AllLinesInOneSequence = []
                    AllLinesInOneSequence.append(line.replace("\n",""))
                    
                else:
                    AllLinesInOneSequence.append(line.replace("\n",""))
            
            HeaderLine = AllLinesInOneSequence[0].replace(">","")
                    
            if HeaderLine in S:
                Ret[HeaderLine] = "".join(AllLinesInOneSequence[1:])
        
        return Ret
    
    def getAlignments(self):
        for GroupKey in self.GroupKey_L[0:1]:
            interestHeader = self.Groups_D[GroupKey][0]
            interestSeq = self.SpeciesToFasta_D[self.interest][interestHeader]
            
            ingroupHeader = self.Groups_D[GroupKey][1]
            ingroupSeq = self.SpeciesToFasta_D[self.ingroup][ingroupHeader]
            
            outgroupHeader = self.Groups_D[GroupKey][2]
            outgroupSeq = self.SpeciesToFasta_D[self.outgroup][outgroupHeader]
            
            incontent = ">%s\n%s\n" % (interestHeader,interestSeq) +\
                        ">%s\n%s\n" % (ingroupHeader,ingroupSeq) +\
                        ">%s\n%s" % (outgroupHeader,outgroupSeq)
            
            inFH = getInputTempFile(incontent)
            os.system("muscle -in %s -out %s/Alignments/Group%s.aln" % (inFH.name , self.ProjDIR,GroupKey))
            inFH.close()
    
    def runMutationToStructure(self):
        P = "%s/" % (self.ProjDIR)
        AlnDIR = P+"Alignments/"
        TblDIR = P+"Tables/"
        RptDIR = P+"Reports/"
        ScmDIR = P+"ScoringMatrices/"
        i = 0
        #i = int(sys.argv[3])
        for Group in self.GroupKey_L:
            
            print "%s %s" % (str(i), Group)
                
            i += 1
            
            inPATH = AlnDIR+"Group"+Group+".aln"
            tblPATH = TblDIR+"Group"+Group+"_table"
            rptPATH = RptDIR+"Group"+Group+"_report.xml"
            scmPATH = ScmDIR+"Group"+Group+"_matrix.xml"
            
            FF = readFasta(inPATH)
            if len(FF[1][FF[0][0]]) > 7000:
                print "skip"
            else:
                
                try:
                    """
                    MTS = MutationToStructure(inPATH , self.Groups_D[Group][0] , self.Groups_D[Group][1] , self.Groups_D[Group][2])
                    ScoringMatrixXML = MTS.RepresentativeSeqScoringMatrix.scoringMatrixXMLPrint()
                    CoverageKeys_D = {re.compile("<ID>(.+?)</ID>").search(C).group(1) : re.compile("<Keys>(.+?)</Keys>").search(C).group(1).split(",") for C in re.findall("<Coverage>.+?</Coverage>",ScoringMatrixXML)}
            
                    #MTSTable = self.getMTSTable(MTS.triplealn.exit['exit_string'],CoverageKeys_D)
                    rptF = open(rptPATH , "w")
                    rptF.write(MTS.triplealn.exit['exit_string'])
                    rptF.close()
                
            
                    
                    #tableF = open(tblPATH,"w")
                    #tableF.write(MTSTable)
                    #tableF.close()
                    
            
                    
                    scmF = open(scmPATH,"w")
                    scmF.write(ScoringMatrixXML)
                    scmF.close()
                    """
                
                    pvaluePATH = self.ProjDIR+"/PValues/Group%s_PValues.txt" % (Group)
                    
                    ScopeAlgorithmTiosSet = ScopeAlgorithmTriosSet(sys.argv[2] , rptPATH , scmPATH , pvaluePATH)
                except Exception as e:
                    pass
            
            
            #i += 1
            
    
    def getMTSTable(self,mutationsXML,CoverageKeys_D):
        #print i
        indRows_L = []
        groupRows_L = []
        groupAccs_D = {}
        s = mutationsXML.split("\n")[:-1]
        #print s[0]
        PDB_L = []
        searchobj = re.compile("<PDBs>(.+?)</PDBs>").search(s[0])
        if searchobj:
            PDB_L = searchobj.group(1).split(";")
        allpdbfiledicts = getAllPDBFileDicts(PDB_L)[1]
        MM = True
        if re.compile("No Mutations").search(s[1]):
            MM = False
        
        if MM:
        
            for mut in s[1:]:
                mutinfo = re.compile("<M>(.+?)</M>").search(mut).group(1).split("|")
                mutpos = mutinfo[0]
                mutstart = mutinfo[1][0]
                mutend = mutinfo[1][1]
                
                PDBRes_L = re.compile("<R>(.+?)</R>").search(mut).group(1).split(",")
                
                if len(PDBRes_L) == 1 and PDBRes_L[0] == "NOCOVERAGE":
                    indRows_L.append(" ".join([mutpos , mutstart , mutend , "none"]))
                else:
                    
                    
                    for PDBRes in PDBRes_L:
                        acc = PDBRes.split("|")[0].lower()
                        chainpos = PDBRes.split("|")[1]
                        
                        if acc in groupAccs_D.keys():
                            pass
                        else:
                            groupAccs_D[acc] = []
                        groupAccs_D[acc].append(chainpos)
                        
                        if "XMLResidue_D" in allpdbfiledicts[acc]:
                            if chainpos in allpdbfiledicts[acc]['XMLResidue_D'].keys():    
                                pdbxml = allpdbfiledicts[acc]['XMLResidue_D'][chainpos]
                                D = self.parsePDBXMLLine(pdbxml)
                                indRows_L.append(" ".join([mutpos , mutstart , mutend , acc , chainpos , D['t'] , D['sas'] , D['rsas'] , D['b'] , D['helix'] , D['sheet'] , D['ssbond'] , D['link'] , D['phi'] , D['psi']]))
                            else:
                                indRows_L.append(" ".join([mutpos , mutstart , mutend , "none"]))
                        else:
                            indRows_L.append(" ".join([mutpos , mutstart , mutend , "none"]))
        
            for acc in groupAccs_D.keys():
                pdbxml_L = []
                for chainpos in groupAccs_D[acc]:
                    
                    
                    if "XMLResidue_D" in allpdbfiledicts[acc]:
                        if chainpos in allpdbfiledicts[acc]["XMLResidue_D"].keys():
                            pdbxml_L.append(allpdbfiledicts[acc]["XMLResidue_D"][chainpos])
                
                if "XMLResidue_D" in allpdbfiledicts[acc]:
                    if len(pdbxml_L) > 1:
                        pdbXMLD_L = [self.parsePDBXMLLine(pdbxml) for pdbxml in pdbxml_L]
                        allPos_L = [D["alphacarb"] for D in pdbXMLD_L]
                        allPosPoint_L = []
                        
                        for Pos in allPos_L:
                            pp = makePoint(Pos)
                            if pp:
                                allPosPoint_L.append(pp)
                        
                        Distances_L = getCombinatorialListOfPairwiseDistances(allPosPoint_L)
                        
                        allV_L = [D["v"] for D in pdbXMLD_L]
                        allVecPoint_L = []
                        for Vec in allV_L:
                            vv = makePoint(Vec)
                            if vv:
                                allVecPoint_L.append(vv)
                        Coplanar_L = getCombinatorialListOfCoplanarities(allVecPoint_L)
                        
                        Dismean = str(round(numpy.mean(Distances_L),3))
                        Disstdev = str(round(numpy.std(Distances_L),3))
                        
                        Anglemean = str(round(numpy.mean(Coplanar_L),3))
                        Anglestdev = str(round(numpy.std(Coplanar_L),3))
                        
                        Pvalue = "NA"
                        RandDist = self.getRelativeDistanceRandDistForIndex(allpdbfiledicts[acc]["XMLResidue_D"],CoverageKeys_D[acc],len(pdbxml_L))
                        if RandDist != None:
                            Pvalue = percentileofscore(RandDist,Dismean) / 100.0
                        
                        groupRows_L.append(" ".join([acc , str(len(pdbxml_L)) , Dismean,Disstdev,Anglemean,Anglestdev,str(Pvalue)]))
        
        ret = "#IndividualRows\n"+"\n".join(indRows_L)+"\n"+"#GroupRows\n"+"\n".join(groupRows_L)
        
        return ret
                    
    
    def parsePDBXMLLine(self,pdbxml):
        t = re.compile("<t>(.+?)</t>").search(pdbxml).group(1)
        sas = re.compile("<s>(.+?)</s>").search(pdbxml).group(1).split(";")[0]
        rsas = str(round(float(float(sas) / self.SAS_D[t]) ,4))
        b = re.compile("<B>(.+?)</B>").search(pdbxml).group(1).split(";")[2]
        
        helix = "N"
        sheet = "N"
        ssbond = "N"
        link = "N"
        
        f = re.compile("<F>(.+?)</F>").search(pdbxml).group(1).split("-")
        if f[0] == "None":
            pass
        else:
            if "HELIX" in set(f):
                helix ="Y"
            if "SHEET" in set(f):
                sheet ="Y"
            if "SSBOND" in set(f):
                ssbond ="Y"
            if "link" in set(f):
                link ="Y"
        
        phipsi = re.compile("<A>(.+?)</A>").search(pdbxml).group(1).split(";")
        phi = phipsi[0]
        psi = phipsi[1]
        
        alphacarb = re.compile("<p>(.+?)</p>").search(pdbxml).group(1).split(";")[1].split(",")
        v = re.compile("<V>(.+?)</V>").search(pdbxml).group(1).split(",")
        
        return {"t":t , "sas":sas , "rsas":rsas , "b":b , "helix":helix,"sheet":sheet,"ssbond":ssbond,"link":link ,"phi":phi,"psi":psi,"alphacarb":alphacarb,"v":v}
    
    "general method for writing/retrieving a random distribution array"
    def getAnyAverageRandomDist(self,AveragedNumbers_L):
        FinalNumbers_L = [AveragedNumber for AveragedNumber in AveragedNumbers_L if math.isnan(AveragedNumber) == False]
        return array(FinalNumbers_L)
    "gets the average of numbers in a list"
    def getAveragedData(self,Numbers_L):
        return numpy.mean([Number for Number in Numbers_L if Number != None])
    "gets a random sample of integers to be used as random indices to draw numbers for the random distributions"
    def getRandomSampleOfIntegers(self,MaxLength,Index):
        return random.sample(range(MaxLength),Index)
    "get distance random distribution for a single PDB ID and index"
    def getRelativeDistanceRandDistForIndex(self,PDBXMLResidue_D,CoverageKeys_L,length):
        Ret = None
        try:
            BuildDist = []
            for i in range(0,10000):
                sample = []
                randomKeys_L = [CoverageKeys_L[i] for i in self.getRandomSampleOfIntegers(len(CoverageKeys_L),length)]
                for Key in randomKeys_L:
                    if Key in PDBXMLResidue_D.keys():
                        sample.append(makePoint(self.parsePDBXMLLine(PDBXMLResidue_D[Key])["alphacarb"]))
                distances_L = getCombinatorialListOfPairwiseDistances(sample)
                BuildDist.append(self.getAveragedData(distances_L))
            Ret = self.getAnyAverageRandomDist(BuildDist)
        except Exception as e:
            pass
        return Ret
    
    def getDatasetSummary(self):
        DIR = self.ProjDIR+"/Tables/"
        IndHeader = "Group SeqPos From To PDB ResPos PDBRes SAS RSAS BFactor Helix Sheet SSBond Link Phi Psi"
        GroupHeader = "Group PDB nMut Dismean Disstdev Anglemean Anglestdev Pvalue"
        AllIndRows = [IndHeader]
        AllGroupRows = [GroupHeader]
        
        for F in sorted(os.listdir(DIR)):
            Group = F.replace("_table","")
            
            ToInd = True
            IndRows = []
            GroupRows = []
            allLines_L = [line.replace("\n","") for line in open(DIR+F,"r").readlines()]
            
            for line in allLines_L[1:]:
                if line.startswith("#"):
                    ToInd = False
                else:
                    if ToInd:
                        IndRows.append(line)
                    else:
                        GroupRows.append(line)
            
            for IndRow in IndRows:
                ls = IndRow.split()
                toadd = ""
                if len(ls) == 0:
                    pass
                else:
                    if len(ls) == 4:
                        toadd = Group+" "+IndRow+" "+" ".join(["NA NA NA NA NA NA NA NA NA NA NA"])
                    else:
                        toadd = Group+" "+IndRow
                    AllIndRows.append(toadd)
            
            for GroupRow in GroupRows:
                AllGroupRows.append(Group+" "+GroupRow)
        
        OutDIR = self.ProjDIR+"/Summary/"
        open(OutDIR+"DatasetInd","w").write("\n".join(AllIndRows))
        open(OutDIR+"DatasetGroup","w").write("\n".join(AllGroupRows))
                
    """
    def getRandomSummary(self):
        ScoringMatrixDIR = self.ProjDIR+"/ScoringMatrices/"
        AllCoverages_L = []
        
        for F in sorted(os.listdir(ScoringMatrixDIR)):
                scmF = open(ScoringMatrixDIR+F,"r")
                scmContent = scmF.read()
                scmF.close()
                Coverages_L = re.findall("<Coverage>.+?</Coverage>" , scmContent)
                for Coverage in Coverages_L:
                    ID = re.compile("<ID>(.+?)</ID>").search(Coverage).group(1)
                    Res_L = re.compile("<Keys>(.+?)</Keys>").search(Coverage).group(1).split(",")
                    AllCoverages_L.append([ID,Res_L])
        
        for z in range(0,1):
            
            mutationOutF = open(self.ProjDIR+"/Summary/RandomMutations"+str(z) , "w")
            mutationOutF.write("From To\n")
            mutationOutF.close()
            mutationOutF = open(self.ProjDIR+"/Summary/RandomMutations"+str(z) , "a")
            
            indOutF = open(self.ProjDIR+"/Summary/RandomInd"+str(z) , "w")
            indOutF.write("PDB ResPos PDBRes SAS RSAS BFactor Helix Sheet SSBond Link Phi Psi\n")
            indOutF.close()
            indOutF = open(self.ProjDIR+"/Summary/RandomInd"+str(z) , "a")
            
            groupOutF = open(self.ProjDIR+"/Summary/RandomGroup"+str(z) , "w")
            groupOutF.write("PDB nMut Dismean Disstdev Anglemean Anglestdev Pvalue\n")
            groupOutF.close()
            groupOutF = open(self.ProjDIR+"/Summary/RandomGroup"+str(z),"a")
            
            
            IndRows_L = []
            GroupRows_L = []
            Blosum62 = getBlosum62_D()
            
            
            DatasetIndF = open(self.ProjDIR+"/Summary/DatasetInd","r")
            for line in DatasetIndF:
                if line.endswith("Psi\n"):
                    pass
                else:
                    ls = line.split()
                    From = ls[2]
                    if From == "X":
                        pass
                    else:
                        To = getMutation(Blosum62,From)[1]
                        mutationOutF.write(From+" "+To+"\n")
                        
                        if ls[4] == "none":
                            pass
                        else:
                            IndRowNotFinished = True
                        
                            while IndRowNotFinished:
                                ThisCoverage = AllCoverages_L[random.randint(0,len(AllCoverages_L)-1)]
                                acc = ThisCoverage[0]
                                res_L = ThisCoverage[1]
                                #print len(res_L)
                                ThisRes = ThisCoverage[1][random.randint(0,len(res_L)-1)]
                                
                                allpdbfiledicts = getAllPDBFileDicts([ThisCoverage[0]])[1]
                                if 'XMLResidue_D' in allpdbfiledicts[acc].keys():
                                    if ThisRes in allpdbfiledicts[acc]['XMLResidue_D']:
                                        pdbxmlline = allpdbfiledicts[acc]['XMLResidue_D'][ThisRes]
                                        D = self.parsePDBXMLLine(pdbxmlline)
                                        outline = " ".join([acc , ThisRes , D['t'] , D['sas'] , D['rsas'] , D['b'] , D['helix'] , D['sheet'] , D['ssbond'] , D['link'] , D['phi'] , D['psi']])+"\n"
                                        indOutF.write(outline)
                                        IndRowNotFinished = False
            
            mutationOutF.close()
            indOutF.close()
            
            b = 0
            DatasetGroupF = open(self.ProjDIR+"/Summary/DatasetGroup","r")
            for line in DatasetGroupF:
                print b
                b += 1
                if line.endswith("Pvalue\n"):
                    pass
                else:
                    ls = line.split()
                    nResidues = int(ls[2])
                    
                    GroupRowNotFinished = True
                    
                    while GroupRowNotFinished:
                        ThisCoverage = AllCoverages_L[random.randint(0,len(AllCoverages_L)-1)]
                        acc = ThisCoverage[0]
                        res_L = ThisCoverage[1]
                        #print len(res_L)
                        
                        try:
                            ThisRes_L = random.sample(res_L , nResidues)
                            
                            allpdbfiledicts = getAllPDBFileDicts([ThisCoverage[0]])[1]
                            if 'XMLResidue_D' in allpdbfiledicts[acc].keys():
                                allThere = True
                                
                                for Res in ThisRes_L:
                                    if Res in allpdbfiledicts[acc]['XMLResidue_D']:
                                        pass
                                    else:
                                        allThere = False
                                
                                if allThere:
                                    pdbxml_L = [allpdbfiledicts[acc]['XMLResidue_D'][pos] for pos in ThisRes_L]
                                    
                                    pdbXMLD_L = [self.parsePDBXMLLine(pdbxml) for pdbxml in pdbxml_L]
                                    allPos_L = [D["alphacarb"] for D in pdbXMLD_L]
                                    allPosPoint_L = []
                                    
                                    for Pos in allPos_L:
                                        pp = makePoint(Pos)
                                        if pp:
                                            allPosPoint_L.append(pp)
                                    
                                    Distances_L = getCombinatorialListOfPairwiseDistances(allPosPoint_L)
                                    
                                    allV_L = [D["v"] for D in pdbXMLD_L]
                                    allVecPoint_L = []
                                    for Vec in allV_L:
                                        vv = makePoint(Vec)
                                        if vv:
                                            allVecPoint_L.append(vv)
                                    Coplanar_L = getCombinatorialListOfCoplanarities(allVecPoint_L)
                                    
                                    Dismean = str(round(numpy.mean(Distances_L),3))
                                    Disstdev = str(round(numpy.std(Distances_L),3))
                                    
                                    Anglemean = str(round(numpy.mean(Coplanar_L),3))
                                    Anglestdev = str(round(numpy.std(Coplanar_L),3))
                                    
                                    Pvalue = "NA"
                                    RandDist = self.getRelativeDistanceRandDistForIndex(allpdbfiledicts[acc]["XMLResidue_D"],res_L,len(pdbxml_L))
                                    if RandDist != None:
                                        Pvalue = percentileofscore(RandDist,Dismean) / 100.0
                                    
                                    outline = " ".join([acc , str(len(pdbxml_L)) , Dismean,Disstdev,Anglemean,Anglestdev,str(Pvalue)])+"\n"
                                    groupOutF.write(outline)
                                    GroupRowNotFinished = False
                        except Exception as e:
                            GroupRowNotFinished = True
            
            groupOutF.close()
    """
                        
                    
                        
                        
                    
        
            
        
 

#WebServerProgram("HUMAN","PANTR","GORGO","oma")
#WebServerProgram("Ncoriiceps","Mzebra","Csemilaevis","custom")