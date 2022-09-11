# MDAnalysis contact analysis
import numpy
import subprocess
from os import listdir
from os.path import isfile, join

##Class types - always very useful to have for speedy documentation lookup
#AG == "MDAnalysis.core.groups.AtomGroup"
#TS == "MDAnalysis.coordinates.base.Timestep"
#Array == "numpy.ndarray"
#DF == pandas dataframe

#Constants
MY_PATH_STR: str  = "/mnt/d/DE_NOVO_PROTEINS/PROTEIN_ACCESSIBILITY/"    

XTC_FILE_STR: str = "xtc"
TPR_FILE_STR: str = "tpr"
SASA_STR: str = "SASA" 
SASA_INPUT_FILE_STR: str = "sasa_input"
    
ALL_MONOMERS: str = "all_monomers"
MONOMERS_DIMER: str = "monomers_dimer"
TWO_DIMERS: str = "two_dimers"
MONOMER_TRIMER: str = "monomer_trimer"
FOURMER: str = "4mer"
    
    
sim0List = [[150,310],[700,1001]]  
sim1List = [[150,300],[380,850]]
sim11List = [[150,300],[430,750]]
sim12List = [[150,600],[900,1001]]
sim13List = [[150,300],[540,1001]]
sim14List = [[280,1001]]
sim15List = [[350,950]]
sim16List = [[150,230],[400,1001]]
sim17List = [[150,670],[720,830],[930,1001]]
sim18List = [[150,400],[600,1001]]
sim19List = [[200,680],[830,1001]]
sim2List = [[150,580],[682,790],[800,1001]]
sim20List = [[150,360],[460,630],[720,1001]]
sim21List = [[260,1001]]
sim22List = [[200,280],[340,1001]]
sim23List = [[250,1001]]
sim24List = [[150,600],[760,1001]]
sim25List = [[150,300],[520,1001]]
sim26List = [[150,800],[910,1001]]
sim27List = [[150,850]]
sim28List = [[250,500],[640,1001]]
sim29List = [[250,700],[860,1001]]
sim3List = [[150,320],[390,600],[880,1001]]
sim30List = [[150,600],[800,1001]]
sim31List = [[150,600],[770,1001]]
sim32List = [[150,270],[320,500],[800,1001]]
sim33List = [[300,700],[880,1001]]
sim34List = [[200,1001]]
sim35List = [[400,800]]
sim36List = [[150,500],[800,1001]]
sim37List = [[240,550],[780,1001]]
sim4List = [[150,500],[800,1001]]
sim5List = [[150,1001]]
sim6List = [[650,1001]]
sim7List = [[150,320],[570,1001]]
sim9List = [[700,1001]]
    
frameList = [sim0List,
            sim1List,
            sim11List,
            sim12List,
            sim13List,
            sim14List,
            sim15List,
            sim16List,
            sim17List,
            sim18List,
            sim19List,
            sim2List,
            sim20List,
            sim21List,
            sim22List,
            sim23List,
            sim24List,
            sim25List,
            sim26List,
            sim27List,
            sim28List,
            sim29List,
            sim3List,
            sim30List,
            sim31List,
            sim32List,
            sim33List,
            sim34List,
            sim35List,
            sim36List,
            sim37List,
            sim4List,
            sim5List,
            sim6List,
            sim7List,
            sim9List]

oligoColumnList = [ALL_MONOMERS, MONOMERS_DIMER, TWO_DIMERS, MONOMER_TRIMER, FOURMER]

oligo0List = [ALL_MONOMERS,FOURMER]
oligo1List = [ALL_MONOMERS,MONOMERS_DIMER]
oligo11List = [ALL_MONOMERS,MONOMERS_DIMER]
oligo12List = [ALL_MONOMERS,MONOMER_TRIMER]
oligo13List = [ALL_MONOMERS,MONOMER_TRIMER]
oligo14List = [MONOMERS_DIMER]
oligo15List = [MONOMERS_DIMER]
oligo16List = [ALL_MONOMERS,TWO_DIMERS]
oligo17List = [ALL_MONOMERS,MONOMERS_DIMER,MONOMER_TRIMER]
oligo18List = [ALL_MONOMERS,TWO_DIMERS]
oligo19List = [MONOMERS_DIMER,MONOMER_TRIMER]
oligo2List = [ALL_MONOMERS,MONOMERS_DIMER,TWO_DIMERS]
oligo20List = [ALL_MONOMERS, MONOMERS_DIMER, MONOMER_TRIMER]
oligo21List = [MONOMERS_DIMER]
oligo22List = [MONOMERS_DIMER, TWO_DIMERS]
oligo23List = [MONOMERS_DIMER]
oligo24List = [ALL_MONOMERS, MONOMERS_DIMER]
oligo25List = [ALL_MONOMERS, MONOMERS_DIMER]
oligo26List = [ALL_MONOMERS, MONOMERS_DIMER]
oligo27List = [ALL_MONOMERS]
oligo28List = [MONOMERS_DIMER, MONOMER_TRIMER]
oligo29List = [MONOMERS_DIMER, MONOMER_TRIMER]
oligo3List = [ALL_MONOMERS,MONOMERS_DIMER,FOURMER]
oligo30List = [ALL_MONOMERS, MONOMERS_DIMER]
oligo31List = [ALL_MONOMERS, MONOMERS_DIMER]
oligo32List = [ALL_MONOMERS, MONOMERS_DIMER, MONOMER_TRIMER]
oligo33List = [MONOMERS_DIMER,FOURMER]
oligo34List = [MONOMERS_DIMER]
oligo35List = [MONOMERS_DIMER]
oligo36List = [ALL_MONOMERS, MONOMERS_DIMER]
oligo37List = [MONOMERS_DIMER, MONOMER_TRIMER]
oligo4List = [ALL_MONOMERS,TWO_DIMERS]
oligo5List = [ALL_MONOMERS]
oligo6List = [MONOMER_TRIMER]
oligo7List = [ALL_MONOMERS,MONOMER_TRIMER]
oligo9List = [TWO_DIMERS]

oligoList = [oligo0List,
            oligo1List,
            oligo11List,
            oligo12List,
            oligo13List,
            oligo14List,
            oligo15List,
            oligo16List,
            oligo17List,
            oligo18List,
            oligo19List,
            oligo2List,
            oligo20List,
            oligo21List,
            oligo22List,
            oligo23List,
            oligo24List,
            oligo25List,
            oligo26List,
            oligo27List,
            oligo28List,
            oligo29List,
            oligo3List,
            oligo30List,
            oligo31List,
            oligo32List,
            oligo33List,
            oligo34List,
            oligo35List,
            oligo36List,
            oligo37List,
            oligo4List,
            oligo5List,
            oligo6List,
            oligo7List,
            oligo9List]

fileList = [f for f in listdir(MY_PATH_STR) if isfile(join(MY_PATH_STR, f))]

xtcFileList = []
tprFileList = []
for fileStr in fileList:
    if fileStr.endswith(XTC_FILE_STR):
        xtcFileList.append(join(MY_PATH_STR, fileStr))
    elif fileStr.endswith(TPR_FILE_STR):
        tprFileList.append(join(MY_PATH_STR, fileStr))
        
for i in range(len(xtcFileList)):
    xtcFileStr: str = xtcFileList[i]
    tprFileStr: str = tprFileList[i]

    oligoSetList = oligoList[i]
    oligoCountInt: int = 0
    for frameSetList in frameList[i]:
        oligoStateStr: str = oligoSetList[oligoCountInt]
        startFrameStr: str = str(frameSetList[0] * 50)
        endFrameStr: str = str(frameSetList[1] * 50)
        outputFileStr: str = MY_PATH_STR + SASA_STR + "_" + oligoStateStr + "_" + str(i)
        oligoCountInt += 1


        commandList = ["gmx", "sasa", "-f", xtcFileStr, "-s", tprFileStr, "-o", outputFileStr, "-b", startFrameStr, "-e", endFrameStr]

        with open(SASA_INPUT_FILE_STR, "r") as sasaInpFileObj:
            result = subprocess.run(commandList, stdin=sasaInpFileObj)
        sasaInpFileObj.close()
