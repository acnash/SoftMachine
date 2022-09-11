#There are 26 amino acids per helix i.e., 104 amino acids/resides
#I need to run the following script first:
#gmx dump -s npt_prod.tpr | ~/miniconda3/envs/soft-machine/bin/gmxdump2pdb.pl test

import loos
import loos.pyloos
import numpy as np
import sys

NUM_PEPTIDES: int = 4

ALL_MONOMERS: str = "1 + 1 + 1 + 1"
MONOMER_TRIMER: str = "3 + 1"
TWO_DIMERS: str = "2 + 2"
DIMER_MONOMER: str = "2 + 1 + 1"
FOURMER: str = "4"

def combineAGs(agList):
    
    for i in range(len(agList)):
        if i == 0:
            combinedAG = agList[0]
        else:
            combinedAG = combinedAG.merge(agList[i])
    return(combinedAG)


#create system and load in trajectory
asaFileNameStr: str = "/mnt/d/DE_NOVO_PROTEINS/PROTEIN_ACCESSIBILITY/all_monomers_test.txt" #asa_npt_eq_29_radius14.txt"
systemFileStr: str = "/mnt/d/DE_NOVO_PROTEINS/PROTEIN_ACCESSIBILITY/npt_eq_26.psf"
systemAG = loos.createSystem(systemFileStr)
traj = loos.pyloos.Trajectory("/mnt/d/DE_NOVO_PROTEINS/PROTEIN_ACCESSIBILITY/npt_eq_26.xtc", systemAG)

#select "domains"
allChainsAG = loos.selectAtoms(systemAG, 'segid=="PRTA"')  #contains all 236 beads
#nonProteinAG = loos.selectAtoms(systemAG, 'segid!="PRTA"')

lipidAG = loos.selectAtoms(systemAG, 'resname=="DPPC" || resname=="POPC"')
lipidSplitAG = lipidAG.splitByMolecule()



#======================================================================================
#This is all to do with dividing the protein mass into each helix and finding centroid
#======================================================================================
chainsList = allChainsAG.splitByMolecule() #this splits into 136 parts

chainsListLengthInt: int = len(chainsList)
noResiduesInt: int = int(chainsListLengthInt/NUM_PEPTIDES)

#constructs N (noPeptidesInt) List entries, each holding a list of AtomicGroups that represent a single helix. 
#For exampes, this might be of length 4 (4 helices), in which each entry is a list of 34 AtomicGroups that represent the amino acids, 
helixAGList = [chainsList[x:x+noResiduesInt] for x in range(0, chainsListLengthInt, noResiduesInt)]

#use LOOs merge to merge the set of amino acid atomic groups per helix into one Atomic group. They are then stored in a list of their own. 
combinedHelixAGList = []
for i in range(NUM_PEPTIDES):
    combinedHelixAGList.append(combineAGs(helixAGList[i]))


with open(asaFileNameStr, "w") as fileObj:
    frameNumberInt: int = 0
    for frame in traj:
        if frameNumberInt >= 150 and frameNumberInt <= 650:
            contacts = 0.0
            stateValueStr = ALL_MONOMERS
            for helixAG in combinedHelixAGList:
                for lipid in lipidSplitAG:
                    contacts += helixAG.logisticContact(lipid, 10, 2, frame.periodicBox())
            fileObj.write(str(frameNumberInt) + "\t" + str(contacts) + "\t" + stateValueStr + "\n")
        frameNumberInt += 1

fileObj.close()


