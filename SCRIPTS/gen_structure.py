from loos import *
import loos.pyloos

from modeller import *
from modeller.optimizers import ConjugateGradients

from typing import List
import random
import subprocess
import sys
import math
import re

import MDAnalysis
from MDAnalysis.analysis import distances

STEEPEST_DESCENT = "steepest descent"
CONJUGATE_GRADIENT = "conjugate gradient"

#SHORT_EQ = "short_eq"
LONG_EQ = "long_eq"
#SHORT_EQ_REST = "short_eq_rest"
#LONG_EQ_REST = "long_eq_rest"

SHORT_EQ_STRONG = "short_eq_strong"
LONG_EQ_STRONG = "long_eq_strong"
LONG_EQ_MEDIUM = "long_eq_medium"
LONG_EQ_WEAK = "long_eq_weak"

EQ_VELOCITIES_TRUE = True
EQ_VELOCITIES_FALSE = False

#imutableAAPositionList should be 0 to num-1
def design_sequence(availableAAList: List[str], scafSequenceStr: str, imutableAAPositionList: List[int], numAAToChange: int):
    numAAInt: int = len(scafSequenceStr)
    numAvailableAAInt: int = len(availableAAList)
    tempImutableAAPositionList: List = imutableAAPositionList

    if (numAAToChange + len(tempImutableAAPositionList)) > len(scafSequenceStr):
        print(len(imutableAAPositionList))
        #print("Error: there is an inbalance between the scaffold length, the number of desired changes and the number of imutable positions.")
        return None

    newScafStructureStr: str = scafSequenceStr
    changesMadeInt: int = 0

    while changesMadeInt <= numAAToChange:
        randPosInt: int = random.randint(0,numAAInt-1)
        if randPosInt in tempImutableAAPositionList:
            continue
        else:
            newAAStr: str = availableAAList[random.randint(0, numAvailableAAInt-1)]
            newScafStructureStr = newScafStructureStr[:randPosInt] + newAAStr + newScafStructureStr[randPosInt+1:]
            tempImutableAAPositionList.append(randPosInt)
            changesMadeInt = changesMadeInt + 1

    return newScafStructureStr

#this returns the charge on ONE peptide (including the terminal residues)
def calculate_peptide_charge(sequenceStr: str):
    sequenceChargeInt: int = 0
    for aa in sequenceStr:
        if aa == "R" or aa == "K":
            sequenceChargeInt = sequenceChargeInt + 1
        elif aa == "D" or aa == "E":
            sequenceChargeInt = sequenceChargeInt - 1
    return sequenceChargeInt


def test_hydrophobicity(sequenceStr: str, supressOutputBool: bool):
    result = subprocess.run(["perl", "analyze.pl", sequenceStr, "-o", "/dev/stdout"], capture_output=True, text=True)
    deltaGStr: str = result.stdout

    if deltaGStr[0] == "+":
        deltaGStr = deltaGStr[1:len(deltaGStr)-1]

    if deltaGStr == "--0.0":
        deltaGStr = "0.0"

    deltaGFlt = float(deltaGStr)

    if isinstance(deltaGFlt, float):
        return float(deltaGStr)
    else:
        raise Exception("Error assigned free energy prediction of TM insertion.")

#Takes a sequence and file name and outputs an alpha helix of that sequence design
def build_helix(sequenceStr: str, dirLocationStr: str):

    fileNameStr: str = dirLocationStr + "helix.pdb"
    extendedChainStr: str = dirLocationStr + "extended-chain.pdb"
    rotatedFileNameStr: str = dirLocationStr + "helix_rotated.pdb"

    environment = Environ()
    environment.libs.topology.read("${LIB}/top_heav.lib")
    environment.libs.parameters.read("${LIB}/par.lib")

    model = Model(environment)
    model.build_sequence(sequenceStr)
    model.write(file=extendedChainStr)

    allatoms = Selection(model)
    model.restraints.make(allatoms, restraint_type="STEREO", spline_on_site=False)
    model.restraints.add(secondary_structure.Alpha(model.residues))

    cg = ConjugateGradients()
    cg.optimize(allatoms, max_iterations=100)

    model.write(file=fileNameStr)

    #use LOOS to orientate the helix
    modelAtomicGroup = loos.createSystem(fileNameStr)
    modelAtomicGroup.rotate(GCoord(0,1,0), 90)
    modelAtomicGroup.rotate(GCoord(1,0,0), -45)
    pdb = PDB.fromAtomicGroup(modelAtomicGroup)
    
    #this needs to write out to a file
    originalSTDOUT = sys.stdout
    fileObj = open(rotatedFileNameStr, "w")
    sys.stdout = fileObj
    print(pdb)
    sys.stdout = originalSTDOUT
    fileObj.close()


#execute martinize to build a CG representation of a helix and fix the dihedral angles
def build_martinize(inputFileStr: str, topFileStr: str, outputPDBFileStr: str, itpFileStr: str, supressOutputBool):
    result = subprocess.run(["/home/ubuntu/DE_NOVO_PROTEINS/SCRIPTS/martini_shell.sh", inputFileStr, topFileStr, outputPDBFileStr, itpFileStr], capture_output=True, text=True)

    if supressOutputBool == False:
        print("stdout", result.stdout)
        print("stderr", result.stderr)

def moveHelix(cgFileStr, outputFileStr, xDestFlt, yDestFlt, zDestFlt, supressOutputBool):
    #Use loos to build the model and get the x, y, z position of the incoming helix
    xCurrentFlt = 0.0
    yCurrentFlt = 0.0
    zCurrentFlt = 0.0

    modelAtomicGroup = loos.createSystem(cgFileStr)
    residuesVec = modelAtomicGroup.splitByResidue()
    res12AG = residuesVec[11]
    for atom in res12AG:
        atomName = atom.name()
        if atomName == "BB":
            coordinates = atom.coords()
            xCurrentFlt = coordinates.x()
            yCurrentFlt = coordinates.y()
            zCurrentFlt = coordinates.z()

    if xCurrentFlt == 0.0 and yCurrentFlt == 0.0 and zCurrentFlt == 0.0:
        raise Exception("Error: could not find the current coordinates in original helix file.")
    else:
        #gromacs requires nanometers, loos reads in the PDB angstrom numbers
        xDistanceFlt = xDestFlt - (xCurrentFlt/10)
        yDistanceFlt = yDestFlt - (yCurrentFlt/10)
        zDistanceFlt = zDestFlt - (zCurrentFlt/10)

        translateValueStr = str(xDistanceFlt) + " " + str(yDistanceFlt) + " " + str(zDistanceFlt)

        translateCommandList = ["gmx", "editconf", "-f", cgFileStr, "-translate", str(xDistanceFlt), str(yDistanceFlt), str(zDistanceFlt), "-o", outputFileStr]
        result = subprocess.run(translateCommandList, capture_output = True, text=True)

        if supressOutputBool == False:
            print("stdout:", result.stdout)
            print("stderr:", result.stderr)

#The chargeInt values will have already been multiplied by 4 (due to 4 peptides being present)
def write_top_file(topFileStr: str, itpFileStr: str, sequenceStr: str, restraintsDirStr: str, chargeInt: int):

    clCountInt: int = 126
    naCountInt: int = 126

    if chargeInt != 16: #+4 charge on each polyleu because of the charged termini
        if chargeInt < 0:
            clCountInt = clCountInt - chargeInt
        elif chargeInt > 0:
            naCountInt = naCountInt - chargeInt

    else:
        naCountInt = naCountInt - 16

    with open(topFileStr, "w") as fileObj:
        fileObj.write("#include \"../../martini_3.0.b.3.2/martini_v3.0.b.3.2.itp\"\n")
        fileObj.write("#include \"../../martini_3.0.b.3.2/DPPC_M3.itp\"\n")
        fileObj.write("#include \"../../martini_3.0.b.3.2/POPC_M3.itp\"\n")
        fileObj.write("#include \"../../martini_3.0.b.3.2/martini_v3.0_solvents.itp\"\n")
        #the next line adds the atoms for a single protein
        fileObj.write("#include \"" + itpFileStr + "\"\n")
        #if withRestraintsBool == True:
        fileObj.write("#ifdef STRONG\n")
        fileObj.write("#include \"" + restraintsDirStr + "/STRONG_posre.itp" + "\"\n")
        fileObj.write("#endif\n")
   
        fileObj.write("#ifdef MEDIUM\n")
        fileObj.write("#include \"" + restraintsDirStr + "/MEDIUM_posre.itp" + "\"\n")
        fileObj.write("#endif\n")

        fileObj.write("#ifdef WEAK\n")
        fileObj.write("#include \"" + restraintsDirStr + "/WEAK_posre.itp" + "\"\n")
        fileObj.write("#endif\n")
        fileObj.write("#include \"../../martini_3.0.b.3.2/martini_v3.0_ions.itp\"\n")
        fileObj.write("\n")
        fileObj.write("[ system ]\n")
        fileObj.write("; name\n")
        fileObj.write(sequenceStr + "\n")
        fileObj.write("\n")
        fileObj.write("[ molecules ]\n")
        fileObj.write("; name  number\n")
        fileObj.write("DPPC             132\n")
        fileObj.write("POPC             102\n")
        fileObj.write("WN               5860\n")
        fileObj.write("NA               " + str(naCountInt) + "\n")
        fileObj.write("CL               " + str(clCountInt) + "\n")
        fileObj.write("Protein_A        4")
    fileObj.close()

    return naCountInt, clCountInt

#this is for the purpose of testing the helix via an energy minimisation
def write_helix_top_file(topFileStr: str, itpFileStr: str, sequenceStr: str):
    with open(topFileStr, "w") as fileObj:
        fileObj.write("#include \"../../martini_3.0.b.3.2/martini_v3.0.b.3.2.itp\"\n")
        fileObj.write("#include \"" + itpFileStr + "\"\n")
        fileObj.write("\n")
        fileObj.write("[ system ]\n")
        fileObj.write("; name\n")
        fileObj.write(sequenceStr + "\n")
        fileObj.write("\n")
        fileObj.write("[ molecules ]\n")
        fileObj.write("; name  number\n")
        fileObj.write("Protein_A        1")
    fileObj.close()


#takes in a gro file and changes the number of ions depending on the charge
#numPositiveChargeInt and numNegativeChargeInt are the number of positive and negative ions to keep
def adjustCharge(inputFileStr: str, outputFileStr:str , numPositiveChargeInt: int, numNegativeChargeInt):

    #I don't expect this to happen, but if the NA/CL count are both 126 then nothing happens (as that's 
    #the base line number of NA/CL ions in the solution at 0.15M)
    if numPositiveChargeInt == 126 and numNegativeChargeInt == 126:
        return None

    titleStr: str = ""
    numAtomsInt: int = 0
    contentList: List = []
    newContentList: List = []
    counterInt: int = 0
    #read the contents into a list and store the title and numberr of atoms
    with open(inputFileStr) as fileObj:
        for lineStr in fileObj:
            if counterInt == 0 :
                titleStr = lineStr
            elif counterInt == 1:
                numAtomsInt = int(lineStr)
            else:
                contentList.append(lineStr)
            counterInt = counterInt + 1
    fileObj.close()

    #go through the contents andadjust the number of ions to store
    clIonCounterInt: int = 0
    naIonCounterInt: int = 0
    for newLine in contentList:
        splitStrList: List = newLine.split()
        clMatchObj = re.search("\d+CL", splitStrList[0])
        naMatchObj = re.search("\d+NA", splitStrList[0])

        if clMatchObj is not None:
            if clIonCounterInt < numNegativeChargeInt:
                newContentList.append(newLine)
                clIonCounterInt = clIonCounterInt + 1
        elif naMatchObj is not None:
            if naIonCounterInt < numPositiveChargeInt:
                newContentList.append(newLine)
                naIonCounterInt = naIonCounterInt + 1
        else:
            #all other compounds (lipid, water, protein, etc) and also the cell dimensions
            newContentList.append(newLine)
    numAtoms = len(newContentList)-1
    with open(outputFileStr, "w") as fileObj:
        fileObj.write(titleStr)
        fileObj.write(str(numAtoms) + "\n")
        for newLineStr in newContentList:
            fileObj.write(newLineStr)
    fileObj.close()
    

def write_gro_file(inputFileList, outputFileStr: str):
    #read in bilayer file
    topCounterInt: int = 0
    bilayerList = []
    with open(inputFileList[0], "r") as fileObj:
        for lineStr in fileObj:
            if topCounterInt < 2:
                topCounterInt = topCounterInt + 1
                continue
            else:
                bilayerList.append(lineStr)
    fileObj.close()

    #read in helix 0
    topCounterInt: int = 0
    helix0List = []
    with open(inputFileList[1], "r") as fileObj:
        for lineStr in fileObj:
            if topCounterInt < 2:
                topCounterInt = topCounterInt + 1
                continue
            else:
                helix0List.append(lineStr)
    fileObj.close()

    #read in helix 1
    topCounterInt: int = 0
    helix1List = []
    with open(inputFileList[2], "r") as fileObj:
        for lineStr in fileObj:
            if topCounterInt < 2:
                topCounterInt = topCounterInt + 1
                continue
            else:
                helix1List.append(lineStr)
    fileObj.close()

    #read in helix 2
    topCounterInt: int = 0
    helix2List = []
    with open(inputFileList[3], "r") as fileObj:
        for lineStr in fileObj:
            if topCounterInt < 2:
                topCounterInt = topCounterInt + 1
                continue
            else:
                helix2List.append(lineStr)
    fileObj.close()

    #read in helix 3
    topCounterInt: int = 0
    helix3List = []
    with open(inputFileList[4], "r") as fileObj:
        for lineStr in fileObj:
            if topCounterInt < 2:
                topCounterInt = topCounterInt + 1
                continue
            else:
                helix3List.append(lineStr)
    fileObj.close()

    numAtomsInt: int = len(bilayerList) + 4*(len(helix0List)) - 5

    with open(outputFileStr, "w") as fileObj:
        fileObj.write("Combined GRO file\n")
        fileObj.write(str(numAtomsInt) + "\n")
        for i in range(len(bilayerList)-1):
            lineStr: str = bilayerList[i]
            fileObj.write(lineStr)

        for i in range(len(helix0List)-1):
            lineStr: str = helix0List[i]
            fileObj.write(lineStr)

        for i in range(len(helix1List)-1):
            lineStr: str = helix1List[i]
            fileObj.write(lineStr)

        for i in range(len(helix2List)-1):
            lineStr: str = helix2List[i]
            fileObj.write(lineStr)

        for i in range(len(helix3List)-1):
            lineStr: str = helix3List[i]
            fileObj.write(lineStr)

        fileObj.write(bilayerList[len(bilayerList)-1])
    fileObj.close()

#to gromp all files #I am not sure if I need to ever strain the peptides?
def run_grompp(mdpFileStr: str, groFileStr: str, topFileStr: str, outputFileStr: str, supressOutputBool: bool, eqVelocitiesBool, epocDirStr: str, velocityFileStr: str, restraintsBool: bool):
    commandList = ["gmx","grompp","-f",mdpFileStr,"-c",groFileStr,"-p",topFileStr,"-o",outputFileStr, "-maxwarn","1"]

    if restraintsBool == True:
        commandList.append("-r")
        commandList.append(epocDirStr + "/minimization.gro")


    if eqVelocitiesBool == EQ_VELOCITIES_TRUE:
        commandList.append("-t")
        commandList.append(velocityFileStr)

    result = subprocess.run(commandList, capture_output=True, text=True)
    if supressOutputBool == False:
        print("stdout:", result.stdout)
        print("stderr:", result.stderr)

#run mdrun 
def run_simulation(deffnmFileStr: str, inputFileTprStr: str, gpuBool: bool, numMPI: int, numOMP: int, supressOutputBool: bool):
    commandList = ["gmx", "mdrun", "-s", inputFileTprStr, "-deffnm", deffnmFileStr, "-ntmpi", str(numMPI), "-ntomp", str(numOMP)]

    if gpuBool == True:
        commandList.append("-nb")
        commandList.append("gpu")

    result = subprocess.run(commandList, capture_output=True, text=True)
    if supressOutputBool == False:
        print("stdout:" + result.stdout)
        print("stderr:" + result.stderr)


#make a basic index file (only save and quit)
def run_make_index(inputFileStr: str, outputFileStr: str, indexInpFileStr: str, supressOutputBool: bool):
    commandList = ["gmx", "make_ndx", "-f", inputFileStr, "-o", outputFileStr]
    with open(indexInpFileStr, "r") as indexInpFileObj:
        result = subprocess.run(commandList, stdin=indexInpFileObj)
        if supressOutputBool == False:
            print("stderr:", result.stderr)
            print("stdout:", result.stdout)
    indexInpFileObj.close()


def convert_gro_to_pdb(inputGROFileStr: str, outputPDBFileStr: str, supressOutputBool: bool, translationList = None):
    commandList = ["gmx", "editconf", "-f", inputGROFileStr, "-o", outputPDBFileStr]

    if translationList != None:
        commandList.append("-translate")
        commandList.append(str(translationList[0]))
        commandList.append(str(translationList[1]))
        commandList.append(str(translationList[2]))

    result = subprocess.run(commandList, capture_output=True, text=True)
    if supressOutputBool == False:
        print("stdout:" + result.stdout)
        print("stderr:" + result.stderr) 


def convert_pdb_to_gro(inputPDBFileStr: str, outputGROFileStr: str, supressOutputBool: bool, translationList = None):
    commandList = ["gmx", "editconf", "-f", inputPDBFileStr, "-o", outputGROFileStr]

    if translationList != None:
        commandList.append("-translate")
        commandList.append(str(translationList[0]))
        commandList.append(str(translationList[1]))
        commandList.append(str(translationList[2]))

    result = subprocess.run(commandList, capture_output=True, text=True)
    if supressOutputBool == False:
        print("stdout:" + result.stdout)
        print("stderr:" + result.stderr)


def calculate_helix_distance(helixActualFileStr: str, helixReferenceFileStr: str):
    u1 = MDAnalysis.Universe(helixActualFileStr)
    u2 = MDAnalysis.Universe(helixReferenceFileStr)

    ca1 = u1.select_atoms("name BB")
    ca2 = u2.select_atoms("name BB")

    ca1COM = ca1.center_of_mass()
    ca2COM = ca2.center_of_mass()

    #distance = distances.distance_array(ca1COM, ca2COM, box=u1.dimensions)
    actualX = ca1COM[0]
    actualY = ca1COM[1]
    actualZ = ca1COM[2]

    referenceX = ca2COM[0]
    referenceY = ca2COM[1]
    referenceZ = ca2COM[2]

    #turn it from angstrom into nano metres
    deltaX = (referenceX - actualX)*10
    deltaY = (referenceY - actualY)*10
    #deltaZ = (referenceZ - actualZ)*10
    deltaZ = 0

    returnList = [deltaX, deltaY, deltaZ]
    return(returnList)
    






