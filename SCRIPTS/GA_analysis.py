#from typing import List
#import random
#import subprocess
#import sys
import numpy
import pandas
import MDAnalysis
from MDAnalysis.analysis import contacts
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy.stats import ttest_ind, shapiro, f_oneway, ranksums

#def calculate_SASA(xtcFileStr: str, tprFileStr: str, beginStr: str, endStr: str, supressOutputBool: bool):
#    #the input file to specify protein is sasa_protein.txt (a value of 6)
#    processCommandList: List = []
#    processCommandList[0] = "gmx"
#    processCommandList[1] = "sasa"
#    processCommandList[2] = "-f"
#    processCommandList[3] = xtcFileStr
#    processCommandList[4] = "-s"
#    processCommandList[5] = tprFileStr
#
#    result = subprocess.run([processCommandList, capture_output=True, text=True)
#
#    if supressOutputBool == False:
#        print("stdout", result.stdout)
#        print("stderr", result.stderr)
#    pass
NORMAL_STR = "normal"
RADIUS = 10
PROTEIN_SELECTION_STR: str = "protein"
LIPID_SELECTION_STR: str = "resname DPPC or resname POPC"
FRAME_COLLECTION_POINT_INT: int = 1 #this must change
PROTEIN_STR: str = "Protein"
LIPID_STR: str = "Lipid"
PROTEIN_LIPID_STR: str = "ProteinLipid"
CONTROL_DIR_STR: str = "CONTROLS/"

#load the control metrics
allMonomerProteinArray = numpy.loadtxt(CONTROL_DIR_STR + "allMonomerProteinArray.txt")
monomerDimerProteinArray = numpy.loadtxt(CONTROL_DIR_STR + "monomerDimerProteinArray.txt")
twoDimersProteinArray = numpy.loadtxt(CONTROL_DIR_STR + "twoDimersProteinArray.txt")
monomerTrimerProteinArray = numpy.loadtxt(CONTROL_DIR_STR + "monomerTrimerProteinArray.txt")
fourmerProteinArray = numpy.loadtxt(CONTROL_DIR_STR + "fourmerProteinArray.txt")

allMonomerLipidArray = numpy.loadtxt(CONTROL_DIR_STR + "allMonomerLipidArray.txt")
monomerDimerLipidArray = numpy.loadtxt(CONTROL_DIR_STR + "monomerDimerLipidArray.txt")
twoDimersLipidArray = numpy.loadtxt(CONTROL_DIR_STR + "twoDimersLipidArray.txt")
monomerTrimerLipidArray = numpy.loadtxt(CONTROL_DIR_STR + "monomerTrimerLipidArray.txt")
fourmerLipidArray = numpy.loadtxt(CONTROL_DIR_STR + "fourmerLipidArray.txt")

allMonomerProteinLipidArray = numpy.loadtxt(CONTROL_DIR_STR + "allMonomerProteinLipidArray.txt")
monomerDimerProteinLipidArray = numpy.loadtxt(CONTROL_DIR_STR + "monomerDimerProteinLipidArray.txt")
twoDimersProteinLipidArray = numpy.loadtxt(CONTROL_DIR_STR + "twoDimersProteinLipidArray.txt")
monomerTrimerProteinLipidArray = numpy.loadtxt(CONTROL_DIR_STR + "monomerTrimerProteinLipidArray.txt")
fourmerProteinLipidArray = numpy.loadtxt(CONTROL_DIR_STR + "fourmerProteinLipidArray.txt")


def gen_final_frames(epocDirStr: str, popSizeInt: int):
    print("Retrieving all frames after the " + str(FRAME_COLLECTION_POINT_INT) + " frame.")
    scoreList = []
    for i in range(popSizeInt):  #0, 1, 2, ... popSizeInt-1
        groFileStr: str = epocDirStr + str(i) + "/npt_prod.gro"
        xtcFileStr: str = epocDirStr + str(i) + "/npt_prod.xtc"
        universe = MDAnalysis.Universe(groFileStr, xtcFileStr)
        proteinAG = universe.select_atoms(PROTEIN_SELECTION_STR)
        lipidAG = universe.select_atoms(LIPID_SELECTION_STR)

        timeSeriesList = []
        for timeStepTS in universe.trajectory:
            currentFrameInt: int = timeStepTS.frame

            if currentFrameInt >= FRAME_COLLECTION_POINT_INT:
                distanceProteinArray = contacts.distance_array(proteinAG.positions, proteinAG.positions, box=timeStepTS.dimensions)
                distanceLipidArray = contacts.distance_array(lipidAG.positions, lipidAG.positions, box=timeStepTS.dimensions)
                distanceProteinLipidArray = contacts.distance_array(proteinAG.positions, lipidAG.positions, box=timeStepTS.dimensions)

                contactsProteinInt: int = contacts.contact_matrix(distanceProteinArray, RADIUS).sum().item()
                contactsLipidInt: int = contacts.contact_matrix(distanceLipidArray, RADIUS).sum().item()
                contactsProteinLipidInt: int = contacts.contact_matrix(distanceProteinLipidArray, RADIUS).sum().item()

                timeSeriesList.append([contactsProteinInt, contactsLipidInt, contactsProteinLipidInt])


        timeSeriesArray = numpy.array(timeSeriesList)
        timeSeriesDF = pandas.DataFrame(timeSeriesArray, columns=[PROTEIN_STR, LIPID_STR, PROTEIN_LIPID_STR])
        scoreList.append(calculate_score(timeSeriesDF, i))

    return scoreList

def calculate_score(timeSeriesDF, individualIDInt: int):
    print("Calculating score for individual " + str(individualIDInt))
    proteinArray = timeSeriesDF[PROTEIN_STR]
    lipidArray = timeSeriesDF[LIPID_STR]
    proteinLipidArray = timeSeriesDF[PROTEIN_LIPID_STR]
    fitnessScoreList = []
   
    proteinList = []
    stat, p = compare_2_groups(proteinArray, allMonomerProteinArray, NORMAL_STR)
    proteinList.append(stat)
    stat, p = compare_2_groups(proteinArray, monomerDimerProteinArray, NORMAL_STR)
    proteinList.append(stat)
    stat, p = compare_2_groups(proteinArray, twoDimersProteinArray, NORMAL_STR)
    proteinList.append(stat)
    fitnessScoreList.append(stat)
    stat, p = compare_2_groups(proteinArray, monomerTrimerProteinArray, NORMAL_STR)
    proteinList.append(stat)
    stat, p = compare_2_groups(proteinArray, fourmerProteinArray, NORMAL_STR)
    proteinList.append(stat)

    lipidList = []
    stat, p = compare_2_groups(lipidArray, allMonomerLipidArray, NORMAL_STR)
    lipidList.append(stat)
    stat, p = compare_2_groups(lipidArray, monomerDimerLipidArray, NORMAL_STR)
    lipidList.append(stat)
    stat, p = compare_2_groups(lipidArray, twoDimersLipidArray, NORMAL_STR)
    lipidList.append(stat)
    fitnessScoreList.append(stat)
    stat, p = compare_2_groups(lipidArray, monomerTrimerLipidArray, NORMAL_STR)
    lipidList.append(stat)
    stat, p = compare_2_groups(lipidArray, fourmerLipidArray, NORMAL_STR)
    lipidList.append(stat)

    proteinLipidList = []
    stat, p = compare_2_groups(proteinLipidArray, allMonomerProteinLipidArray, NORMAL_STR)
    proteinLipidList.append(stat)
    stat, p = compare_2_groups(proteinLipidArray, monomerDimerProteinLipidArray, NORMAL_STR)
    proteinLipidList.append(stat)
    stat, p = compare_2_groups(proteinLipidArray, twoDimersProteinLipidArray, NORMAL_STR)
    proteinLipidList.append(stat)
    fitnessScoreList.append(stat)
    stat, p = compare_2_groups(proteinLipidArray, monomerTrimerProteinLipidArray, NORMAL_STR)
    proteinLipidList.append(stat)
    stat, p = compare_2_groups(proteinLipidArray, fourmerProteinLipidArray, NORMAL_STR)
    proteinLipidList.append(stat)

    return(fitnessScoreList, proteinList, lipidList, proteinLipidList)


def compare_2_groups(array1NP, array2NP, testTypeStr: str):
    if testTypeStr == NORMAL_STR:
        stat, p = ttest_ind(array1NP, array2NP)
    else: 
        stat, p = ranksums(array1NP, array2NP)
        
    return stat, p


def retrieve_score(scoreList):
    returnScoreList = []
    for individualScoresList in scoreList:
        twoDimerScoreList = individualScoresList[0]
        score = sum(twoDimerScoreList)/len(twoDimerScoreList)
        returnScoreList.append(score)
    
    return(returnScoreList)




