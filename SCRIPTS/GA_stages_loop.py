#Contains the 'step' functions of the GA algorithm.
#GA_de_novo.py should call functions found in here.

from gen_structure import *
from record_progress import *
from typing import List
import os
import os.path
from os import path
import sys, getopt
import shutil
import time

MINIMISATION_CG_STR: str = "minimisation_cg"
MINIMISATION_SD_STR: str = "minimisation_sd"

WITH_RESTRAINTS = True

STRONG_STR = "1000"
MEDIUM_STR = "100"
WEAK_STR = "10"

def build_directory_structure_loop(newRunBool: bool, continuationBool: bool, epocDirStr: str, individualInt: int):
    print("0. Preparing file structure.")
    if newRunBool == True and continuationBool == False:

        print("Building new directory structure. This will delete all earlier files.")
        #print("You have 12 seconds to Ctrl-c to exit and preserve any existing data.")
        #time.sleep(12)

        #build the initial epoc directory
        if not path.exists(epocDirStr):
            os.makedirs(epocDirStr)
        #else:
            #it does exists so delete first and then remake
        #    try:
        #        shutil.rmtree(epocDirStr)
        #    except:
        #    	print("Error: Cannot delete directory tree from", epocDirStr, "during a new run.")
        #    	sys.exit()

    	#then each of individual directories (inside the initial epoc)
        #for i in range(popSizeInt):
        dirStr: str = epocDirStr + str(individualInt)
        if not path.exists(dirStr):
            os.makedirs(dirStr)
            create_progress_file(dirStr)
    elif continuationBool == True and newRunBool == False:

        print("Continuation is in progress. Directory structure is retained.")
        #check to see whether the directoryt structure is present
        if not path.exists(epocDirStr):
            print("Error: Cannot perform a continuation run with a missing directory structure.")
            sys.exit()
        else:
            missingList = []
            #for i in range(popSizeInt):
            dirStr: str = epocDirStr + str(individualInt)
            if not path.exists(dirStr):
                missingList.append(dirStr)
            if len(missingList) > 0:
                print("Trying to perform a continuation run. The following population directory were missing:")
                print(missingList)
                sys.exit()
    else:
        print("Error: Settings for new or continuation are incorrect.")
        sys.exit()

#==========================================================================================
def build_population_loop(epocDirStr: str, individualInt: int, newRunBool: bool, supressOutputBool: bool, populationList):
    availableAAList: List = ["R", "H", "D", "K", "E", "S", "T", "N", "Q", "C", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W"]
    scafSequenceStr: str = "LLLLLLGLLLGLLLLLL"
    numChangesToMakeInt: int = 8
    terminalNStr: str = "KRSN"
    terminalCStr: str = "NSRK"
    #populationList: List = []
    n = 0

    #loops 0 to the popSizeInt-1
    print("1. Building initial population.")
    if newRunBool == True:
        print("--New run. \n--Population will be created.")
        while n < 1: #popSizeInt:
            imutableAAPositionList: List = [6,10]
            newDesignStr: str = design_sequence(availableAAList, scafSequenceStr, imutableAAPositionList, numChangesToMakeInt)
            freeEnergyInsertFlt = test_hydrophobicity(newDesignStr, supressOutputBool)
            if freeEnergyInsertFlt < 0.0:
                if newDesignStr not in populationList:
                    newDesignStr = terminalNStr + newDesignStr + terminalCStr
                    #populationList.append(newDesignStr)
                    n = n + 1
                    #print(newDesignStr, " ", freeEnergyInsertFlt)
        #save_population_list(epocDirStr, populationList)
    else:
        print("--Continuation run.")
        if check_population_list(epocDirStr) == True:
            print("--Found the population list for this epoc.")
            populationList = load_population_list(epocDirStr)
            print("--Loaded the population list.")
        else:
            print("Error: Unable to load a population list during a continuation.\nThis is a big problem and requires further dev support.")
            sys.exit()
    #print("Population size: ", len(populationList))    
    return(newDesignStr) #populationList

#======================================================================================
def construct_helix_loop(individualStr: str, individualInt: int, epocDirStr: str, newRunBool: bool):
    print("2. Constructing helix files.")
    if newRunBool == True:
        print("--New helices constructed.")
        #for i in range(len(populationList)):
        #sequenceStr = populationList[i]
        fileLocationStr: str = epocDirStr + str(individualInt) + "/" #str(i) + "/"
        #build_helix(sequenceStr, fileLocationStr)
        build_helix(individualStr, fileLocationStr)
        set_construct_helices(fileLocationStr)
    else:
        print("--Continuation run.")
        print("--Using existing helices.")
        #for i in range(len(populationList)):
        fileLocationStr: str = epocDirStr + str(individualInt) + "/" #str(i) + "/"
        if check_construct_helices(fileLocationStr) == False:
            print("--Constructed helices couldn't be retrieved.\n--Building helices using saved population list.")
            #fileLocationStr: str = epocDirStr + str(i) + "/"
            #sequenceStr = populationList[i]
            build_helix(individualStr, fileLocationStr)
            set_construct_helices(fileLocationStr)
        else:
            print("--Constructed helices exist.")

#======================================================================================
def convert_helices_to_CG_loop(individualStr: str, individualInt: int, epocDirStr: str, newRunBool: bool, supressOutputBool: bool):
    print("3. Convert helices into CG and construct CG model.")
    if newRunBool == True:
        print("--New helices converted into CG model.")
        #for i in range(len(populationList)):
        inputFileStr: str = epocDirStr + str(individualInt) + "/helix_rotated.pdb"
        if path.exists(inputFileStr):
            topFileStr: str = epocDirStr + str(individualInt) + "/helix_cg.top"
            outputPDBFileStr: str = epocDirStr + str(individualInt) + "/helix_cg.pdb"
            itpFileStr: str = epocDirStr + str(individualInt) + "/Protein_A.itp"
            build_martinize(inputFileStr, topFileStr, outputPDBFileStr, itpFileStr, supressOutputBool)
            progressDirStr: str = epocDirStr + str(individualInt)
            set_convert_helices(progressDirStr)
        else:
            raise FileNotFoundError("Error: File not found whilst building CG structures.")
        print("Finished building CG helix topology (top) and coordinate (pdb) files.")
    else:
        print("--Continuation run.")
        print("--Using existing converted helices.")
        #for i in range(len(populationList)):
        progressDirStr: str = epocDirStr + str(individualInt)
        if check_convert_helices(progressDirStr) == False:
            print("--Converted CG helices couldn't be retrieved.\nConverted helices using saved data.")
            inputFileStr: str = epocDirStr + str(individualInt) + "/helix_rotated.pdb"
            if path.exists(inputFileStr):
                topFileStr: str = epocDirStr + str(individualInt) + "/helix_cg.top"
                outputPDBFileStr: str = epocDirStr + str(individualInt) + "/helix_cg.pdb"
                itpFileStr: str = epocDirStr + str(individualInt) + "/Protein_A.itp"
                build_martinize(inputFileStr, topFileStr, outputPDBFileStr, itpFileStr, supressOutputBool)
                progressDirStr: str = epocDirStr + str(individualInt)
                set_convert_helices(progressDirStr)
            else:
                 raise FileNotFoundError("Error: File not found whilst building CG structures.")

#======================================================================================
def prepare_helix(individualInt: int, epocDirStr: str, parentDirStr: str, numMPI, numOMP, supressOutputBool: bool):
    topFileStr: str = epocDirStr + str(individualInt) + "/helix_temp_cg.top"
    inputPDBFileStr: str = epocDirStr + str(individualInt) + "/helix_cg.pdb"
    inputGROFileStr: str = epocDirStr + str(individualInt) + "/helix_cg_temp.gro"

    #I can't use the helix_cg.pdb as it has no unit cell dimensions
    #convert helix_cg.pdb to helix_cg.gro 
    convert_pdb_to_gro(inputPDBFileStr, inputGROFileStr, supressOutputBool)

    #edit the last line to include dimensions (but be done by hand to avoid messing with the coordinates)
    lines = open(inputGROFileStr, 'r').readlines()
    newLastLine = "  10.00000  10.00000  10.00000\n"
    lines[-1] = newLastLine
    open(inputGROFileStr, 'w').writelines(lines)
   
    #run editconf to put the centre of the protein in the centre of the unit cell otherwise everything breaks from here
    copyCommandList = ["gmx", "editconf", "-f", inputGROFileStr, "-o", inputGROFileStr, "-c"]
    result = subprocess.run(copyCommandList, capture_output=True, text=True)
    if supressOutputBool == False:
        print("stderr:", result.stderr)
        print("stdout:", result.stdout)

    #with open(inputGROFileStr, 'a') as fileObj:
    #    fileObj.writelines(newLastLine)
    #fileObj.close()


    mdpSDFileStr: str = parentDirStr + "helix_SD.mdp"
    #mdpCGFileStr: str = parentDirStr + "helix_CG.mdp"

    outputSDFileStr: str = epocDirStr + str(individualInt) + "/helix_SD.tpr"
    #outputCGFileStr: str = epocDirStr + str(individualInt) + "/helix_CG.tpr"

    itpFileStr: str = epocDirStr + str(individualInt) + "/Protein_A.itp"

    #build a temp top file
    write_helix_top_file(topFileStr, itpFileStr, "temp_helix")

    #perform minisation steps using the converted GRO file
    run_grompp(mdpSDFileStr, inputGROFileStr, topFileStr, outputSDFileStr, supressOutputBool, False, epocDirStr, None, False)
    #check how succesful
    deffnmFileStr: str = epocDirStr + str(individualInt) + "/helix_minimisation_SD"
    run_simulation(deffnmFileStr, outputSDFileStr, False, numMPI, numOMP, supressOutputBool)
    #exit()
    successBool = check_helix_minimisation(epocDirStr + str(individualInt), MINIMISATION_SD_STR)
    if successBool == False:
        print("ERROR: Steepest descent minimisation.")
        return(False)

    #perform minisation steps (CG)
    #run_grompp(mdpCGFileStr, epocDirStr + str(individualInt) + "/helix_minimisation_SD.gro", topFileStr, outputCGFileStr, supressOutputBool, False, epocDirStr, None)
    #check how succesful
    #deffnmFileStr: str = epocDirStr + str(individualInt) + "/helix_minimisation_CG"
    #run_simulation(deffnmFileStr, outputCGFileStr, False, numMPI, numOMP, supressOutputBool)
    #successBool = check_helix_minimisation(epocDirStr + str(individualInt), MINIMISATION_CG_STR)
    #if successBool == False:
    #    print("ERROR: Conjugate gradient minimisation.")
    #    return(False)

    #the energy minimisation may have moved the coordinates too much for a viable membrane insertion, therefore move them back into place.
    #helixTranslationList = calculate_helix_distance(epocDirStr + str(individualInt) + "/helix_minimisation_SD.gro", epocDirStr + str(individualInt) + "/helix_cg.pdb")

    #I must convert the output from the CG minimisation back to helix_cg.pdb
    helixTranslationList = None
    convert_gro_to_pdb(epocDirStr + str(individualInt) + "/helix_minimisation_SD.gro", epocDirStr + str(individualInt) + "/helix_cg.pdb", supressOutputBool, helixTranslationList)

    return(True)

#======================================================================================
def adjust_helix_position_loop(individualInt: int, newRunBool: bool, epocDirStr: str, supressOutputBool: bool):
    print("4. Copy the helix files (three times) and move them into place.")
    numHelicesInt: int = 4

    xCoordList = [2.263, 2.166, 6.980, 7.040]
    yCoordList = [2.252, 6.994, 1.876, 6.672]
    zCoordList = [6.623, 6.427, 6.774, 6.868]

    #for i in range(len(populationList)):
    if newRunBool == True:
        print("--Copying and moving helicex for the first time (new run).")
        for n in range(numHelicesInt):
            copyCommandList = ["cp", epocDirStr + str(individualInt) + "/helix_cg.pdb", epocDirStr + str(individualInt) + "/helix_translated_" + str(n) + "_cg.pdb"]
            result = subprocess.run(copyCommandList, capture_output=True, text=True)
            if supressOutputBool == False:
                print("stderr:", result.stderr)
                print("stdout:", result.stdout)

            inputFileStr: str = epocDirStr + str(individualInt) + "/helix_cg.pdb"
            #the next line translate the helix and converts the output from pdb to gro
            outputFileStr: str = epocDirStr + str(individualInt) + "/helix_translated_" + str(n) + "_cg.gro"
            moveHelix(inputFileStr, outputFileStr, xCoordList[n], yCoordList[n], zCoordList[n], supressOutputBool)
        set_copy_move_helices(epocDirStr + str(individualInt))
    else:
        print("--Continuation run.")
        print("--Copying and moving helices using existing data.")
        if check_copy_move_helices(epocDirStr + str(individualInt)) == False:
            print("--Unable to retrieve helices to move. Attempting again.")
            for n in range(numHelicesInt):
                copyCommandList = ["cp", epocDirStr + str(individualInt) + "/helix_cg.pdb", epocDirStr + str(individualInt) + "/helix_translated_" + str(n) + "_cg.pdb"]
                result = subprocess.run(copyCommandList, capture_output=True, text=True)
                if supressOutputBool == False:
                    print("stderr:", result.stderr)
                    print("stdout:", result.stdout)

                inputFileStr: str = epocDirStr + str(individualInt) + "/helix_cg.pdb"
                #the next line translate the helix and converts the output from pdb to gro
                outputFileStr: str = epocDirStr + str(individualInt) + "/helix_translated_" + str(n) + "_cg.gro"
                moveHelix(inputFileStr, outputFileStr, xCoordList[n], yCoordList[n], zCoordList[n], supressOutputBool)
            set_copy_move_helices(epocDirStr + str(individualInt))

#======================================================================================
def write_top_gro_files_loop(individualInt: int, epocDirStr: str, scriptsDirStr, newRunBool: bool, sequenceStr: str, supressOutputBool: bool):
    print("5. Write the TOP file, restraint file and construct the combined GRO file.")
    if newRunBool == True:
        #for i in range(popSizeInt):
            #build the restraint file for the helices

        #build the restraint file for the helices
        #STRONG, MEDIUM, WEAK
        outputRestFileStr: str = epocDirStr + str(individualInt) + "/STRONG_posre.itp" #"/Protein_restr_A.itp"
        groFileStr: str = epocDirStr + str(individualInt) + "/helix_cg.pdb"
        systemInpFileStr: str = scriptsDirStr + "System.inp"
        restrCommandList = ["gmx", "genrestr", "-f", groFileStr, "-o", outputRestFileStr, "-fc", STRONG_STR, STRONG_STR, STRONG_STR]
        with open(systemInpFileStr, "r") as systemInpFileObj:
            resultRestr = subprocess.run(restrCommandList, stdin=systemInpFileObj)
            if supressOutputBool == False:
                print("stderr:", resultRestr.stderr)
                print("stdout:", resultRestr.stdout)
        systemInpFileObj.close()

        outputRestFileStr: str = epocDirStr + str(individualInt) + "/MEDIUM_posre.itp"
        restrCommandList = ["gmx", "genrestr", "-f", groFileStr, "-o", outputRestFileStr, "-fc", MEDIUM_STR, MEDIUM_STR, MEDIUM_STR]
        with open(systemInpFileStr, "r") as systemInpFileObj:
            resultRestr = subprocess.run(restrCommandList, stdin=systemInpFileObj)
            if supressOutputBool == False:
                print("stderr:", resultRestr.stderr)
                print("stdout:", resultRestr.stdout)
        systemInpFileObj.close()

        outputRestFileStr: str = epocDirStr + str(individualInt) + "/WEAK_posre.itp"
        restrCommandList = ["gmx", "genrestr", "-f", groFileStr, "-o", outputRestFileStr, "-fc", WEAK_STR, WEAK_STR, WEAK_STR]
        with open(systemInpFileStr, "r") as systemInpFileObj:
            resultRestr = subprocess.run(restrCommandList, stdin=systemInpFileObj)
            if supressOutputBool == False:
                print("stderr:", resultRestr.stderr)
                print("stdout:", resultRestr.stdout)
        systemInpFileObj.close()

        #get peptide charge (which will be *4)
        chargeInt: int = calculate_peptide_charge(sequenceStr) * 4

        #construct the system level topology file
        topFileStr = epocDirStr + str(individualInt) + "/system.top"
        itpFileStr = epocDirStr + str(individualInt) + "/Protein_A.itp"
        #sequenceStr = populationList[i]
        restraintsDirStr = epocDirStr + str(individualInt) # + "/posre.itp" #"/Protein_restr_A.itp"
        #the returning naCountInt and clCountInt are the number of ions to keep in the final gro file
        naCountInt, clCountInt = write_top_file(topFileStr, itpFileStr, sequenceStr, restraintsDirStr, chargeInt)

        #write the combined gro file
        bilayerGroFileStr: str = scriptsDirStr + "temp_template.gro"
        helix0GroFileStr: str = epocDirStr + str(individualInt) + "/helix_translated_0_cg.gro"
        helix1GroFileStr: str = epocDirStr + str(individualInt) + "/helix_translated_1_cg.gro"
        helix2GroFileStr: str = epocDirStr + str(individualInt) + "/helix_translated_2_cg.gro"
        helix3GroFileStr: str = epocDirStr + str(individualInt) + "/helix_translated_3_cg.gro"
        inputFileList = [bilayerGroFileStr, helix0GroFileStr, helix1GroFileStr, helix2GroFileStr, helix3GroFileStr]
        outputFileStr: str = epocDirStr + str(individualInt) + "/combined.gro"
        write_gro_file(inputFileList, outputFileStr)

        #adjust the charge on the water (these are the number of ions we need)
        adjustCharge(outputFileStr, outputFileStr, naCountInt, clCountInt)

        #reorder the atom number (index) in the gro files
        commandList = ["gmx", "editconf", "-f", outputFileStr, "-o", outputFileStr]
        reorderResult = subprocess.run(commandList, capture_output=True, text=True)
        if supressOutputBool == False:
            print("stdout:", reorderResult.stdout)
            print("stderr:", reorderResult.stderr)

        #build the index file required for all NVT/NPT simulations
        indexFileStr: str = epocDirStr + str(individualInt) + "/index.ndx"
        indexInpFileStr: str = scriptsDirStr + "indexInput.inp"
        run_make_index(outputFileStr, indexFileStr, indexInpFileStr, supressOutputBool)
        set_top_geometry(epocDirStr + str(individualInt)) 
    else:
        print("--Continuation run.")
        #for i in range(popSizeInt):
        if check_top_geometry(epocDirStr + str(individualInt)) == False:
            print("--Reconstructing topology and geometry files.")
            
            #build the restraint file for the helices
            #STRONG, MEDIUM, WEAL
            outputRestFileStr: str = epocDirStr + str(individualInt) + "/STRONG_posre.itp" #"/Protein_restr_A.itp"
            groFileStr: str = epocDirStr + str(individualInt) + "/helix_cg.pdb"
            systemInpFileStr: str = scriptsDirStr + "System.inp"
            restrCommandList = ["gmx", "genrestr", "-f", groFileStr, "-o", outputRestFileStr, "-fc", STONG_STR, STRONG_STR, STRONG_STR]
            with open(systemInpFileStr, "r") as systemInpFileObj:
                resultRestr = subprocess.run(restrCommandList, stdin=systemInpFileObj)
                if supressOutputBool == False:
                    print("stderr:", resultRestr.stderr)
                    print("stdout:", resultRestr.stdout)
            systemInpFileObj.close()

            outputRestFileStr: str = epocDirStr + str(individualInt) + "/MEDIUM_posre.itp"
            restrCommandList = ["gmx", "genrestr", "-f", groFileStr, "-o", outputRestFileStr, "-fc", MEDIUM_STR, MEDIUM_STR, MEDIUM_STR]
            with open(systemInpFileStr, "r") as systemInpFileObj:
                resultRestr = subprocess.run(restrCommandList, stdin=systemInpFileObj)
                if supressOutputBool == False:
                    print("stderr:", resultRestr.stderr)
                    print("stdout:", resultRestr.stdout)
            systemInpFileObj.close()

            outputRestFileStr: str = epocDirStr + str(individualInt) + "/WEAK_posre.itp"
            restrCommandList = ["gmx", "genrestr", "-f", groFileStr, "-o", outputRestFileStr, "-fc", WEAK_STR, WEAK_STR, WEAK_STR]
            with open(systemInpFileStr, "r") as systemInpFileObj:
                resultRestr = subprocess.run(restrCommandList, stdin=systemInpFileObj)
                if supressOutputBool == False:
                    print("stderr:", resultRestr.stderr)
                    print("stdout:", resultRestr.stdout)
            systemInpFileObj.close()



            #get peptide charge (which will be *4)
            chargeInt: int = calculate_peptide_charge(populationList[i]) * 4

            #construct the system level topology file
            topFileStr = epocDirStr + str(individualInt) + "/system.top"
            itpFileStr = epocDirStr + str(individualInt) + "/Protein_A.itp"
            #sequenceStr = populationList[i]
            restraintsDirStr = epocDirStr + str(individualInt) #+ "/Posre.itp" #"/Protein_restr_A.itp"
            #the returning naCountInt and clCountInt are the number of ions to keep in the final gro file
            naCountInt, clCountInt = write_top_file(topFileStr, itpFileStr, sequenceStr, restraintsDirStr, chargeInt)

            #write the combined gro file
            bilayerGroFileStr: str = scriptsDirStr + "temp_template.gro"
            helix0GroFileStr: str = epocDirStr + str(individualInt) + "/helix_translated_0_cg.gro"
            helix1GroFileStr: str = epocDirStr + str(individualInt) + "/helix_translated_1_cg.gro"
            helix2GroFileStr: str = epocDirStr + str(individualInt) + "/helix_translated_2_cg.gro"
            helix3GroFileStr: str = epocDirStr + str(individualInt) + "/helix_translated_3_cg.gro"
            inputFileList = [bilayerGroFileStr, helix0GroFileStr, helix1GroFileStr, helix2GroFileStr, helix3GroFileStr]
            outputFileStr: str = epocDirStr + str(individualInt) + "/combined.gro"
            write_gro_file(inputFileList, outputFileStr)

            #adjust the charge on the water (these are the number of ions we need)
            adjustCharge(outputFileStr, outputFileStr, naCountInt, clCountInt)

            #reorder the atom number (index) in the gro files
            commandList = ["gmx", "editconf", "-f", outputFileStr, "-o", outputFileStr]
            reorderResult = subprocess.run(commandList, capture_output=True, text=True)
            if supressOutputBool == False:
                print("stdout:", reorderResult.stdout)
                print("stderr:", reorderResult.stderr)

            #build the index file required for all NVT/NPT simulations
            indexFileStr: str = epocDirStr + str(individualInt) + "/index.ndx"
            indexInpFileStr: str = scriptsDirStr + "indexInput.inp"
            run_make_index(outputFileStr, indexFileStr, indexInpFileStr, supressOutputBool)
            set_top_geometry(epocDirStr + str(individualInt))

#======================================================================================
def run_energy_minimisation_loop(parentDirStr, epocDirStr, individualInt, newRunBool, supressOutputBool, numMPI, numOMP, minimisationTypeStr):
#This won't be running conjugate gradient anymore. It is to unreliable. 

    print("6. Performing energy minimisation using " + minimisationTypeStr + " over epoc.")

    if minimisationTypeStr == STEEPEST_DESCENT:
        mdpFileStr: str = parentDirStr + "minimization_SD.mdp"
    else:
        mdpFileStr: str = parentDirStr + "minimization_CG.mdp"

    #for i in range(popSizeInt):
    if minimisationTypeStr == STEEPEST_DESCENT:
        groFileStr: str = epocDirStr + str(individualInt) + "/combined.gro"
    else:
        groFileStr: str = epocDirStr + str(individualInt) + "/minimization.gro"
    #restrFileStr: str = epocDirStr + str(individualInt) + "/helix_translated_0_cg.gro"
    topFileStr: str = epocDirStr + str(individualInt) + "/system.top"
    outputFileStr: str = epocDirStr + str(individualInt) + "/minimization.tpr"
    #a new run
    if newRunBool == True:
        run_grompp(mdpFileStr, groFileStr, topFileStr, outputFileStr, supressOutputBool, EQ_VELOCITIES_FALSE, epocDirStr + str(individualInt), None, False)
        #set the grompp
        #set_grompp(epocDirStr + str(individualInt), "minimization.tpr", MINIMISATION_GROMPP_SUCCESS, MINIMISATION_GROMPP_FAIL, RERUN_FALSE)
    else:
    #a continuation run
        print("--Continuation run. Checking whether grompp needs to rerun.")
        #If the check failed
        if check_minimisation_grompp(epocDirStr + str(individualInt)) == False:
            print("--Unsuccessful or missing tpr file. Building again.")
            run_grompp(mdpFileStr, groFileStr, topFileStr, outputFileStr, supressOutputBool, EQ_VELOCITIES_FALSE, epocDirStr + str(individualInt), None, False)
            #set the grompp
            #set_grompp(epocDirStr + str(individualInt), "minimization.tpr", MINIMISATION_GROMPP_SUCCESS, MINIMISATION_GROMPP_FAIL, RERUN_TRUE)
        else:
            print("--Found existing md tpr file.")

    deffnmFileStr: str = epocDirStr + str(individualInt) + "/minimization"
    if newRunBool == True:
        run_simulation(deffnmFileStr, outputFileStr, False, numMPI, numOMP, supressOutputBool)
        print("I am here... check the simulation progress if it is energy minimisation...")

        #set the MD
        #set_mdrun(epocDirStr + str(individualInt), "minimization.log", MINIMISATION_MDRUN_SUCCESS, MINIMISATION_MDRUN_FAIL, MINIMISATION_CONVERGE, RERUN_FALSE)
    else:
        print("--Continuation run. Checking whether mdrun is required.")
        #If the check failed
        if check_minimisation_mdrun(epocDirStr + str(individualInt)) == False:
            print("--Unsuccessful or missing energy minimisation. Performing mdrun.")
            run_simulation(deffnmFileStr, outputFileStr, False, numMPI, numOMP, supressOutputBool)
            #set the MD
            #set_mdrun(epocDirStr + str(individualInt), "minimization.log", MINIMISATION_MDRUN_SUCCESS, MINIMISATION_MDRUN_FAIL, MINIMISATION_CONVERGE, RERUN_TRUE)
        else:
            print("--Found existing md trajectory files.")

#======================================================================================
def run_equilibrium_npt_loop(parentDirStr: str, epocDirStr: str, supressOutputBool: bool, newRunBool: bool, numMPI: int, numOMP: int, individualInt: int, equilibrationTypeStr: str):
    print("7. Performing equilibrium NPT simulations.")
    #if equilibrationTypeStr == SHORT_EQ_REST:
    #   mdpFileStr: str = parentDirStr + "NPT_short_eq_rest.mdp"
    #elif equilibrationTypeStr == LONG_EQ_REST:
    #    mdpFileStr: str = parentDirStr + "NPT_long_eq_rest.mdp"
    if equilibrationTypeStr == SHORT_EQ_STRONG:
        mdpFileStr: str = parentDirStr + "NPT_short_eq_STRONG.mdp"
    elif equilibrationTypeStr == LONG_EQ_STRONG:
        mdpFileStr: str = parentDirStr + "NPT_long_eq_STRONG.mdp"
    elif equilibrationTypeStr == LONG_EQ_MEDIUM:
        mdpFileStr: str = parentDirStr + "NPT_long_eq_MEDIUM.mdp"
    elif equilibrationTypeStr == LONG_EQ_WEAK:
        mdpFileStr: str = parentDirStr + "NPT_long_eq_WEAK.mdp"
    else: #where LONG_EQ
        mdpFileStr: str = parentDirStr + "NPT_long_eq.mdp"

    #for i in range(popSizeInt):

    if equilibrationTypeStr == SHORT_EQ_STRONG:
        groFileStr: str = epocDirStr + str(individualInt) + "/minimization.gro"
        velocitiesStr: str = EQ_VELOCITIES_FALSE
        velocityFileStr = None
    else:
        groFileStr: str = epocDirStr + str(individualInt) + "/npt_eq.gro"
        velocitiesStr: str = EQ_VELOCITIES_TRUE
        velocityFileStr: str = epocDirStr + str(individualInt) + "/npt_eq.cpt"
    topFileStr: str = epocDirStr + str(individualInt) + "/system.top"
    outputFileStr: str = epocDirStr + str(individualInt) + "/npt_eq.tpr"
    if newRunBool == True:
        run_grompp(mdpFileStr, groFileStr, topFileStr, outputFileStr, supressOutputBool, velocitiesStr, epocDirStr + str(individualInt), velocityFileStr, True)
        #set the grompp
        #set_grompp(epocDirStr + str(individualInt), "npt_eq.tpr", NPT_EQ_GROMPP_SUCCESS, NPT_EQ_GROMPP_FAIL, RERUN_FALSE)
    else:
        print("--Continuation run. Checking whether grompp needs to rerun.")
        if check_x_grompp(epocDirStr + str(individualInt)) == False:
            print("--Unsuccessful or missing tpr file. Building again.")
            run_grompp(mdpFileStr, groFileStr, topFileStr, outputFileStr, supressOutputBool, velocitiesStr, epocDirStr + str(individualInt), velocityFileStr, True)
            #set grompp
            #set_grompp(epocDirStr + str(individualInt), "npt_eq.tpr", NPT_EQ_GROMPP_SUCCESS, NPT_EQ_GROMPP_FAIL, RERUN_FALSE)
        else:
            print("--Found existing md tpr file.")

    deffnmFileStr: str = epocDirStr + str(individualInt) + "/npt_eq"
    if newRunBool == True:
        run_simulation(deffnmFileStr, outputFileStr, False, numMPI, numOMP, supressOutputBool)
        #set the MD
        #set_mdrun(epocDirStr + str(individualInt), "npt_eq.log", NPT_EQ_MDRUN_SUCCESS, NPT_EQ_MDRUN_FAIL, NPT_MDRUN_FINISHED, RERUN_FALSE)
    else:
        print("--Continuation run. Checking whether mdrun is required.")
        #If the check failed
        if check_npt_eq_mdrun(epocDirStr + str(individualInt)) == False:
            print("Unsuccessful or missing NPT equilibration. Performing mdrun.")
            run_simulation(deffnmFileStr, outputFileStr, False, numMPI, numOMP, supressOutputBool)
            #set the MD
            #set_mdrun(epocDirStr + str(individualInt), "npt_eq.log", NPT_EQ_MDRUN_SUCCESS, NPT_EQ_MDRUN_FAIL, NPT_MDRUN_FINISHED, RERUN_TRUE)
        else:
            print("--Found existing md trajectory files.")

    #after the long equilibrium check to see whether the simulation ran. If it didn't, then continue creating new simulations. 
    return check_npt_eq_mdrun(epocDirStr + str(individualInt))

#==============================================================================================
def run_production_npt_loop(parentDirStr: str, epocDirStr, supressOutputBool: bool, newRunBool: bool, numMPI: int, numOMP, individualInt: int):
    print("8. Performing production NPT simulations.")
    #for i in range(popSizeInt):
    mdpFileStr: str = parentDirStr + "NPT_prod.mdp"
    groFileStr: str = epocDirStr + str(individualInt) + "/npt_eq.gro"
    topFileStr: str = epocDirStr + str(individualInt) + "/system.top"
    outputFileStr: str = epocDirStr + str(individualInt) + "/npt_prod.tpr"
    cptFileStr: str = epocDirStr + str(individualInt) + "/npt_eq.cpt"
    if newRunBool == True:
        run_grompp(mdpFileStr, groFileStr, topFileStr, outputFileStr, supressOutputBool, EQ_VELOCITIES_TRUE, epocDirStr + str(individualInt), cptFileStr, False)
        #set the grompp
        #set_grompp(epocDirStr + str(individualInt), "npt_prod.tpr", NPT_PRODUCTION_GROMPP_SUCCESS, NPT_PRODUCTION_GROMPP_FAIL, RERUN_FALSE)
    else:
        print("--Continuation run. Checking whether grompp needs to rerun.")
        if check_production_grompp(epocDirStr + str(individualInt)) == False:
            print("--Unsuccessful or missing tpr file. Building again.")
            run_grompp(mdpFileStr, groFileStr, topFileStr, outputFileStr, supressOutputBool, EQ_VELOCITIES_TRUE, epocDirStr + str(individualInt), cptFileStr, False)
            #set grompp
            #set_grompp(epocDirStr + str(individualInt), "npt_prod.tpr", NPT_PRODUCTION_GROMPP_SUCCESS, NPT_PRODUCTION_GROMPP_FAIL, RERUN_FALSE)
        else:
            print("--Found existing md tpr file.")

    deffnmFileStr: str = epocDirStr + str(individualInt) + "/npt_prod"
    if newRunBool == True:
        run_simulation(deffnmFileStr, outputFileStr, False, numMPI, numOMP, supressOutputBool)
        #set the MD
        #set_mdrun(epocDirStr + str(individualInt), "npt_prod.log", NPT_PRODUCTION_SUCCESS, NPT_PRODUCTION_FAIL, NPT_MDRUN_FINISHED, RERUN_FALSE)
    else:
        print("--Continuation run. Checking whether mdrun is required.")
        #If the check failed
        if check_production_mdrun(epocDirStr + str(individualInt)) == False:
            print("Unsuccessful or missing NPT production. Performing mdrun.")
            run_simulation(deffnmFileStr, outputFileStr, False, numMPI, numOMP, supressOutputBool)
            #set the MD
            #set_mdrun(epocDirStr + str(individualInt), "npt_prod.log", NPT_PRODUCTION_SUCCESS, NPT_PRODUCTION_FAIL, NPT_MDRUN_FINISHED, RERUN_TRUE)
        else:
            print("--Found existing md trajectory files.")

    #after the long equilibrium check to see whether the simulation ran. If it didn't, then continue creating new simulations.
    return check_npt_production_mdrun(epocDirStr + str(individualInt))

#==============================================================================================
def sp_crossover(individualOneStr: str, individualTwoStr: str):
    pointInt: int = random.randint(4,19)

    individualOneLeftStr: str = individualOneStr[0:pointInt]
    individualOneRightStr: str = individualOneStr[pointInt:len(individualOneStr)+1]

    individualTwoLeftStr: str = individualTwoStr[0:pointInt]
    individualTwoRightStr: str = individualTwoStr[pointInt:len(individualTwoStr)+1]

    newIndividualOneStr = individualOneLeftStr + individualTwoRightStr
    newIndividualTwoStr = individualTwoLeftStr + individualOneRightStr

    return newIndividualOneStr, newIndividualTwoStr


#================================================================================================
def pair_individuals_for_crossover(populationList: List):
    
    #shuffle the population list to make this look random. 
    random.shuffle(populationList)

    middleIndex =  len(populationList)//2

    firstList: List = populationList[:middleIndex]
    secondList: List = populationList[:middleIndex:]


    #firstList: List = populationList[0:((len(populationList)-1)/2)]
    #secondList: List = populationList[len(populationList)/2:len(populationList)-1]

    if len(firstList) != len(secondList):
        print("Error: selected individuals for crossover could not be divided into two equal lists.")
        sys.exit()

    newPopulationList: List = []
    for i in range(len(firstList)):
        newFirstIndividualStr, newSecondIndividualStr = sp_crossover(firstList[i], secondList[i])
        newPopulationList.append(newFirstIndividualStr)
        newPopulationList.append(newSecondIndividualStr)

    return newPopulationList


#================================================================================================
def check_for_mutation(populationList: List, mutateChanceFlt: float) -> List:

    newPopulationList: List = []
    for individualStr in populationList:
        chanceFlt = random.random()
        if chanceFlt <= mutateChanceFlt:
            newIndividualStr = mutate_individual(individualStr)
            newPopulationList.append(newIndividualStr)
        else:
            newPopulationList.append(individualStr)

    return newPopulationList


#================================================================================================
def mutate_individual(individualStr: str) -> str:
    residueNumbersList = [4,5,6,7,8,9,11,12,13,15,16,17,18,19,20]
    availableAAList: List = ["R", "H", "D", "K", "E", "S", "T", "N", "Q", "C", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W"]

    pointInt: int = random.randint(0, len(residueNumbersList)-1)
    pointInt = residueNumbersList[pointInt]

    isDifferentBool: bool = False
    while isDifferentBool == False:
        toMutateStr: str = individualStr[pointInt]
        aaPointInt: int = random.randint(0, len(availableAAList)-1)
        aaNewStr: str = availableAAList[aaPointInt]
        if aaNewStr != toMutateStr:
            isDifferentBool = True

    individualNewStr = individualStr[0:pointInt] + aaNewStr + individualStr[pointInt+1:len(individualStr)]
    
    return individualNewStr

#==================================================================================================
def select_stochastic_population(populationList: List, populationScoreList: List) -> List: #return List
    totalScoreFlt = sum(populationScoreList)

    #all scores fall between 0 and 1. Total score = 1. 
    adjustedScoreList: List = []
    for score in populationScoreList:
        adjustedScoreList.append(score/totalScoreFlt)



    #line the scores up and adjust them so they form a score line between 0 and 1.
    usScoreList: List = []
    for i in range(len(adjustedScoreList)):
        if i == 0:
            usScoreList.append(adjustedScoreList[i])
        else:
            usScoreList.append(adjustedScoreList[i] + usScoreList[i-1])


    #randomly select the scores
    selectedIndividualList: List = []
    for i in range(len(populationList)):
        randSelectDbl = random.random()
        print(randSelectDbl)
        for j in range(len(usScoreList)):
            if j == 0:
                if randSelectDbl >=0 and randSelectDbl < usScoreList[j]:
                    selectedIndividualList.append(j)
                    break
            else:
                if randSelectDbl >= usScoreList[j-1] and randSelectDbl < usScoreList[j]:
                    selectedIndividualList.append(j)
                    break

    #check that the right number of selections had been made
    if len(selectedIndividualList) != len(populationList):
        print("Error: The stochastic universal sampling failed to select the right number of individuals.")
        sys.exit()
    
    newPopulationList: List = [populationList[i] for i in selectedIndividualList]

    return newPopulationList
