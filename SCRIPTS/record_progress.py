#all code related to the management of restarting and saving progress
#status to consider

#created individual file data
#minimisation generated
#minimisation finished

#npt eq generated
#npt eq finished

#npt production generated
#npt production finished

#I need to ask the question:
#did it grompp successfully?
#did it mdrun successfully?

import os.path

RERUN_TRUE = True
RERUN_FALSE = False

POPULATION_LIST_SUCCESS: str = "Population list creation successful"

CONSTRUCT_HELICES_SUCCESS: str = "Helices construction successful"
CONVERT_HELICES_SUCCESS: str = "Helices conversion to CG successful"
COPY_MOVE_HELICES_SUCCESS: str = "Copying and moving helices successful"
TOP_GEOMETRY_SUCCESS: str = "Generating top and geometry successful"

MINIMISATION_GROMPP_SUCCESS: str = "Minimisation grompp successful"
MINIMISATION_GROMPP_FAIL: str = "Minimisation grompp failed"

MINIMISATION_MDRUN_SUCCESS: str = "Minimisation mdrun successful"
MINIMISATION_MDRUN_FAIL: str = "Minimisation mdrun failed"
MINIMISATION_CONVERGE: str = "Steepest Descents converged to Fmax < 1000"
MINIMISATION_FLT_ERROR: str = "Floating point exception"

NPT_EQ_GROMPP_SUCCESS: str = "NPT eq grompp successful"
NPT_EQ_GROMPP_FAIL: str = "NPT eq grompp failed"

NPT_EQ_MDRUN_SUCCESS: str = "NPT eq mdrun successful"
NPT_EQ_MDRUN_FAIL: str = "NPT eq mdrun failed"

NPT_PRODUCTION_GROMPP_SUCCESS: str = "NPT production grompp successful"
NPT_PRODUCTION_GROMPP_FAIL: str = "NPT production grompp failed"

NPT_PRODUCTION_SUCCESS: str = "NPT production mdrun successful"
NPT_PRODUCTION_FAIL: str = "NPT production mdrun failed"

NPT_MDRUN_FINISHED: str = "Finished mdrun on"

ANALYSIS_SUCCESS: str = "Analysis successful"
ANALYSIS_FAIL: str = "Analysis failed"

EPOC_PROGRESS_FILE_STR: str = "epoc_progress.txt"

MINIMISATION_CG_STR: str = "minimisation_cg"
MINIMISATION_SD_STR: str = "minimisation_sd"

MAXIMUM_FORCE_STR: str = "Maximum force"
HELIX_MAX_FORCE_FLT = 750.0

#=============================================================================
#=============================================================================

#this is the new and improved way of loading a continuation simulation
#read in the file and adjust the loop logic
def save_epoc_progress(parentDirStr: str, epocInt: int, individualInt: int):
    fileStr: str = parentDirStr + EPOC_PROGRESS_FILE_STR

    with open(fileStr, "w") as fileObj:
        fileObj.write("epoc\tindividual\n")
        fileObj.write(str(epocInt) + "\t" + str(individualInt))
    fileObj.close()


#returns the epoc and the individual (as integers)
def load_epoc_progress(parentDirStr: str):
    fileStr: str = parentDirStr + EPOC_PROGRESS_FILE_STR

    with open(fileStr, "r") as fileObj:
        for lineStr in fileObj:
            lineList = lineStr.split()

    return(int(lineList[0]), int(lineList[1]))



#Checking the population list during a continuation run simply involves looking
#to see whether the population_list is present in each EPOC_no directory. 
def save_population_list(epocDirStr: str, populationList):
    if epocDirStr[-1] == "/":
        fileStr: str = epocDirStr + "population_list.txt"
    else:
        fileStr: str = epocDirStr + "/population_list.txt"

    with open(fileStr, "w") as fileObj:
        for individualStr in populationList:
            fileObj.write(individualStr + "\n")
    fileObj.close()


def check_population_list(epocDirStr: str) -> bool:
    if epocDirStr[-1] == "/":
        fileStr: str = epocDirStr + "population_list.txt"
    else:
        fileStr: str = epocDirStr + "/population_list.txt"

    resultBool: bool = es.path.exists(fileStr)
    return(resultBool)


def load_population_list(epocDirStr: str):
    if epocDirStr[-1] == "/":
        fileStr: str = epocDirStr + "population_list.txt"
    else:
        fileStr: str = epocDirStr + "/population_list.txt"

    populationList = []
    with open(fileStr, "r") as fileObj:
        for line in fileObj:
            populationList.append(line)
    fileObj.close()

    return populationList

def save_population_individual(epocDirStr: str, individualStr: str):
    if epocDirStr[-1] == "/":
        fileStr: str = epocDirStr + "population_list.txt"
    else:
        fileStr: str = epocDirStr + "/population_list.txt"

    with open(fileStr, "a") as fileObj:
        for individualStr in populationList:
            fileObj.write(individualStr + "\n")
    fileObj.close()

#=============================================================================
#=============================================================================

def set_construct_helices(dirStr: str):
    if dirStr[-1] == "/":
        fileStr: str = dirStr + "progress.txt"
    else:
        fileStr: str = dirStr + "/progress.txt" 
    with open(fileStr, "a") as fileObj:
        fileObj.write(CONSTRUCT_HELICES_SUCCESS + "\n")
    fileObj.close()


def check_construct_helices(dirStr: str):
    if dirStr[-1] == "/":
        fileStr: str = dirStr + "progress.txt"
    else:
        fileStr: str = dirStr + "/progress.txt"
    foundBool: bool = False
    with open(fileStr, "r") as fileObj:
        for lineStr in fileObj:
            lineStr = lineStr.replace("\n","")
            if lineStr == CONSTRUCT_HELICES_SUCCESS:
                foundBool = True
    fileObj.close()
    return foundBool

#=============================================================================
#=============================================================================

def check_convert_helices(dirStr: str):
    if dirStr[-1] == "/":
        fileStr: str = dirStr + "progress.txt"
    else:
        fileStr: str = dirStr + "/progress.txt"
    foundBool: bool = False
    with open(fileStr, "r") as fileObj:
        for lineStr in fileObj:
            lineStr = lineStr.replace("\n","")
            if lineStr == CONVERT_HELICES_SUCCESS:
                foundBool = True
    fileObj.close()
    return foundBool


def set_convert_helices(dirStr: str) -> bool:
    if dirStr[-1] == "/":
        fileStr: str = dirStr + "progress.txt"
    else:
        fileStr: str = dirStr + "/progress.txt"
    with open(fileStr, "a") as fileObj:
        fileObj.write(CONVERT_HELICES_SUCCESS + "\n")
    fileObj.close()

#=============================================================================
#=============================================================================

def set_copy_move_helices(dirStr: str):
    if dirStr[-1] == "/":
        fileStr: str = dirStr + "progress.txt"
    else:
        fileStr: str = dirStr + "/progress.txt"
    with open(fileStr, "a") as fileObj:
        fileObj.write(COPY_MOVE_HELICES_SUCCESS + "\n")
    fileObj.close()


def check_copy_move_helices(dirStr: str) -> bool:
    if dirStr[-1] == "/":
        fileStr: str = dirStr + "progress.txt"
    else:
        fileStr: str = dirStr + "/progress.txt"
    foundBool: bool = False
    with open(fileStr, "r") as fileObj:
        for lineStr in fileObj:
            lineStr = lineStr.replace("\n","")
            if lineStr == COPY_MOVE_HELICES_SUCCESS:
                foundBool = True
    fileObj.close()
    return foundBool

#=============================================================================
#=============================================================================

def set_top_geometry(dirStr: str):
    if dirStr[-1] == "/":
        fileStr: str = dirStr + "progress.txt"
    else:
        fileStr: str = dirStr + "/progress.txt"
    with open(fileStr, "a") as fileObj:
        fileObj.write(TOP_GEOMETRY_SUCCESS + "\n")
    fileObj.close()


def check_top_geometry(dirStr: str) -> bool:
    if dirStr[-1] == "/":
        fileStr: str = dirStr + "progress.txt"
    else:
        fileStr: str = dirStr + "/progress.txt"
    foundBool: bool = False
    with open(fileStr, "r") as fileObj:
        for lineStr in fileObj:
            lineStr = lineStr.replace("\n","")
            if lineStr == TOP_GEOMETRY_SUCCESS:
                foundBool = True
    fileObj.close()
    return foundBool


#=============================================================================
#=============================================================================

#This is used to create the initial progress file - this must always happen whilst building the directory structure
def create_progress_file(dirStr: str):
    if dirStr[-1] == "/":
        fileStr: str = dirStr + "progress.txt"
    else:
        fileStr: str = dirStr + "/progress.txt"

    with open(fileStr, "w") as fileObj:
        fileObj.write(";Individual progress status file\n")
    fileObj.close()


#checks the progress file
#def check_population_list(dirStr: str) -> bool:
#    return check_progress(dirStr, POPULATION_LIST_SUCCESS)

def check_minimisation_grompp(dirStr: str) -> bool:
    return check_progress(dirStr, MINIMISATION_GROMPP_SUCCESS)

def check_minimisation_mdrun(dirStr: str) -> bool:
    resultSuccessBool = check_progress(dirStr, MINIMISATION_MDRUN_SUCCESS)
    resultFltBool = check_progress(dirStr, MINIMISATION_FLT_ERROR)


def check_npt_eq_grompp(dirStr: str) -> bool:
    return check_progress(dirStr, NPT_EQ_GROMPP_SUCCESS)

def check_npt_eq_mdrun(dirStr: str) -> bool:
    return check_simulation_progress(dirStr, NPT_EQ_MDRUN_SUCCESS)
    #return check_progress(dirStr, NPT_EQ_MDRUN_SUCCESS)

def check_npt_production_grompp(dirStr: str) -> bool:
    return check_progress(dirStr, NPT_PRODUCTION_GROMPP_SUCCESS)

def check_npt_production_mdrun(dirStr: str) -> bool:
    return check_simulation_progress(dirStr, NPT_PRODUCTION_SUCCESS)
    #return check_progress(dirStr, NPT_PRODUCTION_MDRUN_SUCCESS)

def check_analysis(dirStr: str) -> bool:
    return check_progress(dirStr, ANALYSIS_SUCCESS)

def check_simulation_progress(dirStr: str, runTypeStr: str):
    if dirStr[-1] != "/":
        dirStr = dirStr + "/"

    if runTypeStr == NPT_EQ_MDRUN_SUCCESS:
        fileStr: str= dirStr + "npt_eq.gro"
    else:
        fileStr: str = dirStr + "npt_prod.gro"

    return os.path.isfile(fileStr)


#This is for any grompp - checks whether a tpr file exists and then writes that to a log file
def set_grompp(dirStr: str, tprFileStr: str, gromppSuccessStr: str, gromppFailStr: str, rerunBool: bool):
    if dirStr[-1] == "/":
        fileStr: str = dirStr + "progress.txt"
        tprFileStr: str = dirStr + tprFileStr
    else:
        fileStr: str = dirStr + "/progress.txt"
        tprFileStr: str = dirStr + "/" + tprFileStr

    if rerunBool == False:
        with open(fileStr, "a") as fileObj:
            if os.path.isfile(tprFileStr):
                fileObj.write(gromppSuccessStr + "\n")
            else:
                fileObj.write(gromppFailStr + "\n")
        fileObj.close()
    else:
        #if this is a continuation run I've got to clear the progress file of previous grompp efforts
        entryList = []
        with open(fileStr, "r") as fileObj:
            for lineStr in fileObj:
                lineStr = lineStr.replace("\n","")
                if lineStr != gromppSuccessStr or lineStr != gromppFailStr:
                    entryList.append(lineStr + "\n")
                else:
                    if os.path.isfile(tprFileStr):
                        entryList.append(fileObj.write(gromppSuccessStr + "\n"))
                    else:
                        entryList.append(gromppFailStr + "\n")
        fileObj.close()
        with open(fileStr, "w") as fileObj:
            for lineStr in entryList:
                fileObj.write(lineStr)
        fileObj.close()
        

#This is for any mdrun - check whether a particular entry was printed to a log file and whether the log file exists
def set_mdrun(dirStr: str, logFileStr: str, mdrunSuccessStr: str, mdrunFailStr: str, identifierStr: str, rerunBool: bool):
    if dirStr[-1] == "/":
        fileStr: str = dirStr + "progress.txt"
        logFileStr: str = dirStr + logFileStr
    else:
        fileStr: str = dirStr + "/progress.txt"
        logFileStr: str = dirStr + "/" + logFileStr

    mdrunStateBool: bool = False
    #if mdrun file does not exist then set to fail
    if os.path.isfile(logFileStr):
        #if mdrun didn't complete then set to fail
        with open(logFileStr, "r") as fileObj:
            for lineStr in fileObj:
                if lineStr.startswith(identifierStr) == True:
                    mdrunStateBool = True
                    break
        fileObj.close()
    else:
        mdrunStateBool = False

    #check whether the progress file already contains an mdrun entry of this kind
    #else set to pass
    if rerunBool == False:
        with open(fileStr, "a") as fileObj:
            if mdrunStateBool == False:
                fileObj.write(mdrunFailStr + "\n")
            else:
                fileObj.write(mdrunSuccessStr + "\n")
        fileObj.close()
    else:
        #if this is a continuation run I've got to clear the progress file of previous mdrun
        entryList = []
        with open(fileStr, "r") as fileObj:
            for lineStr in fileObj:
                lineStr = lineStr.replace("\n","")
                if lineStr != mdrunSuccessStr or lineStr != mdrunFailStr:
                    entryList.append(lineStr + "\n")
                #else:
                    #I haven't got a clue what I need to do here....
                    #I totally forgot!
		    #this is where a check to see whether the md run was successful
                    #if os.path.isfile(tprFileStr):
                    #    entryList.append(fileObj.write(mdrunSuccessStr + "\n"))
                    #else:
                    #    entryList.append(mdrunFailStr + "\n")
        fileObj.close()
        #reconstruct the progress file - but I need to know that the mdrun worked... I think I look at mdrunStateBool...
        with open(fileStr, "w") as fileObj:
            for lineStr in entryList:
                fileObj.write(lineStr)
            if mdrunStateBool == False:
                fileObj.write(mdrunFailStr + "\n")
            else:
                fileObj.write(mdrunSuccessStr + "\n")
        fileObj.close()


#this is used to check the progress of each part of the calculation by looking at the progress.txt file
#returns True if progressEntryStr is present
def check_progress(dirStr: str, progressEntryStr: str) -> bool:
    if dirStr[-1] == "/":
        fileStr: str = dirStr + "progress.txt"
    else:
        fileStr: str = dirStr + "/progress.txt"
    with open(fileStr, "r") as fileObj:
        for lineStr in fileObj:
            lineStr = lineStr.replace("\n","")
            if lineStr == progressEntryStr:
                fileObj.close()
                return True
    fileObj.close()
    return False


#this is to check whether the minimisation of the helix (after CG and before insertion) was energy minimised
def check_helix_minimisation(dirStr: str, minimisationTypeStr: str):
    filePathStr: str = dirStr + "/helix_minimisation_SD.gro"
    logFilePathStr: str = dirStr + "/helix_minimisation_SD.log"
    if minimisationTypeStr == MINIMISATION_SD_STR:
        resultBool: bool = os.path.exists(filePathStr)
    else:
        resultBool: bool = os.path.exists(filePathStr)

    #if the .gro file has been produced we need to then check what the maxiumum force is. 
    if resultBool == True:
        print("Checking helix max force from minimisation run.")
        with open(logFilePathStr) as fileObj:
            for lineStr in fileObj:
                if lineStr.startswith(MAXIMUM_FORCE_STR):
                    splitList = lineStr.split()
                    forceStr: str = splitList[3]
                    forceFlt: flt = float(forceStr)
        fileObj.close()
        if forceFlt > HELIX_MAX_FORCE_FLT:
            print("Maximum force is too great and exceeds " + str(HELIX_MAX_FORCE_FLT))
            resultBool = False
    else:
        print("helix_minimisation_SD.gro does not exist. Energy minimisation failed. No need to check the maximum force.")

    return(resultBool)

#==================================================
def record_individual_progress(fileStr, epocTime, individualInt, individualStr, fitnessScore, ppScore, plScore, llScore):
    with open(fileStr, "a") as fileObj:
        fileObj.write(str(epocTime) + "\t" +
            str(individualInt) + "\t" + 
            str(individualStr) + "\t" +
            str(fitnessScore) + "\t" +
            str(ppScore) + "\t" +
            str(plScore) + "\t" +
            str(llScore) + "\n")
    fileObj.close()



def save_GA_state(fileStr, mutatedNewPopulationList, t):
    with open(fileStr, "a") as fileObj:
        for entryStr in mutatedNewPopulationList:
            fileObj.write(str(entryStr) + "\t" + str(t) + "\n")
    fileObj.close()

def save_population_list(fileStr, populationList):
    with open(fileStr, "w") as fileObj:
        for entryStr in populationList:
            fileObj.write(str(entryStr) + "\n")
    fileObj.close()
