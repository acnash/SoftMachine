#Make sure you are running the soft-machine conda environment. The scripts will not work otherwise. 

from gen_structure import *
from record_progress import *
from GA_stages import *
from GA_analysis import *
from typing import List
import os
import os.path
from os import path
import sys, getopt, shutil, time, random
import numpy
import pandas
import MDAnalysis
from MDAnalysis.analysis import contacts
from os import listdir
from os.path import isfile, join
from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy.stats import ttest_ind, shapiro, f_oneway, ranksums

divideLineStr: str = "--------------------------------------------------------------------------\n"
supressOutputBool: bool = True
popSizeInt: int = -1
noEpocsInt: int = -1
continuationBool: bool = False
newRunBool: bool = False
numMPI: int = 0
numOMP: int = 0
mutationRateFlt: float = 0

def helpDisplay():
    print("GA_de_novo.py [commands]\n")
    print("A new calculation will delete any existing calculation directories and files. You have 12 seconds to abort the run")
    print("before data is deleted. Always specify new and continuation commands but ensure they aren't the same. A continuation")
    print("run will step through the epocs and population folders providing those values have been specified and the directory")
    print("structure exists.")
    print("\n")
    print("-h\t--help\tHelp display (this).\n")
    print("COMMANDS")
    print("-v\t--verbose\t(Verbose reporting.")
    print("-n\t--new\t[y/n] New calculations starting from epoc 0. Explicitly required to avoid accidentally writing over previous steps.")
    print("-c\t--continue\t[y/n] Continue calculations. Useful if the system crashed/restarted/etc.")
    print("-e\t--epoc\tNumber of epocs. Always required, even for a continuation run.")
    print("-p\t--population\tPopulation size. Always required, even for a continuation run.")
    print("-m\t--ntmpi\tJust like gromacs. ADD MORE.")
    print("-o\t--ntomp\tJust like gromacs. ADD MORE.")
    print("-x\t--mutation\tThe mutation rate per residue. Between 0 and 1.")
    print("Examples:")
    print("python GA_de_novo.py -n y -c n -e 1 -p 1 -m 4 -o 8 -v y -x 0.05")
    print("Run a NEW calculation with a population of one design for one epoc. Set Gromacs ntmpi to 4 and ntomp to 8. Verbose output is on.")

argvList = sys.argv[1:]
try:
    opts, args = getopt.getopt(argvList,"h:v:n:c:e:p:m:o:x:", 
            ["help","verbose","new","continue","epoc","population","ntmpi","ntomp","mutation"])
except getopt.GetoptError:
    print("GA_de_novo.py [commands] -- use -h for command descriptions.")
    helpDisplay()
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-h","--help"):
        helpDisplay()
        sys.exit()
    if opt in ("-v", "--verbose"):
        if str(arg) == "y":
            print("Working in verbose mode.")
            supressOutputBool = False
    if opt in ("-n", "--new"):
        if str(arg) == "y":
            print("Working as a new calculation.")
            newRunBool = True
    if opt in ("-c", "--continue"):
        if str(arg) == "y":
            print("Working as a continuation (restart) run.")
            continuationBool = True
    if opt in ("-e", "--epoc"):
        print("Evolving for", str(arg), "epocs.")
        noEpocsInt = int(arg)
    if opt in ("-p", "--population"):
        print("Population size of", str(arg))
        popSizeInt = int(arg)
    if opt in ("-m","--ntmpi"):
        print("MPI:", str(arg))
        numMPI = int(arg)
    if opt in ("-o","--ntomp"):
        print("OMP:", str(arg))
        numOMP = int(arg)
    if opt in ("-x","--mutation"):
        print("Mutation rate:", str(arg))
        mutationRateFlt = float(arg)

if numMPI == 0 or numOMP == 0:
    print("Error: Multi-threading has not been set up.")
    sys.exit()
if noEpocsInt < 0:
    print("Error: Please provide the number of epocs.\n")
    helpDisplay()
    sys.exit()
if popSizeInt < 0:
    print("Error: Please provide the population size.\n")
    helpDisplay()
    sys.exit()
if continuationBool == newRunBool:
    print("Error: Current continuation and new are set to the same. \nPlease specify whether this is a new calculation or a continuation.\n")
    helpDisplay()
    sys.exit()

print("###################################################################")
print("#                                                                 #")
print("# .▄▄ ·       ·▄▄▄▄▄▄▄▄    • ▌ ▄ ·.  ▄▄▄·  ▄▄·  ▄ .▄▪   ▐ ▄ ▄▄▄ . #")
print("# ▐█ ▀. ▪     ▐▄▄·•██      ·██ ▐███▪▐█ ▀█ ▐█ ▌▪██▪▐███ •█▌▐█▀▄.▀· #")
print("# ▄▀▀▀█▄ ▄█▀▄ ██▪  ▐█.▪    ▐█ ▌▐▌▐█·▄█▀▀█ ██ ▄▄██▀▐█▐█·▐█▐▐▌▐▀▀▪▄ #")
print("# ▐█▄▪▐█▐█▌.▐▌██▌. ▐█▌·    ██ ██▌▐█▌▐█ ▪▐▌▐███▌██▌▐▀▐█▌██▐█▌▐█▄▄▌ #")
print("#  ▀▀▀▀  ▀█▄▀▪▀▀▀  ▀▀▀     ▀▀  █▪▀▀▀ ▀  ▀ ·▀▀▀ ▀▀▀ ·▀▀▀▀▀ █▪ ▀▀▀  #")
print("#                                                                 #")
print("#             By Anthony Nash (2021) - University of Oxford       #")
print("###################################################################")
print("\n")

for t in range(noEpocsInt):
    ###########################################################################################
    #Prepare directory structure
    ###########################################################################################
    parentDirStr: str = "/home/ubuntu/DE_NOVO_PROTEINS/"
    epocDirStr: str = parentDirStr + "EPOC_" + str(t) + "/"
    scriptsDirStr: str = parentDirStr + "SCRIPTS/"

    build_directory_structure(newRunBool, continuationBool, epocDirStr, popSizeInt)
    print(divideLineStr)

    ###########################################################################################
    #Build the initial population
    ###########################################################################################
    if t == 0:
        populationList = build_population(epocDirStr, popSizeInt, newRunBool, supressOutputBool)
        print(divideLineStr)

    ###########################################################################################
    #Construct the helix using Modeller
    ###########################################################################################
    construct_helix(populationList, epocDirStr, newRunBool)
    print(divideLineStr)

    ###########################################################################################
    #Convert to CG using martinize
    ###########################################################################################
    convert_helices_to_CG(populationList, epocDirStr, newRunBool, supressOutputBool)
    print(divideLineStr)

    ###########################################################################################
    #Copy the helix files (three times) and move them into place
    ###########################################################################################
    adjust_helix_position(populationList, newRunBool, epocDirStr, supressOutputBool)
    print(divideLineStr)

    ###########################################################################################
    #Write the TOP file, restraint file and construct the combined GRO file
    ###########################################################################################
    write_top_gro_files(popSizeInt, epocDirStr, scriptsDirStr, newRunBool, populationList, supressOutputBool)
    print(divideLineStr)

    ###########################################################################################
    #Performing steepest descent energy minimisation for epoc
    ###########################################################################################        
    run_energy_minimisation(parentDirStr, epocDirStr, popSizeInt, newRunBool, supressOutputBool, numMPI, numOMP, STEEPEST_DESCENT)
    print(divideLineStr)

    ###########################################################################################
    #Performing conjugate gradient energy minimisation for epoc
    ###########################################################################################
    run_energy_minimisation(parentDirStr, epocDirStr, popSizeInt, newRunBool, supressOutputBool, numMPI, numOMP, CONJUGATE_GRADIENT)
    print(divideLineStr)

    ###########################################################################################
    #Performing the 50ps equilibrium NPT simulations for initial epoc dt = 0.005
    ###########################################################################################
    run_equilibrium_npt(parentDirStr, epocDirStr, supressOutputBool, newRunBool, numMPI, numOMP, popSizeInt, SHORT_EQ)
    print(divideLineStr)

    ###########################################################################################
    #Performing the 50 ns equilibrium NPT simulations for initial epoc dt = 0.025
    ###########################################################################################
    run_equilibrium_npt(parentDirStr, epocDirStr, supressOutputBool, newRunBool, numMPI, numOMP, popSizeInt, LONG_EQ)
    print(divideLineStr)

    ###########################################################################################
    #Performing NPT production run for epoc
    ###########################################################################################
    run_production_npt(parentDirStr, epocDirStr, supressOutputBool, newRunBool, numMPI, numOMP, popSizeInt)
    print(divideLineStr)
    
    ###########################################################################################
    #Generate scores for individuals in epoc
    #The score list is an 1xN dimensional List. Each entry then has 4 entries (each a list).The
    #first of these entries is the t-score according to a comparison with a fourmer control.
    ###########################################################################################
    print("Calculating fitness scores.")
    scoreList = gen_final_frames(epocDirStr, popSizeInt)
    twoDimerScoreList = retrieve_score(scoreList)
    print(divideLineStr)
 
    ###########################################################################################
    #Select individuals to perform crossover
    #selectedIndividualList should only contain a list of integer IDs (0 to populationSizeInt-1)
    ###########################################################################################
    print("Selecting individuals to propogate.")
    selectedIndividualList = select_stochastic_population(populationList, twoDimerScoreList)
    print(divideLineStr)

    ###########################################################################################
    #Perform crossover
    ###########################################################################################
    print("Performing crossover.")
    newPopulationList = pair_individuals_for_crossover(selectedIndividualList)
    print(divideLineStr)

    ###########################################################################################
    #Mutate individuals at a 0.05 chance (5%) and replace the populations
    ###########################################################################################
    mutatedNewPopulationList = check_for_mutation(newPopulationList, mutationRateFlt)
    print(divideLineStr)
    
    print(newPopulationList)
    print(mutatedNewPopulationList)




