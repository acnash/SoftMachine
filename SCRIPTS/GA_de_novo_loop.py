#Make sure you are running the soft-machine conda environment. The scripts will not work otherwise. 

from gen_structure import *
from record_progress import *
from GA_stages_loop import *
from GA_analysis import *
from typing import List
import os
import os.path
from os import path
import sys, getopt, shutil, time, random
import numpy
import pandas
import time
import MDAnalysis
from MDAnalysis.analysis import contacts
from os import listdir
from os.path import isfile, join
from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy.stats import ttest_ind, shapiro, f_oneway, ranksums

PARENT_DIR_STR: str = "/home/ubuntu/DE_NOVO_PROTEINS/"
SCRIPTS_DIR_STR: str = PARENT_DIR_STR + "SCRIPTS/"
DIVIDE_LINE_STR: str = "--------------------------------------------------------------------------\n"

outputRecordFileStr: str = PARENT_DIR_STR + "GA_fitness_analyses.txt"
stateFileStr: str = PARENT_DIR_STR + "GA_state.txt"

supressOutputBool: bool = True
popSizeInt: int = -1
noEpocsInt: int = -1
continuationBool: bool = False
newRunBool: bool = False
numMPI: int = 0
numOMP: int = 0
mutationRateFlt: float = 0
sleepTimeInt: int = 0

def helpDisplay():
    print("GA_de_novo_loop.py [commands]\n")
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
    print("-p\t--population\tPopulation size. Always required, even for a continuation run. This much be divisible by 2.")
    print("-m\t--ntmpi\tJust like gromacs. ADD MORE.")
    print("-o\t--ntomp\tJust like gromacs. ADD MORE.")
    print("-x\t--mutation\tThe mutation rate per residue. Between 0 and 1.")
    print("-d\t--debug\tA wait time of 25 second between each step. This is so one can read the output dialog at each step.")
    print("Examples:")
    print("python GA_de_novo_loop.py -n y -c n -e 3 -p 4 -m 4 -o 8 -v y -x 0.05")
    print("Run a NEW calculation with a population of four designs over two epocs. Set Gromacs ntmpi to 1 and ntomp to 4. Verbose output is on.")

argvList = sys.argv[1:]
try:
    opts, args = getopt.getopt(argvList,"h:v:n:c:e:p:m:o:x:d:", 
            ["help","verbose","new","continue","epoc","population","ntmpi","ntomp","mutation","debug"])
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
    if opt in ("-d","--debug"):
        sleepTimeInt = int(arg)
        print("Debug mode. Sleep time is " + str(sleepTimeInt) +  " seconds.")

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

if popSizeInt % 2 != 0:
    print("Error: Population size must be divisible by 2.\n")
    helpDisplay()
    sys.exit()

if continuationBool == True:
    if os.path.isfile(outputRecordFileStr):
        print("Analysis output file exists. Deleting for a new run.")
        os.remove(outputRecordFileStr)


populationList = [""] * popSizeInt
for t in range(noEpocsInt):
   
    i: int = 0
    #i is the population identifier

    if continuationBool == True:
        #remember to remove the continuation functionalist in future development
        continuationBool = False
        t, i = load_epoc_progress(PARENT_DIR_STR)
        i = i + 1

    
    while i < popSizeInt:
        ###########################################################################################
        #Prepare directory structure
        ###########################################################################################
        epocDirStr: str = PARENT_DIR_STR+ "EPOC_" + str(t) + "/"
        build_directory_structure_loop(newRunBool, continuationBool, epocDirStr, i)
        time.sleep(sleepTimeInt)

        ###########################################################################################
        #Build the initial population
        ###########################################################################################
        if t == 0:
            newDesignStr: str = build_population_loop(epocDirStr, i, newRunBool, supressOutputBool, populationList)
            populationList[i] = newDesignStr
            time.sleep(sleepTimeInt)


        #If t > 1, we need to check whether any entries in the populationList has existed before?
        #if not, continue as usual, else load in their earlier score and don't build or run a simulation



        ###########################################################################################
        #Construct the helix using Modeller
        ###########################################################################################
        construct_helix_loop(populationList[i], i, epocDirStr, newRunBool)
        time.sleep(sleepTimeInt)

        ###########################################################################################
        #Convert to CG using martinize
        ###########################################################################################
        convert_helices_to_CG_loop(populationList[i], i, epocDirStr, newRunBool, supressOutputBool)
        time.sleep(sleepTimeInt)

        ###########################################################################################
        #Aggressively energy minimise the CG single helix before inserting into the bilayer
        #--Check that the coordinates do not change too much such that the helix can't be inserted
        #--into the prepared bilayer.
        ###########################################################################################
        successBool = prepare_helix(i, epocDirStr, PARENT_DIR_STR, numMPI, numOMP, supressOutputBool)
        time.sleep(sleepTimeInt)
        if successBool == False:
            print("Failed to equilibrate the helix structure.")
            continue


        ###########################################################################################
        #Copy the helix files (three times) and move them into place
        ###########################################################################################
        adjust_helix_position_loop(i, newRunBool, epocDirStr, supressOutputBool)
        time.sleep(sleepTimeInt)

        ###########################################################################################
        #Write the TOP file, restraint file and construct the combined GRO file
        ###########################################################################################
        write_top_gro_files_loop(i, epocDirStr, SCRIPTS_DIR_STR, newRunBool, populationList[i], supressOutputBool)
        time.sleep(sleepTimeInt)

        ###########################################################################################
        #Performing steepest descent energy minimisation for epoc
        ###########################################################################################        
        run_energy_minimisation_loop(PARENT_DIR_STR, epocDirStr, i, newRunBool, supressOutputBool, numMPI, numOMP, STEEPEST_DESCENT)
        time.sleep(sleepTimeInt)

        ###########################################################################################
        #Performing conjugate gradient energy minimisation for epoc
        ###########################################################################################
        #run_energy_minimisation_loop(PARENT_DIR_STR, epocDirStr, i, newRunBool, supressOutputBool, numMPI, numOMP, CONJUGATE_GRADIENT)

        ###########################################################################################
        #Performing the 50ps equilibrium NPT simulations for initial epoc dt = 0.005 with restraints
        ###########################################################################################
        successBool = run_equilibrium_npt_loop(PARENT_DIR_STR, epocDirStr, supressOutputBool, newRunBool, numMPI, numOMP, i, SHORT_EQ_STRONG)
        time.sleep(sleepTimeInt)
        if successBool == False:
            print("Failed to run equilibrium with short time step.")
            continue
        

        ###########################################################################################
        #Performing the 50 ns equilibrium NPT simulations for initial epoc dt = 0.025 with restraints
        ###########################################################################################
        successBool = run_equilibrium_npt_loop(PARENT_DIR_STR, epocDirStr, supressOutputBool, newRunBool, numMPI, numOMP, i, LONG_EQ_STRONG)
        time.sleep(sleepTimeInt)
        if successBool == False:
            print("Failed to run equilibrium with long time step and strong restraints.")
            continue
        successBool = run_equilibrium_npt_loop(PARENT_DIR_STR, epocDirStr, supressOutputBool, newRunBool, numMPI, numOMP, i, LONG_EQ_MEDIUM)
        time.sleep(sleepTimeInt)
        if successBool == False:
            print("Failed to run equilibrium with long time step and medium restraints.")
            continue
        successBool = run_equilibrium_npt_loop(PARENT_DIR_STR, epocDirStr, supressOutputBool, newRunBool, numMPI, numOMP, i, LONG_EQ_WEAK)
        time.sleep(sleepTimeInt)
        if successBool == False:
            print("Failed to run equilibrium with long time step and weak restraints.")
            continue


        ###########################################################################################
        #Performing the 50 ns equilibrium NPT simulations for initial epoc dt = 0.025
        ###########################################################################################
        successBool = run_equilibrium_npt_loop(PARENT_DIR_STR, epocDirStr, supressOutputBool, newRunBool, numMPI, numOMP, i, LONG_EQ)
        time.sleep(sleepTimeInt)
        if successBool == False:
            print("Failed to run equilibrium with long time step.")
            continue


        ###########################################################################################
        #Performing NPT production run for epoc
        ###########################################################################################
        successBool = run_production_npt_loop(PARENT_DIR_STR, epocDirStr, supressOutputBool, newRunBool, numMPI, numOMP, i)
        time.sleep(sleepTimeInt)
        if successBool == False:
            print("Failed to run a full production run.")
            continue

        i = i + 1
        #save the current position in the epoc is case of a continuation run
        save_epoc_progress(PARENT_DIR_STR, t, i)
        save_population_list(epocDirStr + "population_list.txt", populationList)
        #END of preperation loop - all the sequences should work

    ###########################################################################################
    #Generate scores for individuals in epoc
    #The score list is an 1xN dimensional List. Each entry then has 4 entries (each a list).The
    #first of these entries is the t-score according to a comparison with a fourmer control.
    ###########################################################################################
    print("Calculating fitness scores.")
    scoreList = gen_final_frames(epocDirStr, popSizeInt)
    twoDimerScoreList = retrieve_score(scoreList)

    ###########################################################################################
    #Record population progress to file for analysis
    ###########################################################################################
    print("Record population progress.")
    for i in range(popSizeInt):
        individualStr: str = populationList[i]
        fitnessScore = twoDimerScoreList[i]
        ppScore = (scoreList[i])[1]
        plScore = (scoreList[i])[2]
        llScore = (scoreList[i])[3]
        record_individual_progress(outputRecordFileStr, t, i, individualStr, fitnessScore, ppScore, plScore, llScore)


    ###########################################################################################
    #Select individuals to perform crossover
    #selectedIndividualList should only contain a list of integer IDs (0 to populationSizeInt-1)
    ###########################################################################################
    print("Selecting individuals to propogate.")
    selectedIndividualList = select_stochastic_population(populationList, twoDimerScoreList)

    #I think this is a bit of a problem (selectedIndividualList)

    ###########################################################################################
    #Perform crossover
    ###########################################################################################
    print("Performing crossover.")
    newPopulationList = pair_individuals_for_crossover(selectedIndividualList)

    ###########################################################################################
    #Mutate individuals at a 0.05 chance (5%) and replace the populations
    ###########################################################################################
    print("Performing mutation")
    mutatedNewPopulationList = check_for_mutation(newPopulationList, mutationRateFlt)
   
    ###########################################################################################
    #Save the current state of the system
    #save_GA_state - append the master_GA document with the new population list 
    #save_population_list - saves the population list responsible for this epoc to the epoc root
    #
    #WARNING: a continuation run will not know how to handle the GA_state.txt - or any of these 
    #
    ###########################################################################################
    print("Saving GA state and population list")
    save_GA_state(stateFileStr, mutatedNewPopulationList, t)
    #save_population_list(epocDirStr + "population_list.txt", populationList)


    ###########################################################################################
    #Swap old population for new population
    ###########################################################################################
    populationList = mutatedNewPopulationList
    print("epoc finished")



