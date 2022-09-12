# SoftMachine

## Download and setup

## Dependencies
Gromacs

Modeller

## Running the software
```
>python GA_de_novo.py [commands]

A new calculation will delete any existing calculation directories and files. You have 12 seconds to abort the run before data is deleted. Always specify new and continuation commands but ensure they aren't the same. A continuation run will step through the epocs and population folders providing those values have been specified and the directory structure exists.

-h  --help  Help display (this).

COMMANDS
-v  --verbose (Verbose reporting).
-n  --new [y/n] New calculations starting from epoc 0. Explicitly required to avoid accidentally writing over previous steps.
-c  --continue  [y/n] Continue calculations. Useful if the system crashed/restarted/etc.
-e  --epoc  Number of epocs. Always required, even for a continuation run.
-p  --population  Population size. Always required, even for a continuation run.
-m  --ntmpi Just like gromacs. ADD MORE.
-o  --ntomp Just like gromacs. ADD MORE.
-x  --mutation  The mutation rate per residue. Between 0 and 1.

Examples:
>python GA_de_novo.py -n y -c n -e 1 -p 1 -m 4 -o 8 -v y -x 0.05
Run a NEW calculation with a population of one design for one epoc. Set Gromacs ntmpi to 4 and ntomp to 8. Verbose output is on.
```
