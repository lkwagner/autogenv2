
# Autogen v2.

Goal: Provide tools that make it easy to manage and store QWalk runs.

#Object definitions

## Level 1: constructing input files and gathering data

 * The Writer object has all options to set up a job as member variables, and can construct input files for a run given those input files.
 * The Reader object can read the output of a run and report on its completion or lack thereof. This object defines what a successful run is. 
 * The Runner object executes a set of runs for a given code on a given computer system.


| Object | Functions                                               | Members                       | Notes                                                                                                                                                              |
|--------|---------------------------------------------------------|-------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Writer | __init__, set\_options(dictionary), is\_consistent(other) |                               | Write the input files and store the options for the calculation. Will typically have some version of write() or qwalk\_input, but we do not have a standard so far. |
| Reader | __init__, collect(list of outputfiles), write\_summary() | output: dictionary of outputs | Analyze the output and store the results of the calculation in its member variable output.                                                                         |
| Runner | __init__, check\_status(), run(inpfn,outfn)              | queue id                      | Execute the job, possibly using a queue system, and report on its completion.                                                                                      |


## Level 2: Manager

The Manager object consists of at least one Writer, Reader, and Runner. The Manager owns the directory for the time of its execution, that is, until it says it's done.
A manager may need multiple nextstep() calls before it's done. For example, one might want to run crystal and then properties using a queue to run both. 

Functions
 * __init__: Will be different for different managers, since they need different input objects.
 * nextstep() : Attempt to complete the next step in the calculation.
 * is\_consistent(other) : Check that the plan in other is consistent with the current plan.
 * write\_summary() : Print out some information about the calculation

## Level 3: Job (or Job Recipe)
 
A Job manages a sequence of managers.
These objects define a simulation protocol. 
For example, one might perform:
 * CRYSTAL DFT calculation
 * convert to QWalk format
 * variance optimize a Jastrow
 * perform diffusion Monte Carlo. 
The Job manages the Managers to attain this recipe. In this example, it would have four Managers to complete each of the steps, and coordinate the inputs and outputs to make sure that the QWalk converter knows where the CRYSTAL outputs are, etc.

The Job also has the responsibility of converting the raw data of the various managers into a coherent data structure. 

Functions:
 * __init__: options are set here. 
 * nextstep() : run the next step
 * status() : return the current status of the calculation
  

## Level 4: Job Ensemble

A job ensemble is a sequence of jobs. It manages the directory structure and recovering state. 


## Storage

Assuming that we coded the Managers, Writers, Readers, and Runners correctly, then we should be able to save these directly to JSON or pickle objects. A saved job ensemble is just a dump of a sequence of sequences of Managers.

#Job management

The above objects can be used in a program to do whatever run you'd like to do. However, typically there are some standard management tasks that can be handled with a job management program. 


## EnsembleRun (run\_plan.py)

 * Input: plan.pickle *[Comment: should we use pickle or another storage mechanism here?]* This file should have a sequence of jobs

EnsembleRun will then loop through jobs and:
 * Check to see if there is a directory called jobid
 * If the directory exists, read in status.pickle from the directory.
 * Update the status of the job.
 * For each element in the job, compare plan to status. This is where discrepancies and changes in the plan get resolved.
 * Determine what the next step in the job is and execute it

# Examples

`make_plan.py` will make a sample plan, and `run_plan.py` will run any plan and report on the results.

# To Do/ideas

To Do:
 * Bundling.
 * PySCF periodic boundary conditions.
 * Model fitting

