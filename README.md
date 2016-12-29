
# Autogen v2.

Goal: Provide tools that make it easy to manage and store QWalk runs.

#Object definitions

## Level 1: constructing input files and gathering data

 * The Writer object has all options to set up a job as member variables, and can construct input files for a run given those input files.
 * The Reader object can read the output of a run and report on its completion or lack thereof. This object defines what a successful run is. 
 * The Runner object executes a set of runs for a given code on a given computer system.

There may be some helper objects below this level.

| Object | Functions                                               | Members                       | Notes                                                                                                                                                              |
|--------|---------------------------------------------------------|-------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Writer | __init__, set\_options(dictionary), is\_consistent(other) |                               | Write the input files and store the options for the calculation. Will typically have some version of write() or qwalk\_input, but we do not have a standard so far. |
| Reader | __init__, collect(list of outputfiles), write\_summary() | output: dictionary of outputs | Analyze the output and store the results of the calculation in its member variable output.                                                                         |
| Runner | __init__, check\_status(), run(inpfn,outfn)              | queue id                      | Execute the job, possibly using a queue system, and report on its completion.                                                                                      |


## Level 2: Managing completion of a job element

 * The Manager object consists of a Writer, Reader, and Runner. Its role is to complete the job element according to the requirements set out in Reader.

The Manager owns the directory for the time of its execution, that is, until it says it's done.

*[Comment: Maybe this should be ElementManager?]*

## Level 3: Job ensemble
 
Each Manager defines a job element, and a sequence of managers defines a job. A job ensemble is a sequence of jobs. 

*[Comment: we may wish to define a Job object that formalizes this and comes with an 'id' variable]*

## Storage

Assuming that we coded the Managers, Writers, Readers, and Runners correctly, then we should be able to save these directly to JSON or pickle objects. A saved job ensemble is just a dump of a sequence of sequences of Managers.

#Job management

The above objects can be used in a program to do whatever run you'd like to do. However, typically there are some standard management tasks that can be handled with a job management program. 


## EnsembleRun

 * Input: plan.pickle *[Comment: should we use pickle or another storage mechanism here?]* This file should have a sequence of jobs

EnsembleRun will then loop through jobs and:
 * Check to see if there is a directory called jobid
 * If the directory exists, read in status.pickle from the directory.
 * Update the status of the job.
 * For each element in the job, compare plan to status. This is where discrepancies and changes in the plan get resolved.
 * Determine what the next step in the job is and execute it

# Examples

`make_plan.py` will make a sample plan with h2 and n2, and `run_plan.py` will run any plan.

# To Do/ideas
Anyone who uses autogen is welcome to comment here.

To Do:
 * QWalk runs.
 * Make it use a queuing system/run parallel.
 * Export to database format.
 * Bundling.

Ideas:
 * Instead of `self.completed` as a T/F flag, have `self.status` as an informative string.
