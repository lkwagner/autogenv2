
# Autogen v2.

Goal: Provide tools that make it easy to manage and store QWalk runs.


## Level 1: constructing input files and gathering data

 * The Writer object has all options to set up a job as member variables, and can construct input files for a run given those input files.
 * The Reader object can read the output of a run and report on its completion or lack thereof. This object defines what a successful run is. 
 * The Runner object executes a set of runs for a given code on a given computer system.

There may be some helper objects slightly below this level.

## Level 2: Managing completion of a job element

 * The Manager object consists of a Writer, Reader, and Runner. It may consist of several of those object. Its job is to complete the job according to the requirements set out in Reader.

## Level 3: Managing completion of a list of job elements


