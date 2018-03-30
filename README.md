
# Autogen v2.

Goal: Provide tools that make it easy to manage and store QWalk runs.

# Getting started.

Check out the intro folder, which constains a set of jupyter notebooks to help you understand how this package works.

Working through all the notebooks should only take ~10 minutes.

# Object definitions

## Level 1: constructing input files and gathering data

 * The Writer object has all options to set up a job as member variables, and can construct input files for a run given those input files.
 * The Reader object can read the output of a run and report on its completion or lack thereof. This object defines what a successful run is. 
 * The Runner object executes a set of runs for a given code on a given computer system.

## Level 2: Manager

The Manager object consists of at least one Writer, Reader, and Runner. The Manager owns the directory for the time of its execution, that is, until it says it's done.
A manager may need multiple nextstep() calls before it's done. For example, one might want to run crystal and then properties using a queue to run both. 

Calling `nextstep` repeatedly will attempt to finish the calculation, one step at a time.
