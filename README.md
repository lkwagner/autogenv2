
# Autogen v2.

Goal: Remove `job_record` variables, and instead use member data in classes.

## Method heirarchy.

Top: `job_control.py` calls appropriate member functions for checking status of
jobs and taking appropriate action.

Mid: Writers and runner for Crystal and QWalk define operations that electronic
structure calculation requires.

Note: It might be useful to have a parent class for these classes.

Bottom: Submission classes are used by Runners to submit jobs. Submission
classes have a heirarchy of ConcreteSubmitter < MachineLevelSubmitter < LocalOrRemoteSubmitter.
