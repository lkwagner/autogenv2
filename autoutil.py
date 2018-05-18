#!/usr/bin/env python3
''' Simple utilities for interacting with autogen runs.'''

import pickle as pkl
import sys
import argparse

def get_info(pickle):
  print("Info about %s..."%pickle)
  man=pkl.load(open(pickle,'rb'))

  print("  Queue id: {}".format(man.runner.queueid))
  try:
    print("  Properties queue id: {}".format(man.prunner.queueid))
  except AttributeError:
    pass

if __name__=='__main__':

  parser=argparse.ArgumentParser("Autogen untilities.")
  parser.add_argument('manager',type=str,help='Pickle file to look at.')
  # Can add more options as needed.

  args=parser.parse_args()
  get_info(args.manager)
  
