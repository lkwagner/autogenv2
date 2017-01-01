import sys
sys.path.append("../")

import pickle as pkl

with open('plan.pickle','rb') as inpf:
  plan=pkl.load(inpf)

plan.nextstep()

print("\n\n\n###############################-- Summary \n")
report=plan.generate_report()

for rep in report:
  print(rep['id'],rep.keys())

import json
json.dump(report,open("report.json",'w'),indent=1)
