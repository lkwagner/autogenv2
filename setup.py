import yaml

print("Setting up paths...")
curpaths=yaml.load(open('paths.yaml','r'))

for key in curpaths:
  path=input("Enter path for %s: "%key)
  curpaths[key]=path

yaml.dump(curpaths,open('paths.yaml','w'),default_flow_style=False)
print("Done. paths.yaml should now be up-to-date.")
