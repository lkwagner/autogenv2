import yaml
import os

print()
print("Setting up paths.")
print()
print("To use this, make sure you fill out 'paths.yaml'")
print("You'll need to know the location of Crystal, Crystal's properties component, QWalk,")
print("and the PySCF library if you want to use any of them in autogen.")
print("Use absolute paths where possible.")
print()

curpaths={}
keys=['crystal','Pcrystal','properties','Pproperties','qwalk','pyscf']

input_paths=yaml.load(open('paths.yaml'))

for key in sorted(keys):
  if key in input_paths and os.path.exists(input_paths[key]):
    path=input_paths[key]
    curpaths[key]=path
  else:
    print("Didn't find path for %s, either update it in paths.json, or avoid using it with autogen."%key)

with open('autopaths.py','w') as outf:
  outf.write('paths={}'.format(curpaths))

print("Done. autopaths.py should now be up-to-date.")
