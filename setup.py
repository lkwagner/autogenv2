import os

print()
print("Setting up paths.")
print()
print("You'll need to know the location of Crystal, Crystal's properties component, QWalk,")
print("and the PySCF library if you want to use any of them in autogen.")
print("Use absolute paths.")
print()

curpaths={}
keys=['Pcrystal','Pproperties','qwalk','pyscf']

for key in sorted(keys):
  path=input("--Enter path for %s (empty to disable): "%key)
  if path=='': 
    path=None
  else:
    assert os.path.exists(path),\
        "Couldn't find '%s', you should check that it exists, or that it's spelled correctly."%path
  curpaths[key]=path

with open('autopaths.py','w') as outf:
  outf.write('paths={}'.format(curpaths))

print("Done. autopaths.py should now be up-to-date.")
