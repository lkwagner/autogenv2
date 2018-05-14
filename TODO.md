- [ ] Bad use case: accidentally putting the same writer for multiple managers
      will screw up the process because the first time the writer writes, it delares itself 'done'
- [ ] Manage path setup with setup script.
  - crystal, properties, pyscf, and qwalk executible names and locations.
- [ ] Comparable defaults between PySCF and Crystal are not the same.
- [ ] Cannot yet do different paths with QMC and crystal because the basis and the orb file needs to be in the same place.
  - May need to have the crystal manager export the results to the same directory (this means there are redundant files then).
- [ ] test case where the trial function gets called before the function itself.
- [ ] PySCF Si run in `simple_test` doesn't make any sense.
