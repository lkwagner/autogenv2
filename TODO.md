- [ ] Model fitting example.
- [ ] How to handle qwfiles better? Interaction between manager and trialfunc is currently not too flexible. For example it depends on kpoints being a list.
  - Possible to move the `export_qwalk` into `trialfunc`?
- [ ] Bad use case: accidentally putting the same writer for multiple managers
      will screw up the process because the first time the writer writes, it delares itself 'done'
- [ ] Comparable defaults between PySCF and Crystal are not the same.
- [ ] Cannot yet do different paths with QMC and crystal because the basis and the orb file needs to be in the same place.
  - May need to have the crystal manager export the results to the same directory (this means there are redundant files then).
- [ ] PySCF Si run in `simple_test` doesn't make any sense.
- [x] Manage path setup with setup script.
  - crystal, properties, pyscf, and qwalk executible names and locations.
 ---------------------------------------
- [ ] Redesign of data storage:
# Orbitals object
 - coeff: 2D array
 - basis: dictionary of labels

# System
 - Atomic positions 
 - Atomic labels
 - ECP 
 - nspin
 - k-point
 - lattice vectors
 
# Slater
 - occupation: [det,s,orb]
 - wt
 - orbital object
 
# Jastrow
 - Cusp
 - coefficients
 - basis
