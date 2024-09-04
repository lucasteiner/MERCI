# MERCI

## Molecular Ensemble Reordering sCrIpt

This Python script reorders a set of molecular structures (referred to as "ensemble") using a reference ensemble of molecular structures. The script aligns, reorders, and matches the molecules based on RMSD (Root Mean Square Deviation) criteria.

### Features

- **Reorder Molecular Ensembles**: The script takes two sets of molecular structures: a target ensemble and a reference ensemble, and attempts to find the best matching reference structure for each target structure.
- **Alignment and Sorting**: Uses centroid sorting, the Kabsch algorithm, and the Hungarian algorithm to reorder the atoms and minimize RMSD.
- **Command-line Interface**: Accepts input and output filenames as arguments, making it easy to integrate into workflows.

### Requirements

- Python 3.x
- Required Python libraries:
  - `numpy`
  - `scipy`


### Usage

To use the script, you need to provide three command-line arguments:

    reference_ensemble.xyz: The filename of the ensemble of target structures (in XYZ format).
    target_ensemble.xyz:    The filename of the reference ensemble of structures (in XYZ format).
    reordered_ensemble.xyz: The filename where the reordered ensemble should be saved.

### Command Example
```bash
python reference_ensemble.xyz target_ensemble.xyz reordered_ensemble.xyz
```

### Script Logic

1.) Parse XYZ Files: The script reads the molecular structures from the provided XYZ files.
2.) Moment of Inertia Calculation: For each molecule, moments of inertia are calculated and used as a preliminary filter.
3.) Kabsch and Hungarian Algorithms: The Kabsch algorithm aligns the structures, while the Hungarian algorithm is used for optimal atom reordering based on RMSD minimization.
4.) Reordering: Once a match with an RMSD < 0.1 is found, the atom order is determined and applied to all structures in the target ensemble.
5.) Output: The reordered target ensemble is written to the specified XYZ file.



### XYZ Format

The input and output files must be in XYZ format, where each structure includes:

- The number of atoms on the first line.
- A comment line containing the energy.
- The atom symbols and coordinates for each atom.

4
-76.123456
C 0.000000 0.000000 0.000000
O 1.200000 0.000000 0.000000
H -0.800000 0.800000 0.000000
H -0.800000 -0.800000 0.000000

