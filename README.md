# MERCI

# Molecular Ensemble Reordering sCrIpt

This Python script reorders a set of molecular structures (referred to as "ensemble") using a reference ensemble of molecular structures. The script aligns, reorders, and matches the molecules based on RMSD (Root Mean Square Deviation) criteria.

## Features

- **Reorder Molecular Ensembles**: The script takes two sets of molecular structures: a target ensemble and a reference ensemble, and attempts to find the best matching reference structure for each target structure.
- **Alignment and Sorting**: Uses centroid sorting, the Kabsch algorithm, and the Hungarian algorithm to reorder the atoms and minimize RMSD.
- **Command-line Interface**: Accepts input and output filenames as arguments, making it easy to integrate into workflows.

## Requirements

- Python 3.x
- Required Python libraries:
  - `numpy`
  - `scipy`

You can install these dependencies using pip:

```bash
pip install numpy scipy
