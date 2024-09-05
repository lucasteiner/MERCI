#!/bin/env python3
import numpy as np
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
import argparse

def parse_xyz(filename):
    """
    Parses an XYZ file and returns a list of molecular structures.

    Parameters:
    filename (str): The path to the XYZ file.

    Returns:
    list: A list of tuples, where each tuple contains the energy of the structure 
          and a list of atoms with their coordinates.
    """
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    molecules = []
    i = 0
    while i < len(lines):
        natoms = int(lines[i].strip())
        energy = float(lines[i + 1].strip().split()[0])
        atoms = [
            (parts[0], np.array([float(parts[1]), float(parts[2]), float(parts[3])]))
            for parts in (lines[i + 2 + j].split() for j in range(natoms))
        ]
        molecules.append((energy, atoms))
        i += 2 + natoms
    return molecules


def calculate_moments_of_inertia(atoms):
    """
    Calculates the moments of inertia for a set of atoms.

    Parameters:
    atoms (list): A list of atoms and their coordinates.

    Returns:
    numpy.ndarray: The moment of inertia tensor.
    """
    coords = np.array([atom[1] for atom in atoms])
    masses = np.array([12.0 if atom[0] == 'C' else 14.0 if atom[0] == 'N'
                       else 16.0 if atom[0] == 'O' else 1.0 if atom[0] == 'H' else 32.0 for atom in atoms])
    center_of_mass = np.sum(masses[:, np.newaxis] * coords, axis=0) / np.sum(masses)
    coords -= center_of_mass
    
    Ixx = np.sum(masses * (coords[:, 1]**2 + coords[:, 2]**2))
    Iyy = np.sum(masses * (coords[:, 0]**2 + coords[:, 2]**2))
    Izz = np.sum(masses * (coords[:, 0]**2 + coords[:, 1]**2))
    Ixy = -np.sum(masses * coords[:, 0] * coords[:, 1])
    Ixz = -np.sum(masses * coords[:, 0] * coords[:, 2])
    Iyz = -np.sum(masses * coords[:, 1] * coords[:, 2])
    return np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])


def calculate_rmsd(V, W):
    """
    Calculates the root mean square deviation (RMSD) between two sets of points.

    Parameters:
    V (numpy.ndarray): First set of coordinates.
    W (numpy.ndarray): Second set of coordinates.

    Returns:
    float: The RMSD value.
    """
    return np.sqrt(np.mean(np.sum((V - W)**2, axis=1)))


def kabsch(P, Q):
    """
    Implements the Kabsch algorithm to find the optimal rotation matrix that minimizes RMSD.

    Parameters:
    P (numpy.ndarray): First set of points.
    Q (numpy.ndarray): Second set of points.

    Returns:
    numpy.ndarray: The optimal rotation matrix.
    """
    C = np.dot(np.transpose(P), Q)
    V, S, W = np.linalg.svd(C)
    if (np.linalg.det(V) * np.linalg.det(W)) < 0.0:
        V[:, -1] = -V[:, -1]
    return np.dot(V, W)


def reorder_by_centroid(atoms):
    """
    Reorders atoms by their distance from the centroid.

    Parameters:
    atoms (list): A list of atoms and their coordinates.

    Returns:
    tuple: A tuple of the reordered atoms and their original indices.
    """
    coords = np.array([atom[1] for atom in atoms])
    centroid = np.mean(coords, axis=0)
    distances = np.linalg.norm(coords - centroid, axis=1)
    sorted_indices = np.argsort(distances)
    return [atoms[i] for i in sorted_indices], sorted_indices


def reorder_atoms_hungarian(reference_atoms, target_atoms):
    """
    Reorders atoms using the Hungarian algorithm based on minimal distance.

    Parameters:
    reference_atoms (list): A list of atoms from the reference structure.
    target_atoms (list): A list of atoms from the target structure.

    Returns:
    tuple: A tuple of reordered target atoms and the order indices.
    """
    ref_coords = np.array([atom[1] for atom in reference_atoms])
    target_coords = np.array([atom[1] for atom in target_atoms])
    cost_matrix = cdist(ref_coords, target_coords)
    row_ind, col_ind = linear_sum_assignment(cost_matrix)
    return [target_atoms[i] for i in col_ind], col_ind


def invert_positions(arr):
    """
    Inverts the positions of an array (i.e., mapping index to its position).

    Parameters:
    arr (numpy.ndarray): The array to invert.

    Returns:
    numpy.ndarray: The inverted array.
    """
    inverted_array = np.zeros_like(arr, dtype=int)
    inverted_array[arr] = np.arange(len(arr))
    return inverted_array


def find_best_sort_order_ensemble(ensemble_filename, reference_ensemble_filename):
    """
    Finds the best atom sort order based on RMSD between reference and target ensembles.

    Parameters:
    ensemble_filename (str): The XYZ file of the target ensemble.
    reference_ensemble_filename (str): The XYZ file of the reference ensemble.

    Returns:
    tuple: The reordered atoms, the best sort order, reference energy, and target energy.
    """
    ensemble_molecules = parse_xyz(ensemble_filename)
    reference_molecules = parse_xyz(reference_ensemble_filename)
    
    for ref_energy, ref_atoms in reference_molecules:
        ref_moi = calculate_moments_of_inertia(ref_atoms)
        
        for energy, atoms in ensemble_molecules:
            if abs(energy - ref_energy) > 1e-5:
                continue
            
            target_moi = calculate_moments_of_inertia(atoms)
            if np.linalg.norm(ref_moi - target_moi) > 1e1:
                continue
            
            ref_atoms_sorted, ref_sort_indices = reorder_by_centroid(ref_atoms)
            target_atoms_sorted, sort_indices = reorder_by_centroid(atoms)
            
            ref_coords_sorted = np.array([atom[1] for atom in ref_atoms_sorted])
            target_coords_sorted = np.array([atom[1] for atom in target_atoms_sorted])
            R = kabsch(ref_coords_sorted, target_coords_sorted)
            target_coords_aligned = np.dot(target_coords_sorted, R)
            
            reordered_atoms, hungarian_indices = reorder_atoms_hungarian(
                ref_atoms_sorted, [(atom[0], coord) for atom, coord in zip(target_atoms_sorted, target_coords_aligned)]
            )
            reordered_coords = np.array([atom[1] for atom in reordered_atoms])
            
            R_final = kabsch(ref_coords_sorted, reordered_coords)
            final_coords_aligned = np.dot(reordered_coords, R_final)
            
            rmsd = calculate_rmsd(ref_coords_sorted, final_coords_aligned)
            
            if rmsd < 0.1:
                combined_order = sort_indices[hungarian_indices[invert_positions(ref_sort_indices)]]
                return reordered_atoms, combined_order, ref_energy, energy

    return None, None, None, None


def reorder_molecule(molecule, combined_order):
    """
    Reorders the atoms of a molecule based on a given sort order.

    Parameters:
    molecule (tuple): A tuple containing the energy and atom list of the molecule.
    combined_order (numpy.ndarray): The atom reorder indices.

    Returns:
    tuple: The reordered molecule.
    """
    _, atoms = molecule
    reordered_atoms = [atoms[i] for i in combined_order]
    return molecule[0], reordered_atoms


def write_xyz(filename, molecules):
    """
    Writes the reordered molecular structures to an XYZ file.

    Parameters:
    filename (str): The output filename.
    molecules (list): A list of reordered molecules.
    """
    with open(filename, 'w') as file:
        for energy, atoms in molecules:
            file.write(f"{len(atoms)}\n")
            file.write(f"{energy}\n")
            for atom in atoms:
                file.write(f"{atom[0]} {atom[1][0]:.6f} {atom[1][1]:.6f} {atom[1][2]:.6f}\n")


def main():
    """
    Main function that handles command-line arguments and executes the reordering process.
    """
    parser = argparse.ArgumentParser(description='Reorder molecular ensembles using RMSD matching.')
    parser.add_argument('ensemble_filename', type=str, help='XYZ file for the ensemble of target structures.')
    parser.add_argument('reference_ensemble_filename', type=str, help='XYZ file for the ensemble of reference structures.')
    parser.add_argument('output_filename', type=str, help='XYZ file to write the reordered structures.')

    args = parser.parse_args()

    best_order_atoms, best_sort_order, ref_energy, target_energy = find_best_sort_order_ensemble(
        args.ensemble_filename, args.reference_ensemble_filename
    )

    if best_sort_order is not None:
        print(f"Match found with RMSD < 0.1. Reference energy: {ref_energy}, Target energy: {target_energy}")
        ensemble_molecules = parse_xyz(args.ensemble_filename)
        reordered_ensemble = [reorder_molecule(molecule, best_sort_order) for molecule in ensemble_molecules]
        write_xyz(args.output_filename, reordered_ensemble)
        print(f"Reordered ensemble written to {args.output_filename}")
    else:
        print("No suitable order found in the reference ensemble.")


if __name__ == '__main__':
    main()
