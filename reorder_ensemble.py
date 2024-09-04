import numpy as np
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
import argparse

# Function to parse XYZ file
def parse_xyz(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    molecules = []
    i = 0
    while i < len(lines):
        natoms = int(lines[i].strip())
        comment = lines[i + 1].strip()
        energy = float(comment.split()[0])
        atoms = []
        for j in range(natoms):
            parts = lines[i + 2 + j].split()
            atoms.append((parts[0], np.array([float(parts[1]), float(parts[2]), float(parts[3])])))
        molecules.append((energy, atoms))
        i += 2 + natoms
    return molecules

# Function to calculate moments of inertia
def calculate_moments_of_inertia(atoms):
    coords = np.array([atom[1] for atom in atoms])
    masses = np.array([12.0 if atom[0] == 'C' else 14.0 if atom[0] == 'N' else 16.0 if atom[0] == 'O' else 1.0 if atom[0] == 'H' else 32.0 for atom in atoms])
    center_of_mass = np.sum(masses[:, np.newaxis] * coords, axis=0) / np.sum(masses)
    coords -= center_of_mass
    Ixx = np.sum(masses * (coords[:, 1]**2 + coords[:, 2]**2))
    Iyy = np.sum(masses * (coords[:, 0]**2 + coords[:, 2]**2))
    Izz = np.sum(masses * (coords[:, 0]**2 + coords[:, 1]**2))
    Ixy = -np.sum(masses * coords[:, 0] * coords[:, 1])
    Ixz = -np.sum(masses * coords[:, 0] * coords[:, 2])
    Iyz = -np.sum(masses * coords[:, 1] * coords[:, 2])
    return np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])

# Function to calculate RMSD
def calculate_rmsd(V, W):
    return np.sqrt(np.mean(np.sum((V - W)**2, axis=1)))

# Kabsch algorithm to align two sets of points
def kabsch(P, Q):
    C = np.dot(np.transpose(P), Q)
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    U = np.dot(V, W)
    return U

# Function to reorder atoms based on centroid distances
def reorder_by_centroid(atoms):
    coords = np.array([atom[1] for atom in atoms])
    centroid = np.mean(coords, axis=0)
    distances = np.linalg.norm(coords - centroid, axis=1)
    sorted_indices = np.argsort(distances)
    return [atoms[i] for i in sorted_indices], sorted_indices

# Function to apply the Hungarian algorithm to reorder atoms
def reorder_atoms_hungarian(reference_atoms, target_atoms):
    ref_coords = np.array([atom[1] for atom in reference_atoms])
    target_coords = np.array([atom[1] for atom in target_atoms])
    cost_matrix = cdist(ref_coords, target_coords)
    row_ind, col_ind = linear_sum_assignment(cost_matrix)
    return [target_atoms[i] for i in col_ind], col_ind

# Modified function to handle an ensemble of reference structures
def find_best_sort_order_ensemble(ensemble_filename, reference_ensemble_filename):
    ensemble_molecules = parse_xyz(ensemble_filename)
    reference_molecules = parse_xyz(reference_ensemble_filename)
    
    for ref_energy, ref_atoms in reference_molecules:
        ref_moi = calculate_moments_of_inertia(ref_atoms)
        
        for energy, atoms in ensemble_molecules:
            if abs(energy - ref_energy) > 1e-5:
                continue
            
            target_moi = calculate_moments_of_inertia(atoms)
            if np.linalg.norm(ref_moi - target_moi) > 1e1:  # Updated threshold
                continue
            
            # Pre-sort by centroid distances
            ref_atoms_sorted, ref_sort_indices = reorder_by_centroid(ref_atoms)
            target_atoms_sorted, _ = reorder_by_centroid(atoms)
            
            # First Kabsch alignment
            ref_coords_sorted = np.array([atom[1] for atom in ref_atoms_sorted])
            target_coords_sorted = np.array([atom[1] for atom in target_atoms_sorted])
            R = kabsch(ref_coords_sorted, target_coords_sorted)
            target_coords_aligned = np.dot(target_coords_sorted, R)
            
            # Reorder using Hungarian algorithm
            reordered_atoms, hungarian_indices = reorder_atoms_hungarian(ref_atoms_sorted, [(atom[0], coord) for atom, coord in zip(target_atoms_sorted, target_coords_aligned)])
            reordered_coords = np.array([atom[1] for atom in reordered_atoms])
            
            # Second Kabsch alignment after reordering
            R_final = kabsch(ref_coords_sorted, reordered_coords)
            final_coords_aligned = np.dot(reordered_coords, R_final)
            
            rmsd = calculate_rmsd(ref_coords_sorted, final_coords_aligned)
            
            if rmsd < 0.1:
                print(f"RMSD: {rmsd}")
                combined_order = ref_sort_indices[hungarian_indices]
                return reordered_atoms, combined_order, ref_energy, energy

    return None, None, None, None

# Function to reorder a molecule based on the combined sort order
def reorder_molecule(molecule, combined_order):
    _, atoms = molecule
    reordered_atoms = [atoms[i] for i in combined_order]
    return molecule[0], reordered_atoms

# Function to write the reordered ensemble to an XYZ file
def write_xyz(filename, molecules):
    with open(filename, 'w') as file:
        for energy, atoms in molecules:
            file.write(f"{len(atoms)}\n")
            file.write(f"{energy}\n")
            for atom in atoms:
                file.write(f"{atom[0]} {atom[1][0]:.6f} {atom[1][1]:.6f} {atom[1][2]:.6f}\n")

# Main function to handle command-line arguments and execute the script
def main():
    parser = argparse.ArgumentParser(description='Reorder molecular ensembles using RMSD matching.')
    parser.add_argument('ensemble_filename', type=str, help='XYZ file for the ensemble of target structures.')
    parser.add_argument('reference_ensemble_filename', type=str, help='XYZ file for the ensemble of reference structures.')
    parser.add_argument('output_filename', type=str, help='XYZ file to write the reordered structures.')

    args = parser.parse_args()

    best_order_atoms, best_sort_order, ref_energy, target_energy = find_best_sort_order_ensemble(args.ensemble_filename, args.reference_ensemble_filename)

    if best_sort_order is not None:
        print(f"Match found with RMSD < 0.1. Reference energy: {ref_energy}, Target energy: {target_energy}")
        ensemble_molecules = parse_xyz(args.ensemble_filename)
        reordered_ensemble = [reorder_molecule(molecule, best_sort_order) for molecule in ensemble_molecules]
        write_xyz(args.output_filename, reordered_ensemble)
        print(f"Reordered ensemble written to {args.output_filename}")
    else:
        print("No suitable order found in the reference ensemble.")

# Entry point for script execution
if __name__ == '__main__':
    main()
