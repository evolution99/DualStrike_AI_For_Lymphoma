#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DualStrike Docking Simulation Script
------------------------------------

This script performs the following tasks:
1. Prepares the NPM-ALK fusion protein sequence and saves it as a FASTA file.
2. Guides the user to submit the sequence to the AlphaFold server manually.
3. Processes the AlphaFold output and converts it to PDB format.
4. Prepares the DualStrike inhibitor structure from its SMILES notation.
5. Uses the SwissDock API to perform molecular docking between the protein and inhibitor.
6. Monitors the docking job and retrieves the results.
7. Analyzes and visualizes the docking results.

Author: [Your Name]
Date: [Date]
"""

# Section 1: Import Necessary Libraries
import os
import time
import requests
import tarfile
import subprocess
from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO, MMCIFParser
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Suppress warnings (optional)
import warnings
warnings.filterwarnings('ignore')

# Section 2: Prepare the NPM-ALK Fusion Protein Sequence

# Define the NPM-ALK fusion protein sequence
npm_alk_sequence = """>NPM_ALK_fusion
MEDSMDMDMSPLRPQNYLFGCELKADKDYHFKVDNDENEHQLSLRTVSLGAGAKDELHIVEAEAMNYEGSPI
KVTLATLKMSVQPTVSLGGFEITPPVVLRLKCGSGPVHISGQHLVAVEEDAESEDEEEEDVKDEVHGGKNKT
PSILPSDGLSRTCQPSNALEGKVTVRPDLSV"""

# Save the sequence to a FASTA file
fasta_file = "npm_alk_fusion.fasta"
with open(fasta_file, "w") as file:
    file.write(npm_alk_sequence)

print(f"FASTA sequence saved to {fasta_file}")

# Section 3: Instructions for Submitting to AlphaFold

print("\n--- AlphaFold Submission ---")
print("Please submit the sequence in 'npm_alk_fusion.fasta' to the AlphaFold Protein Structure Database manually.")
print("1. Go to https://alphafold.ebi.ac.uk/submit.")
print("2. Upload the 'npm_alk_fusion.fasta' file.")
print("3. Submit the job and wait for the prediction to complete.")
print("4. Download the predicted structure in CIF format (e.g., 'AF-NPM_ALK_fusion.cif').")
print("5. Save the file in this script's directory as 'npm_alk_fusion.cif'.")

input("\nPress Enter after you have downloaded 'npm_alk_fusion.cif' and saved it in this directory...")

# Section 4: Process the AlphaFold Output and Convert to PDB Format

# Define file paths
cif_file = "npm_alk_fusion.cif"
pdb_file = "npm_alk_fusion.pdb"

# Check if the CIF file exists
if not os.path.isfile(cif_file):
    print(f"\nError: {cif_file} not found. Please ensure you have downloaded the file from AlphaFold.")
    exit(1)
else:
    # Parse the CIF file and convert to PDB format
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("NPM_ALK_Fusion", cif_file)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)
    print(f"\nStructure converted to PDB format and saved as {pdb_file}")

# Section 5: Prepare the DualStrike Inhibitor Structure

# Define the SMILES string for DualStrike (hypothetical example)
dualstrike_smiles = "CC(C)(C)OC(=O)NC1=CC=C(C=C1)C(=O)NC2=NC=C(S2)C3=CN=C(N)N=C3N"

# Convert SMILES to RDKit molecule
dualstrike_mol = Chem.MolFromSmiles(dualstrike_smiles)
if dualstrike_mol is None:
    print("\nError: Invalid SMILES string for DualStrike.")
    exit(1)
else:
    # Add hydrogens
    dualstrike_mol = Chem.AddHs(dualstrike_mol)
    # Generate 3D coordinates
    AllChem.EmbedMolecule(dualstrike_mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(dualstrike_mol)
    # Save the molecule to SDF format
    ligand_file = "dualstrike.sdf"
    writer = Chem.SDWriter(ligand_file)
    writer.write(dualstrike_mol)
    writer.close()
    print(f"\nDualStrike structure saved to {ligand_file}")

# Section 6: Use the SwissDock Web Interface to Perform Docking

print("\n--- SwissDock Submission ---")
print("SwissDock does not provide a public API suitable for automation. We will use the web interface manually.")
print("1. Go to http://www.swissdock.ch/docking#.")
print("2. Upload the receptor (protein) file 'npm_alk_fusion.pdb'.")
print("3. Upload the ligand file 'dualstrike.sdf'.")
print("4. Set the docking parameters as needed (use default for this example).")
print("5. Submit the docking job.")
print("6. Note the 'Job ID' provided after submission.")

job_id = input("\nEnter the SwissDock Job ID: ")

print("\nPlease wait for the docking job to complete. This may take several hours.")
print("You can monitor the job status on the SwissDock website using your Job ID.")

input("\nPress Enter after the docking job has completed and you have downloaded the results (e.g., 'results.zip')...")

# Section 7: Analyze and Visualize the Docking Results

# Assuming the user has downloaded and extracted the results into a folder named 'SwissDock_results'

results_folder = input("\nEnter the path to the extracted SwissDock results folder (e.g., 'SwissDock_results'): ")

# Check if the folder exists
if not os.path.isdir(results_folder):
    print(f"\nError: Folder '{results_folder}' not found. Please ensure you have entered the correct path.")
    exit(1)

# Locate the 'clusters.dlg' file
clusters_file = os.path.join(results_folder, "clusters.dlg")

if not os.path.isfile(clusters_file):
    print("\nError: 'clusters.dlg' file not found in the results folder.")
    exit(1)
else:
    # Parse the clusters.dlg file to extract binding energies
    binding_energies = []
    with open(clusters_file, 'r') as f:
        for line in f:
            if "Estimated Free Energy of Binding" in line:
                try:
                    energy = float(line.strip().split()[-2])
                    binding_energies.append(energy)
                except ValueError:
                    continue
    if binding_energies:
        # Create a DataFrame
        df = pd.DataFrame(binding_energies, columns=['Binding Energy (kcal/mol)'])
        print("\nBinding Energy Statistics:")
        print(df.describe())
        # Plot a histogram
        plt.figure(figsize=(8,6))
        plt.hist(df['Binding Energy (kcal/mol)'], bins=20, color='blue', edgecolor='black')
        plt.title('Distribution of Binding Energies')
        plt.xlabel('Binding Energy (kcal/mol)')
        plt.ylabel('Frequency')
        plt.show()
    else:
        print("\nNo binding energies found in the 'clusters.dlg' file.")
