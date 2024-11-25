# DualStrike Docking Simulation

This project performs a molecular docking simulation between the NPM-ALK fusion protein and the novel dual-site inhibitor, DualStrike, using AlphaFold for protein structure prediction and SwissDock for docking.

## **Table of Contents**

- [Prerequisites](#prerequisites)
- [Setup Instructions](#setup-instructions)
- [Running the Script](#running-the-script)
- [Detailed Steps](#detailed-steps)
  - [1. Prepare the NPM-ALK Fusion Protein Sequence](#1-prepare-the-npm-alk-fusion-protein-sequence)
  - [2. Submit Sequence to AlphaFold](#2-submit-sequence-to-alphafold)
  - [3. Process AlphaFold Output](#3-process-alphafold-output)
  - [4. Prepare the DualStrike Inhibitor Structure](#4-prepare-the-dualstrike-inhibitor-structure)
  - [5. Perform Docking with SwissDock](#5-perform-docking-with-swissdock)
  - [6. Analyze and Visualize Results](#6-analyze-and-visualize-results)
- [Understanding the Outputs](#understanding-the-outputs)
- [Troubleshooting](#troubleshooting)
- [References](#references)

---

## **Prerequisites**

Before you begin, ensure you have the following installed:

- Python 3.7 or higher
- pip (Python package installer)
- The following Python packages:
  - Biopython
  - RDKit
  - NumPy
  - pandas
  - matplotlib
  - requests

**Note:** RDKit installation can be non-trivial on some systems. Follow the official installation guide: [RDKit Install](https://www.rdkit.org/docs/Install.html)

---

## **Setup Instructions**

1. **Clone or Download the Repository:**

   ```bash
   git clone https://github.com/yourusername/dualstrike_docking.git
   cd dualstrike_docking
   ```
2. Create a Virtual Environment (Optional but Recommended):
  ```bash
   python3 -m venv venv
   source venv/bin/activate  # On Windows use: venv\Scripts\activate
   ```
3. Install the Required Packages:
   ```bash
   pip install biopython rdkit pandas numpy matplotlib requests
   ```
   Note: If you encounter issues installing RDKit via pip, refer to the official installation guide.

## Running the Script

### Making the Script Executable (Optional)
```bash
chmod +x dualstrike_docking.py
```

### Running Options
```bash
# Standard execution
python dualstrike_docking.py

# If made executable
./dualstrike_docking.py
```

## Detailed Steps

### 1. Prepare the NPM-ALK Fusion Protein Sequence
- Script automatically writes the sequence to `npm_alk_fusion.fasta`
- No manual action required

### 2. Submit Sequence to AlphaFold
*Note: Manual submission required due to lack of public API*

1. [Visit AlphaFold Protein Server website](https://alphafoldserver.com/)
2. Upload generated `npm_alk_fusion.fasta`
3. Submit and await completion
4. Download predicted structure (CIF format)
5. Save as `npm_alk_fusion.cif` in script directory

### 3. Process AlphaFold Output
- Script converts CIF to PDB format
- Ensure `npm_alk_fusion.cif` is in script directory
- Generates `npm_alk_fusion.pdb`

### 4. Prepare DualStrike Inhibitor Structure
- Uses RDKit for 3D structure generation
- Creates `dualstrike.sdf` from SMILES notation

### 5. Perform SwissDock Docking
*Note: Manual submission required*

1. Visit SwissDock Web Interface
2. Upload:
   - Receptor: `npm_alk_fusion.pdb`
   - Ligand: `dualstrike.sdf`
3. Configure parameters (defaults acceptable)
4. Submit and note Job ID
5. Wait for completion

### 6. Analyze and Visualize Results
1. Download results archive from SwissDock
2. Extract to folder (e.g., `SwissDock_results`)
3. Provide results folder path when prompted
4. Script generates statistics and visualizations

## Understanding the Outputs

### Binding Energy Statistics
- Provides descriptive statistics:
  - Mean
  - Standard deviation
  - Minimum/Maximum values

### Visualization
- Histogram of binding energy distribution

### Output Files
- `npm_alk_fusion.fasta`: Fusion protein sequence
- `npm_alk_fusion.pdb`: Protein structure
- `dualstrike.sdf`: Inhibitor structure
- Results folder: Contains all SwissDock output files

## Troubleshooting

### RDKit Installation Issues
```bash
# Alternative installation via conda
conda create -c conda-forge -n my-rdkit-env rdkit
conda activate my-rdkit-env
```

### Common Issues
- **Missing Files/Incorrect Paths**
  - Verify file names and locations
  - Check paths when prompted

- **Script Errors**
  - Check for missing dependencies
  - Verify file formats

- **Docking Issues**
  - Review SwissDock error messages
  - Verify file formatting

## References
- [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/)
- [SwissDock Docking Server](http://www.swissdock.ch/)
- [RDKit Documentation](https://www.rdkit.org/docs/)
- [Biopython Documentation](https://biopython.org/wiki/Documentation)
- [Matplotlib Documentation](https://matplotlib.org/stable/contents.html)

## Contact Information
- Author: Siraj Raval, o1-preview, Claude-3.5-Sonnet
- Email: hello@sirajrvaval.com

## Disclaimer
This script and associated instructions are provided for educational and research purposes. The DualStrike inhibitor is a hypothetical molecule used for demonstration. Always ensure compliance with all relevant laws, regulations, and ethical guidelines when conducting research.

## License
MIT LICENSE



   
