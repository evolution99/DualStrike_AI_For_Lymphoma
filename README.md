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
