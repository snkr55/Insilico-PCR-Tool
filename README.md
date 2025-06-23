# In-Silico PCR Simulation Tool

This Python-based tool simulates the Polymerase Chain Reaction (PCR) process *in silico*. It allows users to design, validate, and score primer pairs based on laboratory-like constraints, then simulate amplification from a target DNA region.

## Features

- Input DNA sequence via FASTA file
- Select a target region of interest
- Automatically generate forward and reverse primer candidates
- Validate primers based on:
  - Primer length (18–25 bp)
  - GC content (40–60%)
  - Melting temperature (Tm: 50–65°C)
  - Homopolymers
  - Tm compatibility (ΔTm ≤ 2°C)
  - Amplicon size constraints
- Score and rank primer pairs
- Simulate PCR and output amplicon details

## Requirements

- Python 3.10+
- Biopython
- pandas
- primer3-py
All required packages are listed in ```requirements.txt```.

## Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/snkr55/Insilico-PCR-Tool.git
cd Insilico-PCR-Tool
```

## Usage

### 1. Prepare the Environment:
Make sure Python and the following packages are installed:

```bash
pip install -r requirements.txt
```

This includes:
- biopython
- pandas
- primer3-py
- Python version: 3.10 or above

### 2. Prepare Input Files
Place your input FASTA file in the ```data/``` folder.

Ensure the sequence is in correct FASTA format with a single continuous sequence.

### 3. Run the Pipeline
From the root directory, run the tool via terminal:

```bash
python main.py
```

The script will:
- Prompt you to select a target region within the DNA
- Generate and validate primer candidates
- Score and display the top primer pairs
- Ask you to choose a primer pair
- Simulate in-silico PCR and output the amplicon

### 4. Output
Results will be saved to the ```results/``` directory:

- ```CandidateForwardPrimers_Data.csv```
- ```CandidateReversePrimers_Data.csv```
- ```ValidPrimerPairs_Data.csv```
- ```ScoredPrimerPairs_Data.csv```

PCR product details will be printed in the terminal.

## License
This project is open-source and free to use for academic or personal use.

## Contact and Support

For questions, suggestions, or collaboration, please connect on LinkedIn or email at sneha.050598@gmail.com
