# In-Silico PCR Simulation Tool

This Python script performs an in-silico simulation of the Polymerase Chain Reaction (PCR) using a user-provided DNA template and primer sequences. It identifies potential amplicons based on primer binding, mismatch tolerance, and melting temperature (Tm) constraints.

## Features

* Accepts user input for DNA template and primer sequences.
* Calculates melting temperature (Tm) of primers using nearest-neighbor thermodynamics.
* Identifies primer binding positions allowing for a specified number of mismatches.
* Simulates PCR products (amplicons) based on primer match locations and product size constraints.
* Reports all valid amplicons along with their positions, lengths, and sequences.

## Requirements

* Python 3.x
* Biopython library

## Default Parameters

| Parameter      | Value       | Description                                   |
| -------------- | ----------- | --------------------------------------------- |
| Tm range       | 50–65°C     | Acceptable melting temperature for primers    |
| Product size   | 20–50 bases | Allowed size range for predicted PCR products |
| Max mismatches | 1           | Maximum mismatches allowed in primer binding  |

## Output

For each predicted PCR product, the script prints:

* Primer melting temperatures
* Start and end positions
* Length of the amplicon
* The PCR product (amplicon) sequence

## Example Output
![image](https://github.com/user-attachments/assets/24dbbb06-7d8f-4a8c-952f-f705af8bc245)

