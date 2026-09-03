# Heterocycle DFT Analysis

Computational workflow for the electronic-structure analysis of pyrrole, furan, thiophene, and pyridine.

## Method

Calculations use B3LYP-D3(BJ)/def2-TZVP with TightSCF, TightOpt, and DEFGRID3 in the gas phase. Each structure is optimized and followed by a harmonic frequency calculation.

The analysis reports HOMO and LUMO energies, HOMO-LUMO gaps, approximate ionization potentials and electron affinities derived from the frontier orbital energies, chemical hardness, and electronegativity.

## Requirements

- MATLAB
- ORCA 6.x

## Running the workflow

Set the ORCA executable path in `heterocycle_dft_workflow.m`:

```matlab
ORCA_EXE = 'C:\path\to\orca.exe';
```

Then run:

```matlab
heterocycle_dft_workflow
```

ORCA input and XYZ files are stored in `inputs/`. Raw calculation files and processed numerical results are written to `outputs/`. Figures are written to `figures/`.
