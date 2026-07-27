# Scripts

This directory contains the scripts used to execute the computationally intensive parts of the Bayesian analysis performed during this MSc research project.

The gravity-darkened transit modelling was developed using Jupyter notebooks; however, the full MCMC simulations were executed outside Jupyter on the University College London Observatory (UCLO) computing environment due to their long runtime.

## Contents

### `run_30_12_25.py`

Python script exported from the 05b_mcmc_cosine_parameterisation notebook.

This script:

- Loads the contamination-corrected light curve (`Data/processed/hat-p-70b_normalised.csv`).
- Defines the gravity-darkened transit model using **STARRY**.
- Implements the Brewer MCMC algorithm.
- Generates individual MCMC trace files for each simulation run.

### `run_mcmc.cmd`

Windows command script used to automate the MCMC execution.

The script repeatedly launches the Python program, allowing long simulations to be performed from the command line rather than inside Jupyter. During the original research, this workflow was executed on the UCLO computing environment, producing multiple trace files that were later combined into a single final MCMC chain.

## Computational Workflow

```text
Notebook 1
        │
        ▼
Processed light curve
        │
        ▼
Contamination correction
        │
        ▼
run_30_12_25.py
        │
        ▼
run_mcmc.cmd
        │
        ▼
Individual trace files
        │
        ▼
FINAL_30_12_25.csv
        │
        ▼
Trace analysis notebook
```

> **Note:** These scripts were developed for the software environment available at the University College London Observatory (UCLO). Reproducing the complete workflow may require recreating the original Python environment and installing compatible versions of the **STARRY** package and its dependencies.