# Data

This directory contains the datasets produced and used throughout the analysis.

## Processed

Contains the processed TESS light curves.

- `hat-p-70b.csv` — Output of Notebook 1 after downloading, cleaning, detrending, folding, and binning the TESS observations.
- `hat-p-70b_normalised.csv` — Revised light curve with contamination from a nearby star removed, following the methodology described in the dissertation. This dataset is used for the gravity-darkened transit modelling and MCMC analysis.

## Traces

Contains the individual MCMC trace files generated during the Bayesian parameter estimation. Each file represents one segment of the Markov chain Monte Carlo.

## Final

Contains the concatenated MCMC chain produced from the individual trace files. This dataset is used to generate the posterior distributions and corner plots presented in the dissertation.