<p align="center">
  <img src="Assets/Banner.png" width="100%">
</p>

![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=white)
![TESS](https://img.shields.io/badge/Data-TESS-red?-informational)
![MCMC](https://img.shields.io/badge/Method-MCMC-orange)
![Research](https://img.shields.io/badge/Type-MSc%20Research-purple)
![Status](https://img.shields.io/badge/Status-Completed-brightgreen)

# Gravity Darkening in Exoplanet Transits

A Python-based Astrophysics research project investigating how gravity-darkened transit photometry can be used to constrain the true three-dimensional spin–orbit geometry of the hot Jupiter system **HAT-P-70.**

This repository contains the code, analysis, and results developed during my MSc Astrophysics research project at University College London (UCL).

The project combines real TESS observations, statistical inference, and gravity-darkened transit modelling to investigate the orbital architecture of rapidly rotating exoplanet systems.


## Technical Skills

- Python programming
- Scientific data analysis
- Statistical modelling
- Bayesian parameter estimation (MCMC)
- Time-series analysis
- Data visualisation
- Model fitting and optimisation
- Scientific computing
- Research communication

## Technologies

- Python
- NumPy
- SciPy
- Matplotlib
- Lightkurve
- STARRY
- Jupyter Notebook
- TESS Space Telescope Data

## Why this project?

Although this research was conducted in Astrophysics, the analytical workflow is broadly applicable to data science. The project involved working with real observational data, developing reproducible Python workflows, performing statistical inference, validating competing models, and communicating results through scientific visualisation and technical documentation.

Many of these skills are directly transferable to data analysis, scientific computing, and machine learning roles.

## Project Status

**MSc project completed**

Repository currently being reorganised and documented for public release.

Additional documentation, cleaned notebooks and reproducibility notes will be added over time.

## Methodology

### 1. TESS Light Curve Preprocessing

The first stage of the project focuses on preparing high-quality TESS photometric observations for gravity-darkened transit modelling.

This notebook:

- Searches and downloads all publicly available TESS observations of HAT-P-70.
- Inspects the quality of each observing sector.
- Removes poor-quality measurements using the TESS quality flags.
- Stitches observations from multiple sectors into a single light curve.
- Detrends the light curve using the **Lightkurve** package.
- Folds and bins the transit signal.
- Exports clean datasets for the subsequent modelling and parameter estimation.

> **Note:** At the time this analysis was performed, TESS Sector 98 had not yet been released. The preprocessing therefore uses all publicly available observations available during the MSc project.

**Notebook:** [`01_tess_lightcurve_preprocessing.ipynb`](Notebooks/01_tess_lightcurve_preprocessing.ipynb)

### 2. Transit Light Curve Fitting

This notebook implements an analytical transit model based on the formulation presented by Haswell (*Exoplanet Transits*), incorporating quadratic limb darkening.

The workflow:

- Defines an analytical transit light curve model.
- Applies quadratic limb-darkening coefficients.
- Imports the preprocessed TESS light curve.
- Fits the analytical model to the observations.
- Optimises the transit parameters.
- Produces the best-fitting reference model used in the subsequent gravity-darkened analysis.

**Notebook:** [`02_transit_model_fitting.ipynb`](Notebooks/02_transit_model_fitting.ipynb)

### 3. Stellar Oblateness

This notebook derives the stellar parameters required for gravity-darkened transit modelling.

The workflow:

- Defines the stellar physical parameters.
- Computes the allowable stellar inclination range.
- Calculates the stellar rotation rate and oblateness.
- Investigates how the stellar shape varies with inclination.
- Estimates the polar and equatorial temperatures using gravity-darkening relations.
- Provides the stellar parameters used as inputs for the STARRY transit simulations.

**Notebook:** [`03_stellar_oblateness.ipynb`](Notebooks/03_stellar_oblateness.ipynb)

4. **Parameter Estimation**
   - Applied Markov Chain Monte Carlo (MCMC) methods to estimate the best-fitting system parameters.

5. **Model Evaluation**
   - Compared simulated and observed light curves using statistical goodness-of-fit metrics.

6. **Scientific Interpretation**
   - Derived the three-dimensional spin–orbit geometry of the planetary system and compared the results with published literature.