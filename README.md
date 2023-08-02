# Modeling Earth's Top-of-Atmosphere Radiance Spectrum

This project involves modeling the top-of-atmosphere radiance spectrum measured by a nadir-viewing satellite sensor for various atmospheric conditions and surface types. The modeling incorporates radiative transfer processes in a layered model atmosphere.

## Model Overview

The project aims to simulate the spectral radiance observed by a satellite instrument pointing at Earth, considering:

- Atmospheric absorption, scattering, and transmission
- Surface reflection modeled with spectral albedo data  
- Impacts of cloud cover fraction
- Degradation of spectral resolution

The output is the modeled top-of-atmosphere albedo or radiance spectrum from 280 nm to 1000 nm.

## Code Overview

The core script is `cp3.py` which implements the radiative transfer model.

Key functions:

- `fracrad()` - Computes radiance based on input process flags 
- `light_scattered()` - Calculates Rayleigh scattering 
- `light_absorbed()` - Models gaseous absorption

Workflow:

1. Read in atmospheric profile, gas concentrations, surface albedo data
2. Calculate optical path length through model atmosphere layers 
3. Loop over wavelengths to compute radiative transfer
4. Output radiance spectrum under different assumptions

## Data

Input data files in the `data` folder:

- `ztp_eff.txt` - Layered atmospheric profile (height, temperature, pressure)
- `conc.txt` - Gas volume mixing ratio profiles (H2O, O3, O2) 
- `albedo.txt` - Spectral surface albedo data for different surfaces
- `solar.txt` - Solar spectral irradiance 
- `h2o.txt` - H2O absorption cross sections
- `o2.txt` - O2 absorption cross sections
- `o3.txt` - O3 absorption cross sections

## Methods

The radiative transfer model uses a layered plane-parallel atmosphere and iterates through the layers computing transmission, absorption, and scattering.

The key steps are:

- Set up model atmosphere layers (height, temperature, pressure) 
- Calculate layer optical depth using absorption/scattering coefficients
- Compute transmission of solar beam through atmosphere layers
- Calculate Rayleigh scattering source term in each layer
- Compute surface reflection using spectral Lambertian albedo 
- Propagate reflected radiance back through atmosphere
- Output top-of-atmosphere radiance for each wavelength

This is done for a range of wavelengths to produce the full spectral radiance.

Key assumptions:

- Atmosphere is in hydrostatic equilibrium
- Local thermodynamic equilibrium 
- No polarization effects
- Lambertian surface reflection

## Analysis

The script generates various plots analyzing the modeled radiance spectrum:

- Effect of atmospheric absorption and scattering
- Impact of spectral resolution 
- Effect of cloud cover fraction
- Variations for different surface types

These analyses provide insights into the factors that shape the top-of-atmosphere radiance spectrum of Earth.


