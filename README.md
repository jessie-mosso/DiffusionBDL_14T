# DiffusionBDL_14T
Processing codes associated with: In vivo diffusion MRS and MRI reveal cerebellar microstructure alterations in the rat developing brain during hepatic encephalopathy, Imaging Neuroscience, 2026

Please cite the above should you use any material from this repository.

Jessie Mosso, NYU, 11/02/2026

## Third party softwares:
- FID-A, Jamie Near: https://github.com/CIC-methods/FID-A
- LCModel

## Files description: 
> 1_Processing

is used to process STEAM/STELASER data and includes sub functions needed. 

Raw data can be found here:

> 2_LCModel_quantification
 
contains the control files and basis sets for LCModel quantification of STEAM/STE-LASER. The 4 excel tables contain summaries. 


> 3_Diffusion_fitting
> 
contains code to read the output of LCModel, fit diffusion models and compute statistics
