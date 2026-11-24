
# Info 
All code and data needed to replicate the analyses done in the manuscript: "Pipistrelle bat diet in urban and natural habitat is equally diverse but differs in composition", 
by Sander W. J. Buddendorf, Nina van Boxem, Judith E. Risse, Agata S. Pijl, Marcel E. Visser, Kamiel Spoelstra.

This repository used data after the 'bioinformatics' steps (QC, cleaning, classification of ASV's, training of classifier). All code of these steps can be found in /code/01_bioinformatics, 
but the raw sequences are stored elsewhere. These are uploaded to a ENA BioProject with ID: PRJEB116066.

# Installation
- Clone this repository
- Open `Urban_bat_diet_Rproject.Rproj` in RStudio (requires R version 4.3.2)
- Run the following to install all required packages at the exact versions used in this analysis:

```r
install.packages("renv")
renv::restore()
```

# Data
Before being able to run any code, the data needs to be downloaded from Zenodo (https://doi.org/10.34894/BLY8GN). Make sure all data is placed in an directory 'data' in the root.

# How to run the code
- Open the RProject (Urban_bat_diet_Rproject.Rproj)
- From the project, run /code/02_cleaning/complete_clean_pipeline.R 
- Then, run all scripts in /code/03_analyses_and_visualisations/
- All model outputs and visualisations needed, for both the main manuscript and the supplementary materials, will be stored in the /output directory
- ! Be aware that some scripts have a revised version (ending in '_R1')

# Metadata
## Explanation of variables / data format:
Information on how sample names are structured based on two examples:
AMV50828 - AMV: site | 5: roost ID within site | 0828: date, August 28th
A1000925 - A: site | 100: roost ID within site | 0925: date, September 25th

An exception to this rule are 
- the samples starting with LV, these are samples of a colleague scientist, standing for Laatvlieger (Serotine bat). These samples are not part of the BioClock project
- U1d0828: The 'd' stands for dak (which stands for roof) because this was a sample taken from a slightly different location within the study site
- L820815.1 and L820815.2: these were two DNA extractions from two differen subsamples of the same sample. Different pellets were used per extraction.
