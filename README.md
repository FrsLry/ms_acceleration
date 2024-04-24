# Acceleration and demographic rates behind bird decline in North America

## Authors

*François Leroy, Petr Keil* - Department of Spatial Sciences, Faculty of Environmental Sciences, Czech University of Life Sciences Prague, Kamýcká 129, 16500 Praha-Suchdol, Czech Republic

*Marta Jarzyna* - Department of Evolution, Ecology and Organismal Biology, The Ohio State University, Columbus, Ohio, 43210, USA - Translational Data Analytics Institute, The Ohio State University, Columbus, Ohio, 43210, USA

## Description

This repository contains the scripts to reproduce the statistical analyses, figures, and a reproducible example to fit the dynamic *N*-mixture model for the Blue Jay. This is how it is structured:

* **Statistical analyses and figures:** the script `spatial_per_ species_analysis.R` contains the spatial analyses (maps) and the per species analyses. The script `relative_growth_rate_analysis.R` contains the spatial analysis for the growth rate relative to time 1 (Fig. 2c,d). The scripts `per_family_analysis.R` and `per_habitat_analysis.R` contain the per family and per 'preferred habitat' analyses, respectively. 

* **Reproducible example:** in order to fit the dynamic *N*-mixture model, simply run the script `example_blue_jay.R`. The species can be changed by changing the variable `species` in the script. Please note that it may take *ca.* 2 days for the model to fit. The data provided contains the input data for 564 species. 

## System requirements

All the dynamic *N*-mixture models were fit on the [Ohio Supercomputer Center ](https://www.oh-tech.org/) under Linux operating system with JAGS ver. 4.3.0. This repository was created on Windows 10, with R ver. 4.2.1.   

## License

Code and figures in this repository are under [CC-BY license](https://creativecommons.org/share-your-work/cclicenses/). 

If you use it, please give us credit by citing: https://github.com/FrsLry/ms_acceleration