# Hotspots of acceleration and demographic processes behind decline of North American birds

## Authors

*François Leroy, Petr Keil* - Department of Spatial Sciences, Faculty of Environmental Sciences, Czech University of Life Sciences Prague, Kamýcká 129, 16500 Praha-Suchdol, Czech Republic

*François Leroy, Marta Jarzyna* - Department of Evolution, Ecology and Organismal Biology, The Ohio State University, Columbus, Ohio, 43210, USA

*Marta Jarzyna* - Translational Data Analytics Institute, The Ohio State University, Columbus, Ohio, 43210, USA

## Description

This repository contains the scripts to reproduce the statistical analyses, figures, a reproducible example to fit the dynamic *N*-mixture model for the Blue Jay, a simulation framework to demonstrate how this model can identify survival and recruitment from abundance data, and a correlative analysis between patterns of accelerating decline and covariates. This is how it is structured:

* **Dynamic N-mixture models:** The dynamic *N*-mixture model can be found in `.\model\jagsModel.txt` and can be used to reproduce the results. We created a reproducible example with the Blue Jay in `.\scripts\example_blue_jay.R`. The species can be changed by changing the variable `species` in the script. Please note that it may take *ca.* 2 days for the model to fit. The data provided contains the input data for 564 species. The summary of the models can be found in `.\summary_models\`. 

* **Propagating uncertainty in mixed models:** The second part of the analysis can be found in `.\mixed_model\`. There are 3 scripts to samples the chains at the routes/species, family and habitat levels. The mixed models can be fit by using the scripts in `.\mixed_model\script_fitting\`. The commented model is: `.\mixed_model\script_fitting\trend_model_N_spatial.R`.

* **Figures:** The figures can be recreated using the script `.\mixed_model\figures.R`. 

* **Simulations:** In `.\simulations\`, users can use a simulation framework, with known recruitment and survival parameters, and then use the dynamic N-mixture model (Dail & Madsen 2011) to correctly estimate the recruitment and survival parameters. Users are invited to test several values of parameters. Also accessible [here](https://frslry.github.io/ms_acceleration/simulations/simulations.html).

* **Coincidence between hotspots of acceleration and other variables:** We explored how environmental and anthropogenic variables correlate with the patterns of acceleration, change in abundance, recruitment rate, and loss rate using Random Forests and variable importance analyses. These analyses can be found in `.\interpreting_patterns\`. 

## System requirements

All the dynamic *N*-mixture models were fit on the [Ohio Supercomputer Center ](https://www.oh-tech.org/) under Linux operating system with JAGS ver. 4.3.0. This repository was created on Windows 10, with R ver. 4.2.1.   

## License

Code and figures in this repository are under [CC-BY license](https://creativecommons.org/share-your-work/cclicenses/). 

If you use it, please give us credit by citing our manuscript: https://doi.org/10.32942/X21032