# Railway Localization and Mapping

This repository contains: 

1. a localization filter implementation based on the fusion of IMU and GPS data, especially tailored to railway vehicles, and 
2. a method to create track-maps based on the results of the before mentioned localization filter. 

All functions have been realized with Matlab/Simulink 2019b.

A detailed description of the implementation can be found in my PhD thesis. As soon as it will be published, I will provide a link here. For the moment I have to refer you to my publications: https://orcid.org/0000-0002-0429-1787. Unfortunately, they don't cover all the implementation details found in this repository.

**Descriptions of the image below:** Localization result of the approach presented here "new approach" compared to a "common" EKF approach. The new approach is characterized by smaller uncertainties in the cross-track direction of the track when the traveled track geometry has been detected.

![Railway Localization: Example Image](https://github.com/hans42/localization/blob/master/plots/err_ellip_exam_standalone_C.png)

# Usage

- Matlab integration
  - Clone repository and add it to your Matlab path
- Loading recorded IMU and GPS data
  - You can use any recorded IMU and GPS data which follows the [LRT Data Sharing Principles v1.2](LrtDataSetGuidelines-v1_2.pdf)
  - In my thesis I used the datasets mentioned below: see [Datasets](#data-sets)
  - Put the data somewhere on your Matlab path
  - Adjust the paths to the data in the import scripts under: _./_helperScripts/import/_ (see also [Notes.1](#notes))
- Configuration
  1. Adjust the localization filter and mapping settings in the config files (or use the initial setup as it is): _./_config/_ (see also [Notes.1](#notes))
  2. Only if you use a dataset not mentioned below: Adjust switch case statement in initLocalization.m script to load your import scripts and configuration files
- Localization:
  - Open simulink model localization.slx
  - Set simulation time manually
  - During loading the simulink model the initialization script "initLocalization.m" is executed automatically via a callback
  - Run the simulation
  - There are several "to workspace" blocks in the simulink model which transfer the simulation results into the Matlab workspace for further processing
- Mapping
  - After running the localization (by running the simulink model "localization.slx") you can start the "runMapOptimization.m" script to create a track map from the positioning data calculated from the simulink model
  - After finishing the mapping process the results are available in the Matlab workspace and you are asked, if you want to store the results in a \*.mat to be able to easily load the data again later on 
- Plots
  - In _./plots/_ you find some scripts to plot the localization and mapping results (see also [Notes.2](#notes))

## Notes
1. For the [Datasets](#data-sets) mentioned below, there are already corresponding import and config scripts you can use as a template. The datasets are abbreviated in the filenames via the letters "BS", "C" and "NT".
2. Before the localization and mapping results can be plotted some processing steps are performed. They are carried out by the scripts under _./_helperScripts/preparation_. The plot script automatically starts these processing scripts. At the end of some processing steps you are asked, if you want to store the processed data in a \*.mat file, so that you can easily load them later on again.
3. Some general functionalities to [work with datasets](https://github.com/hans42/railway-dataset-functions) and my [Railway Toolbox](https://github.com/hans42/railway-toolbox) are integrated via git [submodules](https://git-scm.com/book/de/v2/Git-Tools-Submodule). Details on these submodules can be found in the corresponding repositories.
 
# Datasets
  
The implemenation has been used with these datasets so far:
  
- https://doi.org/10.25534/tudatalib-360
- https://doi.org/10.48328/tudatalib-166.3
- https://doi.org/10.25534/tudatalib-359

But it should be possible to use any recorded IMU and GPS data which follows the [LRT Data Sharing Principles v1.2](LrtDataSetGuidelines-v1_2.pdf)
