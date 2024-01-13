# Summary
These are MCMC samplers, data simulators, and processing and run scripts for a occupancy-abundance model for abundance-mediated species interactions with a full example case study. This model framework models detection/non-detection data to estimated occupancy, abundance, and interactions between the species of interest. This model enables the user to apply detection/non-detection data collected over repeat surveys to model interactions as a function of abundance. These samplers are presented in Twining et al. submitted, and are based on adaptations of the [Waddle et al. (2010)](https://www.jstor.org/stable/25680391) formulation for modelling species interactions within an occupancy model, but instead of modelling the state model of a subordinate species a function of the occupancy states, it is modelled as a function of abundance (N). We provide a range of MCMC samplers for different ecological scenarios between two or more species including disease- and predator- mediated competition, intraguild predation, mesopredator release, and tri-trophic cascades.

# The working directory

Below you fill find descriptions of each folder in this repository and files contained within them.

# Models and data simulators (./occupancy_abundance_model/models_and_simulations)
This folder has seven R scripts. Within each script are MCMC samplers and a data simulator for different iterations of the occupancy-abundance model. Model code and implentation is conducted in the nimble package [(de Valpine et al. 2022).](https://cran.r-project.org/web/packages/nimble/index.html) 

**1. two_species_occupancy_abundance_data_simulator_and_model.R**

This script contains a data simulator and model for simulating a two species occupancy abundance model with spatially varying interaction terms. The code is commented throughout to describe each part of the simulator and script.

**2. three_species_occupancy_abundance_abundance_model_and_simulator_pinemarten_greysquirrel_redsquirrel_example.R**

This script contains a data simulator and model for simulating a three species occupancy abundance model - modeling the occupancy of the dominant species on abundance of intermediate, and abundance of subordinate, with additional interactions between abundance of intermediae with the subordinate. The code is commented throughout to describe each part of the simulator and script.

**3. three_species_abundance_abundance_abundance_model_and_simulator_otter_urchin_kelp_example.R**

This script contains a data simulator and model for simulating a three species abundance model - modeling the abundance of the dominant species on abundance of intermediate, and abundnace of subordinate, with additional interactions between abundance of intermediae with the subordinate. The code is commented throughout to describe each part of the simulator and script.

**4. three_species_abundance_occupancy_occupancy_model_and_simulator_deer_liverfluke_moose_example.R**

This script contains a data simulator and model for simulating a three species abundance model - modeling the abundance of the dominant species on occupancy of intermediate, and occupancy of subordinate, with additional interactions between occupancy of intermediae with the subordinate. The code is commented throughout to describe each part of the simulator and script.

**5. occupancy_vs_abundance_mediated_interactions_model_and_simulator.R**

This script contains a data simulator and model for simulating a two species occupancy abundance model with and fits models with interactions as a function of occupancy and as a function of abundance. The code is commented throughout to describe each part of the simulator and script.

**6. occupancy_abundance_model_three_species_simulator_and_model_variable_detection.R**

This script contains a data simulator and model for simulating a three species abundance model - modeling the abundance of the dominant species on abundance of intermediate, and occupancy of subordinate, with additional interactions between abundance of intermediae with the subordinate. The code is commented throughout to describe each part of the simulator and script.

**7. two_species_occupancy_abundance_simulator_spatially_varying_interaction_terms** 

This script contains a data simulator and model for simulating a two species occupancy abundance model with spatially varying interaction terms. The code is commented throughout to describe each part of the simulator and script.


# The coyote-fisher-marten case study folder (./occupancy_abundance_model/case_study)

This folder has nine files.

**1. alNY_2013-2021_6kmbuffer_allsitecovs.csv**

This file contains all of the summarized spatial covariate data at a 6km scale used in the analysis. Each row is a different 6km 2 pixel in New York State, each column is a covariate, and each cell is a value.

The covariates used in the analysis.

| **Covariate**    | **Description**                                                                              | **Source**                                                                                                                                       |
|------------------|----------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------|
| Deciduous        | Proportion of a 6 km2 buffer around the detector made up of deciduous forest                 | NLCD, 2019 (https://www.mrlc.gov/data/nlcd-2019-land-cover-conus)                                                                                |
| Coniferous       | Proportion of a 6 km2 buffer around the detector made up of coniferous forest                | NLCD, 2019 (https://www.mrlc.gov/data/nlcd-2019-land-cover-conus)                                                                                |
| Mixed            | Proportion of a 6 km2 buffer around the detector made up of mixed forest                     | NLCD, 2019 (https://www.mrlc.gov/data/nlcd-2019-land-cover-conus)                                                                                |
| Pasture          | Proportion of a 6 km2 buffer around the detector made up of pasture                          | NLCD, 2019 (https://www.mrlc.gov/data/nlcd-2019-land-cover-conus)                                                                                |
| Cultivated.Crops | Proportion of a 6 km2 buffer around the detector made up of cultivated crops                 | NLCD, 2019 (https://www.mrlc.gov/data/nlcd-2019-land-cover-conus)                                                                                |
| road_density     | Mean number of km of road per km2 in each 6 km2 buffer                                       | calculated from primary and secondary roads raster provided by the NYSDEC, hosted on the github                                                   |
| snow_depth       | Mean daily snow depth(m) of the 6 km2 buffer around each detector across the sampling period | National Operational Hydrologic Remote Sensing Centre, 2004. Snow data assimilation system (SNODAS) products (https://doi.org/10.7265/N5TB14TC). |
| forest_edge      | Edge density of combined class of all forest in each 6 km2 buffer around detectors           | NLCD, 2019 (https://www.mrlc.gov/data/nlcd-2019-land-cover-conus)                                                                                |
| GPP              | The amount of carbon captured by plants (kg C MJ-1) in at the detectors                      | MODIS Land Satellite, 2017 (https://lpdaac.usgs.gov/products/mod17a2hv006/)                                                                      |
| Camera           | The camera trap model used at each site (Bushnell, Recoynx, Browning)                        | Fieldwork datasheets                                                                                                                             |

**2.allNY__2013-2021_15kmbuffer_allsitecovs.csv**

This file contains all of the summarized spatial covariate data at a 15km scale used in the analysis. Each row is a different 15km 2 pixel in New York State, each column is a covariate, and each cell is a value.

The covariates used in the analysis.

| **Covariate**    | **Description**                                                                               | **Source**                                                                                                                                       |
|------------------|-----------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------|
| Deciduous        | Proportion of a 15 km2 buffer around the detector made up of deciduous forest                 | NLCD, 2019 (https://www.mrlc.gov/data/nlcd-2019-land-cover-conus)                                                                                |
| Coniferous       | Proportion of a 15 km2 buffer around the detector made up of coniferous forest                | NLCD, 2019 (https://www.mrlc.gov/data/nlcd-2019-land-cover-conus)                                                                                |
| Mixed            | Proportion of a 15 km2 buffer around the detector made up of mixed forest                     | NLCD, 2019 (https://www.mrlc.gov/data/nlcd-2019-land-cover-conus)                                                                                |
| Pasture          | Proportion of a 15 km2 buffer around the detector made up of pasture                          | NLCD, 2019 (https://www.mrlc.gov/data/nlcd-2019-land-cover-conus)                                                                                |
| Cultivated.Crops | Proportion of a 15 km2 buffer around the detector made up of cultivated crops                 | NLCD, 2019 (https://www.mrlc.gov/data/nlcd-2019-land-cover-conus)                                                                                |
| road_density     | Mean number of km of road per km2 in each 15 km2 buffer                                       | calculated from primary and secondary roads raster provided by the NYSDEC, hosted on the githu                                                   |
| snow_depth       | Mean daily snow depth(m) of the 15 km2 buffer around each detector across the sampling period | National Operational Hydrologic Remote Sensing Centre, 2004. Snow data assimilation system (SNODAS) products (https://doi.org/10.7265/N5TB14TC). |
| forest_edge      | Edge density of combined class of all forest in each 15 km2 buffer around detectors           | NLCD, 2019 (https://www.mrlc.gov/data/nlcd-2019-land-cover-conus)                                                                                |
| GPP              | The amount of carbon captured by plants (kg C MJ-1) in at the detectors                       | MODIS Land Satellite, 2017 (https://lpdaac.usgs.gov/products/mod17a2hv006/)                                                                      |
| Deer             | The probability of occupancy (ψ) of white-tailed deer in a 15 km2 buffer around each detector | Calculated from detection/non-detection data and covariate hosted in the repository (see below)                                                  |
| Camera           | The camera trap model used at each site (Bushnell, Recoynx, Browning)                         | Fieldwork datasheets                                                                                                                             |

**3. allNY_2013-2021_7dayocc_americanmarten_detection_nondetection.csv**

This file contains detection/non-detection data for American marten between the years of 2013-2021 in New York state. Each row is a site, each column is an occasion, each cell is a detection/non-detection record.

**4. allNY_2013_2021_7dayocc_coyote_counts.csv**

This file contains the count data for coyote between the years of 2013-2021 in New York state. Each row is a site, each column is an occasion, each cell is a detection/non-detection record.

**5. allNY_2013-2021_7dayocc_fisher_counts.csv**

This file contains the count data for fisher between the years of 2013-2021 in New York state. Each row is a site, each column is an occasion, each cell is a detection/non-detection record.

**6. 7dayoccDetectionswhitetaileddeerNZ.csv**

This file contains the detection/non-detection data for white-tailed deer between the years of 2016-2018 in the northern zone of New York state. Each row is a site, each column is an occasion, each cell is a detection/non-detection record.

**7. juliandays_allNY_2013-2021.csv**

This file contains the sampling dates (in ordinal format) for each sampling occasion for the winter surveys from 2013-2021. Each row is a site, each column is an occasion.

**8. 'coyote_fisher_marten_abu_ocu_binomialobsmodel_nimble_model.R'**

This is the nimble occupancy-abundance model that is fit to the coyote-fisher-marten data files above. The code is commented out to describe each part of the model.

**9. nimble_occu_abu_coyote_fisher_marten_binomialversion_processingandruncode.R**

This is the data loading, formatting, and run script for the case study three species occupancy-abundance model. he code is commented out to describe each stage of the process. 

