Reproducible statistical modelling framework
=========================================================================================

<!-- badges: start -->

[![Ask Us Anything][0a]][0b]
![Open Source Love][0c]
[![license](https://img.shields.io/badge/license-MIT%20+%20file%20LICENSE-lightgrey.svg)](https://choosealicense.com/)

[0a]: https://img.shields.io/badge/Ask%20us-anything-1abc9c.svg
[0b]: https://github.com/open-aims/bcs_mixing_model/issues/new
[0c]: https://badges.frapsoft.com/os/v2/open-source.svg?v=103

## Content 

This repository contains the results presented in the paper “Predicting coral trends and attributing drivers of change from local to global scales” by Vercelloni et al. 

The files `Appendix_A.html` and `Appendix_B.html` correspond to the manuscript appendices. They can be reproduced following these simple instructions:

In the `Appendix_A` folder, the code is stored in the `scripts/` directory. The scripts should be run in sequence. The final script, `8.report.qmd`, automatically generates the visualizations in HTML format (e.g., `quarto render 8.report.qmd`). Note that the folder `figures` needs to be manually created. 

In the `Appendix_B` folder, the `appendix_B.qmd` file contains the code to run the statistical modelling pipeline for the four scenarios described in the paper. Once the pipeline has finished running, this file can be rendered to automatically visualize the results (e.g., `quarto render appendix_B.qmd`).

## Data tables 

Monitoring and disturbance data are provided in the `Appendix_A/data` folder, which contains three files:

- reef_data_aggregated.csv: observations of coral cover at data-tier locations. The table contains eight variables - `Tier5` (data-tier identifier), `fYEAR` (survey year), `COUNT` (count of coral cover), `TOTAL` (total point count), `n_obs` (number of monitored sites within the data-tier), `COVER` (percent coral cover), and geolocation coordinates of each data-tier. 
- hexpred.shp: a shapefile containing disturbance information across all tiers. The table contains 18 variables with `Tier5` (all-tier identifier), `fYEAR` (survey year), `svrty_d*` (severity of heat stress and associated time lags), `max_dhw*` (maximum values of degree heating weeks for the survey year and associated time lags), `svrty_c*` (severity of cyclone exposure and associated time lags), `max_cyc*` (maximum cyclone exposure for the survey year expressed as hours of exposure to cyclonic waves and associated time lags), `As_Data` (yes if `Tier5` is a data-tier, no is `Tier5` is a predictive-tier), `reefid` (reef(s) identifier within the `Tier5`), `reef_ar` (portion of the `Tier5` that contains coral reefs in square metres), `geometry` (geolocation of the `Tier5`).  
- MarineWaterBodiesV2_4: a shapefile indicating inshore and offshore reef zones along the Great Barrier Reef.

To cite the data, please use:
Australian Institute of Marine Science (AIMS). (2025). ReefCloud public data (Australia, 2006–2024). Retrieved from https://apps.aims.gov.au/metadata/view/1997de64-b274-427f-8cdc-d01b54e623f9
 (accessed 01-Dec-2025).

## License

This repository is provided by the authors under the MIT License ([MIT](http://opensource.org/licenses/MIT)).

## Bug reporting
* Please [report any issues or bugs](https://github.com/open-aims/namma_methods_comparison/issues).
