Reproducible study results
=========================================================================================

<!-- badges: start -->

[![Ask Us Anything][0a]][0b]
![Open Source Love][0c]
[![license](https://img.shields.io/badge/license-MIT%20+%20file%20LICENSE-lightgrey.svg)](https://choosealicense.com/)

[0a]: https://img.shields.io/badge/Ask%20us-anything-1abc9c.svg
[0b]: https://github.com/open-aims/bcs_mixing_model/issues/new
[0c]: https://badges.frapsoft.com/os/v2/open-source.svg?v=103

## Content 

This repository contains the results presented in the paper “Predicting coral trends and attributing drivers of change from local to broad spatial scales” by Vercelloni et al. 

The files `Appendix_A.html` and `Appendix_B.html` correspond to the manuscript appendices. They can be reproduced by rendering the quarto documents located in the `scripts/` directories. 

`Appendix_C.html` can be also reproduced by rendering the quarto document located in `Appendix_C` folder. 

## Data tables 

Monitoring and disturbance data are provided in the `Appendix_A/data` and `Appendix_B/data` folder, which contains these files:

- `reef_data_aggregated.csv`: observations of coral cover at data-tier locations. The table contains eight variables - `Tier4` (marine ecoregion id),`Tier5` (data-tier identifier), `fYEAR` (survey year), `COUNT` (count of coral cover), `TOTAL` (total point count), `LONGITUDE` (associated with Tier5 centroid), `LATITUDE` (associated with Tier5 centroid). 
- `covariates_full_tier5.RData`: a shapefile containing disturbance information across all tiers. The table contains 22 variables with `Tier5` (all-tier identifier), `year`, `severity_dhw*` (severity of heat stress and associated time lags), `max_dhw*` (maximum values of degree heating weeks for the survey year and associated time lags), `severity_c*` (severity of cyclone exposure and associated time lags), `max_cyc*` (maximum cyclone exposure for the survey year expressed as hours of exposure to cyclonic waves and associated time lags).  
- `tier5.sf.RData`: a shapefile indicating Tier5 locations.
- `tiers.lookup.RData`: a lookup table linking Tier5 with Tier4 levels. 
- `reef_layer.sf.RData`: a shapefile of the Tropical Coral Reefs of the World 

To cite the data, please use:
- `Appendix_A`: Australian Institute of Marine Science (AIMS). (2025). ReefCloud public data (Australia, 2006–2024). Retrieved from https://apps.aims.gov.au/metadata/view/1997de64-b274-427f-8cdc-d01b54e623f9
 (accessed 01-Dec-2025).
- `Appendix_B`: National Oceanographic and Atmospheric Association (NOAA). U.S. Pacific Islands National Coral Reef Monitoring Program (NCRMP). Retrived from https://www.fisheries.noaa.gov/inport/item/25274 (accessed 01-Dec-2025).

## License

This repository is provided by the authors under the MIT License ([MIT](http://opensource.org/licenses/MIT)).

## Bug reporting
* Please [report any issues or bugs](https://github.com/open-aims/namma_methods_comparison/issues).
