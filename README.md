Reproducible modelling framework
=========================================================================================

<!-- badges: start -->

[![Ask Us Anything][0a]][0b]
![Open Source Love][0c]
[![license](https://img.shields.io/badge/license-MIT%20+%20file%20LICENSE-lightgrey.svg)](https://choosealicense.com/)

[0a]: https://img.shields.io/badge/Ask%20us-anything-1abc9c.svg
[0b]: https://github.com/open-aims/bcs_mixing_model/issues/new
[0c]: https://badges.frapsoft.com/os/v2/open-source.svg?v=103

## Content 

This repository contains the results presented in the paper “Predicting coral trends and attributing drivers of change from local to global scales”

The files `Appendix_A.html` and `Appendix_B.html` correspond to the manuscript appendices. They can be reproduced following these simple instructions:

In the `Appendix_A` folder, the code is stored in the `scripts/` directory. The scripts should be run in sequence. The final script, `8.report.qmd`, automatically generates the visualizations in HTML format (e.g., `quarto render 8.report.qmd`). Note that the folder `figures` needs to be manually created. 

In the `Appendix_B` folder, the `appendix_B.qmd` file contains the code to run the modelling pipeline for the four scenarios described in the paper. Once the pipeline has finished running, this file can be rendered to automatically visualize the results (e.g., `quarto render appendix_B.qmd`).


## License

This repository is provided by the authors under the MIT License ([MIT](http://opensource.org/licenses/MIT)).

## Bug reporting
* Please [report any issues or bugs](https://github.com/open-aims/namma_methods_comparison/issues).
