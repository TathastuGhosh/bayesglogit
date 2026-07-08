[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-3.5%2B-blue.svg)](https://www.r-project.org/)
[![Status](https://img.shields.io/badge/Status-Active-success.svg)]()

# bayesglogit

## Bayesian Grouped Logistic Regression with Shrinkage Priors

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![R](https://img.shields.io/badge/R-3.5%2B-blue.svg)](https://www.r-project.org/)

`bayesglogit` implements three grouped shrinkage priors for logistic regression using Pólya-Gamma data augmentation:

- **BGHS** — Bayesian Grouped Horseshoe
- **GRASP** — Grouped Regression with Adaptive Shrinkage Priors\
- **GIGG** — Group Inverse-Gamma Gamma

## Installation

``` r
# Install from GitHub
remotes::install_github("TathastuGhosh/bayesglogit")
```

## Quick Example

``` r
library(bayesglogit)

# Fit BGHS model
fit <- bghs(y ~ ., data = my_data, groups = my_groups)
summary(fit)
```

## Methods

| Method    | Description                                       |
|-----------|---------------------------------------------------|
| `bghs()`  | Bayesian Grouped Horseshoe                        |
| `grasp()` | Grouped Regression with Adaptive Shrinkage Priors |
| `gigg()`  | Group Inverse-Gamma Gamma                         |

## Reference

Ghosh, T. (2026). *Extending the Idea of Bayesian Grouped Regression under Logistic Setup*. MSc Dissertation, University of Calcutta.

## License

MIT License © 2026 Tathastu Ghosh
