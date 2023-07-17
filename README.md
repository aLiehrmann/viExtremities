# What is vizExtremities ?

The `vizExtremities` package laverages the Shiny framework to offer an interactive way to visualize the extremities (3' and 5' ends) of long read sequecing. Consequently, it enables researchers to gain insights into the structure of RNA isoforms with ease.


## How can I get vizExtremities ?

Make sure that `remotes` is installed by running
`install.packages("remotes")`, then type

``` r
remotes::install_github("aLiehrmann/vizExtremities")
```

## Quick start

Launch the Shiny application with 5Go of dedicated RAM.
``` r
library(vizExtremities)
vizExtremities(maxMemory=5)
```
