# MetID
## Description
This R package implements Baysian approach with pathway and network information to prioritize putative identifications of metabolites.

## Installation
You can install MetID package from GitHub with:
```r
# install.packages("devtools")
devtools::install_github("ressomlab/INDEED")
```

## Example
```r
library(MetID)
## check if colnames of dataset meet requirement
names(demo1)
## change colnames
colnames(demo1) <- c('query_m.z','name','formula','exact_m.z','pubchem_cid','kegg_id')
## get scores
out <- get_scores_for_LC_MS(demo1, type = 'data.frame', na='-', mode='POS')
```
