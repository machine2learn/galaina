
# Galaina

## Description
Within the course of the EU project Aggressotype, Machine2Learn BV developed a platform called Galaina to extract causal relationships from biological and medical databases that are in tabular format. Using the causal relationships from Galaina, we show that it is possible to derive subtypes, i.e., meaningful grouping of variables (here also called features) present in those databases. Galaina is a platform to read data, find causal associations, and generate a diagram that shows the causal relationships between features and subtypes. 

## Installation
### Linux (Ubuntu)

#### Terminal
Install v8 and some R packages

```bash
sudo apt install libv8-dev 
sudo apt install r-bioc-rbgl r-cran-devtools r-cran-dplyr r-cran-shiny r-cran-ggplot2 r-cran-tidyr r-cran-caret r-cran-nnet r-cran-mvtnorm r-bioc-graph r-cran-gridextra r-cran-psych
```

#### R console

```r
install.packages(c("daggity", "infotheo", "sbgcop"))

library(devtools)  # to use install_github
install_github(“cran/BDgraph@2.44”)  # Install version 2.44 of BDgraph

source("http://bioconductor.org/biocLite.R") 
biocLite("Rgraphviz")  # for visualizing causal discovery results (a graph) 

install.packages("latex2exp")  # for plot_consistency.R
install.packages("polycor")
install.packages("pcalg")
install.packages("gRain")
install.packages("bnlearn")
install.packages("ConfigParser")
install.packages("stringi")
install.packages("ggplotify")
```

## Usage

### Activation
In a terminal, execute the folloeing command:
```python
python main.py
```
Now the internal web-server is active and it can be used by a web browser (preferably Chrome or Firefox). The default URL is `http://127.0.0.1:5000/`.

### File formats
#### Data
Data must be privided as a CSV file

#### Factor models
We assume that each variable has one and only one factor with non-zero loading. 
    <!-- Variable_value = Loading_matrix(Variable, Factor) * Factor_value -->

Factor models are provided via a CSV file. The variables must be indicated with the column names of the input data file. 
A variable can be a factor only if no other variables have it as a factor.
The factor model file must satisfy one of the following formats:
1. Factor Loading Matrix
    CSV file s.t. each cell contain the factor loading for the related variable/factor combination.
    <!-- `cell_value = Loading_matrix(Variable, Factor)` -->
    First must contain the string `Variable`

    Example:
    | *Variable* | Factor_01 | Factor_02 | Factor_02 | ... |
    |---|---|---|---|---|
    |Var_01    |0.5 |0 | 0 | ...|
    |Var_02    |0.5 |0 | 0 | ...|
    |Var_03    |0 |0 | 0.25 | ...|

1. Factor Table
    CSV file with 3 columns listing the loading for each variable/factor combination. 
    Header is fixed and must containg the column names `Variable`, `Factor`, `Loading`.

    Example:
    | *Variable* | *Factor* | *Loading* |
    |---|---|---|
    |Var_01    | Factor_01 | 0.5 |
    |Var_02    | Factor_01 | 0.5 |
    |Var_03    | Factor_02 | 0.25 |
    |...| ... |...|

1. Factor Variable List
    CSV file with 2 columns listing for each factor the associated linear combinations of variables.
    Header is fixed and must containg the column names `Variable`, `Factor`, `Loading`.

    Example:
    | *Factor*| *Variable_list* |
    |---|---|---|
    | Factor_01 | Var_01 + Var_02 |
    | Factor_02 | Var_03  |
    |...| ... |



## Acknowledgement and Disclaimer 
This project has received funding from the European Union's Seventh Framework Programme for research, technological development and demonstration under grant agreement no 602805.
The research leading to these results has received funding from the European Community's Seventh Framework Programme (FP7/2007-2013) under Grant Agreement no 278948.
This project is co-funded by the Ambient Assisted Living (AAL) Joint programme, by the german BMBF, the French ANR, the Austrian BMVIT.
