# Galaina
Galaina is a platform to read data, find probabilistic causal associations, and generate a diagram that shows the causal relationships among the data variables or factors (i.e. group of variables expressing the same hypothetical constructs).
 
## Description
Within the course of the EU project Aggressotype, Machine2Learn BV developed a platform called Galaina to extract probabilistic causal relationships from biological and medical databases that are in tabular format. Using the causal relationships extracted by Galaina, we could identify the variable subtypes, i.e., meaningful grouping of variables (here also called features) present in those databases.

## Installation

### Python and R Path Settings
This part is common to all the operating systems.

#### Python

We assume that Python 3 is the default python interpreter. By running from the terminal the following command, all the required packages will be installed:
```bash
pip install -r requirements.txt
```

#### Local setup
The path to the `Rscript` executable is read from the file `meta_config.ini`:
```ini
[r_front_end]
path_r_binary_command = /usr/local/bin/Rscript
```
Please assign the path to your local `Rscript` file to the key `path_r_binary_command` in the `r_front_end` section  of the file `meta_config.ini`.

### Linux (Ubuntu)

#### Terminal
Install v8 and some R packages

```bash
sudo apt install libv8-dev
sudo apt install r-bioc-rbgl r-cran-devtools r-cran-dplyr r-cran-shiny r-cran-ggplot2 r-cran-tidyr r-cran-caret r-cran-nnet r-cran-mvtnorm r-bioc-graph r-cran-gridextra r-cran-psych
```

#### R console

```r
install.packages(c("dagitty", "infotheo", "sbgcop"))
library(devtools)  # to use install_github
install_github("cran/BDgraph@2.44")  # Install version 2.44 of BDgraph
source("http://bioconductor.org/biocLite.R") 
biocLite("Rgraphviz")  # for visualizing causal discovery results (a graph) 
install.packages(c("latex2exp", "polycor", "pcalg", "gRain", "bnlearn", "ConfigParser", "stringi", "ggplotify"))
install.packages(c("here"))
```


### MacOS 

#### Terminal
Install v8 via Homebrew
```bash
brew install v8
```
<!-- brew install v8@3.15 -->

#### R console

```r
install.packages("V8")
source("http://bioconductor.org/biocLite.R") 
biocLite("RBGL")
install.packages(c("devtools", "dplyr", "shiny", "ggplot2", "tidyr", "caret", "nnet"), dependencies=TRUE)
install.packages("mvtnorm")
install.packages(c("dagitty", "infotheo", "sbgcop"))
library(devtools)  # to use install_github
install_github("cran/BDgraph@2.44")  # Install version 2.44 of BDgraph
source("http://bioconductor.org/biocLite.R") 
biocLite("Rgraphviz")  # for visualizing causal discovery results (a graph) 
install.packages(c("latex2exp", "polycor", "pcalg", "gRain", "bnlearn", "ConfigParser", "stringi", "ggplotify"))
install.packages("here")
```


## Usage

### Activation of the Internal Web Server

In a terminal, execute the following command (Python 3 is the assumed python interpreter):
```python
python main.py
```
### Access the Internal Web Server

1. Now the internal web-server is active and it can be accessed by a web browser (preferably Chrome or Firefox). The default URL is `http://127.0.0.1:5000/`. 
1. Sign in with the following credentials:
   * Username: `test`
   * Password: `test12345`
1. Upload the data and factor model files, edit the configurations accordingly, and run.


## File formats

### Input Data
Galaina accepts as input one or multiple datasets.

#### Single dataset
Each dataset must be provided as CSV files where:
* The first column contains the instances (i.e. rows) identifiers.
* The first rows contains a header with the names of the variables (also known as indicators). Basically they are the column names.
* The following cell values are interpreted as NaN (empty string '' is the recommended choice): '', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN', '-NaN', '-nan', '1.#IND', '1.#QNAN', 'N/A', 'NA', 'NULL', 'NaN', 'n/a', 'nan', 'null'.

##### Example:

We refer to this example dataset as **Dataset_01**

| *Identifier* | *Variable1_01* | *Variable1_02* | *Variable1_03* | ... |
| ------------ | -------------- | -------------- | -------------- | --- |
| id01         | 45             | 0              | 1              | ... |
| id02         | 0.5            | 42             | 0              | ... |
| id03         | -8             | 0              | 0.25           | ... |
| id04         | 48             | 10             | 5.25           | ... |

#### Multiple datasets
If  multiple datasets are given as input, they will be inner-joined with respect to the first column. 
This means that Galaina will create internally a new database made by only the instances present in all the input datasets and with all the columns of the input datasets.

##### Example:

Together with Dataset_01, we give as input two more datasets; the value of `Variable3_02` for `id02` in `Dataset_03` is missing:

**Dataset_02**

| *Identifier* | *Variable2_01* | *Variable2_02* | *Variable2_03* | ... |
| ------------ | -------------- | -------------- | -------------- | --- |
| id02         | 4.5            | 442            | 4              | ... |
| id03         | 48             | 4              | 4.25           | ... |
| id04         | 448            | 40             | 9.25           | ... |

**Dataset_03**

| *Identifier* | *Variable3_01* | *Variable3_02* | *Variable3_03* | ... |
| ------------ | -------------- | -------------- | -------------- | --- |
| id01         | 75             | 3              | 3              | ... |
| id02         | 3.5            |                | 0.3            | ... |
| id04         | 38             | 30             | 3.25           | ... |
| id05         | 58             | 35             | 5.25           | ... |

Then Galaina will process the resulting dataset, made only by `id02` and `id04`:

| *Identifier* | *Variable1_01* | *Variable1_02* | *Variable1_03* | *Variable2_01* | *Variable2_02* | *Variable2_03* | *Variable3_01* | *Variable3_02* | *Variable3_03* | ... |
| ------------ | -------------- | -------------- | -------------- | -------------- | -------------- | -------------- | -------------- | -------------- | -------------- | --- |
| id02         | 0.5            | 42             | 0              | 4.5            | 442            | 4              | 3.5            |                | 0.3            | ... |
| id04         | 48             | 10             | 5.25           | 448            | 40             | 9.25           | 38             | 30             | 3.25           | ... |


### Input Factor models
We assume that each variable has one and only one factor with non-zero loading. 
    <!-- Variable_value = Loading_matrix(Variable, Factor) * Factor_value -->

Factor models are provided via a CSV file. The variables must be indicated with the column names of the input data file. 
If a variable has no factor associated, a dummy factor must be present in the factor model: that factor must be called as the variable itself and will have just that variable associated; the factor loading will be 1.

<!-- **Example:** -->
##### Factor Model Example

All the example tables of this section are based on the following factor model. This factor model for six variables `Var01`,...`Var06` (indicators) and two factors `Factor02` and `Factor03`:
* `Var01` has no factor associated.
* `Var02` and `Var03` have factor `Factor02` with loadings 0.75 and 0.25, respectively.
* `Var04`, `Var05`, and `Var06` have factor `Factor03` with loadings 5, 25, and 25, respectively.

<!-- A variable can be a factor only if no other variables have it as a factor. -->
The factor model file must satisfy one of the following formats:


#### 1. Factor Loading Matrix

CSV file such that each cell contains the factor loading for the related variable/factor combination.
<!-- `cell_value = Loading_matrix(Variable, Factor)` -->
The first cell must contain the string `Variable`.

<!-- **Example:** -->
##### Example

Factor Loading Matrix representation of [factor model example](#factor-model-example) above.

| *Variable* | Var01 | Factor02 | Factor03 |
| ---------- | ----- | -------- | -------- |
| Var01      | 1     | 0        | 0        |
| Var02      | 0     | 0.75     | 0        |
| Var03      | 0     | 0.25     | 0        |
| Var04      | 0     | 0        | 5        |
| Var05      | 0     | 0        | 25       |
| Var06      | 0     | 0        | 25       |


#### 2. Factor Table

CSV file with 3 columns listing the loading for each variable/factor combination. 
Header is fixed and must contain the column names `Variable`, `Factor`, `Loading`.

##### Example
<!-- **Example:** -->

Factor Table representation of [factor model example](#factor-model-example) above.

| *Variable* | *Factor* | *Loading* |
| ---------- | -------- | --------- |
| Var01      | Var01    | 1         |
| Var02      | Factor02 | 0.75      |
| Var03      | Factor02 | 0.25      |
| Var04      | Factor03 | 5         |
| Var05      | Factor03 | 25        |
| Var06      | Factor03 | 25        |


#### 3. Factor Variable List

CSV file with 2 columns listing for each factor the associated linear combinations of variables; the coefficients are the inverse of the factor loadings.
Header is fixed and must contain the column names `Factor` and `Variable_list` (or `Factor` and `Variable_set`).

<!-- **Example:** -->
##### Example

Factor Variable List representation of [factor model example](#factor-model-example) above.

| *Factor* | *Variable_list*                           |
| -------- | ----------------------------------------- |
| Var01    | Var01                                     |
| Factor02 | 1.33 * Var02 + 4 * Var03                  |
| Factor03 | 0.2 * Var04 + 0.04 * Var05 + 0.04 * Var06 |

### Output Files

#### Path Diagram

The structure model derived by Galaina is available as a path diagram; the width of a path indicates the likelihood of that path. The diagram is available as a PDF file

#### Structure Model Data

The data used for generating the output diagram is available in multiple format:
* bn.strength - This is the format used by the [bn learn](http://www.bnlearn.com/documentation/man/bn.strength-class.html) package. The output path is determined by the `output_path_bn_strength_obj` parameter value in the `config.ini`.
A `bn.strength` object is an R data frame with the following four columns (one row for each path):

    * `from`, `to`: the factors (or single-factor indicators) connected by the path.
    * `strength`: the likelihood of the path.
    * `direction`: the likelihood that the path is directed from `from` towards `to`. A value equal to 0.5 means that no direction is more likely, so both path edges will have no arrows; a value higher than 0.5 indicates that the path is directed from `from` towards `to`.

    Among the additional attributes of this object, there is the `threshold`: a numeric value, the threshold used to determine if a path likelihood (strength) is significantly high, and thus the path is reliable. 
    Only paths with likelihood above this threshold are displayed in the final structure model.

* CSV (2 files) - The model is saved to a table with the data of the dataframe above. A second version of this CSV file contains only the paths with significantly high likelihood (i.e. `strength >= threshold`); e refer to this CSV file as **filtered**.

    ##### Example

    Let us suppose that by running Galaina on a given dataset and the [example factor model](#factor-model-example) we obtain the following structure model:
    *  A path from `Var01` to `Factor02` with path likelihood 0.8 and direction likelihood 0.6.
    *  A directionless path connecting `Factor02` and `Factor03` with path likelihood 0.6.

    Then the filtered CSV with the structure model is as follow:

| from     | to       | strength | direction |
| -------- | -------- | -------- | --------- |
| Var01    | Factor02 | 0.8      | 0.6       |
| Factor02 | Factor03 | 0.6      | 0.5       |
| Factor03 | Factor02 | 0.6      | 0.5       |


## Acknowledgement and Disclaimer 

This project has received funding from the European Union's Seventh Framework Programme for research, technological development and demonstration under grant agreement no 602805.
The research leading to these results has received funding from the European Community's Seventh Framework Programme (FP7/2007-2013) under Grant Agreement no 278948.
This project is co-funded by the Ambient Assisted Living (AAL) Joint programme, by the German BMBF, the French ANR, the Austrian BMVIT.
