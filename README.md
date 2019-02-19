# Galaina

* This folder `galaina` contains the Galaina software with user interface
* To run this software:
```
conda activate galaina_env
python main.py
```
Now the internal web-server is active and it can be used.
The server uses by default the `R` code in this folder, but in theory the config file can be changed to use other files

* I believe that the folder `software` was used for some testing and playing with different configuration
* The R code used is
    * `galaina/R_code/20180725_use_config_ini_final_part.R` :: run all the remaining part of the pipeline, that is in `R`. Inside it uses the file below.
    * `galaina/R_code/my_inferCopulaFactorModel.R` :: core `R` file computing the causal graph


## Other code

* `galaina/R_code/just_for_playing_with_graph.R` :: used to test writing on `dot` file, currently not implemented