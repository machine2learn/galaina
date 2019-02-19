library(utils)
library(utils)
library(tools)  # to get file extension
library(dplyr)  # actually I load it only for filter
library(ggplot2)  # to save in SVG on R installed via homebrew  https://stackoverflow.com/questions/50767445/how-to-fix-failed-to-load-cairo-dll-in-r
library(ggplotify)  # to convert plot to grob objects
library(grid)
library(stringi)  # to check if strings are empty
library(sbgcop)
library(infotheo)
library(BDgraph)
library(pcalg)  # depends on bioconductor
library(bnlearn)
library(ConfigParser)

library(Rgraphviz)  # for toDot

args <- c("/Users/fabio/projects/aggressotype/code/software/software_demo/config_estonian_restricted_demo_relaxed_less_null.ini")
config <- ConfigParser$new(Sys.getenv(), optionxform=identity)
config$read(filepath = args[1])
# path_to_bn_strength_obj = "/Users/fabio/bn_strength_obj.rds"  # TODO put in INI
path_to_bn_strength_obj <- config$get(option = "output_path_bn_strength_obj", fallback = "", section = "output_paths")
path_to_avg_bn_obj <- config$get( option = "output_path_avg_bn_obj", fallback = "", section = "output_paths" )
bn_strength_obj <- readRDS(path_to_bn_strength_obj)
avg_bn_obj <- readRDS(path_to_avg_bn_obj)
graph_nel <- strength.plot(x = avg_bn_obj, strength = bn_strength_obj, threshold = attributes(bn_strength_obj)$threshold)

# # #
tmp_key_v <- names(graph_nel@renderInfo@edges[["lwd"]])

lwd_key2data_key <- as.list(setNames(object = gsub("[~]", "|", tmp_key_v), nm = tmp_key_v ))


for (iName in names(lwd_key2data_key)) {
  graph_nel@edgeData@data[[lwd_key2data_key[[iName]]]][["penwidth"]] <- graph_nel@renderInfo@edges[["lwd"]][[iName]]
  cat(c(graph_nel@edgeData@data[[lwd_key2data_key[[iName]]]][["penwidth"]], "\n"))
}
toDot(graph = graph_nel, filename = "~/Documents/figures/estonian_experiment_penwidth.gv")

# prova <- as.bn(graph_nel)
# colnames(prova[['arcs']]) <- c('from', 'to')