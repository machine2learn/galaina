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

# dumb_script_for_computing_edge
# args <- c("/Users/fabio/projects/aggressotype/code/software/software_demo/config_estonian_restricted_demo_relaxed_less_null.ini")
args <- c("/Users/fabio/projects/aggressotype/code/galaina/user_data/test/rpq/config_1/config.ini")

config <- ConfigParser$new(Sys.getenv(), optionxform=identity)
config$read(filepath = args[1])

path_pc_algo_obj <- config$get(option = "output_path_pc_algo_obj", section = "output_paths" )
run2pcalgo_ls <- readRDS(path_pc_algo_obj)

path_to_bn_strength_obj <- config$get(option = "output_path_bn_strength_obj", fallback = "", section = "output_paths")
bn_strength_obj <- readRDS(path_to_bn_strength_obj)

path_to_avg_bn_obj <- config$get(option = "output_path_avg_bn_obj", fallback = "", section = "output_paths")
avg_bn_obj <- readRDS(path_to_avg_bn_obj)

# graph_nel <- strength.plot(x = avg_bn_obj, strength = bn_strength_obj, 
#                            # threshold = attributes(bn_strength_obj)$threshold, main = first_title
#                            )


causal_discovery_algorithm_run_n <- length(run2pcalgo_ls)
# run2directed_amat_ls <- vector(mode = "list", length = causal_discovery_algorithm_run_n)
# run2undirected_amat_ls <- vector(mode = "list", length = causal_discovery_algorithm_run_n)

# run2undirected_amat_ls <- run2directed_amat_ls
# run2single_arrow_amat <- run2directed_amat_ls
# run2alt_single_arrow_amat <- run2directed_amat_ls

# Maybe could encode it in all one matrix instead of using 6 binary matrices. But it makes it more readable
run2undirected_bool_amat <- vector(mode = "list", length = causal_discovery_algorithm_run_n)
run2double_arrow_bool_amat <- vector(mode = "list", length = causal_discovery_algorithm_run_n)
run2double_circle_bool_amat <- vector(mode = "list", length = causal_discovery_algorithm_run_n)
#
run2single_arrow_bool_amat <- vector(mode = "list", length = causal_discovery_algorithm_run_n)
run2alt_single_arrow_bool_amat <- vector(mode = "list", length = causal_discovery_algorithm_run_n)
run2single_circle_bool_amat <- vector(mode = "list", length = causal_discovery_algorithm_run_n)

for (iRun_n in 1:causal_discovery_algorithm_run_n) {
  tmp_amat <- as(run2pcalgo_ls[[iRun_n]], "amat")
  unusual_values_v <- !(tmp_amat %in% c(0,1))
  if (any(unusual_values_v)) {
    for (jUnique_n in unique(tmp_amat[unusual_values_v])) {
      # cat(iRun_n, 'has value', jUnique_n, 'not symmetrically for', sum(tmp_amat[unusual_values_v] == jUnique_n), 'times', '\n')
      # cat(iRun_n, 'has value', jUnique_n, 'not symmetrically for', sum(tmp_amat[unusual_values_v] == jUnique_n), 'times', '\n')
      # if (any((tmp_amat != t(tmp_amat)) & tmp_amat == jUnique_n)) {
        cat(iRun_n, 'has value', jUnique_n, 'non symmetrically for a total of', sum((tmp_amat == jUnique_n) & (tmp_amat != t(tmp_amat))), 'values', '\n')
      # }
      # if (any((tmp_amat == jUnique_n) & (tmp_amat != t(tmp_amat)))) {
        cat(iRun_n, 'has value', jUnique_n, 'symmetrically for a total of', sum((tmp_amat == jUnique_n) & (tmp_amat == t(tmp_amat))), 'values', '\n')
      # }
    }
    # TODO
    # Convert the cells with 2 and 3
    # [arrow] a,b = 2 ; b,a = 3 => a,b=0 ;  a,b = 0 ; b,a = 1
    # [circle] a,b = 1 ; b,a = 3 => a,b=0 ;  a,b = 0 ; b,a = 0.5
    # [double arrow] 2 ; 2 -> 1 ; 1
    
    # Wait... why therare 2?
    # tmp_amat <- 0 + (tmp_amat != 0)  
    
    # Reduce to binary, not sure it is the right step... for 2,2 edges it makes sense, but if there is a 3 no
  }
  
  no_or_undirected_edges_in_amat_mask <- (tmp_amat == t(tmp_amat))  # this captures when transposed elemetns are both 1 or both 0
  #
  undirected_mask <- no_or_undirected_edges_in_amat_mask & (tmp_amat == 1) # & (t(tmp_amat) == 1) 
  double_arrow_mask <- no_or_undirected_edges_in_amat_mask & (tmp_amat == 2) # & (t(tmp_amat) == 2) 
  # This should not exist as PC algorithm output, but I put it for completeness
  double_circle_mask <- no_or_undirected_edges_in_amat_mask & (tmp_amat == 3) # & (t(tmp_amat) == 3) 
  #
  single_arrow_mask <- !(no_or_undirected_edges_in_amat_mask) & (tmp_amat == 1) & (t(tmp_amat) == 0) 
  alt_single_arrow_mask <- !(no_or_undirected_edges_in_amat_mask) & (tmp_amat == 3) & (t(tmp_amat) == 2) 
  single_circle_mask <- !(no_or_undirected_edges_in_amat_mask) & (tmp_amat == 3) & (t(tmp_amat) == 1) 
  #
  zero_amat <- tmp_amat
  zero_amat[] <- 0L
  run2undirected_bool_amat[[iRun_n]] <- zero_amat
  run2undirected_bool_amat[[iRun_n]][undirected_mask] <- 1
  
  run2double_arrow_bool_amat[[iRun_n]] <- zero_amat
  run2double_arrow_bool_amat[[iRun_n]][double_arrow_mask] <- 1
  
  run2double_circle_bool_amat[[iRun_n]] <- zero_amat
  run2double_circle_bool_amat[[iRun_n]][double_circle_mask] <- 1
  
  run2undirected_bool_amat[[iRun_n]] <- zero_amat
  run2undirected_bool_amat[[iRun_n]][undirected_mask] <- 1
  
  #
  run2single_arrow_bool_amat[[iRun_n]] <- zero_amat
  run2single_arrow_bool_amat[[iRun_n]][single_arrow_mask] <- 1
  run2alt_single_arrow_bool_amat[[iRun_n]] <- zero_amat
  run2alt_single_arrow_bool_amat[[iRun_n]][alt_single_arrow_mask] <- 1
  run2single_circle_bool_amat[[iRun_n]] <- zero_amat
  run2single_circle_bool_amat[[iRun_n]][single_circle_mask] <- 1
  
  # I guess it also captures i,j =0 and j,i = 0
  # binarized_version_tmp_amat <- 0 + (tmp_amat != 0)
  # directed_amat <- binarized_version_tmp_amat
  # undirected_amat <- binarized_version_tmp_amat
  # 
  # no_or_undirected_edges_in_amat_mask <- (tmp_amat == t(tmp_amat))  # this captures when transposed elemetns are both 1 or both 0
  
  # directed_amat[no_or_undirected_edges_in_amat_mask] <- 0
  # undirected_amat[!(no_or_undirected_edges_in_amat_mask)] <- 0
  
  # if (any(unusual_values_v)) {
  #   # Those 3 masks are symmetric, so if used for counting edges (i.e. for determining edge likelihood) might have to count half in some cases
  #   undirected_mask <- no_or_undirected_edges_in_amat_mask & (tmp_amat == 1) # & (t(tmp_amat) == 1) 
  #   double_arrow_mask <- no_or_undirected_edges_in_amat_mask & (tmp_amat == 2) # & (t(tmp_amat) == 2) 
  #   # Tdis should not exist as PC algorithm output, but I put it for completeness
  #   double_circle_mask <- no_or_undirected_edges_in_amat_mask & (tmp_amat == 3) # & (t(tmp_amat) == 3) 
  #   #
  #   single_arrow_mask <- !(no_or_undirected_edges_in_amat_mask) & (tmp_amat == 1) & (t(tmp_amat) == 0) 
  #   alt_single_arrow_mask <- !(no_or_undirected_edges_in_amat_mask) & (tmp_amat == 3) & (t(tmp_amat) == 2) 
  #   single_circle_mask <- !(no_or_undirected_edges_in_amat_mask) & (tmp_amat == 3) & (t(tmp_amat) == 1) 
  #   
  #   # Wait... better do it later
  #   directed_amat[double_arrow_mask] <- 1  # both edges are directed
  #   undirected_amat[undirected_mask] <- 1  # should be already there
  #   directed_amat[single_arrow_mask] <- 1  # should be already there
  #   directed_amat[alt_single_arrow_mask] <- 1  # PAG amat notation
  #   directed_amat[single_circle_mask] <- 0.75  # but circle counts half for direction
  #   
  #   
  # }
  # #
  
  # run2directed_amat_ls[[iRun_n]] <- directed_amat
  # run2undirected_amat_ls[[iRun_n]] <- undirected_amat
}

inv_causal_discovery_algorithm_run_n <- causal_discovery_algorithm_run_n ^ -1

avg_single_circle_amat <- inv_causal_discovery_algorithm_run_n * Reduce('+', run2single_circle_bool_amat)
avg_single_arrow_amat <- inv_causal_discovery_algorithm_run_n * Reduce('+', run2single_arrow_bool_amat)
avg_alt_single_arrow_amat <- inv_causal_discovery_algorithm_run_n * Reduce('+', run2alt_single_arrow_bool_amat)
avg_undirected_amat <- inv_causal_discovery_algorithm_run_n * Reduce('+', run2undirected_bool_amat)
avg_double_circle_amat <- inv_causal_discovery_algorithm_run_n * Reduce('+', run2double_circle_bool_amat)
avg_double_arrow_amat <- inv_causal_discovery_algorithm_run_n * Reduce('+', run2double_arrow_bool_amat)

avg_single_arrow_amat <- avg_alt_single_arrow_amat + avg_single_arrow_amat

# Used just to see if it is ok
asym_score_mat <- avg_single_arrow_amat + 0.75 * avg_single_circle_amat
numerator_dir_score_mat <- asym_score_mat + 0.5 * avg_undirected_amat + avg_double_arrow_amat
denominator_dir_score_mat <-  asym_score_mat + t(asym_score_mat) + avg_undirected_amat + avg_double_arrow_amat
dir_score_mat <- numerator_dir_score_mat
dir_score_mat[numerator_dir_score_mat != 0] <- dir_score_mat[numerator_dir_score_mat != 0] / denominator_dir_score_mat[numerator_dir_score_mat != 0]


numerator_mark_score_mat <- avg_single_arrow_amat + 0.5 * avg_single_circle_amat + 0.5 * avg_double_arrow_amat + 0.25 * avg_double_circle_amat
denominator_mark_score_mat <- numerator_mark_score_mat + avg_undirected_amat
mark_score_mat <- numerator_mark_score_mat
mark_score_mat[numerator_mark_score_mat != 0] <- mark_score_mat[numerator_mark_score_mat != 0] / denominator_mark_score_mat[numerator_mark_score_mat != 0]

# Implemet my arc strength. But don't use it now because otherwise maybe we'll have to find a way to create a a new bn_strength_obj with:
arc_strength_mat <- 0.5 * (avg_undirected_amat + avg_double_circle_amat + avg_double_arrow_amat) + avg_single_arrow_amat + avg_single_circle_amat
arc_strength_mat <- arc_strength_mat + t(arc_strength_mat)

# amat2arcs

# - this new strength
# - newly def direction
# And new strength means that the threshold would have to be computed again

# # Take average among all networks
# UnArc_mat <- inv_causal_discovery_algorithm_run_n * Reduce('+', run2undirected_amat_ls)
# DirArc_mat <- inv_causal_discovery_algorithm_run_n * Reduce('+', run2directed_amat_ls)
# 
# # compute Mark (arrow, circle scores)
# Mark_score_mat <- UnArc_mat
# Mark_score_mat[] <- 0L
# nonzero_position_mask <- (UnArc_mat != 0)
# Mark_score_mat[nonzero_position_mask] <- UnArc_mat[nonzero_position_mask] / (UnArc_mat[nonzero_position_mask] + DirArc_mat[nonzero_position_mask])


# Mark_score_mat[UnArc_mat != 0] <- UnArc_mat[UnArc_mat != 0] /(UnArc_mat[UnArc_mat != 0] + DirArc_mat[UnArc_mat != 0])


# path_to_bn_strength_obj <- config$get(option = "output_path_bn_strength_obj", fallback = "", section = "output_paths")
# bn_strength_obj <- readRDS(path_to_bn_strength_obj)
