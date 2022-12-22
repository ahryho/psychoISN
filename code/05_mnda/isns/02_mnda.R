library(data.table)

library(mnda)
library(keras)
# install_keras()
library(aggregation)

library(tictoc)

# rslt_dir   <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-isns/results/5_fold_cv/all/"

# args        <- commandArgs(T)
# input_dir <- as.character(args[1]) 
# out_dir     <- as.character(args[2]) 

rslt_dir    <- "~/bio/code/kul/dex-stim-human-array-isns/results/data/5_fold_cv/all/mnda/"
input_dir   <- paste0(rslt_dir, "isns/input/")
out_dir     <- paste0(rslt_dir, "isns/results/")

system(paste0("mkdir -p ", out_dir))

isns_fn_lst <- list.files(input_dir, full.names = T)

## Load data

lapply(isns_fn_lst[-c(1:2)], function(input_graph_fn){
  
  sample_id  <- gsub(".*_input_(.+).csv", "\\1", input_graph_fn) 
  out_fn     <- paste0(out_dir, "isn_mnda_distances_", sample_id, ".rds")
  graph_data <- fread(input_graph_fn, sep = ";")
  
  ## Create embedding space
  
  print("Start creating embedding space ...", quote = F)
  tic("Create embedding space")
  
  embedding_space_list <- mnda_embedding_2layer(graph_data, train.rep = 50, 
                                                walk.rep = 50,
                                                n.steps = 10,
                                                embedding.size = 10, 
                                                batch.size = 32,
                                                epochs = 10,
                                                random.walk = T, null.perm = TRUE)
  
  toc()
  
  ## Calculate the node-pair distances and p-values
  
  print("Start calculation of the node-pair distances and p-values ...", quote = F)
  tic("Calculate the node-pair distances and p-values")
  
  mnda_dist <- mnda_node_detection_2layer(embedding_space_list, p.adjust.method = "BH")
  
  toc()
  
  ## Save results
  
  tic("Save results")
  saveRDS(list(embedding_space = embedding_space_list, mnda_dist = mnda_dist), out_fn)
  toc()
})
