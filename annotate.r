rm(list = ls())

library(data.table)
library(parallel)
source("annotate_fun.r")

# Your existing data.table
chiResults <- readRDS("data/input/5feb_tem_ezchi_c3.rds")

# path to gene annotation file SQLite database
dbPath <- "data/input/annotations/chrom3_1.db" # Replace with the path to your database

# Setup parallel processing
no_cores <- detectCores() - 2  # Leave one core free for system processes
cl <- makeCluster(no_cores)

# Load RSQLite on each worker
clusterEvalQ(cl, library(RSQLite))

# Export necessary variables to the cluster
clusterExport(cl, varlist = c("get_info", "dbPath"))

# Apply the function in parallel and create a new column
chiResults[, info := parSapply(cl, nucPosition, get_info, dbPath)]

# Stop the cluster
stopCluster(cl)

saveRDS(chiResults, "data/output/annotated_5feb_chrom3.rds")
