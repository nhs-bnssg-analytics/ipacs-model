# Visit-based model (for P2 + P3 pathways) ------------------------------------

# Set up ----------------------------------------------------------------------
# Create objects for simulation using setup_all(), returning list of 2 where
# 1 contains object names and 2 contains the objects. Extract from list into
# workspace using assign()
setup_bed <- setup_all("bed")
for (i in seq_along(setup_bed[[1]])){
  assign(setup_bed[[1]][i], setup_bed[[2]][[i]])
}

# Record time taken
start_time_bed <- Sys.time()

# Run simulation ---------------------------------------------------------------
# For each pathway/location/scenario, run the simulation (nruns times)
# Saves results as sim_res, a list
sim_res <- lapply(1:(ncol(arr_rates_bed) - 1), function(node) {
  # Set parameters -------------------------------------------------------------
  # For that pathway and simulation, assign parameters (initial occupancy,
  # initial queue, arrival rates, LOS distribution parameters, capacity, loss
  # (balking condition))
  node_init_occ <- as.numeric(init_occ_bed[[node]])
  node_init_niq <- as.numeric(init_niq_bed[[node]])
  node_arr_rates <- arr_rates_bed %>% select(date, pathway_vector_bed[node])
  node_srv_dist <- srv_dist_bed[[node]]
  node_srv_params <- srv_params_bed[[node]]
  node_cap <- cap_bed[[node]]
  node_loss <- loss_bed[[node]]

  # Run using parallel processing ----------------------------------------------
  # Intialisation for parallel processing
  # detectCores() but -1 as want to make you you have one left to do other
  # stuff on, then makecluster() to set the amount of clusters you want your
  # code to run on
  cl <- parallel::makeCluster(detectCores() - 1)

  # Create a cluster with all parameters needed for running simfn
  clusterExport(cl = cl,
                varlist = c("node", "node_init_occ", "node_init_niq",
                            "node_arr_rates", "node_srv_dist",
                            "node_srv_params", "node_cap", "node_loss",
                            "rtdist"),
                envir = environment())
  clusterEvalQ(cl = cl, c(library(tidyr), library(dplyr)))

  # Apply using parallel processing simfn for nruns time using information in cl
  tres <- parLapply(cl, 1:nruns, simfn,
                    node = node,
                    node_init_occ = node_init_occ,
                    node_init_niq = node_init_niq,
                    node_arr_rates = node_arr_rates,
                    node_srv_params = node_srv_params,
                    node_cap = node_cap,
                    node_loss = node_loss)
  stopCluster(cl)

  # Bind together results -----------------------------------------------------
  # Bind together results from each run into one dataframe
  # For tres1 (res) and tres2 (res_arr_neg)
  tres1 <- do.call("bind_rows",
                   lapply(seq_along(tres), function(x) tres[[x]]))

  return(tres1)
})

# Bind together results from each pathway and scenario
res1 <- do.call("bind_rows", lapply(seq_along(sim_res),
                                    function(x) sim_res[[x]]))

# Replace NA mean wait with 0
res1$mean_wait[is.nan(res1$mean_wait)] <- 0

# Print processing time
print(difftime(Sys.time(), start_time_bed), quote = FALSE)


# Find average results for report ----------------------------------------------
# Pivot dataframe so instead of columns for arr_admit, arr_no_admit, mean_wait,
# niq and occ, these are now rows (so from wide to long format)
res1mean <- res1 %>%
  pivot_longer(
    cols = c(occ, niq, arr_admit, arr_no_admit, mean_wait),
    names_to = "measure",
    values_to = "value") %>%
  group_by(node, time, measure) %>%
  summarise(mean = mean(value, na.rm = TRUE))

# Add column to res1mean with capacity, based on node number referencing to cap
cap_size <- as.numeric()
for (x in seq_along(res1mean$node)) {
  cap_size[x] <- cap_bed[[res1mean$node[x]]]
}
res1mean["capacity"] <- cap_size

# Find mean for given measure for each node for each day
beds_required <- find_mean("occ")
niq_result <- find_mean("niq")
wait_result <- find_mean("mean_wait")

# Create cost column ----------------------------------------------------------
# Calculate costs using find_costs()
p2_beds_cost <- find_costs(oc = beds_required[grep("P2", pathway_vector_bed)],
                           node = "P2_B", cost_type = "community_cost")
p3_beds_cost <- find_costs(oc = beds_required[grep("P3", pathway_vector_bed)],
                           node = "P3_B", cost_type = "community_cost")
niq_cost <- find_costs(oc = niq_result, node = "P2_B", cost_type = "acute_dtoc")

# Add together NIQ costs and P2 and P3 bed costs
beds_cost <- mapply(
  "+", niq_cost, c(p2_beds_cost, p3_beds_cost), SIMPLIFY = FALSE)


# Save results in an output dataframe ------------------------------------------
# Create output dataframe with mean result for each measure for each
# pathway/scenario/location, with rows for each day (1-181)
meansoutput <- cbind(
  data.frame(arr_rates_bed$date),
  sapply(c(beds_required, niq_result, wait_result, beds_cost),
         function(x) data.frame(round(data.frame(x))))
)
colnames <- cbind(
  c("date", sapply(c("occ", "niq", "wait", "cost"),
                   function(x) paste0(pathway_vector_bed, "__", x))))
colnames(meansoutput) <- colnames

# Save to excel, with filename based on the input file
write.csv(meansoutput,
          paste0("outputs/report_data/bed_output_using_",
                 gsub(".xlsx", "", input_filename), ".csv"),
          row.names = FALSE)

# Create and save stochastic outputs dataframe --------------------------------
# Quantile outputs for optional statistic report
res1q <- res1 %>%
  pivot_longer(cols = c(occ, niq),
               names_to = "measure",
               values_to = "value") %>%
  group_by(node, time, measure) %>%
  summarise(q05 = quantile(value, 0.05, na.rm = TRUE),
            q5 = quantile(value, 0.5, na.rm = TRUE),
            q95 = quantile(value, 0.95, na.rm = TRUE),
            mean = mean(value, na.rm = TRUE),
            q10 = quantile(value, 0.1, na.rm = TRUE),
            q25 = quantile(value, 0.25, na.rm = TRUE),
            q75 = quantile(value, 0.75, na.rm = TRUE),
            q90 = quantile(value, 0.9, na.rm = TRUE),
            q025 = quantile(value, 0.025, na.rm = TRUE),
            q975 = quantile(value, 0.975, na.rm = TRUE)) %>%
  ungroup()

# Add column to res1q with capacity, based on node number referencing to cap
cap_size_q <- as.numeric()
for (x in seq_along(res1q$node)) {
  cap_size_q[x] <- cap_bed[[res1q$node[x]]]
}
res1q["capacity"] <- cap_size_q

# Loop through nodes and apply transformation with scenarios
pathway_vector_bed <- as.list(pathway_vector_bed)
node_values <- sapply(res1q$node, function(x) pathway_vector_bed[[x]])
res1q$node_name <- node_values

# Save to excel, with filename based on the input file
write.csv(res1q,
          paste0("outputs/stochastic_data/stochastic_bed_output_using_",
                 gsub(".xlsx", "", input_filename), ".csv"),
          row.names = FALSE)
