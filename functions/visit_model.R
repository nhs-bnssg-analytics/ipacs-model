# Visit-based model (for P1 pathways) -----------------------------------------

# Set up ----------------------------------------------------------------------

# Create objects for simulation using setup_all(), returning list of 2 where
# 1 contains object names and 2 contains the objects. Extract from list into
# workspace using assign()
setup_visit <- setup_all("visit")
for (i in seq_along(setup_visit[[1]])){
  assign(setup_visit[[1]][i], setup_visit[[2]][[i]])
}

# Initial service rate and end service rate, and their standard deviation
# Create lists with sd_isr or sd_esr repeated for number of visit scenarios
isr <- as.integer(scenarios_visit$IVR)
end_sr <- as.integer(scenarios_visit$FVR)
sd_isr <- as.double(rep(sd_isr, nrow(scenarios_visit)))
sd_esr <- as.double(rep(sd_esr, nrow(scenarios_visit)))

# Create n_slots, the number of visit slots available per day. This is based
# on an average visit rate (as from mean of isr and end_sr) multiplied by the
# capacity for P1 (cap_visit)
n_slots  <- cap_visit * mean(c(isr, end_sr))

# Create object to store outputs from each scenario
visits_based_output <- NULL
visits_based_output_q <- NULL


# Model -----------------------------------------------------------------------

# Repeat for each scenario in visit-based pathways
for (z in seq_along(pathway_vector_visit)) {
  # detectCores() but -1 as want to make you you have one left to do other
  # stuff on. Cores are processors that can work on tasks. Then makecluster()
  # to set the amount of clusters you want your code to run on
  cl <- parallel::makeCluster(detectCores() - 1)
  registerDoSNOW(cl)
  run <- 1

  # Use foreach() to repeat operation for each run
  results <- foreach(run = 1:nruns, .combine = "rbind") %dopar% {
    set.seed(nruns * (temp_seed - 1) + run)

    # Output variables
    ent_sys <- 0 # number of entities that entered the system
    left_sys <- 0 # number of entities that left the system

    # Create necessary data structures
    # Captures output after warmup
    output <- create_output_df(nrow = sim_length)
    patients <- create_patient_df(nrow = (sim_length + warmup) * 10)

    # Stores wait times for patients who leave system
    waittime_vec <- create_wait_df()

    # List with required visit vectors for each patient
    req_visits <- list()

    # Create resources, *10 to make it sufficiently large
    resources <- matrix(data = n_slots[z], nrow = (sim_length + warmup) * 10)

    # Initialising counter for patients dataframe, plus ID and t (day)
    npat <- 0
    id <- 0
    t <- 1

    ### Initial conditions (already in P1) ------------------------------------
    # For each patient...
    for (j in 1:init_occ_visit[[z]]) {
      add_patient_output <- add_patient(in_system = TRUE)
      id <- add_patient_output[[1]]
      npat <- add_patient_output[[2]]
      req_visits <- add_patient_output[[3]]
      patients <- add_patient_output[[4]]
      resources <- add_patient_output[[5]]
    }

    ### Initial conditions (waiting to go to P1) ------------------------------
    # For each patient...
    for (j in 1:init_niq_visit[[z]]) {
      add_patient_output <- add_patient(in_system = FALSE)
      id <- add_patient_output[[1]]
      npat <- add_patient_output[[2]]
      req_visits <- add_patient_output[[3]]
      patients <- add_patient_output[[4]]
      resources <- add_patient_output[[5]]
    }

    # Set ent_sys as npat
    ent_sys <- npat

    ### Simulation ------------------------------------------------------------
    for (t in 1:(sim_length + warmup)) {
      # Sample from poisson distribution to get number of arrivals
      # t is the day, and z+1 is the appropriate pathway/location/scenario
      narr <- round(rpois(n = 1,
                          lambda = as.numeric(arr_rates_visit[t, z + 1])))

      # If there are arrivals...
      if (narr > 0) {
        ent_sys <- ent_sys + narr

        # For each arrived patient...
        for (j in 1:narr) {
          add_patient_output <- add_patient(in_system = FALSE)
          id <- add_patient_output[[1]]
          npat <- add_patient_output[[2]]
          req_visits <- add_patient_output[[3]]
          patients <- add_patient_output[[4]]
          resources <- add_patient_output[[5]]
        }
      }

      # Find patients in queue, increment wait time column by one day
      in_q <- which((patients$start_service > t) & (patients$id > 0))
      if (length(in_q) > 0) {
        patients[in_q, "wait_time"] <- patients[in_q, "wait_time"] + 1
      }

      # Recording output from the day warm up period has finished
      if (t > warmup) {
        output[t - warmup, ] <- c(
          RUNX = run,
          node = pathway_vector_visit[z],
          day = t,
          q_length = length(in_q),
          n_slots_used = n_slots[z] - (resources[t, ]),
          patients_in_service = (n_slots[z] - (resources[t, ])) /
            (mean(c(isr[z], end_sr[z]))),
          res_used = 1 - (resources[t, ] / n_slots[z]),
          res_idle = resources[t, ] / n_slots[z],
          in_sys = (ent_sys - left_sys))

        # Remove patients whose service has ended from the patients table
        remove <- which(patients$end_service == t)
        if (length(remove) > 0) {
          if (t >= warmup) {
            df <- data.frame(
              RUNX = run,
              day_ = t,
              scen_ = pathway_vector_visit[z],
              start_service = patients$start_service[remove],
              waittime = patients[remove, 6])
            waittime_vec <- rbind(waittime_vec, df) #keeping waiting time
          }
          patients <- patients[-remove, ] #remove from patient list
          npat <- npat - length(remove)
          left_sys <- left_sys + length(remove)
        }
      }
    }
    list <- list(output, resources, waittime_vec)

    return(list)
  }
  stopCluster(cl)

  # Extract results from above (contains results from each run)...

  # results[,1] contains "output"
  out <- do.call(rbind, results[, 1]) %>%
    mutate_at(c("n_slots_used", "patients_in_service", "res_used",
                "res_idle", "in_sys"), as.numeric) %>%
    mutate_at(c("RUNX", "day", "q_length"), as.integer)

  # results[,2] contains "resources"
  res <- do.call(cbind, results[, 2])
  colnames(res) <- c(1:nruns)

  # results[,3] contains "waittimes"
  wait <- do.call(rbind, results[, 3])

  # Create dataframe for summary information from each run
  summary <- create_summary_df(nruns)
  summary$LOS <- 1 / mean_los_visit[[z]]
  summary$ISR <- isr[z]
  summary$nruns <- nruns
  summary$sim_length <- sim_length
  summary$warm_up <- warmup
  summary$capacity <- n_slots[z]

  for (k in 1:nruns) {
    # Extract results for that run
    r_out <- which(out[, "RUNX"] == k)
    k_wait <- which(wait[, "RUNX"] == k)
    # Add results for that run (with k_wait or r_out extracting the rows)
    summary[k, "mean_wait"] <- round(mean(wait[k_wait, "waittime"]), 2)
    for (x in c("q_length", "res_used", "res_idle", "in_sys")){
      summary[k, x] <- round(mean(out[r_out, x]), 2)
    }
  }

  # Normal outputs ------------------------------------------------------------
  # Groups by day (e.g. day 1) and node (e.g. P1_B_BCap_Blos_Barr)
  # Finds average results for each day (for niq + occ, plus other measures
  # if needed)
  ts_output <- out %>%
    group_by(day, node) %>%
    summarise(
      niq = mean(q_length),
      in_sys = mean(in_sys),
      n_slots_used = mean(n_slots_used),
      occ = mean(patients_in_service),
      mean_res_idle = mean(res_idle),
      mean_res_used = mean(res_used)
    ) %>%
    ungroup()

  # Create node column by extracting first two parts of the scenario
  # (e.g. "P1_LocB", dropping "BCap_Bloc_BArr")
  ts_output$node <- sapply(ts_output$node, function(x) {
    paste(unlist(str_split(x, "_"))[1:2], collapse = "_")})

  # Create costs column
  ts_output <- left_join(ts_output, costs_visit, by = "node") %>%
    mutate(cost = (niq * acute_dtoc) + (n_slots_used * community_cost))

  # Waits by day
  ts_waits <- wait %>%
    group_by(day_, scen_) %>%
    summarise(wait = mean(waittime)) %>%
    ungroup()

  # Combine ts_output and ts_waits, and then save the scenario results
  ts_output <- cbind(ts_output, ts_waits)
  visits_based_output <- rbind(visits_based_output, ts_output)

  # Stochastic outputs --------------------------------------------------------
  # Quantiles for niq, occ and wait for optional stochastic report
  q_niq <- out %>%
    group_by(day, node) %>%
    summarise(
      q05 = quantile(q_length, 0.05), q5 = quantile(q_length, 0.5),
      q95 = quantile(q_length, 0.95), mean = mean(q_length),
      q10 = quantile(q_length, 0.1), q25 = quantile(q_length, 0.25),
      q75 = quantile(q_length, 0.75), q90 = quantile(q_length, 0.9),
      q025 = quantile(q_length, 0.025), q975 = quantile(q_length, 0.975)
    ) %>%
    ungroup()
  q_niq$measure <- rep("niq", nrow(q_niq))

  q_occ <- out %>%
    group_by(day, node) %>%
    summarise(
      q05 = quantile(patients_in_service, 0.05),
      q5 = quantile(patients_in_service, 0.5),
      q95 = quantile(patients_in_service, 0.95),
      mean = mean(patients_in_service),
      q10 = quantile(patients_in_service, 0.1),
      q25 = quantile(patients_in_service, 0.25),
      q75 = quantile(patients_in_service, 0.75),
      q90 = quantile(patients_in_service, 0.9),
      q025 = quantile(patients_in_service, 0.025),
      q975 = quantile(patients_in_service, 0.975)
    ) %>%
    ungroup()
  q_occ$measure <- rep("occ", nrow(q_occ))

  q_waits <- wait %>%
    group_by(day_, scen_) %>%
    summarise(
      q05 = quantile(waittime, 0.05), q5 = quantile(waittime, 0.5),
      q95 = quantile(waittime, 0.95), mean = mean(waittime),
      q10 = quantile(waittime, 0.1), q25 = quantile(waittime, 0.25),
      q75 = quantile(waittime, 0.75), q90 = quantile(waittime, 0.9),
      q025 = quantile(waittime, 0.025), q975 = quantile(waittime, 0.975)
    ) %>%
    ungroup() %>%
    rename_at(1, ~ "day") %>%
    rename_at(2, ~ "node")
  q_waits$measure <- rep("waits", nrow(q_waits))

  # For each scenario for optional stochastic report
  ts_output_q <- rbind(q_niq, q_occ, q_waits)
  ts_output_q <- cbind(
    ts_output_q,
    data.frame(arr_rates_visit$date[seq_along(arr_rates_visit$date)])) %>%
    rename_at(14, ~"date")

  # Rowbind each scenario for optional stochastic report
  visits_based_output_q <- rbind(visits_based_output_q, ts_output_q)
}


# Save model results ----------------------------------------------------------

# Extract and pivot results for NIQ, OCC, wait and costs, then append together,
# along with the dates
outcomes <- c("niq", "occ", "wait", "cost")
ptvisits_list <- list(arr_rates_visit["date"])
for (i in seq_along(outcomes)){
  ptvisits_list[[i + 1]] <- visits_based_output[c("day", "scen_",
                                                  outcomes[i])] %>%
    pivot_wider(names_from = "scen_", values_from = as.name(outcomes[i])) %>%
    select(-day) %>%
    setNames(paste0(names(.), "__", outcomes[i]))
}
means_output_v <- do.call("cbind", ptvisits_list)

# Save output to csv with output filename based on input filename
output_filename <- paste0("outputs/report_data/visit_output_using_",
                          gsub(".xlsx", "", input_filename),
                          ".csv")
write.csv(means_output_v, output_filename, row.names = FALSE)

# Correct date formatting, and save the quantiles to csv for optional
# stochastic report. Saved in long format for plotting
visits_based_output_q$date <- as.Date(visits_based_output_q$date)
stoch_filename <- paste0("outputs/stochastic_data/stochastic_visit_output_using_",
                         gsub(".xlsx", "", input_filename),
                         ".csv")
write.csv(visits_based_output_q, stoch_filename, row.names = FALSE)
