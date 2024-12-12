input_path <- system.file("extdata/capeverde.zip", package = "tidygenclust")
clump_res <- input_path %>% gt_clumppling()
get_modes <- tidy(clump_res,matrix="modes")
get_major_modes <- tidy(clump_res,matrix="major_modes")
get_q_modes <- tidy(clump_res,matrix="q_modes")
get_q_major_modes <- tidy(clump_res,matrix="q_major_modes")

# Plotting the results
autoplot(clump_res, type = "modes_within_k", k=3)
autoplot(clump_res, type = "major_modes")
autoplot(clump_res, type = "all_modes")
autoplot(clump_res, type = "modes")

