## This script is to sample the values of abundance N, recruitment R, Survival S, loss L, growth rate g, recruitment rate r, and loss rate l from the MCMC chains
## in order to propagate the uncertainties
library(dplyr)

## First, let's filter the models for which all Rhats < 1.1
## Analyse summaries
file.list <-list.files(path = "summary_models/", full.names = T)
## Get the species names
species <- gsub("^summary_|\\.rds$","",list.files(path = "summary_models/"))

for(i in 1:length(file.list)){
  assign(species[i],
         readRDS(file.list[i]))
  # print(i)
}

## Filter species that have all parameters Rhat < 1.1
rhats <- data.frame(species = as.character(), param = as.character(), Rhat = as.numeric())
for(sp in species){
  d <- get(sp)
  rhats <- rbind(rhats,
                 data.frame(sp, rownames(d[1:20,]), d[1:20,"Rhat"]))
  # print(sp)
}
colnames(rhats) <- c("species", "param", "rhat")
rhats$rhat <- as.numeric(rhats$rhat)

### Take off the species with at least 1 Rhat > 1.1
ok_sp <-
  rhats %>%
  group_by(species) %>%
  filter(all(rhat <= 1.1, na.rm = T)) %>%
  distinct(species) %>% pull()

## Cleaning
rm(list = ls()[!(ls() %in% c("ok_sp"))])
gc()

##############################################
## List files with MCMC chains
files <- list.files("model_output/", full.names = T)

## Filter the files that have good Rhats
files <- files[gsub("^model_output/|\\.rds$","",files) %in% ok_sp]

## Number of samples
s <- 500
## Create the datasets that will have s samples of the overall N, R, and L at each route
N_overall <- array(0, dim = c(1033, 35, s))
R_overall <- array(0, dim = c(1033, 34, s))
L_overall <- array(0, dim = c(1033, 34, s))
## For the per-species analysis: create dataset that will have mean and sd of the sampled species abundance at the North American scale
mean_sd_N_sp <- data.frame(species = as.character(), year = as.numeric(), mean_N = as.numeric(), sd_N = as.numeric())
mean_sd_R_sp <- data.frame(species = as.character(), year = as.numeric(), mean_R = as.numeric(), sd_R = as.numeric())
mean_sd_L_sp <- data.frame(species = as.character(), year = as.numeric(), mean_L = as.numeric(), sd_L = as.numeric())

mean_sd_g_sp <- data.frame(species = as.character(), year = as.numeric(), mean_g = as.numeric(), sd_g = as.numeric())
mean_sd_r_sp <- data.frame(species = as.character(), year = as.numeric(), mean_r = as.numeric(), sd_r = as.numeric())
mean_sd_l_sp <- data.frame(species = as.character(), year = as.numeric(), mean_l = as.numeric(), sd_l = as.numeric())
mean_sd_absg_sp <- data.frame(species = as.character(), year = as.numeric(), mean_absg = as.numeric(), sd_absg = as.numeric())

## Iterate over each species
for(file in files){

  ## Keep track of time
  start <- Sys.time()

  ## Load MCMC chains
  sp <- readRDS(file)

  ## Number of samples in the original MCMC chains
  it = dim(sp[[1]])[1]

  ## First index is the chain (1:3), e.g. sp[[1]]
  ## Second index: sp[[1]][chain_iteration,parameter*route*year]

  # Create the arrays that will take the dataset with dimensions 1033 x 35 x 2500 x 3 being route x year x iteration x chain
  N_final <- array(NA, dim = c(1033, 35, it, 3))
  R_final <- array(NA, dim = c(1033, 34, it, 3))
  S_final <- array(NA, dim = c(1033, 34, it, 3))
  L_final <- array(NA, dim = c(1033, 34, it, 3))

  ## Loop over each chain
  for(chain in 1:3){

    ## Extract the abundance N, recruitment R and Survival S.
    ## The rows are the number of iterations saved in the MCMC algo (2500), then the columns are the values for each route and each year.
    ## For the columns, you have first the values for a given year for all the routes (1033), then the values for the year after etc
    ## Note that there are 35 year for N, and 34 for R and S
    N <- sp[[chain]][,dimnames(sp[[chain]])[[2]][grepl("^N", dimnames(sp[[chain]])[[2]])]]
    R <- sp[[chain]][,dimnames(sp[[chain]])[[2]][grepl("^R", dimnames(sp[[chain]])[[2]])]]
    S <- sp[[chain]][,dimnames(sp[[chain]])[[2]][grepl("^S", dimnames(sp[[chain]])[[2]])]]

    ## Reshaping N, R, S into [route,year,sample]
    N <- array(t(N), dim = c(1033, 35, it))
    R <- array(t(R), dim = c(1033, 34, it))
    S <- array(t(S), dim = c(1033, 34, it))

    ## Insert in the final dataset
    N_final[,,,chain] <- N
    R_final[,,,chain] <- R
    S_final[,,,chain] <- S
    L_final <- N_final[,-35,,] - S_final

    rm(N,R,S)

  }
  rm(S_final, sp)

  # Pick the chain randomly
  chain_s <- sample(1:3, s, replace = T)
  ## Pick the iteration randomly
  it_s <- sample(1:it, s)

  ## For N
  N_sampled <- array(NA, dim = c(1033, 35, s))
  for(route in 1:1033){
    for(year in 1:35){
      for(sample in 1:s){

        N_sampled[route, year, sample] <- N_final[route, year, it_s[sample], chain_s[sample]]

      }
    }
  }
  rm(N_final)

  ## For R
  R_sampled <- array(NA, dim = c(1033, 34, s))
  for(route in 1:1033){
    for(year in 1:34){
      for(sample in 1:s){

        R_sampled[route, year, sample] <- R_final[route, year, it_s[sample], chain_s[sample]]

      }
    }
  }
  rm(R_final)

  ## For L
  L_sampled <- array(NA, dim = c(1033, 34, s))
  for(route in 1:1033){
    for(year in 1:34){
      for(sample in 1:s){

        L_sampled[route, year, sample] <- L_final[route, year, it_s[sample], chain_s[sample]]

      }
    }
  }
  rm(L_final)

  ## Add N, R, L to get the overall N, R, L at each route
  N_overall <- N_overall + N_sampled
  R_overall <- R_overall + R_sampled
  L_overall <- L_overall + L_sampled

  ## Per species analysis: sum N,R,L over each route for each sample in order to have a distribution of overall abundance at the national level for this species
  N_sp <- t(colSums(N_sampled, dims = 1)) # row = sample / column = year / value = abundance summed over all rows for each year
  R_sp <- t(colSums(R_sampled, dims = 1)) # row = sample / column = year / value = recruitment summed over all rows for each year
  L_sp <- t(colSums(L_sampled, dims = 1)) # row = sample / column = year / value = loss summed over all rows for each year

  ## Compute absolute growth rate, per capita growth rate, recruitment and loss rates for the species at the scale of North America
  g_sp <- (N_sp[,-1] - N_sp[,-35])/N_sp[,-35] # row = sample / column = year
  r_sp <- R_sp/N_sp[,-35]                     # row = sample / column = year
  l_sp <- L_sp/N_sp[,-35]                     # row = sample / column = year
  absg_sp <- N_sp[,-1] - N_sp[,-35]

  ## Per species analysis: compute mean and sd over the s samples
  mean_sd_N_sp <-
  rbind(mean_sd_N_sp,
        data.frame(species = gsub("^model_output/|\\.rds$","",file),
                   year = 1:35,
                   mean_N = apply(N_sp, 2, mean),
                   sd_N = apply(N_sp, 2, sd)))

  mean_sd_R_sp <-
    rbind(mean_sd_R_sp,
          data.frame(species = gsub("^model_output/|\\.rds$","",file),
                     year = 1:34,
                     mean_R = apply(R_sp, 2, mean),
                     sd_R = apply(R_sp, 2, sd)))

  mean_sd_L_sp <-
    rbind(mean_sd_L_sp,
          data.frame(species = gsub("^model_output/|\\.rds$","",file),
                     year = 1:34,
                     mean_L = apply(L_sp, 2, mean),
                     sd_L = apply(L_sp, 2, sd)))

  mean_sd_g_sp <-
    rbind(mean_sd_g_sp,
          data.frame(species = gsub("^model_output/|\\.rds$","",file),
                     year = 1:34,
                     mean_g = apply(g_sp, 2, mean),
                     sd_g = apply(g_sp, 2, sd)))

  mean_sd_r_sp <-
    rbind(mean_sd_r_sp,
          data.frame(species = gsub("^model_output/|\\.rds$","",file),
                     year = 1:34,
                     mean_R = apply(r_sp, 2, mean),
                     sd_R = apply(r_sp, 2, sd)))

  mean_sd_l_sp <-
    rbind(mean_sd_l_sp,
          data.frame(species = gsub("^model_output/|\\.rds$","",file),
                     year = 1:34,
                     mean_l = apply(l_sp, 2, mean),
                     sd_l = apply(l_sp, 2, sd)))

  mean_sd_absg_sp <-
    rbind(mean_sd_absg_sp,
          data.frame(species = gsub("^model_output/|\\.rds$","",file),
                     year = 1:34,
                     mean_absg = apply(absg_sp, 2, mean),
                     sd_absg = apply(absg_sp, 2, sd)))


  rm(N_sampled, R_sampled, L_sampled)
  gc()

  print(paste0("Iteration ", file," took ", Sys.time() - start))

}

## Save the output of this loop
saveRDS(N_overall, "mixed_model/save_samples/N_overall.rds")
# saveRDS(R_overall, "mixed_model/save_samples/R_overall.rds")
# saveRDS(L_overall, "mixed_model/save_samples/L_overall.rds")
saveRDS(mean_sd_N_sp, "mixed_model/save_samples/mean_sd_N_sp.rds")
# saveRDS(mean_sd_R_sp, "mixed_model/save_samples/mean_sd_R_sp.rds")
# saveRDS(mean_sd_L_sp, "mixed_model/save_samples/mean_sd_L_sp.rds")
saveRDS(mean_sd_g_sp, "mixed_model/save_samples/mean_sd_g_sp.rds")
# saveRDS(mean_sd_r_sp, "mixed_model/save_samples/mean_sd_rr_sp.rds")
# saveRDS(mean_sd_l_sp, "mixed_model/save_samples/mean_sd_lr_sp.rds")
saveRDS(mean_sd_absg_sp, "mixed_model/save_samples/mean_sd_absg_sp.rds")

## Compute rates g,r,l of the overall abundance at each route
# First load the samples
N_overall <- readRDS("mixed_model/save_samples/N_overall.rds")
# R_overall <- readRDS("mixed_model/save_samples/R_overall.rds")
# L_overall <- readRDS("mixed_model/save_samples/L_overall.rds")

g_overall <- (N_overall[,-1,] - N_overall[,-35,]) / N_overall[,-35,]
absg_overall <- N_overall[,-1,] - N_overall[,-35,]
# r_overall <- R_overall / N_overall[,-35,]
# l_overall <- L_overall / N_overall[,-35,]

## Compute growth rate relative to time 1
# gt0_overall <- array(NA, dim = c(1033, 35, 500))
# for(k in 1:dim(gt0_overall)[3]){
#   gt0_overall[,,k] <- (N_overall[,,k] - N_overall[,1,k]) / N_overall[,1,k]
# }


## Compute the mean and sd of N,R,L,g,r,l, for each year and route over the s samples
mean_N_overall <- apply(N_overall, c(1,2), mean)
sd_N_overall <- apply(N_overall, c(1,2), sd)
# mean_R_overall <- apply(R_overall, c(1,2), mean)
# sd_R_overall <- apply(R_overall, c(1,2), sd)
# mean_L_overall <- apply(L_overall, c(1,2), mean)
# sd_L_overall <- apply(L_overall, c(1,2), sd)

mean_g_overall <- apply(g_overall, c(1,2), mean)
sd_g_overall <- apply(g_overall, c(1,2), sd)
# mean_r_overall <- apply(r_overall, c(1,2), mean)
# sd_r_overall <- apply(r_overall, c(1,2), sd)
# mean_l_overall <- apply(l_overall, c(1,2), mean)
# sd_l_overall <- apply(l_overall, c(1,2), sd)

mean_absg_overall <- apply(absg_overall, c(1,2), mean)
sd_absg_overall <- apply(absg_overall, c(1,2), sd)

# mean_gt0_overall <- apply(gt0_overall, c(1,2), mean)
# sd_gt0_overall <- apply(gt0_overall, c(1,2), sd)

saveRDS(mean_N_overall, "mixed_model/save_samples/mean_N_overall.rds")
saveRDS(sd_N_overall, "mixed_model/save_samples/sd_N_overall.rds")
# saveRDS(mean_R_overall, "mixed_model/save_samples/mean_R_overall.rds")
# saveRDS(sd_R_overall, "mixed_model/save_samples/sd_R_overall.rds")
# saveRDS(mean_L_overall, "mixed_model/save_samples/mean_L_overall.rds")
# saveRDS(sd_L_overall, "mixed_model/save_samples/sd_L_overall.rds")

saveRDS(mean_g_overall, "mixed_model/save_samples/mean_g_overall.rds")
saveRDS(sd_g_overall, "mixed_model/save_samples/sd_g_overall.rds")
# saveRDS(mean_r_overall, "mixed_model/save_samples/mean_rr_overall.rds")
# saveRDS(sd_r_overall, "mixed_model/save_samples/sd_rr_overall.rds")
# saveRDS(mean_l_overall, "mixed_model/save_samples/mean_lr_overall.rds")
# saveRDS(sd_l_overall, "mixed_model/save_samples/sd_lr_overall.rds")

# saveRDS(mean_gt0_overall, "mixed_model/save_samples/mean_gt0_overall.rds")
# saveRDS(sd_gt0_overall, "mixed_model/save_samples/sd_gt0_overall.rds")

saveRDS(mean_absg_overall, "mixed_model/save_samples/mean_absg_overall.rds")
saveRDS(sd_absg_overall, "mixed_model/save_samples/sd_absg_overall.rds")

## Run check over time series
# route = 1000
# sample = 30
# chain = 2
# all.equal(N_final[route,,sample,chain],
#           unname(sp[[chain]][,dimnames(sp[[chain]])[[2]][grepl(paste0("^N\\[",route,","), dimnames(sp[[chain]])[[2]])]][sample,]))

