## This script is to sample the values of abundance N, recruitment R, Survival S, loss L, growth rate g, recruitment rate r, and loss rate l from the MCMC chains
## in order to propagate the uncertainties at the family level
library(dplyr)
library(tidyr)
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

## Let's check the Rhats
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

## Create the lists of each family ##########
family <- readr::read_csv("data/phylo.csv")
hab <- readRDS("data/selected_species.rds")

family <- family %>%
  right_join(hab %>% rename(English_Common_Name = COM_NAME_updated)) %>%
  select(-AOU) %>%
  distinct() %>%
  drop_na()

## Filter the species with Rhat < 1.1
family = family[family$English_Common_Name %in% ok_sp,]

## Number of samples
s <- 500
## Create the datasets that will have s samples of the overall N, R, and L at each route
N_family <- array(0, dim = c(length(unique(family$Family)), 35, s))
R_family <- array(0, dim = c(length(unique(family$Family)), 34, s))
L_family <- array(0, dim = c(length(unique(family$Family)), 34, s))

## Give name of the first dimension as the family
dimnames(N_family) <- list(unique(family$Family))
dimnames(R_family) <- list(unique(family$Family))
dimnames(L_family) <- list(unique(family$Family))

## For the per-species analysis: create dataset that will have mean and sd of the sampled species abundance at the North American scale
mean_sd_N_fam <- data.frame(family = as.character(), year = as.numeric(), mean_N = as.numeric(), sd_N = as.numeric())
mean_sd_R_fam <- data.frame(family = as.character(), year = as.numeric(), mean_R = as.numeric(), sd_R = as.numeric())
mean_sd_L_fam <- data.frame(family = as.character(), year = as.numeric(), mean_L = as.numeric(), sd_L = as.numeric())

mean_sd_g_fam <- data.frame(family = as.character(), year = as.numeric(), mean_g = as.numeric(), sd_g = as.numeric())
mean_sd_r_fam <- data.frame(family = as.character(), year = as.numeric(), mean_r = as.numeric(), sd_r = as.numeric())
mean_sd_l_fam <- data.frame(family = as.character(), year = as.numeric(), mean_l = as.numeric(), sd_l = as.numeric())
mean_sd_absg_fam <- data.frame(family = as.character(), year = as.numeric(), mean_absg = as.numeric(), sd_absg = as.numeric())

## Iterate over each family
for(fam in unique(family$Family)){

  tmp_fam <- family[family$Family == fam,]

  for(species in tmp_fam$English_Common_Name){

    ## Keep track of time
    start <- Sys.time()

    sp <- readRDS(files[grepl(species, files)])

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

    ## Sum N,R,L over each route for each sample in order to have a distribution of overall abundance at the national level for this species
    N_sp <- t(colSums(N_sampled, dims = 1)) # row = sample / column = year / value = abundance summed over all rows for each year
    R_sp <- t(colSums(R_sampled, dims = 1)) # row = sample / column = year / value = recruitment summed over all rows for each year
    L_sp <- t(colSums(L_sampled, dims = 1)) # row = sample / column = year / value = loss summed over all rows for each year


    N_family[which(fam == unique(family$Family)),,] <- N_family[which(fam == unique(family$Family)),,] + t(N_sp)
    R_family[which(fam == unique(family$Family)),,] <- R_family[which(fam == unique(family$Family)),,] + t(R_sp)
    L_family[which(fam == unique(family$Family)),,] <- L_family[which(fam == unique(family$Family)),,] + t(L_sp)

    print(paste0("Iteration ", species," took ", Sys.time() - start))

  }

}

## Save the output of this loop
saveRDS(N_family, "mixed_model/save_samples/N_family.rds")
saveRDS(R_family, "mixed_model/save_samples/R_family.rds")
saveRDS(L_family, "mixed_model/save_samples/L_family.rds")


## Compute rates g,r,l of families
# First load the samples
N_families <- readRDS("mixed_model/save_samples/N_family.rds")
# R_families <- readRDS("mixed_model/save_samples/R_family.rds")
# L_families <- readRDS("mixed_model/save_samples/L_family.rds")

g_families <- (N_families[,-1,] - N_families[,-35,]) / N_families[,-35,]
# r_families <- R_families / N_families[,-35,]
# l_families <- L_families / N_families[,-35,]

absg_families <- N_families[,-1,] - N_families[,-35,]

# ## Compute growth rate relative to time 1
# gt0_overall <- array(NA, dim = c(1033, 35, 500))
# for(k in 1:dim(gt0_overall)[3]){
#   gt0_overall[,,k] <- (N_overall[,,k] - N_overall[,1,k]) / N_overall[,1,k]
# }

# ## Compute the mean and sd of N,R,L,g,r,l, for each year and route over the s samples
mean_N_families <- apply(N_families, c(1,2), mean)
sd_N_families <- apply(N_families, c(1,2), sd)
# mean_R_overall <- apply(R_overall, c(1,2), mean)
# sd_R_overall <- apply(R_overall, c(1,2), sd)
# mean_L_overall <- apply(L_overall, c(1,2), mean)
# sd_L_overall <- apply(L_overall, c(1,2), sd)

mean_g_families <- apply(g_families, c(1,2), mean)
sd_g_families <- apply(g_families, c(1,2), sd)
# mean_r_families <- apply(r_families, c(1,2), mean)
# sd_r_families <- apply(r_families, c(1,2), sd)
# mean_l_families <- apply(l_families, c(1,2), mean)
# sd_l_families <- apply(l_families, c(1,2), sd)

# mean_gt0_overall <- apply(gt0_overall, c(1,2), mean)
# sd_gt0_overall <- apply(gt0_overall, c(1,2), sd)

mean_absg_families <- apply(absg_families, c(1,2), mean)
sd_absg_families <- apply(absg_families, c(1,2), sd)

saveRDS(mean_N_families, "mixed_model/save_samples/mean_N_families.rds")
saveRDS(sd_N_families, "mixed_model/save_samples/sd_N_families.rds")
# saveRDS(mean_R_overall, "mixed_model/save_samples/mean_R_overall.rds")
# saveRDS(sd_R_overall, "mixed_model/save_samples/sd_R_overall.rds")
# saveRDS(mean_L_overall, "mixed_model/save_samples/mean_L_overall.rds")
# saveRDS(sd_L_overall, "mixed_model/save_samples/sd_L_overall.rds")

saveRDS(mean_g_families, "mixed_model/save_samples/mean_g_families.rds")
saveRDS(sd_g_families, "mixed_model/save_samples/sd_g_families.rds")
# saveRDS(mean_r_families, "mixed_model/save_samples/mean_rr_families.rds")
# saveRDS(sd_r_families, "mixed_model/save_samples/sd_rr_families.rds")
# saveRDS(mean_l_families, "mixed_model/save_samples/mean_lr_families.rds")
# saveRDS(sd_l_families, "mixed_model/save_samples/sd_lr_families.rds")

# saveRDS(mean_gt0_overall, "mixed_model/save_samples/mean_gt0_overall.rds")
# saveRDS(sd_gt0_overall, "mixed_model/save_samples/sd_gt0_overall.rds")

saveRDS(mean_absg_families, "mixed_model/save_samples/mean_absg_families.rds")
saveRDS(sd_absg_families, "mixed_model/save_samples/sd_absg_families.rds")

## Run check over time series
# route = 1000
# sample = 30
# chain = 2
# all.equal(N_final[route,,sample,chain],
#           unname(sp[[chain]][,dimnames(sp[[chain]])[[2]][grepl(paste0("^N\\[",route,","), dimnames(sp[[chain]])[[2]])]][sample,]))

