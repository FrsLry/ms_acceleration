# Script that processes all of the summary files for all species, and plots
# the estimated S, R, N in time
# Author: Petr Keil
# ------------------------------------------------------------------------------

rm(list=ls())

library(plotrix) # to plot data frame as a figure

# ------------------------------------------------------------------------------

# read the raw JAGS data
dat.all <- readRDS("data/jags_data_par.rds")

# read list of summary files for all species
filelist <- list.files(path = "summary_models")

# example species
filename <- filelist[10]
filename

# ------------------------------------------------------------------------------

# function that plots a polynomial regression of two temporal matrices
trend.plotter <- function(mat, poly.deg = 3, family = "gaussian", reg.col=reg.col)
{
  # matrix for time values
  year <- mat
  year[] <- rep(1:ncol(mat), each=nrow(mat))
  # fit a polynomial
  dat <- na.omit(data.frame(mat = as.vector(mat), year = as.vector(year)))

  mod <- glm(mat ~ poly(year, poly.deg), data = dat, family = family)
  # make predictions
  prds <- predict(mod, newdata=data.frame(year = min(year):max(year)))
  # plot the line
  lines(1:ncol(mat), prds, col = reg.col, lwd=2)
}


# ------------------------------------------------------------------------------
# function that takes a matrix of means and matrix of SDs, aggregates time series
# across sites, and plots the aggregated (summed) values together with their
# Bayesian unceratinty
aggregate.uncertainty <- function(mean.mat, sd.mat, Nrep = 500)
{
  N <- nrow(mean.mat)*ncol(mean.mat) # N of all cells in the matrix
  smp.sums <- matrix(ncol = ncol(mean.mat),
                     nrow = Nrep) # empty matrix for sample storage

  for(i in 1:Nrep)
  {
    smp.mat <- matrix(rnorm(N, as.vector(mean.mat), as.vector(sd.mat)), # take a draw from the normal
                      nrow = nrow(mean.mat),
                      ncol = ncol(mean.mat),
                      byrow = FALSE)
    #smp.mat[smp.mat < 0] <- 0 # replace negative values with 0
    smp.sums[i,] <- colSums(smp.mat, na.rm = TRUE) # aggregate across sites
  }
  res <- apply(X=smp.sums, MARGIN=2, FUN=quantile, probs = c(0.025, 0.975))

  return(res)
}
# x <- aggregate.uncertainty(Smat, Smat.sd, Nrep=100)


# ------------------------------------------------------------------------------


# main loop through all species
for(filename in filelist)
{
  # read the summary file
  summ <- readRDS(paste('summary_models/', filename, sep=''))
  mean.rhat <- mean(summ$Rhat[1:20])

  # clean species name for nice file output
  species <- gsub(".rds", replacement = "", x = filename)
  species <- gsub("summary_", replacement = "", x = species)
  print(species)

  # extract matrix of C (observed counts) from the data
  spec <- dat.all[species][[species]]
  Cmat <- spec$C#[,-1]

  # extract S, N and R from the summary data frame
  N <- summ[grepl("N", row.names(summ)),]

  # abundance N
  Nmat.mean <- matrix(N$mean, ncol = 35, nrow = 1033, byrow=FALSE)
  Nmat <- matrix(N$"50%", ncol = 35, nrow = 1033, byrow=FALSE)
  Nmat.sd <- matrix(N$sd, ncol = 35, nrow = 1033, byrow=FALSE)


  # index of time series with min and max mean N values
  # (this will be helpful to identify how S and R match N)
  max.ind<-which( rowMeans(Nmat) %in% max(rowMeans(Nmat)) )[1]
  min.ind<-which( rowMeans(Nmat) %in% min(rowMeans(Nmat)) )[1]

  # summary table of the main  parameters
  sm <- summ[c(1:20),c("mean", "sd", "2.5%", "50%", "97.5%", "Rhat")]
  sm <- data.frame(parameter = rownames(sm),
                   mean = round(sm$mean,2),
                   SD = round(sm$sd,2),
                   q2.5 = round(sm$"2.5%",2),
                   q50 = round(sm$"50%", 2),
                   q97.5 = round(sm$"97.5%",2),
                   Rhat = sm$Rhat)
  converged <- max(sm$Rhat, na.rm=T) < 1.1 # has the model converged

  # parameters of the graphical output
  max.col <- "#4daf4a" # marking time series with max mean abundance
  min.col <- "#377eb8" # marking time series with min mean abundance
  reg.col <- "#ff7f00" # marking the polynomial regression line
  alpha <- 0.7 # transparency of points or lines


  # ----------------------------------
  # Graphical output
  png(filename=paste("summary_graphics_per-species/per_species_curves/", species, ".png", sep=""),
      width=1900, height=450, res=120)

    par(mfrow=c(1, 5))

    # plot the parameters as a table in the plot
    par(mar = c(3,0.6,0,0))
    plot(c(0, 0.15), c(0, 1), ann = F, type = 'n', xaxt = 'n', yaxt = 'n',
         frame = F, main = species)
    plotrix::addtable2plot(0.01, 0, sm, xjust = 0, bty = 'o',
                           cex=0.8, vlines=T, xpad=0.2)
    # mark species that haven't converged
    if(converged == FALSE)
    {
      text(x=0.06, y=0.5, "NOT\nCONVERGED", col="red", cex = 3)
    }

    # reset figure margins
    par(mai=c(0.7,0.6, 0.3, 0.1))

    # estimated abundance N
    matplot(t(Nmat), col=gray(0.3, alpha = alpha),
            type = "l",lty = 1, ylab = "Nt", xlab="t",
            main = "Est. N", ylim=c(0, max(Nmat)))
    lines(Nmat[max.ind,], col = max.col, lwd=2)
    lines(Nmat[min.ind,], col = min.col, lwd=2)
    trend.plotter(Nmat, reg.col=reg.col)


    # raw observed counts in the data C
    matplot(t(Cmat), col=gray(0.3, alpha = alpha),
            type = "l",lty = 1, ylab = "Ct", xlab="t",
            main = "Obs. counts",
            ylim=c(0, max(Nmat)))
    trend.plotter(Cmat, reg.col=reg.col)
    lines(Cmat[max.ind,], col = max.col, lwd=2)
    lines(Cmat[min.ind,], col = min.col, lwd=2)

    # plot N against C
    plot(Cmat, Nmat, col=gray(0.3, alpha = alpha),
         xlab = "Ct", ylab = "Nt", main = "Est. N vs obs. counts")
    abline(a=0, b=1)


    # ================ AGGREGATION ACCROSS TIME SERIES =========================
    # total estimated abundance summed accross all time series:
    Ntot <- colSums(Nmat, na.rm = TRUE)
    # uncertaint about the estimated abundance
    Nunc <- aggregate.uncertainty(Nmat, Nmat.sd, Nrep = 500)

    # total observed abundance summed across all time series
    ObsTot <- colSums(Cmat, na.rm = TRUE)


    plot(1:35, Ntot, type = "l", lty = 1, lwd = 2, ylim = c(0, max(Nunc[2,])),
         ylab = "Number of individuals", xlab="t",
         main = "Total N")
      lines(1:35, Nunc[1,], lwd = 1, lty = 1)
      lines(1:35, Nunc[2,], lwd = 1, lty = 1)
      lines(1:35, ObsTot, lwd=2, col="red")

  dev.off()
}



