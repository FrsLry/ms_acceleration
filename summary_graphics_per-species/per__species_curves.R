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
filename <- filelist[1]
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
    smp.mat <- matrix(rnorm(N, mean.mat, sd.mat), # take a draw from the normal
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

  # extract matrix of C from the data
  spec <- dat.all[species][[species]]
  Cmat <- spec$C#[,-1]

  # extract S, N and R from the summary data frame
  S <- summ[grepl("S", row.names(summ)),]
  N <- summ[grepl("N", row.names(summ)),]
  R <- summ[grepl("R", row.names(summ)),]

  # extract global parameters
  gamma50 <- summ[grepl("gamma", row.names(summ)),]$"50%"
  gamma25 <- summ[grepl("gamma", row.names(summ)),]$"2.5%"
  gamma975 <- summ[grepl("gamma", row.names(summ)),]$"97.5%"
  phi50 <- summ[grepl("phi", row.names(summ)),]$"50%"
  phi25 <- summ[grepl("phi", row.names(summ)),]$"2.5%"
  phi975 <- summ[grepl("phi", row.names(summ)),]$"97.5%"

  # estimated survival S (number of surviving individuals)
  Smat.mean <- matrix(S$mean, ncol = 34, nrow = 1033, byrow=FALSE)
  Smat <- matrix(S$"50%", ncol = 34, nrow = 1033, byrow=FALSE)
  Smat.sd <- matrix(S$mean, ncol = 34, nrow = 1033, byrow=FALSE)

  # recruitment R (number of recruited individuals)
  Rmat.mean <- matrix(R$mean, ncol = 34, nrow = 1033, byrow=FALSE)
  Rmat <- matrix(R$"50%", ncol = 34, nrow = 1033, byrow=FALSE)
  Rmat.sd <- matrix(R$sd, ncol = 34, nrow = 1033, byrow=FALSE)

  # abundance N
  Nmat.mean <- matrix(N$mean, ncol = 35, nrow = 1033, byrow=FALSE)
  Nmat <- matrix(N$"50%", ncol = 35, nrow = 1033, byrow=FALSE)
  Nmat.sd <- matrix(N$sd, ncol = 35, nrow = 1033, byrow=FALSE)

  # check if they the conversion to matrix went well
  Smat[1:5, 1:5]
  S[1:5,]

  # index of time series with min and max mean N values
  # (this will be helpful to identify how S and R match N)
  max.ind<-which( rowMeans(Nmat) %in% max(rowMeans(Nmat)) )
  min.ind<-which( rowMeans(Nmat) %in% min(rowMeans(Nmat)) )

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
  par.col <- "#e41a1c" # marking estimated parameters of the model (phi, gamma)
  reg.col <- "#ff7f00" # marking the polynomial regression line
  alpha <- 0.7 # transparency of points or lines


  # ----------------------------------
  # Graphical output
  png(filename=paste("summary_graphics_per-species/per_species_curves/", species, ".png", sep=""),
      width=1500, height=1500, res=120)

    par(mfrow=c(4,4))

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
            main = "Est. N")
    lines(Nmat[max.ind,], col = max.col, lwd=2)
    lines(Nmat[min.ind,], col = min.col, lwd=2)
    trend.plotter(Nmat, reg.col=reg.col)


    # raw observed counts in the data C
    matplot(t(Cmat), col=gray(0.3, alpha = alpha),
            type = "l",lty = 1, ylab = "Ct", xlab="t",
            main = "Obs. counts",
            ylim=c(0, max(Nmat)))
    trend.plotter(Cmat, reg.col=reg.col)

    # plot N against C
    plot(Cmat, Nmat, col=gray(0.3, alpha = alpha),
         xlab = "Ct", ylab = "Nt", main = "Est. N vs obs. counts")
    abline(a=0, b=1)

    # survival S vs N
    plot(Nmat[,-ncol(Nmat)], Smat, col=gray(0.3, alpha = 0.5),
         xlab = "Nt", ylab = "St+1",
         main = "Survival vs N")
    abline(a=0, b=1)

    # survival S
    matplot(t(Smat), col=gray(0.3, alpha = alpha),
            type = "l",lty = 1, ylab = "St+1", xlab="t+1",
            main = "Survival", ylim=c(0, max(Nmat)))
    lines(Smat[max.ind,], col = max.col, lwd=2)
    lines(Smat[min.ind,], col = min.col, lwd=2)
    trend.plotter(Smat, reg.col=reg.col)

    # survival rate (red line is phi)
    SdivNmat <- Smat/Nmat[,-ncol(Nmat)]
    matplot(t(SdivNmat), col=gray(0.3, alpha = 0.5), ylim = c(0,1),
            type = "l",lty = 1, ylab = "St+1 / Nt", xlab="t+1",
            main = "Survival rate")
    lines(SdivNmat[max.ind,], col = max.col, lwd=2)
    lines(SdivNmat[min.ind,], col = min.col, lwd=2)
    abline(h = phi50, lwd= 2, col=par.col)
    abline(h = phi25, lwd= 1, col=par.col, lty=2)
    abline(h = phi975, lwd= 1, col=par.col, lty=2)
    trend.plotter(SdivNmat, reg.col=reg.col)

    # loss rate
    Lmat <- Nmat[,-ncol(Nmat)] - Smat
    LdivNmat <- Lmat/Nmat[,-ncol(Nmat)]
    matplot(t(LdivNmat), col=gray(0.3, alpha = alpha), ylim = c(0,1),
            type = "l",lty = 1, ylab = "Lt+1 / Nt", xlab="t+1",
            main = "Loss rate")
    lines(LdivNmat[max.ind,], col = max.col, lwd=2)
    lines(LdivNmat[min.ind,], col = min.col, lwd=2)
    trend.plotter(LdivNmat, reg.col=reg.col)

    # artificial recruitment as a random draw from poission
    RmatArti <- Rmat
    RmatArti[] <- rpois(n=ncol(Rmat)*nrow(Rmat), lambda=gamma50)
    matplot(t(RmatArti), col=gray(0.3, alpha = alpha),
            type = "l",lty = 1, ylab = "Rt+1", xlab="t+1",
            ylim=c(0, max(Rmat)),
            main = "Simul. Poisson recruitment")
    trend.plotter(RmatArti, reg.col=reg.col)

    # recruitment R
    matplot(t(Rmat), col=gray(0.3, alpha = alpha),
            type = "l",lty = 1, ylab = "Rt+1", xlab="t+1",
            main = "Recruitment")
    lines(Rmat[max.ind,], col = max.col, lwd=2)
    lines(Rmat[min.ind,], col = min.col, lwd=2)
    abline(h = gamma50, lwd= 2, col=par.col)
    abline(h = gamma25, lwd= 1, col=par.col, lty=2)
    abline(h = gamma975, lwd= 1, col=par.col, lty=2)
    trend.plotter(Rmat, reg.col=reg.col)

    # recruitment rate
    RdivNmat <- Rmat/Nmat[,-ncol(Nmat)]
    RdivNmat[is.infinite(RdivNmat)] <- NA
    RdivNmat[is.nan(RdivNmat)] <- NA
    matplot(t(RdivNmat), col=gray(0.3, alpha = alpha),
            type = "l",lty = 1, ylab = "Rt+1 / Nt", xlab="t+1",
            main = "Recrutiment rate")
    lines(RdivNmat[max.ind,], col = max.col, lwd=2)
    lines(RdivNmat[min.ind,], col = min.col, lwd=2)
    trend.plotter(RdivNmat, reg.col=reg.col)

    # N against R/N
    plot(Nmat[,-ncol(Nmat)], RdivNmat, col=gray(0.3, alpha = alpha),
         xlab = "Nt", ylab = "Rt+1 / Nt",
         main = "N vs recruitment rate")

    # ================ AGGREGATION ACCROSS TIME SERIES =========================
    Ctot <- colSums(Cmat, na.rm = TRUE)
    Ntot <- colSums(Nmat, na.rm = TRUE)
    Stot <- colSums(Smat, na.rm = TRUE)
    Rtot <- colSums(Rmat, na.rm = TRUE)

    # all the count-based variables
    plot(1:35, Ntot, type = "l", lty = 1, lwd = 2, ylim = c(0, max(Ntot)),
         ylab = "Number of individuals", xlab="t",
         main = "Total sums of count variables")
      Nunc <- aggregate.uncertainty(Nmat, Nmat.sd, Nrep = 500)
      lines(1:35, Nunc[1,], lwd = 1, lty = 1)
      lines(1:35, Nunc[2,], lwd = 1, lty = 1)

      lines(1:35, Ctot, lwd = 2, lty = 1, col = "grey")
      # total survival
      lines(2:35, Stot, lwd = 2, lty = 3)
      Sunc <- aggregate.uncertainty(Smat, Smat.sd, Nrep = 500)
      lines(2:35, Sunc[1,], lwd = 1, lty = 3)
      lines(2:35, Sunc[2,], lwd = 1, lty = 3)
      # total recruitment
      lines(2:35, Rtot, lwd = 2, lty = 4)
      Runc <- aggregate.uncertainty(Rmat, Rmat.sd, Nrep = 500)
      lines(2:35, Runc[1,], lwd = 1, lty = 4)
      lines(2:35, Runc[2,], lwd = 1, lty = 4)

    legend("topright", legend = c("Nt", "Ct", "St", "Rt"),
           lty = c(1,1,3,4),
           lwd = c(2,2,2,2),
           col = c("black","grey","black","black"))


    # total recruitment rate
    plot(2:35, Rtot/Ntot[-length(Ntot)], type = "l", lwd = 2,
         ylab = "Rt+1 / Nt", xlab = "t+1", lty = 4,
         main = "Total recruitment rate")

    # total survival rate
    plot(2:35, Stot/Ntot[-length(Ntot)], type = "l", lwd = 2,
         ylab = "St+1 / Nt", xlab = "t+1", lty = 3,
         main = "Total survival rate")

    # total loss rate
    plot(2:35, (Ntot[-length(Ntot)]-Stot)/Ntot[-length(Ntot)], type = "l", lwd = 2,
         ylab = "Lt+1 / Nt", xlab = "t+1", lty = 1,
         main = "Total loss rate")

  dev.off()
}



