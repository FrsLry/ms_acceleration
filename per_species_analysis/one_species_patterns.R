# Script that processes all of the summary files for all species, and plots
# the estimated S, R, N in time
# Author: Petr Keil
# ------------------------------------------------------------------------------

rm(list=ls())


library(plotrix)
# ------------------------------------------------------------------------------

# read the raw JAGS data
dat.all <- readRDS("data/jags_data_par.rds")

# read list of summary files for all species
splist <- list.files(path = "summary_models")


# ------------------------------------------------------------------------------

# function that plots a polynomial regression of two temporal matrices
trend.plotter <- function(mat)
{
  # matrix for time values
  year <- mat
  year[] <- rep(1:ncol(mat), each=nrow(mat))
  # fit a polynomial
  dat <- na.omit(data.frame(mat = as.vector(mat), year = as.vector(year)))

  mod <- lm(mat ~ poly(year,3), data = dat)
  # make predictions
  prds <- predict(mod, newdata=data.frame(year = min(year):max(year)))
  # plot the line
  lines(1:ncol(mat), prds, col = "blue", lwd=2)
}


# ------------------------------------------------------------------------------


# main loop through all species
for(species in splist)
{
  # read the summary file
  summ <- readRDS(paste('summary_models/', species, sep=''))
  mean.rhat <- mean(summ$Rhat[1:20])

  # clean species name for nice file output
  species <- gsub(".rds", replacement = "", x = species)
  species <- gsub("summary_", replacement = "", x = species)

  # extract matrix of C from the data
  spec <- dat.all[species][[species]]
  Cmat <- spec$C[,-1]

  # extract S, N and R from the summary data frame
  S <- summ[grepl("S", row.names(summ)),]
  N <- summ[grepl("N", row.names(summ)),]
  R <- summ[grepl("R", row.names(summ)),]

  # extract global parameters
  gamma <- summ[grepl("gamma", row.names(summ)),]$mean
  gamma25 <- summ[grepl("gamma", row.names(summ)),]$"2.5%"
  gamma975 <- summ[grepl("gamma", row.names(summ)),]$"97.5%"
  phi <- summ[grepl("phi", row.names(summ)),]$mean
  phi25 <- summ[grepl("phi", row.names(summ)),]$"2.5%"
  phi975 <- summ[grepl("phi", row.names(summ)),]$"97.5%"

  # estimated survival S (number of surviving individuals)
  Smat <- matrix( ncol = 34, nrow = 1033, byrow=TRUE)
  Smat[] <- S$mean

  # recruitment R (number of recruited individuals)
  Rmat <- matrix( ncol = 34, nrow = 1033, byrow=TRUE)
  Rmat[] <- R$mean

  # artificial recruitment as a random draw from poission
  RmatArti <- Rmat
  RmatArti[] <- rpois(n=ncol(Rmat)*nrow(Rmat), lambda=gamma)

  # abundance N
  Nmat <- matrix( ncol = 35, nrow = 1033, byrow=TRUE)
  Nmat[] <- N$mean
  Nmat <- Nmat[,-1]

  # proportion of surviving individuals
  SdivNmat <- Smat/Nmat

  # proportion of recruited individuals
  RdivNmat <- Rmat[,-ncol(Rmat)] / Nmat[,-1]
  RdivNmat[is.infinite(RdivNmat)] <- NA
  RdivNmat[is.nan(RdivNmat)] <- NA

  # ----------------------------------
  # Graphical output
  png(filename=paste("per_species_curves/", species, ".png", sep=""),
      width=1800, height=820, res=120)

    par(mfrow=c(2,5))

    # plot the parameters as a table in the plot
    sm <- summ[c(1:20),c("mean", "sd", "2.5%", "97.5%", "Rhat")]
    sm <- data.frame(parameter = rownames(sm),
                     mean = round(sm$mean,2),
                     SD = round(sm$sd,2),
                     q2.5 = round(sm$"2.5%",2),
                     q97.5 = round(sm$"97.5%",2),
                     Rhat = sm$Rhat)
    par(mar = c(5,0.6,0,0))
    plot(c(0, 0.1), c(0, 1), ann = F, type = 'n', xaxt = 'n', yaxt = 'n',
         frame = F, main = species)
    plotrix::addtable2plot(0.01, 0, sm, xjust = 0, bty = 'o',
                           cex=0.8, vlines=T, xpad=0.2)


    # reset figure margins
    par(mai=c(0.7,0.6, 0.3, 0.1))

    # raw observed counts in the data C
    matplot(t(Cmat), col=gray(0.3, alpha = 0.5),
            type = "l",lty = 1, ylab = "Ct", xlab="year",
            main = "Obs. counts",
            ylim=c(0, max(Nmat)))
    trend.plotter(Cmat)


    # estimated abundance N
    matplot(t(Nmat), col=gray(0.3, alpha = 0.5),
            type = "l",lty = 1, ylab = "Nt", xlab="year",
            main = "Est. N")
    trend.plotter(Nmat)


    # plot N against C
    plot(Cmat, Nmat, col=gray(0.3, alpha = 0.5),
         xlab = "Ct", ylab = "Nt", main = "Est. N vs obs. counts")
    abline(a=0, b=1)


    # N against R/N
    plot(Nmat[,-1], RdivNmat, col=gray(0.3, alpha = 0.5),
         xlab = "Nt+1", ylab = "Rt / Nt+1",
         main = "N vs prop. recruited")

    # artificial recruitment
    matplot(t(RmatArti), col=gray(0.3, alpha = 0.5),
            type = "l",lty = 1, ylab = "Rt", xlab="year",
            ylim=c(0, max(Rmat)),
            main = "Simul. Poisson recruitment")
    trend.plotter(RmatArti)


    # recruitment R (red line is gamma)
    matplot(t(Rmat), col=gray(0.3, alpha = 0.5),
            type = "l",lty = 1, ylab = "Rt", xlab="year",
            main = "Recruitment")
    abline(h = gamma, lwd= 2, col="red")
    abline(h = gamma25, lwd= 1, col="red", lty=2)
    abline(h = gamma975, lwd= 1, col="red", lty=2)
    trend.plotter(Rmat)


    # survival S
    matplot(t(Smat), col=gray(0.3, alpha = 0.5),
            type = "l",lty = 1, ylab = "St", xlab="year",
            main = "Survival", ylim=c(0, max(Nmat)))
    trend.plotter(Smat)



    # proportion of R in N
    matplot(t(RdivNmat), col=gray(0.3, alpha = 0.5),
            type = "l",lty = 1, ylab = "Rt / Nt+1", xlab="year",
            main = "Prop. recruited", ylim=c(0,1))
    trend.plotter(RdivNmat)


    # proportion of surviving individuals (red line is phi)
    matplot(t(SdivNmat), col=gray(0.3, alpha = 0.5), ylim = c(0,1),
            type = "l",lty = 1, ylab = "St / Nt", xlab="year",
            main = "Prop. survived")
    abline(h = phi, lwd= 2, col="red")
    abline(h = phi25, lwd= 1, col="red", lty=2)
    abline(h = phi975, lwd= 1, col="red", lty=2)
    trend.plotter(SdivNmat)

   dev.off()
}



