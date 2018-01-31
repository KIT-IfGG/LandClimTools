\name{plot_gradient}
\alias{plot_gradient}
\title{
Gradient plot
}
\description{
Plot succession or elevation gradient based on Landclim decadal output.
}
\usage{
plot_gradient(x, y, col = NULL, ...)
}

\arguments{
  \item{x}{
x values, either elevation or decades from e.g. elevation output.
}
  \item{y}{
y values as matrix. Species columns from e.g. elevation output. 
}
  \item{col}{
Color vector of length ncol(y). Use e.g. landclim_colors().
}
  \item{\dots}{
Additional arguments for function plot().
}
}

\author{
Klara Dolos
}

\examples{

dat <- read.table(system.file("elevation_biomass_out.csv", package = "LandClimTools"), sep = ",", dec = ".", header = TRUE)

species <- c("abiealba" , "piceabie", "fagusilv", "pinusilv", "querpetr")

plot_gradient(dat$decade[dat$elevation==823], dat[dat$elevation==823, c(species)], col=landclim_colors(length(species)), xlab="Decade", ylab="Biomass", ylim=c(0,400))
legend("topleft", legend = species, col = landclim_colors(length(species)), pch=15, pt.cex=2)

plot_gradient(dat$elevation[dat$decade==50], dat[dat$decade==50, c(species)], xlab="Elevation", ylab="Biomass", col=landclim_colors(length(species)), ylim=c(0,400))
legend("topleft", legend = species, col = landclim_colors(length(species)), pch=15, pt.cex=2)



}

