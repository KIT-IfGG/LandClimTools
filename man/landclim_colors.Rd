\name{landclim_colors}
\alias{landclim_colors}

\title{
LandClim color palette
}
\description{
Create a color vector for plotting.
}
\usage{
landclim_colors(n)
}

\arguments{
  \item{n}{
Number of different tree species = number of colors.
}
}

\value{
Character vector with colors.
}

\author{
Klara Dolos
}
\note{
Please feel free to enhance this palette!
}

\examples{
n <- 10
pie(rep(1, n), col=landclim_colors(n))

### Get number of available colors.
landclim_colors(99)

}

