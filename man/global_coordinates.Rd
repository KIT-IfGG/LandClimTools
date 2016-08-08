\name{global_coordinates}
\alias{global_coordinates}

\title{Convert local coordinates to global ones}
\description{takes local coordinates, the numbers of rows/columns and the width/height of each part of the grid. From that the global coordinates are calculated}

\usage{global_coordinates(x.local, y.local, row, col, a)}

\arguments{
  \item{x.local}{local x coordinate}
  \item{y.local}{local y coordinate}
  \item{row}{row index}
  \item{col}{column index}
  \item{a}{size (both width and height) of each local cell}
}

\value{the calculated global coordinates}
