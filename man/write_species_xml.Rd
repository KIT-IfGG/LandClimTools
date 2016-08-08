\name{write_species_xml}
\alias{write_species_xml}

\title{Write LandClim species parameter file}
\description{Writes the given data into the given file in the XML-format for species of LandClim}
\usage{
  write_species_xml(x, file)
}

\arguments{
  \item{x}{
    %%     ~~Describe \code{x} here~~
  }
  \item{file}{
    %%     ~~Describe \code{file} here~~
  }
}

\examples{
## The function is currently defined as
function (x, file)
{
    names <- colnames(x)
    suppressWarnings(tr <- xmlTree("species"))
    for (i in 1:NROW(x)) {
        tr$addTag("set", close = FALSE)
        for (j in names) {
            tr$addTag(j, as.character(x[i, j]))
        }
        tr$closeTag()
    }
    invisible(saveXML(tr$value(), file = file))
  }
}
