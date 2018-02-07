\name{run_landclim}
\alias{run_landclim}

\title{
Run LandClim
}
\description{

}
\usage{
run_landclim(control_file)
}

\arguments{
  \item{control_file}{Name of the control file (e.g. control.xml).

}
}
\details{
The working directory needs to be set to the level above "Input" and "Output" folders of LandClim, thus at the "site" level. Until now, the folder names are hard-coded, thus need to be "Input" and "Output", not e.g. input! Feel free to fix this!

It is assumed that LandClim is installed in a way, that the command `landclim` is available in the command line (regardless of the current working directory).

If that is not the case, you can do the following in Ubuntu:
Add the following line to the end of the file \code{.bashrc} in your home directory: \code{PATH="$PATH:/path/to/landclim"} (replace \code{/path/to/landclim} with the real path of the directory in which LandClim can be found).
Then run the following command in the command line: \code{source ~/.bashrc}

Now you should be able to simply enter \code{landclim} in the command line to start LandClim. And calling \code{run_landclim()} in R should also work.
}

\references{
### Add LandClim refrences
}
\author{
Klara Dolos
}
