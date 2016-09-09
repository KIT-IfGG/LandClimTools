# Windows

*we currently don't have a guide on setting up the environment on Windows machines*

# Ubuntu 16.04 LTS (Xenial Xerus)

## Setting up the development environment

First install the required software:
* git version control system
* PDFLaTeX (part of the `texlive` LaTeX-distribution)
* GNU R statistical computation and graphics system
* GNOME XML library `libxml2` ≥ 2.6.3 (required by the R-package `XML`)
* Cartographic projection library `proj.4` ≥ 4.4.9 (required by the R-package `rgdal`)
* Geospatial Data Abstraction Library `libgdal` ≥ 1.6.3 (required by the R-package `rgdal`)

All of those are in the main Ubuntu repositories, so they can simply be installed by executing the following command in the command line:
```bash
$ sudo apt install git texlive-latex-base r-base libxml2-dev libproj-dev libgdal-dev
```

---

Then you'll need to install the R-Packages required.

For this start an R workspace with the following command:
```bash
$ R
```

In the R command line, type the following commands. When you are asked, if you want to create a personal library, you can answer with `y` for yes. And if you are asked to select a mirror, choose the one nearest to you (e.g. in Germany the one at the university of Münster):
```r
> install.packages("rgdal", dependencies = TRUE)
> install.packages("raster", dependencies = TRUE)
> install.packages("XML", dependencies = TRUE)
```

Quit the R command line with the following command:
```r
> q()
```
When it asks you if you want to save the workspace image, you can answer with `n` for no.

---

Now you have everything to start building and developing the `LandClimTools` package.

You can download the package from Github with the following commands. In the first command replace the dummy-path with an existing directory where you want the package to be:
```bash
$ cd /home/your-username/directory-containing-the-package/
$ git clone https://github.com/KIT-IfGG/LandClimTools.git
```

---

That's it! You are ready to build and develop. :tada: :thumbsup:

If you want, you can use an IDE (integrated development environment) like RStudio, VisualStudio or the StatET plugin for Eclipse.

These provide additional editing functionality like autocompletion or syntax highlighting, but they are not necessary to work with the code. Any text editor will work.

We recommend RStudio, a project file for it is included in the source code.

## Building the R-package
Do the following in the command line. In the first command replace the dummy-path with the correct path to the directory where you have placed the source code of `LandClimTools`:
```bash
$ cd /home/your-username/directory-containing-the-package/LandClimTools
$ mkdir build
$ cd build
$ R CMD build ..
```
After you have done that, you'll find a file called similar to `LandClimTools_1.2.3.tar.gz` (with the current version number) in the `build/` directory.

## Checking the package for errors

After building the package, you can check for errors in it by running the following command inside the `build/` directory:
```bash
$ R CMD check LandClimTools_1.2.3.tar.gz
```
This might take a while and some windows may appear on your screen.
