`fda.clust` (Version 0.1.1) Clustering Methods for Functional Data
Analysis
================

<!-- README.Rmd is compiled to README.md. Please edit README.Rmd and then knit to update README.md -->
<!-- Badges -->

[![Licence](https://img.shields.io/badge/licence-GPL--2-blue.svg)](https://www.gnu.org/licenses/gpl-2.0.en.html)

<!-- Logo -->

<img src="inst/figures/fda.clust.png" alt="Logo fda.clust" align="right" width="140"/>

`fda.clust` is a development package that provides specialized methods
for clustering functional data. Inspired by the functional data analysis
framework, this package applies and adapts clustering techniques—such as
functional k-means, hierarchical clustering, and mixture models—to the
continuous nature of functional data. Additionally, it includes tools
for cluster visualization and result validation.

Version `0.1.1` is the first development release, focusing on
implementing the main methods and establishing the basic framework for
functional clustering.

## Package Overview

The `fda.clust` package offers:

- **Functional clustering methods**:
  - functional k-means,
  - hierarchical clustering,
  - functional dbscan.
- **Visualization tools**: Plots and graphical representations to
  explore and understand the resulting clusters.
- **Cluster validation**: Metrics for evaluating the quality and
  stability of the obtained partitions.

The goal is to provide a comprehensive toolbox for users interested in
segmenting functional data—be it for research or for applying advanced
statistical methodologies in real-world problems.

## Installation

Once available, you can install the development version of `fda.clust`
from GitHub. For example:

``` r
# install.packages("devtools") # if you don't have it already
library(devtools)
devtools::install_github("username/fda.clust")
# Or:
# remotes::install_github("username/fda.clust")
```

\##Example

Below is a simple example that shows how you might load the package and
apply a clustering method (note: this is a placeholder and may need
adjustment once the functions are implemented):

``` r
library(fda.clust)

# Suppose 'fd_data' is a functional data object
# cluster_result <- fda_clust_kmeans(fd_data, centers = 3)
# plot(cluster_result)
```

## Issues & Feature Requests

For reporting issues, bugs, feature requests, etc., please use the
[Github Issues](https://github.com/moviedo5/fda.clust/issues). Your
contributions and feedback are always welcome.

## Documentation

A hands on introduction to can be found in the reference
[vignette](https://www.jstatsoft.org/index.php/jss/article/view/v051i04/v51i04.pdf).

## References

- Febrero-Bande, M. and Oviedo de la Fuente, M. (2012). Statistical
  Computing in Functional Data Analysis: The R Package fda.usc. *Journal
  of Statistical Software*, **51**(4):1-28,
  [DOI](http://www.jstatsoft.org/v51/i04/)

<!-- 
# data
# generador de datos usados en classif.DD/TFM
fnt.sim() parabola ojo
Simulacion DF.R estan los tres modelos,
&#10;# growth
# 
Medidas de bondad del ajuste
FV2006
fpc:::cluster.stats()$dunn
fpc:::dbscan vs dbscan
som kohonen
cluster:::silhouette
Abre el archivo DESCRIPTION con un editor de texto sin formato (como VS Code, Notepad++ o RStudio).
Revisa cada linea y asegurate de que todas las lineas tengan el formato correcto, especialmente las que contienen texto multilinea.
Verifica errores comunes:
Asegurate de que cada linea este separada por una nueva linea (\n).
Revisa los campos multilinea. Por ejemplo, la Description debe estar sangrada si ocupa varias lineas, por ejemplo:
dcf
Copiar codigo
Description: This package provides tools for clustering functional data. 
  The clustering methods are based on the use of distances 
  between curves.
&#10;# Regenerar un archivo DESCRIPTIO usethis::use_description()
Detectar caracteres no visibles:
lines <- readLines("DESCRIPTION")
print(lines)
# codificado en UTF-8.
tools::showNonASCIIfile("DESCRIPTION")
rm(list = c("rproc2mu"))
rm(list = c("rproc2clust"))
&#10;library(roxygen2)
# setwd("C:/Users/Manuel Oviedo/github/fda.clust")
getwd()
#pkgbuild::compile_dll()
roxygenize()
devtools::document() 
library(tools)
tools::checkRd("man/fdbscan.Rd")
tools::checkRd("man/fhcust.Rd")
tools::checkRd("man/fmeanshift.Rd")
tools::checkRd("man/fhcust.Rd")
tools::checkRd("man/fkmeans.Rd")
tools::checkRd("man/rproc2mu.Rd")
tools::checkRd("man/rproc2clust.Rd")
library(devtools)
devtools::build()
devtools::check(manual = TRUE)  # problemas con miktex
devtools::install()
&#10;#  setwd("C:/Users/Manuel Oviedo/github")
build_manual(pkg = "fda.clust", path = NULL)
devtools::install_github("moviedo5/fda.usc",auth_user="moviedo5")
# devtools::install_github("moviedo5/fda.usc",auth_user="moviedo5")
R CMD build fda.clust
R CMD check fda.clust_0.1.1.tar.gz --as-cran
R CMD INSTALL fda.clust_0.1.1.tar.gz --build
R-wind-builder fda.clust_0.1.1.tar.gz --as-cran
R CMD Rd2pdf fda.clust
library(pkgdown)
# usethis::use_pkgdown()
# Build website:
#pkgdown::build_site()
build_site(new_process = TRUE)
-->
