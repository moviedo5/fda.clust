# fda.clust: Clustering Methods for Functional Data Analysis

[![Licence](https://img.shields.io/badge/licence-GPL--2-blue.svg)](https://www.gnu.org/licenses/gpl-2.0.en.html)

The `fda.clust` package provides specialized methods for clustering functional data, inspired by the functional data analysis (FDA) framework. This package offers tools for clustering, validation, and visualization of functional data. It allows users to work with real and simulated datasets for performance evaluation of clustering methods.

## **Features**

- **Clustering Methods**:
  - `fkmeans`: Functional k-means clustering.
  - `fdbscan`: Functional DBSCAN clustering.
  - `fmeanshift`: Functional mean-shift clustering.
  - `fhclust`: Functional hierarchical clustering.

- **Validation Metrics**:
  - **Silhouette**: Measure of cohesion and separation.
  - **Dunn**: Ratio of the smallest inter-cluster distance to the largest intra-cluster distance.
  - **Davies-Bouldin**: Average similarity between clusters.
  - **Calinski-Harabasz**: Ratio of between-cluster dispersion to within-cluster dispersion.

- **Data Simulation**:
  - `rprocKclust`: Simulate functional data for a predefined number of clusters.
  - `rprocKmu`: Simulate mean functions for multiple clusters.

- **Datasets**:
  - **ECG200**: 200 heartbeats classified as normal or myocardial infarction.
  - **ECG5000**: 5000 heartbeats classified into four groups.
  - **growth_ldata**: Longitudinal growth data from the Berkeley Growth Study.

---

## **Installation**

To install the development version from GitHub, run the following command in R:

```r
# Install the development version from GitHub
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("moviedo5/fda.clust")
```

---

## **Usage**

Below is an example demonstrating how to use the **fkmeans** function for clustering functional data.

```r
library(fda.clust)

# Load the example dataset ECG200
data(ECG200)

# Perform k-means clustering on functional data with 2 clusters
set.seed(123)
result <- fkmeans(ECG200$x, ncl = 2)

# Plot the functional data with cluster assignments
plot(ECG200$x, col = result$cluster, main = "ECG200 Clustered with fkmeans")
```

---

## **Datasets**

The package includes three functional datasets to test clustering methods:

1. **ECG200**: Electrical signals from heartbeats (2 classes: normal and myocardial infarction).
2. **ECG5000**: A larger dataset of 5000 heartbeats (4 classes).
3. **growth_ldata**: Longitudinal growth data from the Berkeley Growth Study, including the heights of boys and girls at 31 ages.

---

## **Available Functions**

### **Clustering Methods**
- **`fkmeans`**: Perform functional k-means clustering.
- **`fdbscan`**: Perform functional DBSCAN clustering.
- **`fmeanshift`**: Perform functional mean-shift clustering.
- **`fhclust`**: Perform functional hierarchical clustering.

### **Cluster Validation and Measures**
- **`fclust.measures`**: Evaluate the quality of clusters using internal indices such as silhouette, Dunn, Davies-Bouldin, and Calinski-Harabasz indices.

### **Data Simulation**
- **`rprocKclust`**: Generate functional data with known clusters for testing clustering methods.
- **`rprocKmu`**: Generate mean functions for multiple clusters.

### **Utility Functions**
- **`kmeans.assig.groups`**: Assign functional data to clusters based on distances.
- **`kmeans.centers.update`**: Update cluster centers during k-means clustering.
- **`kmeans.fd.dist`**: Calculate distances for functional k-means clustering.

---

## **Documentation**

To learn more about the **fda.clust** package, check the vignettes for a more comprehensive overview of its functionalities. You can access the vignettes directly in R:

```r
vignette("Introduction", package = "fda.clust")
vignette("Simulations", package = "fda.clust")
```


Details on specific functions are in the [reference
manual](https://moviedo5.github.io/fda.clust/index.html).
To learn more about the functions and their usage, you can refer to the
**pkgdown** documentation site. 

Vignettes:

- [**Pkg introduction**](https://moviedo5.github.io/fda.clust/articles/Introduction.html)

- [**Generation of functional clusters**](https://moviedo5.github.io/fda.clust/articles/Simulations.html).

- [**Example using Canadian Weather**](https://moviedo5.github.io/fda.clust/articles/ExampleCanadianWeather.html).

- [**Example using Spanish Weather**](https://moviedo5.github.io/fda.clust/articles/ExampleSpanishWeather.html).

- [**Example multivariate FDA usingh Spanish Weather**](https://moviedo5.gmithub.io/fda.clust/articles/MultivariateSpanishWeather.html).


---

## **Issues & Feature Requests**

For reporting issues, bugs, feature requests, etc., please use the [Github Issues](https://github.com/moviedo5/fda.clust/issues) page. Contributions and feedback are always welcome.

---

## **References**

- Febrero-Bande, M. and Oviedo de la Fuente, M. (2012). Statistical Computing in Functional Data Analysis: The R Package fda.usc. *Journal of Statistical Software*, **51**(4):1-28, [DOI](http://www.jstatsoft.org/v51/i04/).


<!-- 
# data
# generador de datos usados en classif.DD/TFM
fnt.sim() parabola ojo
Simulacion DF.R estan los tres modelos,
# growth
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
# Regenerar un archivo DESCRIPTIO usethis::use_description()
Detectar caracteres no visibles:
lines <- readLines("DESCRIPTION")
print(lines)
# codificado en UTF-8.
tools::showNonASCIIfile("DESCRIPTION")
rm(list = c("rproc2mu"))
rm(list = c("rproc2clust"))

# Limpia los archivos previos
unlink("NAMESPACE")
unlink("man", recursive = TRUE)


library(roxygen2)
# setwd("C:/Users/Manuel Oviedo/github/fda.clust")
getwd()
#pkgbuild::compile_dll()
roxygenize()
devtools::document() 
library(tools)
tools::checkRd("man/fdbscan.Rd")
tools::checkRd("man/fmeanshift.Rd")
tools::checkRd("man/mfkmeans.Rd")

tools::checkRd("man/fhclust.Rd")
tools::checkRd("man/fkmeans.Rd")
#tools::checkRd("man/rproc2mu.Rd")
#tools::checkRd("man/rproc2clust.Rd")
tools::checkRd("man/fclust.measures.Rd")
tools::checkRd("man/ECG5000.Rd")
tools::checkRd("man/ECG200.Rd")
tools::checkRd("man/growth_ldata.Rd")


remove.packages("fda.clust") 
devtools::install()
# library(fda.clust)
# exists("mfkmeans", where = asNamespace("fda.clust"))
#system.file(package = "fda.clust")

# 1. Inicializar la web del paquete
library(pkgdown)
#pkgdown::init_site()

# 2. Construir la web
pkgdown::build_site()

# Step 1: Update documentation and vignettes
devtools::document()
devtools::build_vignettes()

# Step 2: Build pkgdown site
pkgdown::build_site()

# Step 3: Rebuild package
devtools::build()

# Step 4: Check package
devtools::check()
# devtools::check(manual = TRUE) 

devtools::install()


#  setwd("C:/Users/Manuel Oviedo/github")
build_manual(pkg = "fda.clust", path = NULL)

unlink(file.path(tempdir(), "lastMiKTeXException"))
unlink(tempdir(), recursive = TRUE)

devtools::install_github("moviedo5/fda.usc",auth_user="moviedo5")
# devtools::install_github("moviedo5/fda.usc",auth_user="moviedo5")
R CMD build fda.clust
R CMD check fda.clust_0.1.1.tar.gz --as-cran
R CMD INSTALL fda.clust_0.1.1.tar.gz --build
R-wind-builder fda.clust_0.1.1.tar.gz --as-cran
R CMD build --resave-data fda.clust
     
R CMD Rd2pdf fda.clust

library(pkgdown)
# usethis::use_pkgdown()
# Build website:
#pkgdown::build_site()
build_site(new_process = TRUE)

devtools::build_vignettes()
#unlink("inst/doc", recursive = TRUE)
#devtools::build_vignettes()

-->

