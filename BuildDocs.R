knitr::knit("README.Rmd") # Knit README
devtools::document() # Rebuild Documentation
devtools::build_vignettes() #Build Vignettes
pkgdown::build_site() #Build pkgdown site
