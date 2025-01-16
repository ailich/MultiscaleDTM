r<- erupt()
set.seed(5)
r[sample(1:ncell(r), size = 500, replace = FALSE)]<- NA

# groundhog.library("MultiscaleDTM", date = "2024-7-1") # Surfaces except for SAPA generated with this


test_that("Test SlpAsp_queen", {
  test<- SlpAsp(r = r, w = c(5,7), unit = "degrees", method = "queen", metrics = c("slope", "aspect", "eastness", "northness"), na.rm=TRUE) |> values(mat=TRUE)
  expect<- readRDS(system.file("testdata", "slp_asp_queen.RDS", package= "MultiscaleDTM"))
  expect_equal(test, expect)
})

test_that("Test SlpAsp_rook", {
  test<- SlpAsp(r = r, w = c(5,7), unit = "degrees", method = "queen", metrics = c("slope", "aspect", "eastness", "northness"), na.rm=TRUE) |> values(mat=TRUE)
  expect<- readRDS(system.file("testdata", "slp_asp_rook.RDS", package= "MultiscaleDTM"))
  expect_equal(test, expect)
})

test_that("Test Qfit", {
  test<- Qfit(r, w = c(5,7), unit = "degrees", metrics = c("elev", "qslope", "qaspect", "qeastness", "qnorthness", "profc", "planc", "twistc", "meanc", "maxc", "minc", "features"), na.rm = TRUE) |> values(mat=TRUE)
  expect<- readRDS(system.file("testdata", "qmetrics.RDS", package= "MultiscaleDTM"))
  expect_equal(test, expect)
})

test_that("Test VRM", {
  test<- VRM(r, w=c(5,7), na.rm = TRUE) |> values(mat=TRUE)
  expect<- readRDS(system.file("testdata", "vrm.RDS", package= "MultiscaleDTM"))
  expect_equal(test, expect)
})


test_that("Test SAPA", {
  test<- SAPA(r, w=c(5,7), slope_correction = FALSE, na.rm=TRUE) |> values(mat=TRUE)
  expect<- readRDS(system.file("testdata", "sapa.RDS", package= "MultiscaleDTM"))
  expect_equal(test, expect)
})

test_that("Test AdjSD", {
  test<-  AdjSD(r, w=c(5,7), na.rm = TRUE) |> values(mat=TRUE)
  expect<- readRDS(system.file("testdata", "adj_sd.RDS", package= "MultiscaleDTM"))
  expect_equal(test,expect)
})

test_that("Test RIE", {
  test<-  RIE(r, w=c(5,7), na.rm = TRUE) |> values(mat=TRUE)
  expect<- readRDS(system.file("testdata", "rie.RDS", package= "MultiscaleDTM"))
  expect_equal(test, expect)
})

test_that("Test RelPos", {
  test<-  RelPos(r, w=matrix(data = c(1,NA,1), nrow = 3, ncol=3), shape = "custom", fun = "median", na.rm = TRUE) |> values(mat=TRUE)
  expect<- readRDS(system.file("testdata", "rp.RDS", package= "MultiscaleDTM"))
  expect_equal(test, expect)
})

test_that("Test TPI", {
  test<-  TPI(r, w=c(5,5), shape= "rectangle", na.rm = TRUE) |> values(mat=TRUE)
  expect<- readRDS(system.file("testdata", "tpi.RDS", package= "MultiscaleDTM"))
  expect_equal(test, expect)
})

test_that("Test DMV", {
  test<- DMV(r, w=2, shape= "circle", na.rm = TRUE, stand="range") |> values(mat=TRUE)
  expect<- readRDS(system.file("testdata", "dmv.RDS", package= "MultiscaleDTM"))
  expect_equal(test, expect)
})

test_that("Test BPI", {
  test<- BPI(r, w = c(4,6), unit = "cell", stand= "sd", na.rm = TRUE) |> values(mat=TRUE)
  expect<- readRDS(system.file("testdata", "bpi.RDS", package= "MultiscaleDTM"))
  expect_equal(test, expect)
})

