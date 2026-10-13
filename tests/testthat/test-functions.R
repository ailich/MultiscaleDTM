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

test_that("Test Qfit queen", {
  test<- Qfit(r, w = c(5,7), unit = "degrees", method = "queen", metrics = c("elev", "qslope", "qaspect", "qeastness", "qnorthness", "profc", "planc", "twistc", "meanc", "maxc", "minc", "features"), na.rm = TRUE) |> values(mat=TRUE)
  expect<- readRDS(system.file("testdata", "qmetrics_queen.RDS", package= "MultiscaleDTM"))
  expect_equal(test, expect)
})

test_that("Test Qfit boundary", {
  test<- Qfit(r, w = c(5,7), unit = "degrees", method = "boundary", metrics = c("elev", "qslope", "qaspect", "qeastness", "qnorthness", "profc", "planc", "twistc", "meanc", "maxc", "minc", "features"), na.rm = TRUE) |> values(mat=TRUE)
  expect<- readRDS(system.file("testdata", "qmetrics_boundary.RDS", package= "MultiscaleDTM"))
  expect_equal(test, expect)
})

test_that("Qfit methods agree on a 3x3 window", {
  # Every cell of a 3x3 window is on its outer edge, so all methods select all 9 cells
  metrics<- c("elev", "qslope", "qaspect", "qeastness", "qnorthness", "profc", "planc", "twistc", "meanc", "maxc", "minc", "features")
  full<- Qfit(r, w = c(3,3), metrics = metrics, na.rm = TRUE, method = "full") |> values(mat=TRUE)
  expect_equal(Qfit(r, w = c(3,3), metrics = metrics, na.rm = TRUE, method = "queen") |> values(mat=TRUE), full)
  expect_equal(Qfit(r, w = c(3,3), metrics = metrics, na.rm = TRUE, method = "boundary") |> values(mat=TRUE), full)
})

test_that("Qfit methods agree on a 3x3 window with force_center", {
  # The central cell must stay among the selected cells even though it does not
  # affect the fit, otherwise the flat window shortcut in the C++ kernels zeroes
  # the coefficients wherever the outer edge is flat but the centre differs
  # (e.g. a single cell pit), erasing the feature force_center exists to detect.
  full<- Qfit(r, w = c(3,3), method = "full", force_center = TRUE, na.rm = TRUE, metrics = c(), return_params = TRUE, outlier_quantile = c(0,1)) |> values(mat=TRUE)
  expect_equal(Qfit(r, w = c(3,3), method = "queen", force_center = TRUE, na.rm = TRUE, metrics = c(), return_params = TRUE, outlier_quantile = c(0,1)) |> values(mat=TRUE), full)
  expect_equal(Qfit(r, w = c(3,3), method = "boundary", force_center = TRUE, na.rm = TRUE, metrics = c(), return_params = TRUE, outlier_quantile = c(0,1)) |> values(mat=TRUE), full)

  # A flat plateau with a one cell pit at its centre: curvature must be detected
  pit_r<- terra::rast(nrows=11, ncols=11, xmin=0, xmax=11, ymin=0, ymax=11, crs="EPSG:32610")
  terra::values(pit_r)<- 100
  pit_r[6,6]<- 50
  pit_cell<- terra::cellFromRowCol(pit_r, 6, 6)
  for(m in c("full", "queen", "boundary")){
    curv<- Qfit(pit_r, w = c(3,3), method = m, force_center = TRUE, metrics = "maxc", outlier_quantile = c(0,1))
    expect_equal(terra::extract(curv, pit_cell)[1,1], -60)
  }
})

test_that("Qfit recovers an exact quadratic surface under every method", {
  qr_rast<- terra::rast(nrows=40, ncols=40, xmin=0, xmax=40, ymin=0, ymax=40, crs="EPSG:32610")
  xy<- terra::xyFromCell(qr_rast, 1:terra::ncell(qr_rast))
  A<- 2; B<- -3; C<- 1.5
  terra::values(qr_rast)<- A*xy[,1]^2 + B*xy[,2]^2 + C*xy[,1]*xy[,2] + 4*xy[,1] - 5*xy[,2] + 10

  full<- Qfit(qr_rast, w = c(7,7), method = "full", metrics = c(), return_params = TRUE, outlier_quantile = c(0,1)) |> values(mat=TRUE)
  for(m in c("queen", "boundary")){
    # An exact quadratic is fit exactly by any full rank subset of points
    expect_equal(Qfit(qr_rast, w = c(7,7), method = m, metrics = c(), return_params = TRUE, outlier_quantile = c(0,1)) |> values(mat=TRUE), full)
  }
  keep<- stats::complete.cases(full)
  expect_equal(unname(full[keep, "a"]), rep(A, sum(keep)))
  expect_equal(unname(full[keep, "b"]), rep(B, sum(keep)))
  expect_equal(unname(full[keep, "c"]), rep(C, sum(keep)))
})

test_that("Qfit rejects invalid methods and rank deficient windows", {
  expect_error(Qfit(r, w = c(5,7), method = "rook", metrics = "qslope"), "rook")
  expect_error(Qfit(r, w = c(5,7), method = "bogus", metrics = "qslope"), "must be")
  expect_error(Qfit(r, w = c(1,7), method = "full", metrics = "qslope"), "quadratic surface")
})

test_that("Test equivalence of Pfit and Qfit", {
  test<- Pfit(r, w = c(5,7), unit = "degrees", metrics = c("pslope", "paspect", "peastness", "pnorthness"), na.rm = FALSE) |> values(mat=TRUE) |> `dimnames<-`(NULL)
  expect<- Qfit(r, w = c(5,7), unit = "degrees", metrics = c("qslope", "qaspect", "qeastness", "qnorthness"), na.rm = FALSE, outlier_quantile = c(0,1)) |> values(mat=TRUE) |> `dimnames<-`(NULL)
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

