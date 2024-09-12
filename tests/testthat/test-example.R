# run example from function documentation
example(spatial_simulate, echo = FALSE)

test_that("example object has correct class", {
  expect_s4_class(spe, "SpatialExperiment")
})

test_that("example objects have correct dimensions", {
  expect_equal(dim(spe), c(1, 748))
})
