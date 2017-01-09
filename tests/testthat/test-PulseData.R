context("PulseData creation and normalisation")

test_that("normalisation works", {
  pd <- new.env()
  pd$user_conditions <- data.frame(conditions=rep(c("a","b"), each=3))
  pd$count_data <- matrix(rep(c(1,2,4),2),ncol=6, nrow=20, byrow=TRUE)
  colnames(pd$count_data) <- rownames(pd$user_conditions)
  rownames(pd$count_data) <- 1:20
  class(pd) <- "PulseData"    
  normalise(pd)
  expect_equivalent(pd$norm_factors, rep(c(.5,1,2),2))
})
