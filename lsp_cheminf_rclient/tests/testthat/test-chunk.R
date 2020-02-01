test_that("dataframe is chunked", {
  df <- data.frame(x = 1:10)
  chunked <- chunk_df(df, 3)
  expect_length(chunked, 3)
  expect_equal(sum(as.integer(lapply(chunked, nrow))), 10)
})
