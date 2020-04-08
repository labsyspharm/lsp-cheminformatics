test_that("morgan similarities are returned", {
  sim <- chemical_similarity(
    c("resveratrol" = "InChI=1S/C14H12O3/c15-12-5-3-10(4-6-12)1-2-11-7-13(16)9-14(17)8-11/h1-9,15-17H/b2-1+"),
    c(
      "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
      "aspirin" = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
    )
  )
  expect_equal(nrow(sim), 2L)
})

test_that("substructure matches are returned", {
  sim <- match_substructure(
    c("secondary_amine" = "[H]N(C)C"),
    c(
      "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
      "aspirin" = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
    ),
    query_identifier = "smiles"
  )
  expect_equal(nrow(sim), 1L)
  expect_equal(length(sim$match[[1]]), 6)
})

test_that("morgan similarities with threshold are returned", {
  sim <- chemical_similarity_threshold(
    c("resveratrol" = "InChI=1S/C14H12O3/c15-12-5-3-10(4-6-12)1-2-11-7-13(16)9-14(17)8-11/h1-9,15-17H/b2-1+"),
    c(
      "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
      "aspirin" = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
    ),
    threshold = 0.1
  )
  expect_equal(nrow(sim), 1L)
  expect_equal(sim[["query"]], "resveratrol")
  expect_equal(sim[["target"]], "aspirin")
})

test_that("morgan similarities with threshold using precalculated fingerprints are returned", {
  sim <- chemical_similarity_threshold(
    c("a" = "880DF", "b" = "881DF", "c" = "870E0"),
    threshold = 0.8,
    precalculated = TRUE
  )
  expect_equal(nrow(sim), 1L)
  expect_equal(sim[["query"]], "a")
  expect_equal(sim[["target"]], "b")
})

test_that("identical compounds are found", {
  sim <- compound_identity(
    c("resveratrol" = "InChI=1S/C14H12O3/c15-12-5-3-10(4-6-12)1-2-11-7-13(16)9-14(17)8-11/h1-9,15-17H/b2-1+"),
    c(
      "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
      "aspirin" = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
      "resveratrol2" = "InChI=1S/C14H12O3/c15-12-5-3-10(4-6-12)1-2-11-7-13(16)9-14(17)8-11/h1-9,15-17H/b2-1+"
    )
  )
  expect_equal(nrow(sim), 1L)
  expect_equal(sim[["query"]], "resveratrol")
  expect_equal(sim[["target"]], "resveratrol2")
})

test_that("chemical fingerprints are returned", {
  fp <- calculate_fingerprints(
    c(
      "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
      "aspirin" = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
    )
  )
  expect_equal(nrow(fp), 2L)
})
