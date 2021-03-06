test_that("inchis are converted in smiles", {
  sim <- convert_compound_descriptor(
    c(
      "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
      "aspirin" = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
    ),
    target_descriptor = "smiles"
  )
  expect_equal(nrow(sim), 2L)
})

test_that("molecular masses are calculated", {
  mass <- molecular_mass(
    c(
      "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
      "aspirin" = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
    )
  )
  expect_equal(nrow(mass), 2L)
})

test_that("organic compounds are detected", {
  res <- is_organic(
    c(
      "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
      "nitrous oxide" = "InChI=1S/N2O/c1-2-3"
    )
  )
  expect_equal(nrow(res), 2L)
  res_map <- magrittr::set_names(
    res[["is_organic"]],
    res[["compound"]]
  )
  expect_false(res_map[["nitrous oxide"]])
  expect_true(res_map[["tofacitnib"]])
})
