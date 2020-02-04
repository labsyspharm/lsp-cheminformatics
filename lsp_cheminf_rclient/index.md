# LSP Cheminformatics R client

The R client provides easy access to the functionality provided by the LSP Cheminformatics server.

## Installing the R client

``` r
devtools::install_github("labsyspharm/lsp-cheminformatics", subdir = "lsp_cheminf_rclient")
```

The client only functions while the lspcheminf server is running.

## Query for similarities

The R client can be used to query the chemical similarity between compounds. By default the Morgan fingerprinting algorithm is used.

``` r
library(lspcheminf)
chemical_similarity(
  c("resveratrol" = "InChI=1S/C14H12O3/c15-12-5-3-10(4-6-12)1-2-11-7-13(16)9-14(17)8-11/h1-9,15-17H/b2-1+"),
  c(
    "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
    "aspirin" = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
  )
)
#> # A tibble: 2 x 3
#>   query        score target
#>   <chr>        <dbl> <chr>
#> 1 resveratrol 0.0571 tofacitnib
#> 2 resveratrol 0.132  aspirin
```

If we want to compute the similarity based on topological fingerprints instead and pass some custom arguments to the RDKit function
we can do this:

  ``` r
chemical_similarity(
  c("resveratrol" = "InChI=1S/C14H12O3/c15-12-5-3-10(4-6-12)1-2-11-7-13(16)9-14(17)8-11/h1-9,15-17H/b2-1+"),
  c(
    "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
    "aspirin" = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
  ),
  fingerprint_type = "topological",
  fingerprint_args = list("minPath" =  2, "useHs" = FALSE)
)
#> # A tibble: 2 x 3
#>   query        score target
#>   <chr>        <dbl> <chr>
#> 1 resveratrol 0.0967 tofacitnib
#> 2 resveratrol 0.198  aspirin
```

## Query for substructure matches

``` r
match_substructure(
  c("secondary_amine" = "[H]N(C)C"),
  c(
    "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
    "aspirin" = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
  ),
  query_identifier = "smiles"
)
#> # A tibble: 1 x 3
#>   match      query           target
#>   <list>     <chr>           <chr>
#> 1 <list [6]> secondary_amine tofacitnib
```

Six matches for secondary amine groups where found in tofacitnib. The atom indices for each match are stored in the `match`
list column.

## Additional functionality

All additional functionality is documented in the `Reference`
