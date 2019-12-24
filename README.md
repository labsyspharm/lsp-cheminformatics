# LSP Cheminformatics Tools

Set of chemoinformatics tools to find canonical representations of compounds
and fingerprint them. Implemented as python package with commandline interface
and webserver with JSON API.

## Installation

### Installing Anaconda

We recommend installing this package in a separate Anaconda environment. If Ananconda is not
available, first [install Anaconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

### Creating the environment and installing dependencies

Once anaconda is installed, create a new environment for this package and install RDKit
([see here for more help](https://www.rdkit.org/docs/Install.html)).

``` bash
conda create -c rdkit -n lspcheminf_env python=3.7 rdkit click flask pandas gunicorn marshmallow apispec
conda activate lspcheminf_env
conda install -c conda-forge molvs
```

### Installing the LSP Cheminformatics tools

``` bash
conda activate lspcheminf_env
pip install --no-deps 'git+https://github.com/labsyspharm/lsp-cheminformatics.git#egg=lspcheminf&subdirectory=lsp_cheminf_server'
```

## Running the server

The JSON API is exposed by running the server in the background. The following command runs the server on port 8000.

``` bash
conda activate lspcheminf_env
gunicorn --workers=4 -b 127.0.0.1:8000 -t 600 lspcheminf
```

Once the server is running the documentation for the JSON API is available at [http://127.0.0.1:8000/doc](http://127.0.0.1:8000/doc).

## Querying the server

Any software capable of sending and receiving JSON can be used to query the server. For convenience, there is an R client
that implements the JSON API in simple functions. In this example we use the R client to query for the chemical similarities between a number of compounds.

### Installing the R client

``` r
devtools::install_github("labsyspharm/lsp-cheminformatics", subdir = "lsp_cheminf_rclient")
```

### Query for similarities

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

### Query for substructure matches

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

At the moment there is no R client code yet for the additional functionalities of the JSON API, like the molecule drawing,
ID conversion etc.

## Funding

This work was supported by NIH grants U54-HL127365, U24-DK116204 and U54-HL127624.
