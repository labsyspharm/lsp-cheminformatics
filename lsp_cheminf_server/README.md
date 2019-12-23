# TAS Cheminformatics Tools

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

```bash
conda create -c rdkit -n tas_cheminformatics python=3.7 rdkit click flask pandas gunicorn marshmallow apispec
conda activate tas_cheminformatics
conda install -c conda-forge molvs
```

### Installing the TAS Cheminformatics tools

```bash
conda activate tas_cheminformatics
pip install --no-deps 'git+https://github.com/labsyspharm/small-molecule-suite-maintenance.git#egg=tas_cheminformatics&subdirectory=cheminformatics'
```

## Running the server

The JSON API is exposed by running the server in the background. The following command runs the server on port 8000.

```bash
conda activate tas_cheminformatics
gunicorn --workers=4 -b 127.0.0.1:8000 -t 600 tas_cheminformatics
```

Once the server is running the documentation for the JSON API is available at [http://127.0.0.1:8000/doc](http://127.0.0.1:8000/doc).

## Querying the server

In this example we use the JSON API from R to query for the chemical similarities between a number of compounds.

### Installing the R client

```r
devtools::install_github("labsyspharm/small-molecule-suite-maintanance", subdir = "")
```
