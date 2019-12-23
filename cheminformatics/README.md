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
conda create -c rdkit conda-forge -n tas_cheminformatics rdkit python=3.7
conda install click flask pandas gunicorn marshmallow apispec molvs
```

### Installing the TAS Cheminformatics tools

```bash
pip install --no-deps 'git+https://github.com/labsyspharm/small-molecule-suite-maintenance.git#egg=tas_cheminformatics&subdirectory=cheminformatics'
```

## Running the server

The JSON API

```bash
conda activate tas_cheminformatics
gunicorn --workers=4 -b 127.0.0.1:8000 -t 600 tas_cheminformatics
```
