# TAS Chemoinformatics Tools

Set of chemoinformatics tools to find canonical representations of compounds
and fingerprint them. Implemented as python package with commandline interface
and webserver with JSON API.

## Installation

TODO

## Using JSON API

The JSON API

```bash
gunicorn --workers=4 -b 127.0.0.1:8000 -t 600 tas_chemoinformatics
```
