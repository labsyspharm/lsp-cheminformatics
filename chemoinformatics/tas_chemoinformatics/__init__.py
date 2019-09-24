from flask import Flask

app = Flask(__name__)

import tas_chemoinformatics.tautomer_server
import tas_chemoinformatics.fingerprints
