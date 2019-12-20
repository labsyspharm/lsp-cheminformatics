import sys

from flask import Flask


app = Flask(__name__)
# To enable flask to run in directly
application = app

@app.route("/doc")
def doc_route():
    return app.send_static_file("apidoc.html")

if sys.version_info[0] >= 3:
    import tas_cheminformatics.tautomer_server
    import tas_cheminformatics.properties_server
    import tas_cheminformatics.draw_server
import tas_cheminformatics.fingerprint_server
