import sys

from flask import Flask


app = Flask(__name__)
# To enable flask to run in directly
application = app

@app.route("/doc")
def doc_route():
    return app.send_static_file("apidoc.html")

if sys.version_info[0] >= 3:
    import tas_chemoinformatics.tautomer_server
    import tas_chemoinformatics.convert_server
    import tas_chemoinformatics.draw_server
import tas_chemoinformatics.fingerprint_server
