import sys

from flask import Flask


app = Flask(__name__)
# To enable flask to run in directly
application = app


@app.route("/doc")
def doc_route():
    return app.send_static_file("apidoc.html")


import lspcheminf.tautomer_server
import lspcheminf.properties_server
import lspcheminf.draw_server
import lspcheminf.fingerprint_server
