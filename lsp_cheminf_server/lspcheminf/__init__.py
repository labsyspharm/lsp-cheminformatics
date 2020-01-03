import sys
import traceback

from flask import Flask, jsonify


app = Flask(__name__)
# To enable flask to run in directly
application = app


@app.route("/doc")
def doc_route():
    return app.send_static_file("apidoc.html")


@app.errorhandler(500)
def internal_server_error(e):
    return (
        jsonify(
            {
                "error": e.__class__.__name__,
                "traceback": "".join(traceback.format_tb(e.__traceback__)),
            }
        ),
        500,
    )


import lspcheminf.tautomer_server
import lspcheminf.properties_server
import lspcheminf.draw_server
import lspcheminf.fingerprint_server
