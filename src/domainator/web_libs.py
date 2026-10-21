"""The JavaScript libraries the HTML reports load from a CDN, at pinned versions.

Each tag carries a Subresource Integrity hash, so a browser refuses a file that differs from
the one these hashes were taken from. The reports need a network connection to show their
charts; the SSN viewer (ssn_viewer_html.py) needs none.

To change a version, change the URL and recompute its hash:
    curl -sSf URL | openssl dgst -sha384 -binary | openssl base64 -A
"""

def _script(url: str, integrity: str) -> str:
    return f'<script src="{url}" integrity="{integrity}" crossorigin="anonymous"></script>'


D3 = _script("https://cdn.jsdelivr.net/npm/d3@7.9.0/dist/d3.min.js",
             "sha384-CjloA8y00+1SDAUkjs099PVfnY2KmDC2BZnws9kh8D/lX1s46w6EPhpXdqMfjK6i")
PLOTLY = _script("https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.27.0/plotly.min.js",
                 "sha384-Hl48Kq2HifOWdXEjMsKo6qxqvRLTYqIGbvlENBmkHAxZKIGCXv43H6W1jA671RzC")
TABULATOR = (
    '<link href="https://cdn.jsdelivr.net/npm/tabulator-tables@5.4.3/dist/css/tabulator.min.css" rel="stylesheet" '
    'integrity="sha384-kG8f0sPG81gBEhEJEjpzgamvakONxaSBS7clMM1kRIVKCe8WiotTGWC9rQrEr7hJ" crossorigin="anonymous">\n'
    + _script("https://cdn.jsdelivr.net/npm/tabulator-tables@5.4.3/dist/js/tabulator.min.js",
              "sha384-RP2rMM8yN1UN1x7O19WtHXv/JiLe9JhuP048uaKBFGgmwhKCxWJ93WUoup94FQGB")
)
