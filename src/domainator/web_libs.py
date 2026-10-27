"""The JavaScript libraries the HTML reports load from a CDN, at pinned versions.

Each tag carries a Subresource Integrity hash, so a browser refuses a file that differs from
the one these hashes were taken from. The reports that use them need a network connection to
show their charts; the pages built by report_html.py and the SSN viewer need none.

To change a version, change the URL and recompute its hash:
    curl -sSf URL | openssl dgst -sha384 -binary | openssl base64 -A
"""

def _script(url: str, integrity: str) -> str:
    return f'<script src="{url}" integrity="{integrity}" crossorigin="anonymous"></script>'


D3 = _script("https://cdn.jsdelivr.net/npm/d3@7.9.0/dist/d3.min.js",
             "sha384-CjloA8y00+1SDAUkjs099PVfnY2KmDC2BZnws9kh8D/lX1s46w6EPhpXdqMfjK6i")
PLOTLY = _script("https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.27.0/plotly.min.js",
                 "sha384-Hl48Kq2HifOWdXEjMsKo6qxqvRLTYqIGbvlENBmkHAxZKIGCXv43H6W1jA671RzC")
