"""
Local server for the step-by-step visualizer.

Run from anywhere:   python web/server.py      then open http://localhost:8000

GET /                  serves web/index.html
GET /api/trace?read=N&mutations=..&coverage=..&error_rate=..&seed=..&min_depth=..&max_edit=..
                       runs the C trace program for read N (-1 = auto pick) with
                       those settings (each optional, range-checked) and returns its JSON

The C program (trace.exe) is built with gcc the first time, and again
whenever a .c or .h file under src/ is newer than it. Python standard
library only; nothing to install.
"""

import http.server
import json
import os
import subprocess
import sys
import urllib.parse

WEB_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(WEB_DIR)
EXE = os.path.join(ROOT, "trace.exe" if os.name == "nt" else "trace")
PORT = int(os.environ.get("PORT", "8000"))   # set PORT to use another port

STAGE_DIRS = ["stage1_load", "stage2_minimizers", "stage3_index", "stage4_bloom",
              "stage5_vote", "stage6_top3", "stage7_verify", "stage8_pileup",
              "stage9_score", "simulate"]
# Settings the page may change: name -> (type, default, lowest, highest), in trace.c's argument order.
SETTINGS = [
    ("mutations", int, 20, 1, 500),
    ("coverage", float, 30.0, 1.0, 100.0),
    ("error_rate", float, 0.01, 0.0, 0.10),
    ("seed", int, 42, 0, 1_000_000_000),
    ("min_depth", int, 3, 1, 50),
    ("max_edit", int, 3, 0, 10),
]

SOURCES = ["stage1_load/seq.c", "stage2_minimizers/minimizer.c", "stage3_index/index.c",
           "stage4_bloom/bloom.c", "stage5_vote/vote.c", "stage6_top3/top3.c",
           "stage7_verify/verify.c", "stage8_pileup/pileup.c", "stage9_score/score.c",
           "simulate/simulate.c", "demo/trace.c"]


def newest_source_time():
    """Latest modification time of any .c or .h file under src/."""
    newest = 0.0
    for folder, _, files in os.walk(os.path.join(ROOT, "src")):
        for name in files:
            if name.endswith((".c", ".h")):
                newest = max(newest, os.path.getmtime(os.path.join(folder, name)))
    return newest


def build_if_needed():
    """Compiles trace.exe when it is missing or older than the sources. Returns an error string or None."""
    if os.path.exists(EXE) and os.path.getmtime(EXE) >= newest_source_time():
        return None
    cmd = ["gcc", "-std=c99", "-O2", "-Wall", "-Wextra"]
    cmd += ["-I" + os.path.join(ROOT, "src", d) for d in STAGE_DIRS]
    cmd += [os.path.join(ROOT, "src", s) for s in SOURCES]
    cmd += ["-lm", "-o", EXE]
    print("building trace program ...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return "gcc failed:\n" + result.stderr
    return None


class Handler(http.server.SimpleHTTPRequestHandler):
    """Serves web/ as static files plus the /api/trace endpoint."""

    def __init__(self, *args, **kwargs):
        """Serve files from the web folder regardless of where the server was started."""
        super().__init__(*args, directory=WEB_DIR, **kwargs)

    def do_GET(self):
        """Routes /api/trace to the C program; everything else is a static file."""
        url = urllib.parse.urlparse(self.path)
        if url.path != "/api/trace":
            return super().do_GET()

        query = urllib.parse.parse_qs(url.query)
        try:
            read = int(query.get("read", ["-1"])[0])
        except ValueError:
            return self.send_json(400, {"error": "read must be a whole number"})

        args = []
        for name, kind, default, low, high in SETTINGS:
            raw = query.get(name, [str(default)])[0]
            try:
                value = kind(raw)
            except ValueError:
                return self.send_json(400, {"error": f"{name} must be a number"})
            if not low <= value <= high:
                return self.send_json(400, {"error": f"{name} must be between {low} and {high}"})
            args.append(str(value))

        error = build_if_needed()
        if error:
            return self.send_json(500, {"error": error})

        result = subprocess.run([EXE, os.path.join(ROOT, "sequence.fasta"), str(read)] + args,
                                capture_output=True, text=True, cwd=ROOT)
        if result.returncode != 0:
            return self.send_json(500, {"error": result.stdout or result.stderr})
        body = result.stdout.encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def send_json(self, status, obj):
        """Sends a small JSON response with the given HTTP status."""
        body = json.dumps(obj).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


def main():
    """Builds the C program once, then serves until Ctrl+C."""
    error = build_if_needed()
    if error:
        print(error)
        sys.exit(1)
    server = http.server.ThreadingHTTPServer(("127.0.0.1", PORT), Handler)
    print(f"visualizer running at http://localhost:{PORT}  (Ctrl+C to stop)")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
