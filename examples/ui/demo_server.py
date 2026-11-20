#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Optional loopback-only reference UI. No external packages; not a hosted service."""
import argparse
from http.server import BaseHTTPRequestHandler, HTTPServer
import json
from pathlib import Path
import secrets
import sys
from paacs_bridge import Bridge, BridgeError, MAX_LINE, parse_json, response


def preview(bridge, job):
    """Bounded display geometry. The solver mesh is neither copied nor modified.

    50k source nodes / 200k scanned cells / 4k surface faces. Above the node
    budget ask a real frontend to use the paged C API or its own LOD pipeline.
    """
    summary = bridge.summary(job)
    nn, nc = int(summary['nodes']), int(summary['cells'])
    if nn > 50_000:
        return dict(positions=[], faces=[], message='Preview node budget exceeded; use paged mesh access.', truncated=True)
    positions = []
    for offset in range(0, nn, 4096):
        data, count = bridge._page(job, str(offset), min(4096, nn-offset), False)
        positions.extend(list(data[i].position) for i in range(count))
    faces = []
    truncated = nc > 200_000
    for offset in range(0, min(nc, 200_000), 4096):
        data, count = bridge._page(job, str(offset), min(4096, min(nc, 200_000)-offset), True)
        for i in range(count):
            v = data[i]
            if v.type not in (2, 3, 4):
                continue
            if len(faces) == 4000:
                truncated = True
                break
            # Tri6 is displayed by its corners; result reads keep all six nodes.
            corners = 4 if v.type == 4 else 3
            faces.append(dict(nodes=[int(v.nodes[k]) for k in range(corners)],
                              dimension=v.dimension, tag=v.tag, cell=str(offset+i)))
        if len(faces) == 4000:
            truncated = truncated or offset+i+1 < nc
            break
    return dict(positions=positions, faces=faces, truncated=truncated,
                message='Surface preview is capped at 4,000 faces.' if truncated else 'Full emitted surface; drag to rotate, click to inspect a face tag.')


def make_server(bridge, port=0):
    token = secrets.token_urlsafe(32)
    html = Path(__file__).with_name('index.html').read_text().replace('__PAACS_TOKEN__', token).encode()

    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *_):
            pass

        def send(self, status, body, mime):
            self.send_response(status)
            self.send_header('Content-Type', mime)
            self.send_header('Content-Length', str(len(body)))
            self.send_header('Cache-Control', 'no-store')
            self.send_header('X-Content-Type-Options', 'nosniff')
            self.send_header('Content-Security-Policy', "default-src 'self'; script-src 'self' 'unsafe-inline'; style-src 'self' 'unsafe-inline'; frame-ancestors 'none'; connect-src 'self'")
            self.end_headers()
            self.wfile.write(body)

        def valid_host(self):
            return self.headers.get('Host') == f'127.0.0.1:{self.server.server_port}'

        def do_GET(self):
            if not self.valid_host():
                return self.send(403, b'Invalid host', 'text/plain')
            if self.path != '/':
                return self.send(404, b'Not found', 'text/plain')
            self.send(200, html, 'text/html; charset=utf-8')

        def do_POST(self):
            origin = f'http://127.0.0.1:{self.server.server_port}'
            if (not self.valid_host() or self.headers.get('X-PAACS-Token') != token or
                    self.headers.get('Origin', origin) != origin):
                return self.send(403, b'Forbidden', 'text/plain')
            if self.path != '/api':
                return self.send(404, b'Not found', 'text/plain')
            try:
                count = int(self.headers.get('Content-Length', '-1'))
                if not 0 <= count <= MAX_LINE:
                    raise BridgeError('request exceeds limit')
                payload = self.rfile.read(count)
                command = parse_json(payload)
                if not isinstance(command, dict):
                    raise BridgeError('expected command object')
                # HTTP demo never accepts arbitrary filesystem writes.
                if command.get('method') == 'export':
                    raise BridgeError('use the native/stdio bridge for archive export')
                if command.get('method') == 'preview':
                    from paacs_bridge import object_keys
                    object_keys(command, ('version', 'id', 'method', 'params'), ('version', 'method', 'params'))
                    if type(command['version']) is not int or command['version'] != 1:
                        raise BridgeError('unsupported protocol version')
                    object_keys(command['params'], ('job',), ('job',))
                    out = dict(version=1, ok=True, result=preview(bridge, command['params']['job']))
                else:
                    out = response(bridge, payload)
                self.send(200, json.dumps(out, allow_nan=False).encode(), 'application/json')
            except (ValueError, TypeError, OverflowError, RecursionError) as exc:
                out = dict(ok=False, error=dict(message=str(exc)[:512]))
                self.send(400, json.dumps(out).encode(), 'application/json')

    class LocalServer(HTTPServer):
        # A partial HTTP request cannot occupy the single control thread forever.
        def get_request(self):
            sock, address = super().get_request()
            sock.settimeout(5)
            return sock, address

    return LocalServer(('127.0.0.1', port), Handler)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--library', required=True)
    parser.add_argument('--port', type=int, default=0)
    args = parser.parse_args()
    with Bridge(args.library) as bridge:
        with make_server(bridge, args.port) as server:
            print(f'PAACS reference UI: http://127.0.0.1:{server.server_port}', flush=True)
            print('Keep this terminal open. Ctrl+C cancels work and exits.', file=sys.stderr)
            try:
                server.serve_forever()
            except KeyboardInterrupt:
                pass


if __name__ == '__main__':
    main()
