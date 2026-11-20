#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""JSON/ctypes contract tests. Standard library only; no browser/network needed."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'examples/ui'))
from paacs_bridge import Bridge, BridgeError, MAX_LINE, response
LIBRARY = str(Path(sys.argv.pop(1)).resolve())


class ProtocolTests(unittest.TestCase):
    def setUp(self):
        self.b = Bridge(LIBRARY, maximum_result_bytes=64*1024**2, maximum_retained_bytes=128*1024**2, compute_hash=True)

    def tearDown(self):
        self.b.close()

    def run_job(self, spec):
        j = self.b.submit(spec)['job']
        self.assertTrue(self.b.wait(j, 20_000))
        self.assertEqual(self.b.status(j)['state'], 'succeeded', self.b.status(j))
        return j

    def test_schema(self):
        for route, count in [('primitive', 49), ('adaptive_box', 19)]:
            s = self.b.schema(route)
            self.assertEqual(s['abi'], 1)
            self.assertEqual(len(s['options']), count)
            self.assertEqual(len({x['key'] for x in s['options']}), count)

    def test_all_primitive_families_and_elements(self):
        for kind in ('cylinder', 'cone', 'frustum', 'tube', 'torus', 'elbow'):
            for elements in ('hex_dominant', 'tet4', 'tet10'):
                options = dict(kind=kind, elements=elements, target_size=.6,
                               maximum_chord_error=.03, maximum_relative_volume_error=.03)
                if kind == 'tube':
                    options['inner_radius'] = .5
                j = self.run_job(dict(route='primitive', options=options))
                s = self.b.summary(j)
                self.assertGreater(int(s['nodes']), 0)
                self.assertGreater(s['minimum_quality'], 0)
                self.b.release(j)

    def test_profile_field_and_export(self):
        j = self.run_job(dict(route='primitive', options=dict(height=1, profile_global_size=.4,
            target_size=.4, elements='tet4', maximum_chord_error=.04),
            profile_controls=[dict(begin=0, end=0, size=.2, growth=1.3)]))
        with tempfile.TemporaryDirectory() as d:
            p = str(Path(d)/'profile.paacsmesh')
            self.b.export(j, p)
            self.assertGreater(Path(p).stat().st_size, 100)
        self.b.release(j)
        j = self.run_job(dict(route='adaptive_box', options=dict(target_size=.6), field=dict(
            nodes=[dict(operation='distance_ball', a=[.5,.5,.5], parameters=[.1,0,0,0]),
                   dict(operation='threshold_linear', first=0, parameters=[0,.35,.2,.6])], output=1)))
        self.assertGreater(int(self.b.summary(j)['cell_types'][5]), 0)
        self.b.release(j)

    def test_exact_ids_and_pages(self):
        j = self.run_job(dict(route='adaptive_box', options=dict(target_size=.6)))
        s = self.b.summary(j)
        self.assertIsInstance(j, str)
        self.assertEqual(len(s['hash']), 16)
        first = self.b.page(j, '0', 17, cells=True)
        self.assertEqual(first['next_offset'], '17')
        self.assertEqual([x['index'] for x in first['values']], [str(x) for x in range(17)])
        self.assertTrue(all(isinstance(n, str) for c in first['values'] for n in c['nodes']))
        self.assertEqual(self.b.page(j, s['nodes'], 10)['values'], [])
        for offset, count in [('18446744073709551616',1), ('-1',1), ('1.2',1), ('0',4097), ('0',-1), (0,1)]:
            with self.assertRaises(BridgeError):
                self.b.page(j, offset, count)
        self.b.release(j)
        with self.assertRaises(BridgeError):
            self.b.status(j)

    def test_explicit_profile_disable(self):
        spec = dict(route='primitive', options=dict(target_size=.5, maximum_chord_error=.03))
        a = self.run_job(spec)
        spec['options']['profile_enabled'] = False
        # This profile length does not match height. Disabled controls must not apply.
        spec['profile_controls'] = [dict(begin=0, end=0, size=.2, growth=1.3)]
        b = self.run_job(spec)
        self.assertEqual(self.b.summary(a)['hash'], self.b.summary(b)['hash'])
        self.b.release(a)
        self.b.release(b)

    def test_unknown_options_and_nonfinite_inputs(self):
        for options in [dict(typo=1), dict(target_size=float('nan')), dict(target_size=True),
                        dict(threads=.5), dict(elements='invalid'), dict(radius=10**1000), dict(threads=-1)]:
            with self.assertRaises(BridgeError):
                self.b.submit(dict(route='primitive', options=options))
        with self.assertRaises(BridgeError):
            self.b.submit(dict(route='primitive', field=None))
        with self.assertRaises(BridgeError):
            self.b.submit(dict(route='adaptive_box', profile_controls=[]))
        with self.assertRaises(BridgeError):
            self.b.submit(dict(route='adaptive_box', field=dict(nodes=[dict(operation='scale', first=5)], output=0)))
        self.assertEqual(len(self.b.jobs), 0)

    def test_failure_preserves_previous_result(self):
        first = self.run_job(dict(route='primitive', options=dict(target_size=.5, maximum_chord_error=.03)))
        old = self.b.summary(first)
        bad = self.b.submit(dict(route='primitive', options=dict(radius=0)))['job']
        self.assertTrue(self.b.wait(bad, 10_000))
        self.assertEqual(self.b.status(bad)['state'], 'failed')
        with self.assertRaises(BridgeError):
            self.b.summary(bad)
        self.assertEqual(self.b.summary(first), old)
        self.b.release(bad)
        self.b.release(first)

    def test_cancellation(self):
        for _ in range(10):
            j = self.b.submit(dict(route='primitive', options=dict(elements='tet4', target_size=.1)))['job']
            accepted = self.b.cancel(j)['accepted']
            self.assertTrue(self.b.wait(j, 10_000))
            if accepted:
                self.assertEqual(self.b.status(j)['state'], 'cancelled')
                with self.assertRaises(BridgeError):
                    self.b.summary(j)
            self.assertFalse(self.b.cancel(j)['accepted'])
            self.b.release(j)

    def test_request_schema_refusals(self):
        for text in [b'{bad', b'\xff', '[]', '{"version":1,"version":1,"method":"schema"}',
                     '{"version":NaN}', '{"version":1,"method":[]}',
                     '{"version":true,"method":"schema","params":{"route":"primitive"}}',
                     '{"version":1,"method":"status","params":{"job":1}}',
                     '{"version":1,"method":"schema","params":{"route":"primitive","extra":1}}']:
            self.assertFalse(response(self.b, text)['ok'])
        valid = response(self.b, '{"version":1,"id":"hello","method":"schema","params":{"route":"primitive"}}')
        self.assertTrue(valid['ok'])
        self.assertEqual(valid['id'], 'hello')

    def test_stdio_recovers_after_oversize_and_bad_json(self):
        schema = json.dumps(dict(version=1, id='recover', method='schema', params=dict(route='adaptive_box')))
        payload = b'x'*(MAX_LINE+17)+b'\n'+b'{broken\n'+schema.encode()+b'\n'
        p = subprocess.run([sys.executable, str(ROOT/'examples/ui/paacs_bridge.py'), '--library', LIBRARY],
                           input=payload, capture_output=True, timeout=30)
        self.assertEqual(p.returncode, 0, p.stderr.decode())
        lines = [json.loads(x) for x in p.stdout.splitlines()]
        self.assertEqual(len(lines), 3)
        self.assertEqual(lines[0]['error']['code'], 'request_limit')
        self.assertFalse(lines[1]['ok'])
        self.assertTrue(lines[2]['ok'])
        self.assertEqual(lines[2]['id'], 'recover')


if __name__ == '__main__':
    unittest.main(verbosity=2)
