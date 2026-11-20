#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Standard-library ctypes adapter + versioned JSON Lines worker.

Run --library /absolute/path/to/libpaacs_ui.so (or .dylib / paacs_ui.dll).
One input line -> one response line. No HTTP listener or dependency is needed.
See docs/UI_INTEGRATION.md for the protocol and native route coverage.
"""
import argparse
import ctypes as C
import json
import math
from pathlib import Path
import sys

ABI = 1
PROTOCOL = 1
MAX_LINE = 1_048_576
MAX_PAGE = 4096
U32, U64, SIZE, PTR = C.c_uint32, C.c_uint64, C.c_size_t, C.c_void_p
DOUBLE3 = C.c_double * 3
STATES = ('queued', 'waiting_for_engine', 'running', 'succeeded', 'failed', 'cancelled')
STAGES = ('queued', 'preparing', 'generating', 'finalizing', 'complete')
FIELD_OPS = ('constant', 'distance_point', 'distance_segment', 'distance_triangle',
             'distance_box', 'distance_ball', 'threshold_linear', 'threshold_smooth',
             'minimum', 'maximum', 'scale', 'clamp')


class Option(C.Structure):
    _fields_ = [('key', C.c_char_p), ('kind', U32), ('default', C.c_double),
                ('minimum', C.c_double), ('maximum', C.c_double), ('choices', C.c_char_p)]


class Snapshot(C.Structure):
    _fields_ = [(x, U64) for x in ('id', 'completed', 'total')] + [
        ('elapsed_seconds', C.c_double)] + [(x, U32) for x in (
            'state', 'stage', 'cancellation_requested', 'error', 'domain', 'engine_reason')] + [
        ('diagnostic', C.c_char * 192)]


class Summary(C.Structure):
    _fields_ = [(x, U64) for x in ('nodes', 'cells', 'connectivity', 'logical_bytes', 'hash')] + [
        ('cell_types', U64 * 10), ('lower', DOUBLE3), ('upper', DOUBLE3),
        ('minimum_quality', C.c_double), ('mesh_volume', C.c_double),
        ('quality_warning', U32), ('hash_available', U32)]


class Node(C.Structure):
    _fields_ = [('position', DOUBLE3), ('dimension', C.c_int32), ('tag', C.c_int32)]


class Cell(C.Structure):
    _fields_ = [('type', U32), ('node_count', U32), ('dimension', C.c_int32),
                ('tag', C.c_int32), ('nodes', U64 * 10)]


class ProfileControl(C.Structure):
    _fields_ = [(x, C.c_double) for x in ('begin', 'end', 'size', 'growth')]


class FieldNode(C.Structure):
    _fields_ = [(x, U32) for x in ('operation', 'first', 'second')] + [
        ('a', DOUBLE3), ('b', DOUBLE3), ('c', DOUBLE3), ('parameters', C.c_double * 4)]


class BridgeError(ValueError):
    def __init__(self, message, code='invalid_request'):
        super().__init__(message)
        self.code = code


def check(status, message='native operation rejected'):
    if status:
        raise BridgeError(f'{message} (status {status})', f'native_{status}')


def object_keys(value, allowed, required=()):
    if not isinstance(value, dict):
        raise BridgeError('expected an object')
    extra, missing = set(value) - set(allowed), set(required) - set(value)
    if extra or missing:
        raise BridgeError(f'unknown keys {sorted(extra)}; missing keys {sorted(missing)}')


def integer(value, minimum=0, maximum=2**32-1):
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise BridgeError(f'expected integer in [{minimum}, {maximum}]')
    return value


def number(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise BridgeError('expected a finite number')
    try:
        v = float(value)
    except OverflowError as exc:
        raise BridgeError('number too large') from exc
    if not math.isfinite(v):
        raise BridgeError('expected a finite number')
    return v


def vector(value, n):
    if not isinstance(value, list) or len(value) != n:
        raise BridgeError(f'expected {n} numbers')
    return [number(x) for x in value]


class Bridge:
    """Thread-confined control object. Native meshing happens on its own worker.

    Use close() / a context manager. Jobs/results are explicitly released;
    no finalizers call into a DLL during interpreter teardown.
    """
    def __init__(self, library, maximum_jobs=8, maximum_result_bytes=512*1024**2,
                 maximum_retained_bytes=512*1024**2, compute_hash=False):
        path = Path(library).resolve(strict=True)
        self.lib = C.CDLL(str(path))
        self.jobs, self.results = {}, {}
        self.service = None
        if not isinstance(compute_hash, bool):
            raise BridgeError('compute_hash must be boolean')
        self.compute_hash = compute_hash
        def bind(name, result, *args):
            f = getattr(self.lib, 'paacs_ui_' + name)
            f.restype, f.argtypes = result, list(args)
            setattr(self, '_' + name, f)
        bind('abi_version', U32)
        bind('release_version', C.c_char_p)
        bind('abi_size', SIZE, U32)
        bind('service_create_with_flags', PTR, U32, U32, U64, U64, U32)
        bind('service_destroy', None, PTR)
        bind('service_retained_bytes', U64, PTR)
        bind('request_create', PTR, U32, C.c_char_p)
        bind('request_destroy', None, PTR)
        bind('option_count', SIZE, PTR)
        bind('option_at', U32, PTR, SIZE, C.POINTER(Option))
        bind('request_set', U32, PTR, C.c_char_p, C.c_double)
        bind('request_axial_positions', U32, PTR, C.POINTER(C.c_double), SIZE)
        bind('request_profile_controls', U32, PTR, C.POINTER(ProfileControl), SIZE)
        bind('request_field', U32, PTR, C.POINTER(FieldNode), SIZE, U32, C.c_double, C.c_double, C.c_double)
        bind('submit', U32, PTR, PTR, C.POINTER(PTR), C.c_char_p, SIZE)
        bind('job_release', None, PTR)
        bind('job_snapshot', U32, PTR, C.POINTER(Snapshot))
        bind('job_cancel', U32, PTR, C.POINTER(U32))
        bind('job_wait', U32, PTR, U32, C.POINTER(U32))
        bind('job_result', U32, PTR, C.POINTER(PTR))
        bind('result_release', None, PTR)
        bind('result_summary', U32, PTR, C.POINTER(Summary))
        bind('result_nodes', U32, PTR, U64, C.POINTER(Node), SIZE, C.POINTER(SIZE))
        bind('result_cells', U32, PTR, U64, C.POINTER(Cell), SIZE, C.POINTER(SIZE))
        bind('result_write_archive', U32, PTR, C.c_char_p)
        if self._abi_version() != ABI:
            raise BridgeError('UI ABI version mismatch')
        for i, t in enumerate((Option, Snapshot, Summary, Node, Cell, ProfileControl, FieldNode)):
            if self._abi_size(i) != C.sizeof(t):
                raise BridgeError(f'ABI layout mismatch: {t.__name__}')
        self.service = self._service_create_with_flags(ABI, integer(maximum_jobs, 1),
            integer(maximum_result_bytes, 1, 2**64-1), integer(maximum_retained_bytes, 1, 2**64-1), int(bool(compute_hash)))
        if not self.service:
            raise BridgeError('could not create native job service')

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        if self.service:
            self._service_destroy(self.service)
            self.service = None
        for r in self.results.values():
            self._result_release(r)
        self.results.clear()
        for j in self.jobs.values():
            self._job_release(j)
        self.jobs.clear()

    def _request(self, route):
        if route not in ('primitive', 'adaptive_box'):
            raise BridgeError('unknown route')
        r = self._request_create(ABI, route.encode())
        if not r:
            raise BridgeError('request allocation failed')
        return r

    def schema(self, route):
        r = self._request(route)
        try:
            options = []
            for i in range(self._option_count(r)):
                o = Option()
                check(self._option_at(r, i, C.byref(o)))
                options.append(dict(key=o.key.decode(), kind=('number', 'integer', 'boolean')[o.kind],
                    default=o.default, minimum=o.minimum, maximum=o.maximum,
                    choices=o.choices.decode().split(',') if o.choices else []))
            return dict(protocol=PROTOCOL, abi=ABI, release=self._release_version().decode(),
                route=route, options=options, field_operations=list(FIELD_OPS) if route == 'adaptive_box' else [],
                supports_profile=route == 'primitive', supports_axial_positions=route == 'primitive',
                maximum_page=MAX_PAGE, hash_enabled=self.compute_hash, cancellation='cooperative', progress='stages; generation indeterminate')
        finally:
            self._request_destroy(r)

    def submit(self, spec):
        object_keys(spec, ('route', 'options', 'profile_controls', 'axial_positions', 'field'), ('route',))
        r = self._request(spec['route'])
        try:
            schema = {o['key']: o for o in self.schema(spec['route'])['options']}
            options = spec.get('options', {})
            object_keys(options, schema)
            for key, raw in options.items():
                desc = schema[key]
                if desc['choices'] and isinstance(raw, str):
                    if raw not in desc['choices']:
                        raise BridgeError(f'unknown choice for {key}')
                    v = desc['choices'].index(raw)
                elif desc['kind'] == 'boolean' and isinstance(raw, bool):
                    v = int(raw)
                else:
                    v = number(raw)
                check(self._request_set(r, key.encode(), v), f'invalid option {key}')
            if 'axial_positions' in spec:
                raw = spec['axial_positions']
                if not isinstance(raw, list) or len(raw) > 100001:
                    raise BridgeError('axial_positions needs at most 100001 entries')
                data = (C.c_double * len(raw))(*(number(x) for x in raw))
                check(self._request_axial_positions(r, data, len(raw)))
            if 'profile_controls' in spec:
                raw = spec['profile_controls']
                if not isinstance(raw, list) or len(raw) > 1024:
                    raise BridgeError('profile_controls needs at most 1024 entries')
                data = (ProfileControl * len(raw))()
                for i, control in enumerate(raw):
                    object_keys(control, ('begin', 'end', 'size', 'growth'), ('begin', 'end', 'size', 'growth'))
                    data[i] = ProfileControl(*(number(control[k]) for k in ('begin', 'end', 'size', 'growth')))
                check(self._request_profile_controls(r, data, len(raw)))
                # Presence enables a profile unless explicitly disabled by the user.
                check(self._request_set(r, b'profile_enabled', float(options.get('profile_enabled', True))))
            if 'field' in spec:
                self._set_field(r, spec['field'])
            job, error = PTR(), C.create_string_buffer(192)
            check(self._submit(self.service, r, C.byref(job), error, len(error)), error.value.decode() or 'submission rejected')
            # The expression's error text is evaluated after the native call above.
            s = Snapshot()
            try:
                check(self._job_snapshot(job, C.byref(s)))
                job_id = str(s.id)
                self.jobs[job_id] = job
                return dict(job=job_id)
            except BaseException:
                self._job_release(job)
                raise
        finally:
            self._request_destroy(r)

    def _set_field(self, request, field):
        if field is None:
            check(self._request_field(request, None, 0, 0, 0, 0, 0))
            return
        object_keys(field, ('nodes', 'output', 'minimum_size', 'maximum_size', 'size_factor'), ('nodes', 'output'))
        raw = field['nodes']
        if not isinstance(raw, list) or not 1 <= len(raw) <= 256:
            raise BridgeError('field requires 1..256 nodes')
        nodes = (FieldNode * len(raw))()
        for i, node in enumerate(raw):
            object_keys(node, ('operation', 'first', 'second', 'a', 'b', 'c', 'parameters'), ('operation',))
            if node['operation'] not in FIELD_OPS:
                raise BridgeError('unknown field operation')
            n = nodes[i]
            n.operation = FIELD_OPS.index(node['operation'])
            n.first = integer(node.get('first', 2**32-1))
            n.second = integer(node.get('second', 2**32-1))
            n.a[:] = vector(node.get('a', [0, 0, 0]), 3)
            n.b[:] = vector(node.get('b', [0, 0, 0]), 3)
            n.c[:] = vector(node.get('c', [0, 0, 0]), 3)
            n.parameters[:] = vector(node.get('parameters', [1, 0, 0, 0]), 4)
        check(self._request_field(request, nodes, len(nodes), integer(field['output']),
            number(field.get('minimum_size', .001)), number(field.get('maximum_size', 1)),
            number(field.get('size_factor', 1))))

    def _job(self, job):
        if not isinstance(job, str) or job not in self.jobs:
            raise BridgeError('unknown job handle')
        return self.jobs[job]

    def status(self, job):
        s = Snapshot()
        check(self._job_snapshot(self._job(job), C.byref(s)))
        return dict(job=str(s.id), state=STATES[s.state], stage=STAGES[s.stage],
            completed=str(s.completed), total=str(s.total), elapsed_seconds=s.elapsed_seconds,
            cancellation_requested=bool(s.cancellation_requested), error=s.error, domain=s.domain,
            engine_reason=s.engine_reason, diagnostic=s.diagnostic.decode(),
            retained_result_bytes=str(self._service_retained_bytes(self.service)))

    def cancel(self, job):
        accepted = U32()
        check(self._job_cancel(self._job(job), C.byref(accepted)))
        return dict(accepted=bool(accepted.value))

    def wait(self, job, milliseconds=100):
        done = U32()
        check(self._job_wait(self._job(job), integer(milliseconds, 0, 60_000), C.byref(done)))
        return bool(done.value)

    def _result(self, job):
        self._job(job)
        if job not in self.results:
            r = PTR()
            check(self._job_result(self.jobs[job], C.byref(r)), 'no published result')
            try:
                self.results[job] = r
            except BaseException:
                self._result_release(r)
                raise
        return self.results[job]

    def summary(self, job):
        s = Summary()
        check(self._result_summary(self._result(job), C.byref(s)))
        return dict(job=job, nodes=str(s.nodes), cells=str(s.cells), connectivity=str(s.connectivity),
            logical_bytes=str(s.logical_bytes), hash=f'{s.hash:016x}' if s.hash_available else None, cell_types=[str(x) for x in s.cell_types],
            lower=list(s.lower), upper=list(s.upper), minimum_quality=s.minimum_quality,
            mesh_volume=s.mesh_volume, quality_warning=bool(s.quality_warning))

    def _page(self, job, offset, count, cells):
        if not isinstance(offset, str) or not offset.isascii() or not offset.isdecimal() or len(offset) > 20:
            raise BridgeError('offset must be a decimal uint64 string')
        start = integer(int(offset), 0, 2**64-1)
        n = integer(count, 0, MAX_PAGE)
        data, written = ((Cell if cells else Node) * n)(), SIZE()
        f = self._result_cells if cells else self._result_nodes
        check(f(self._result(job), start, data, n, C.byref(written)))
        return data, written.value

    def page(self, job, offset='0', count=1024, cells=False):
        data, n = self._page(job, offset, count, cells)
        values = []
        for i in range(n):
            v = data[i]
            item = dict(index=str(int(offset)+i), dimension=v.dimension, tag=v.tag)
            if cells:
                item.update(type=v.type, nodes=[str(v.nodes[k]) for k in range(v.node_count)])
            else:
                item['position'] = list(v.position)
            values.append(item)
        return dict(offset=offset, next_offset=str(int(offset)+n), values=values)

    def export(self, job, path):
        if not isinstance(path, str) or not path or '\0' in path:
            raise BridgeError('expected output path without NUL')
        check(self._result_write_archive(self._result(job), str(Path(path).resolve()).encode()), 'archive write failed')
        return dict(path=str(Path(path).resolve()))

    def release(self, job):
        handle = self._job(job)
        if self.status(job)['state'] not in STATES[3:]:
            raise BridgeError('cancel and wait for termination before releasing a job')
        if job in self.results:
            self._result_release(self.results.pop(job))
        self._job_release(handle)
        del self.jobs[job]
        return dict(released=True)

    def dispatch(self, command):
        object_keys(command, ('version', 'id', 'method', 'params'), ('version', 'method'))
        if type(command['version']) is not int or command['version'] != PROTOCOL:
            raise BridgeError('unsupported protocol version')
        if 'id' in command and (not isinstance(command['id'], str) or len(command['id']) > 128):
            raise BridgeError('id must be a string of at most 128 characters')
        method, p = command['method'], command.get('params', {})
        signatures = {
            'schema': (('route',), ('route',)), 'submit': (('spec',), ('spec',)),
            'status': (('job',), ('job',)), 'cancel': (('job',), ('job',)),
            'summary': (('job',), ('job',)), 'release': (('job',), ('job',)),
            'nodes': (('job', 'offset', 'count'), ('job',)),
            'cells': (('job', 'offset', 'count'), ('job',)), 'export': (('job', 'path'), ('job', 'path'))}
        if not isinstance(method, str) or method not in signatures:
            raise BridgeError('unknown method')
        object_keys(p, *signatures[method])
        if method in ('nodes', 'cells'):
            return self.page(**p, cells=method == 'cells')
        return getattr(self, method)(**p)


def parse_json(text):
    def pairs(values):
        out = {}
        for k, v in values:
            if k in out:
                raise BridgeError(f'duplicate key: {k}')
            out[k] = v
        return out
    def bad_constant(value):
        raise BridgeError(f'nonfinite JSON number: {value}')
    return json.loads(text, object_pairs_hook=pairs, parse_constant=bad_constant)


def response(bridge, text):
    command = None
    try:
        command = parse_json(text)
        result = bridge.dispatch(command)
        return dict(version=PROTOCOL, id=command.get('id'), ok=True, result=result)
    except (ValueError, TypeError, OverflowError, RecursionError, MemoryError) as exc:
        ident = command.get('id') if isinstance(command, dict) else None
        if not isinstance(ident, str) or len(ident) > 128:
            ident = None
        return dict(version=PROTOCOL, id=ident, ok=False,
                    error=dict(code=getattr(exc, 'code', 'invalid_request'), message=str(exc)[:512]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--library', required=True)
    parser.add_argument('--hash', action='store_true', help='compute the optional strict result hash')
    args = parser.parse_args()
    with Bridge(args.library, compute_hash=args.hash) as bridge:
        while True:
            line = sys.stdin.buffer.readline(MAX_LINE+1)
            if not line:
                break
            if len(line) > MAX_LINE:
                # Reject and drain exactly this line using bounded reads.
                while line and not line.endswith(b'\n'):
                    line = sys.stdin.buffer.readline(MAX_LINE+1)
                out = dict(version=PROTOCOL, id=None, ok=False,
                           error=dict(code='request_limit', message='JSON line exceeds 1 MiB'))
            else:
                out = response(bridge, line)
            print(json.dumps(out, allow_nan=False, separators=(',', ':')), flush=True)


if __name__ == '__main__':
    main()
