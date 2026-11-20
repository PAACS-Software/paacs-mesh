#!/usr/bin/env python3
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

"""Optional live loopback test. Requires permission to open a local socket."""
import argparse
import json
from pathlib import Path
import re
import sys
import threading
import time
import urllib.error
import urllib.request

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'examples/ui'))
from paacs_bridge import Bridge
from demo_server import make_server


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--library',required=True)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    with Bridge(args.library) as bridge:
        with make_server(bridge) as server:
            thread=threading.Thread(target=server.serve_forever)
            thread.start()
            try:
                base=f'http://127.0.0.1:{server.server_port}'
                with urllib.request.urlopen(base,timeout=10) as r:
                    html=r.read().decode()
                token=re.search("token='([^']+)'",html).group(1)
                def post(command,extra=None):
                    headers={'Content-Type':'application/json','X-PAACS-Token':token}
                    headers.update(extra or {})
                    request=urllib.request.Request(base+'/api',data=json.dumps(command).encode(),headers=headers)
                    with urllib.request.urlopen(request,timeout=10) as r:
                        return json.load(r)
                def api(method,params):
                    r=post(dict(version=1,method=method,params=params))
                    if not r['ok']:raise RuntimeError(r)
                    return r['result']
                assert len(api('schema',{'route':'primitive'})['options'])==49
                assert len(api('schema',{'route':'adaptive_box'})['options'])==19
                for header in [{'X-PAACS-Token':'bad'}, {'Origin':'http://example.invalid'}, {'Host':'example.invalid'}]:
                    try:
                        post(dict(version=1,method='schema',params={'route':'primitive'}),header)
                        raise AssertionError('foreign token/origin/host accepted')
                    except urllib.error.HTTPError as e:
                        assert e.code==403
                j=api('submit',dict(spec=dict(route='primitive',options=dict(elements='tet4',target_size=.4,maximum_chord_error=.03))))['job']
                deadline=time.monotonic()+20
                while time.monotonic()<deadline:
                    status=api('status',dict(job=j))
                    if status['state'] in ('succeeded','failed','cancelled'):break
                    time.sleep(.01)
                assert status['state']=='succeeded',status
                summary=api('summary',dict(job=j));preview=api('preview',dict(job=j))
                assert summary['hash'] is None
                assert len(preview['positions'])==int(summary['nodes'])
                assert len(preview['faces'])==int(summary['cell_types'][2])
                assert not preview['truncated']
                assert all(f['dimension']==2 and f['tag'] in (1,2,3) for f in preview['faces'])
                assert all(0<=i<len(preview['positions']) for f in preview['faces'] for i in f['nodes'])
                try:
                    post(dict(version=1,method='export',params=dict(job=j,path='forbidden.paacsmesh')))
                    raise AssertionError('HTTP export allowed')
                except urllib.error.HTTPError as e:
                    assert e.code==400
                    assert 'native/stdio' in e.read().decode()
                api('release',dict(job=j))
                record=dict(passed=True,browser_rendered=False,nodes=len(preview['positions']),surface_faces=len(preview['faces']),
                    checks=['schema','host_origin_token','live_submission','poll','summary','actual_surface_preview','entity_tags','hash_default_off','http_export_disabled'])
                if args.output:
                    args.output.mkdir(parents=True,exist_ok=True)
                    for name,data in [('http-test',record),('http-summary',summary),('http-preview',preview)]:
                        (args.output/(name+'.json')).write_text(json.dumps(data,indent=2))
                print(json.dumps(record))
            finally:
                server.shutdown();thread.join()


if __name__=='__main__':
    main()
